----------------------------------
-- Switching model of a half-bridge LLC resonant DC-DC converter, 400 V -> 51 V.
--
--   Vdc -+-[Qh]-+          Lr        Cr        ideal xfmr 1 : n     +--[D1]--+- + Vout -+- load
--        |      | sw o--~~~~~~--||------+----------)||(----------+  |        |          |
--       [Ql]    |                       |        primary  sec    +--[D2]--+  Cout       R
--        |      |                     [Lm]                          center-tapped      |
--   0V --+------+---------------------- 0V ---------------------------- 0V -------------+
--
-- Series-resonant tank Lr, Cr with the magnetizing inductance Lm across the
-- transformer primary. The half-bridge runs a fixed 50 % square wave; the
-- output is regulated by *frequency* (a PI loop moves f_sw around the
-- series resonance f_r = 1/(2*pi*sqrt(Lr*Cr))).
--
-- Rectifier :
--   conducting : the secondary clamps the primary to +/- n*Vout. i_Lr
--                resonates (Lr-Cr), i_Lm ramps linearly, the output gets
--                n*|i_Lr - i_Lm|. Diode D1 while i_Lr > i_Lm, D2 while
--                i_Lr < i_Lm.
--   freewheel  : below resonance i_Lr falls to meet i_Lm before the
--                half-period ends; the rectifier then blocks, the primary
--                is unclamped and Lr + Lm resonate together with Cr while
--                the output coasts on Cout. The turn-off instant is found
--                by a short bisection of the half-period (like the boost
--                PFC's DCM search), after which i_Lm is locked to i_Lr.
--
-- state vector : (0 => i_Lr resonant inductor current
--                 1 => i_Lm magnetizing current
--                 2 => v_Cr resonant capacitor voltage (DC bias = Vdc/2)
--                 3 => Vout output capacitor voltage)
----------------------------------
LIBRARY ieee  ;
    USE ieee.NUMERIC_STD.all  ;
    USE ieee.std_logic_1164.all  ;
    use ieee.math_real.all;
    use std.textio.all;

library vunit_lib;
context vunit_lib.vunit_context;

    LIBRARY ode;
    use ode.write_pkg.all;
    use ode.ode_pkg.all;

entity llc_converter_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of llc_converter_tb is

    constant clock_period : time := 1 ns;

    signal simulator_clock    : std_logic := '0';
    signal simulation_counter : natural   := 0;
    -----------------------------------
    -- simulation specific signals ----

    constant n_segments : natural := 2;   -- Qh half, Ql half

    constant minimum_time_step : real := 1.0e-9;   -- seconds per integer time tick

    function to_seconds(ticks : integer) return real is
    begin
        return real(ticks) * minimum_time_step;
    end function;

    function to_ticks(seconds : real) return natural is
    begin
        return natural(round(seconds / minimum_time_step));
    end function;

    constant stoptime       : real    := 4.0e-3;
    constant stoptime_ticks : natural := to_ticks(stoptime);

    signal realtime_ticks : natural := 0;

    constant vdc     : real := 400.0;      -- input DC-link voltage
    constant n_turns : real := 4.2;        -- transformer turns ratio
    constant vref    : real := 51.0;       -- output setpoint

    -- frequency-control loop (PI on Vout -> switching frequency). The tank
    -- (n, Lr, Cr) is sized so the loop settles around 0.78*f_r : the
    -- converter runs *below* resonance, so the resonant half-cycle finishes
    -- before the switching half-period ends and the secondary current sits
    -- at zero (rectifier freewheeling) for the rest of the half - the LLC
    -- analogue of the boost PFC's discontinuous conduction.
    constant f_centre : real := 100.0e3;
    constant f_min    : real := 55.0e3;
    constant f_max    : real := 200.0e3;
    constant kf_p     : real := 1.5e3;     -- Hz per V
    constant kf_i     : real := 1.6e6;     -- Hz per V-s

    -- modulator -> stimulus : the half-period length for this cycle
    signal llc_half_ticks : natural := to_ticks(1.0/(2.0*f_centre));
    signal vout_meas      : real    := vref;

    signal sim_ready        : boolean := false;
    signal modulation_ready : boolean := false;

    ----------------------
    function hb_modulator(sw_state : bit) return real is
    begin
        if sw_state = '1' then
            return 1.0;
        else
            return 0.0;
        end if;
    end hb_modulator;
    ----------------------

begin

------------------------------------------------------------------------
    simtime : process
    begin
        test_runner_setup(runner, runner_cfg);
        wait until realtime_ticks >= stoptime_ticks;
        test_runner_cleanup(runner); -- Simulation ends here
        wait;
    end process simtime;

    simulator_clock <= not simulator_clock after clock_period/2.0;
------------------------------------------------------------------------

    stimulus : process(simulator_clock)

        constant lr    : real := 40.0e-6;    -- series resonant inductor
        constant cr    : real := 38.0e-9;    -- series resonant capacitor  (f_r ~ 129 kHz)
        constant lm    : real := 400.0e-6;   -- magnetizing inductance
        constant r_tank: real := 0.3;        -- tank loss (Lr ESR + winding)
        constant r_mag : real := 0.5;        -- bleeds the magnetizing DC offset (tau = Lm/r_mag)
        constant cout  : real := 470.0e-6;   -- output capacitor
        variable rload : real := 8.5;        -- output load resistance (steps to 6.5 Ohm)

        -- rk5 sub-steps per half switching period (resolve the resonant arc)
        constant steps_per_half : positive := 12;
        constant fw_substeps    : positive := 4;    -- rk5 steps for the freewheel tail
        constant fw_iters       : natural  := 12;   -- freewheel-search bisections
        constant i_rect_thr     : real     := 0.01; -- secondary current at rectifier turn-off

        -- (i_Lr, i_Lm, v_Cr, Vout). v_Cr starts discharged (it charges to its
        -- Vdc/2 DC bias over the first few cycles); Vout pre-charged.
        constant init_state_vector : real_vector(0 to 3) :=
            (0 => 0.0, 1 => 0.0, 2 => 0.0, 3 => vref);

        constant sw_seq : bit_vector(0 to n_segments-1) := "10";

        variable sw_state           : bit := '0';
        variable simulation_started : boolean := false;

        variable seg_index : natural range 0 to n_segments-1 := 0;
        variable half_ticks : natural := 0;
        variable rect_mode  : integer range -1 to 1 := 1;   -- +1 D1, -1 D2, 0 freewheel

        impure function deriv_llc(t : real; states : real_vector) return real_vector is
            variable retval : states'subtype := (others => 0.0);
            variable v_sw   : real;
            variable v_pri  : real;
            variable i_sec  : real := 0.0;
            variable di_fw  : real;
            alias i_Lr is states(0);
            alias i_Lm is states(1);
            alias v_Cr is states(2);
            alias vout is states(3);
        begin
            v_sw := hb_modulator(sw_state) * vdc;   -- 0 or Vdc

            if rect_mode = 0 then
                -- freewheel : Lr + Lm in series carry the same current
                di_fw     := (v_sw - v_Cr - i_Lr*(r_tank + r_mag)) / (lr + lm);
                retval(0) := di_fw;
                retval(1) := di_fw;
            else
                v_pri     := real(rect_mode) * n_turns * vout;
                retval(0) := (v_sw - v_Cr - v_pri - i_Lr*r_tank) / lr;
                retval(1) := (v_pri - i_Lm*r_mag) / lm;
                i_sec     := real(rect_mode) * n_turns * (i_Lr - i_Lm);   -- >= 0 while conducting
            end if;

            retval(2) := i_Lr / cr;
            retval(3) := (i_sec - vout/rload) / cout;

            return retval;
        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_llc);

        variable llc_rk5 : init_state_vector'subtype := init_state_vector;

        file file_handler : text open write_mode is "llc_converter_tb.dat";

        variable now_ticks  : natural := 0;
        variable trial      : init_state_vector'subtype := init_state_vector;
        variable f_lo, f_hi, f_mid : real := 0.0;
        variable cond_ticks : natural := 0;

        ----------------------------------------------------------------
        -- silent forward integration of `st` by `dticks` ticks in `nsub`
        -- equal rk5 steps (uses the current rect_mode)
        procedure run(variable st : inout real_vector;
                      t0_ticks : natural; dticks : natural; nsub : positive) is
            variable hs : real;
        begin
            if dticks < 1 then return; end if;
            hs := to_seconds(dticks) / real(nsub);
            for k in 0 to nsub-1 loop
                rk5(to_seconds(t0_ticks) + real(k)*hs, st, hs);
            end loop;
        end procedure;
        ----------------------------------------------------------------
        procedure log_row(tn : real) is
        begin
            write_to(file_handler,(tn
                    , llc_rk5(0)                                      -- T_i0 : resonant current i_Lr
                    , llc_rk5(1)                                      -- T_i1 : magnetizing current i_Lm
                    , n_turns * abs(llc_rk5(0) - llc_rk5(1))          -- T_i2 : rectified output current
                    , llc_rk5(2)                                      -- B_u0 : resonant cap voltage v_Cr
                    , hb_modulator(sw_state) * vdc                    -- B_u1 : half-bridge node
                    , llc_rk5(3)                                      -- B_u2 : output voltage
                    , vref                                           -- B_u3 : setpoint
                ));
        end procedure;
        ----------------------------------------------------------------
        -- forward integration of the real state, matching run's stepping,
        -- logging one row per step and advancing now_ticks by `dticks`
        procedure run_logged(dticks : natural; nsub : positive) is
            variable hs : real;
            variable t0 : real := to_seconds(now_ticks);
        begin
            if dticks < 1 then return; end if;
            hs := to_seconds(dticks) / real(nsub);
            for k in 0 to nsub-1 loop
                log_row(t0 + real(k)*hs);
                rk5(t0 + real(k)*hs, llc_rk5, hs);
            end loop;
            now_ticks := now_ticks + dticks;
        end procedure;
        ----------------------------------------------------------------
        -- signed secondary current for the current half : > 0 while the
        -- half's diode conducts, <= 0 once the rectifier turns off
        function i_sec_signed(st : real_vector; hm : integer) return real is
        begin
            return real(hm) * (st(0) - st(1));
        end function;
        ----------------------------------------------------------------

    begin
        if rising_edge(simulator_clock) then
            sim_ready <= false;

            if not simulation_started then
                write_plot_config(file_handler, "title", "Half-bridge LLC resonant converter, 400 V -> 51 V");
                write_plot_config(file_handler, "T_title", "Tank and rectifier currents");
                write_plot_config(file_handler, "T_ylabel", "Current [A]");
                write_plot_config(file_handler, "B_title", "Resonant cap, half-bridge node, output");
                write_plot_config(file_handler, "B_ylabel", "Voltage [V]");
                write_plot_config(file_handler, "label_T_i0", "Resonant current i_Lr");
                write_plot_config(file_handler, "label_T_i1", "Magnetizing current i_Lm");
                write_plot_config(file_handler, "label_T_i2", "Rectified output current (pri-referred)");
                write_plot_config(file_handler, "label_B_u0", "Resonant cap voltage v_Cr");
                write_plot_config(file_handler, "label_B_u1", "Half-bridge node voltage");
                write_plot_config(file_handler, "label_B_u2", "Output voltage");
                write_plot_config(file_handler, "label_B_u3", "Output setpoint");
                write_plot_config(file_handler, "drawstyle_B_u1", "steps-pre");

                init_simfile(file_handler, ("time"
                ,"T_i0" ,"T_i1" ,"T_i2"
                ,"B_u0" ,"B_u1" ,"B_u2" ,"B_u3"
                ));

                sim_ready          <= true;
                simulation_started := true;

        elsif modulation_ready then
            simulation_counter <= simulation_counter + 1;

            if seg_index = 0 then
                half_ticks := llc_half_ticks;
                vout_meas  <= llc_rk5(3);
            end if;

            sw_state := sw_seq(seg_index);

            if to_seconds(realtime_ticks) >= 2.0e-3 then
                rload := 6.5;
            end if;

            now_ticks := realtime_ticks;

            -- this half-period's rectifier diode : D1 for the Qh half, D2 for Ql
            if seg_index = 0 then rect_mode := 1; else rect_mode := -1; end if;

            -- test-integrate the whole half : does the rectifier turn off?
            trial := llc_rk5;
            run(trial, now_ticks, half_ticks, steps_per_half);

            if i_sec_signed(trial, rect_mode) > i_rect_thr then
                -- rectifier conducts for the whole half (at / above resonance)
                run_logged(half_ticks, steps_per_half);
            else
                -- below resonance : bisect for the rectifier turn-off instant
                f_lo := 0.0;
                f_hi := 1.0;
                for iter in 1 to fw_iters loop
                    f_mid := 0.5*(f_lo + f_hi);
                    trial := llc_rk5;
                    run(trial, now_ticks, natural(round(real(half_ticks)*f_mid)), steps_per_half);
                    if i_sec_signed(trial, rect_mode) > i_rect_thr then
                        f_lo := f_mid;
                    else
                        f_hi := f_mid;
                    end if;
                end loop;

                cond_ticks := natural(round(real(half_ticks) * f_lo));
                if cond_ticks < 1          then cond_ticks := 1;          end if;
                if cond_ticks > half_ticks then cond_ticks := half_ticks; end if;

                -- power-delivery portion, up to the turn-off
                run_logged(cond_ticks, steps_per_half);

                -- freewheel : lock i_Lm to i_Lr, primary unclamped, output coasts
                llc_rk5(1) := llc_rk5(0);
                rect_mode  := 0;
                if half_ticks - cond_ticks > 0 then
                    run_logged(half_ticks - cond_ticks, fw_substeps);
                end if;
            end if;

            realtime_ticks <= now_ticks;

            if seg_index < n_segments-1 then
                seg_index := seg_index + 1;
            else
                seg_index := 0;
            end if;

            sim_ready <= true;

            end if; -- handshake
        end if; -- rising_edge
    end process stimulus;
------------------------------------------------------------------------
    -- Frequency modulator : PI on the output voltage moves f_sw around the
    -- series resonance (lower f_sw = more gain).
    p_modulation : process(simulator_clock)
        variable t, dt, t_prev : real := 0.0;
        variable verr, f_sw : real := 0.0;
        variable f_integ : real := 0.0;
    begin
        if rising_edge(simulator_clock) then
            modulation_ready <= false;

            if sim_ready then
                t      := to_seconds(realtime_ticks);
                dt     := t - t_prev;
                t_prev := t;

                verr := vref - vout_meas;
                f_sw := f_centre - kf_p*verr - f_integ;
                if f_sw > f_min and f_sw < f_max then
                    f_integ := f_integ + kf_i*verr*dt;   -- conditional integration
                end if;
                f_sw := f_centre - kf_p*verr - f_integ;
                if    f_sw < f_min then f_sw := f_min;
                elsif f_sw > f_max then f_sw := f_max;
                end if;

                llc_half_ticks   <= to_ticks(1.0/(2.0*f_sw));
                modulation_ready <= true;
            end if;
        end if;
    end process p_modulation;
------------------------------------------------------------------------
end vunit_simulation;
