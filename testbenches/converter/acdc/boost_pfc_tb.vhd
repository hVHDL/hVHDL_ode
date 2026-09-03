----------------------------------
-- Switching model of a single-phase boost power-factor-correction (PFC)
-- front end : a diode-bridge rectifier feeding a boost stage whose inductor
-- current is shaped to follow the rectified line voltage (unity power
-- factor) while a slow outer loop regulates the DC output. Same simulation
-- technique as the other converter testbenches : the power stage is a set
-- of ODEs integrated with a fixed-step rk5, and each switching period is
-- walked as two sub-intervals (Q on, Q off) whose lengths come from the
-- modulator.
--
--   v_ac ~ --[bridge]-- v_rect --L,RL--+--[D]--+-- + Vout --+-- load
--                                      |       |            |
--                                     [Q]     Cout           R
--                                      |       |             |
--   0V ---------------------------------+-------+---- 0V -----+
--
--   Q on  : diL/dt = (v_rect - iL*RL)/L                  (inductor charges)
--   Q off : diL/dt = (v_rect - Vout - iL*RL)/L , iL > 0  (boost diode conducts)
--           diL/dt = 0                          , iL <= 0 (DCM, diode blocks)
--
-- Discontinuous conduction (light load, near the line zero crossings) is
-- resolved rather than clamped : when a Q-off sub-interval would drive iL
-- below 0, a short bisection of the off-time (dcm_iters trial integrations
-- of a state copy) locates the negative-slope length, the real state is
-- integrated only up to that zero crossing, and iL is then held at 0
-- (vL = 0, Vout decaying into the load) for the rest of the sub-interval.
--
-- Modulation (p_modulation), two nested loops :
--   outer : slow PI on (Vref - LPF(Vout)) -> line-current amplitude i_amp
--           (the low-pass keeps the double-line-frequency ripple out of the
--            current shape)
--   inner : per switching period, PI current loop trimming the boost feed-
--           forward duty :   d = (1 - v_rect/Vout) + PI(i_ref - iL)
--           with i_ref = i_amp * |sin(w_line * t)|  (rectified-sine shape)
--   sub-interval lengths are d*Tsw (Q on) and (1-d)*Tsw (Q off).
--
-- state vector : (0 => iL   boost inductor current
--                 1 => Vout output capacitor voltage)
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

entity boost_pfc_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of boost_pfc_tb is

    constant clock_period : time := 1 ns;

    signal simulator_clock    : std_logic := '0';
    signal simulation_counter : natural   := 0;
    -----------------------------------
    -- simulation specific signals ----

    constant n_segments : natural := 2;   -- Q on, Q off

    -- Time is kept as an integer number of ticks. One tick lasts
    -- minimum_time_step seconds; convert to real seconds only where one is
    -- genuinely needed (the rk5 step size and the plot file).
    constant minimum_time_step : real := 10.0e-9;   -- seconds per integer time tick

    function to_seconds(ticks : integer) return real is
    begin
        return real(ticks) * minimum_time_step;
    end function;

    function to_ticks(seconds : real) return natural is
    begin
        return natural(round(seconds / minimum_time_step));
    end function;

    constant stoptime       : real    := 200.0e-3;  -- ten line cycles
    constant stoptime_ticks : natural := to_ticks(stoptime);

    signal realtime_ticks : natural := 0;

    constant sw_frequency : real := 65.0e3;
    constant t_sw         : real := 1.0/sw_frequency;

    constant f_line : real := 50.0;                 -- mains frequency [Hz]
    constant v_pk   : real := 230.0 * sqrt(2.0);    -- mains peak voltage (230 Vrms)
    constant vref   : real := 400.0;                -- DC output setpoint

    -- control gains
    constant kv_p     : real := 0.10;     -- voltage loop P  (A per V)
    constant kv_i     : real := 4.0;      -- voltage loop I  (A per V-s)
    constant i_amp_max: real := 15.0;     -- line-current amplitude clamp [A]
    constant kip      : real := 0.05;     -- current loop P  (duty per A)
    constant kii      : real := 200.0;    -- current loop I  (duty per A-s)

    -- modulator -> stimulus : the two sub-interval lengths for this period
    type seg_tick_array is array (0 to n_segments-1) of natural;
    signal pfc_seg_ticks : seg_tick_array := (others => 0);
    signal il_meas       : real := 0.0;   -- inductor current, period-sampled
    signal vout_meas     : real := vref;  -- output voltage,  period-sampled
    signal i_amp_cmd     : real := 0.0;   -- line-current amplitude command (for logging)

    -- Handshake : stimulus pulses sim_ready to ask for a fresh modulator
    -- solution; p_modulation pulses modulation_ready; stimulus then
    -- integrates the next sub-interval.
    signal sim_ready        : boolean := false;
    signal modulation_ready : boolean := false;

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

        constant l     : real := 1.0e-3;     -- boost inductor
        constant rl    : real := 100.0e-3;   -- inductor series resistance
        constant cout  : real := 470.0e-6;   -- output capacitor
        variable rload : real := 940.0;      -- output load resistance (steps to 190 Ohm)

        -- each sub-interval is integrated as this many equal rk5 steps
        constant steps_per_segment : positive := 2;

        -- discontinuous-conduction handling : when the Q-off inductor current
        -- would reach 0 mid sub-interval, bisect the off-time (dcm_iters
        -- trial integrations) to locate the zero crossing, integrate the real
        -- state only up to it, then hold iL = 0 (vL = 0) for the remainder.
        constant dcm_enable : boolean := true;
        constant dcm_iters  : natural := 12;

        -- (iL, Vout). Output cap starts pre-charged to the rectified peak
        -- (the value it reaches through the bridge + boost diode before the
        -- boost stage starts switching).
        constant init_state_vector : real_vector(0 to 1) := (0 => 0.0, 1 => 325.0);

        -- the "switch matrix" : Q on for sub-interval 0, off for 1
        constant sw_seq : bit_vector(0 to n_segments-1) := "10";

        variable sw_state           : bit := '0';
        variable simulation_started : boolean := false;

        variable seg_index : natural range 0 to n_segments-1 := 0;
        variable cur_seg   : seg_tick_array := (others => 0);

        impure function q_on return boolean is
        begin return sw_state = '1'; end function;

        impure function deriv_pfc(t : real; states : real_vector) return real_vector is
            variable retval  : states'subtype := (others => 0.0);
            variable v_rect  : real;
            variable i_diode : real := 0.0;
            alias iL   is states(0);
            alias vout is states(1);
        begin
            v_rect := abs(v_pk * sin(2.0*MATH_PI*f_line*t));

            if q_on then
                retval(0) := (v_rect - iL*rl) / l;                 -- inductor charges
            elsif iL > 0.0 then
                retval(0) := (v_rect - vout - iL*rl) / l;          -- boost diode conducts
                i_diode   := iL;
            else
                retval(0) := 0.0;                                  -- DCM : diode blocks
            end if;

            retval(1) := (i_diode - vout/rload) / cout;            -- d(Vout)/dt

            return retval;
        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_pfc);

        variable pfc_rk5 : init_state_vector'subtype := init_state_vector;

        file file_handler : text open write_mode is "boost_pfc_tb.dat";

        variable seg_ticks  : natural := 0;
        variable now_ticks  : natural := 0;

        -- DCM zero-crossing search
        variable trial      : init_state_vector'subtype := init_state_vector;
        variable f_lo, f_hi, f_mid : real := 0.0;
        variable cond_ticks, idle_ticks : natural := 0;

        ----------------------------------------------------------------
        -- silent forward integration of `st` by `dticks` ticks in `nsub`
        -- equal rk5 steps, starting the line-voltage clock at t0_ticks
        -- (used for the DCM trial integrations)
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
        -- forward integration of the real state (pfc_rk5), advancing
        -- now_ticks by exactly `dticks` and logging one row per step
        procedure run_logged(dticks : natural; nsub : positive) is
            variable step : natural := 1;
            variable remaining  : natural := dticks;
            variable tn   : real;
        begin
            if dticks < 1 then return; end if;
            step := dticks / nsub;
            if step < 1 then step := 1; end if;
            while remaining > 0 loop
                if remaining < step then step := remaining; end if;
                tn := to_seconds(now_ticks);
                write_to(file_handler,(tn
                        , abs(v_pk * sin(2.0*MATH_PI*f_line*tn))                 -- T_i0 : rectified line voltage
                        , pfc_rk5(0) * 40.0                                      -- T_i1 : inductor current x40
                        , i_amp_cmd * abs(sin(2.0*MATH_PI*f_line*tn)) * 40.0     -- T_i2 : current reference x40
                        , pfc_rk5(1)                                             -- B_u0 : output voltage
                        , vref                                                  -- B_u1 : setpoint
                    ));
                rk5(tn, pfc_rk5, to_seconds(step));
                now_ticks := now_ticks + step;
                remaining := remaining - step;
            end loop;
        end procedure;
        ----------------------------------------------------------------

    begin
        if rising_edge(simulator_clock) then
            sim_ready <= false;   -- default: sim_ready is a single-cycle pulse

            if not simulation_started then
                write_plot_config(file_handler, "title", "Single-phase boost PFC");
                write_plot_config(file_handler, "T_title", "Line voltage and inductor current");
                write_plot_config(file_handler, "T_ylabel", "V [V] , A*40");
                write_plot_config(file_handler, "B_title", "DC output voltage");
                write_plot_config(file_handler, "B_ylabel", "Voltage [V]");
                write_plot_config(file_handler, "label_T_i0", "Rectified line voltage");
                write_plot_config(file_handler, "label_T_i1", "Inductor current x 40");
                write_plot_config(file_handler, "label_T_i2", "Shaped current reference x 40");
                write_plot_config(file_handler, "label_B_u0", "Output voltage");
                write_plot_config(file_handler, "label_B_u1", "Output setpoint");

                init_simfile(file_handler, ("time"
                ,"T_i0"
                ,"T_i1"
                ,"T_i2"
                ,"B_u0"
                ,"B_u1"
                ));

                sim_ready          <= true;   -- kick off the first handshake
                simulation_started := true;

        -- one sub-interval per completed modulation handshake
        elsif modulation_ready then
            simulation_counter <= simulation_counter + 1;

            -- once per period : latch the sub-interval lengths and take a
            -- period-synchronised sample of iL and Vout for the modulator
            if seg_index = 0 then
                cur_seg   := pfc_seg_ticks;
                il_meas   <= pfc_rk5(0);
                vout_meas <= pfc_rk5(1);
            end if;

            sw_state := sw_seq(seg_index);

            -- load step at mid-run
            if to_seconds(realtime_ticks) >= 100.0e-3 then
                rload := 190.0;
            end if;

            seg_ticks := cur_seg(seg_index);
            if seg_ticks < 1 then
                seg_ticks := 1;
            end if;
            now_ticks := realtime_ticks;

            if (seg_index = 0) or (not dcm_enable) or (pfc_rk5(0) <= 0.0) then
                -- Q on, or DCM disabled, or the inductor is already idle
                -- (iL = 0, deriv holds it there and Vout decays) : integrate
                -- the whole sub-interval directly.
                run_logged(seg_ticks, steps_per_segment);
                if pfc_rk5(0) < 0.0 then
                    pfc_rk5(0) := 0.0;   -- safety clamp
                end if;

            else
                -- Q off with the boost diode conducting : the inductor
                -- current has a negative slope. Test whether it reaches 0
                -- within the off-time.
                trial := pfc_rk5;
                run(trial, now_ticks, seg_ticks, steps_per_segment);

                if trial(0) > 0.0 then
                    -- stays positive : continuous conduction this period
                    run_logged(seg_ticks, steps_per_segment);
                else
                    -- discontinuous : bisect the off-time for the length of
                    -- the negative-slope (diode-conducting) sub-interval
                    f_lo := 0.0;
                    f_hi := 1.0;
                    for iter in 1 to dcm_iters loop
                        f_mid := 0.5*(f_lo + f_hi);
                        trial := pfc_rk5;
                        run(trial, now_ticks,
                            natural(round(real(seg_ticks) * f_mid)), steps_per_segment);
                        if trial(0) > 0.0 then
                            f_lo := f_mid;
                        else
                            f_hi := f_mid;
                        end if;
                    end loop;

                    -- f_lo is the largest tested fraction with iL still > 0
                    cond_ticks := natural(round(real(seg_ticks) * f_lo));
                    if cond_ticks < 1         then cond_ticks := 1;         end if;
                    if cond_ticks > seg_ticks then cond_ticks := seg_ticks; end if;
                    idle_ticks := seg_ticks - cond_ticks;

                    -- diode-conducting portion, up to the zero crossing
                    run_logged(cond_ticks, steps_per_segment);
                    pfc_rk5(0) := 0.0;   -- snap iL to exactly zero at the crossing

                    -- idle portion : vL = 0, iL = 0, Vout decays into the load
                    run_logged(idle_ticks, 1);
                end if;
            end if;

            if pfc_rk5(0) < 0.0 then
                pfc_rk5(0) := 0.0;   -- the bridge + boost diode block reverse iL
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
    -- Boost PFC modulator : a slow outer PI on the DC output sets the line-
    -- current amplitude, an inner per-period law shapes the duty so the
    -- inductor current tracks a rectified sine.
    p_modulation : process(simulator_clock)
        variable t, dt, t_prev : real := 0.0;
        variable v_rect, verr, i_amp, i_ref, ierr, d_corr, d : real := 0.0;
        variable v_integ  : real := 4.1;        -- pre-charge : steady i_amp for the initial load
        variable i_integ  : real := 0.0;        -- current-loop integrator
        variable vout_filt : real := vref;      -- low-passed output voltage
        variable seg      : real_vector(0 to n_segments-1);
        constant d_corr_max : real := 0.35;    -- current-loop duty authority
        constant tau_v      : real := 8.0e-3; -- output-voltage filter time constant
                                               -- (>> 1/2f_line so the loop ignores
                                               --  the double-line-frequency ripple)
    begin
        if rising_edge(simulator_clock) then
            modulation_ready <= false;   -- default: single-cycle pulse

            if sim_ready then
                t      := to_seconds(realtime_ticks);
                dt     := t - t_prev;
                t_prev := t;

                v_rect := abs(v_pk * sin(2.0*MATH_PI*f_line*t));

                -- outer voltage loop -> line-current amplitude, acting on a
                -- low-passed Vout so the double-line-frequency ripple does
                -- not fold straight back into the current amplitude
                if dt > 0.0 then
                    vout_filt := vout_filt + (vout_meas - vout_filt) * (dt / tau_v);
                end if;
                verr := vref - vout_filt;
                i_amp := kv_p*verr + v_integ;
                if i_amp > 0.0 and i_amp < i_amp_max then
                    v_integ := v_integ + kv_i*verr*dt;   -- conditional integration
                end if;
                i_amp := kv_p*verr + v_integ;
                if    i_amp < 0.0        then i_amp := 0.0;
                elsif i_amp > i_amp_max  then i_amp := i_amp_max;
                end if;

                -- rectified-sine current reference, PI current loop trimming
                -- the boost feed-forward duty (1 - v_rect/Vout)
                i_amp_cmd <= i_amp;
                i_ref  := i_amp * abs(sin(2.0*MATH_PI*f_line*t));
                ierr   := i_ref - il_meas;
                d_corr := kip*ierr + i_integ;
                if d_corr > -d_corr_max and d_corr < d_corr_max then
                    i_integ := i_integ + kii*ierr*dt;   -- conditional integration
                end if;
                d_corr := kip*ierr + i_integ;
                if    d_corr >  d_corr_max then d_corr :=  d_corr_max;
                elsif d_corr < -d_corr_max then d_corr := -d_corr_max;
                end if;

                d := (1.0 - v_rect/vout_meas) + d_corr;
                if    d < 0.02 then d := 0.02;
                elsif d > 0.98 then d := 0.98;
                end if;

                seg(0) := d       * t_sw;   -- Q on
                seg(1) := (1.0-d) * t_sw;   -- Q off
                for i in 0 to n_segments-1 loop
                    pfc_seg_ticks(i) <= to_ticks(seg(i));
                end loop;
                modulation_ready <= true;
            end if;
        end if;
    end process p_modulation;
------------------------------------------------------------------------
end vunit_simulation;
