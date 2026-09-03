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
    constant minimum_time_step : real := 1.0e-9;   -- 1 ns per integer time tick

    function to_seconds(ticks : integer) return real is
    begin
        return real(ticks) * minimum_time_step;
    end function;

    function to_ticks(seconds : real) return natural is
    begin
        return natural(round(seconds / minimum_time_step));
    end function;

    constant stoptime       : real    := 100.0e-3;  -- five line cycles
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

        constant l     : real := 2.0e-3;     -- boost inductor
        constant rl    : real := 100.0e-3;   -- inductor series resistance
        constant cout  : real := 470.0e-6;   -- output capacitor
        variable rload : real := 240.0;      -- output load resistance (steps to 190 Ohm)

        -- each sub-interval is integrated as this many equal rk5 steps
        constant steps_per_segment : positive := 2;

        -- (iL, Vout). Output cap starts pre-charged to the rectified peak
        -- (the value it reaches through the bridge + boost diode before the
        -- boost stage starts switching).
        constant init_state_vector : real_vector(0 to 1) := (0 => 0.0, 1 => vref);

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
        variable step_ticks : natural := 1;
        variable now_ticks  : natural := 0;
        variable h          : real    := 0.0;   -- rk5 step size, seconds
        variable t_now      : real    := 0.0;

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
            if to_seconds(realtime_ticks) >= 40.0e-3 then
                rload := 190.0;
            end if;

            seg_ticks := cur_seg(seg_index);
            if seg_ticks < 1 then
                seg_ticks := 1;
            end if;

            step_ticks := seg_ticks / steps_per_segment;
            if step_ticks < 1 then
                step_ticks := 1;
            end if;
            h := to_seconds(step_ticks);
            now_ticks := realtime_ticks;

            for s in 1 to steps_per_segment loop
                t_now := to_seconds(now_ticks);
                write_to(file_handler,(t_now
                        , abs(v_pk * sin(2.0*MATH_PI*f_line*t_now))    -- T_i0 : rectified line voltage
                        , pfc_rk5(0) * 40.0                            -- T_i1 : inductor current (scaled)
                        , i_amp_cmd * abs(sin(2.0*MATH_PI*f_line*t_now)) * 40.0  -- T_i2 : current reference (scaled)
                        , pfc_rk5(1)                                   -- B_u0 : output voltage
                        , vref                                        -- B_u1 : setpoint
                    ));

                rk5(t_now, pfc_rk5, h);
                now_ticks := now_ticks + step_ticks;
            end loop;

            -- the boost inductor current cannot go negative (bridge + boost
            -- diode); clamp any zero-crossing overshoot from the fixed step
            if pfc_rk5(0) < 0.0 then
                pfc_rk5(0) := 0.0;
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
