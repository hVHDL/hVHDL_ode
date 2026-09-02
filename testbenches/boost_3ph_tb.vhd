----------------------------------
-- Switching model of three interleaved (parallel) synchronous boost
-- converters sharing one input source and one output capacitor / load.
--
-- Same technique as boost_converter_tb, but the three phases switch
-- independently, 120 deg apart, so the simulation can no longer step by
-- whole on/off intervals. Instead it advances to the next switching edge
-- among all phases, integrates that segment with a fixed-step rk5, then
-- toggles whichever phase(s) are due.
--
--        Vin o--+--L1--+--[ Q2a ]--+----------+---> Vout
--               |      |           |          |
--               |    [ Q1a ]       |          |
--               |      |           |          |
--               +--L2--+--[ Q2b ]--+          C   Rload
--               |      |           |          |
--               |    [ Q1b ]       |          |
--               |      |           |          |
--               +--L3--+--[ Q2c ]--+----------+---> GND
--               |      |           |
--               |    [ Q1c ]       |
--        GND o--+------+-----------+
--
-- Q1x on (sw_state(x) = '0') : phase x inductor charges from Vin.
-- Q2x on (sw_state(x) = '1') : phase x inductor drives the shared output.
--
-- state vector : (0 => iL1 , 1 => iL2 , 2 => iL3 , 3 => vC)
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

entity boost_3ph_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of boost_3ph_tb is

    constant clock_period : time := 1 ns;

    signal simulator_clock    : std_logic := '0';
    signal simulation_counter : natural   := 0;
    -----------------------------------
    -- simulation specific signals ----

    constant n_phases : natural := 3;

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

    constant stoptime       : real    := 4.0e-3;
    constant stoptime_ticks : natural := to_ticks(stoptime);

    signal realtime_ticks : natural := 0;

    ----------------------
    -- half-bridge modulator : 1.0 while a phase's switch node is tied to
    -- the high rail (high-side Q2x conducting, that inductor feeding the
    -- output), 0.0 while it is tied to GND (low-side Q1x conducting, that
    -- inductor charging). Multiplying a rail voltage or a branch current by
    -- this gives the switched quantity, e.g.
    --   v_sw  = hb_modulator(sw_state) * vc
    --   i_out = hb_modulator(sw_state) * il
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

        variable udc    : real := 12.0;      -- input voltage
        constant l      : real := 8.0e-6;   -- per-phase inductance
        constant c      : real := 100.0e-6;  -- shared output capacitance
        constant rl     : real := 15.0e-3;   -- per-phase inductor series resistance
        variable rload  : real := 8.0;       -- output load resistance

        variable sw_frequency : real := 100.0e3;
        variable t_sw : real := 1.0/sw_frequency;
        variable duty : real := 0.5;         -- low-side (charge) on-time fraction

        -- each inter-edge segment is integrated as this many equal rk5 steps
        constant steps_per_segment : positive := 1;

        variable seed1, seed2 : positive := 1;
        variable rand : real;

        -- iL1, iL2, iL3, vC. vC starts at Vin (output cap pre-charged
        -- through the Q2x body diodes before switching begins).
        constant init_state_vector : real_vector :=
            (0 => 0.0, 1 => 0.0, 2 => 0.0, n_phases => 12.0);

        type sw_array   is array (0 to n_phases-1) of bit;
        type tick_array is array (0 to n_phases-1) of natural;

        -- '0' -> charging (Q1x on), '1' -> transferring (Q2x on)
        variable sw_state  : sw_array   := (others => '0');
        -- absolute tick at which each phase next toggles
        variable next_edge : tick_array := (others => 0);
        variable phases_initialised : boolean := false;

        impure function deriv_boost(t : real; states : real_vector) return real_vector is
            variable retval      : states'subtype := (others => 0.0);
            variable cap_current : real := 0.0;
            alias vc is states(n_phases);
        begin
            for k in 0 to n_phases-1 loop
                -- states(k) is phase k inductor current
                retval(k)   := (udc - states(k) * rl - hb_modulator(sw_state(k)) * vc) * (1.0/l);
                cap_current := cap_current + hb_modulator(sw_state(k)) * states(k);
            end loop;
            retval(n_phases) := (cap_current - vc/rload) * (1.0/c);

            return retval;
        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_boost);

        variable boost_rk5 : init_state_vector'subtype := init_state_vector;

        file file_handler : text open write_mode is "boost_3ph_tb.dat";

        variable charge_ticks   : natural := 0;
        variable transfer_ticks : natural := 0;
        variable period_ticks   : natural := 0;
        variable offset         : natural := 0;
        variable target_ticks   : natural := 0;
        variable segment_ticks  : natural := 0;
        variable step_ticks     : natural := 1;
        variable now_ticks      : natural := 0;
        variable h              : real    := 0.0;   -- rk5 step size, seconds

    begin
        if rising_edge(simulator_clock) then
            if simulation_counter = 0 then
                write_plot_config(file_handler, "title", "Three interleaved boost converters");
                write_plot_config(file_handler, "T_title", "Phase and total input currents");
                write_plot_config(file_handler, "T_ylabel", "Current [A]");
                write_plot_config(file_handler, "B_title", "Output and input voltage");
                write_plot_config(file_handler, "B_ylabel", "Voltage [V]");
                write_plot_config(file_handler, "label_T_i0", "Phase 1 inductor current");
                write_plot_config(file_handler, "label_T_i1", "Phase 2 inductor current");
                write_plot_config(file_handler, "label_T_i2", "Phase 3 inductor current");
                write_plot_config(file_handler, "label_T_i3", "Total input current");
                write_plot_config(file_handler, "label_B_u0", "Output voltage");
                write_plot_config(file_handler, "label_B_u1", "Input voltage");
                write_plot_config(file_handler, "label_B_u2", "Phase 1 switch-node voltage");
                -- keep the switch node square with one rk5 step per segment
                write_plot_config(file_handler, "drawstyle_B_u2", "steps-pre");

                init_simfile(file_handler, ("time"
                ,"T_i0"
                ,"T_i1"
                ,"T_i2"
                ,"T_i3"
                ,"B_u0"
                ,"B_u1"
                ,"B_u2"
                ));
            end if;
            simulation_counter <= simulation_counter + 1;

            -------------------------
            -- small random dither on the duty cycle, same trick as
            -- fc_4level_tb, keeps the solver from locking onto a perfectly
            -- periodic orbit and broadens the spectral content.
            uniform(seed1, seed2, rand);
            duty := 0.5 + ((rand - 0.5) * 2.0) * 0.002;

            charge_ticks   := to_ticks(t_sw * duty);
            transfer_ticks := to_ticks(t_sw * (1.0 - duty));
            period_ticks   := charge_ticks + transfer_ticks;

            -- load step half way through the run
            if to_seconds(realtime_ticks) >= 2.0e-3 then
                rload := 4.0;
            end if;

            -- stagger the phases 360/n_phases degrees apart on the first pass
            if not phases_initialised then
                for k in 0 to n_phases-1 loop
                    offset := (k * period_ticks) / n_phases;
                    if offset < charge_ticks then
                        sw_state(k)  := '0';
                        next_edge(k) := charge_ticks - offset;
                    else
                        sw_state(k)  := '1';
                        next_edge(k) := period_ticks - offset;
                    end if;
                end loop;
                phases_initialised := true;
            end if;

            now_ticks := realtime_ticks;

            -- advance only as far as the next switching edge among all phases
            target_ticks := next_edge(0);
            for k in 1 to n_phases-1 loop
                if next_edge(k) < target_ticks then
                    target_ticks := next_edge(k);
                end if;
            end loop;

            segment_ticks := target_ticks - now_ticks;
            if segment_ticks < 1 then
                segment_ticks := 1;   -- guard against coincident edges
            end if;

            step_ticks := segment_ticks / steps_per_segment;
            if step_ticks < 1 then
                step_ticks := 1;
            end if;
            h := to_seconds(step_ticks);

            for s in 1 to steps_per_segment loop
                write_to(file_handler,(to_seconds(now_ticks)
                        , boost_rk5(0)                                          -- T_i0 : phase 1 current
                        , boost_rk5(1)                                          -- T_i1 : phase 2 current
                        , boost_rk5(2)                                          -- T_i2 : phase 3 current
                        , boost_rk5(0) + boost_rk5(1) + boost_rk5(2)            -- T_i3 : total input current
                        , boost_rk5(n_phases)                                   -- B_u0 : output voltage
                        , udc                                                  -- B_u1 : input voltage
                        , hb_modulator(sw_state(0)) * boost_rk5(n_phases)       -- B_u2 : phase 1 switch node
                    ));

                rk5(to_seconds(now_ticks), boost_rk5, h);
                now_ticks := now_ticks + step_ticks;
            end loop;

            -- toggle every phase whose edge has now been reached
            for k in 0 to n_phases-1 loop
                if now_ticks >= next_edge(k) then
                    if sw_state(k) = '0' then
                        sw_state(k)  := '1';
                        next_edge(k) := next_edge(k) + transfer_ticks;
                    else
                        sw_state(k)  := '0';
                        next_edge(k) := next_edge(k) + charge_ticks;
                    end if;
                end if;
            end loop;

            realtime_ticks <= now_ticks;

        end if; -- rising_edge
    end process stimulus;
------------------------------------------------------------------------
end vunit_simulation;
