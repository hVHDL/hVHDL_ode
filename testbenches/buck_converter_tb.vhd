----------------------------------
-- Switching model of a synchronous buck converter.
--
-- Same simulation technique as fc_4level_tb : the power stage is a set of
-- ODEs integrated with a fixed-step Dormand-Prince rk5, and the switching
-- is modelled by stepping through the converter's switch states, using a
-- step length equal to that sub-interval's on- or off-time. No averaging,
-- so the inductor current ripple and switch-node waveform are visible.
--
--            Vin o--+--[ Q1 ]--+---L---+----+---> Vout
--                   |          |       |    |
--                            [ Q2 ]  RL(L)  C   Rload
--                   |          |       |    |
--            GND o--+----------+-------+----+---> GND
--
-- state vector : (0 => iL , 1 => vC)
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

entity buck_converter_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of buck_converter_tb is

    constant clock_period : time := 1 ns;

    signal simulator_clock    : std_logic := '0';
    signal simulation_counter : natural   := 0;
    -----------------------------------
    -- simulation specific signals ----

    signal realtime : real := 0.0;
    constant stoptime : real := 4.0e-3;

    ----------------------
    -- half-bridge modulator : 1.0 while the switch node is tied to the
    -- high rail (high-side device conducting), 0.0 while it is tied to GND
    -- (low-side device conducting). Multiplying a rail voltage or a branch
    -- current by this gives the switched quantity, e.g.
    --   v_sw  = hb_modulator(sw_state) * udc
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
        wait until realtime >= stoptime;
        test_runner_cleanup(runner); -- Simulation ends here
        wait;
    end process simtime;

    simulator_clock <= not simulator_clock after clock_period/2.0;
------------------------------------------------------------------------

    stimulus : process(simulator_clock)

        variable udc    : real := 24.0;
        constant l      : real := 47.0e-6;
        constant c      : real := 100.0e-6;
        constant rl     : real := 15.0e-3;   -- inductor series resistance
        variable rload  : real := 6.0;       -- output load resistance

        variable sw_frequency : real := 100.0e3;
        variable t_sw : real := 1.0/sw_frequency;
        variable duty : real := 0.5;

        -- number of fixed rk5 steps used to integrate each conduction
        -- interval. 1 => one step spanning the whole on- (or off-) time,
        -- as fc_4level_tb does. A larger value subdivides the interval into
        -- that many equal rk5 steps, giving a finer waveform and a smaller
        -- per-step truncation error at the cost of more solver calls.
        constant steps_per_on_time  : positive := 1;
        constant steps_per_off_time : positive := 1;

        variable seed1, seed2 : positive := 1;
        variable rand : real;

        -- iL, vC
        constant init_state_vector : real_vector := (0 => 0.0, 1 => 0.0);

        -- '1' -> high-side conduction interval, '0' -> low-side interval
        variable sw_state : bit := '1';

        ----------
        -- total length of the conduction interval about to be integrated
        impure function get_interval_length return real is
        begin
            if sw_state = '1' then
                return t_sw * duty;
            else
                return t_sw * (1.0 - duty);
            end if;
        end get_interval_length;
        ----------
        -- how many rk5 steps that interval is split into
        impure function get_substeps return positive is
        begin
            if sw_state = '1' then
                return steps_per_on_time;
            else
                return steps_per_off_time;
            end if;
        end get_substeps;
        ----------

        impure function deriv_buck(t : real; states : real_vector) return real_vector is
            variable retval : states'subtype := (others => 0.0);
            variable v_sw : real := 0.0;
            alias il is states(0);
            alias vc is states(1);
        begin
            v_sw := hb_modulator(sw_state) * udc;

            retval(0) := (v_sw - il * rl - vc) * (1.0/l);
            retval(1) := (il - vc/rload)       * (1.0/c);

            return retval;
        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_buck);

        variable buck_rk5 : init_state_vector'subtype := init_state_vector;

        file file_handler : text open write_mode is "buck_converter_tb.dat";

        variable interval_length : real := 0.0;
        variable substeps        : positive := 1;
        variable h               : real := 0.0;
        variable t_now           : real := 0.0;

    begin
        if rising_edge(simulator_clock) then
            if simulation_counter = 0 then
                write_plot_config(file_handler, "title", "Synchronous buck converter switching model");
                write_plot_config(file_handler, "T_title", "Inductor current");
                write_plot_config(file_handler, "T_ylabel", "Current [A]");
                write_plot_config(file_handler, "B_title", "Switch-node and output voltage");
                write_plot_config(file_handler, "B_ylabel", "Voltage [V]");
                write_plot_config(file_handler, "label_T_i0", "Inductor current");
                write_plot_config(file_handler, "label_B_u0", "Switch-node voltage");
                write_plot_config(file_handler, "label_B_u1", "Output voltage");
                write_plot_config(file_handler, "label_B_u2", "Input voltage");
                -- draw the switch-node and input voltages as stepped lines so
                -- they stay square even with one rk5 step per conduction
                -- interval (steps_per_on_time = steps_per_off_time = 1).
                write_plot_config(file_handler, "drawstyle_B_u0", "steps-pre");

                init_simfile(file_handler, ("time"
                ,"T_i0"
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
            rand := ((rand - 0.5) * 2.0) * 0.002;
            duty := 0.5 + rand;

            -- load step half way through the run
            if realtime >= 2.0e-3 then
                rload := 3.0;
            end if;

            -- integrate the conduction interval as N equal rk5 steps,
            -- logging a sample before each one.
            interval_length := get_interval_length;
            substeps        := get_substeps;
            h               := interval_length / real(substeps);
            t_now           := realtime;

            for step in 1 to substeps loop
                write_to(file_handler,(t_now
                        , buck_rk5(0)                          -- T_i0 : inductor current
                        , hb_modulator(sw_state) * udc         -- B_u0 : switch-node voltage
                        , buck_rk5(1)                          -- B_u1 : output voltage
                        , udc                                 -- B_u2 : input voltage
                    ));

                rk5(t_now, buck_rk5, h);
                t_now := t_now + h;
            end loop;

            realtime <= realtime + interval_length;

            -- advance to the other half of the switching period
            if sw_state = '1' then
                sw_state := '0';
            else
                sw_state := '1';
            end if;

        end if; -- rising_edge
    end process stimulus;
------------------------------------------------------------------------
end vunit_simulation;
