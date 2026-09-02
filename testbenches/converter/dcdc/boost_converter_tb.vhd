----------------------------------
-- Switching model of a synchronous boost converter.
--
-- Same simulation technique as fc_4level_tb / buck_converter_tb : the power
-- stage is a set of ODEs integrated with a fixed-step Dormand-Prince rk5,
-- and the switching is modelled by stepping through the converter's switch
-- states with a step length equal to that sub-interval's on- or off-time.
-- No averaging, so the inductor current ripple and switch-node waveform are
-- visible.
--
--            Vin o---L---+---[ Q2 ]---+----+---> Vout
--                        |            |    |
--                      [ Q1 ]         C   Rload
--                        |            |    |
--            GND o-------+------------+----+---> GND
--
-- Q1 on  (sw_state = '0') : inductor charges from Vin, switch node at GND,
--                           load is supplied only by C.
-- Q2 on  (sw_state = '1') : inductor drives the output, switch node at vC.
--
-- Ideal (CCM) conversion ratio : Vout = Vin / (1 - duty)
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

entity boost_converter_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of boost_converter_tb is

    constant clock_period : time := 1 ns;

    signal simulator_clock    : std_logic := '0';
    signal simulation_counter : natural   := 0;
    -----------------------------------
    -- simulation specific signals ----

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

    -- low-side (charge) duty ratio (Vout = Vin / (1 - duty)), produced by
    -- p_modulation.
    constant base_duty : real := 0.5;
    signal   duty      : real := base_duty;

    -- Handshake between the simulation step (stimulus) and the duty-cycle
    -- modulator (p_modulation): stimulus pulses sim_ready to ask for a fresh
    -- duty ratio; p_modulation computes it and pulses modulation_ready;
    -- stimulus then runs the next conduction interval.
    signal sim_ready        : boolean := false;
    signal modulation_ready : boolean := false;

    ----------------------
    -- half-bridge modulator : 1.0 while the switch node is tied to the
    -- high rail (high-side Q2 conducting, inductor feeding the output),
    -- 0.0 while it is tied to GND (low-side Q1 conducting, inductor
    -- charging). Multiplying a rail voltage or a branch current by this
    -- gives the switched quantity, e.g.
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
        constant l      : real := 47.0e-6;
        constant c      : real := 100.0e-6;
        constant rl     : real := 15.0e-3;   -- inductor series resistance
        variable rload  : real := 24.0;      -- output load resistance

        variable sw_frequency : real := 100.0e3;
        variable t_sw : real := 1.0/sw_frequency;

        -- number of fixed rk5 steps used to integrate each conduction
        -- interval. 1 => one step spanning the whole on- (or off-) time,
        -- as fc_4level_tb does. A larger value subdivides the interval into
        -- that many equal rk5 steps, giving a finer waveform and a smaller
        -- per-step truncation error at the cost of more solver calls.
        constant steps_per_on_time  : positive := 1;
        constant steps_per_off_time : positive := 1;

        -- iL, vC. vC starts at Vin, as the output cap would be pre-charged
        -- through Q2's body diode before switching begins.
        constant init_state_vector : real_vector := (0 => 0.0, 1 => 12.0);

        -- '0' -> low-side Q1 on (charge), '1' -> high-side Q2 on (transfer).
        -- duty is the low-side on-time fraction, so Vout = Vin/(1-duty).
        variable sw_state : bit := '0';
        variable simulation_started : boolean := false;

        ----------
        -- total length of the conduction interval about to be integrated
        impure function get_interval_length return real is
        begin
            if sw_state = '1' then
                return t_sw * (1.0 - duty);   -- transfer interval
            else
                return t_sw * duty;           -- charge interval
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

        impure function deriv_boost(t : real; states : real_vector) return real_vector is
            variable retval : states'subtype := (others => 0.0);
            variable v_sw   : real := 0.0;
            variable i_out  : real := 0.0;   -- inductor current delivered to the output
            alias il is states(0);
            alias vc is states(1);
        begin
            v_sw  := hb_modulator(sw_state) * vc;   -- switch-node voltage
            i_out := hb_modulator(sw_state) * il;   -- inductor current into the output

            retval(0) := (udc - il * rl - v_sw) * (1.0/l);
            retval(1) := (i_out - vc/rload)      * (1.0/c);

            return retval;
        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_boost);

        variable boost_rk5 : init_state_vector'subtype := init_state_vector;

        file file_handler : text open write_mode is "boost_converter_tb.dat";

        variable substeps       : positive := 1;
        variable interval_ticks : natural  := 0;
        variable step_ticks     : natural  := 1;
        variable now_ticks      : natural  := 0;
        variable h              : real     := 0.0;   -- rk5 step size, seconds

    begin
        if rising_edge(simulator_clock) then
            sim_ready <= false;   -- default: sim_ready is a single-cycle pulse

            if not simulation_started then
                write_plot_config(file_handler, "title", "Synchronous boost converter switching model");
                write_plot_config(file_handler, "T_title", "Inductor current");
                write_plot_config(file_handler, "T_ylabel", "Current [A]");
                write_plot_config(file_handler, "B_title", "Switch-node and output voltage");
                write_plot_config(file_handler, "B_ylabel", "Voltage [V]");
                write_plot_config(file_handler, "label_T_i0", "Inductor current");
                write_plot_config(file_handler, "label_B_u0", "Switch-node voltage");
                write_plot_config(file_handler, "label_B_u1", "Output voltage");
                write_plot_config(file_handler, "label_B_u2", "Input voltage");
                -- draw the switch-node voltage as a stepped line so it stays
                -- square even with one rk5 step per conduction interval
                -- (steps_per_on_time = steps_per_off_time = 1).
                write_plot_config(file_handler, "drawstyle_B_u0", "steps-pre");

                init_simfile(file_handler, ("time"
                ,"T_i0"
                ,"B_u0"
                ,"B_u1"
                ,"B_u2"
                ));

                sim_ready          <= true;   -- kick off the first handshake
                simulation_started := true;

        -- one conduction interval per completed duty-cycle handshake
        elsif modulation_ready then
            simulation_counter <= simulation_counter + 1;

            -- load step half way through the run
            if to_seconds(realtime_ticks) >= 2.0e-3 then
                rload := 12.0;
            end if;

            -- quantise the conduction interval to an integer number of
            -- ticks, then integrate it as N equal rk5 steps, logging a
            -- sample before each one.
            interval_ticks := to_ticks(get_interval_length);
            substeps       := get_substeps;
            step_ticks     := interval_ticks / substeps;   -- integer ticks per rk5 step
            if step_ticks < 1 then
                step_ticks := 1;
            end if;
            h              := to_seconds(step_ticks);
            now_ticks      := realtime_ticks;

            for step in 1 to substeps loop
                write_to(file_handler,(to_seconds(now_ticks)
                        , boost_rk5(0)                          -- T_i0 : inductor current
                        , hb_modulator(sw_state) * boost_rk5(1) -- B_u0 : switch-node voltage
                        , boost_rk5(1)                          -- B_u1 : output voltage
                        , udc                                   -- B_u2 : input voltage
                    ));

                rk5(to_seconds(now_ticks), boost_rk5, h);
                now_ticks := now_ticks + step_ticks;
            end loop;

            realtime_ticks <= now_ticks;

            -- advance to the other half of the switching period
            if sw_state = '1' then
                sw_state := '0';
            else
                sw_state := '1';
            end if;

            -- interval done: ask p_modulation for the next duty
            sim_ready <= true;

            end if; -- handshake
        end if; -- rising_edge
    end process stimulus;
------------------------------------------------------------------------
    -- Duty-cycle modulator : the duty command in its own process so it can
    -- be produced independently (a control loop, a schedule, ...). Runs only
    -- when stimulus pulses sim_ready, and pulses modulation_ready when the
    -- new duty is ready. Here it is just base_duty with a small dither
    -- (spectral-broadening trick from fc_4level_tb).
    p_modulation : process(simulator_clock)
        variable seed1 : positive := 7;
        variable seed2 : positive := 13;
        variable rand  : real;
    begin
        if rising_edge(simulator_clock) then
            modulation_ready <= false;   -- default: single-cycle pulse

            if sim_ready then
                uniform(seed1, seed2, rand);
                duty <= base_duty + ((rand - 0.5) * 2.0) * 0.002;
                modulation_ready <= true;
            end if;
        end if;
    end process p_modulation;
------------------------------------------------------------------------
end vunit_simulation;
