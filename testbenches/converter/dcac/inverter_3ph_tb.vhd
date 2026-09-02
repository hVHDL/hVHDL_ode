----------------------------------
-- Switching model of a three-phase voltage-source inverter with an LC
-- output filter.
--
-- Same power stage as buck_3ph_tb (three half-bridge legs, each with a
-- filter inductor), and the same 3-phase LC filter as lcr_3ph_tb : three
-- separate capacitors, one per phase, star-connected to a floating
-- neutral. The neutral voltage un is not a state; it is solved
-- algebraically each derivative evaluation from the constraint that the
-- three inductor currents sum to zero (3-wire system), using
-- get_neutral_voltage() from lcr_models_pkg.
--
--   +Vdc/2 o--+--[Q1a]--+--La,Ra--o a --Ca--+           Rload (wye) from
--             |         |                   |           each phase node to
--             |       [Q2a]                 |           0 V; the caps and
--             |         |                   |           load return through
--   -Vdc/2 o--+--[Q1b]--+--Lb,Rb--o b --Cb--+-- n       the floating neutral
--             |         |                   |   (un)
--             |       [Q2b]                 |
--             |         |                   |
--             +--[Q1c]--+--Lc,Rc--o c --Cc--+
--             |         |
--             |       [Q2c]
--   -Vdc/2 o--+---------+
--
-- Q1x on (sw_state(x) = '1') : leg x switch node at +Vdc/2
-- Q2x on (sw_state(x) = '0') : leg x switch node at -Vdc/2
--
-- state vector : (0..2 => iLa,iLb,iLc , 3..5 => uc_a,uc_b,uc_c)
-- The cap voltage states uc(k) are already referenced to the floating
-- neutral, so uc(k) is the line-to-neutral phase voltage directly; un only
-- enters the inductor equations.
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
    use ode.lcr_models_pkg.all;   -- get_neutral_voltage()

entity inverter_3ph_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of inverter_3ph_tb is

    constant clock_period : time := 1 ns;

    signal simulator_clock    : std_logic := '0';
    signal simulation_counter : natural   := 0;
    -----------------------------------
    -- simulation specific signals ----

    constant n_phases : natural := 3;
    constant n_states : natural := 2*n_phases;   -- 3 inductor currents + 3 cap voltages

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

    constant stoptime       : real    := 6.0e-3;
    constant stoptime_ticks : natural := to_ticks(stoptime);

    signal realtime_ticks : natural := 0;

    -- modulator command : sinusoidal per-leg duty, Vphase = (duty-0.5)*Vdc
    constant mod_index : real := 0.9;       -- 0..1
    constant f_out     : real := 500.0;     -- output fundamental frequency [Hz]

    type real_phase_array is array (0 to n_phases-1) of real;
    signal phase_duty : real_phase_array := (others => 0.5);

    -- Handshake between the simulation step (stimulus) and the modulator
    -- (p_modulation): stimulus pulses sim_ready to ask for a fresh set of
    -- duty references; p_modulation computes them and pulses
    -- modulation_ready; stimulus then runs the next integration segment.
    signal sim_ready        : boolean := false;
    signal modulation_ready : boolean := false;

    ----------------------
    -- half-bridge modulator : 1.0 when a leg's high-side device conducts
    -- (switch node at +Vdc/2), 0.0 when the low-side conducts (-Vdc/2).
    --   v_sw = (hb_modulator(sw_state) - 0.5) * vdc
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

        variable vdc    : real := 400.0;     -- DC-link voltage (legs swing +/- vdc/2)
        constant l      : real_vector(1 to 3) := (others => 300.0e-6);  -- per-phase filter inductance
        constant c      : real_vector(1 to 3) := (others => 10.0e-6);   -- per-phase filter capacitance
        constant r      : real_vector(1 to 3) := (others => 100.0e-3);  -- per-phase inductor series resistance
        variable rload  : real := 12.0;      -- wye load resistance, phase node to 0 V

        variable sw_frequency : real := 40.0e3;
        variable t_sw : real := 1.0/sw_frequency;

        -- each inter-edge segment is integrated as this many equal rk5 steps
        constant steps_per_segment : positive := 2;

        constant init_state_vector : real_vector(0 to n_states-1) := (others => 0.0);

        type sw_array   is array (0 to n_phases-1) of bit;
        type tick_array is array (0 to n_phases-1) of natural;

        -- '1' -> high-side on (+Vdc/2), '0' -> low-side on (-Vdc/2)
        variable sw_state  : sw_array   := (others => '1');
        variable next_edge : tick_array := (others => 0);
        variable legs_initialised    : boolean := false;
        variable simulation_started  : boolean := false;

        ----------
        -- inductor voltages before the neutral shift : ul(k) = vsw(k) - uc(k) - iL(k)*R(k)
        impure function get_inductor_voltages(states : real_vector) return real_vector is
            variable ul : real_vector(1 to 3);
            alias il : real_vector(1 to 3) is states(0 to 2);
            alias uc : real_vector(1 to 3) is states(3 to 5);
        begin
            for k in 1 to 3 loop
                ul(k) := (hb_modulator(sw_state(k-1)) - 0.5)*vdc - uc(k) - il(k)*r(k);
            end loop;
            return ul;
        end function;
        ----------
        impure function deriv_inv(t : real; states : real_vector) return real_vector is
            variable retval : states'subtype := (others => 0.0);
            variable ul     : real_vector(1 to 3);
            variable un     : real;
            alias il : real_vector(1 to 3) is states(0 to 2);
            alias uc : real_vector(1 to 3) is states(3 to 5);
        begin
            ul := get_inductor_voltages(states);
            un := get_neutral_voltage(ul, l);

            for k in 1 to 3 loop
                retval(k-1) := (ul(k) - un) / l(k);              -- d(iL)/dt
                retval(k+2) := (il(k) - uc(k)/rload) / c(k);     -- d(uc)/dt  (states 3..5)
            end loop;

            return retval;
        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_inv);

        variable inv_rk5 : init_state_vector'subtype := init_state_vector;

        file file_handler : text open write_mode is "inverter_3ph_tb.dat";

        variable on_ticks      : tick_array := (others => 0);
        variable cycle_start   : tick_array := (others => 0);
        variable period_ticks  : natural := 0;   -- fixed switching period in ticks
        variable target_ticks  : natural := 0;
        variable segment_ticks : natural := 0;
        variable step_ticks    : natural := 1;
        variable now_ticks     : natural := 0;
        variable h             : real    := 0.0;   -- rk5 step size, seconds

    begin
        if rising_edge(simulator_clock) then
            sim_ready <= false;   -- default: sim_ready is a single-cycle pulse

            if not simulation_started then
                write_plot_config(file_handler, "title", "Three-phase voltage-source inverter");
                write_plot_config(file_handler, "T_title", "Phase inductor currents");
                write_plot_config(file_handler, "T_ylabel", "Current [A]");
                write_plot_config(file_handler, "B_title", "Leg-a switch node and phase output voltages");
                write_plot_config(file_handler, "B_ylabel", "Voltage [V]");
                write_plot_config(file_handler, "label_T_i0", "Phase a inductor current");
                write_plot_config(file_handler, "label_T_i1", "Phase b inductor current");
                write_plot_config(file_handler, "label_T_i2", "Phase c inductor current");
                write_plot_config(file_handler, "label_B_u0", "Leg a switch-node voltage");
                write_plot_config(file_handler, "label_B_u1", "Phase a voltage (line-to-neutral)");
                write_plot_config(file_handler, "label_B_u2", "Phase b voltage (line-to-neutral)");
                write_plot_config(file_handler, "label_B_u3", "Phase c voltage (line-to-neutral)");
                -- keep the switch node square with one rk5 step per segment
                write_plot_config(file_handler, "drawstyle_B_u0", "steps-pre");

                init_simfile(file_handler, ("time"
                ,"T_i0"
                ,"T_i1"
                ,"T_i2"
                ,"B_u0"
                ,"B_u1"
                ,"B_u2"
                ,"B_u3"
                ));

                sim_ready          <= true;   -- kick off the first handshake
                simulation_started := true;

        -- one integration segment per completed modulation handshake
        elsif modulation_ready then
            simulation_counter <= simulation_counter + 1;

            -------------------------
            -- All three legs share one PWM carrier (no interleaving). The
            -- switching period is the constant period_ticks; only each leg's
            -- in-cycle on-time follows its (sinusoidal) duty, so the period
            -- never accumulates rounding error.
            period_ticks := to_ticks(t_sw);
            for k in 0 to n_phases-1 loop
                on_ticks(k) := to_ticks(t_sw * phase_duty(k));
            end loop;

            -- load step half way through the run
            if to_seconds(realtime_ticks) >= 2.0e-3 then
                rload := 6.0;
            end if;

            if not legs_initialised then
                for k in 0 to n_phases-1 loop
                    sw_state(k)    := '1';
                    next_edge(k)   := on_ticks(k);      -- '1'->'0' edge
                    cycle_start(k) := period_ticks;     -- next cycle boundary
                end loop;
                legs_initialised := true;
            end if;

            now_ticks := realtime_ticks;

            -- advance only as far as the next switching edge among all legs
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
                        , inv_rk5(0)                               -- T_i0 : phase a current
                        , inv_rk5(1)                               -- T_i1 : phase b current
                        , inv_rk5(2)                               -- T_i2 : phase c current
                        , (hb_modulator(sw_state(0)) - 0.5)*vdc    -- B_u0 : leg a switch node
                        , inv_rk5(3)                               -- B_u1 : phase a line-to-neutral
                        , inv_rk5(4)                               -- B_u2 : phase b line-to-neutral
                        , inv_rk5(5)                               -- B_u3 : phase c line-to-neutral
                    ));

                rk5(to_seconds(now_ticks), inv_rk5, h);
                now_ticks := now_ticks + step_ticks;
            end loop;

            -- toggle every leg whose edge has now been reached
            for k in 0 to n_phases-1 loop
                if now_ticks >= next_edge(k) then
                    if sw_state(k) = '1' then
                        -- high-side off; low-side holds until the fixed cycle boundary
                        sw_state(k)  := '0';
                        next_edge(k) := cycle_start(k);
                    else
                        -- new cycle: advance the grid by exactly one period
                        sw_state(k)    := '1';
                        next_edge(k)   := cycle_start(k) + on_ticks(k);
                        cycle_start(k) := cycle_start(k) + period_ticks;
                    end if;
                end if;
            end loop;

            realtime_ticks <= now_ticks;

            -- segment done: ask p_modulation for the next duty references
            sim_ready <= true;

            end if; -- handshake
        end if; -- rising_edge
    end process stimulus;
------------------------------------------------------------------------
    -- Modulator : sinusoidal, 120 deg-shifted duty reference per leg, in its
    -- own process so it can be replaced by a closed-loop current/voltage
    -- controller. Runs only when stimulus pulses sim_ready and pulses
    -- modulation_ready when the new references are ready.
    p_modulation : process(simulator_clock)
        variable t : real;
    begin
        if rising_edge(simulator_clock) then
            modulation_ready <= false;   -- default: single-cycle pulse

            if sim_ready then
                t := to_seconds(realtime_ticks);
                for k in 0 to n_phases-1 loop
                    phase_duty(k) <= 0.5 + 0.5 * mod_index
                        * sin(2.0*MATH_PI*f_out*t + real(k) * (2.0*MATH_PI/3.0));
                end loop;
                modulation_ready <= true;
            end if;
        end if;
    end process p_modulation;
------------------------------------------------------------------------
end vunit_simulation;
