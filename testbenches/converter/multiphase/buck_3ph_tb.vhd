----------------------------------
-- Switching model of a three-phase interleaved synchronous buck converter :
-- three buck phases sharing one input source and one output capacitor /
-- load, switched 120 deg apart.
--
-- Same technique as boost_3ph_tb : the three phases switch independently,
-- so the simulation advances to the next switching edge among all phases,
-- integrates that segment with a fixed-step rk5, then toggles whichever
-- phase(s) are due.
--
--        Vin o--+--[ Q1a ]--+--L1--+----------+---> Vout
--               |           |      |          |
--               |         [ Q2a ]  |          |
--               |           |      |          |
--               +--[ Q1b ]--+--L2--+          C   Rload
--               |           |      |          |
--               |         [ Q2b ]  |          |
--               |           |      |          |
--               +--[ Q1c ]--+--L3--+----------+---> GND
--               |           |      |
--               |         [ Q2c ]  |
--        GND o--+-----------+------+
--
-- Q1x on (sw_state(x) = '1') : phase x switch node at Vin  (inductor rising)
-- Q2x on (sw_state(x) = '0') : phase x switch node at GND  (inductor falling)
--
-- Ideal (CCM) conversion ratio : Vout = duty * Vin
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

entity buck_3ph_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of buck_3ph_tb is

    constant clock_period : time := 1 ns;

    signal simulator_clock    : std_logic := '0';
    signal simulation_counter : natural   := 0;
    -----------------------------------
    -- simulation specific signals ----

    constant n_phases : natural := 3;

    -- Time is kept as an integer number of ticks. One tick lasts
    -- minimum_time_step seconds; convert to real seconds only where one is
    -- genuinely needed (the rk5 step size and the plot file).
    constant modulator_frequency : real := 125.0e6;
    constant minimum_time_step : real := 1.0/modulator_frequency;   -- 1 ns per integer time tick

    function to_seconds(ticks : integer) return real is
    begin
        return real(ticks) * minimum_time_step;
    end function;

    function to_ticks(seconds : real) return natural is
    begin
        return natural(round(seconds / minimum_time_step));
    end function;

    constant stoptime       : real    := 500.0e-3;
    constant stoptime_ticks : natural := to_ticks(stoptime);

    signal realtime_ticks : natural := 0;

    type real_phase_array is array (0 to n_phases-1) of real;

    -- per-phase high-side duty ratio (Vout = duty * Vin). Produced by the
    -- modulator process below so each phase can be commanded on its own -
    -- set base_duty per index for deliberate current sharing / phase
    -- trimming, or replace the process body with a real control loop.
    constant base_duty  : real_phase_array := (others => 0.40);
    signal   phase_duty : real_phase_array := base_duty;

    ----------------------
    -- half-bridge modulator : 1.0 while a phase's switch node is tied to
    -- the high rail (high-side Q1x conducting), 0.0 while it is tied to GND
    -- (low-side Q2x conducting). Multiplying a rail voltage by this gives
    -- the switched node voltage, e.g.
    --   v_sw = hb_modulator(sw_state) * udc
    function hb_modulator(sw_state : bit) return real is
    begin
        if sw_state = '1' then
            return 1.0;
        else
            return 0.0;
        end if;
    end hb_modulator;
    ----------------------
    -- Handshake between the simulation step (stimulus) and the duty-cycle
    -- modulator (p_modulation). stimulus pulses sim_ready to ask for a fresh
    -- set of duty ratios; p_modulation computes them and pulses
    -- modulation_ready; stimulus then runs the next integration segment.
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

        variable udc    : real := 12.0;      -- input voltage
        constant l      : real := 4.7e-6;    -- per-phase inductance
        constant c      : real := 100.0e-6;  -- shared output capacitance
        constant rl     : real := 5.0e-3;    -- per-phase inductor series resistance
        variable rload  : real := 1.2;       -- output load resistance

        variable sw_frequency : real := 100.0e3;
        variable t_sw : real := 1.0/sw_frequency;

        -- each inter-edge segment is integrated as this many equal rk5 steps
        constant steps_per_segment : positive := 1;

        -- iL1, iL2, iL3, vC. Output starts fully discharged.
        constant init_state_vector : real_vector :=
            (0 => 0.0, 1 => 0.0, 2 => 0.0, n_phases => 0.0);

        type sw_array   is array (0 to n_phases-1) of bit;
        type tick_array is array (0 to n_phases-1) of natural;

        -- '1' -> high-side on (node at Vin), '0' -> low-side on (node at GND)
        variable sw_state  : sw_array   := (others => '1');
        -- absolute tick at which each phase next toggles
        variable next_edge : tick_array := (others => 0);
        variable phases_initialised  : boolean := false;
        variable simulation_started  : boolean := false;

        impure function deriv_buck(t : real; states : real_vector) return real_vector is
            variable retval      : states'subtype := (others => 0.0);
            variable cap_current : real := 0.0;
            alias vc is states(n_phases);
        begin
            for k in 0 to n_phases-1 loop
                -- states(k) is phase k inductor current; in a buck it always
                -- flows to the output (freewheeling through the low-side
                -- device when the phase is off), so every phase feeds the
                -- output node regardless of switch state.
                retval(k)   := (hb_modulator(sw_state(k)) * udc - states(k) * rl - vc) * (1.0/l);
                cap_current := cap_current + states(k);
            end loop;
            retval(n_phases) := (cap_current - vc/rload) * (1.0/c);

            return retval;
        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_buck);

        variable buck_rk5 : init_state_vector'subtype := init_state_vector;

        file file_handler : text open write_mode is "buck_3ph_tb.dat";

        variable on_ticks      : tick_array := (others => 0);
        -- start tick of each phase's *next* switching cycle. Edges are
        -- scheduled off this fixed grid, so a phase's period is always
        -- exactly period_ticks regardless of duty dither and the phases
        -- never drift relative to each other.
        variable cycle_start   : tick_array := (others => 0);
        variable period_ticks  : natural := 0;   -- fixed switching period in ticks
        variable offset        : natural := 0;
        variable target_ticks  : natural := 0;
        variable segment_ticks : natural := 0;
        variable step_ticks    : natural := 1;
        variable now_ticks     : natural := 0;
        variable h             : real    := 0.0;   -- rk5 step size, seconds

    begin
        if rising_edge(simulator_clock) then
            sim_ready <= false;   -- default: sim_ready is a single-cycle pulse

            if not simulation_started then
                -- write the plot header once, then kick off the first
                -- handshake so p_modulation fills in phase_duty before the
                -- first integration segment runs.
                write_plot_config(file_handler, "title", "Three-phase interleaved buck converter");
                write_plot_config(file_handler, "T_title", "Phase and total output currents");
                write_plot_config(file_handler, "T_ylabel", "Current [A]");
                write_plot_config(file_handler, "B_title", "Output and input voltage");
                write_plot_config(file_handler, "B_ylabel", "Voltage [V]");
                write_plot_config(file_handler, "label_T_i0", "Phase 1 inductor current");
                write_plot_config(file_handler, "label_T_i1", "Phase 2 inductor current");
                write_plot_config(file_handler, "label_T_i2", "Phase 3 inductor current");
                write_plot_config(file_handler, "label_T_i3", "Total output current");
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

                sim_ready          <= true;
                simulation_started := true;

        -- one integration segment per completed duty-cycle handshake
        elsif modulation_ready then
            simulation_counter <= simulation_counter + 1;

            -------------------------
            -- per-phase high-side on-time from the modulator's duty command.
            -- The switching period is the constant period_ticks; only the
            -- in-cycle on/off split follows the (dithered) duty, so the
            -- period never accumulates rounding error.
            period_ticks := to_ticks(t_sw);
            for k in 0 to n_phases-1 loop
                on_ticks(k) := to_ticks(t_sw * phase_duty(k));
            end loop;

            -- load step half way through the run
            if to_seconds(realtime_ticks) >= 2.0e-3 then
                rload := 0.6;
            end if;

            -- stagger the phases 360/n_phases degrees apart on the first pass
            if not phases_initialised then
                for k in 0 to n_phases-1 loop
                    offset := (k * period_ticks) / n_phases;
                    -- cycle_start holds this phase's next cycle boundary
                    cycle_start(k) := period_ticks - offset;
                    if offset < on_ticks(k) then
                        sw_state(k)  := '1';
                        next_edge(k) := on_ticks(k) - offset;   -- '1'->'0' edge
                    else
                        sw_state(k)  := '0';
                        next_edge(k) := cycle_start(k);         -- '0'->'1' edge
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
                        , buck_rk5(0)                                     -- T_i0 : phase 1 current
                        , buck_rk5(1)                                     -- T_i1 : phase 2 current
                        , buck_rk5(2)                                     -- T_i2 : phase 3 current
                        , buck_rk5(0) + buck_rk5(1) + buck_rk5(2)         -- T_i3 : total output current
                        , buck_rk5(n_phases)                              -- B_u0 : output voltage
                        , udc                                            -- B_u1 : input voltage
                        , hb_modulator(sw_state(0)) * udc                 -- B_u2 : phase 1 switch node
                    ));

                rk5(to_seconds(now_ticks), buck_rk5, h);
                now_ticks := now_ticks + step_ticks;
            end loop;

            -- toggle every phase whose edge has now been reached
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

            -- segment done: ask p_modulation for the next duty set
            sim_ready <= true;

            end if; -- handshake
        end if; -- rising_edge
    end process stimulus;
------------------------------------------------------------------------
    -- Duty-cycle modulator : one duty ratio per phase, in its own process so
    -- each phase is driven independently. Runs only when stimulus pulses
    -- sim_ready, and pulses modulation_ready when the new duties are ready.
    -- Here it is just each phase's base duty with a small independent dither
    -- (spectral-broadening trick from fc_4level_tb); swap in per-phase
    -- references or a control loop as needed.
    p_modulation : process(simulator_clock)
        variable seed1 : positive := 7;
        variable seed2 : positive := 13;
        variable rand  : real;
    begin
        if rising_edge(simulator_clock) then
            modulation_ready <= false;   -- default: single-cycle pulse

            if sim_ready then
                for k in 0 to n_phases-1 loop
                    uniform(seed1, seed2, rand);
                    phase_duty(k) <= base_duty(k) + ((rand - 0.5) * 2.0) * 0.002;
                end loop;
                modulation_ready <= true;
            end if;
        end if;
    end process p_modulation;
------------------------------------------------------------------------
end vunit_simulation;
