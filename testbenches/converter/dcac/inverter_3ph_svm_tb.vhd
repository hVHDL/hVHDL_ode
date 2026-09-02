----------------------------------
-- Switching model of a three-phase voltage-source inverter driven by
-- space-vector modulation from a switch-state matrix, in the same style as
-- the flying-capacitor testbenches (fc_4level_tb / fc_5level_tb) : a
-- constant matrix holds the switch state for every sub-interval, the
-- modulator picks a row (here the SVM sector) and the sub-interval dwell
-- times, and the stimulus process just walks the matrix, integrating each
-- sub-interval with a fixed-step rk5.
--
-- Power stage is identical to inverter_3ph_tb : three half-bridge legs,
-- three filter inductors, three per-phase capacitors star-connected to a
-- floating neutral whose voltage un is solved algebraically each
-- derivative evaluation (get_neutral_voltage, from lcr_models_pkg).
--
-- Two switch-state matrices are provided, selected by `scheme` :
--   svpwm_7seg   - classic 7-segment symmetric SVM, V0..V7..V0,
--                  6 commutations per PWM period.
--   dpwmmin_5seg - discontinuous PWM (DPWMMIN), 5-segment V0..V0 using
--                  only the 000 zero vector : 4 commutations per period,
--                  one leg clamped to -Vdc/2 for the whole 60 deg sector,
--                  and no commutation at sector boundaries.
-- In both, T1 = m*sin(60deg-alpha), T2 = m*sin(alpha), Tz = 1 - T1 - T2,
-- with (sector, alpha) from the rotating reference vector angle.
--
-- state vector : (0..2 => iLa,iLb,iLc , 3..5 => uc_a,uc_b,uc_c)
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

entity inverter_3ph_svm_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of inverter_3ph_svm_tb is

    constant clock_period : time := 1 ns;

    signal simulator_clock    : std_logic := '0';
    signal simulation_counter : natural   := 0;
    -----------------------------------
    -- simulation specific signals ----

    constant n_phases : natural := 3;
    constant n_states : natural := 2*n_phases;   -- 3 inductor currents + 3 cap voltages

    -- Modulation / switch-state matrix selection :
    --   svpwm_7seg   : classic 7-segment symmetric SVM,  6 commutations / period
    --   dpwmmin_5seg : discontinuous PWM (DPWMMIN), only the V0 zero vector,
    --                  4 commutations / period, one leg clamped to the
    --                  negative rail for the whole 60 deg sector
    type modulation_scheme is (svpwm_7seg, dpwmmin_5seg);
    constant scheme : modulation_scheme := dpwmmin_5seg;

    function seg_count(s : modulation_scheme) return natural is
    begin
        case s is
            when svpwm_7seg   => return 7;
            when dpwmmin_5seg => return 5;
        end case;
    end function;
    constant n_segments : natural := seg_count(scheme);   -- sub-intervals per PWM period

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

    constant sw_frequency : real := 40.0e3;
    constant t_sw         : real := 1.0/sw_frequency;

    -- space-vector modulation command
    constant mod_index : real := 0.9;      -- 0..1 (linear range), m = sqrt(3)*|Vref|/Vdc
    constant f_out     : real := 500.0;    -- output fundamental frequency [Hz]

    -- modulator -> stimulus : the SVM sector for the current PWM period and
    -- the dwell time (in ticks) of each of its 7 sub-intervals.
    type seg_tick_array is array (0 to n_segments-1) of natural;
    signal svm_sector    : natural range 0 to 5 := 0;
    signal svm_seg_ticks : seg_tick_array := (others => 0);

    -- Handshake : stimulus pulses sim_ready to ask for a fresh modulator
    -- solution; p_modulation pulses modulation_ready when it is available;
    -- stimulus then integrates the next sub-interval.
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

        -- each sub-interval is integrated as this many equal rk5 steps
        constant steps_per_segment : positive := 2;

        constant init_state_vector : real_vector(0 to n_states-1) := (others => 0.0);

        type sw_array is array (0 to n_phases-1) of bit;

        -- switch-state matrices : row = SVM sector, entry = "abc" leg pattern
        -- for that sub-interval. The active process only indexes the one
        -- selected by `scheme` (see sw_pattern below).

        -- 7-segment symmetric SVM :  V0 - Va - Vb - V7 - Vb - Va - V0
        type sw_row_7 is array (0 to 6) of bit_vector(0 to n_phases-1);
        type sw_mtx_7 is array (0 to 5) of sw_row_7;
        constant sw_matrix_7seg : sw_mtx_7 := (
            0 => ("000", "100", "110", "111", "110", "100", "000"),
            1 => ("000", "110", "010", "111", "010", "110", "000"),
            2 => ("000", "010", "011", "111", "011", "010", "000"),
            3 => ("000", "011", "001", "111", "001", "011", "000"),
            4 => ("000", "001", "101", "111", "101", "001", "000"),
            5 => ("000", "101", "100", "111", "100", "101", "000"));

        -- 5-segment DPWMMIN :  V0 - Vx - Vy - Vx - V0  (only the 000 zero
        -- vector). Every transition is a single leg commutation, every row
        -- begins and ends at 000 so sector changes cost nothing, and the
        -- leg that stays '0' across the whole row is clamped to -Vdc/2 for
        -- that sector : c,c,a,a,b,b for sectors 0..5.
        type sw_row_5 is array (0 to 4) of bit_vector(0 to n_phases-1);
        type sw_mtx_5 is array (0 to 5) of sw_row_5;
        constant sw_matrix_dpwmmin : sw_mtx_5 := (
            0 => ("000", "100", "110", "100", "000"),   -- V0 V1 V2 V1 V0
            1 => ("000", "010", "110", "010", "000"),   -- V0 V3 V2 V3 V0
            2 => ("000", "010", "011", "010", "000"),   -- V0 V3 V4 V3 V0
            3 => ("000", "001", "011", "001", "000"),   -- V0 V5 V4 V5 V0
            4 => ("000", "001", "101", "001", "000"),   -- V0 V5 V6 V5 V0
            5 => ("000", "100", "101", "100", "000"));   -- V0 V1 V6 V1 V0

        impure function sw_pattern(sector : natural; seg : natural) return bit_vector is
        begin
            case scheme is
                when dpwmmin_5seg => return sw_matrix_dpwmmin(sector)(seg);
                when others       => return sw_matrix_7seg(sector)(seg);
            end case;
        end function;

        variable sw_state           : sw_array := (others => '0');
        variable simulation_started : boolean := false;

        -- current PWM period, latched at segment 0 so it holds for the
        -- whole 7-segment sequence
        variable seg_index  : natural range 0 to n_segments-1 := 0;
        variable cur_sector : natural range 0 to 5 := 0;
        variable cur_seg    : seg_tick_array := (others => 0);

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

        file file_handler : text open write_mode is "inverter_3ph_svm_tb.dat";

        variable seg_ticks  : natural := 0;
        variable step_ticks : natural := 1;
        variable now_ticks  : natural := 0;
        variable h          : real    := 0.0;   -- rk5 step size, seconds

    begin
        if rising_edge(simulator_clock) then
            sim_ready <= false;   -- default: sim_ready is a single-cycle pulse

            if not simulation_started then
                write_plot_config(file_handler, "title", "Three-phase inverter, space-vector modulation");
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
                -- keep the switch node square with one rk5 step per sub-interval
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

        -- one sub-interval per completed modulation handshake
        elsif modulation_ready then
            simulation_counter <= simulation_counter + 1;

            -- latch the sector / dwell times for the whole PWM period
            if seg_index = 0 then
                cur_sector := svm_sector;
                cur_seg    := svm_seg_ticks;
            end if;

            -- switch state for this sub-interval, straight from the matrix
            for k in 0 to n_phases-1 loop
                sw_state(k) := sw_pattern(cur_sector, seg_index)(k);
            end loop;

            -- load step half way through the run
            if to_seconds(realtime_ticks) >= 2.0e-3 then
                rload := 6.0;
            end if;

            seg_ticks := cur_seg(seg_index);
            if seg_ticks < 1 then
                seg_ticks := 1;   -- guard against a zero-dwell sub-interval
            end if;

            step_ticks := seg_ticks / steps_per_segment;
            if step_ticks < 1 then
                step_ticks := 1;
            end if;
            h := to_seconds(step_ticks);
            now_ticks := realtime_ticks;

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

            realtime_ticks <= now_ticks;

            -- advance to the next sub-interval of the sequence
            if seg_index < n_segments-1 then
                seg_index := seg_index + 1;
            else
                seg_index := 0;
            end if;

            -- sub-interval done: ask p_modulation for the next solution
            sim_ready <= true;

            end if; -- handshake
        end if; -- rising_edge
    end process stimulus;
------------------------------------------------------------------------
    -- Space-vector modulator : from the rotating reference angle it picks
    -- the sector (row of the switch-state matrix) and computes the
    -- sub-interval dwell times for the selected scheme. Replace with a
    -- closed-loop current/voltage regulator by driving |Vref| and its angle
    -- from a controller instead of a fixed magnitude and a free-running
    -- phase.
    p_modulation : process(simulator_clock)
        variable t      : real;
        variable theta  : real;
        variable alpha  : real;
        variable sector : natural range 0 to 5;
        variable t1f, t2f, tzf, scale : real;
        variable seg    : real_vector(0 to n_segments-1);
    begin
        if rising_edge(simulator_clock) then
            modulation_ready <= false;   -- default: single-cycle pulse

            if sim_ready then
                t     := to_seconds(realtime_ticks);
                theta := 2.0*MATH_PI*f_out*t;
                theta := theta - 2.0*MATH_PI*floor(theta / (2.0*MATH_PI));   -- wrap to [0, 2pi)

                sector := integer(floor(theta / (MATH_PI/3.0)));
                if sector > 5 then
                    sector := 5;
                end if;
                alpha := theta - real(sector) * (MATH_PI/3.0);              -- 0 .. pi/3

                -- dwell time fractions of t_sw for the two active vectors
                t1f := mod_index * sin(MATH_PI/3.0 - alpha);
                t2f := mod_index * sin(alpha);
                tzf := 1.0 - t1f - t2f;
                if tzf < 0.0 then
                    -- overmodulation : scale the active vectors to fill the period
                    scale := 1.0 / (t1f + t2f);
                    t1f := t1f * scale;
                    t2f := t2f * scale;
                    tzf := 0.0;
                end if;

                case scheme is
                    when svpwm_7seg =>
                        -- V0 Va Vb V7 Vb Va V0
                        seg(0) := tzf * 0.25;
                        seg(1) := t1f * 0.5;
                        seg(2) := t2f * 0.5;
                        seg(3) := tzf * 0.5;
                        seg(4) := t2f * 0.5;
                        seg(5) := t1f * 0.5;
                        seg(6) := tzf * 0.25;

                    when dpwmmin_5seg =>
                        -- V0 - Vx - Vy - Vx - V0, one zero vector only. The
                        -- V0-adjacent (odd) active vector sits on the
                        -- outside so every step is a single commutation; in
                        -- even sectors that is the first active vector, in
                        -- odd sectors the second.
                        seg(0) := tzf * 0.5;
                        seg(4) := tzf * 0.5;
                        if (sector mod 2) = 0 then
                            seg(1) := t1f * 0.5;
                            seg(2) := t2f;
                            seg(3) := t1f * 0.5;
                        else
                            seg(1) := t2f * 0.5;
                            seg(2) := t1f;
                            seg(3) := t2f * 0.5;
                        end if;
                end case;

                for i in 0 to n_segments-1 loop
                    svm_seg_ticks(i) <= to_ticks(t_sw * seg(i));
                end loop;
                svm_sector       <= sector;
                modulation_ready <= true;
            end if;
        end if;
    end process p_modulation;
------------------------------------------------------------------------
end vunit_simulation;
