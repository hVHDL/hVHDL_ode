----------------------------------
-- Switching model of a dual-active-half-bridge (DHB) DC-DC converter with
-- single-phase-shift (SPS) modulation, built from a switch-state matrix in
-- the flying-capacitor style : a constant matrix holds the two half-bridge
-- states for every sub-interval of the switching period, the modulator
-- produces the sub-interval lengths (from the phase shift) and which matrix
-- row to use (power direction), and the stimulus process just walks the
-- matrix, integrating each sub-interval with a fixed-step rk5.
--
-- Same as dab_converter_tb but each side is a single half-bridge leg
-- working against a split-capacitor divider (Cd1a/Cd1b, Cd2a/Cd2b), so the
-- transformer sees a +/-Vdc/2 square wave (half the full-bridge amplitude)
-- and each cap voltage is a state.
--
--   Vdc1 -+--Cd1a--+   Lk1,R1   m   Lk2,R2   +--Cd2a--+- + Vdc2 -+- load
--         |      M1|            |            |M2      |          |
--        [Ph]------+-a o-~~~~-+-o-+-~~~~-o b--+------[Sh]         R
--        [Pl]      |          |Lmag,Rmag      |      [Sl]         |
--   0V  -+--Cd1b---+          0V              +--Cd2b--+--- 0V ---+
--
--   v_p     = hb_p * vc1a - (1 - hb_p) * vc1b          -> +/- Vdc1/2
--   v_s_pri = n * (hb_s * vc2a - (1 - hb_s) * vc2b)    -> +/- n*Vdc2/2
--
-- A stiff Vdc1 source pins vc1a + vc1b = Vdc1, so on the primary the
-- transformer current i1 only shifts the midpoint (+/- i1/2 into the two
-- caps). On the secondary vc2a + vc2b = Vdc2 is free : the half-bridge
-- routes n*i2 to Vdc2+ (hs) or 0 (ls) and the load sits across the pair.
--
-- Transformer T model : leakage split into Lk1 (primary) and Lk2
-- (secondary, primary-referred) with a magnetizing branch Lmag,Rmag from
-- the midpoint node m to ground. The three inductor currents (i1, i2, im)
-- are state variables and the midpoint voltage vm is solved algebraically
-- each derivative evaluation from KCL at m (d/dt(i1 - i2 - im) = 0), the
-- same way the 3-phase LC models solve the floating-neutral voltage
-- (get_neutral_voltage from lcr_models_pkg).
--
-- SPS period, four sub-intervals (forward power, d > 0) :
--   [0, phi)        p=+ s=-   dur d*(Tsw/2)
--   [phi, Tsw/2)    p=+ s=+   dur (1-d)*(Tsw/2)
--   [Tsw/2, +phi)   p=- s=+   dur d*(Tsw/2)
--   [.., Tsw)       p=- s=-   dur (1-d)*(Tsw/2)
--
-- state vector : (0 => i1   primary  leakage current
--                 1 => i2   secondary leakage current (primary-referred)
--                 2 => im   magnetizing current
--                 3 => vc1a primary   split cap, Vdc1+ to M1
--                 4 => vc1b primary   split cap, M1 to 0
--                 5 => vc2a secondary split cap, Vdc2+ to M2
--                 6 => vc2b secondary split cap, M2 to 0   (Vdc2 = vc2a + vc2b))
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
    use ode.lcr_models_pkg.all;   -- get_neutral_voltage() for the T-model midpoint

entity dhb_converter_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of dhb_converter_tb is

    constant clock_period : time := 1 ns;

    signal simulator_clock    : std_logic := '0';
    signal simulation_counter : natural   := 0;
    -----------------------------------
    -- simulation specific signals ----

    constant n_legs     : natural := 2;   -- primary + secondary half-bridge high-side switch
    constant n_segments : natural := 4;   -- sub-intervals per switching period

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

    constant sw_frequency : real := 50.0e3;
    constant t_sw         : real := 1.0/sw_frequency;

    constant vdc1     : real := 400.0;   -- primary DC-link voltage (fixed source)
    constant n_turns  : real := 1.0;     -- transformer turns ratio Np:Ns
    constant vref     : real := 400.0;   -- secondary DC-link setpoint

    -- output-voltage PI regulator gains (d per volt, d per volt-second)
    constant reg_kp   : real := 1.0e-2;
    constant reg_ki   : real := 25.0;
    constant d_max    : real := 0.45;    -- phase-shift fraction clamp (|phi| < Tsw/2)

    -- modulator -> stimulus : power direction (matrix row) and the tick
    -- length of each of the four sub-intervals for this period.
    type seg_tick_array is array (0 to n_segments-1) of natural;
    signal dhb_dir       : natural range 0 to 1 := 0;   -- 0 = forward, 1 = reverse
    signal dhb_seg_ticks : seg_tick_array := (others => 0);
    signal vout          : real := vref;                -- measured Vdc2 fed back

    -- Handshake : stimulus pulses sim_ready to ask for a fresh modulator
    -- solution; p_modulation pulses modulation_ready; stimulus then
    -- integrates the next sub-interval.
    signal sim_ready        : boolean := false;
    signal modulation_ready : boolean := false;

    ----------------------
    -- half-bridge modulator : 1.0 for a conducting high-side switch, 0.0
    -- for its (complementary) low-side. Against a split-cap divider a leg
    -- output is  v = (hb_modulator(sw) - 0.5) * Vdc  (i.e. +/- Vdc/2).
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

        -- transformer T model (primary-referred)
        constant lk1   : real := 3.0e-6;     -- primary leakage inductance
        constant lk2   : real := 3.0e-6;     -- secondary leakage inductance
        constant lmag  : real := 1.0e-3;     -- magnetizing inductance (midpoint to ground)
        constant r1    : real := 30.0e-3;    -- primary winding resistance
        constant r2    : real := 30.0e-3;    -- secondary winding resistance
        -- lumped magnetizing-branch resistance (core loss + winding); also
        -- sets the flux-balancing time constant Lmag/Rmag ~ 1 ms
        constant rmag  : real := 1.0;

        -- split-capacitor dividers : two caps per half-bridge. On the
        -- primary a stiff Vdc1 source keeps vc1a + vc1b = Vdc1, so only the
        -- midpoint moves; on the secondary vc2a + vc2b = Vdc2 is free and
        -- carries the output (load across the pair).
        constant cd1   : real := 470.0e-6;   -- each primary split cap
        constant cd2   : real := 470.0e-6;   -- each secondary split cap
        variable rload : real := 40.0;       -- secondary load resistance

        -- each sub-interval is integrated as this many equal rk5 steps
        constant steps_per_segment : positive := 2;

        -- (i1, i2, im, vc1a, vc1b, vc2a, vc2b). Split caps start balanced,
        -- the secondary pair pre-charged so vc2a + vc2b = vref.
        constant init_state_vector : real_vector(0 to 6) :=
            (0 => 0.0, 1 => 0.0, 2 => 0.0,
             3 => vdc1/2.0, 4 => vdc1/2.0, 5 => vref/2.0, 6 => vref/2.0);

        type sw_array is array (0 to n_legs-1) of bit;

        -- switch-state matrix : sw_matrix(dir)(segment) is the "P S" high-side
        -- pattern for that sub-interval. Row 0 forward (primary leads), row 1
        -- reverse (secondary leads). Each half-bridge is always in a +/-Vdc/2
        -- state : plain SPS.
        type sw_seg_row    is array (0 to n_segments-1) of bit_vector(0 to n_legs-1);
        type sw_dir_matrix is array (0 to 1) of sw_seg_row;
        constant sw_matrix : sw_dir_matrix := (
            0 => ("10", "11", "01", "00"),   -- forward : primary leads
            1 => ("01", "11", "10", "00"));   -- reverse : secondary leads

        variable sw_state           : sw_array := (others => '0');
        variable simulation_started : boolean := false;

        variable seg_index : natural range 0 to n_segments-1 := 0;
        variable cur_dir   : natural range 0 to 1 := 0;
        variable cur_seg   : seg_tick_array := (others => 0);

        -- half-bridge high-side states hb_p, hb_s in {0.0, 1.0}
        impure function hb_p return real is
        begin return hb_modulator(sw_state(0)); end function;
        impure function hb_s return real is
        begin return hb_modulator(sw_state(1)); end function;

        impure function deriv_dhb(t : real; states : real_vector) return real_vector is
            variable retval  : states'subtype := (others => 0.0);
            variable v_p     : real;
            variable v_s_pri : real;   -- secondary bridge voltage, primary-referred
            variable vm      : real;   -- T-model midpoint node voltage
            variable i_load  : real;
            variable ubranch : real_vector(1 to 3);
            constant lbranch : real_vector(1 to 3) := (lk1, lk2, lmag);
            alias i1   is states(0);
            alias i2   is states(1);
            alias im   is states(2);
            alias vc1a is states(3);   -- primary   : Vdc1+ to M1
            alias vc1b is states(4);   -- primary   : M1 to 0
            alias vc2a is states(5);   -- secondary : Vdc2+ to M2
            alias vc2b is states(6);   -- secondary : M2 to 0
        begin
            -- half-bridge output = switch node minus its split-cap midpoint
            v_p     := hb_p * vc1a - (1.0 - hb_p) * vc1b;
            v_s_pri := n_turns * (hb_s * vc2a - (1.0 - hb_s) * vc2b);

            -- T-model midpoint vm from KCL at m : d/dt(i1 - i2 - im) = 0,
            --   vm = (u1/Lk1 + u2/Lk2 + u3/Lmag) / (1/Lk1 + 1/Lk2 + 1/Lmag)
            ubranch := (v_p - i1*r1, v_s_pri + i2*r2, im*rmag);
            vm := get_neutral_voltage(ubranch, lbranch);

            retval(0) := (v_p - vm - i1*r1)     / lk1;              -- d(i1)/dt   primary leakage
            retval(1) := (vm - v_s_pri - i2*r2) / lk2;              -- d(i2)/dt   secondary leakage
            retval(2) := (vm - im*rmag)         / lmag;             -- d(im)/dt   magnetizing

            -- primary split caps : i1 always enters the midpoint M1 and the
            -- stiff source holds vc1a + vc1b = Vdc1, so it splits +/- i1/2.
            retval(3) := -i1 / (2.0*cd1);                           -- d(vc1a)/dt
            retval(4) :=  i1 / (2.0*cd1);                           -- d(vc1b)/dt

            -- secondary split caps : the transformer current n*i2 is routed
            -- to Vdc2+ (hs) or 0 (ls) by the secondary half-bridge, the load
            -- sits across vc2a + vc2b.
            i_load := (vc2a + vc2b) / rload;
            retval(5) := (hb_s*n_turns*i2 - i_load)         / cd2;   -- d(vc2a)/dt
            retval(6) := (-(1.0 - hb_s)*n_turns*i2 - i_load) / cd2;  -- d(vc2b)/dt

            return retval;
        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_dhb);

        variable dhb_rk5 : init_state_vector'subtype := init_state_vector;

        file file_handler : text open write_mode is "dhb_converter_tb.dat";

        variable seg_ticks  : natural := 0;
        variable step_ticks : natural := 1;
        variable now_ticks  : natural := 0;
        variable h          : real    := 0.0;   -- rk5 step size, seconds

    begin
        if rising_edge(simulator_clock) then
            sim_ready <= false;   -- default: sim_ready is a single-cycle pulse

            if not simulation_started then
                write_plot_config(file_handler, "title", "Dual-active-half-bridge : SPS, T-model transformer, split-cap dividers");
                write_plot_config(file_handler, "T_title", "Transformer branch currents (T model)");
                write_plot_config(file_handler, "T_ylabel", "Current [A]");
                write_plot_config(file_handler, "B_title", "Half-bridge voltages and secondary DC link");
                write_plot_config(file_handler, "B_ylabel", "Voltage [V]");
                write_plot_config(file_handler, "label_T_i0", "Primary leakage current i1");
                write_plot_config(file_handler, "label_T_i1", "Secondary leakage current i2");
                write_plot_config(file_handler, "label_T_i2", "Magnetizing current im");
                write_plot_config(file_handler, "label_B_u0", "Primary half-bridge voltage (+/- Vdc1/2)");
                write_plot_config(file_handler, "label_B_u1", "Secondary half-bridge voltage (primary-referred)");
                write_plot_config(file_handler, "label_B_u2", "Secondary DC-link voltage");
                write_plot_config(file_handler, "label_B_u3", "Secondary DC-link setpoint");
                write_plot_config(file_handler, "label_B_u4", "Primary split-cap midpoint");
                write_plot_config(file_handler, "label_B_u5", "Secondary split-cap midpoint");
                -- keep the square bridge voltages square with one rk5 step per sub-interval
                write_plot_config(file_handler, "drawstyle_B_u0", "steps-pre");
                write_plot_config(file_handler, "drawstyle_B_u1", "steps-pre");

                init_simfile(file_handler, ("time"
                ,"T_i0"
                ,"T_i1"
                ,"T_i2"
                ,"B_u0"
                ,"B_u1"
                ,"B_u2"
                ,"B_u3"
                ,"B_u4"
                ,"B_u5"
                ));

                sim_ready          <= true;   -- kick off the first handshake
                simulation_started := true;

        -- one sub-interval per completed modulation handshake
        elsif modulation_ready then
            simulation_counter <= simulation_counter + 1;

            -- once per period : latch the direction / dwell times, and take
            -- a period-synchronised sample of Vdc2 for the regulator (keeps
            -- the switching ripple out of the phase-shift command)
            if seg_index = 0 then
                cur_dir := dhb_dir;
                cur_seg := dhb_seg_ticks;
                vout    <= dhb_rk5(5) + dhb_rk5(6);   -- Vdc2 = vc2a + vc2b
            end if;

            -- switch state for this sub-interval, straight from the matrix
            for k in 0 to n_legs-1 loop
                sw_state(k) := sw_matrix(cur_dir)(seg_index)(k);
            end loop;

            -- load step half way through the run
            if to_seconds(realtime_ticks) >= 2.0e-3 then
                rload := 20.0;
            end if;

            seg_ticks := cur_seg(seg_index);
            if seg_ticks < 1 then
                seg_ticks := 1;   -- guard against a zero-length sub-interval (d -> 0 or 0.5)
            end if;

            step_ticks := seg_ticks / steps_per_segment;
            if step_ticks < 1 then
                step_ticks := 1;
            end if;
            h := to_seconds(step_ticks);
            now_ticks := realtime_ticks;

            for s in 1 to steps_per_segment loop
                write_to(file_handler,(to_seconds(now_ticks)
                        , dhb_rk5(0)                                                       -- T_i0 : i1
                        , dhb_rk5(1)                                                       -- T_i1 : i2
                        , dhb_rk5(2)                                                       -- T_i2 : im
                        , hb_p*dhb_rk5(3) - (1.0-hb_p)*dhb_rk5(4)                          -- B_u0 : primary HB voltage
                        , n_turns*(hb_s*dhb_rk5(5) - (1.0-hb_s)*dhb_rk5(6))                -- B_u1 : secondary HB voltage (pri-ref)
                        , dhb_rk5(5) + dhb_rk5(6)                                          -- B_u2 : Vdc2 = vc2a + vc2b
                        , vref                                                            -- B_u3 : setpoint
                        , dhb_rk5(4)                                                       -- B_u4 : primary split-cap midpoint (vc1b)
                        , dhb_rk5(6)                                                       -- B_u5 : secondary split-cap midpoint (vc2b)
                    ));

                rk5(to_seconds(now_ticks), dhb_rk5, h);
                now_ticks := now_ticks + step_ticks;
            end loop;

            realtime_ticks <= now_ticks;

            -- advance to the next sub-interval of the SPS sequence
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
    -- Single-phase-shift modulator with an output-voltage PI loop. Produces
    -- the phase-shift fraction d (|d| < d_max), then the four sub-interval
    -- lengths  d*(Tsw/2), (1-d)*(Tsw/2), d*(Tsw/2), (1-d)*(Tsw/2)  (always
    -- summing to exactly Tsw), and the matrix row from sign(d). Replace the
    -- PI with a fixed d for an open-loop run.
    p_modulation : process(simulator_clock)
        variable t, dt, t_prev : real := 0.0;
        variable err, d, d_free, dabs, half : real := 0.0;
        variable integ : real := 0.0;
        variable seg   : real_vector(0 to n_segments-1);
    begin
        if rising_edge(simulator_clock) then
            modulation_ready <= false;   -- default: single-cycle pulse

            if sim_ready then
                t      := to_seconds(realtime_ticks);
                dt     := t - t_prev;
                t_prev := t;

                err    := vref - vout;
                d_free := reg_kp*err + integ;
                if d_free < d_max and d_free > -d_max then
                    integ := integ + reg_ki*err*dt;   -- conditional integration (anti-windup)
                end if;
                d := reg_kp*err + integ;
                if    d >  d_max then d :=  d_max;
                elsif d < -d_max then d := -d_max;
                end if;

                dabs := abs(d);
                half := t_sw * 0.5;
                seg(0) := dabs * half;
                seg(1) := (1.0 - dabs) * half;
                seg(2) := dabs * half;
                seg(3) := (1.0 - dabs) * half;

                for i in 0 to n_segments-1 loop
                    dhb_seg_ticks(i) <= to_ticks(seg(i));
                end loop;
                if d >= 0.0 then
                    dhb_dir <= 0;
                else
                    dhb_dir <= 1;
                end if;
                modulation_ready <= true;
            end if;
        end if;
    end process p_modulation;
------------------------------------------------------------------------
end vunit_simulation;
