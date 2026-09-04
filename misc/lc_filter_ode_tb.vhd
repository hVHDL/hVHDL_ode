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

-- Time-domain model of the LC EMI filter in misc/emi_filter_model.qsch.
-- This is the VHDL leg of the worked model-verification example in
-- misc/README.md : the same filter checked against a bench measurement
-- three ways (measurement / QSPICE .ac / this hVHDL_ode RK4 run).
--
-- Built the same way as ode/testbenches/template_tb.vhd :
--
--   * deriv_lc() returns d/dt of the state vector for the circuit ODE
--   * a fixed-step Runge-Kutta integrator (generic_rk4, instantiated with
--     that deriv) advances the state one timestep per clock
--   * the stimulus is a DC level plus uniform dither -- the dither
--     broadens the input spectrum so a transfer function can be
--     estimated from the recorded time series (Welch H1 estimator in
--     ode/python/freq_response.py, driven by the #CONFIG freq_pair_* line
--     below), exactly the trick template_tb uses
--   * write_to() logs columns to lc_filter_ode_tb.dat ; the #CONFIG lines
--     tell ode/python/test_plot.py how to draw the time-domain traces and
--     the Bode plot
--
-- State vector (6) : the three L/C ladder sections of the schematic,
--   0: i(L1)   1: v(C1)     (L1 = 1.2 mH, R1 = 0.5, C1 = 680 nF, esr 0.1)
--   2: i(L2)   3: v(C2)     (L2 = 19 uH,  R3 = 0.2, C2 = 680 nF, esr 0.1)
--   4: i(L3)   5: v(C3)     (L3 = 19 uH,  R5 = 0.2, C3 = 680 nF, esr 0.1)
-- The probed output is the L1/C1 node V(N04) = v(C1) + esr1*(i(L1)-i(L2)),
-- the node the measured sweep in misc/L1C1_resp.csv was taken at.
entity lc_filter_ode_tb is
  generic (
    runner_cfg : string;
    -- 2 us / 0.5 s = 250k samples, split into 5 Welch windows -> ~10 Hz
    -- resolution at fs = 500 kHz, same budget as template_tb
    g_timestep : real := 2.0e-6;
    g_stoptime : real := 0.5
  );
end;

architecture vunit_simulation of lc_filter_ode_tb is

    constant clock_period : time := 1 ns;

    signal simulator_clock    : std_logic := '0';
    signal simulation_counter : natural   := 0;

    signal realtime : real := 0.0;

begin

------------------------------------------------------------------------
    simtime : process
    begin
        test_runner_setup(runner, runner_cfg);
        wait until realtime >= g_stoptime;
        test_runner_cleanup(runner); -- Simulation ends here
        wait;
    end process simtime;

    simulator_clock <= not simulator_clock after clock_period/2.0;
------------------------------------------------------------------------

    stimulus : process(simulator_clock)

        -- component values, misc/emi_filter_model.qsch
        constant l1 : real := 1.2e-3;   constant rl1 : real := 0.5;
        constant l2 : real := 19.0e-6;  constant rl2 : real := 0.2;
        constant l3 : real := 19.0e-6;  constant rl3 : real := 0.2;
        constant c1 : real := 680.0e-9; constant esr1 : real := 0.1;
        constant c2 : real := 680.0e-9; constant esr2 : real := 0.1;
        constant c3 : real := 680.0e-9; constant esr3 : real := 0.1;

        constant timestep : real := g_timestep;

        variable seed1, seed2 : positive := 1;
        variable rand : real;
        variable input_voltage : real := 1.0;

        constant init_state : real_vector := (0 to 5 => 0.0);

        ----------
        impure function deriv_lc(t : real; states : real_vector) return real_vector is
            variable d : states'subtype := (others => 0.0);
            alias il1 is states(0); alias uc1 is states(1);
            alias il2 is states(2); alias uc2 is states(3);
            alias il3 is states(4); alias uc3 is states(5);
            variable v1, v2, v3 : real;
        begin
            -- node voltages : capacitor voltage plus its esr drop
            v1 := uc1 + esr1 * (il1 - il2);
            v2 := uc2 + esr2 * (il2 - il3);
            v3 := uc3 + esr3 *  il3;

            d(0) := (input_voltage - il1*rl1 - v1) / l1; -- d i(L1)
            d(1) := (il1 - il2) / c1;                    -- d v(C1)
            d(2) := (v1 - il2*rl2 - v2) / l2;            -- d i(L2)
            d(3) := (il2 - il3) / c2;                    -- d v(C2)
            d(4) := (v2 - il3*rl3 - v3) / l3;            -- d i(L3)
            d(5) := il3 / c3;                            -- d v(C3)
            return d;
        end function;

        procedure rk4 is new generic_rk4 generic map(deriv_lc);

        variable state : init_state'subtype := init_state;

        impure function n04 return real is  -- V(N04), the L1/C1 node
        begin
            return state(1) + esr1 * (state(0) - state(2));
        end function;

        file file_handler : text open write_mode is "lc_filter_ode_tb.dat";
        use ode.real_vector_pkg.all;
    begin
        if rising_edge(simulator_clock) then
            if simulation_counter = 0 then
                write_plot_config(file_handler, "title", "LC EMI filter (misc/emi_filter_model.qsch) - hVHDL_ode RK4");
                write_plot_config(file_handler, "T_title", "L1 inductor current");
                write_plot_config(file_handler, "T_ylabel", "Current [A]");
                write_plot_config(file_handler, "B_title", "Input and L1/C1 node voltage");
                write_plot_config(file_handler, "B_ylabel", "Voltage [V]");
                write_plot_config(file_handler, "label_T_i0", "i(L1)");
                write_plot_config(file_handler, "label_B_u0", "input V1");
                write_plot_config(file_handler, "label_B_u1", "V(N04)  [L1/C1 node]");
                write_plot_config(file_handler, "xlim", "0,5e-3");

                -- Bode plot : H(f) = V(N04) / Vin, estimated from the dithered run
                write_plot_config(file_handler, "combined_layout", "true");
                write_plot_config(file_handler, "freq_unwrap_phase", "true");
                -- the stimulus is sampled and held once per timestep, so the
                -- Welch estimate carries a half-sample (T/2) zero-order-hold
                -- group delay the real circuit has not; deskew it so the phase
                -- lines up with the QSPICE .ac / measured sweep
                write_plot_config(file_handler, "freq_deskew", "zoh");
                write_plot_config(file_handler, "freq_fs", real'image(1.0/timestep));
                write_plot_config(file_handler, "freq_num_windows", "5");
                write_plot_config(file_handler, "freq_xlim", "5e2,1e5");
                write_plot_config(file_handler, "mag_ylim", "-60,40");
                write_plot_config(file_handler, "phase_ylim", "-220,20");
                write_plot_config(file_handler, "freq_title",
                    "LC EMI filter V(N04)/Vin : hVHDL_ode vs QSPICE vs measured");
                write_plot_config(file_handler, "freq_pair_L1C1", "B_u0,B_u1");
                write_plot_config(file_handler, "label_L1C1", "V(N04)/Vin  (hVHDL_ode)");
                -- dump the RK4 response so it can be overlaid without re-running
                -- the sim ; misc/README.md walks through the full three-way
                -- (measured / QSPICE / hVHDL_ode) comparison
                write_plot_config(file_handler, "freq_save_L1C1", "misc/lc_vhdl_resp.csv");

                init_simfile(file_handler, ("time", "T_i0", "B_u0", "B_u1"));
            end if;
            simulation_counter <= simulation_counter + 1;

            -- DC + 10% uniform dither, same broadband-excitation trick as
            -- template_tb ; the small DC step part-way through adds low
            -- frequency content
            uniform(seed1, seed2, rand);
            input_voltage := 1.0 + ((rand - 0.5) * 2.0) * 0.1;
            if realtime > g_stoptime/2.0 then
                input_voltage := input_voltage + 1.0;
            end if;

            realtime <= realtime + timestep;

            write_to(file_handler, (realtime, state(0), input_voltage, n04));

            rk4(realtime, state, timestep);

        end if; -- rising_edge
    end process stimulus;
------------------------------------------------------------------------
end vunit_simulation;
