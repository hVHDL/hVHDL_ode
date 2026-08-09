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

entity template_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of template_tb is

    constant clock_period : time := 1 ns;

    signal simulator_clock    : std_logic := '0';
    signal simulation_counter : natural   := 0;
    -----------------------------------
    -- simulation specific signals ----

    signal realtime : real := 0.0;
    constant stoptime : real := 100.0e-3;

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

        variable udc    : real := 10.0;
        variable i_load : real := 0.0;
        constant l      : real := 10.0e-6;
        constant c      : real := 100.0e-6;
        constant rl     : real := 1.0e-3;

        variable sw_frequency : real := 200.0e3;
        constant timestep : real := 1.0/sw_frequency;
        variable base_duty : real := 0.6;
        variable duty : real := 0.6;

        variable seed1, seed2 : positive := 1;
        variable rand : real;

        -- i_l, uc
        constant init_state_vector : real_vector := (0 => 0.0, 1 => 0.0);

        ----------
        -- average bridge voltage over a switching cycle, no switching ripple
        function get_bridge_voltage(duty : real; udc : real) return real is
        begin

            return duty * udc;

        end get_bridge_voltage;

        ----------
        impure function deriv_lcr(t : real; states : real_vector) return real_vector is
            variable retval : states'subtype := (others => 0.0);
            alias il is states(0);
            alias uc is states(1);
        begin

            if t > 250.0e-6 then i_load := 10.0; end if;
            if t > 600.0e-6 then base_duty := 0.4; end if;

            retval(0) := (get_bridge_voltage(duty, udc) - il * rl - uc) * (1.0/l);
            retval(1) := (il - i_load) * (1.0/c);

            return retval;

        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_lcr);

        variable lcr_rk5 : init_state_vector'subtype := init_state_vector;

        file file_handler : text open write_mode is "template_tb.dat";
        use ode.real_vector_pkg.all;
    begin
        if rising_edge(simulator_clock) then
            simulation_counter <= simulation_counter + 1;
            if simulation_counter = 0 then
                write_plot_config(file_handler, "title", "LCR filter testbench");
                write_plot_config(file_handler, "T_title", "Inductor current");
                write_plot_config(file_handler, "T_ylabel", "Current [A]");
                write_plot_config(file_handler, "B_title", "Bridge and capacitor voltages");
                write_plot_config(file_handler, "B_ylabel", "Voltage [V]");
                write_plot_config(file_handler, "label_T_i0", "Inductor current");
                write_plot_config(file_handler, "label_B_u0", "Bridge voltage");
                write_plot_config(file_handler, "label_B_u1", "Capacitor voltage");

                write_plot_config(file_handler, "combined_layout", "true");
                write_plot_config(file_handler, "freq_fs", real'image(sw_frequency));
                write_plot_config(file_handler, "freq_nperseg", "5000");
                -- write_plot_config(file_handler, "freq_xlim", "1000,80000");
                write_plot_config(file_handler, "freq_pair_iL", "B_u0,T_i0");
                write_plot_config(file_handler, "freq_pair_uC", "B_u0,B_u1");
                write_plot_config(file_handler, "label_iL", "Bridge voltage to inductor current");
                write_plot_config(file_handler, "label_uC", "Bridge voltage to capacitor voltage");

                init_simfile(file_handler, ("time"
                ,"T_i0"
                ,"B_u0"
                ,"B_u1"
                ));
            end if;

            -- 5% random dither on duty, same technique as fc_4level_tb,
            -- broadens the bridge voltage's spectral content so the
            -- frequency response above can be estimated from it.
            uniform(seed1, seed2, rand);
            rand := ((rand - 0.5) * 2.0) * 0.01;
            duty := base_duty + rand;

            realtime <= realtime + timestep;

            write_to(file_handler,(realtime
                    ,lcr_rk5(0)
                    ,get_bridge_voltage(duty, udc)
                    ,lcr_rk5(1)
                ));

            rk5(realtime, lcr_rk5, timestep);

        end if; -- rising_edge
    end process stimulus;
------------------------------------------------------------------------
end vunit_simulation;
