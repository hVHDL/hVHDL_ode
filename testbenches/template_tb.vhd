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
    constant stoptime : real := 500.0e-3;

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
        constant l      : real := 2.0e-6;
        constant c      : real := 100.0e-6;
        constant rl     : real := 2.0e-3;

        variable sw_frequency : real := 500.0e3;
        constant timestep : real := 1.0/sw_frequency;
        variable base_duty : real := 0.6;
        variable duty : real := 0.6;

        variable seed1, seed2 : positive := 1;
        variable rand : real;
        variable input_voltage : real := 0.0;

        -- i_l, uc
        constant init_state_vector : real_vector := (0 => 0.0, 1 => 0.0, 2 => 0.0, 3 => 0.0);

        ----------
        -- average bridge voltage over a switching cycle, no switching ripple
        function get_bridge_voltage(duty : real; udc : real) return real is
        begin

            return duty * udc;

        end get_bridge_voltage;

        ----------
        impure function deriv_lcr(t : real; states : real_vector) return real_vector is
            variable retval : states'subtype := (others => 0.0);
            alias il  is states(0);
            alias uc  is states(1);
            alias i2  is states(2);
            alias uc2 is states(3);
        begin

            -- if t > 250.0e-6 then i_load := 10.0; end if;
            if t > 600.0e-6 then base_duty := 0.4; end if;

            retval(0) := (input_voltage - il * rl - uc) * (1.0/l);
            retval(1) := (il - i2) * (1.0/c);
            retval(2) := (uc - uc2) * (1.0/c);
            retval(3) := (i2- i_load) * (1.0/c);

            return retval;

        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_lcr);

        variable lcr_rk5 : init_state_vector'subtype := init_state_vector;

        file file_handler : text open write_mode is "template_tb.dat";
        use ode.real_vector_pkg.all;
    begin
        if rising_edge(simulator_clock) then
            if simulation_counter = 0 then
                write_plot_config(file_handler, "title", "LCR filter testbench");
                write_plot_config(file_handler, "T_title", "Inductor current");
                write_plot_config(file_handler, "T_ylabel", "Current [A]");
                write_plot_config(file_handler, "B_title", "Bridge and capacitor voltages");
                write_plot_config(file_handler, "B_ylabel", "Voltage [V]");
                write_plot_config(file_handler, "label_T_i0", "Inductor current");
                write_plot_config(file_handler, "label_B_u0", "Bridge voltage");
                write_plot_config(file_handler, "label_B_u1", "Capacitor voltage");
                write_plot_config(file_handler, "xlim", "0,5e-3");
                -- ylim to fit the large-signal excitation below (input_voltage
                -- steps 1 V -> 2 V); leave commented to autoscale
                write_plot_config(file_handler, "T_ylim", "-8,8");
                write_plot_config(file_handler, "B_ylim", "-0.5,4.5");

                -- bode plot
                write_plot_config(file_handler, "combined_layout", "true");
                -- this LCR stage is lightly damped (rl only), so the two
                -- resonance peaks reach ~+40..+55 dB and above the cutoff the
                -- capacitor-voltage estimate is mostly noise; keep the phase
                -- wrapped (unwrap would run that noise off the axis) and give
                -- the magnitude/phase enough range to show the peaks.
                -- write_plot_config(file_handler, "freq_unwrap_phase", "true");
                write_plot_config(file_handler, "freq_fs", real'image(sw_frequency));
                write_plot_config(file_handler, "freq_num_windows", "5");
                write_plot_config(file_handler, "freq_xlim", "2e2,100e3");
                write_plot_config(file_handler, "mag_ylim", "-60,60");
                write_plot_config(file_handler, "phase_ylim", "-200,200");
                write_plot_config(file_handler, "freq_pair_iL", "B_u0,T_i0");
                write_plot_config(file_handler, "freq_pair_uC", "B_u0,B_u1");
                write_plot_config(file_handler, "label_iL", "Bridge voltage to inductor current");
                write_plot_config(file_handler, "label_uC", "Bridge voltage to capacitor voltage");
                -- save the computed response so it can be overlaid on a
                -- later run's plot, e.g.: python test_plot.py other.dat template_tb_uC.csv
                -- write_plot_config(file_handler, "freq_save_uC", "template_tb_uC.csv");

                init_simfile(file_handler, ("time"
                ,"T_i0"
                ,"B_u0"
                ,"B_u1"
                ));
            end if;
            simulation_counter <= simulation_counter + 1;

            -- 5% random dither on duty, same technique as fc_4level_tb,
            -- broadens the bridge voltage's spectral content so the
            -- frequency response above can be estimated from it.
            uniform(seed1, seed2, rand);
            rand := ((rand - 0.5) * 2.0) * 0.1;

            input_voltage := rand + 1.0;

            if realtime > 2.5e-3 then
                input_voltage := rand + 2.0;
            end if;

            realtime <= realtime + timestep;

            write_to(file_handler,(realtime
                    ,lcr_rk5(0)
                    ,input_voltage
                    ,lcr_rk5(3)
                ));

            rk5(realtime, lcr_rk5, timestep);

        end if; -- rising_edge
    end process stimulus;
------------------------------------------------------------------------
end vunit_simulation;
