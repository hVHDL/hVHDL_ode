----------------------------------
-- 4-level flying-capacitor converter, switching model.
--
-- Split into a stimulus process, a reference process and a modulator entity,
-- joined by the usual sim_ready / modulation_ready handshake :
--   stimulus            - integrates the LCR + flying-cap plant with a
--                          fixed-step rk5, one switch sub-interval per
--                          handshake, and feeds the inductor current back.
--   p_reference         - builds the (dithered, scheduled) modulator
--                          reference and the DC-link voltage, and hence the
--                          0..1 duty_ratio command.
--   fc_4level_modulator - takes that duty_ratio and, on each request, emits
--                          the next switch pattern and the time to apply it.
--
-- The dither on the modulator reference broadens its spectrum so the
-- reference-to-output / reference-to-current frequency responses can be
-- estimated from the single run (see test_plot.py).
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
    use ode.fc_modulator_common_pkg.all;
    use ode.fc_4level_modulator_pkg.all;

entity fc_4level_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of fc_4level_tb is

    constant clock_period : time := 1 ns;

    signal simulator_clock    : std_logic := '0';
    signal simulation_counter : natural   := 0;

    signal realtime : real := 0.0;
    constant stoptime : real := 500.0e-3;

    constant initial_dc_link : real := 200.0;

    constant sw_frequency : real := 1000.0e3;
    constant t_sw         : real := 1.0/sw_frequency;

    constant pr : real := 0.0;   -- proportional current droop on the reference

    -- stimulus -> p_reference -> modulator -> stimulus
    signal sim_ready       : boolean := false;
    signal reference_ready : boolean := false;
    signal modulator_in    : fc_modulator_input_record  := fc_modulator_input_init;
    signal modulator_out   : fc_modulator_output_record := fc_modulator_output_init;

    -- p_reference outputs
    signal udc_sig            : real := initial_dc_link;   -- DC-link voltage
    signal modulator_reference : real := 66.666/2.0;       -- command / logging
    signal duty_ratio         : real := (66.666/2.0)/initial_dc_link;
    -- stimulus -> p_reference
    signal il_meas            : real := 0.0;               -- inductor current feedback

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

        -- flat names for the modulator outputs this process consumes
        alias mod_sw_state     is modulator_out.next_switch_pattern;
        alias mod_step         is modulator_out.switching_time;
        alias modulation_ready is modulator_out.modulation_ready;

        variable i_load : real := -10.0;
        constant l      : real := 10.0e-6;
        constant c      : real := 10.0e-6;
        constant rl     : real := 50.0e-3;
        constant cfc    : real := 4.0e-6;

        -- i_l, uc, ufc1, ufc2
        constant init_state_vector : real_vector :=
            (0 => 0.0, 1 => 150.0, 2 => 66.0, 3 => 132.0);

        impure function deriv_lcr(t : real; states : real_vector) return real_vector is
            variable retval : states'subtype := (others => 0.0);
            variable bridge_voltage : real := 0.0;
            alias il   is states(0);
            alias uc   is states(1);
            alias ufc1 is states(2);
            alias ufc2 is states(3);
        begin
            bridge_voltage := get_fc_bridge_voltage(mod_sw_state, udc_sig, (ufc1, ufc2));

            retval(0) := (bridge_voltage - il * rl - uc) * (1.0/l);
            retval(1) := (il - i_load) * (1.0/c);
            retval(2) := -fc_modulator(mod_sw_state(1 downto 0)) * il / cfc;
            retval(3) := -fc_modulator(mod_sw_state(2 downto 1)) * il / cfc;

            return retval;
        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_lcr);

        variable lcr_rk5 : init_state_vector'subtype := init_state_vector;
        file file_handler : text open write_mode is "fc_4level_tb.dat";
        variable simulation_started : boolean := false;
        variable steplen : real := 0.0;

    begin
        if rising_edge(simulator_clock) then
            sim_ready <= false;

            if not simulation_started then
                write_plot_config(file_handler, "title", "4-level flying-capacitor converter");
                write_plot_config(file_handler, "T_title", "Inductor current");
                write_plot_config(file_handler, "T_ylabel", "Current [A]");
                write_plot_config(file_handler, "B_title", "Bridge, output and flying-cap voltages");
                write_plot_config(file_handler, "B_ylabel", "Voltage [V]");
                write_plot_config(file_handler, "label_T_i0", "Inductor current");
                write_plot_config(file_handler, "label_B_u0", "Output voltage");
                write_plot_config(file_handler, "label_B_u1", "Flying-cap 1 voltage");
                write_plot_config(file_handler, "label_B_u2", "Flying-cap 2 voltage");
                write_plot_config(file_handler, "label_B_u3", "DC-link voltage");
                write_plot_config(file_handler, "label_B_u4", "Modulator reference");

                -- frequency response : the dither on the modulator reference
                -- broadens its spectrum, so the reference-to-output and
                -- reference-to-current transfer functions can be estimated
                -- from a single run (see write_pkg / test_plot.py).
                write_plot_config(file_handler, "combined_layout", "true");
                write_plot_config(file_handler, "freq_unwrap_phase", "true");
                write_plot_config(file_handler, "freq_fs", real'image(sw_frequency));
                write_plot_config(file_handler, "freq_num_windows", "20");
                write_plot_config(file_handler, "freq_xlim", "1000,300000");
                write_plot_config(file_handler, "mag_ylim", "-60,20");
                write_plot_config(file_handler, "phase_ylim", "-360,90");
                write_plot_config(file_handler, "freq_title", "Modulator reference frequency response");
                write_plot_config(file_handler, "freq_pair_uC", "B_u4,B_u0");
                write_plot_config(file_handler, "freq_pair_iL", "B_u4,T_i0");
                write_plot_config(file_handler, "label_uC", "Reference -> output voltage");
                write_plot_config(file_handler, "label_iL", "Reference -> inductor current");

                init_simfile(file_handler, ("time"
                ,"T_i0"
                ,"B_u0" ,"B_u1" ,"B_u2" ,"B_u3" ,"B_u4"
                ));

                sim_ready          <= true;   -- kick off the first handshake
                simulation_started := true;

            elsif modulation_ready then
                simulation_counter <= simulation_counter + 1;
                steplen := mod_step;

                write_to(file_handler,(realtime
                        , lcr_rk5(0)          -- T_i0 : inductor current
                        , lcr_rk5(1)          -- B_u0 : output voltage
                        , lcr_rk5(2)          -- B_u1 : flying-cap 1
                        , lcr_rk5(3)          -- B_u2 : flying-cap 2
                        , udc_sig             -- B_u3 : DC link
                        , modulator_reference -- B_u4 : modulator reference
                    ));

                rk5(realtime, lcr_rk5, steplen);
                realtime <= realtime + steplen;
                il_meas  <= lcr_rk5(0);

                -- output-voltage perturbation
                if realtime < 400.0e-3 and realtime + steplen >= 400.0e-3 then
                    lcr_rk5(1) := lcr_rk5(1) - 50.0;
                end if;

                sim_ready <= true;

            end if; -- handshake
        end if; -- rising_edge
    end process stimulus;
------------------------------------------------------------------------
    p_reference : process(simulator_clock)

        constant voltage_offset : real := 66.666/2.0;

        variable udc : real := initial_dc_link;
        variable seed1, seed2 : positive := 1;
        variable rand : real;
        variable reference : real := voltage_offset;
    begin
        if rising_edge(simulator_clock) then
            reference_ready <= false;

            if sim_ready then
                -- dithered modulator reference, with scheduled steps
                uniform(seed1, seed2, rand);
                rand := ((rand - 0.5) * 2.0) * 1.0;
                reference := voltage_offset + rand;
                if realtime > 150.0e-3 then reference := 66.6666;  end if;
                if realtime > 250.0e-3 then reference := 140.9999; end if;

                -- scheduled DC-link steps
                if realtime > 100.0e-3 then udc := 150.0; end if;
                if realtime > 300.0e-3 then udc := 300.0; end if;

                reference := reference - il_meas*pr;

                modulator_reference <= reference;
                udc_sig             <= udc;
                duty_ratio          <= reference / udc;
                reference_ready     <= true;
            end if;
        end if;
    end process p_reference;
------------------------------------------------------------------------
    -- fc_4level_modulator : one 0..1 duty over the whole 0..Udc output range,
    -- no flying-cap balancing (switching_time_trim left at 0). The nominal
    -- switching period is t_sw*2 so a full 6-state sequence spans one carrier.
    modulator_in.modulation_requested <= reference_ready;
    modulator_in.duty_ratio           <= duty_ratio;
    modulator_in.switching_time_trim  <= (others => 0.0);
    modulator_in.t_sw                 <= t_sw*2.0;

    u_modulator : entity ode.fc_4level_modulator
        port map (
            clock         => simulator_clock,
            modulator_in  => modulator_in,
            modulator_out => modulator_out);
------------------------------------------------------------------------
end vunit_simulation;
