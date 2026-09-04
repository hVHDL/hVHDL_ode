----------------------------------
-- 5-level flying-capacitor converter, switching model.
--
-- Split into a stimulus process and a modulator entity, joined by the
-- usual sim_ready / modulation_ready handshake :
--   stimulus             - integrates the LCR + flying-cap plant with a
--                           fixed-step rk5, one switch sub-interval per
--                           handshake, and feeds the inductor current back.
--   fc_5level_modulator  - takes a single 0..1 duty_ratio over the whole
--                           0..Udc range and, on each request, emits the
--                           next switch pattern and the time to apply it.
--
-- Open-loop : no flying-cap balancing (switching_time_trim left at 0). The
-- modulator still carries that hook, and fc_modulator_common_pkg.get_fc_trims
-- can fill it from the measured cap errors if a balancing run is wanted.
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
    use ode.fc_5level_modulator_pkg.all;

entity fc_5level_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of fc_5level_tb is

    constant clock_period : time := 1 ns;

    signal simulator_clock    : std_logic := '0';
    signal simulation_counter : natural   := 0;

    signal realtime : real := 0.0;
    constant stoptime : real := 400.0e-3;

    constant initial_dc_link     : real := 200.0;
    constant initial_voltage_ref : real := 149.0;

    constant sw_frequency : real := 500.0e3;
    constant t_sw         : real := 1.0/sw_frequency;

    constant pr : real := 0.0;   -- proportional current droop on the reference

    -- stimulus <-> modulator : request strobe and the wrapped port records
    signal sim_ready     : boolean := false;
    signal modulator_in  : fc_modulator_input_record  := fc_modulator_input_init;
    signal modulator_out : fc_modulator_output_record := fc_modulator_output_init;

    signal udc_sig      : real         := initial_dc_link;   -- DC-link (constant here)
    -- stimulus -> modulator
    signal il_meas      : real         := 0.0;
    -- modulator command : one 0..1 duty over the whole 0..Udc output range
    signal modulator_reference : real := initial_voltage_ref;   -- logging / command
    signal duty_ratio          : real := initial_voltage_ref/initial_dc_link;

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

        variable i_load : real := 10.1111;
        constant l      : real := 20.0e-6;
        constant c      : real := 10.0e-6;
        constant rl     : real := 100.0e-3;
        constant cfc    : real := 2.0e-6;

        -- il, uc, ufc1, ufc2, ufc3
        constant init_state_vector : real_vector := (
              0 => 0.0
            , 1 => initial_voltage_ref
            , 2 => initial_dc_link*1.0/4.0
            , 3 => initial_dc_link*2.0/4.0
            , 4 => initial_dc_link*3.0/4.0);

        impure function deriv_lcr(t : real; states : real_vector) return real_vector is
            variable retval : states'subtype := (others => 0.0);
            variable bridge_voltage : real := 0.0;
            alias il   is states(0);
            alias uc   is states(1);
            alias ufc1 is states(2);
            alias ufc2 is states(3);
            alias ufc3 is states(4);
        begin
            bridge_voltage := get_fc_bridge_voltage(mod_sw_state, udc_sig, (ufc1, ufc2, ufc3));

            retval(0) := (bridge_voltage - il * rl - uc) * (1.0/l);
            retval(1) := (il - i_load) * (1.0/c);
            retval(2) := -fc_modulator(mod_sw_state(1 downto 0)) * il / cfc;
            retval(3) := -fc_modulator(mod_sw_state(2 downto 1)) * il / cfc;
            retval(4) := -fc_modulator(mod_sw_state(3 downto 2)) * il / cfc;

            return retval;
        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_lcr);

        variable lcr_rk5 : init_state_vector'subtype := init_state_vector;
        file file_handler : text open write_mode is "fc_5level_tb.dat";
        variable simulation_started : boolean := false;
        variable steplen : real := 0.0;

    begin
        if rising_edge(simulator_clock) then
            sim_ready <= false;

            if not simulation_started then
                write_plot_config(file_handler, "title", "5-level flying-capacitor converter");
                write_plot_config(file_handler, "T_title", "Flying-cap, DC-link, reference and output voltages");
                write_plot_config(file_handler, "T_ylabel", "Voltage [V]");
                write_plot_config(file_handler, "B_title", "Inductor current");
                write_plot_config(file_handler, "B_ylabel", "Current [A]");
                write_plot_config(file_handler, "label_T_i0", "Flying-cap 1 voltage");
                write_plot_config(file_handler, "label_T_i1", "Flying-cap 2 voltage");
                write_plot_config(file_handler, "label_T_i2", "Flying-cap 3 voltage");
                write_plot_config(file_handler, "label_T_i3", "DC-link voltage");
                write_plot_config(file_handler, "label_T_i4", "Modulator reference");
                write_plot_config(file_handler, "label_T_i5", "Output voltage");
                write_plot_config(file_handler, "label_B_u0", "Inductor current");

                init_simfile(file_handler, ("time"
                ,"T_i0" ,"T_i1" ,"T_i2" ,"T_i3" ,"T_i4" ,"T_i5"
                ,"B_u0"
                ));

                sim_ready          <= true;   -- kick off the first handshake
                simulation_started := true;

            elsif modulation_ready then
                simulation_counter <= simulation_counter + 1;
                steplen := mod_step;

                write_to(file_handler,(realtime
                        , lcr_rk5(2)   -- T_i0 : flying-cap 1
                        , lcr_rk5(3)   -- T_i1 : flying-cap 2
                        , lcr_rk5(4)   -- T_i2 : flying-cap 3
                        , udc_sig             -- T_i3 : DC link
                        , modulator_reference -- T_i4 : modulator reference
                        , lcr_rk5(1)   -- T_i5 : output voltage
                        , lcr_rk5(0)   -- B_u0 : inductor current
                    ));

                rk5(realtime, lcr_rk5, steplen);
                realtime <= realtime + steplen;
                il_meas  <= lcr_rk5(0);

                sim_ready <= true;

            end if; -- handshake
        end if; -- rising_edge
    end process stimulus;
------------------------------------------------------------------------
    -- modulator command : a single 0..1 duty over the whole 0..Udc output
    -- range, no flying-cap balancing (switching_time_trim left at 0)
    modulator_reference <= initial_voltage_ref - il_meas*pr;
    duty_ratio          <= modulator_reference / udc_sig;

    modulator_in.modulation_requested <= sim_ready;
    modulator_in.duty_ratio           <= duty_ratio;
    modulator_in.switching_time_trim  <= (others => 0.0);
    modulator_in.t_sw                 <= t_sw;

    u_modulator : entity ode.fc_5level_modulator
        port map (
            clock         => simulator_clock,
            modulator_in  => modulator_in,
            modulator_out => modulator_out);
------------------------------------------------------------------------
end vunit_simulation;
