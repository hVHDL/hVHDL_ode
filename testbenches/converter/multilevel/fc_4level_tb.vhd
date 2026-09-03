----------------------------------
-- 4-level flying-capacitor converter, switching model.
--
-- Split into two processes like the other converter testbenches :
--   stimulus     - integrates the LCR + flying-cap plant with a fixed-step
--                  rk5, one switch sub-interval per handshake, and feeds
--                  the inductor current back to the modulator.
--   p_modulation - from the (dithered) modulator reference picks the row of
--                  fc_4_sw_matrix and the sub-interval length, runs the
--                  flying-cap balancing integrators, and steps Vdc / the
--                  reference on a schedule.
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

entity fc_4level_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of fc_4level_tb is

    constant clock_period : time := 1 ns;

    signal simulator_clock    : std_logic := '0';
    signal simulation_counter : natural   := 0;

    signal realtime : real := 0.0;
    constant stoptime : real := 500.0e-3;

    constant sw_frequency : real := 1000.0e3;
    constant t_sw         : real := 1.0/sw_frequency;

    subtype sw_states is bit_vector(2 downto 0);

    -- stimulus <-> p_modulation handshake
    signal sim_ready        : boolean := false;
    signal modulation_ready : boolean := false;
    -- p_modulation -> stimulus
    signal mod_sw_state : sw_states := "111";
    signal mod_step     : real      := t_sw*0.5;   -- sub-interval length [s]
    signal udc_sig      : real      := 200.0;      -- DC-link voltage
    signal mod_ref      : real      := 0.0;        -- modulator reference (for logging)
    signal fc1_balance  : real      := 0.0;        -- flying-cap balancing integrators
    signal fc2_balance  : real      := 0.0;
    -- stimulus -> p_modulation
    signal il_meas      : real      := 0.0;        -- inductor current feedback

    ----------------------
    function fc_modulator(gate_signals : bit_vector) return real is
        variable retval : real;
    begin
        CASE gate_signals is
            WHEN "10" => retval := -1.0;
            WHEN "01" => retval :=  1.0;
            WHEN others => retval := 0.0;
        end CASE;
        return retval;
    end fc_modulator;
    ----------------------
    function number_of_ones(vector : bit_vector) return natural is
        variable retval : natural := 0;
    begin
        for i in vector'range loop
            if vector(i) = '1' then retval := retval + 1; end if;
        end loop;
        return retval;
    end number_of_ones;
    ----------------------
    function get_fc_bridge_voltage(sw_state : sw_states ; udc : real; ufc : real_vector) return real is
        variable bridge_voltage : real := 0.0;
    begin
        for i in ufc'range loop
            bridge_voltage := bridge_voltage + fc_modulator(sw_state(i+1 downto i)) * ufc(i);
        end loop;
        bridge_voltage := bridge_voltage + fc_modulator('0' & sw_state(sw_state'high)) * udc;
        return bridge_voltage;
    end get_fc_bridge_voltage;
    ----------------------
    function get_fc_duty(vref : real; udc : real ; level_bits : bit_vector) return real is
        variable retval : real := 0.0;
        variable imax : natural := level_bits'high;
        constant fc_vdiv : real := udc/real(imax+1);
    begin
        retval := vref/fc_vdiv;
        for i in 1 to imax loop
            if vref > real(i)*fc_vdiv then
                retval := (vref - real(i)*fc_vdiv)/(fc_vdiv);
            end if;
        end loop;
        return retval;
    end get_fc_duty;
    ----------------------
    function get_next_step_length(t_sw : real; pwm : bit; duty : real) return real is
    begin
        if pwm = '1' then
            return t_sw * duty;
        else
            return t_sw * (1.0 - duty);
        end if;
    end get_next_step_length;
    ----------------------

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
                write_plot_config(file_handler, "T_title", "Inductor current and flying-cap balancing");
                write_plot_config(file_handler, "T_ylabel", "Current [A]");
                write_plot_config(file_handler, "B_title", "Bridge, output and flying-cap voltages");
                write_plot_config(file_handler, "B_ylabel", "Voltage [V]");
                write_plot_config(file_handler, "label_T_i0", "Inductor current");
                write_plot_config(file_handler, "label_T_i1", "FC1 balancing integrator");
                write_plot_config(file_handler, "label_T_i2", "FC2 balancing integrator");
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
                ,"T_i0" ,"T_i1" ,"T_i2"
                ,"B_u0" ,"B_u1" ,"B_u2" ,"B_u3" ,"B_u4"
                ));

                sim_ready          <= true;   -- kick off the first handshake
                simulation_started := true;

        elsif modulation_ready then
            simulation_counter <= simulation_counter + 1;
            steplen := mod_step;

            write_to(file_handler,(realtime
                    , lcr_rk5(0)         -- T_i0 : inductor current
                    , fc1_balance        -- T_i1 : FC1 balancing integrator
                    , fc2_balance        -- T_i2 : FC2 balancing integrator
                    , lcr_rk5(1)         -- B_u0 : output voltage
                    , lcr_rk5(2)         -- B_u1 : flying-cap 1
                    , lcr_rk5(3)         -- B_u2 : flying-cap 2
                    , udc_sig            -- B_u3 : DC link
                    , mod_ref            -- B_u4 : modulator reference
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
    p_modulation : process(simulator_clock)

        constant voltage_offset : real := 66.666/2.0;

        variable udc : real := 200.0;
        variable seed1, seed2 : positive := 1;
        variable rand : real;
        variable modulator_reference : real := 0.0;

        variable level_bits        : bit_vector(2 downto 0) := (others => '0');
        variable ones_in_low_state : natural range 0 to 2 := 0;

        variable avg_current, avg_current_diff : real := 0.0;
        variable prev_current, prev_avg_current : real := 0.0;
        variable sw_integ : real_vector(0 to 1) := (others => 0.0);

        variable pwm     : bit := '1';
        variable fc_duty : real := 0.5;

        type sw_vector    is array (natural range <>) of bit_vector;
        type sw_seq_matrix is array (natural range <>) of sw_vector;
        variable state_index : natural range 0 to 5 := 0;
        constant fc_4_sw_matrix : sw_seq_matrix(0 to 2)(0 to 5)(0 to 2) := (
            0 => ("001", "000", "010", "000", "100", "000"),
            1 => ("011", "010", "110", "100", "101", "001"),
            2 => ("111", "101", "111", "110", "111", "011"));
        variable next_sw_state : sw_states;

    begin
        if rising_edge(simulator_clock) then
            modulation_ready <= false;

            if sim_ready then
                -- dithered modulator reference, with scheduled steps
                uniform(seed1, seed2, rand);
                rand := ((rand - 0.5) * 2.0) * 1.0;
                modulator_reference := voltage_offset + rand;
                if realtime > 150.0e-3 then modulator_reference := 66.6666;  end if;
                if realtime > 250.0e-3 then modulator_reference := 140.9999; end if;

                -- scheduled DC-link steps
                if realtime > 100.0e-3 then udc := 150.0; end if;
                if realtime > 300.0e-3 then udc := 300.0; end if;

                -- output voltage level of the reference
                level_bits := (others => '0');
                if modulator_reference >= udc*0.0/3.0 then level_bits(0) := '1'; end if;
                if modulator_reference >= udc*1.0/3.0 then level_bits(1) := '1'; end if;
                if modulator_reference >= udc*2.0/3.0 then level_bits(2) := '1'; end if;
                ones_in_low_state := number_of_ones(level_bits) - 1;

                -- fed-back inductor current -> average current change
                avg_current      := (il_meas - prev_current)/2.0;
                avg_current_diff := avg_current - prev_avg_current;
                prev_avg_current := avg_current;
                prev_current     := il_meas;

                -- flying-cap balancing integrators (leaky), one per cell
                CASE mod_sw_state(2 downto 1) is
                    WHEN "01" => sw_integ(1) := sw_integ(1) - sw_integ(1)*80.0e-3 + avg_current_diff;
                    WHEN "10" => sw_integ(1) := sw_integ(1) - sw_integ(1)*80.0e-3 - avg_current_diff;
                    WHEN others => null;
                end CASE;
                CASE mod_sw_state(1 downto 0) is
                    WHEN "01" => sw_integ(0) := sw_integ(0) - sw_integ(0)*80.0e-3 + avg_current_diff;
                    WHEN "10" => sw_integ(0) := sw_integ(0) - sw_integ(0)*80.0e-3 - avg_current_diff;
                    WHEN others => null;
                end CASE;

                -- walk the switch-state matrix
                if state_index mod 2 = 0 then pwm := '1'; else pwm := '0'; end if;
                next_sw_state := fc_4_sw_matrix(ones_in_low_state)(state_index);
                if state_index < 5 then
                    state_index := state_index + 1;
                else
                    state_index := 0;
                end if;

                fc_duty := get_fc_duty(modulator_reference, udc, level_bits);

                mod_sw_state     <= next_sw_state;
                mod_step         <= get_next_step_length(t_sw*2.0, pwm, fc_duty);
                udc_sig          <= udc;
                mod_ref          <= modulator_reference;
                fc1_balance      <= sw_integ(0);
                fc2_balance      <= sw_integ(1);
                modulation_ready <= true;
            end if;
        end if;
    end process p_modulation;
------------------------------------------------------------------------
end vunit_simulation;
