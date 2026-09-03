----------------------------------
-- 5-level flying-capacitor converter, switching model.
--
-- Split into two processes like the other converter testbenches :
--   stimulus     - integrates the LCR + flying-cap plant with a fixed-step
--                  rk5, one switch sub-interval per handshake, and feeds
--                  the inductor current and the three flying-cap voltages
--                  back to the modulator.
--   p_modulation - from the modulator reference picks the row of
--                  fc_5_sw_matrix and computes a switching time
--                  *separately for each of the 8 states* in the sequence :
--                  the nominal level dwell plus a per-cell trim that
--                  lengthens the states draining an over-charged flying
--                  cap and shortens the ones charging it. fc_modulator
--                  sums to zero over the sequence for every cell, so the
--                  level average and the sequence period are unchanged.
--
-- A flying-cap disturbance at 200 ms exercises that per-state balancing
-- (it settles far quicker than the topology's own natural balancing).
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

    subtype sw_states is bit_vector(3 downto 0);

    -- stimulus <-> p_modulation handshake
    signal sim_ready        : boolean := false;
    signal modulation_ready : boolean := false;
    -- p_modulation -> stimulus
    signal mod_sw_state : sw_states := "0001";
    signal mod_step     : real      := t_sw*0.5;
    signal udc_sig      : real      := initial_dc_link;
    signal mod_ref      : real      := 0.0;
    -- stimulus -> p_modulation
    signal il_meas      : real         := 0.0;
    signal fc_meas      : real_vector(0 to 2) :=
        (0 => initial_dc_link*1.0/4.0,
         1 => initial_dc_link*2.0/4.0,
         2 => initial_dc_link*3.0/4.0);

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
            if vref >= real(i)*fc_vdiv then
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
                    , udc_sig      -- T_i3 : DC link
                    , mod_ref      -- T_i4 : modulator reference
                    , lcr_rk5(1)   -- T_i5 : output voltage
                    , lcr_rk5(0)   -- B_u0 : inductor current
                ));

            rk5(realtime, lcr_rk5, steplen);
            realtime <= realtime + steplen;
            il_meas    <= lcr_rk5(0);
            fc_meas(0) <= lcr_rk5(2);
            fc_meas(1) <= lcr_rk5(3);
            fc_meas(2) <= lcr_rk5(4);

            -- flying-cap disturbance : bump cell 2 off its udc/2 setpoint
            -- so the per-state switching times have something to correct
            if realtime < 200.0e-3 and realtime + steplen >= 200.0e-3 then
                lcr_rk5(3) := lcr_rk5(3) + 15.0;
            end if;

            sim_ready <= true;

            end if; -- handshake
        end if; -- rising_edge
    end process stimulus;
------------------------------------------------------------------------
    p_modulation : process(simulator_clock)

        constant pr           : real := 0.0;
        constant fc_kt        : real := 1.0e-8;   -- cap error [V] -> per-state time trim [s]
        constant fc_trim_max  : real := 0.30;     -- trim clamp, as a fraction of t_sw

        variable udc  : real := initial_dc_link;
        variable modulator_reference : real := 0.0;

        variable level_bits        : bit_vector(3 downto 0) := (others => '0');
        variable ones_in_low_state : natural range 0 to 3 := 0;

        variable pwm     : bit := '1';
        variable fc_duty : real := 0.5;

        type sw_vector     is array (natural range <>) of bit_vector;
        type sw_seq_matrix is array (natural range <>) of sw_vector;
        variable state_index : natural range 0 to 7 := 0;
        constant fc_5_sw_matrix : sw_seq_matrix(0 to 3)(0 to 7)(0 to 3) := (
            0 => ("0001", "0000", "0100", "0000", "0010", "0000", "1000", "0000"),
            1 => ("0011", "0001", "1001", "1000", "1100", "0100", "0110", "0010"),
            2 => ("0111", "0011", "1011", "1001", "1101", "1100", "1110", "0110"),
            3 => ("1111", "1110", "1111", "1101", "1111", "1011", "1111", "0111"));
        variable next_sw_state : sw_states;

        -- one switching time per state of the 8-state sequence
        variable sw_times : real_vector(0 to 7) := (others => t_sw*0.5);
        variable vfc_err  : real_vector(0 to 2) := (others => 0.0);
        variable trim     : real := 0.0;

    begin
        if rising_edge(simulator_clock) then
            modulation_ready <= false;

            if sim_ready then
                modulator_reference := initial_voltage_ref - il_meas*pr;

                -- output voltage level of the reference
                level_bits := (others => '0');
                if modulator_reference >= udc*0.0/4.0 then level_bits(0) := '1'; end if;
                if modulator_reference >= udc*1.0/4.0 then level_bits(1) := '1'; end if;
                if modulator_reference >= udc*2.0/4.0 then level_bits(2) := '1'; end if;
                if modulator_reference >= udc*3.0/4.0 then level_bits(3) := '1'; end if;
                ones_in_low_state := number_of_ones(level_bits) - 1;

                -- flying-cap voltage errors against the k*udc/4 setpoints
                for k in 0 to 2 loop
                    vfc_err(k) := fc_meas(k) - real(k+1)*udc/4.0;
                end loop;

                if state_index mod 2 = 0 then pwm := '1'; else pwm := '0'; end if;
                fc_duty       := get_fc_duty(modulator_reference, udc, level_bits);
                next_sw_state := fc_5_sw_matrix(ones_in_low_state)(state_index);

                -- switching time for THIS state : nominal level dwell, then a
                -- per-cell trim. A cell driven "01" discharges its flying cap
                -- (positive inductor current), "10" charges it; so lengthen the
                -- discharging states of an over-charged cap and shorten the
                -- charging ones. fc_modulator sums to zero over the 8-state
                -- sequence for every cell -> level average and period hold.
                sw_times(state_index) := get_next_step_length(t_sw, pwm, fc_duty);
                for k in 0 to 2 loop
                    trim := fc_modulator(next_sw_state(k+1 downto k)) * vfc_err(k) * fc_kt;
                    if    trim >  t_sw*fc_trim_max then trim :=  t_sw*fc_trim_max;
                    elsif trim < -t_sw*fc_trim_max then trim := -t_sw*fc_trim_max;
                    end if;
                    sw_times(state_index) := sw_times(state_index) + trim;
                end loop;
                if sw_times(state_index) < t_sw*0.01 then
                    sw_times(state_index) := t_sw*0.01;
                end if;

                mod_sw_state     <= next_sw_state;
                mod_step         <= sw_times(state_index);
                udc_sig          <= udc;
                mod_ref          <= modulator_reference;
                modulation_ready <= true;

                if state_index < 7 then
                    state_index := state_index + 1;
                else
                    state_index := 0;
                end if;
            end if;
        end if;
    end process p_modulation;
------------------------------------------------------------------------
end vunit_simulation;
