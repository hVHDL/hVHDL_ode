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
    -----------------------------------
    -- simulation specific signals ----

    signal realtime : real := 0.0;
    constant stoptime : real := 400.0e-3;

    ----------------------
    function fc_modulator
    (
        gate_signals : bit_vector
    )
    return real is
        variable retval : real;
    begin
        CASE gate_signals is
            WHEN "10" => retval := -1.0;
            WHEN "01" => retval := 1.0;
            WHEN others => retval := 0.0;
        end CASE;
        
        return retval;
    end fc_modulator;
    ----------------------
    function number_of_ones(vector : bit_vector) return natural is
        variable retval : natural := 0;
    begin
        for i in vector'range loop
            if vector(i) = '1'
            then
                retval := retval + 1;
            end if;
        end loop;
        return retval;
    end number_of_ones;
    ----------

    subtype sw_states is bit_vector(3 downto 0);
    ----------
    function get_fc_bridge_voltage(sw_state : sw_states ; udc : real; ufc : real_vector) return real is
        variable bridge_voltage : real := 0.0;
    begin

        for i in ufc'range loop
            bridge_voltage := bridge_voltage + fc_modulator(sw_state(i+1 downto i)) * ufc(i);
        end loop;
        bridge_voltage := bridge_voltage + fc_modulator('0' & sw_state(sw_state'high)) * udc;

        return bridge_voltage;

    end get_fc_bridge_voltage;

    ----------
    function get_fc_duty(vref : real; udc : real ; level_bits : bit_vector) return real is
        variable retval : real := 0.0;
        variable imax : natural := level_bits'high;
        constant fc_vdiv : real := udc/real(imax+1);
    begin
        retval := vref/fc_vdiv;

        -- find voltage level of vref
        for i in 1 to imax loop
            if vref >= real(i)*fc_vdiv
            then
                retval := (vref - real(i)*fc_vdiv)/(fc_vdiv);
            end if;
        end loop;

        return retval;

    end get_fc_duty;


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

        constant initial_dc_link : real := 200.0;
        constant initial_voltage_ref : real := 149.0;
        variable udc    : real := initial_dc_link;
        variable i_load : real := 10.1111;
        constant l      : real := 20.0e-6;
        constant c      : real := 10.0e-6;
        constant rl     : real := 100.0e-3;
        constant cfc    : real := 2.0e-6;

        variable sw_frequency : real := 500.0e3;
        variable t_sw : real := 1.0/sw_frequency;
        variable duty : real := 0.5;

        variable seed1, seed2 : positive := 1;
        variable rand : real;

        -- i_l, uc, ufc
        constant init_state_vector : real_vector := (
              0 => 0.0
            , 1 => initial_voltage_ref -- udc
            , 2 => initial_dc_link*1.0/4.0   -- fc1
            , 3 => initial_dc_link*2.0/4.0   -- fc2
            , 4 => initial_dc_link*3.0/4.0); -- fc3

        variable sw_state      : sw_states := "0001";
        variable next_sw_state : sw_states := "1110";
        variable prev_sw_state : sw_states := "1111";

        ----------
        impure function deriv_lcr(t : real; states : real_vector) return real_vector is
            variable retval : states'subtype := (others => 0.0);
            variable bridge_voltage : real := 0.0;
            alias il is states(0);
            alias uc is states(1);
            alias ufc1 is states(2);
            alias ufc2 is states(3);
            alias ufc3 is states(4);
        begin

            bridge_voltage :=  get_fc_bridge_voltage(sw_state, udc, (ufc1, ufc2, ufc3));

            retval(0) := (bridge_voltage - il * rl - uc) * (1.0/l);
            retval(1) := (il - i_load) * (1.0/c);
            retval(2) := -fc_modulator(sw_state(1 downto 0)) * il / cfc;
            retval(3) := -fc_modulator(sw_state(2 downto 1)) * il / cfc;
            retval(4) := -fc_modulator(sw_state(3 downto 2)) * il / cfc;

            return retval;

        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_lcr);

        variable lcr_rk5 : init_state_vector'subtype := init_state_vector;

        file file_handler : text open write_mode is "fc_5level_tb.dat";

        ------------------- modulator variables ----------------------
        variable steplength : real := t_sw * (duty);

        variable pwm : bit := '1';
        variable modulator_reference : real := 0.0;
        variable fc_duty : real := 0.5;

        ---------------- end modulator variables ---------------------
        ---------------------------------------------------------------
        function get_next_step_length(t_sw : real; pwm : bit; duty : real) return real
        is
            variable retval : real := 1.0e-9;
            variable high_time : real := 1.0e-9;
            variable low_time : real := 1.0e-9;
        begin
            high_time := t_sw * duty;
            low_time  := t_sw * (1.0-duty);

            if pwm = '1'
            then
                retval := high_time;
            else
                retval := low_time;
            end if;

            return retval;

        end get_next_step_length;
        ---------------------------------------------------------------

        variable level_bits : bit_vector(3 downto 0) := (others => '0');
        variable ones_in_high_state : natural := 0;
        variable ones_in_low_state : natural := 0;

        type sw_vector is array (natural range <>) of bit_vector;
        type sw_matrix is array (natural range <>) of sw_vector;

        variable state_index : natural range 0 to 7 := 0;
        constant fc_5_sw_matrix : sw_matrix(0 to 3)(0 to 7)(0 to 3) := (
            0 =>(
                "0001",
                "0000",
                "0100",
                "0000",
                "0010",
                "0000",
                "1000",
                "0000"),
            1 =>(
                "0011",  -- 1-c2    | -c2
                "0001",  -- c3-c2   | c3-2c2
                "1001",  -- c3-c1   | 2c3-c2c-c1
                "1000",  -- c2-c1   | 2c3-c2-2c1
                "1100",  -- c2      | 2c3-2c1
                "0100",  -- c1      | 2c3-c1
                "0110",  -- 1-c3+c1 | c3
                "0010"), -- 1-c3    | 0
            2 =>(
                "0111",  -- c3      | c3
                "0011",  -- c2      | c3+c2
                "1011",  -- 1-c3+c2 | 2c2
                "1001",  -- 1-c3+c1 | -c3+2c2+c1
                "1101",  -- 1-c2+c1 | -c3+c2+2c1
                "1100",  -- 1-c2    | -c3+2c1
                "1110",  -- 1-c1    | -c3+c1
                "0110"), -- c3-c1   | 0
            3 =>(
                "1111",
                "1110",
                "1111",
                "1101",
                "1111",
                "1011",
                "1111",
                "0111"));
        variable avg_current : real := 10.0;
        variable avg_current_diff : real := 10.0;
        variable prev_current : real := -10.0;
        variable prev_avg_current : real := 10.0;
        variable voltage_offset : real := 200.0/4.0/2.0;

        variable sw_integ : real_vector(0 to 7) := (others => 0.0);

        variable sw_times : real_vector(0 to 7) := (1 => 0.00, others => 0.0);
        function sign (a : real) return real is
            variable retval :real := 0.0;
        begin
            if a >= 0.0 then
                retval := 1.0;
            else
                retval := -1.0;
            end if;
            return retval;
        end sign;

        constant fc_kp : real := 0.0;
        constant fc_ki : real := 0.00;
        variable pr : real := 0.0;

    begin
        if rising_edge(simulator_clock) then
            simulation_counter <= simulation_counter + 1;
            if simulation_counter = 0 then
                init_simfile(file_handler, ("time"
                ,"T_i0"
                ,"T_i1"
                ,"T_i2"
                ,"T_i3"
                ,"T_i4"
                ,"T_i5"
                ,"B_u0"
                -- ,"B_u1"
                -- ,"B_u2"
                -- ,"B_u3"
                -- ,"B_u4"
                ));
            end if;
            -------------------------

            write_to(file_handler,(realtime
                    -- ,lcr_rk5(0)          -- ,"T_i0"
                    -- ,avg_current_diff
                    -- ,sw_integ(0) + udc/3.0    
                    -- ,sw_integ(1) + udc*2.0/3.0
                    ,lcr_rk5(2) -- ,"B_u1"
                    ,lcr_rk5(3) -- ,"B_u1"
                    ,lcr_rk5(4)          -- ,"B_u2"
                    ,udc
                    ,modulator_reference -- ,"B_u4"
                    ,lcr_rk5(1)          -- ,"B_u0"
                    -- ,get_fc_bridge_voltage(sw_state, udc, ufc => (0 => lcr_rk5(2), 1 => lcr_rk5(3), 2 => lcr_rk5(4)))
                    ,lcr_rk5(0)          -- ,"T_i0"
                    -- ,avg_current
                    -- ,avg_current_diff
                    -- ,lcr_rk5(3)          -- ,"B_u2"
                    -- ,lcr_rk5(1)          -- ,"B_u0"
                    -- ,lcr_rk5(2)          -- ,"B_u1"
                    -- ,lcr_rk5(3)          -- ,"B_u2"
                    -- ,modulator_reference -- ,"B_u4"
                    -- ,sw_integ(2)
                    -- ,sw_integ(3)
                    -- ,sw_integ(4)
                ));

            rk5(realtime, lcr_rk5, steplength);
            realtime <= realtime + steplength;

            uniform(seed1, seed2, rand);
            rand := ((rand - 0.5) * 2.0) * 1.0;

                modulator_reference := initial_voltage_ref - lcr_rk5(0)*pr;
                -- if realtime > 150.0e-3 then
                --     modulator_reference := initial_voltage_ref*2.0/3.0 - lcr_rk5(0)*pr;
                --     modulator_reference := initial_voltage_ref+60.0 - lcr_rk5(0)*pr;
                    -- modulator_reference := (udc/2.0-abs((realtime mod 100.0e-3)/100.0e-3 * udc-udc/2.0))*2.0;
                -- end if;
                -- if realtime > 250.0e-3 then
                    -- modulator_reference := initial_voltage_ref;
                -- end if;
                -- if realtime > 350.0e-3 then
                --     modulator_reference := 190.0;
                --     modulator_reference := initial_voltage_ref - lcr_rk5(0)*pr;
                -- end if;

            ------- modulator -----------

            --
            level_bits := (others => '0');
            if modulator_reference >= udc*0.0/4.0 then level_bits(0) := '1'; end if;
            if modulator_reference >= udc*1.0/4.0 then level_bits(1) := '1'; end if;
            if modulator_reference >= udc*2.0/4.0 then level_bits(2) := '1'; end if;
            if modulator_reference >= udc*3.0/4.0 then level_bits(3) := '1'; end if;
            --
            ones_in_low_state  := number_of_ones(level_bits)-1;

            avg_current := (lcr_rk5(0))/2.0;
            avg_current_diff := -(avg_current - prev_avg_current);
            prev_avg_current := avg_current;
            -- sw_integ(state_index) :=sw_integ(state_index) + avg_current;
            prev_current := lcr_rk5(0);

            -- if realtime > 250.0e-3 then udc := 190.0; end if;
            -- if realtime > 300.0e-3 then udc := 300.0; end if;

            if state_index mod 2 = 0 then
                sw_times((state_index)) := -(sign(avg_current_diff)) *avg_current_diff* fc_ki;
                pwm := '1';
            else
                sw_times(state_index) := (sign(avg_current_diff)) *avg_current_diff* fc_ki;
                pwm := '0';
            end if;
            if sw_times(state_index) <= 0.0 then
                sw_times(state_index) := 0.0;
            end if;

            fc_duty       := get_fc_duty(modulator_reference, udc, level_bits);
            if (ones_in_low_state = 1) or (ones_in_low_state = 2) then
                steplength    := get_next_step_length(t_sw, pwm, fc_duty) + (avg_current_diff*fc_kp + sw_times(state_index))*t_sw;
            else
                if state_index mod 2 = 1 then
                    steplength    := get_next_step_length(t_sw, pwm, fc_duty);
                else
                    steplength    := get_next_step_length(t_sw, pwm, fc_duty);
                end if;
            end if;

            next_sw_state := fc_5_sw_matrix(ones_in_low_state)(state_index);
            if state_index < 7 then
                state_index := state_index + 1;
            else
                state_index := 0;
            end if;

            prev_sw_state := sw_state; -- not needed at the moment
            sw_state      := next_sw_state;

            -- if realtime < 400.0e-3 and realtime + steplength >= 400.0e-3
            -- then
            --     lcr_rk5(1) := lcr_rk5(1) - 50.0;
            -- end if;

            -----------------------------

        end if; -- rising_edge
    end process stimulus;	

------------------------------------------------------------------------
end vunit_simulation;
