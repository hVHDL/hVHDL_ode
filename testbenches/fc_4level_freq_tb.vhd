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

entity fc_4level_freq_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of fc_4level_freq_tb is

    constant clock_period : time := 1 ns;
    
    signal simulator_clock    : std_logic := '0';
    signal simulation_counter : natural   := 0;
    -----------------------------------
    -- simulation specific signals ----

    signal realtime : real := 0.0;
    constant stoptime : real := 500.0e-3;

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

    subtype sw_states is bit_vector(2 downto 0);
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
    function get_fc_duty(vref : real; udc : real ; imax : natural := 2) return real is
        variable retval : real := 0.0;
        constant fc_vdiv : real := udc/real(imax+1);
    begin
        retval := vref/fc_vdiv;

        -- find voltage level of vref
        for i in 1 to imax loop
            if vref > real(i)*fc_vdiv
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

        variable udc    : real := 200.0;
        variable i_load : real := -10.0;
        constant l      : real := 10.0e-6;
        constant c      : real := 10.0e-6;
        constant rl     : real := 10.0e-3;
        constant cfc    : real := 4.0e-6;

        variable sw_frequency : real := 500.0e3;
        variable t_sw : real := 1.0/sw_frequency;
        variable duty : real := 0.5;

        variable seed1, seed2 : positive := 1;
        variable rand : real;

        -- i_l, uc, ufc
        constant init_state_vector : real_vector := (
              0 => 0.0
            , 1 => 150.0
            , 2 => 66.0    -- fc1
            , 3 => 132.0); -- fc2

        variable sw_state      : sw_states := "111";
        variable next_sw_state : sw_states := "110";
        variable prev_sw_state : sw_states := "101";

        ----------
        impure function deriv_lcr(t : real; states : real_vector) return real_vector is
            variable retval : states'subtype := (others => 0.0);
            variable bridge_voltage : real := 0.0;
            alias il is states(0);
            alias uc is states(1);
            alias ufc1 is states(2);
            alias ufc2 is states(3);
        begin

            bridge_voltage :=  get_fc_bridge_voltage(sw_state, udc, (ufc1, ufc2));

            retval(0) := (bridge_voltage - il * rl - uc) * (1.0/l);
            retval(1) := (il ) * (1.0/c);
            retval(2) := -fc_modulator(sw_state(1 downto 0)) * il / cfc;
            retval(3) := -fc_modulator(sw_state(2 downto 1)) * il / cfc;

            return retval;

        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_lcr);

        variable lcr_rk5 : init_state_vector'subtype := init_state_vector;

        file file_handler : text open write_mode is "fc_4level_tb.dat";

        ------------------- modulator variables ----------------------
        variable high_time  : real := t_sw * (duty);
        variable low_time   : real := t_sw * (1.0-duty);

        variable steplength : real := t_sw * (duty);

        variable pwm : bit := '1';
        variable modulator_reference : real := 167.0;
        variable fc_duty : real := 0.5;

        variable next_high_state : BIT_VECTOR(2 downto 0) := "111";
        variable next_low_state  : BIT_VECTOR(2 downto 0) := "110";

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

        variable level_bits : bit_vector(2 downto 0) := (others => '0');
        variable ones_in_high_state : natural := 0;
        variable ones_in_low_state : natural := 0;

    begin
        if rising_edge(simulator_clock) then
            simulation_counter <= simulation_counter + 1;
            if simulation_counter = 0 then
                init_simfile(file_handler, ("time"
                ,"T_i0"
                ,"B_u0"
                ,"B_u1"
                ,"B_u2"
                ,"B_u3"
                ,"B_u4"
                ));
            end if;
            -------------------------

            write_to(file_handler,(realtime
                    ,lcr_rk5(0)          -- ,"T_i0"
                    ,lcr_rk5(1)          -- ,"B_u0"
                    ,lcr_rk5(2)          -- ,"B_u1"
                    ,lcr_rk5(3)          -- ,"B_u2"
                    ,udc                 -- ,"B_u3"
                    ,modulator_reference -- ,"B_u4"
                ));

            -- write_to(file_handler,(realtime
            --         ,lcr_rk5(0) 
            --         -- ,get_fc_bridge_voltage(sw_state, udc, ufc => (0 => lcr_rk5(2), 1 => lcr_rk5(3)))
            --         ,lcr_rk5(1) 
            --         ,lcr_rk5(2) 
            --         ,lcr_rk5(3) 
            --         ,udc
            --     ));

            rk5(realtime, lcr_rk5, steplength);
            realtime <= realtime + steplength;

            uniform(seed1, seed2, rand);
            rand := ((rand - 0.5) * 2.0) * 1.0;

            if simulation_counter mod 1 = 0
            then
                modulator_reference := 167.0 + rand;
            end if;

            ------- modulator -----------
            pwm := not pwm;

            --
            if modulator_reference >= udc*0.0/3.0 then level_bits(0) := '1'; end if;
            if modulator_reference >= udc*1.0/3.0 then level_bits(1) := '1'; end if;
            if modulator_reference >= udc*2.0/3.0 then level_bits(2) := '1'; end if;
            --
            ones_in_high_state := number_of_ones(level_bits);
            ones_in_low_state  := number_of_ones(level_bits)-1;

            fc_duty    := get_fc_duty(modulator_reference, udc);
            steplength := get_next_step_length(t_sw, pwm, fc_duty);

            prev_sw_state := sw_state;
            sw_state      := next_sw_state;

            -- ones_in_high_state
            -- ones_in_low_state 

            case number_of_ones(sw_state) is
                WHEN 3 => 
                    CASE prev_sw_state is
                        WHEN "110"  => next_sw_state := "101";
                        WHEN "101"  => next_sw_state := "011";
                        WHEN "011"  => next_sw_state := "110";
                        WHEN others => next_sw_state := "111";
                    end CASE;
                WHEN others   => next_sw_state := "111";
            end CASE;

            -----------------------------


        end if; -- rising_edge
    end process stimulus;	

    control : process(simulator_clock) is
    begin
        if rising_edge(simulator_clock)
        then 
        end if;
    end process;
------------------------------------------------------------------------
end vunit_simulation;
