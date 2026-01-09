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
    -----------------------------------
    -- simulation specific signals ----

    signal realtime : real := 0.0;
    constant stoptime : real := 10.0e-3;

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
        constant rl     : real := 50.0e-3;
        constant cfc    : real := 4.25e-6;

        variable sw_frequency : real := 250.0e3;
        variable t_sw : real := 1.0/sw_frequency;
        variable duty : real := 0.6;

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

        -- i_l, uc, ufc
        constant init_state_vector : real_vector := (0 => 0.0, 1 => 150.0,  2 => 66.0, 3 => 132.0);

        subtype sw_states is bit_vector(2 downto 0);

        variable sw_state      : sw_states := "111";
        variable next_sw_state : sw_states := "110";
        variable prev_sw_state : sw_states := "101";

        ----------
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
        function get_next_sw_state(sw_state : sw_states; prev_state : sw_states) return sw_states is
            variable next_sw_state : sw_states;
        begin

            case number_of_ones(sw_state) is
                WHEN 3 => 
                    CASE prev_state is
                        WHEN "110" => next_sw_state := "101";
                        WHEN "101" => next_sw_state := "011";
                        WHEN "011" => next_sw_state := "110";
                        WHEN others => next_sw_state := "111";
                    end CASE;
                WHEN others   => next_sw_state := "111";
            end CASE;

            return next_sw_state;
        end get_next_sw_state;

        ----------
        impure function get_step_length return real is
            variable step_length : real := 1.0e-9;
        begin

            case sw_state is
                WHEN "111"  => step_length := t_sw * duty;
                WHEN others => step_length := t_sw * (1.0-duty);
            end CASE;

            return step_length;

        end get_step_length;
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
        impure function deriv_lcr(t : real; states : real_vector) return real_vector is
            variable retval : states'subtype := (others => 0.0);
            variable bridge_voltage : real := 0.0;
            alias il is states(0);
            alias uc is states(1);
            alias ufc1 is states(2);
            alias ufc2 is states(3);
        begin

            -- if t > 250.0e-6 then i_load := 10.0; end if;
            -- if t > 600.0e-6 then duty := 0.4; end if;
            if t > 5.0e-3 then udc := 170.0; end if;
            -- if t > 10.0e-3 then udc := 10.0; end if;

            bridge_voltage :=  get_fc_bridge_voltage(sw_state, udc, (ufc1, ufc2));

            retval(0) := (bridge_voltage - il * rl - uc) * (1.0/l);
            retval(1) := (il - i_load - uc/15.0) * (1.0/c);
            retval(2) := -fc_modulator(sw_state(1 downto 0)) * il / cfc;
            retval(3) := -fc_modulator(sw_state(2 downto 1)) * il / cfc;

            return retval;

        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_lcr);

        variable lcr_rk5 : init_state_vector'subtype := init_state_vector;

        file file_handler : text open write_mode is "fc_4level_tb.dat";

    begin
        if rising_edge(simulator_clock) then
            simulation_counter <= simulation_counter + 1;

            if simulation_counter = 0 then
                init_simfile(file_handler, ("time"
                ,"T_i0"
                -- ,"B_u0"
                ,"B_u1"
                ,"B_u2"
                ,"B_u3"
                ,"B_u4"
                ));
            end if;

            write_to(file_handler,(realtime
                    ,lcr_rk5(0) 
                    -- ,get_fc_bridge_voltage(prev_sw_state, udc, ufc => (0 => lcr_rk5(2), 1 => lcr_rk5(3)))
                    ,lcr_rk5(1) 
                    ,lcr_rk5(2) 
                    ,lcr_rk5(3) 
                    ,udc
                ));

            write_to(file_handler,(realtime
                    ,lcr_rk5(0) 
                    -- ,get_fc_bridge_voltage(sw_state, udc, ufc => (0 => lcr_rk5(2), 1 => lcr_rk5(3)))
                    ,lcr_rk5(1) 
                    ,lcr_rk5(2) 
                    ,lcr_rk5(3) 
                    ,udc
                ));

            rk5(realtime, lcr_rk5, get_step_length);

            realtime <= realtime + get_step_length;
            prev_sw_state := sw_state;
            sw_state      := next_sw_state;
            next_sw_state := get_next_sw_state(sw_state, prev_sw_state);

        end if; -- rising_edge
    end process stimulus;	
------------------------------------------------------------------------
end vunit_simulation;
