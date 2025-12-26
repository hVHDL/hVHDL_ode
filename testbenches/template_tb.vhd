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
    constant timestep : real := 10.0e-6;
    constant stoptime : real := 5.0e-3;

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
        variable u_in   : real := 10.0;
        variable i_load : real := 0.0;
        constant l      : real := 100.0e-6;
        constant c      : real := 100.0e-6;
        constant rl      : real := 300.0e-3;

        variable sw_frequency : real := 16.0e3;
        variable t_sw : real := 1.0/sw_frequency;
        variable duty : real := 0.5;

        constant init_state_vector : real_vector := (0.0, 0.0);

        type sw_states is (dc, zero);
        variable sw_state      : sw_states := dc;
        variable next_sw_state : sw_states := zero;

        ----------
        impure function get_step_length return real is
            variable step_length : real := 1.0e-9;
        begin
            case sw_state is
                WHEN dc => 
                    step_length := t_sw * duty;
                WHEN zero => 
                    step_length := t_sw * (1.0-duty);
            end CASE;

            return step_length;

        end get_step_length;
        ----------
        impure function get_bridge_voltage(sw_state : sw_states) return real is
            variable bridge_voltage : real := 0.0;
        begin
            CASE sw_state is
                WHEN dc   => bridge_voltage := u_in;
                WHEN zero => bridge_voltage := 0.0;
            end CASE;

            return bridge_voltage;
        end get_bridge_voltage;

        ----------
        impure function deriv_lcr(t : real; states : real_vector) return real_vector is
            variable retval : init_state_vector'subtype := init_state_vector;
            variable bridge_voltage : real := 0.0;
        begin


            if realtime > 2.5e-3 then duty := 0.7; end if;

            bridge_voltage := get_bridge_voltage(sw_state);
            retval(0) := (bridge_voltage - states(0) * rl - states(1)) * (1.0/l);
            retval(1) := (states(0) - i_load) * (1.0/c);

            return retval;

        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_lcr);

        variable lcr_rk5 : init_state_vector'subtype := init_state_vector;

        file file_handler : text open write_mode is "template_tb.dat";
    begin
        if rising_edge(simulator_clock) then
            simulation_counter <= simulation_counter + 1;
            if simulation_counter = 0 then
                init_simfile(file_handler, ("time"
                ,"T_i0"
                ,"B_u0"
                ,"B_u1"
                ));
            end if;

            if simulation_counter > 0 then


                realtime <= realtime + get_step_length;

                write_to(file_handler,(realtime
                        ,lcr_rk5(0) 
                        ,get_bridge_voltage(next_sw_state)
                        ,lcr_rk5(1) 
                    ));

                write_to(file_handler,(realtime
                        ,lcr_rk5(0) 
                        ,get_bridge_voltage(sw_state)
                        ,lcr_rk5(1) 
                    ));

                rk5(realtime, lcr_rk5, get_step_length);

                sw_state := next_sw_state;
                if sw_state = dc
                then
                    next_sw_state := zero;
                else
                    next_sw_state := dc;
                end if;

            end if;

        end if; -- rising_edge
    end process stimulus;	
------------------------------------------------------------------------
end vunit_simulation;
