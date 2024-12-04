LIBRARY ieee  ; 
    USE ieee.NUMERIC_STD.all  ; 
    USE ieee.std_logic_1164.all  ; 
    use ieee.math_real.all;
    use std.textio.all;

library vunit_lib;
context vunit_lib.vunit_context;

    use work.write_pkg.all;
    use work.ode_pkg.all;
    use work.lcr_models_pkg.all;

entity lcr_3ph_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of lcr_3ph_tb is

    constant clock_period      : time    := 1 ns;
    
    signal simulator_clock     : std_logic := '0';
    signal simulation_counter  : natural   := 0;
    -----------------------------------
    -- simulation specific signals ----

    signal realtime   : real := 0.0;
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

        variable timestep : real := 1.0e-6;
        variable simtime : real := 0.0;

        variable i_load : real_vector (0 to 1) := (others => 0.0);
        variable uin    : real_vector (1 to 3) := (
                        sin(simtime*1000.0*math_pi*2.0)
                        ,sin((simtime*1000.0+1.0/3.0)*math_pi*2.0)
                        ,sin((simtime*1000.0 + 2.0/3.0)*math_pi*2.0));

        constant l      : real_vector (1 to 3) := (1 => 80.0e-6, others => 80.0e-6);
        constant c      : real_vector (1 to 3) := (others => 60.0e-6);
        constant r      : real_vector (1 to 3) := (others => 100.0e-3);

        ------------
        impure function deriv(states : real_vector) return real_vector is
        begin
            return deriv_lcr(states, i_load, uin, l, c, r);
        end deriv;

        procedure rk1 is new generic_rk1 generic map(deriv);
        procedure rk2 is new generic_rk2 generic map(deriv);
        procedure rk4 is new generic_rk2 generic map(deriv);

        variable k2 : am_state_array(1 to 4)(0 to 5) := (others => (others => 0.0));
        procedure am2 is new am2_generic generic map(deriv);

        variable k4 : am_state_array(1 to 4)(0 to 5) := (others => (others => 0.0));
        procedure am4 is new am4_generic generic map(deriv);

        variable lcr_rk1 : lcr_model_3ph_record := init_lcr_model;
        variable lcr_rk2 : lcr_model_3ph_record := init_lcr_model;
        variable lcr_rk4 : lcr_model_3ph_record := init_lcr_model;

        variable lcr_am2 : lcr_model_3ph_record := init_lcr_model;
        variable lcr_am4 : lcr_model_3ph_record := init_lcr_model;

        file file_handler : text open write_mode is "lcr_3ph_tb.dat";

    begin
        if rising_edge(simulator_clock) then
            simulation_counter <= simulation_counter + 1;
            if simulation_counter = 0 then
                init_simfile(file_handler, ("time", 
                "T_u0",
                "T_u1",
                "T_u2",
                "B_i0",
                "B_i1",
                "B_i2",
                "B_st"
                ));
            end if;

            if simulation_counter > 0 then

                simtime := realtime;
                write_to(file_handler,(realtime
                        ,get_capacitor_voltage(lcr_rk4)(0)
                        ,get_capacitor_voltage(lcr_rk4)(1)
                        ,get_capacitor_voltage(lcr_rk4)(2)
                        ,lcr_rk4.states(0)
                        ,lcr_rk4.states(1)
                        ,lcr_rk4.states(2) 
                        ,timestep
                    ));
                uin := (sin(simtime*1000.0*math_pi*2.0)
                       ,sin((simtime*1000.0+1.0/3.0)*math_pi*2.0)
                       ,sin((simtime*1000.0 + 2.0/3.0)*math_pi*2.0));

                rk1(lcr_rk1.states, timestep);
                rk2(lcr_rk2.states, timestep);
                rk4(lcr_rk4.states, timestep);

                am2(k2,lcr_am2.states, timestep);
                am4(k4,lcr_am4.states, timestep);

                realtime <= realtime + timestep;

            end if;

        end if; -- rising_edge
    end process stimulus;	
------------------------------------------------------------------------
end vunit_simulation;
