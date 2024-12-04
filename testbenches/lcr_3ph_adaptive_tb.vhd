LIBRARY ieee  ; 
    USE ieee.NUMERIC_STD.all  ; 
    USE ieee.std_logic_1164.all  ; 
    use ieee.math_real.all;
    use std.textio.all;

library vunit_lib;
context vunit_lib.vunit_context;

    use work.write_pkg.all;
    use work.adaptive_ode_pkg.all;
    use work.ode_pkg.all;
    use work.lcr_models_pkg.all;

entity lcr_3ph_adaptive_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of lcr_3ph_adaptive_tb is

    constant clock_period      : time    := 1 ns;
    
    signal simulator_clock     : std_logic := '0';
    signal simulation_counter  : natural   := 0;
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

        variable timestep : real := 10.0e-9;

        variable i_load : real_vector (0 to 1) := (others => 0.0);
        variable uin    : real_vector (1 to 3) := (0.0 , -0.0 , 0.0);
        constant l      : real_vector (1 to 3) := (2 => 80.0e-6, others => 80.0e-6);
        constant c      : real_vector (1 to 3) := (others => 60.0e-6);
        constant r      : real_vector (1 to 3) := (others => 100.0e-3);

        ------------
        impure function deriv(states : real_vector) return real_vector is
        begin
            return deriv_lcr(states, i_load, uin, l, c, r);
        end deriv;

        procedure rk23 is new generic_adaptive_rk23 generic map(maxstep => 10.0e-3, deriv => deriv);
        ------------

        variable lcr_rk23 : real_vector(0 to 5) := (others => 0.0);

        file file_handler : text open write_mode is "lcr_3ph_adaptive_tb.dat";
        variable simtime : real := 0.0;

        variable err  : real ;
        variable z_n1 : real_vector(lcr_rk23'range);

    begin
        if rising_edge(simulator_clock) then
            simulation_counter <= simulation_counter + 1;
            if simulation_counter = 0 then
                init_simfile(file_handler, ("time" 
                ,"T_u0"
                ,"T_u1"
                ,"T_u2"
                ,"B_i0"
                ,"B_i1"
                ,"B_i2"
                ,"B_st"
                ));
            end if;

            simtime := realtime;
            uin := (sin((simtime*1000.0-1.0/3.0)*math_pi*2.0)
                    , sin(simtime*1000.0*math_pi*2.0)
                    , sin((simtime*1000.0 + 1.0/3.0)*math_pi*2.0));

            z_n1 := deriv(lcr_rk23);


            if simulation_counter > 0 then


                simtime := realtime;
                uin := (sin((simtime*1000.0-1.0/3.0)*math_pi*2.0)
                        , sin(simtime*1000.0*math_pi*2.0)
                        , sin((simtime*1000.0 + 1.0/3.0)*math_pi*2.0));

                rk23(lcr_rk23, z_n1 , simtime, err , timestep);

                -- if realtime > 5.0e-3 then i_load := (2.0, -1.0); end if;
                -- if realtime > 5.0 then i_load := (-20.0, 10.0); end if;

                realtime <= realtime + timestep;
                write_to(file_handler,(realtime
                        ,lcr_rk23(3) 
                        ,lcr_rk23(4) 
                        ,lcr_rk23(5) 
                        ,lcr_rk23(0) 
                        ,lcr_rk23(1) 
                        ,lcr_rk23(2) 
                        ,timestep * 100.0
                    ));

            end if;

        end if; -- rising_edge
    end process stimulus;	
------------------------------------------------------------------------
end vunit_simulation;
