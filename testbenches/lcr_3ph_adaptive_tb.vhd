LIBRARY ieee  ; 
    USE ieee.NUMERIC_STD.all  ; 
    USE ieee.std_logic_1164.all  ; 
    use ieee.math_real.all;

package lcr_models_pkg is

    impure function deriv_lcr (
        states : real_vector
        ; i_load : real_vector
        ; uin : real_vector
        ; c : real_vector
        ; l : real_vector) 
        return real_vector;

end package;

package body lcr_models_pkg is

    impure function deriv_lcr (
        states : real_vector
        ; i_load : real_vector
        ; uin : real_vector
        ; c : real_vector
        ; l : real_vector) 
        return real_vector is

        variable retval : real_vector(0 to 5);

        variable ul : real_vector(1 to 3) := (0.0 , 0.0 , 0.0);
        constant r  : real := 0.11;
        alias il    : real_vector(1 to 3) is states(0 to 2);
        alias uc    : real_vector(1 to 3) is states(3 to 5);

        variable un  : real := 0.0;

        constant div : real                := 1.0/(l(1)*l(2) + l(1)*l(3) + l(2)*l(3));
        constant a   : real_vector(1 to 3) := (l(2)*l(3)/div, l(1)*l(3)/div, l(1)*l(2)/div);


        variable dil : real_vector(1 to 3);
        variable duc : real_vector(1 to 3);

    begin
        ul(1) := uin(1) - uc(1) - il(1) * r;
        ul(2) := uin(2) - uc(2) - il(2) * r;
        ul(3) := uin(3) - uc(3) - il(3) * r;
        un := a(1)*ul(1) + a(2)*ul(2) + a(3)*ul(3);

        dil(1) := (ul(1)-un)/l(1);
        dil(2) := (ul(2)-un)/l(2);
        dil(3) := (ul(3)-un)/l(3);

        duc(1) := (il(1) - i_load(0))/c(1);
        duc(2) := (il(2) - i_load(1))/c(2);
        duc(3) := (il(3) - (i_load(1) - i_load(0)))/c(3);

        retval := (dil(1), dil(2), dil(3), duc(1), duc(2), duc(3));

        return retval;

    end deriv_lcr;
end package body;
-----------------
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
    constant stoptime : real := 10.0;

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

        variable timestep : real := 10.0e-6;


        variable i_load : real_vector(0 to 1) := (others => 0.0);
        variable uin : real_vector(1 to 3) := (1.0 , -0.5 , 0.5);
        constant l : real_vector(1 to 3) := (2 => 50.0e-6, others => 10.0e-6);
        constant c : real_vector(1 to 3) := (others => 100.0e-6);

        ------------
        impure function deriv(states : real_vector) return real_vector is
        begin
            return deriv_lcr(states, i_load, uin, l, c);
        end deriv;

        procedure rk23 is new generic_adaptive_rk23 generic map(deriv);
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

            z_n1 := deriv(lcr_rk23);

            if simulation_counter > 0 then

                simtime := realtime;
                rk23(lcr_rk23, z_n1 , simtime, err , timestep);

                if realtime > 5.0e-3 then i_load := (2.0, -1.0); end if;
                if realtime > 5.0 then i_load := (-20.0, 10.0); end if;

                realtime <= realtime + timestep;
                write_to(file_handler,(realtime,
                        lcr_rk23(0) ,
                        lcr_rk23(1) ,
                        lcr_rk23(2) ,
                        lcr_rk23(3) ,
                        lcr_rk23(4) ,
                        lcr_rk23(5) ,
                        timestep
                    ));

            end if;

        end if; -- rising_edge
    end process stimulus;	
------------------------------------------------------------------------
end vunit_simulation;
