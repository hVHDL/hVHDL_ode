LIBRARY ieee  ; 
    USE ieee.NUMERIC_STD.all  ; 
    USE ieee.std_logic_1164.all  ; 
    use ieee.math_real.all;
    use std.textio.all;

library vunit_lib;
context vunit_lib.vunit_context;

    use work.write_pkg.all;
    use work.ode_pkg.all;

entity lcr_3ph_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of lcr_3ph_tb is

    constant clock_period      : time    := 1 ns;
    
    signal simulator_clock     : std_logic := '0';
    signal simulation_counter  : natural   := 0;
    -----------------------------------
    -- simulation specific signals ----

    signal realtime : real := 0.0;
    constant timestep : real := 10.0e-6;
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

        variable u_in : real := 10.0;
        variable i_load : real := 0.0;

        variable uin          : real_vector(1 to 3) := (0.0 , 0.0 , 0.0);
        variable state_vector : real_vector(0 to 5) := (others => 0.0);

        constant l : real_vector(1 to 3) := (others => 10.0e-6);

        impure function deriv_lcr (states : real_vector) return real_vector is

            variable retval : real_vector(0 to 1) := (0.0, 0.0);

            variable ul : real_vector(1 to 3) := (0.0 , 0.0 , 0.0);
            constant r  : real := 0.01;
            alias il    : real_vector(1 to 3) is state_vector(0 to 2);
            alias uc    : real_vector(1 to 3) is state_vector(3 to 5);

            variable un : real := 0.0;
            constant div : real := 1.0/(l(1)*l(2) + l(1)*l(3) + l(2)*l(3));
            constant a : real_vector(1 to 3) := (l(2)*l(3)/div, l(1)*l(3)/div, l(1)*l(2)/div);

            constant l      : real := 100.0e-6;
            constant c      : real := 100.0e-6;

            variable dil : real_vector(1 to 3);
            variable duc : real_vector(1 to 3);

        begin
            ul(1) := uin(1) - uc(1) - il(1) * r;
            ul(2) := uin(2) - uc(2) - il(2) * r;
            ul(3) := uin(3) - uc(3) - il(3) * r;
            un := a(1)*ul(1) + a(2)*ul(2) + a(3)*ul(3);

            dil(1) := (ul(1)-un)/l;
            dil(2) := (ul(2)-un)/l;
            dil(3) := (ul(3)-un)/l;

            duc(1) := (il(1))/c;
            duc(2) := (il(2))/c;
            duc(3) := (il(3))/c;

            retval(0) := (u_in - states(0) * 0.1 - states(1)) * (1.0/l);
            retval(1) := (states(0) - i_load) * (1.0/c);

            return retval;
        end function;

        procedure rk1 is new generic_rk1 generic map(deriv_lcr);
        procedure rk2 is new generic_rk2 generic map(deriv_lcr);

        variable k2 : am_array := (others => (others => 0.0));
        procedure am2 is new am2_generic generic map(deriv_lcr);

        variable k4 : am_array := (others => (others => 0.0));
        procedure am4 is new am4_generic generic map(deriv_lcr);

        variable lcr_rk1 : real_vector(0 to 1) := (0.0, 0.0);
        variable lcr_rk2 : real_vector(0 to 1) := (0.0, 0.0);

        variable lcr_am2     : real_vector(0 to 1) := (0.0, 0.0);
        variable lcr_am4     : real_vector(0 to 1) := (0.0, 0.0);

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
                "B_i2"
                ));
            end if;

            if simulation_counter > 0 then

                rk1(lcr_rk1, timestep);
                rk2(lcr_rk2, timestep);

                am2(k2,lcr_am2, timestep);
                am4(k4,lcr_am4, timestep);

                if realtime > 5.0e-3 then i_load := 2.0; end if;

                realtime <= realtime + timestep;
                write_to(file_handler,(realtime,
                        lcr_rk2(0) ,
                        lcr_am2(0) ,
                        lcr_am4(0) ,
                        lcr_rk2(1) ,
                        lcr_am2(1) ,
                        lcr_am4(1)
                    ));

            end if;

        end if; -- rising_edge
    end process stimulus;	
------------------------------------------------------------------------
end vunit_simulation;
