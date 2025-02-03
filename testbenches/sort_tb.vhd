
LIBRARY ieee  ; 
    USE ieee.NUMERIC_STD.all  ; 
    USE ieee.std_logic_1164.all  ; 
    use ieee.math_real.all;

library vunit_lib;
context vunit_lib.vunit_context;

    use work.sort_pkg.all;

entity sort_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of sort_tb is

    constant clock_period      : time    := 1 ns;

    signal simulator_clock     : std_logic := '0';
    signal simulation_counter  : natural   := 0;

    constant mixed_vector  : real_vector(0 to 4) := (2.0, 0.0, 1.0, 4.0, 3.0);
    constant sorted_vector : real_vector(0 to 4) := (0.0, 1.0, 2.0, 3.0, 4.0);

    signal test_vector : real_vector(mixed_vector'range) := mixed_vector;

begin

    simtime : process
    begin
        test_runner_setup(runner, runner_cfg);
        wait until simulation_counter = 5;
        check(test_vector = sorted_vector);
        test_runner_cleanup(runner); -- Simulation ends here
        wait;
    end process simtime;	

    simulator_clock <= not simulator_clock after clock_period / 2.0;

    stimulus : process(simulator_clock)

    begin
        if rising_edge(simulator_clock) then
            simulation_counter <= simulation_counter + 1;

            test_vector <= insertion_sort(mixed_vector);

        end if; -- rising_edge
    end process stimulus;	

end vunit_simulation;
