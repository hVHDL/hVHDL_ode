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
    constant stoptime : real := 1.0e-3;

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
        variable i_load : real := 0.0;
        constant l      : real := 1.0e-6;
        constant c      : real := 100.0e-6;
        constant rl     : real := 20.0e-3;
        constant cfc    : real := 40.0e-6;

        variable sw_frequency : real := 200.0e3;
        variable t_sw : real := 1.0/sw_frequency;
        variable duty : real := 0.5;

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
        constant init_state_vector : real_vector := (0 => 0.0, 1 => 0.0,  2 => udc*0.45);

        subtype sw_states is bit_vector(1 downto 0);

        constant dc     : bit_vector(1 downto 0) := "11";
        constant p_half : bit_vector(1 downto 0) := "10";
        constant n_half : bit_vector(1 downto 0) := "01";
        constant zero   : bit_vector(1 downto 0) := "00";

        type sw_state_record is record
            sw_state : sw_states;
            next_sw_state : sw_states;
        end record;

        constant init_sw_state : sw_state_record := (dc, p_half);

        variable sw_state      : sw_states := dc;
        variable next_sw_state : sw_states := p_half;
        variable prev_sw_state : sw_states := n_half;

        ----------
        function get_next_sw_state(sw_state : sw_states; prev_state : sw_states) return sw_states is
            variable next_sw_state : sw_states;
        begin

            case sw_state is
                WHEN dc => 
                    if prev_state = n_half then 
                        next_sw_state := p_half;
                    else
                        next_sw_state := n_half;
                    end if;

                WHEN p_half => next_sw_state := dc;
                WHEN n_half => next_sw_state := dc;
                WHEN zero   => next_sw_state := n_half;

            end CASE;

            return next_sw_state;
        end get_next_sw_state;

        variable fc_duty : real := 1.0;

        ----------
        impure function get_step_length return real is
            variable step_length : real := 1.0e-9;
        begin
            case sw_state is
                WHEN dc     => step_length := t_sw * duty;
                WHEN p_half => step_length := t_sw * (1.0-duty);
                WHEN n_half => step_length := t_sw * (1.0-duty);
                WHEN zero   => step_length := t_sw * (1.0-duty);
            end CASE;

            return step_length;

        end get_step_length;
        ----------
        impure function get_bridge_voltage(sw_state : sw_states ; udc : real; ufc : real_vector) return real is
            variable bridge_voltage : real := 0.0;
        begin
            CASE sw_state is
                WHEN dc     => bridge_voltage := udc;
                WHEN p_half => bridge_voltage := udc-ufc(0);
                WHEN n_half => bridge_voltage := ufc(0);
                WHEN zero   => bridge_voltage := 0.0;
            end CASE;

            return bridge_voltage;
        end get_bridge_voltage;

        ----------
        impure function deriv_lcr(t : real; states : real_vector) return real_vector is
            variable retval : init_state_vector'subtype := init_state_vector;
            variable bridge_voltage : real := 0.0;
        begin

            if t > 250.0e-6 then i_load := 10.0; end if;
            if t > 600.0e-6 then duty := 0.8; end if;

            bridge_voltage := 
                fc_modulator(('0', sw_state(1))) * udc
              + fc_modulator((sw_state(1), sw_state(0))) * states(2);

            retval(0) := (bridge_voltage - states(0) * rl - states(1)) * (1.0/l);
            retval(1) := (states(0)      - i_load) * (1.0/c);
            retval(2) := fc_modulator((sw_state(1), sw_state(0))) * states(0);

            return retval;

        end function;

        procedure rk5 is new generic_rk5 generic map(deriv_lcr);

        variable lcr_rk5 : init_state_vector'subtype := init_state_vector;

        file file_handler : text open write_mode is "template_tb.dat";
        use ode.real_vector_pkg.all;
    begin
        if rising_edge(simulator_clock) then
            simulation_counter <= simulation_counter + 1;
            if simulation_counter = 0 then
                init_simfile(file_handler, ("time"
                ,"T_i0"
                ,"B_u0"
                ,"B_u1"
                ,"B_u2"
                ));
            end if;


            realtime <= realtime + get_step_length;

            write_to(file_handler,(realtime
                    ,lcr_rk5(0) 
                    ,get_bridge_voltage(prev_sw_state, udc, ufc => (0 => lcr_rk5(2)))
                    ,lcr_rk5(1) 
                    ,lcr_rk5(2) 
                ));

            write_to(file_handler,(realtime
                    ,lcr_rk5(0) 
                    ,get_bridge_voltage(sw_state, udc, ufc => (0 => lcr_rk5(2)))
                    ,lcr_rk5(1) 
                    ,lcr_rk5(2) 
                ));

            rk5(realtime, lcr_rk5, get_step_length);

            prev_sw_state := sw_state;
            sw_state      := next_sw_state;
            next_sw_state := get_next_sw_state(sw_state, prev_sw_state);

        end if; -- rising_edge
    end process stimulus;	
------------------------------------------------------------------------
end vunit_simulation;
