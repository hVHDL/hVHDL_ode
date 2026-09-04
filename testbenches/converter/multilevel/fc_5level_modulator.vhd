----------------------------------
-- 5-level flying-capacitor modulator.
--
-- Pulled out of fc_5level_tb.p_modulation into a clocked entity. Its
-- signals are wrapped into records (fc_modulator_input_record /
-- fc_modulator_output_record, in fc_5level_modulator_pkg); clock is a
-- separate port. There are no generics : the switching period is carried on
-- modulator_in.t_sw so it can be changed during simulation.
--
--   in   modulator_in.modulation_requested - strobe : solve the next sub-interval
--   in   modulator_in.duty_ratio           - one 0.0 .. 1.0 command over the whole
--                                            0 .. Udc output range. to_level() selects
--                                            the output level (row of fc_5_sw_matrix);
--                                            the fractional part is the within-level
--                                            duty that sets the nominal dwell.
--   in   modulator_in.switching_time_trim  - one added dwell time per state of the
--                                            8-state sequence. Supply flying-cap
--                                            balancing here (see get_fc_trims); leave
--                                            at 0 for an open-loop run.
--   in   modulator_in.t_sw                 - nominal switching period [s]
--
--   out  modulator_out.modulation_ready    - pulses the cycle the outputs are valid
--   out  modulator_out.next_switch_pattern - "abcd" gate pattern for the coming sub-interval
--   out  modulator_out.switching_time      - how long to apply it [s]
--   out  modulator_out.switching_times     - the dwell time of every state of the sequence
--   out  modulator_out.level_switch_vector - the active level's 8 gate patterns
--
-- Each request emits one state of the 8-state sequence : its switch pattern
-- from fc_5_sw_matrix and a dwell of the nominal level time
-- (get_next_step_length) plus switching_time_trim for that state, floored at
-- 1 % of t_sw so a large trim cannot collapse the period.
----------------------------------
LIBRARY ieee  ;
    USE ieee.std_logic_1164.all  ;

LIBRARY ode;
    use ode.fc_5level_modulator_pkg.all;

entity fc_5level_modulator is
    port (
        clock         : in  std_logic;
        modulator_in  : in  fc_modulator_input_record;
        modulator_out : out fc_modulator_output_record := fc_modulator_output_init
    );
end entity fc_5level_modulator;

architecture rtl of fc_5level_modulator is
begin

    modulator : process(clock)
        variable state_index : natural range 0 to 7 := 0;
        variable d           : real := 0.0;
        variable level       : natural range 0 to 3 := 0;
        variable within_duty : real := 0.0;
        variable pwm         : bit  := '1';
        variable t_sw        : real := 0.0;
        variable times       : sw_time_vector := (others => 0.0);
    begin
        if rising_edge(clock) then
            modulator_out.modulation_ready <= false;

            if modulator_in.modulation_requested then
                t_sw := modulator_in.t_sw;

                -- split the 0..1 command into a level (matrix row) and the
                -- within-level duty ratio that sets the nominal dwell
                d := modulator_in.duty_ratio;
                if d < 0.0 then d := 0.0; end if;
                if d > 1.0 then d := 1.0; end if;
                level := to_level(d);
                within_duty := d * 4.0 - real(level);
                if within_duty < 0.0 then within_duty := 0.0; end if;
                if within_duty > 1.0 then within_duty := 1.0; end if;

                -- nominal level dwell + the caller's trim, for every state of
                -- the sequence (even states are the PWM high half, odd the low)
                for s in 0 to 7 loop
                    if s mod 2 = 0 then pwm := '1'; else pwm := '0'; end if;
                    times(s) := get_next_step_length(t_sw, pwm, within_duty)
                                + modulator_in.switching_time_trim(s);
                    if times(s) < t_sw*0.01 then
                        times(s) := t_sw*0.01;
                    end if;
                end loop;

                modulator_out.next_switch_pattern <= fc_5_sw_matrix(level)(state_index);
                modulator_out.switching_time      <= times(state_index);
                modulator_out.switching_times     <= times;
                modulator_out.level_switch_vector <= fc_5_sw_matrix(level);
                modulator_out.modulation_ready    <= true;

                if state_index < 7 then
                    state_index := state_index + 1;
                else
                    state_index := 0;
                end if;
            end if;
        end if;
    end process modulator;

end architecture rtl;
