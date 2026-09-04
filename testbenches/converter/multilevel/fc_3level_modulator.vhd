----------------------------------
-- 3-level flying-capacitor modulator (2-cell stack).
--
-- Thin per-level wrapper : the record ports come from fc_3level_modulator_pkg,
-- the actual work (to_level, get_sequence_times) from fc_modulator_common_pkg.
-- Clock is a separate port; there are no generics, the switching period is
-- carried on modulator_in.t_sw so it can change during simulation.
--
--   in   modulator_in.modulation_requested - strobe : solve the next sub-interval
--   in   modulator_in.duty_ratio           - 0.0 .. 1.0 over the whole 0..Udc range ;
--                                            to_level() picks the output level (row of
--                                            fc_sw_matrix), the fractional part is the
--                                            within-level duty that sets the nominal dwell
--   in   modulator_in.switching_time_trim  - one added dwell time per sequence state
--                                            (flying-cap balancing, see get_fc_trims) ;
--                                            leave at 0 for an open-loop run
--   in   modulator_in.t_sw                 - nominal switching period [s]
--
--   out  modulator_out.modulation_ready    - pulses the cycle the outputs are valid
--   out  modulator_out.next_switch_pattern - gate pattern for the coming sub-interval
--   out  modulator_out.switching_time      - how long to apply it [s]
--   out  modulator_out.switching_times     - the dwell time of every state of the sequence
--   out  modulator_out.level_switch_vector - the active level's gate patterns
----------------------------------
LIBRARY ieee  ;
    USE ieee.std_logic_1164.all  ;

LIBRARY ode;
    use ode.fc_modulator_common_pkg.all;
    use ode.fc_3level_modulator_pkg.all;

entity fc_3level_modulator is
    port (
        clock         : in  std_logic;
        modulator_in  : in  fc_modulator_input_record;
        modulator_out : out fc_modulator_output_record := fc_modulator_output_init
    );
end entity fc_3level_modulator;

architecture rtl of fc_3level_modulator is
begin

    modulator : process(clock)
        variable state_index : natural range 0 to n_states-1 := 0;
        variable level       : natural range 0 to n_cells-1  := 0;
        variable times       : sw_time_vector := (others => 0.0);
    begin
        if rising_edge(clock) then
            modulator_out.modulation_ready <= false;

            if modulator_in.modulation_requested then
                level := to_level(modulator_in.duty_ratio, n_cells);
                times := get_sequence_times(modulator_in.duty_ratio, n_cells,
                                            modulator_in.t_sw,
                                            modulator_in.switching_time_trim);

                modulator_out.next_switch_pattern <= fc_sw_matrix(level)(state_index);
                modulator_out.switching_time      <= times(state_index);
                modulator_out.switching_times     <= times;
                modulator_out.level_switch_vector <= fc_sw_matrix(level);
                modulator_out.modulation_ready    <= true;

                if state_index < n_states-1 then
                    state_index := state_index + 1;
                else
                    state_index := 0;
                end if;
            end if;
        end if;
    end process modulator;

end architecture rtl;
