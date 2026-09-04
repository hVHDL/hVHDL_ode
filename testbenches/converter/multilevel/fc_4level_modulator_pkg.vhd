----------------------------------
-- 4-level flying-capacitor topology data. The modulation code is shared,
-- in fc_modulator_common_pkg ; this package holds only what depends on the
-- 3-cell stack : the level-sized subtypes, the switch-state matrix and the
-- wrapped port records.
----------------------------------
LIBRARY ieee  ;
LIBRARY ode;
    use ode.fc_modulator_common_pkg.all;

package fc_4level_modulator_pkg is

    constant n_cells  : positive := 3;
    constant n_states : positive := 2*n_cells;   -- states per switch sequence

    -- gate pattern of the cells ; one dwell value per sequence state ;
    -- one output level's row of the sequence
    subtype sw_states      is bit_vector(n_cells-1 downto 0);
    subtype sw_time_vector is real_vector(0 to n_states-1);
    subtype sw_sequence    is sw_vector(0 to n_states-1)(n_cells-1 downto 0);

    -- one row per output level (0 .. n_cells-1), one column per sequence state
    constant fc_sw_matrix : sw_seq_matrix(0 to n_cells-1)(0 to n_states-1)(n_cells-1 downto 0) := (
        0 => ("001", "000", "010", "000", "100", "000"),
        1 => ("011", "010", "110", "100", "101", "001"),
        2 => ("111", "101", "111", "110", "111", "011"));

    -- modulator ports, wrapped. clock stays a separate port.
    type fc_modulator_input_record is record
        modulation_requested : boolean;          -- strobe : solve the next sub-interval
        duty_ratio           : real;             -- 0.0 .. 1.0 over the whole 0..Udc range
        switching_time_trim  : sw_time_vector;   -- per-state dwell trim (flying-cap balancing)
        t_sw                 : real;             -- nominal switching period [s], settable at run time
    end record;

    type fc_modulator_output_record is record
        modulation_ready    : boolean;           -- pulses the cycle the fields below are valid
        next_switch_pattern : sw_states;         -- gate pattern for the coming sub-interval
        switching_time      : real;              -- how long to apply it [s]
        switching_times     : sw_time_vector;    -- the whole sequence's dwell times
        level_switch_vector : sw_sequence;       -- the active level's gate patterns
    end record;

    constant fc_modulator_input_init : fc_modulator_input_record := (
        modulation_requested => false,
        duty_ratio           => 0.0,
        switching_time_trim  => (others => 0.0),
        t_sw                 => 1.0e-6);

    constant fc_modulator_output_init : fc_modulator_output_record := (
        modulation_ready    => false,
        next_switch_pattern => (others => '0'),
        switching_time      => 0.0,
        switching_times     => (others => 0.0),
        level_switch_vector => fc_sw_matrix(0));

end package fc_4level_modulator_pkg;
