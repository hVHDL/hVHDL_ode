----------------------------------
-- Shared bits of the 4-level flying-capacitor modulation scheme, used by
-- both fc_4level_modulator (the entity) and fc_4level_tb (which models the
-- bridge). Mirrors fc_5level_modulator_pkg for the 3-cell stack.
----------------------------------
LIBRARY ieee  ;
    use ieee.math_real.all;

package fc_4level_modulator_pkg is

    -- gate pattern of the 3 flying-cap cells, "abc"
    subtype sw_states is bit_vector(2 downto 0);

    -- switch-state sequence : one row per output level (0..2), one column
    -- per state of the 6-state sequence
    type sw_vector     is array (natural range <>) of bit_vector;
    type sw_seq_matrix is array (natural range <>) of sw_vector;
    constant fc_4_sw_matrix : sw_seq_matrix(0 to 2)(0 to 5)(0 to 2) := (
        0 => ("001", "000", "010", "000", "100", "000"),
        1 => ("011", "010", "110", "100", "101", "001"),
        2 => ("111", "101", "111", "110", "111", "011"));

    -- one output level's row of the sequence : its 6 gate patterns
    subtype sw_sequence    is sw_vector(0 to 5)(0 to 2);
    -- one value per state of the 6-state sequence
    subtype sw_time_vector is real_vector(0 to 5);

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
        level_switch_vector : sw_sequence;       -- the active level's 6 gate patterns
    end record;

    constant fc_modulator_input_init : fc_modulator_input_record := (
        modulation_requested => false,
        duty_ratio           => 0.0,
        switching_time_trim  => (others => 0.0),
        t_sw                 => 1.0e-6);

    constant fc_modulator_output_init : fc_modulator_output_record := (
        modulation_ready    => false,
        next_switch_pattern => "001",
        switching_time      => 0.0,
        switching_times     => (others => 0.0),
        level_switch_vector => fc_4_sw_matrix(0));

    -- +1 when a cell drains its flying cap, -1 when it charges it, 0 otherwise
    function fc_modulator(gate_signals : bit_vector) return real;

    -- output level (matrix row) for a 0..1 duty over the whole 0..Udc range
    function to_level(duty : real) return natural;

    -- nominal dwell of a state : t_sw*duty on the PWM high half,
    -- t_sw*(1-duty) on the low half
    function get_next_step_length(t_sw : real; pwm : bit; duty : real) return real;

    -- flying-cap balancing : one dwell-time trim per state of the 6-state
    -- sequence. A cell driven "01" discharges its flying cap (positive
    -- inductor current), "10" charges it; so lengthen the discharging states
    -- of an over-charged cap and shorten the charging ones. fc_modulator sums
    -- to zero over the sequence for every cell -> level average and period hold.
    function get_fc_trims(
        level       : natural;
        vfc_err     : real_vector;
        t_sw        : real;
        fc_kt       : real := 1.0e-8;    -- cap error [V] -> per-state time trim [s]
        fc_trim_max : real := 0.30       -- trim clamp, as a fraction of t_sw
    ) return real_vector;

end package fc_4level_modulator_pkg;

package body fc_4level_modulator_pkg is

    function fc_modulator(gate_signals : bit_vector) return real is
        variable retval : real;
    begin
        CASE gate_signals is
            WHEN "10" => retval := -1.0;
            WHEN "01" => retval :=  1.0;
            WHEN others => retval := 0.0;
        end CASE;
        return retval;
    end fc_modulator;

    function to_level(duty : real) return natural is
        variable d : real := duty;
        variable r : integer;
    begin
        if d < 0.0 then d := 0.0; end if;
        if d > 1.0 then d := 1.0; end if;
        r := integer(floor(d * 3.0));
        if r > 2 then r := 2; end if;
        return r;
    end to_level;

    function get_next_step_length(t_sw : real; pwm : bit; duty : real) return real is
    begin
        if pwm = '1' then
            return t_sw * duty;
        else
            return t_sw * (1.0 - duty);
        end if;
    end get_next_step_length;

    function get_fc_trims(
        level       : natural;
        vfc_err     : real_vector;
        t_sw        : real;
        fc_kt       : real := 1.0e-8;
        fc_trim_max : real := 0.30
    ) return real_vector is
        variable pat    : sw_states;
        variable trim   : real;
        variable retval : real_vector(0 to 5) := (others => 0.0);
    begin
        for s in 0 to 5 loop
            pat := fc_4_sw_matrix(level)(s);
            for k in 0 to 1 loop
                trim := fc_modulator(pat(k+1 downto k)) * vfc_err(k) * fc_kt;
                if    trim >  t_sw*fc_trim_max then trim :=  t_sw*fc_trim_max;
                elsif trim < -t_sw*fc_trim_max then trim := -t_sw*fc_trim_max;
                end if;
                retval(s) := retval(s) + trim;
            end loop;
        end loop;
        return retval;
    end get_fc_trims;

end package body fc_4level_modulator_pkg;
