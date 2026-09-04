----------------------------------
-- Common flying-capacitor modulation code, shared by the 3-, 4- and
-- 5-level testbenches. The level-specific packages
-- (fc_3level_modulator_pkg / fc_4level_modulator_pkg /
-- fc_5level_modulator_pkg) keep only the topology data : the cell count,
-- the level-sized subtypes, the switch-state matrix and the wrapped port
-- records. Everything that does not depend on the level lives here.
--
-- Conventions for an n-cell (n+1 level) stack :
--   * gate pattern       : bit_vector(n-1 downto 0), one bit per cell
--   * flying caps         : n-1 of them, ufc indexed 0 to n-2
--   * a switch sequence   : 2*n states, even states the PWM high half
--   * fc_sw_matrix        : row per output level 0..n-1, column per state
----------------------------------
LIBRARY ieee  ;
    use ieee.math_real.all;

package fc_modulator_common_pkg is

    -- a switch-state sequence : rows of gate-pattern vectors, and the matrix
    -- of them (one row per output level)
    type sw_vector     is array (natural range <>) of bit_vector;
    type sw_seq_matrix is array (natural range <>) of sw_vector;

    -- +1 when a cell drains its flying cap, -1 when it charges it, 0 otherwise
    function fc_modulator(gate_signals : bit_vector) return real;

    -- multilevel bridge voltage for a gate pattern : each cell's +/- flying-cap
    -- contribution plus the top cell's DC-link term. sw_state is "n-1 downto 0",
    -- ufc is "0 to n-2" (one entry per flying cap).
    function get_fc_bridge_voltage(
        sw_state : bit_vector;
        udc      : real;
        ufc      : real_vector) return real;

    -- nominal dwell of a state : t_sw*duty on the PWM high half,
    -- t_sw*(1-duty) on the low half
    function get_next_step_length(t_sw : real; pwm : bit; duty : real) return real;

    -- output level (matrix row, 0 .. n_cells-1) for a 0..1 duty over the
    -- whole 0..Udc range
    function to_level(duty : real; n_cells : positive) return natural;

    -- dwell time of every state of a 2*n_cells-state sequence : the nominal
    -- level time for that state plus the caller's per-state trim, floored at
    -- 1 % of t_sw so a large trim cannot collapse the period. The number of
    -- states is taken from trim'length.
    function get_sequence_times(
        duty_ratio : real;
        n_cells    : positive;
        t_sw       : real;
        trim       : real_vector) return real_vector;

    -- flying-cap balancing : one dwell trim per state of the sequence row
    -- `seq`. A cell driven "01" discharges its flying cap (positive inductor
    -- current), "10" charges it; so lengthen the discharging states of an
    -- over-charged cap and shorten the charging ones. fc_modulator sums to
    -- zero over the sequence for every cell -> level average and period hold.
    -- One cell per vfc_err entry.
    function get_fc_trims(
        seq         : sw_vector;
        vfc_err     : real_vector;
        t_sw        : real;
        fc_kt       : real := 1.0e-8;    -- cap error [V] -> per-state time trim [s]
        fc_trim_max : real := 0.30       -- trim clamp, as a fraction of t_sw
    ) return real_vector;

end package fc_modulator_common_pkg;

package body fc_modulator_common_pkg is

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

    function get_fc_bridge_voltage(
        sw_state : bit_vector;
        udc      : real;
        ufc      : real_vector) return real is
        variable bridge_voltage : real := 0.0;
    begin
        for i in ufc'range loop
            bridge_voltage := bridge_voltage + fc_modulator(sw_state(i+1 downto i)) * ufc(i);
        end loop;
        bridge_voltage := bridge_voltage + fc_modulator('0' & sw_state(sw_state'high)) * udc;
        return bridge_voltage;
    end get_fc_bridge_voltage;

    function get_next_step_length(t_sw : real; pwm : bit; duty : real) return real is
    begin
        if pwm = '1' then
            return t_sw * duty;
        else
            return t_sw * (1.0 - duty);
        end if;
    end get_next_step_length;

    function to_level(duty : real; n_cells : positive) return natural is
        variable d : real := duty;
        variable r : integer;
    begin
        if d < 0.0 then d := 0.0; end if;
        if d > 1.0 then d := 1.0; end if;
        r := integer(floor(d * real(n_cells)));
        if r > n_cells-1 then r := n_cells-1; end if;
        return r;
    end to_level;

    function get_sequence_times(
        duty_ratio : real;
        n_cells    : positive;
        t_sw       : real;
        trim       : real_vector) return real_vector is
        variable d, within : real;
        variable lvl       : natural;
        variable pwm       : bit;
        variable times     : real_vector(0 to trim'length-1) := (others => 0.0);
    begin
        d := duty_ratio;
        if d < 0.0 then d := 0.0; end if;
        if d > 1.0 then d := 1.0; end if;
        lvl    := to_level(d, n_cells);
        within := d * real(n_cells) - real(lvl);
        if within < 0.0 then within := 0.0; end if;
        if within > 1.0 then within := 1.0; end if;

        for s in times'range loop
            if s mod 2 = 0 then pwm := '1'; else pwm := '0'; end if;
            times(s) := get_next_step_length(t_sw, pwm, within) + trim(trim'low + s);
            if times(s) < t_sw*0.01 then
                times(s) := t_sw*0.01;
            end if;
        end loop;
        return times;
    end get_sequence_times;

    function get_fc_trims(
        seq         : sw_vector;
        vfc_err     : real_vector;
        t_sw        : real;
        fc_kt       : real := 1.0e-8;
        fc_trim_max : real := 0.30
    ) return real_vector is
        variable trim   : real;
        variable retval : real_vector(0 to seq'length-1) := (others => 0.0);
    begin
        for s in 0 to seq'length-1 loop
            for k in 0 to vfc_err'length-1 loop
                trim := fc_modulator(seq(seq'low + s)(k+1 downto k))
                        * vfc_err(vfc_err'low + k) * fc_kt;
                if    trim >  t_sw*fc_trim_max then trim :=  t_sw*fc_trim_max;
                elsif trim < -t_sw*fc_trim_max then trim := -t_sw*fc_trim_max;
                end if;
                retval(s) := retval(s) + trim;
            end loop;
        end loop;
        return retval;
    end get_fc_trims;

end package body fc_modulator_common_pkg;
