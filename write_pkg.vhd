library ieee;
    use ieee.std_logic_1164.all;
    use ieee.math_real.all;
    use std.textio.all;
    use ieee.numeric_std.all;

package write_pkg is

    subtype real_number_array is real_vector;
    type stringarray is array (natural range <>) of string;
------------------------------------------------------------------------
    procedure write_to (
        file filee : text;
        data_to_be_written : real_number_array);
------------------------------------------------------------------------
    procedure init_simfile (
        file filee : text;
        data_to_be_written : stringarray);
------------------------------------------------------------------------
    -- writes a "#CONFIG key=value" line that test_plot.py reads to
    -- configure plot titles/labels. Recognized keys: title, xlabel,
    -- <prefix>_title, <prefix>_ylabel (e.g. T_title, B_ylabel).
    procedure write_plot_config (
        file filee : text;
        key   : string;
        value : string);
------------------------------------------------------------------------
end package write_pkg;

package body write_pkg is

------------------------------------------------------------------------
    procedure write_to
    (
        file filee : text;
        data_to_be_written : real_number_array
        
    ) is
        variable row : line;
        constant number_of_characters_between_columns : integer := 30;
    begin
        
        for i in data_to_be_written'range loop
            write(row , data_to_be_written(i) , left , number_of_characters_between_columns);
        end loop;

        writeline(filee , row);
    end write_to;
------------------------------------------------------------------------
    procedure init_simfile
    (
        file filee : text;
        data_to_be_written : stringarray
    ) is
        variable row : line;
        constant number_of_characters_between_columns : integer := 30;
    begin
        
        for i in data_to_be_written'range loop
            write(row , data_to_be_written(i) , left , number_of_characters_between_columns);
        end loop;

        writeline(filee , row);
    end init_simfile;
------------------------------------------------------------------------
    procedure write_plot_config
    (
        file filee : text;
        key   : string;
        value : string
    ) is
        variable row : line;
    begin
        write(row, string'("#CONFIG "));
        write(row, key);
        write(row, string'("="));
        write(row, value);
        writeline(filee, row);
    end write_plot_config;
------------------------------------------------------------------------
end package body write_pkg;
