library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.MATH_REAL.ALL;

package sort_generic_pkg is
    generic(type g_sorted_type);

    ------------------
    type event_record is record
        time_until_event : real;
        event_number : g_sorted_type;
    end record;

    type event_array is array (natural range <>) of event_record;
    ------------------
    function insertion_sort(arr : event_array) return event_array;
    ------------------
    function insertion_sort(arr : real_vector) return real_vector;
    ------------------
    function "=" (left : event_array; right : real_vector) return boolean;
    ------------------

end sort_generic_pkg;

package body sort_generic_pkg is

    ------------------
    function ">" (left : event_record; right : event_record) return boolean is
    begin
        return left.time_until_event > right.time_until_event;
    end function;
    ------------------
    function ">=" (left : event_record; right : event_record) return boolean is
    begin
        return left.time_until_event > right.time_until_event;
    end function;
    ------------------
    function "=" (left : event_array; right : real_vector) return boolean is
        variable compared_vector : real_vector(left'range);
    begin
        for i in left'range loop
            compared_vector(i) := left(i).time_until_event;
        end loop;

        return compared_vector = right;

    end function;
    ------------------
    function insertion_sort(arr : event_array) return event_array is
        variable sorted_arr : event_array(arr'range) := arr; -- Copy input array
        variable key : event_record;
        variable j : integer;
    begin

        for i in sorted_arr'left + 1 to sorted_arr'right loop
            key := sorted_arr(i);
            j := i - 1;
            
            -- Shift elements until correct position is found
            while j >= sorted_arr'left and sorted_arr(j).time_until_event > key.time_until_event loop
                sorted_arr(j + 1) := sorted_arr(j);
                j := j - 1;
            end loop;
            
            -- Insert key at the correct position
            sorted_arr(j + 1) := key;
        end loop;

        return sorted_arr;
    end function insertion_sort;
    ------------------
    function insertion_sort(arr : real_vector) return real_vector is
        variable sorted_arr : real_vector(arr'range) := arr; -- Copy input array
        variable key : real;
        variable j : integer;
    begin

        for i in sorted_arr'left + 1 to sorted_arr'right loop
            key := sorted_arr(i);
            j := i - 1;
            
            -- Shift elements until correct position is found
            while j >= sorted_arr'left and sorted_arr(j) > key loop
                sorted_arr(j + 1) := sorted_arr(j);
                j := j - 1;
            end loop;
            
            -- Insert key at the correct position
            sorted_arr(j + 1) := key;
        end loop;

        return sorted_arr;
    end function insertion_sort;
    ------------------

end sort_generic_pkg;

