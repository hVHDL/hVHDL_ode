library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.MATH_REAL.ALL;

package sort_pkg is

    ------------------
    function insertion_sort(arr : real_vector) return real_vector;
    ------------------

end sort_pkg;

package body sort_pkg is

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

end sort_pkg;
