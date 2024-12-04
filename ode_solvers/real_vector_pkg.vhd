LIBRARY ieee  ; 
    USE ieee.NUMERIC_STD.all  ; 
    USE ieee.std_logic_1164.all  ; 
    use ieee.math_real.all;

package real_vector_pkg is

    function "+" (left : real_vector; right : real_vector) return real_vector;
    function "-" (left : real_vector; right : real_vector) return real_vector;
    function "-" (left : real_vector; right : real) return real_vector;
    function "/" (left : real_vector; right : real) return real_vector;
    function "*" (left : real_vector; right : real) return real_vector;

end package;


package body real_vector_pkg is
------------------------------------------
    function "+" (left : real_vector; right : real_vector) return real_vector is
        variable retval : left'subtype;
    begin

        for i in left'range loop
            retval(i) := left(i) + right(i);
        end loop;

        return retval;
    end function "+";
------------------------------------------
    function "-" (left : real_vector; right : real_vector) return real_vector is
        variable retval : left'subtype;
    begin

        for i in left'range loop
            retval(i) := left(i) - right(i);
        end loop;

        return retval;
    end function "-";
------------------------------------------
    function "-" (left : real_vector; right : real) return real_vector is
        variable retval : left'subtype;
    begin

        for i in left'range loop
            retval(i) := left(i) - right;
        end loop;

        return retval;
    end function "-";
------------------------------------------
    function "/" (left : real_vector; right : real) return real_vector is
        variable retval : left'subtype;
    begin

        for i in left'range loop
            retval(i) := left(i) / right;
        end loop;

        return retval;
    end function "/";
------------------------------------------
    function "*" (left : real_vector; right : real) return real_vector is
        variable retval : left'subtype;
    begin

        for i in left'range loop
            retval(i) := left(i) * right;
        end loop;

        return retval;
    end function "*";

end package body;
