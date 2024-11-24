
LIBRARY ieee  ; 
    USE ieee.NUMERIC_STD.all  ; 
    USE ieee.std_logic_1164.all  ; 
    use ieee.math_real.all;

package ode_pkg is 

------------------------------------------
    procedure generic_rk1
    generic(impure function deriv (input : real_vector) return real_vector is <>)
    (
        state    : inout real_vector;
        stepsize : real
    );
------------------------------------------
    procedure generic_rk2
    generic(impure function deriv (input : real_vector) return real_vector is <>)
    (
        state    : inout real_vector;
        stepsize : real);
------------------------------------------
    type rk_adaptive_record is record
        previous_step : real;
    end record;
------------------------------------------
    procedure generic_adaptive_rk23
    generic(impure function deriv (input : real_vector) return real_vector is <>)
    (
        state    : inout real_vector;
        z_n1     : inout real_vector;
        simtime  : inout real;
        err      : inout real;
        stepsize : inout real);
------------------------------------------
    procedure generic_rk4
    generic(impure function deriv (input : real_vector) return real_vector is <>)
    (
        state    : inout real_vector;
        stepsize : real);
------------------------------------------
    type am_state_array is array(natural range <>) of REAL_VECTOR;
    subtype am_array is am_state_array(1 to 4)(0 to 1);

    procedure am2_generic
    generic(impure function deriv (input : real_vector) return real_vector is <>)
    (
        variable adams_steps : inout am_state_array;
        variable state       : inout real_vector;
        stepsize             : real);
------------------------------------------
    procedure am4_generic
    generic(impure function deriv (input : real_vector) return real_vector is <>)
    (
        variable adams_steps : inout am_state_array;
        variable state       : inout real_vector;
        stepsize             : real);
------------------------------------------

end package ode_pkg;

package body ode_pkg is

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
------------------------------------------
    procedure generic_rk1
    generic(impure function deriv (input : real_vector) return real_vector is <>)
    (
        state    : inout real_vector;
        stepsize : real

    ) is
    begin
        state := state + deriv(state)*stepsize;
    end generic_rk1;

------------------------------------------
    procedure generic_rk2
    generic(impure function deriv (input : real_vector) return real_vector is <>)
    (
        state    : inout real_vector;
        stepsize : real
    ) is
        type state_array is array(1 to 2) of real_vector(state'range);
        variable k : state_array;
    begin
        k(1) := deriv(state);
        k(2) := deriv(state + k(1) * stepsize/ 2.0);

        state := state + k(2)*stepsize;

    end generic_rk2;

------------------------------------------
    procedure generic_adaptive_rk23
    generic(impure function deriv (input : real_vector) return real_vector is <>)
    (
        state    : inout real_vector;
        z_n1     : inout real_vector;
        simtime  : inout real;
        err      : inout real;
        stepsize : inout real
    ) is
        type state_array is array(1 to 4) of real_vector(state'range);
        variable k : state_array;
        variable y_n1 : real_vector(state'range);

        variable tolerance : real := 1.0e-4;
        variable h         : real := stepsize;
        variable h_new     : real ;

    begin
        k(1) := z_n1;
        k(2) := deriv(state + k(1) * stepsize * 1.0/2.0);
        k(3) := deriv(state + k(2) * stepsize * 3.0/4.0);

        y_n1 := state + (k(1)*2.0/9.0 + k(2)*1.0/3.0 + k(3)*4.0/9.0) * stepsize;

        k(4) := deriv(y_n1);

        z_n1 := state + (k(1)*7.0/24.0 + k(2)*1.0/4.0 + k(3)*1.0/3.0 + k(4) * 1.0/8.0) * stepsize;

        err := abs((z_n1(0) - y_n1(0)));

        state := y_n1;

        if err < tolerance then
            simtime := simtime + stepsize;
        end if;

        -- cubic root
        h_new := h*cbrt(tolerance/err);
        if h_new < 100.0e-9 then
            h_new := 100.0e-9;
        end if;
        if h_new > 100.0e-6 then
            h_new := 100.0e-6;
        end if;
        stepsize := h_new;

    end generic_adaptive_rk23;

------------------------------------------
    procedure generic_rk4
    generic(impure function deriv (input : real_vector) return real_vector is <>)
    (
        state    : inout real_vector;
        stepsize : real
    ) is
        type state_array is array(1 to 4) of real_vector(state'range);
        variable k : state_array;
    begin
        k(1) := deriv(state);
        k(2) := deriv(state + k(1) * stepsize/ 2.0);
        k(3) := deriv(state + k(2) * stepsize/ 2.0);
        k(4) := deriv(state + k(3) * stepsize);

        state := state + (k(1) + k(2) * 2.0 + k(3) * 2.0 + k(4)) * stepsize/6.0;

    end generic_rk4;
------------------------------------------
------------------------------------------
    procedure am2_generic
    generic(impure function deriv (input : real_vector) return real_vector is <>)
    (
        variable adams_steps : inout am_state_array;
        variable state       : inout real_vector;
        stepsize             : real
    ) is
        alias k is adams_steps;
    begin
        k(2) := k(1);
        k(1) := deriv(state);

        state := state + (k(1)*3.0 - k(2)) * stepsize/2.0;
    end am2_generic;
------------------------------------------
    procedure am4_generic
    generic(impure function deriv (input : real_vector) return real_vector is <>)
    (
        variable adams_steps : inout am_state_array;
        variable state       : inout real_vector;
        stepsize             : real
    ) is
        alias k is adams_steps;
    begin
        k(4) := k(3);
        k(3) := k(2);
        k(2) := k(1);
        k(1) := deriv(state);

        state := state + (k(1)*55.0 - k(2)*59.0 + k(3)*37.0 - k(4)*9.0) * stepsize/24.0;
    end am4_generic;
------------------------------------------
end package body;
