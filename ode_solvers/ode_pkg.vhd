LIBRARY ieee  ; 
    USE ieee.NUMERIC_STD.all  ; 
    USE ieee.std_logic_1164.all  ; 
    use ieee.math_real.all;

package ode_pkg is 

------------------------------------------
    procedure generic_rk1
    generic(impure function deriv (t : real ; input : real_vector) return real_vector is <>)
    (
        t : real;
        state    : inout real_vector;
        stepsize : real
    );
------------------------------------------
    procedure generic_rk2
    generic(impure function deriv (t : real; input : real_vector) return real_vector is <>)
    (
        t : real;
        state    : inout real_vector;
        stepsize : real);
------------------------------------------
    type rk_adaptive_record is record
        previous_step : real;
    end record;
------------------------------------------
    procedure generic_rk4
    generic(impure function deriv (t : real; input : real_vector) return real_vector is <>)
    (
        t : real;
        state    : inout real_vector;
        stepsize : real);
------------------------------------------
    type am_state_array is array(natural range <>) of REAL_VECTOR;
    subtype am_array is am_state_array(1 to 4)(0 to 1);

    procedure am2_generic
    generic(impure function deriv (t : real; input : real_vector) return real_vector is <>)
    (
        t : real;
        variable adams_steps : inout am_state_array;
        variable state       : inout real_vector;
        stepsize             : real);
------------------------------------------
    procedure am4_generic
    generic(impure function deriv (t : real; input : real_vector) return real_vector is <>)
    (
        t : real;
        variable adams_steps : inout am_state_array;
        variable state       : inout real_vector;
        stepsize             : real);
------------------------------------------

end package ode_pkg;

package body ode_pkg is

    use work.real_vector_pkg.all;

------------------------------------------
    procedure generic_rk1
    generic(impure function deriv (t : real; input : real_vector) return real_vector is <>)
    (
        t : real;
        state    : inout real_vector;
        stepsize : real

    ) is
    begin
        state := state + deriv(t, state)*stepsize;
    end generic_rk1;

------------------------------------------
    procedure generic_rk2
    generic(impure function deriv (t : real; input : real_vector) return real_vector is <>)
    (
        t : real;
        state    : inout real_vector;
        stepsize : real
    ) is
        type state_array is array(1 to 2) of real_vector(state'range);
        variable k : state_array;
    begin
        k(1) := deriv(t, state);
        k(2) := deriv(t + stepsize/2.0, state + k(1) * stepsize/ 2.0);

        state := state + k(2)*stepsize;

    end generic_rk2;

------------------------------------------
    procedure generic_rk4
    generic(impure function deriv (t : real; input : real_vector) return real_vector is <>)
    (
        t : real;
        state    : inout real_vector;
        stepsize : real
    ) is
        type state_array is array(1 to 4) of real_vector(state'range);
        variable k : state_array;
    begin
        k(1) := deriv(t, state);
        k(2) := deriv(t + stepsize/2.0 , state + k(1) * stepsize/ 2.0);
        k(3) := deriv(t + stepsize/2.0 , state + k(2) * stepsize/ 2.0);
        k(4) := deriv(t + stepsize     , state + k(3) * stepsize);

        state := state + (k(1) + k(2) * 2.0 + k(3) * 2.0 + k(4)) * stepsize/6.0;

    end generic_rk4;
------------------------------------------
------------------------------------------
    procedure am2_generic
    generic(impure function deriv (t : real; input : real_vector) return real_vector is <>)
    (
        t : real;
        variable adams_steps : inout am_state_array;
        variable state       : inout real_vector;
        stepsize             : real
    ) is
        alias k is adams_steps;
    begin
        k(2) := k(1);
        k(1) := deriv(t, state);

        state := state + (k(1)*3.0 - k(2)) * stepsize/2.0;
    end am2_generic;
------------------------------------------
    procedure am4_generic
    generic(impure function deriv (t : real; input : real_vector) return real_vector is <>)
    (
        t : real;
        variable adams_steps : inout am_state_array;
        variable state       : inout real_vector;
        stepsize             : real
    ) is
        alias k is adams_steps;
    begin
        k(4) := k(3);
        k(3) := k(2);
        k(2) := k(1);
        k(1) := deriv(t, state);

        state := state + (k(1)*55.0 - k(2)*59.0 + k(3)*37.0 - k(4)*9.0) * stepsize/24.0;
    end am4_generic;
------------------------------------------
end package body;
