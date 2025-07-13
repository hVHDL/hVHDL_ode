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
    procedure generic_rk5
    generic(impure function deriv (t : real; input : real_vector) return real_vector is <>)
    (
        t          : in real
        ; state    : inout real_vector
        ; stepsize : in real);
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
    procedure generic_rk5
    generic(impure function deriv (t : real; input : real_vector) return real_vector is <>)
    (
        t          : in real
        ; state    : inout real_vector
        ; stepsize : in real
        -- ; y_n1     : inout real_vector
        -- ; fsal     : inout real_vector
        -- ; vErr     : out real_vector
    ) 
    is
        constant dop2 : real_vector := (0 => 1.0/5.0);
        constant dop3 : real_vector := (3.0/40.0       , 9.0/40.0);
        constant dop4 : real_vector := (44.0/45.0      , -56.0/15.0      , 32.0/9.0);
        constant dop5 : real_vector := (19372.0/6561.0 , -25360.0/2187.0 , 64448.0/6561.0 , -212.0/729.0);
        constant dop6 : real_vector := (9017.0/3168.0  , -355.0/33.0     , 46732.0/5247.0 , 49.0/176.0     , -5103.0/18656.0);
        constant dop7 : real_vector := (35.0/384.0     , 0.0             , 500.0/1113.0   , 125.0/192.0    , -2187.0/6784.0    , 11.0/84.0);

        constant dop8 : real_vector := (5179.0/57600.0 , 0.0     , 7571.0/16695.0 , 393.0/640.0 , -92097.0/339200.0 , 187.0/2100.0 , 1.0/40.0);

        constant tdop : real_vector := (0.0 , 1.0/5.0 , 3.0/10.0 , 4.0/5.0 , 8.0/9.0 , 1.0 , 1.0);
        type state_array is array(1 to 8) of real_vector(state'range);
        variable k : state_array;
        alias h is stepsize;

    begin
        -- k(1) := z;
        k(1) := deriv(t, state);

        k(2) := deriv(t + stepsize*tdop(1), state +
            ( k(1) * dop2(0) 
            ) * stepsize);

        k(3) := deriv(t + stepsize*tdop(2), state +
            ( k(1) * dop3(0)
            + k(2) * dop3(1)
            ) * stepsize);

        k(4) := deriv(t + stepsize*tdop(3), state +
            ( k(1) * dop4(0)
            + k(2) * dop4(1)
            + k(3) * dop4(2)
            ) * stepsize);

        k(5) := deriv(t + stepsize*tdop(4), state +
            ( k(1) * dop5(0)
            + k(2) * dop5(1)
            + k(3) * dop5(2)
            + k(4) * dop5(3)
            ) * stepsize);

        k(6) := deriv(t + stepsize*tdop(5), state +
            ( k(1) * dop6(0)
            + k(2) * dop6(1)
            + k(3) * dop6(2)
            + k(4) * dop6(3)
            + k(5) * dop6(4)
            ) * stepsize);

        state := 
            state +
            ( k(1) * dop7(0)
            + k(2) * dop7(1)
            + k(3) * dop7(2)
            + k(4) * dop7(3)
            + k(5) * dop7(4)
            + k(6) * dop7(5)
            ) * stepsize;

        -- k(7) := deriv(t + stepsize, y_n1);
        -- z := k(7);

        -- vErr := y_n1 - z;
        -- vErr := 
        --     ( k(1) * (dop8(0) - dop7(0))
        --     + k(2) * (dop8(1) - dop7(1))
        --     + k(3) * (dop8(2) - dop7(2))
        --     + k(4) * (dop8(3) - dop7(3))
        --     + k(5) * (dop8(4) - dop7(4))
        --     + k(6) * (dop8(5) - dop7(5))
        --     + k(7) *  dop8(6)
        --     ) * stepsize;
    end generic_rk5;
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
