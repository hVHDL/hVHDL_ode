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
    -- generic_rk5 is the classic Dormand-Prince 5(4) method (DOPRI5), the pair
    -- behind MATLAB's ode45 / SciPy's RK45. It is the default fixed-step
    -- 5th-order solver here.
    procedure generic_rk5
    generic(impure function deriv (t : real; input : real_vector) return real_vector is <>)
    (
        t          : in real
        ; state    : inout real_vector
        ; stepsize : in real);
------------------------------------------
    -- generic_tsit5 is the Tsitouras 5(4) method (Tsit5), offered as an option.
    -- It has noticeably smaller truncation-error constants than DOPRI5, but its
    -- stability polynomial is slightly more amplifying on the imaginary axis;
    -- for very lightly damped resonant plants sampled near the explicit-RK step
    -- limit (e.g. the grid-inverter example in this repository) that difference
    -- can accumulate and diverge, so DOPRI5 is the safer default.
    procedure generic_tsit5
    generic(impure function deriv (t : real; input : real_vector) return real_vector is <>)
    (
        t          : in real
        ; state    : inout real_vector
        ; stepsize : in real);
------------------------------------------
    -- generic_vern7 is Verner's "most efficient" 7(6) method (Vern7), offered
    -- as an option for smooth, high-accuracy problems (e.g. continuous-time
    -- model verification). 9 derivative evaluations per step. It is NOT worth
    -- using on the PWM-switched converter models: their right-hand side is
    -- discontinuous at every switching instant, which drops any high-order
    -- method to first order across the step.
    procedure generic_vern7
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
    -- Tsitouras 5(4) method ("Tsit5"). Option; see the header comment for when
    -- to prefer it over generic_rk5 (DOPRI5).
    -- Coefficients: C. Tsitouras, "Runge-Kutta pairs of order 5(4) satisfying
    -- only the first column simplifying assumption", Comput. Math. Appl. 62
    -- (2011), as used by DifferentialEquations.jl. 6 effective derivative
    -- evaluations per step (same cost as Dormand-Prince); the 7th, FSAL stage
    -- is only needed for the embedded error estimate and is omitted here.
    procedure generic_tsit5
    generic(impure function deriv (t : real; input : real_vector) return real_vector is <>)
    (
        t          : in real
        ; state    : inout real_vector
        ; stepsize : in real
    )
    is
        constant a2 : real_vector := (0 => 0.161);
        constant a3 : real_vector := (-0.008480655492356989 , 0.335480655492357);
        constant a4 : real_vector := (2.8971530571054935    , -6.359448489975075  , 4.3622954328695815);
        constant a5 : real_vector := (5.325864828439257     , -11.748883564062828 , 7.4955393428898365 , -0.09249506636175525);
        constant a6 : real_vector := (5.86145544294642      , -12.92096931784711  , 8.159367898576159  , -0.071584973281401 , -0.028269050394068383);
        -- 5th-order weights (== Tsit5 a7 row, FSAL)
        constant b  : real_vector := (0.09646076681806523   , 0.01                , 0.4798896504144996 , 1.379008574103742  , -3.290069515436081 , 2.324710524099774);
        constant c  : real_vector := (0.0 , 0.161 , 0.327 , 0.9 , 0.9800255409045097 , 1.0 , 1.0);
        type state_array is array(1 to 6) of real_vector(state'range);
        variable k : state_array;
        alias h is stepsize;
    begin
        k(1) := deriv(t, state);

        k(2) := deriv(t + h*c(1), state +
            ( k(1)*a2(0)
            ) * h);

        k(3) := deriv(t + h*c(2), state +
            ( k(1)*a3(0)
            + k(2)*a3(1)
            ) * h);

        k(4) := deriv(t + h*c(3), state +
            ( k(1)*a4(0)
            + k(2)*a4(1)
            + k(3)*a4(2)
            ) * h);

        k(5) := deriv(t + h*c(4), state +
            ( k(1)*a5(0)
            + k(2)*a5(1)
            + k(3)*a5(2)
            + k(4)*a5(3)
            ) * h);

        k(6) := deriv(t + h*c(5), state +
            ( k(1)*a6(0)
            + k(2)*a6(1)
            + k(3)*a6(2)
            + k(4)*a6(3)
            + k(5)*a6(4)
            ) * h);

        state := state +
            ( k(1)*b(0)
            + k(2)*b(1)
            + k(3)*b(2)
            + k(4)*b(3)
            + k(5)*b(4)
            + k(6)*b(5)
            ) * h;

    end generic_tsit5;
------------------------------------------
    -- Verner "most efficient" 7(6) method ("Vern7"). Option; see the header
    -- comment for when to prefer it.
    -- Coefficients: J. H. Verner, "Numerically optimal Runge-Kutta pairs with
    -- interpolants", Numer. Algorithms 53 (2010); values taken from the Float64
    -- ("CompiledFloats") Vern7Tableau in OrdinaryDiffEq.jl. Verified here: every
    -- A-row sums to its node, and b.A^q.1 = 1/(q+1)! holds through q = 6.
    -- 10-stage pair; the fixed-step 7th-order propagator needs stages 1..9
    -- (stage 10 is FSAL / error-estimate only and is omitted).
    procedure generic_vern7
    generic(impure function deriv (t : real; input : real_vector) return real_vector is <>)
    (
        t          : in real
        ; state    : inout real_vector
        ; stepsize : in real
    )
    is
        constant c2 : real := 0.005;
        constant c3 : real := 0.10888888888888888;
        constant c4 : real := 0.16333333333333333;
        constant c5 : real := 0.4555;
        constant c6 : real := 0.6095094489978381;
        constant c7 : real := 0.884;
        constant c8 : real := 0.925;
        constant c9 : real := 1.0;

        constant a21 : real := 0.005;
        constant a31 : real := -1.07679012345679;
        constant a32 : real := 1.185679012345679;
        constant a41 : real := 0.04083333333333333;
        constant a43 : real := 0.1225;
        constant a51 : real := 0.6389139236255726;
        constant a53 : real := -2.455672638223657;
        constant a54 : real := 2.272258714598084;
        constant a61 : real := -2.6615773750187572;
        constant a63 : real := 10.804513886456137;
        constant a64 : real := -8.3539146573962;
        constant a65 : real := 0.820487594956657;
        constant a71 : real := 6.067741434696772;
        constant a73 : real := -24.711273635911088;
        constant a74 : real := 20.427517930788895;
        constant a75 : real := -1.9061579788166472;
        constant a76 : real := 1.006172249242068;
        constant a81 : real := 12.054670076253203;
        constant a83 : real := -49.75478495046899;
        constant a84 : real := 41.142888638604674;
        constant a85 : real := -4.461760149974004;
        constant a86 : real := 2.042334822239175;
        constant a87 : real := -0.09834843665406107;
        constant a91 : real := 10.138146522881808;
        constant a93 : real := -42.6411360317175;
        constant a94 : real := 35.76384003992257;
        constant a95 : real := -4.3480228403929075;
        constant a96 : real := 2.0098622683770357;
        constant a97 : real := 0.3487490460338272;
        constant a98 : real := -0.27143900510483127;

        constant b1 : real := 0.04715561848627222;
        constant b4 : real := 0.25750564298434153;
        constant b5 : real := 0.26216653977412624;
        constant b6 : real := 0.15216092656738558;
        constant b7 : real := 0.4939969170032485;
        constant b8 : real := -0.29430311714032503;
        constant b9 : real := 0.08131747232495111;

        type state_array is array(1 to 9) of real_vector(state'range);
        variable k : state_array;
        alias h is stepsize;
    begin
        k(1) := deriv(t, state);

        k(2) := deriv(t + c2*h, state +
            ( k(1)*a21
            ) * h);

        k(3) := deriv(t + c3*h, state +
            ( k(1)*a31
            + k(2)*a32
            ) * h);

        k(4) := deriv(t + c4*h, state +
            ( k(1)*a41
            + k(3)*a43
            ) * h);

        k(5) := deriv(t + c5*h, state +
            ( k(1)*a51
            + k(3)*a53
            + k(4)*a54
            ) * h);

        k(6) := deriv(t + c6*h, state +
            ( k(1)*a61
            + k(3)*a63
            + k(4)*a64
            + k(5)*a65
            ) * h);

        k(7) := deriv(t + c7*h, state +
            ( k(1)*a71
            + k(3)*a73
            + k(4)*a74
            + k(5)*a75
            + k(6)*a76
            ) * h);

        k(8) := deriv(t + c8*h, state +
            ( k(1)*a81
            + k(3)*a83
            + k(4)*a84
            + k(5)*a85
            + k(6)*a86
            + k(7)*a87
            ) * h);

        k(9) := deriv(t + c9*h, state +
            ( k(1)*a91
            + k(3)*a93
            + k(4)*a94
            + k(5)*a95
            + k(6)*a96
            + k(7)*a97
            + k(8)*a98
            ) * h);

        state := state +
            ( k(1)*b1
            + k(4)*b4
            + k(5)*b5
            + k(6)*b6
            + k(7)*b7
            + k(8)*b8
            + k(9)*b9
            ) * h;

    end generic_vern7;
------------------------------------------
    -- Dormand-Prince 5(4) ("DOPRI5"), the 5th-order formula of the pair - the
    -- method behind MATLAB's ode45 / SciPy's RK45. Default fixed-step 5th-order
    -- solver; see generic_tsit5 for a lower-error alternative.
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
