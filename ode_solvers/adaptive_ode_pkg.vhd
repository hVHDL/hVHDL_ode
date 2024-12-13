LIBRARY ieee  ; 
    USE ieee.NUMERIC_STD.all  ; 
    USE ieee.std_logic_1164.all  ; 
    use ieee.math_real.all;

package adaptive_ode_pkg is 
    constant default_minstep : real := 10.0e-9;
    constant default_maxstep : real := 100.0e-4;
    constant default_tolerance : real := 1.0e-7;

------------------------------------------
    procedure generic_adaptive_rk23
    generic(
        impure function deriv (t : real; input : real_vector) return real_vector is <>
        ; tolerance : real := default_tolerance
        ; minstep : real := default_minstep
        ; maxstep : real := default_maxstep)
    (
        t : real;
        state    : inout real_vector;
        z_n1     : inout real_vector;
        err      : inout real;
        stepsize : inout real);

------------------------------------------
    procedure generic_adaptive_dopri54
    generic(impure function deriv (input : real_vector) return real_vector is <>)
    (
        state    : inout real_vector;
        z_n1     : inout real_vector;
        simtime  : inout real;
        err      : inout real;
        stepsize : inout real);
------------------------------------------

end package adaptive_ode_pkg;

package body adaptive_ode_pkg is 

    use work.real_vector_pkg.all;
------------------------------------------
    function norm(input : real_vector) return real is
        variable retval : real := 0.0;
    begin
        for i in input'range loop
            retval := retval + input(i)**2.0;
        end loop;

        return sqrt(retval);
    end norm;

------------------------------------------
    procedure generic_stepper
    generic(
        minstep : real := default_minstep
        ;maxstep : real := default_maxstep)
    (
        prev_stepsize : in real
        ; vErr : in real_vector
        ; h_new : inout real
        ; err : inout real)
    is
    begin
        err := norm(vErr); 

        if abs(err) > 1.0e-15 then
            h_new := prev_stepsize*cbrt(default_tolerance/err); -- cbrt() is cubic root
            if h_new < minstep then
                h_new := minstep;
            end if;
            if h_new > maxstep then
                h_new := maxstep;
            end if;
        else
            h_new := maxstep;
        end if;
    end generic_stepper;
------------------------------------------
    procedure generic_adaptive_rk23
    generic(
        impure function deriv (t : real; input : real_vector) return real_vector is <>;
        tolerance : real := default_tolerance;
        minstep : real := default_minstep;
        maxstep : real := default_maxstep
    )
    (
        t : real;
        state    : inout real_vector;
        z_n1     : inout real_vector;
        err      : inout real;
        stepsize : inout real
    ) is
        type st_array is array(natural range <>) of real_vector(state'range);
        subtype state_array is st_array(1 to 4);

        procedure rk(y_n1 : inout real_vector; state : real_vector; h: real; z : inout real_vector; vErr : inout real_vector; k : inout state_array) is
        begin
            k(1) := z;
            k(2) := deriv(t + h/2.0, state + k(1) * h * 1.0/2.0);
            k(3) := deriv(t + h*3.0/4.0, state + k(2) * h * 3.0/4.0);

            y_n1 := state + (k(1)*2.0/9.0 + k(2)*1.0/3.0 + k(3)*4.0/9.0) * h;

            k(4) := deriv(t + h, y_n1);

            z  := state + k(4) * h;

            vErr := (k(1)*(-5.0/72.0) + k(2)*1.0/12.0 + k(3)*1.0/9.0 - k(4)*1.0/8.0) * h;
        end rk;
        variable k : state_array;
        variable y_n1 : real_vector(state'range);

        variable h         : real := stepsize;
        variable h_new     : real ;
        variable vErr       : real_vector(state'range);
        variable z : real_vector(z_n1'range) := z_n1;
        procedure stepper is new generic_stepper;


        variable run : boolean := true;
        variable loop_count : natural range 0 to 7 := 0;

    begin

        while(run) loop
            loop_count := loop_count + 1;
            rk(y_n1 => y_n1, state => state, h => h, z => z, vErr => vErr, k => k);
            stepper(prev_stepsize => h, vErr => vErr, h_new => h_new, err => err);
            if err < 1.0e-6 or loop_count >= 7 then
                run := false;
            else
                h := h/4.0;
            end if;
        end loop;
        z_n1     := z;
        state    := y_n1;
        stepsize := h_new;

    end generic_adaptive_rk23;

------------------------------------------
    procedure generic_adaptive_dopri54
    generic(
        impure function deriv (input : real_vector) return real_vector is <>
    )
    (
        state    : inout real_vector;
        z_n1     : inout real_vector;
        simtime  : inout real;
        err      : inout real;
        stepsize : inout real
    ) is
        type state_array is array(1 to 7) of real_vector(state'range);
        variable k : state_array;
        variable y_n1 : real_vector(state'range);

        variable tolerance : real := 1.0e-3;
        variable h         : real := stepsize;
        variable h_new     : real ;
        variable vErr      : real_vector(state'range);

        constant dop2 : real_vector := (0 => 1.0/5.0);
        constant dop3 : real_vector := (3.0/40.0       , 9.0/40.0);
        constant dop4 : real_vector := (3.0/40.0       , 9.0/40.0        , 32.0/9.0);
        constant dop5 : real_vector := (19372.0/6561.0 , -25360.0/2187.0 , 64448.0/6561.0 , -212.0/729.0);
        constant dop6 : real_vector := (9017.0/3168.0  , -355.0/33.0     , 46732.0/5247.0 , 49.0/176.0     , -5103.0/18656.0);
        constant dop7 : real_vector := (35.0/384.0     , 0.0             , 500.0/1113.0   , 125.0/192.0    , -2187.0/6784.0    , 11.0/84.0);

        constant dop8 : real_vector := (5179.0/57600.0 , 0.0             , 7571.0/16695.0 , 393.0/640.0    , -92097.0/339200.0 , 187.0/2100.0 , 1.0/40.0);

    begin
        k(1) := z_n1;

        k(2) := deriv(state +
            ( k(1) * dop2(0) 
            ) * stepsize);

        k(3) := deriv(state +
            ( k(1) * dop3(0)
            + k(2) * dop3(1)
            ) * stepsize);

        k(4) := deriv(state +
            ( k(1) * dop4(0)
            + k(2) * dop4(1)
            + k(3) * dop4(2)
            ) * stepsize);

        k(5) := deriv(state +
            ( k(1) * dop5(0)
            + k(2) * dop5(1)
            + k(3) * dop5(2)
            + k(4) * dop5(3)
            ) * stepsize);

        k(6) := deriv(state +
            ( k(1) * dop6(0)
            + k(2) * dop6(1)
            + k(3) * dop6(2)
            + k(4) * dop6(3)
            + k(5) * dop6(4)
            ) * stepsize);

        k(7) := deriv(state +
            ( k(1) * dop7(0)
            + k(2) * dop7(1)
            + k(3) * dop7(2)
            + k(4) * dop7(3)
            + k(5) * dop7(4)
            + k(6) * dop7(5)
            ) * stepsize);

        z_n1 := deriv(state +
            ( k(1) * dop8(0)
            + k(2) * dop8(1)
            + k(3) * dop8(2)
            + k(4) * dop8(3)
            + k(5) * dop8(4)
            + k(6) * dop8(5)
            + k(7) * dop8(6)
            ) * stepsize);

        y_n1  := state + k(7) * stepsize;
        state := y_n1;

        vErr := 
            ( k(1) * (dop8(0) - dop7(0))
            + k(2) * (dop8(1) - dop7(1))
            + k(3) * (dop8(2) - dop7(2))
            + k(4) * (dop8(3) - dop7(3))
            + k(5) * (dop8(4) - dop7(4))
            + k(6) * (dop8(5) - dop7(5))
            + k(7) *  dop8(6)
            ) * stepsize;

        err := norm(vErr); -- use max value

        -- if err < tolerance then
        -- end if;

        -- cubic root
        h_new := h*cbrt(tolerance/err);
        if h_new < default_minstep then
            h_new := default_minstep;
        end if;
        if h_new > 100.0e-6 then
            h_new := 100.0e-6;
        end if;

        simtime  := simtime + stepsize;
        stepsize := h_new;

    end generic_adaptive_dopri54;


end package body adaptive_ode_pkg;
