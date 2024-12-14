LIBRARY ieee  ; 
    USE ieee.NUMERIC_STD.all  ; 
    USE ieee.std_logic_1164.all  ; 
    use ieee.math_real.all;

package adaptive_ode_pkg is 
    constant default_minstep : real := 1.0e-9;
    constant default_maxstep : real := 10.0e-3;
    constant default_tolerance : real := 1.0e-7;

------------------------------------------
    procedure generic_adaptive_rk23
    generic(
        impure function deriv (t : real; input : real_vector) return real_vector is <>
        ;minstep : real := default_minstep
        ;maxstep : real := default_maxstep
    )
    (
        t : real;
        state    : inout real_vector;
        z_n1     : inout real_vector;
        err      : inout real;
        stepsize : inout real);

------------------------------------------
    procedure generic_adaptive_dopri54
    generic(
        impure function deriv (t : real; input : real_vector) return real_vector is <>
        ;minstep : real := default_minstep
        ;maxstep : real := default_maxstep
    )
    (
        t        : in real;
        state    : inout real_vector;
        z_n1     : inout real_vector;
        err      : inout real;
        stepsize : inout real);
------------------------------------------

end package adaptive_ode_pkg;

package body adaptive_ode_pkg is 

    use work.real_vector_pkg.all;
    type st_array is array(natural range <>) of real_vector;
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
        function error_norm(X : real) return real is <>
        ;minstep : real := default_minstep
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
            h_new := prev_stepsize*error_norm(default_tolerance/err); -- cbrt() is cubic root
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
    --------
    procedure generic_adaptive_rk23
    generic(
        impure function deriv (t : real; input : real_vector) return real_vector is <>
        ;minstep : real := default_minstep
        ;maxstep : real := default_maxstep
    )
    (
        t : real;
        state    : inout real_vector;
        z_n1     : inout real_vector;
        err      : inout real;
        stepsize : inout real
    ) is

        procedure rk
        -- generic(
        --     impure function deriv(t : real; input : real_vector) return real_vector is <>
        -- )
        (t : real; y_n1 : inout real_vector; state : real_vector; h: real; z : inout real_vector; vErr : inout real_vector; k : inout st_array) 
        is
            constant dop2 : real_vector := (0 => 1.0/2.0);
            constant dop3 : real_vector := (0.0       , 3.0/4.0);
            constant dop4 : real_vector := (2.0/9.0      , 1.0/3.0      , 4.0/9.0);
            constant dop5 : real_vector := (7.0/24.0, 1.0/4.0, 1.0/3.0, 1.0/8.0);
        begin
            k(1) := z;
            k(2) := deriv(t + h/2.0, state + k(1) * h * dop2(0));
            k(3) := deriv(t + h*3.0/4.0, state + k(2) * h * dop3(1));

            y_n1 := state + (k(1)*dop4(0) + k(2)*dop4(1) + k(3)*dop4(2)) * h;

            k(4) := deriv(t + h, y_n1);

            z := state + k(4) * h;

            vErr := (
                k(1)*(dop4(0) - dop5(0)) 
                +k(2)*(dop4(1) - dop5(1)) 
                +k(3)*(dop4(2) - dop5(2)) 
                +k(4)*( - dop5(3)) 
            ) * h;
        end rk;
        --------

        -- procedure rk is new rk_generic generic map(deriv);
        procedure stepper is new generic_stepper generic map (cbrt, minstep => minstep, maxstep => maxstep);

        subtype state_array is st_array(1 to 4)(state'range);
        variable k : state_array;
        variable y_n1 : real_vector(state'range);

        variable h     : real := stepsize;
        variable h_new : real ;
        variable vErr  : real_vector(state'range);
        variable z     : real_vector(z_n1'range) := z_n1;

        variable run : boolean := true;
        variable loop_count : natural range 0 to 7 := 0;

    begin

        while(run) loop
            loop_count := loop_count + 1;
            rk(t => t, y_n1 => y_n1, state => state, h => h, z => z, vErr => vErr, k => k);
            stepper(prev_stepsize => h, vErr => vErr, h_new => h_new, err => err);
            if err < 1.0e-4 or loop_count >= 7 then
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
        impure function deriv (t : real; input : real_vector) return real_vector is <>
        ;minstep : real := default_minstep
        ;maxstep : real := default_maxstep
    )
    (
        t        : in real;
        state    : inout real_vector ;
        z_n1     : inout real_vector ;
        err      : inout real        ;
        stepsize : inout real
    ) is

        procedure rk
            (t : real; y_n1 : inout real_vector; state : real_vector; h: real; z : inout real_vector; vErr : inout real_vector; k : inout st_array) 
        is
            constant dop2 : real_vector := (0 => 1.0/5.0);
            constant dop3 : real_vector := (3.0/40.0       , 9.0/40.0);
            constant dop4 : real_vector := (44.0/45.0      , -56.0/15.0      , 32.0/9.0);
            constant dop5 : real_vector := (19372.0/6561.0 , -25360.0/2187.0 , 64448.0/6561.0 , -212.0/729.0);
            constant dop6 : real_vector := (9017.0/3168.0  , -355.0/33.0     , 46732.0/5247.0 , 49.0/176.0     , -5103.0/18656.0);
            constant dop7 : real_vector := (35.0/384.0     , 0.0             , 500.0/1113.0   , 125.0/192.0    , -2187.0/6784.0    , 11.0/84.0);

            constant dop8 : real_vector := (5179.0/57600.0 , 0.0     , 7571.0/16695.0 , 393.0/640.0 , -92097.0/339200.0 , 187.0/2100.0 , 1.0/40.0);

            constant tdop : real_vector := (0.0 , 1.0/5.0 , 3.0/10.0 , 4.0/5.0 , 8.0/9.0 , 1.0 , 1.0);
        begin
            k(1) := z;

            k(2) := deriv(t + h*tdop(1), state +
                ( k(1) * dop2(0) 
                ) * h);

            k(3) := deriv(t + h*tdop(2), state +
                ( k(1) * dop3(0)
                + k(2) * dop3(1)
                ) * h);

            k(4) := deriv(t + h*tdop(3), state +
                ( k(1) * dop4(0)
                + k(2) * dop4(1)
                + k(3) * dop4(2)
                ) * h);

            k(5) := deriv(t + h*tdop(4), state +
                ( k(1) * dop5(0)
                + k(2) * dop5(1)
                + k(3) * dop5(2)
                + k(4) * dop5(3)
                ) * h);

            k(6) := deriv(t + h*tdop(5), state +
                ( k(1) * dop6(0)
                + k(2) * dop6(1)
                + k(3) * dop6(2)
                + k(4) * dop6(3)
                + k(5) * dop6(4)
                ) * h);

            y_n1 := 
                state +
                ( k(1) * dop7(0)
                + k(2) * dop7(1)
                + k(3) * dop7(2)
                + k(4) * dop7(3)
                + k(5) * dop7(4)
                + k(6) * dop7(5)
                ) * h;

            k(7) := deriv(t + h, y_n1);

            z := state +
                ( k(1) * dop8(0)
                + k(2) * dop8(1)
                + k(3) * dop8(2)
                + k(4) * dop8(3)
                + k(5) * dop8(4)
                + k(6) * dop8(5)
                + k(7) * dop8(6)
                ) * h;


            vErr := 
                ( k(1) * (dop8(0) - dop7(0))
                + k(2) * (dop8(1) - dop7(1))
                + k(3) * (dop8(2) - dop7(2))
                + k(4) * (dop8(3) - dop7(3))
                + k(5) * (dop8(4) - dop7(4))
                + k(6) * (dop8(5) - dop7(5))
                + k(7) *  dop8(6)
                ) * h;
            end rk;
        --------
        subtype state_array is st_array(1 to 7)(state'range);
        variable k : state_array;
        variable y_n1 : real_vector(state'range);

        function fifth_root(X : real) return real is
        begin
            return X**(1.0/5.0);
        end fifth_root;

        procedure stepper is new generic_stepper generic map(fifth_root, minstep, maxstep);

        variable h     : real := stepsize;
        variable h_new : real ;
        variable vErr  : real_vector(state'range);
        variable z     : real_vector(z_n1'range) := z_n1;

        variable run : boolean := true;
        variable loop_count : natural range 0 to 7 := 0;


    begin

        -- while(run) loop
        --     loop_count := loop_count + 1;
            rk(t => t, y_n1 => y_n1, state => state, h => h, z => z, vErr => vErr, k => k);
            -- stepper(prev_stepsize => h, vErr => vErr, h_new => h_new, err => err);
            -- if err < 1.0e-6 or loop_count >= 7 then
            --     run := false;
            -- else
            --     h := h/4.0;
            --     run := true;
            -- end if;
        -- end loop;
        z_n1     := y_n1;
        state    := y_n1;
        -- stepsize := h_new;
        stepsize := h;

    end generic_adaptive_dopri54;


end package body adaptive_ode_pkg;
