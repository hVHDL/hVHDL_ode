LIBRARY ieee  ; 
    USE ieee.NUMERIC_STD.all  ; 
    USE ieee.std_logic_1164.all  ; 
    use ieee.math_real.all;

package adaptive_ode_pkg is
    constant default_minstep : real := 1.0e-9;
    constant default_maxstep : real := 10.0e-3;
    constant default_tolerance : real := 1.0e-5;

    -- array of state vectors: stage values, dense-output coefficients, ...
    type st_array is array(natural range <>) of real_vector;

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
    -- Dormand-Prince 5(4) ("DOPRI5", the ode45 pair) - one adaptive step with
    -- a continuous extension. Cleaner than generic_adaptive_dopri54: proper
    -- FSAL (k7 is carried into the next step as k1), the standard embedded
    -- 4th-order error estimate, and `tolerance` exposed as a generic.
    --
    -- Usage: seed k_fsal with deriv(t, state) once before the first call, then
    -- reuse it. After the call, state has advanced by h_accepted (from t to
    -- t + h_accepted), stepsize holds the recommended next step, and `cont`
    -- holds the dense-output coefficients for that step - pass them to
    -- dopri5_dense to evaluate the solution anywhere in [t, t + h_accepted].
    procedure generic_adaptive_dopri5
    generic(
        impure function deriv (t : real; input : real_vector) return real_vector is <>
        ; minstep   : real := default_minstep
        ; maxstep   : real := default_maxstep
        ; tolerance : real := default_tolerance
    )
    (
        t          : in    real;          -- time at the start of the step
        state      : inout real_vector;   -- y(t) in, y(t + h_accepted) out
        k_fsal     : inout real_vector;   -- f(t, y(t)) in, f(t+h, y_new) out
        err        : inout real;          -- error norm of the accepted step
        stepsize   : inout real;          -- step to attempt in, next step out
        h_accepted : out   real;          -- step length actually taken
        cont       : out   st_array       -- (1 to 5)(state'range) dense-output coeffs
    );
------------------------------------------
    -- Evaluate the DOPRI5 continuous extension produced by
    -- generic_adaptive_dopri5. theta = (t_query - t_step_start) / h_accepted,
    -- in [0, 1]. Fourth-order accurate, matches the step ends exactly.
    function dopri5_dense(theta : real; cont : st_array) return real_vector;
------------------------------------------
    -- Adaptive DOPRI5 driver that logs on an arbitrary (uniform) output grid
    -- rather than at the adaptive step points: it steps internally with error
    -- control and calls log_sample(t, state) at t_start, t_start + log_step,
    -- ... up to t_end, interpolating each sample with dopri5_dense. Pass
    -- log_step <= 0.0 to log every accepted step instead. On return `state`
    -- holds y(t_end).
    procedure generic_dopri5_uniform_log
    generic(
        impure function deriv (t : real; input : real_vector) return real_vector is <>
        ; procedure log_sample (t : real; state : real_vector)
        ; minstep   : real := default_minstep
        ; maxstep   : real := default_maxstep
        ; tolerance : real := default_tolerance
    )
    (
        t_start  : in    real;
        t_end    : in    real;
        log_step : in    real;
        state    : inout real_vector;
        stepsize : inout real
    );
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
            h_new := 0.9*prev_stepsize*error_norm(default_tolerance/err); -- cbrt() is cubic root
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

            z := k(4);

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
            -- k(1) := z;
            k(1) := deriv(t, state);

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

            z := k(7);

            -- vErr := y_n1 - z;
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

        while(run) loop
            loop_count := loop_count + 1;
            rk(t => t, y_n1 => y_n1, state => state, h => h, z => z, vErr => vErr, k => k);
            stepper(prev_stepsize => h, vErr => vErr, h_new => h_new, err => err);
            if err < 1.0e-4 or loop_count >= 7 then
                run := false;
            else
                h := h/4.0;
                run := true;
            end if;
        end loop;
        z_n1     := z;
        state    := y_n1;
        stepsize := h_new;
        -- stepsize := h;

    end generic_adaptive_dopri54;

------------------------------------------
    procedure generic_adaptive_dopri5
    generic(
        impure function deriv (t : real; input : real_vector) return real_vector is <>
        ; minstep   : real := default_minstep
        ; maxstep   : real := default_maxstep
        ; tolerance : real := default_tolerance
    )
    (
        t          : in    real;
        state      : inout real_vector;
        k_fsal     : inout real_vector;
        err        : inout real;
        stepsize   : inout real;
        h_accepted : out   real;
        cont       : out   st_array
    ) is
        -- nodes
        constant c2 : real := 1.0/5.0;
        constant c3 : real := 3.0/10.0;
        constant c4 : real := 4.0/5.0;
        constant c5 : real := 8.0/9.0;
        -- Runge-Kutta matrix
        constant a21 : real := 1.0/5.0;
        constant a31 : real := 3.0/40.0;        constant a32 : real := 9.0/40.0;
        constant a41 : real := 44.0/45.0;       constant a42 : real := -56.0/15.0;       constant a43 : real := 32.0/9.0;
        constant a51 : real := 19372.0/6561.0;  constant a52 : real := -25360.0/2187.0;  constant a53 : real := 64448.0/6561.0;  constant a54 : real := -212.0/729.0;
        constant a61 : real := 9017.0/3168.0;   constant a62 : real := -355.0/33.0;      constant a63 : real := 46732.0/5247.0;  constant a64 : real := 49.0/176.0;      constant a65 : real := -5103.0/18656.0;
        -- 5th-order weights (== a7 row, FSAL); b(1) = b(6) = 0 so k2, k7 drop out
        constant b1 : real := 35.0/384.0;
        constant b3 : real := 500.0/1113.0;
        constant b4 : real := 125.0/192.0;
        constant b5 : real := -2187.0/6784.0;
        constant b6 : real := 11.0/84.0;
        -- ei = b_i - bhat_i (5th-order minus embedded 4th-order weights).
        -- e7 is written as -(x) rather than 0.0 - x: nvc <= 1.18.2 folds the
        -- static expression "0.0 - x" to "+x" (constant-folding bug, fixed in
        -- nvc 1.22).
        constant e1 : real := 35.0/384.0     - 5179.0/57600.0;
        constant e3 : real := 500.0/1113.0   - 7571.0/16695.0;
        constant e4 : real := 125.0/192.0    - 393.0/640.0;
        constant e5 : real := -2187.0/6784.0 - (-92097.0/339200.0);
        constant e6 : real := 11.0/84.0      - 187.0/2100.0;
        constant e7 : real := -(1.0/40.0);
        -- dense-output (continuous extension) coefficients, Hairer CONTD5
        constant d1 : real := -12715105075.0 / 11282082432.0;
        constant d3 : real :=  87487479700.0 / 32700410799.0;
        constant d4 : real := -10690763975.0 / 1880347072.0;
        constant d5 : real :=  701980252875.0 / 199316789632.0;
        constant d6 : real := -1453857185.0 / 822651844.0;
        constant d7 : real :=  69997945.0 / 29380423.0;

        subtype vec is real_vector(state'range);
        subtype stage_array is st_array(1 to 7)(state'range);

        variable k    : stage_array;
        variable y0   : vec := state;
        variable y1   : vec;
        variable evec : vec;
        variable h    : real := stepsize;
        variable e    : real := 0.0;
        variable fac  : real := 1.0;
        variable tries : natural := 0;
    begin
        loop
            tries := tries + 1;

            k(1) := k_fsal;                                        -- FSAL: f(t, y0)
            k(2) := deriv(t + c2*h, y0 + (k(1)*a21) * h);
            k(3) := deriv(t + c3*h, y0 + (k(1)*a31 + k(2)*a32) * h);
            k(4) := deriv(t + c4*h, y0 + (k(1)*a41 + k(2)*a42 + k(3)*a43) * h);
            k(5) := deriv(t + c5*h, y0 + (k(1)*a51 + k(2)*a52 + k(3)*a53 + k(4)*a54) * h);
            k(6) := deriv(t + h,    y0 + (k(1)*a61 + k(2)*a62 + k(3)*a63 + k(4)*a64 + k(5)*a65) * h);

            y1 := y0 + (k(1)*b1 + k(3)*b3 + k(4)*b4 + k(5)*b5 + k(6)*b6) * h;

            k(7) := deriv(t + h, y1);                              -- FSAL stage

            evec := (k(1)*e1 + k(3)*e3 + k(4)*e4 + k(5)*e5 + k(6)*e6 + k(7)*e7) * h;
            e    := norm(evec);

            if e > 1.0e-15 then
                fac := 0.9 * exp(0.2 * log(tolerance / e));
                if fac < 0.2 then fac := 0.2; end if;
                if fac > 5.0 then fac := 5.0; end if;
            else
                fac := 5.0;
            end if;

            exit when e <= tolerance or h <= minstep or tries >= 12;

            h := h * fac;                                          -- rejected: shrink and retry
            if h < minstep then h := minstep; end if;
        end loop;

        -- continuous extension over [t, t + h]
        cont(1) := y0;
        cont(2) := y1 - y0;
        cont(3) := k(1)*h - cont(2);
        cont(4) := cont(2) - k(7)*h - cont(3);
        cont(5) := (k(1)*d1 + k(3)*d3 + k(4)*d4 + k(5)*d5 + k(6)*d6 + k(7)*d7) * h;

        state      := y1;
        k_fsal     := k(7);                                        -- -> k(1) of the next step
        err        := e;
        h_accepted := h;

        h := h * fac;                                              -- recommend the next step
        if h < minstep then h := minstep; end if;
        if h > maxstep then h := maxstep; end if;
        stepsize := h;
    end generic_adaptive_dopri5;

------------------------------------------
    function dopri5_dense(theta : real; cont : st_array) return real_vector is
        constant lo : natural := cont'low;
        constant s  : real := theta;
        constant s1 : real := 1.0 - theta;
        variable a  : real_vector(cont(lo)'range);
        variable b  : real_vector(cont(lo)'range);
        variable cc : real_vector(cont(lo)'range);
    begin
        a  := cont(lo + 3) + cont(lo + 4) * s1;   -- cont4 + s1*cont5
        b  := cont(lo + 2) + a * s;               -- cont3 + s *(...)
        cc := cont(lo + 1) + b * s1;              -- cont2 + s1*(...)
        return cont(lo) + cc * s;                 -- cont1 + s *(...)
    end function dopri5_dense;

------------------------------------------
    procedure generic_dopri5_uniform_log
    generic(
        impure function deriv (t : real; input : real_vector) return real_vector is <>
        ; procedure log_sample (t : real; state : real_vector)
        ; minstep   : real := default_minstep
        ; maxstep   : real := default_maxstep
        ; tolerance : real := default_tolerance
    )
    (
        t_start  : in    real;
        t_end    : in    real;
        log_step : in    real;
        state    : inout real_vector;
        stepsize : inout real
    ) is
        -- forward the generic deriv through a concrete local function so it can
        -- be mapped onto the nested generic_adaptive_dopri5 instantiation
        impure function deriv_fwd (t : real; input : real_vector) return real_vector is
        begin
            return deriv(t, input);
        end function;

        procedure dopri5_step is new generic_adaptive_dopri5
            generic map(deriv => deriv_fwd, minstep => minstep, maxstep => maxstep, tolerance => tolerance);

        subtype vec is real_vector(state'range);
        variable y     : vec := state;
        variable kf    : vec;
        variable cont  : st_array(1 to 5)(state'range);
        variable t     : real := t_start;
        variable h_acc : real := 0.0;
        variable e     : real := 0.0;
        variable t_log : real := t_start;
        variable theta : real;
        variable saved_step : real;
        variable clamped    : boolean;
        constant eps   : real := 1.0e-12;
    begin
        kf := deriv(t_start, y);

        if log_step > 0.0 then
            log_sample(t_start, y);
            t_log := t_start + log_step;
        end if;

        while t < t_end - eps loop
            clamped := false;
            if t + stepsize > t_end then
                saved_step := stepsize;
                stepsize   := t_end - t;
                clamped    := true;
            end if;

            dopri5_step(t => t, state => y, k_fsal => kf, err => e,
                        stepsize => stepsize, h_accepted => h_acc, cont => cont);

            if clamped then
                stepsize := saved_step;      -- interval done; keep the natural guess
            end if;

            if log_step > 0.0 then
                while t_log <= t + h_acc + eps and t_log <= t_end + eps loop
                    theta := (t_log - t) / h_acc;
                    if theta < 0.0 then theta := 0.0; end if;
                    if theta > 1.0 then theta := 1.0; end if;
                    log_sample(t_log, dopri5_dense(theta, cont));
                    t_log := t_log + log_step;
                end loop;
            else
                log_sample(t + h_acc, y);
            end if;

            t := t + h_acc;
        end loop;

        state := y;
    end generic_dopri5_uniform_log;


end package body adaptive_ode_pkg;
