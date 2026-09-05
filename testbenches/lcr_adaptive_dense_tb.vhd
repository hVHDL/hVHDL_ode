----------------------------------
-- Adaptive DOPRI5 with dense output.
--
-- The 3-phase LCR filter of lcr_3ph_adaptive_tb, but integrated with
-- generic_dopri5_uniform_log: the solver takes error-controlled adaptive
-- steps internally, while the .dat file is written on a fixed uniform
-- grid (log_step) using the DOPRI5 continuous extension. This decouples
-- the output sampling from the solver's step points.
--
-- state vector : (0..2 => inductor currents, 3..5 => capacitor voltages)
----------------------------------
LIBRARY ieee  ;
    USE ieee.NUMERIC_STD.all  ;
    USE ieee.std_logic_1164.all  ;
    use ieee.math_real.all;
    use std.textio.all;

library vunit_lib;
context vunit_lib.vunit_context;

    use work.write_pkg.all;
    use work.adaptive_ode_pkg.all;
    use work.lcr_models_pkg.all;

entity lcr_adaptive_dense_tb is
  generic (runner_cfg : string);
end;

architecture vunit_simulation of lcr_adaptive_dense_tb is

    constant stoptime : real := 10.0e-3;
    constant log_step : real := 20.0e-6;   -- uniform output grid

begin

------------------------------------------------------------------------
    stimulus : process

        constant i_load : real_vector (0 to 1) := (others => 0.0);
        constant l      : real_vector (1 to 3) := (others => 80.0e-6);
        constant c      : real_vector (1 to 3) := (1 => 200.0e-6, others => 60.0e-6);
        constant r      : real_vector (1 to 3) := (others => 100.0e-3);

        impure function deriv(t : real; states : real_vector) return real_vector is
            variable uin : real_vector(1 to 3);
        begin
            -- 3-phase 1 kHz excitation, close to the filter resonance so the
            -- transient is long enough to make the adaptive step vary
            uin := ( sin( t*1000.0*math_pi*2.0)
                   , sin((t*1000.0 + 1.0/3.0)*math_pi*2.0)
                   , sin((t*1000.0 + 2.0/3.0)*math_pi*2.0) );

            return deriv_lcr(states, i_load, uin, l, c, r);
        end deriv;

        variable state    : real_vector(0 to 5) := (others => 0.0);
        variable stepsize : real := 50.0e-6;

        file file_handler : text open write_mode is "lcr_adaptive_dense_tb.dat";

        procedure log_sample(t : real; s : real_vector) is
        begin
            write_to(file_handler, ( t, s(3), s(4), s(5), s(0), s(1), s(2) ));
        end procedure;

        procedure run is new generic_dopri5_uniform_log
            generic map(deriv => deriv, log_sample => log_sample, maxstep => 200.0e-6);

    begin
        test_runner_setup(runner, runner_cfg);

        init_simfile(file_handler, ("time"
        ,"T_u0"
        ,"T_u1"
        ,"T_u2"
        ,"B_i0"
        ,"B_i1"
        ,"B_i2"
        ));

        run(t_start  => 0.0,
            t_end    => stoptime,
            log_step => log_step,
            state    => state,
            stepsize => stepsize);

        test_runner_cleanup(runner);
        wait;
    end process stimulus;
------------------------------------------------------------------------
end vunit_simulation;
