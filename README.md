# hVHDL_ode

Numerical ODE integrators (Runge-Kutta and Adams-Moulton) as reusable VHDL
packages, plus switching-model testbenches for power-electronic converters
that use them for the plant simulation.

Developed and tested with the open-source [NVC](https://www.nickg.me.uk/nvc/)
simulator, driven by [VUnit](https://vunit.github.io/).

The `testbenches/converter/` switching models, the reorganization, and the
`test_plot.py` plotting helpers were written with [Claude Code](https://claude.com/claude-code)
(Claude Sonnet 5); see the co-authored commits.

## Layout

```
ode_solvers/       integrator packages
  ode_pkg.vhd            fixed-step  generic_rk1/rk2/rk4/rk5, am2/am4 (generic over `deriv`)
  adaptive_ode_pkg.vhd   adaptive    generic_adaptive_rk23 / _dopri54
  real_vector_pkg.vhd    real_vector arithmetic used by the solvers
  sort_pkg.vhd
write_pkg.vhd      write_to / init_simfile / write_plot_config -> whitespace .dat + "#CONFIG" lines
python/            test_plot.py (time-domain + Bode from a .dat), freq_response.py, qspice helpers
qspice_ref_models/ QSPICE schematics used as reference for some models
ode_examples/      grid_inverter_control_tb  (worked control-loop example)
testbenches/
  lcr_*            LC / LCR filter integration tests (fixed vs adaptive, single- and 3-phase)
  sort_tb, template_tb
  converter/       power-electronic converter switching models
    multilevel/    fc_3level_tb, fc_4level_tb, fc_5level_tb        (flying-capacitor buck-type)
    dcdc/          buck_converter_tb, boost_converter_tb          (synchronous 2-switch)
                   dab_converter_tb, dhb_converter_tb             (dual active [half-]bridge, SPS)
                   llc_converter_tb                               (half-bridge LLC resonant, 400 -> 51 V)
    multiphase/    buck_3ph_tb, boost_3ph_tb                      (3-phase interleaved)
    dcac/          inverter_3ph_tb, inverter_3ph_svm_tb           (3-phase VSI + LC filter)
    acdc/          boost_pfc_tb                                   (single-phase boost PFC)
```

## Solvers

Each integrator is a `procedure` generic over an
`impure function deriv(t : real; state : real_vector) return real_vector`.
Instantiate it for a model's derivative and call it once per step:

```vhdl
procedure rk5 is new generic_rk5 generic map(deriv_lcr);
...
rk5(t, state, stepsize);          -- advances `state` in place
```

`generic_rk5` is the Dormand-Prince 5th-order pair; `am2_generic` /
`am4_generic` are multi-step and take a history array. The adaptive
variants in `adaptive_ode_pkg` size the step from an error estimate.

## Converter testbenches

All of them model the power stage as an ODE integrated with a fixed-step
`rk5` and represent the switching by stepping through the converter's
switch states, using a step length equal to that sub-interval's on/off
time (no averaging, so ripple and switch-node waveforms are visible).

Shared conventions in the `converter/` testbenches:

- `hb_modulator(sw_state)` returns `1.0` / `0.0` for a half-bridge leg;
  switched quantities are formed as `hb_modulator(sw) * rail`.
- Time is an integer tick count; `minimum_time_step` sets the seconds per
  tick, and conversions happen only where a real time is needed (rk5 step
  size, the plot file).
- A `stimulus` / `p_modulation` process pair with a
  `sim_ready` / `modulation_ready` handshake: `stimulus` integrates one
  sub-interval, then asks `p_modulation` for the next duty / modulator
  solution.
- `multiphase/` and `dcac/` use a grid-locked edge scheduler so the phase
  relationship does not drift over long runs.
- `inverter_3ph_svm_tb` builds the switching pattern from a switch-state
  matrix in the flying-capacitor style; `scheme` selects 7-segment SVPWM
  or the 5-segment minimal-commutation DPWMMIN matrix.
- `dab_converter_tb` / `dhb_converter_tb` also use a switch-state matrix:
  single-phase-shift modulation of the two (half-)bridges, with a PI loop
  on the phase-shift fraction regulating the secondary DC-link voltage.
  The transformer is a T model - the three branch currents are states and
  the midpoint voltage is solved algebraically like the 3-phase neutral.
- `boost_pfc_tb` is a single-phase boost PFC: diode-bridge rectifier +
  boost stage, with a slow outer voltage PI setting the amplitude of a
  rectified-sine current reference that an inner current loop tracks
  (unity power factor, ~2 % THD in the model). Discontinuous conduction
  near the line zero crossings is resolved with a bisection search for
  the inductor-current zero crossing, not just clamped.
- `llc_converter_tb` is a half-bridge LLC resonant converter regulated by
  frequency (a PI moves f_sw around the series resonance). The Lr-Cr-Lm
  tank is a 3-state ODE; the rectifier freewheel interval below resonance
  is found with the same bisection search as the PFC's DCM.

## Running

```bash
python vunit_run_ode.py            # run everything
python vunit_run_ode.py "*buck*"   # filter by test name
python vunit_run_ode.py -l         # list tests
```

Each testbench writes `<entity>.dat` (whitespace-separated columns with a
`time` column, preceded by `#CONFIG key=value` lines).

## Plotting

```bash
python python/test_plot.py buck_converter_tb.dat
python python/test_plot.py --step buck_converter_tb.dat        # stepped switch-node lines
python python/test_plot.py --freq-only template_tb.dat         # Bode only
```

Column-name prefixes route data to sub-plots: `T_*` to the top axis, `B_*`
to the bottom. Testbenches steer titles, labels, limits, drawstyle and
frequency-response pairs by writing `#CONFIG` lines (see
`write_plot_config` and `read_plot_config` in `python/test_plot.py`);
`--drawstyle` / `--step` and `--freq-pair` override them from the command
line.
