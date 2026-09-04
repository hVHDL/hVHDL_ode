# Model verification example — LC EMI filter

A worked example of checking one circuit model three ways and putting the
results on a single Bode plot:

1. **Bench measurement** — a swept-sine network-analyser run on the real board.
2. **QSPICE** `.ac` — the schematic simulated in the frequency domain.
3. **hVHDL_ode** — the same circuit as an ODE, integrated with a fixed-step
   RK4 in a VUnit testbench, its transfer function recovered from the
   dithered time series with the Welch H1 estimator.

The device under test is the 3-section LC EMI filter in
[`emi_filter_model.qsch`](emi_filter_model.qsch) (L1 = 1.2 mH, then two
19 µH sections, 680 nF shunt caps, small series R / ESR). The probed output
is the L1/C1 node `V(N04)`, relative to the source.

## Files

| file | what it is |
|------|------------|
| `emi_filter_model.qsch` / `.cir` | QSPICE schematic and netlist, `V1 = 1 V AC`, `.ac` sweep |
| `L1C1_resp.csv` | bench measurement, 500 Hz – 1 MHz (Digilent WaveForms Network Analyzer, ADP3450), `V(N04)` relative to source |
| `L1C1_resp_50k1M.csv` | a second measured sweep, 50 kHz – 1 MHz, for the high-frequency detail |
| `lc_filter_ode_tb.vhd` | the hVHDL_ode RK4 testbench (registered in `vunit_run_ode.py`) |
| `plot_lc_ac_response.py` | runs the QSPICE `.ac` and overlays it on a measurement; can also dump either to a `test_plot.py` saved-response CSV |
| `lc_ac_response.png` | saved QSPICE-vs-measured figure |

The response CSVs the steps below generate (`lc_vhdl_resp.csv`,
`lc_meas_resp.csv`, `lc_qspice_resp.csv`) are build artifacts and are not
checked in.

## Run it

All commands are run from the repository root.

### 1. hVHDL_ode

```bash
python vunit_run_ode.py "*lc_filter_ode*"
python python/test_plot.py lc_filter_ode_tb.dat --freq-only
```

The first line writes `lc_filter_ode_tb.dat`. The second draws the RK4 Bode
estimate and, because of the `#CONFIG freq_save_L1C1=misc/lc_vhdl_resp.csv`
line in the testbench, also dumps that curve to `misc/lc_vhdl_resp.csv`.

The stimulus is sampled and held once per timestep, so the Welch estimate
carries a half-sample (T/2) zero-order-hold group delay the real circuit
does not have. The testbench sets `freq_deskew=zoh` to remove exactly that,
so the phase lines up with the `.ac` and measured sweeps.

### 2. measurement and QSPICE

```bash
python misc/plot_lc_ac_response.py --save-measured misc/lc_meas_resp.csv \
                                   --save-sim misc/lc_qspice_resp.csv --no-plot
```

`--save-measured` just reformats `L1C1_resp.csv` and needs nothing extra.
`--save-sim` runs the QSPICE `.ac` analysis and needs QSPICE (QUX.exe /
QSPICE64.exe) plus the `PyQSPICE` package installed; drop it if you only
want the measurement.

Without any `--save-*` flags the script draws its own QSPICE-vs-measured
figure (`lc_ac_response.png`).

### 3. overlay all three

```bash
python python/test_plot.py lc_filter_ode_tb.dat \
       misc/lc_meas_resp.csv misc/lc_qspice_resp.csv --freq-only
```

`lc_filter_ode_tb.dat` gives the hVHDL_ode curve (recomputed live); the two
CSVs are drawn dashed on top. To plot without re-deriving the hVHDL_ode
response, swap the `.dat` for the saved `misc/lc_vhdl_resp.csv` from step 1.

The three curves land on top of each other through the ~3 kHz main
resonance and the higher-frequency ladder features, which is the point:
the ODE model matches both the reference simulator and the real board.
