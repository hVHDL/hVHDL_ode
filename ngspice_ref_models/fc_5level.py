#!/usr/bin/env python3
"""
ngspice reference model of the 5-level flying-capacitor converter in
testbenches/converter/multilevel/fc_5level_tb.vhd.

The VHDL testbench integrates the plant equations directly; this builds the
same converter as an actual switch circuit, so the two can be checked
against each other:

    Udc ─ S3 ─ p3 ─ S2 ─ p2 ─ S1 ─ p1 ─ S0 ─┐
               │          │          │       sw ─ L ─ RL ─ out ─┬─ C
              Cfc2       Cfc1       Cfc0     │                  └─ i_load
               │          │          │       │
    gnd ─ S3'─ n3 ─ S2'─ n2 ─ S1'─ n1 ─ S0'─┘

Cell k is the complementary pair Sk / Sk', gated by bit k of the switch
pattern; flying cap k sits across cell k+1's inputs (p_{k+1} - n_{k+1}),
which gives the same bridge voltage and cap currents as the VHDL
get_fc_bridge_voltage / fc_modulator equations.

The switch-pattern sequence is read from fc_sw_matrix in
fc_5level_modulator_pkg.vhd, and the dwell times are computed the same way
as fc_modulator_common_pkg.get_sequence_times, so the gate waveforms match
what the VHDL modulator emits (open loop, constant duty, no trim).

usage:
    python3 ngspice_ref_models/fc_5level.py               # write netlist, run ngspice
    python3 ngspice_ref_models/fc_5level.py --no-run      # netlist only
    python3 ngspice_ref_models/fc_5level.py --compare fc_5level_tb.dat
    python3 ngspice_ref_models/fc_5level.py --no-run --compare fc_5level_tb.dat   # re-plot
    python3 ngspice_ref_models/fc_5level.py --simulator qspice --compare fc_5level_tb.dat

--simulator qspice runs the same circuit in QSPICE80.exe (Windows only ; under
WSL the install in C:\Program Files\QSPICE is called directly), and reads
the .qraw back.

Results go to ngspice_ref_models/out/ (--outdir). With the default 20 ns
max step the full 400 ms takes ~50 s in ngspice 36 and ~140 s in QSPICE.
The open-loop flying caps are barely damped, so a coarser --tmax drifts
visibly, and so do slow gate edges : both simulators' switches overlap top
and bottom over the edge and bleed the caps (1 ns edges : ~2.7 V cap error
by 400 ms in ngspice, against 0.4 V at the default 0.1 ns).
"""

import argparse
import math
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PKG = ROOT / "testbenches/converter/multilevel/fc_5level_modulator_pkg.vhd"

# plant and modulator parameters, as in fc_5level_tb.vhd
PARAMS = dict(
    udc=200.0,
    vref=149.0,
    f_sw=500.0e3,
    l=20.0e-6,
    c=10.0e-6,
    rl=100.0e-3,
    cfc=2.0e-6,
    i_load=10.1111,
)


def read_sw_matrix(pkg=PKG):
    """fc_sw_matrix rows as lists of bit strings ("3210" order)."""
    text = pkg.read_text()
    body = text[text.index("fc_sw_matrix"):]
    body = body[:body.index(");\n")]
    rows = {}
    for m in re.finditer(r"(\d+)\s*=>\s*\(([^)]*)\)", body):
        rows[int(m.group(1))] = re.findall(r'"([01]+)"', m.group(2))
    return [rows[k] for k in sorted(rows)]


def to_level(duty, n_cells):
    d = min(max(duty, 0.0), 1.0)
    return min(int(math.floor(d * n_cells)), n_cells - 1)


def sequence_times(duty, n_cells, t_sw, n_states):
    """fc_modulator_common_pkg.get_sequence_times with zero trim."""
    d = min(max(duty, 0.0), 1.0)
    within = min(max(d * n_cells - to_level(d, n_cells), 0.0), 1.0)
    times = []
    for s in range(n_states):
        t = t_sw * within if s % 2 == 0 else t_sw * (1.0 - within)
        times.append(max(t, 0.01 * t_sw))
    return times


def gate_pulse(seq, times, cell, t_rise, inverted=False):
    """One cell's gate over the switch sequence as a PULSE source.

    Every cell turns on and off once per sequence, so one pulse per period
    covers it. The pulse starts from the state-0 value so t < tdelay is right.
    (ngspice's repeating PWL, "pwl(...) r=0", was tried first : its run time
    grows with the square of tstop, 20 ms took 7 min against 2.6 s here.)
    """
    vals = [int(st[-1 - cell]) for st in seq]
    if inverted:
        vals = [1 - v for v in vals]
    if len(set(vals)) == 1:
        return f"{vals[0]}"
    starts = [sum(times[:s]) for s in range(len(seq))]
    period = sum(times)
    v0 = vals[0]
    # first edge away from the state-0 value, then the edge back to it
    away = next(s for s in range(1, len(seq)) if vals[s] != v0)
    back = next((s for s in range(away + 1, len(seq)) if vals[s] == v0), len(seq))
    t_away = starts[away]
    t_back = starts[back] if back < len(seq) else period
    width = t_back - t_away - t_rise
    return (f"PULSE({v0} {1 - v0} {t_away:.9g} {t_rise} {t_rise} {width:.9g} {period:.9g})")


def netlist(p, tstop, tmax, tstep, t_rise=1e-10, simulator="ngspice"):
    matrix = read_sw_matrix()
    n_cells = len(matrix)
    n_states = len(matrix[0])
    t_sw = 1.0 / p["f_sw"]
    duty = p["vref"] / p["udc"]
    level = to_level(duty, n_cells)
    seq = matrix[level]
    times = sequence_times(duty, n_cells, t_sw, n_states)

    lines = [
        f"* {n_cells + 1}-level flying-capacitor converter, switching model",
        f"* generated by {Path(__file__).name} from {PKG.name}",
        f"* duty {duty:.6g} -> level {level}, sequence {' '.join(seq)}",
        f"* dwell times {' '.join(f'{t:.4g}' for t in times)}",
        "",
        f"Vdc p{n_cells} 0 {p['udc']}",
        "",
        "* gates, top switch of cell k on when g<k> is high, bottom on when gb<k> is high",
    ]
    for k in range(n_cells):
        lines.append(f"Vg{k}  g{k}  0 {gate_pulse(seq, times, k, t_rise)}")
        lines.append(f"Vgb{k} gb{k} 0 {gate_pulse(seq, times, k, t_rise, inverted=True)}")
    lines += ["", "* switch cells, p0 = n0 = sw"]
    for k in range(n_cells):
        hi_p = f"p{k + 1}"
        hi_n = f"n{k + 1}" if k + 1 < n_cells else "0"
        lo = "sw" if k == 0 else f"p{k}"
        lo_n = "sw" if k == 0 else f"n{k}"
        lines.append(f"S{k}t {hi_p} {lo}   g{k}  0 swmod")
        lines.append(f"S{k}b {lo_n} {hi_n}   gb{k} 0 swmod")
    lines += ["", "* flying caps, balanced start at k/n * Udc"]
    for k in range(n_cells - 1):
        u0 = p["udc"] * (k + 1) / n_cells
        lines.append(f"Cfc{k} p{k + 1} n{k + 1} {p['cfc']} ic={u0}")
    lines += [
        "",
        "* output filter and constant-current load",
        f"L1 sw nl {p['l']} ic=0",
        f"RL nl out {p['rl']}",
        f"C1 out 0 {p['c']} ic={p['vref']}",
        f"Iload out 0 {p['i_load']}",
        "",
    ]
    if simulator == "qspice":
        lines += [
            ".model swmod SW(Ron=1m Roff=10Meg Vt=0.5 Vh=0)",
            f".tran 0 {tstop} 0 {tmax} uic",
            ".save V(p1) V(n1) V(p2) V(n2) V(p3) V(n3) V(out) V(sw) I(L1)",
            ".end",
        ]
    else:
        lines += [
            ".model swmod sw(vt=0.5 vh=0 ron=1m roff=10Meg)",
            ".options method=gear reltol=1e-4 abstol=1e-9 vntol=1e-6 interp",
            "",
            ".control",
            f"  tran {tstep} {tstop} 0 {tmax} uic",
            "  let ufc0 = v(p1) - v(n1)",
            "  let ufc1 = v(p2) - v(n2)",
            "  let ufc2 = v(p3) - v(n3)",
            "  let il = i(L1)",
            "  let uout = v(out)",
            "  let usw = v(sw)",
            "  set wr_singlescale",
            "  set wr_vecnames",
            "  wrdata fc_5level_ngspice.txt ufc0 ufc1 ufc2 uout il usw",
            "  quit",
            ".endc",
            ".end",
        ]
    return "\n".join(lines) + "\n"


def find_ngspice(explicit=None):
    for cand in (explicit, os.environ.get("NGSPICE"), shutil.which("ngspice")):
        if cand and Path(cand).exists():
            return cand
    raise SystemExit("ngspice not found: put it on PATH, set NGSPICE=..., or pass --ngspice")


QSPICE_EXE = "/mnt/c/Program Files/QSPICE/QSPICE80.exe"


def find_qspice(explicit=None):
    for cand in (explicit, os.environ.get("QSPICE"), QSPICE_EXE, shutil.which("QSPICE80.exe")):
        if cand and Path(cand).exists():
            return cand
    raise SystemExit("QSPICE80.exe not found: set QSPICE=..., or pass --qspice")


def read_qraw(path):
    """QSPICE .qraw : text header up to 'Binary:\n', then all-double rows."""
    import numpy as np

    raw = Path(path).read_bytes()
    marker = b"Binary:\n"
    split = raw.index(marker) + len(marker)
    header = raw[:split].decode("latin-1")
    names = re.findall(r"^\t\d+\t(\S+)\t", header, re.M)
    data = np.frombuffer(raw[split:], dtype="<f8")
    return dict(zip(names, data.reshape(-1, len(names)).T))


def qspice_to_txt(qraw, txt, tstep):
    """Resample a .qraw to tstep and write it in the ngspice wrdata layout."""
    import numpy as np

    d = read_qraw(qraw)
    t = np.arange(0.0, d["Time"][-1], tstep)
    cols = {
        "ufc0": d["V(p1)"] - d["V(n1)"],
        "ufc1": d["V(p2)"] - d["V(n2)"],
        "ufc2": d["V(p3)"] - d["V(n3)"],
        "uout": d["V(out)"],
        "il": d["I(L1)"],
        "usw": d["V(sw)"],
    }
    out = np.column_stack([t] + [np.interp(t, d["Time"], c) for c in cols.values()])
    np.savetxt(txt, out, fmt="%.9e", header=" ".join(["time", *cols]), comments="")


def compare(vhdl_dat, spice_txt, out_png, zoom=16e-6, simulator="ngspice"):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    v = pd.read_csv(vhdl_dat, sep=r"\s+", comment="#")
    s = pd.read_csv(spice_txt, sep=r"\s+")

    pairs = [("T_i0", "ufc0", "Flying-cap 1"), ("T_i1", "ufc1", "Flying-cap 2"),
             ("T_i2", "ufc2", "Flying-cap 3"), ("T_i5", "uout", "Output voltage"),
             ("B_u0", "il", "Inductor current")]

    t_end = min(v["time"].iloc[-1], s["time"].iloc[-1])
    v = v[v["time"] <= t_end]
    print(f"{'signal':<18}{'vhdl end':>12}{simulator + ' end':>14}{'max |diff|':>13}")
    for vc, sc, name in pairs:
        interp = np.interp(v["time"], s["time"], s[sc])
        print(f"{name:<18}{v[vc].iloc[-1]:12.4f}{interp[-1]:14.4f}"
              f"{np.max(np.abs(interp - v[vc])):13.4f}")

    fig, axes = plt.subplots(3, 1, figsize=(11, 11))
    for vc, sc, name in pairs[:4]:
        l, = axes[0].plot(s["time"], s[sc], lw=0.8, label=f"{name} {simulator}")
        axes[0].plot(v["time"], v[vc], "--", lw=0.8, color=l.get_color(), label=f"{name} VHDL")
    axes[1].plot(s["time"], s["il"], lw=0.8, label=simulator)
    axes[1].plot(v["time"], v["B_u0"], "--", lw=0.8, label="VHDL")
    axes[1].sharex(axes[0])

    # last few switching periods : the VHDL samples land on the switching
    # instants, so draw them as markers over the continuous ngspice traces
    t0 = t_end - zoom
    sz = s[s["time"] >= t0]
    vz = v[v["time"] >= t0]
    for vc, sc, name in pairs[:3]:
        l, = axes[2].plot(sz["time"], sz[sc], lw=0.8, label=f"{name} {simulator}")
        axes[2].plot(vz["time"], vz[vc], "o", ms=3, color=l.get_color(), label=f"{name} VHDL")

    axes[0].set_ylabel("Voltage [V]")
    axes[1].set_ylabel("Inductor current [A]")
    axes[2].set_ylabel("Flying-cap voltage [V]")
    axes[2].set_xlabel("Time [s]")
    axes[0].set_title(f"5-level flying-capacitor converter, {simulator} vs VHDL switching model")
    axes[2].set_title(f"Last {zoom*1e6:g} us")
    for a in axes:
        a.legend(fontsize=7, ncol=2)
        a.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    print(f"saved {out_png}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tstop", type=float, default=400e-3, help="simulated time [s], default 400 ms as the VHDL tb")
    ap.add_argument("--tmax", type=float, default=20e-9, help="max internal time step [s]")
    ap.add_argument("--tstep", type=float, default=1e-6, help="output sample interval [s] (results are interpolated to it)")
    ap.add_argument("--outdir", default=str(Path(__file__).parent / "out"), help="where the netlist and results go")
    ap.add_argument("--simulator", choices=("ngspice", "qspice"), default="ngspice")
    ap.add_argument("--t-rise", type=float, default=1e-10,
                    help="gate edge time [s] ; slower edges overlap the cell switches "
                         "and bleed the flying caps")
    ap.add_argument("--ngspice", help="ngspice executable")
    ap.add_argument("--qspice", help="QSPICE80.exe (default: the Windows install, run through WSL)")
    ap.add_argument("--no-run", action="store_true", help="don't run the simulator ; with --compare, re-plot the previous results")
    ap.add_argument("--compare", metavar="VHDL_DAT", help="compare against fc_5level_tb.dat and plot")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sim = args.simulator
    cir = outdir / f"fc_5level_{sim}.cir"
    result = outdir / f"fc_5level_{sim}.txt"
    if not (args.no_run and args.compare):
        cir.write_text(netlist(PARAMS, args.tstop, args.tmax, args.tstep, args.t_rise, simulator=sim))
        print(f"wrote {cir}")
    if not args.no_run:
        t0 = time.monotonic()
        if sim == "ngspice":
            subprocess.run([find_ngspice(args.ngspice), "-b", "-o", "fc_5level.log", cir.name],
                           cwd=outdir, check=True)
        else:
            # a Windows exe under WSL : run it in outdir so it finds the
            # netlist by its bare name, wherever outdir is
            subprocess.run([find_qspice(args.qspice), "-binary", cir.name],
                           cwd=outdir, check=True, stdout=subprocess.DEVNULL)
        print(f"{sim} took {time.monotonic() - t0:.1f} s")
        if sim == "qspice":
            qspice_to_txt(cir.with_suffix(".qraw"), result, args.tstep)
        print(f"wrote {result}")

    if args.compare:
        compare(args.compare, result, outdir / f"fc_5level_{sim}_vs_vhdl.png", simulator=sim)


if __name__ == "__main__":
    main()
