#!/usr/bin/env python
"""Run the QSPICE LC-filter AC analysis and overlay it on a measurement.

Netlists and simulates misc/emi_filter_model.qsch (a `.ac` sweep of the
LC EMI filter, V1 = 1 V AC), then plots the simulated transfer function
V(node)/V(source) against the measured Bode response in
misc/L1C1_resp.csv (Digilent WaveForms Network Analyzer export, the
DUT-relative "Channel 2" magnitude/phase).

Needs QSPICE installed (QUX.exe / QSPICE64.exe) and the PyQSPICE package.

    python misc/plot_lc_ac_response.py
    python misc/plot_lc_ac_response.py --node N06 --measured misc/L1C1_resp_50k1M.csv

To compare against the VHDL simulation of the same circuit instead, dump
the frequency responses to CSV in the format python/test_plot.py reads as
a saved frequency response, then overlay them on the VUnit run's output:

    python misc/plot_lc_ac_response.py --save-measured misc/lc_meas_resp.csv \\
                                       --save-sim misc/lc_qspice_resp.csv --no-plot
    python python/test_plot.py <vunit_lc_filter_ode_tb>.dat \\
           misc/lc_meas_resp.csv misc/lc_qspice_resp.csv --freq-only

--save-measured on its own does not need QSPICE.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PyQSPICE import clsQSPICE

HERE = Path(__file__).resolve().parent


def nodes_in_netlist(cir: Path) -> set[str]:
    nodes: set[str] = set()
    for line in cir.read_text(encoding="latin-1").splitlines():
        p = line.split()
        if len(p) >= 3 and p[0][:1].upper() in "RLCVID":
            nodes.update(p[1:3])
    return nodes


def ac_source_node(cir: Path) -> str:
    for line in cir.read_text(encoding="latin-1").splitlines():
        p = line.split()
        if p and p[0][:1].upper() == "V" and any("AC" in t.upper() for t in p[3:]):
            return p[1]
    return "N10"


def run_qspice_ac(qsch: Path, out_node: str):
    """Return (freq[Hz], complex transfer function V(out_node)/V(source))."""
    q = clsQSPICE(str(qsch))
    q.qsch2cir()
    cir = Path(q.path["cir"])

    available = nodes_in_netlist(cir)
    if out_node not in available:
        sys.exit(
            f"node {out_node!r} is not in the netlist.\n"
            f"available nodes: {', '.join(sorted(available))}"
        )
    src = ac_source_node(cir)

    q.cir2qraw()
    q.tstime(["qraw"])
    if not q.ts["qraw"]:
        sys.exit(f"QSPICE did not produce {cir.with_suffix('.qraw').name}")

    df = q.LoadQRAW([f"V({out_node})", f"V({src})"])
    if not str(q.sim.get("Type", "")).startswith("AC"):
        sys.exit(f"expected an AC analysis, got {q.sim.get('Type')!r} "
                 f"-- check the .ac directive in {qsch.name}")

    freq = df["Freq"].astype(float).to_numpy()
    tf = df[f"V({out_node})"].to_numpy() / df[f"V({src})"].to_numpy()
    return freq, tf


def load_measured(csv: Path):
    """Return (freq[Hz], magnitude[dB], phase[deg]) for the DUT channel."""
    lines = csv.read_text(encoding="latin-1").splitlines()
    hdr = next(i for i, l in enumerate(lines) if l.lower().startswith("frequency"))
    m = pd.read_csv(csv, skiprows=hdr, encoding="latin-1")
    mag = next(c for c in m.columns if "Channel 2 Magnitude" in c)
    pha = next(c for c in m.columns if "Channel 2 Phase" in c)
    return m[m.columns[0]].to_numpy(), m[mag].to_numpy(), m[pha].to_numpy()


def save_freq_response_csv(path: Path, label: str, freq, mag_db, phase_deg):
    """Write a frequency response in the format python/test_plot.py reads
    back as a saved response (its save_freq_response / is_saved_freq_response):
    a #FREQRESPONSE marker, a #CONFIG label line, then a
    frequency,magnitude_db,phase_deg table."""
    with open(path, "w", encoding="utf-8") as out:
        out.write("#FREQRESPONSE\n")
        out.write(f"#CONFIG label={label}\n")
        out.write("frequency,magnitude_db,phase_deg\n")
        for f, m, p in zip(freq, mag_db, phase_deg):
            out.write(f"{float(f)},{float(m)},{float(p)}\n")
    print(f"wrote {path}")


def resonance(freq, mag_db):
    i = int(np.argmax(mag_db))
    return float(freq[i]), float(mag_db[i])


def break_wraps(phase_deg):
    """NaN out the samples straddling a +/-180 wrap so the plotted phase
    line does not draw a near-vertical jump across the discontinuity."""
    p = np.asarray(phase_deg, dtype=float).copy()
    jump = np.abs(np.diff(p)) > 300.0
    p[1:][jump] = np.nan
    return p


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qsch", type=Path, default=HERE / "emi_filter_model.qsch")
    ap.add_argument("--measured", type=Path, default=HERE / "L1C1_resp.csv")
    ap.add_argument("--node", default="N04",
                    help="schematic node to probe (default: N04, the L1/C1 output)")
    ap.add_argument("--out", type=Path, default=HERE / "lc_ac_response.png")
    ap.add_argument("--no-show", action="store_true", help="save the figure but do not open a window")
    ap.add_argument("--no-plot", action="store_true",
                    help="skip the matplotlib figure entirely; only write the requested --save-* CSVs")
    ap.add_argument("--save-sim", type=Path, metavar="CSV",
                    help="write the QSPICE response to CSV in python/test_plot.py's "
                         "saved-frequency-response format")
    ap.add_argument("--save-measured", type=Path, metavar="CSV",
                    help="write the measured response to CSV in python/test_plot.py's "
                         "saved-frequency-response format (does not need QSPICE)")
    a = ap.parse_args()

    want_plot = not a.no_plot
    want_sim = want_plot or a.save_sim is not None
    want_meas = want_plot or a.save_measured is not None

    if a.no_plot and a.save_sim is None and a.save_measured is None:
        sys.exit("--no-plot with nothing to do: pass --save-sim and/or --save-measured")

    needed = [a.qsch] if want_sim else []
    if want_meas:
        needed.append(a.measured)
    for f in needed:
        if not f.is_file():
            sys.exit(f"not found: {f}")

    f_sim = tf = mag_sim = pha_sim = None
    if want_sim:
        f_sim, tf = run_qspice_ac(a.qsch, a.node)
        mag_sim = 20.0 * np.log10(np.abs(tf))
        pha_sim = np.degrees(np.angle(tf))
        if a.save_sim:
            save_freq_response_csv(a.save_sim, f"QSPICE V({a.node})/Vin",
                                   f_sim, mag_sim, pha_sim)

    f_meas = mag_meas = pha_meas = None
    if want_meas:
        f_meas, mag_meas, pha_meas = load_measured(a.measured)
        if a.save_measured:
            save_freq_response_csv(a.save_measured, f"measured {a.measured.name}",
                                   f_meas, mag_meas, pha_meas)

    if not want_plot:
        return

    fr_s, mg_s = resonance(f_sim, mag_sim)
    fr_m, mg_m = resonance(f_meas, mag_meas)
    print(f"simulated  V({a.node})/Vin  peak : {fr_s:9.1f} Hz   {mg_s:7.2f} dB")
    print(f"measured   {a.measured.name:<20} peak : {fr_m:9.1f} Hz   {mg_m:7.2f} dB")

    fig, (axm, axp) = plt.subplots(2, 1, sharex=True, figsize=(9, 7))

    axm.semilogx(f_sim, mag_sim, "-", lw=1.6, color="C0", label=f"QSPICE  V({a.node})/Vin")
    axm.semilogx(f_meas, mag_meas, "--", lw=1.3, color="C1", label=f"measured  {a.measured.name}")
    axm.axvline(fr_s, color="C0", ls=":", lw=0.8)
    axm.axvline(fr_m, color="C1", ls=":", lw=0.8)
    axm.set_ylabel("magnitude [dB]")
    axm.set_title("LC filter AC response — QSPICE vs measured")
    axm.grid(True, which="both", alpha=0.3)
    axm.legend()

    axp.semilogx(f_sim, break_wraps(pha_sim), "-", lw=1.6, color="C0")
    axp.semilogx(f_meas, break_wraps(pha_meas), "--", lw=1.3, color="C1")
    axp.set_ylabel("phase [deg]")
    axp.set_xlabel("frequency [Hz]")
    axp.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print(f"wrote {a.out}")
    if not a.no_show:
        plt.show()


if __name__ == "__main__":
    main()
