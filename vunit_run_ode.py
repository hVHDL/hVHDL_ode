#!/usr/bin/env python3

from pathlib import Path
from vunit import VUnit

# ROOT
ROOT = Path(__file__).resolve().parent
VU = VUnit.from_argv()


#this is obsolete and will be rewritten
ode = VU.add_library("ode")
ode.add_source_files(ROOT / "write_pkg.vhd")
ode.add_source_files(ROOT / "ode_solvers/ode_pkg.vhd")
ode.add_source_files(ROOT / "testbenches/lcr_simulation_rk4_tb.vhd")

VU.main()
