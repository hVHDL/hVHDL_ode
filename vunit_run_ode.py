#!/usr/bin/env python3

from pathlib import Path
from vunit import VUnit

# ROOT
ROOT = Path(__file__).resolve().parent
VU = VUnit.from_argv()

ode = VU.add_library("ode")
ode.add_source_files(ROOT / "write_pkg.vhd")
ode.add_source_files(ROOT / "ode_solvers/real_vector_pkg.vhd")
ode.add_source_files(ROOT / "ode_solvers/ode_pkg.vhd")
ode.add_source_files(ROOT / "ode_solvers/adaptive_ode_pkg.vhd")

ode.add_source_files(ROOT / "testbenches/lcr_simulation_rk4_tb.vhd")
ode.add_source_files(ROOT / "testbenches/lcr_3ph_tb.vhd")

VU.main()
