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

ode.add_source_files(ROOT / "ode_solvers/sort_pkg.vhd")

ode.add_source_files(ROOT / "testbenches/lcr_models_pkg.vhd")
ode.add_source_files(ROOT / "testbenches/lcr_simulation_rk4_tb.vhd")
ode.add_source_files(ROOT / "testbenches/lcr_3ph_tb.vhd")
ode.add_source_files(ROOT / "testbenches/lcr_3ph_adaptive_tb.vhd")

ode.add_source_files(ROOT / "ode_examples/grid_inverter_model_pkg.vhd")
ode.add_source_files(ROOT / "ode_examples/grid_inverter_control_tb.vhd")

ode.add_source_files(ROOT / "testbenches/sort_tb.vhd")
ode.add_source_files(ROOT / "testbenches/template_tb.vhd")

# converter switching-model testbenches
ode.add_source_files(ROOT / "testbenches/converter/multilevel/fc_4level_tb.vhd")
ode.add_source_files(ROOT / "testbenches/converter/multilevel/fc_4level_freq_tb.vhd")
ode.add_source_files(ROOT / "testbenches/converter/multilevel/fc_5level_tb.vhd")
ode.add_source_files(ROOT / "testbenches/converter/dcdc/buck_converter_tb.vhd")
ode.add_source_files(ROOT / "testbenches/converter/dcdc/boost_converter_tb.vhd")
ode.add_source_files(ROOT / "testbenches/converter/dcdc/dab_converter_tb.vhd")
ode.add_source_files(ROOT / "testbenches/converter/multiphase/boost_3ph_tb.vhd")
ode.add_source_files(ROOT / "testbenches/converter/multiphase/buck_3ph_tb.vhd")
ode.add_source_files(ROOT / "testbenches/converter/dcac/inverter_3ph_tb.vhd")
ode.add_source_files(ROOT / "testbenches/converter/dcac/inverter_3ph_svm_tb.vhd")

VU.main()
