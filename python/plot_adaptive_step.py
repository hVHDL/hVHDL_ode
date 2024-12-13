

import pandas as pd
import matplotlib.pyplot as plt

vhdl_data = pd.read_csv('./lcr_3ph_adaptive_tb.dat', delim_whitespace=True)
qspice_data = pd.read_csv('./qspice_ref_models/3ph_lc.csv')

fig1, (axT, axB) = plt.subplots(2,1,sharex=True,constrained_layout=True)
vhdl_data.plot(ax=axT, x="time", y="T_u0", label="vhdl u0")
vhdl_data.plot(ax=axT, x="time", y="T_u1", label="vhdl u1")
vhdl_data.plot(ax=axT, x="time", y="T_u2", label="vhdl u2")

vhdl_data.plot(ax=axB, x="time", y="B_i0", label="vhdl i0")
vhdl_data.plot(ax=axB, x="time", y="B_i1", label="vhdl i1")
vhdl_data.plot(ax=axB, x="time", y="B_i2", label="vhdl i2")

qspice_data.plot(ax=axT, x="Time", y="V(uc1)", label="qspice u0")
qspice_data.plot(ax=axT, x="Time", y="V(uc2)", label="qspice u1")
qspice_data.plot(ax=axT, x="Time", y="V(uc3)", label="qspice u2")

qspice_data.plot(ax=axB, x="Time", y="I(L1)", label="qspice i0")
qspice_data.plot(ax=axB, x="Time", y="I(L2)", label="qspice i1")
qspice_data.plot(ax=axB, x="Time", y="I(L3)", label="qspice i2")

plt.show()
plt.close('all')
