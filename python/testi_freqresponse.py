from freq_response import freq_response
import pandas as pd
import matplotlib.pyplot as pyplot
import numpy as np
from scipy.signal import StateSpace, freqresp

C = 10e-6
L = 10e-6
r = 10e-3

ssA = np.array([[0, 1/C],
              [-1/L, -r/L]])

ssB = np.array([[0],
              [1/L]])

ssUC = np.array([[1, 0]])
ssIL = np.array([[0, 1]])

iL_model = StateSpace(ssA, ssB, ssIL, np.array([0]))
uC_model = StateSpace(ssA, ssB, ssUC, np.array([0]))


ode_data = pd.read_csv("./fc_4level_tb.dat", sep='\s+') 
f,uC,_ = freq_response(ode_data["B_u4"], ode_data["B_u0"], fs = 1000e3, nperseg = 10000) 
f,iL,_ = freq_response(ode_data["B_u4"], ode_data["T_i0"], fs = 1000e3, nperseg = 10000) 

frequencies = f*2*np.pi

fig, (ax1, ax2) = pyplot.subplots(2, 1
                                  ,sharex=True
                                  ,figsize=(8, 4))

w, ssIL = freqresp(iL_model, frequencies)
w, ssUC = freqresp(uC_model, frequencies)

ax1.semilogx(f,20*np.log10((np.abs (iL))))
ax1.semilogx(f,20*np.log10((np.abs (ssIL))))
ax1.semilogx(f,20*np.log10((np.abs (uC))))
ax1.semilogx(f,20*np.log10((np.abs (ssUC))))
ax1.set_xlim(1000,500e3)
ax1.set_ylim(-60,50)
ax1.grid(True)
ax2.semilogx(f,np.degrees((np.angle (iL))))
ax2.semilogx(f,np.degrees((np.angle (ssIL))))
ax2.semilogx(f,np.degrees(np.angle (uC*1j))-90)
ax2.semilogx(f,np.degrees((np.angle (ssUC))))
ax2.grid(True)
pyplot.tight_layout()
pyplot.show()
