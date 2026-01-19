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

ssC = np.array([[0, 1]])  # Output: capacitor voltage

ss_model = StateSpace(ssA, ssB, ssC, np.array([0]))


ode_data = pd.read_csv("./fc_4level_tb.dat", sep='\s+') 
f,y,_ = freq_response(ode_data["B_u4"], ode_data["T_i0"], fs = 1000e3, nperseg = 40000) 

frequencies = f*2*np.pi

fig, (ax1, ax2) = pyplot.subplots(2, 1
                                  ,sharex=True
                                  ,figsize=(8, 4))

w, H = freqresp(ss_model, frequencies)

ax1.semilogx(f,20*np.log10((np.abs (y))))
ax1.semilogx(f,36+20*np.log10((np.abs (H))))
ax1.set_xlim(1000,500e3)
ax1.set_ylim(0,80)
ax1.grid(True)
ax2.semilogx(f,np.degrees((np.angle (y))))
ax2.semilogx(f,np.degrees((np.angle (H))))
ax2.grid(True)
pyplot.tight_layout()
pyplot.show()
