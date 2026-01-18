from freq_response import frequency_response
import pandas as pd
import matplotlib.pyplot as pyplot
import numpy as np


ode_data = pd.read_csv("./fc_4level_tb.dat", sep='\s+') 
f,y,_ = frequency_response(ode_data["B_u5"], ode_data["T_i0"] , fs = 500e3, nperseg = 50000) 

fig, (ax1, ax2) = pyplot.subplots(2, 1
                                  ,sharex=True
                                  , figsize=(8, 4))  # figsize=(width, height)

ax1.semilogx(f,20*np.log10((np.abs (y))))
ax1.set_xlim(1000,250e3)
ax1.grid(True)
ax2.semilogx(f,np.degrees((np.angle (y))))
ax2.grid(True)
pyplot.tight_layout()
pyplot.show()
