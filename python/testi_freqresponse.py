from freq_response import frequency_response
import pandas as pd
import matplotlib.pyplot as pyplot
import numpy as np


ode_data = pd.read_csv("./fc_4level_tb.dat", sep='\s+') 
f,y,_ = frequency_response(ode_data["B_u3"], ode_data["B_u4"] , fs = 500e3, nperseg = 10000) 
pyplot.plot(f,20*np.log10((np.abs (y))))
pyplot.show()
