import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import models
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))


df = pd.read_table(parent+"/data/meyer.txt",delimiter=" ",header=None,skiprows=1)
init = np.array([0.02,4000,250])
popt, pcov = curve_fit(models.meyer,df[0].values,df[1].values,p0=init)

with open(parent+"/output/sol_ls_meyer.txt","w") as f:
    for i in range(3):
        f.write("%f\n" % popt[i])

t = np.linspace(df[0].values[0],df[0].values[-1],1000)

plt.plot(df[0].values,df[1].values,"ko",ms=2)
plt.plot(t,models.meyer(t,*popt),lw=1)
plt.show()    
print(popt)

