import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import models
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

size_img = 0.6
plt.rcParams.update({'font.size': 11})
plt.rcParams['figure.figsize'] = [size_img * 6.4,size_img * 4.8]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

df_data = pd.read_table(parent+"/data/andreani.txt",delimiter=" ",header=None,skiprows=1,skipinitialspace=True)
df_sol = pd.read_table(parent+"/output/solution_andreani.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
df_outliers = pd.read_table(parent+"/output/outliers_andreani.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

t = np.linspace(df_data[0].values[0],df_data[0].values[-1],1000)
noutliers = df_outliers[0].values[0]
outliers = np.empty((2,noutliers))

for i in range(noutliers):
    outliers[0,i] = df_data[0].values[df_outliers[0].values[i+1]-1]
    outliers[1,i] = df_data[1].values[df_outliers[0].values[i+1]-1]
    # print(outliers[1][i])

plt.plot(df_data[0].values,df_data[1].values,"ko",ms=2)
plt.plot(t,models.andreani(t,*df_sol.values[0]),lw=1)
plt.plot(outliers[0],outliers[1],'ro',mfc='none',ms=5,mew=0.5)
plt.tick_params(axis='both',direction='in')
plt.xticks(np.arange(-1,3.1,1))
plt.yticks(np.arange(-6,12.1,6))
plt.ylim(-6.5,12)
plt.xlim(-1.1,3.6)
plt.savefig(parent+"/images/andreani_fitting.pdf",bbox_inches = "tight")
# plt.show()