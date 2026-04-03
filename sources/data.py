import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

def plot_seropositive(sero,x,y):
    # t = np.linspace(0,70,1000)
    size_img = 0.6
    plt.rcParams.update({'font.size': 11})
    plt.rcParams['figure.figsize'] = [size_img * 6.4,size_img * 4.8]
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(x,y,"ko",ls=":",ms=3,lw=1)
    plt.tick_params(axis='both',direction='in')
    plt.xticks(np.arange(0,71,10))
    plt.yticks(np.arange(0.2,1.1,0.2))
    plt.ylim(0.0,1.05)
    plt.xlim(-4.0,70)

    if sero == "measles":
        plt.savefig(parent+"/images/sero_measles.pdf",bbox_inches = "tight") 
    elif sero == "mumps":
        plt.savefig(parent+"/images/sero_mumps.pdf",bbox_inches = "tight")
    else:
        plt.savefig(parent+"/images/sero_rubella.pdf",bbox_inches = "tight")

    # plt.show()
    plt.close()

age = np.array([
    1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,19,21,23,25,27,29,31,33,35,40,45,55,65
])

sero_measles = np.array([
    0.207,0.301,0.409,0.589,0.757,0.669,0.797,0.818,0.866,0.859,0.908,0.923,0.889,0.936,0.889,\
    0.898,0.959,0.957,0.937,0.918,0.939,0.967,0.973,0.943,0.967,0.946,0.961,0.968,0.968
])

sero_mumps = np.array([
    0.115,0.147,0.389,0.516,0.669,0.768,0.786,0.798,0.878,0.861,0.844,0.881,0.895,0.882,0.869,\
    0.895,0.911,0.920,0.915,0.950,0.909,0.873,0.880,0.915,0.906,0.933,0.917,0.898,0.839
])

sero_rubella = np.array([
    0.126,0.171,0.184,0.286,0.400,0.503,0.524,0.634,0.742,0.664,0.735,0.815,0.768,0.842,0.760,\
    0.869,0.844,0.852,0.907,0.935,0.921,0.896,0.890,0.949,0.899,0.955,0.937,0.933,0.917
])

samples = len(age)

# Adding outliers manually
deviation = 0.5

sero_measles[16] = deviation
sero_measles[17] = deviation
sero_measles[18] = deviation
sero_measles[19] = deviation

sero_mumps[16] = deviation
sero_mumps[17] = deviation
sero_mumps[18] = deviation
sero_mumps[19] = deviation

sero_rubella[16] = deviation
sero_rubella[17] = deviation
sero_rubella[18] = deviation
sero_rubella[19] = deviation

age_midpoint = np.empty(samples)
age_midpoint[:-1] = (age[:-1] + age[1:]) / 2
age_midpoint[-1]  = 70

with open(parent+"/data/seropositives.txt","w") as f:
    f.write("%i\n" % samples)
    for i in range(samples):
        f.write("%i %f %f %f %f\n" % (age[i],sero_measles[i],sero_mumps[i],sero_rubella[i],age_midpoint[i]))


plot_seropositive("measles",age,sero_measles)
plot_seropositive("mumps",age,sero_mumps)
plot_seropositive("rubella",age,sero_rubella)
