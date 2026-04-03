import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error
import models
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

def main(problem):
    if problem == 1:
        df = pd.read_table(parent+"/data/seropositives.txt",delimiter=" ",header=None,skiprows=1)

        # init_measles = np.array([0.197,0.287,0.021])
        # init_mumps = np.array([0.156,0.250,0.0])
        # init_rubella = np.array([0.0628,0.178,0.020])

        init_measles = np.ones(3)
        init_mumps   = np.ones(3)
        init_rubella = np.ones(3)

        popt_measles, pcov_measles  = curve_fit(models.F,df[0].values,df[1].values,p0=init_measles,bounds=(np.zeros(3),np.ones(3)))
        popt_mumps, pcov_mumps      = curve_fit(models.F,df[0].values,df[2].values,p0=init_mumps,bounds=(np.zeros(3),np.ones(3)))
        popt_rubella, pcov_rubella  = curve_fit(models.F,df[0].values,df[3].values,p0=init_rubella,bounds=(np.zeros(3),np.ones(3)))

        with open(parent+"/output/sol_ls_farrington.txt","w") as f:
            f.write("%f %f %f\n" % (popt_measles[0],popt_measles[1],popt_measles[2]))
            f.write("%f %f %f\n" % (popt_mumps[0],popt_mumps[1],popt_mumps[2]))
            f.write("%f %f %f\n" % (popt_rubella[0],popt_rubella[1],popt_rubella[2])) 

        # y_pred = np.empty((3,29))

        # y_pred[0,:] = models.F(df[0].values,*popt_measles)
        # y_pred[1,:] = models.F(df[0].values,*popt_mumps)
        # y_pred[2,:] = models.F(df[0].values,*popt_rubella)

        # error_measles =  mean_squared_error(df[1].values,y_pred[0,:])
        # error_mumps   =  mean_squared_error(df[2].values,y_pred[1,:])
        # error_rubella =  mean_squared_error(df[3].values,y_pred[2,:])
            
        # print("Mean squared error for Measles:","{:.3e}".format(error_measles))
        # print("Mean squared error for Mumps:","{:.3e}".format(error_mumps))
        # print("Mean squared error for Rubella:","{:.3e}".format(error_rubella))  

        print(popt_measles)
        print(popt_mumps)
        print(popt_rubella)

        t = np.linspace(df[0].values[0],df[0].values[-1],1000)

        plt.plot(df[0].values,df[1].values,"ko",ms=2)
        plt.plot(t,models.F(t,*popt_mumps),lw=1)
        plt.savefig(parent+"/images/farrington.pdf",bbox_inches="tight")

    elif problem == 2:
        df = pd.read_table(parent+"/data/osborne2.txt",delimiter=" ",header=None,skiprows=1)
        init = np.array([1.3,0.65,0.65,0.7,0.6,3.,5.,7.,2.,4.5,5.5])
        popt, pcov = curve_fit(models.osborne2,df[0].values,df[1].values,p0=init)
        
        with open(parent+"/output/sol_ls_osborne2.txt","w") as f:
            for i in range(11):
                f.write("%f\n" % popt[i])

        t = np.linspace(df[0].values[0],df[0].values[-1],1000)

        plt.plot(df[0].values,df[1].values,"ko",ms=2)
        plt.plot(t,models.osborne2(t,*popt),lw=1)
        plt.show()    
        print(popt)

    elif problem == 3:

        df = pd.read_table(parent+"/data/andreani.txt",delimiter=" ",header=None,skiprows=1)
        init = np.array([-1.0,-2.0,1.0,-1.0])
        popt, pcov = curve_fit(models.andreani,df[0].values,df[1].values,p0=init)
        
        with open(parent+"/output/sol_ls_andreani.txt","w") as f:
            for i in range(4):
                f.write("%f\n" % popt[i])

        t = np.linspace(df[0].values[0],df[0].values[-1],1000)

        plt.plot(df[0].values,df[1].values,"ko",ms=2)
        plt.plot(t,models.andreani(t,*popt),lw=1)
        plt.show()    
        print(popt)

    else:
        df = pd.read_table(parent+"/data/andreani1000000.txt",delimiter=" ",header=None,skiprows=1)
        init = np.array([-1.0,-2.0,1.0,-1.0])
        popt, pcov = curve_fit(models.andreani,df[0].values,df[1].values,p0=init)
        
        with open(parent+"/output/sol_ls_andreani1000000.txt","w") as f:
            for i in range(4):
                f.write("%f\n" % popt[i])

        t = np.linspace(df[0].values[0],df[0].values[-1],1000)

        plt.plot(df[0].values,df[1].values,"ko",ms=2)
        plt.plot(t,models.andreani(t,*popt),lw=1)
        # plt.show()    
        print(popt)

main(4)

