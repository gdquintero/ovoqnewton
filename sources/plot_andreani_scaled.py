import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import models
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

size_img = 0.6
# plt.rcParams.update({'font.size': 11})
plt.rcParams['figure.figsize'] = [size_img * 6.4,size_img * 4.8]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def plot_fit(n):

    data = parent+"/data/andreani"+str(n)+".txt"
    outliers_true = parent+"/data/outliers_andreani_scaled_true"+str(n)+".txt"
    sol_ls = parent+"/output/sol_ls_andreani"+str(n)+".txt"

    df_data = pd.read_table(data,delimiter=" ",header=None,skiprows=1,skipinitialspace=True)
    df_sol = pd.read_table(parent+"/output/solution_andreani_scaled.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
    df_outliers = pd.read_table(parent+"/output/outliers_andreani_scaled.txt",delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
    df_outliers_true = pd.read_table(outliers_true,delimiter=" ",header=None,skiprows=0,skipinitialspace=True)
    df_sol_ls = pd.read_table(sol_ls,delimiter=" ",header=None,skiprows=0,skipinitialspace=True)

    t = np.linspace(df_data[0].values[0],df_data[0].values[-1],1000)
    noutliers = df_outliers[0].values[0]
    outliers = np.empty((2,noutliers))

    for i in range(noutliers):
        outliers[0,i] = df_data[0].values[df_outliers[0].values[i+1]-1]
        outliers[1,i] = df_data[1].values[df_outliers[0].values[i+1]-1]

    n_data = len(df_data)

    for i in range(1,n_data+1):
        if i in df_outliers_true[0].values[1:]:
            plt.plot(df_data[0].values[i-1],df_data[1].values[i-1],"ro",ms=0.1)
        else:
            plt.plot(df_data[0].values[i-1],df_data[1].values[i-1],"ko",ms=0.1)
        

    plt.plot(outliers[0],outliers[1],'go',mfc='none',ms=2,mew=0.2)
    plt.plot(t,models.andreani(t,*df_sol.values[0]),lw=1,label="OVO")
    plt.plot(t,models.andreani(t,*df_sol_ls[0].values),lw=1,label="LS")
    plt.tick_params(axis='both',direction='in')
    plt.xticks(np.arange(-1,3.1,1))
    plt.xlim(-1.1,3.6)

    plt.yticks(np.arange(-6,16.1,3))
    plt.ylim(-7.5,16)
    
    # plt.legend(loc='best', fontsize='small', fancybox=True, framealpha=0.5)
    plt.savefig(parent+"/images/andreani_fitting"+str(n)+".pdf",bbox_inches = "tight")
    # plt.show()
    plt.close()

def plot_log(n):

    # fig,ax = plt.subplots(1, 1)
    if n == 100: 
        file = "andreani_scaled_log_100.txt"
    elif n == 1000:
        file = "andreani_scaled_log_1000.txt"
    elif n == 10000:
        file = "andreani_scaled_log_10000.txt"
    elif n == 100000:
        file = "andreani_scaled_log_100000.txt"
    else:
        file = "andreani_scaled_log_1000000.txt"

    df_data = pd.read_table(file,delimiter=" ",header=None,skipinitialspace=True)
    # ax.tick_params(axis='both',direction='in',which='both')
    # ax.loglog(df_data[1].values,df_data[2].values,"-o",color="darkgreen",lw=0.5,ms=1)
    # ax.set_xscale('linear')

    plt.plot(df_data[1].values,df_data[2].values,"-ko",lw=0.5,ms=1)


    if n == 100:
        plt.xticks(np.arange(5,15.1,2))
        plt.yticks(np.arange(0,8.1,1))
        plt.ylim(-0.5,8)
    elif n == 1000:
        plt.xticks(np.arange(50,150.1,20))
        plt.yticks(np.arange(0,6.1,1))
        plt.ylim(-0.5,6)
    elif n == 10000:
        plt.xticks(np.arange(500,1500.1,200))
        plt.yticks(np.arange(0,10.1,2))
        plt.ylim(-0.5,10)
    elif n == 100000:
        plt.xticks(np.arange(5000,15000.1,2000))
        plt.yticks(np.arange(0,10.1,2))
        plt.ylim(-0.5,10)
    else:
        plt.xticks(np.arange(50000,150000.1,20000))
        plt.yticks(np.arange(0,10.1,2))
        plt.ylim(-0.5,10)


    # plt.yscale("log")
    # plt.yticks(np.arange(10**(-1),10,4))
    plt.xlabel("Number of outliers $o$")
    plt.ylabel("$f(x^*)$")
    plt.tick_params(axis='both',direction='in')
    plt.savefig(parent+"/images/andreani_scaled"+str(n)+".pdf",bbox_inches="tight")
    # plt.show()
    plt.close()


# for n in [100,1000,10000,100000,1000000]:
#     plot_log(n)

plot_fit(100000)

