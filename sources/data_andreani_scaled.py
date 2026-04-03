import numpy as np
import matplotlib.pyplot as plt
import os
import models
import random

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

size_img = 0.6
# plt.rcParams.update({'font.size': 11})
plt.rcParams['figure.figsize'] = [size_img * 6.4,size_img * 4.8]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def gen_data(m):
    t = np.linspace(-1,3.5,m)
    y = models.andreani(t,*xsol)
    random.seed(123456)
    np.random.seed(12)
    r = 0.5
    ind = []
    noutliers = 0

    for i in range(m):
        y[i] = y[i] + random.uniform(-r,r)
        
        if random.random() <= 0.1:
            noutliers += 1
            ind.append(i+1)
            operacion = np.random.choice([0,1],p=[0.2, 0.8])

            if operacion == 1:
                y[i] = random.uniform(y[i],15)
            else:
                y[i] = random.uniform(-6,y[i])     
       
    with open(parent+"/data/outliers_andreani_scaled_true"+str(m)+".txt","w") as f:
        f.write("%i\n" % noutliers)
        for i in range(noutliers):
            f.write("%i\n" % ind[i])

    with open(parent+"/data/andreani"+str(m)+".txt","w") as f:
        f.write("%i\n" % m)
        for i in range(m):
            f.write("%f %f\n" % (t[i],y[i]))


xsol = np.array([0,2,-3,1])

for n in [100,1000,10000,100000,1000000]:
    gen_data(n)

# gen_data(100)
# gen_data(1000)
# gen_data(10000)
# gen_data(100000)
# gen_data(1000000)



