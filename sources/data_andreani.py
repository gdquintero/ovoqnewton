import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

m = 46

t = np.empty(m)

for i in range(m):
    t[i] = -1 + 0.1 * i

y = np.array([
    -5.8,-5.159,-4.232,-3.413,-2.696,-1.675,10,10,10,10,10,10,10,10,10,10,0.536,0.473,\
    -0.008,0.299,0.2,-0.299,0.008,-0.073,-0.536,-0.575,-0.584,-0.557,-0.488,0.029,-0.2,\
    0.031,0.328,1.097,1.144,1.675,2.296,3.413,4.232,5.159,6.2,6.961,8.648,9.667,11.224,12.925
])

print(len(y))

with open(parent+"/data/andreani.txt","w") as f:
    f.write("%i\n" % m)
    for i in range(m):
        f.write("%f %f\n" % (t[i],y[i]))

# plt.plot(t,y,"o")
# plt.show()