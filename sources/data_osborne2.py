import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

y = np.array([
    1.366,1.191,1.112,1.013,0.991,0.885,0.831,0.847,0.786,0.725,\
    0.746,0.679,0.608,0.655,0.616,0.606,0.602,0.626,0.651,0.724,\
    0.649,0.649,0.694,0.644,0.624,0.661,0.612,0.558,0.533,0.495,\
    0.500,0.423,0.395,0.375,0.372,0.391,0.396,0.405,0.428,0.429,\
    0.523,0.562,0.607,0.653,0.672,0.708,0.633,0.668,0.645,0.632,\
    0.591,0.559,0.597,0.625,0.739,0.710,0.729,0.720,0.636,0.581,\
    0.428,0.292,0.162,0.098,0.054
])

sup = 6.4

told = np.linspace(0,sup,len(y))

plt.plot(told,y,"ko")


y = np.insert(y,4,1.418)
y = np.insert(y,9,1.146)
y = np.insert(y,15,1.132)
y = np.insert(y,20,1.071)
y = np.insert(y,22,1.175)
y = np.insert(y,24,1.266)
y = np.insert(y,33,0.979)
y = np.insert(y,40,0.867)
y = np.insert(y,42,1.109)
y = np.insert(y,46,1.211)
y = np.insert(y,47,0.942)
y = np.insert(y,67,1.338)
y = np.insert(y,70,1.148)


samples = len(y)
t = np.linspace(0,sup,samples)

plt.plot(t[4],y[4],"ro")
plt.plot(t[9],y[9],"ro")
plt.plot(t[15],y[15],"ro")
plt.plot(t[20],y[20],"ro")
plt.plot(t[22],y[22],"ro")
plt.plot(t[24],y[24],"ro")
plt.plot(t[33],y[33],"ro")
plt.plot(t[40],y[40],"ro")
plt.plot(t[42],y[42],"ro")
plt.plot(t[46],y[46],"ro")
plt.plot(t[47],y[47],"ro")
plt.plot(t[67],y[67],"ro")
plt.plot(t[70],y[70],"ro")

plt.xticks([])
plt.yticks([])

plt.savefig(parent+"/images/osborne-data.pdf",bbox_inches="tight")



with open(parent+"/data/osborne2.txt","w") as f:
    f.write("%i\n" % samples)
    for i in range(samples):
        f.write("%f %f\n" % (t[i],y[i]))

# plt.plot(t,y,"o")
# plt.show()