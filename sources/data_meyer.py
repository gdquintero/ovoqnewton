import numpy as np
import matplotlib.pyplot as plt
import os

cwd = os.getcwd()
parent =  os.path.abspath(os.path.join(cwd,os.pardir))

m = 16

t = np.empty(m)

for i in range(m):
    t[i] = 45 + 5 * (i+1)

y = np.array([
    34780.0, 28610.0, 23650.0, 19630.0, 16370.0, 13720.0, 11540.0, 9744.0,\
     8261.0,  7030.0,  6005.0,  5147.0,  4427.0,  3820.0,  3307.0, 2872.0
])

# Inject 4 outliers by corrupting the observations at indices 4, 8, 12, 16
# (1-based) with a multiplicative factor of 0.5.
outlier_idx = np.array([4, 8, 12, 16]) - 1   # 0-based
y[outlier_idx] = 0.5 * y[outlier_idx]

with open(parent+"/data/meyer.txt","w") as f:
    f.write("%i\n" % m)
    for i in range(m):
        f.write("%f %f\n" % (t[i],y[i]))

plt.plot(t,y,"o")
plt.show()