import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt("cs8_100000_40")
plt.yscale("log")
t = np.array([1, 2, 3, 4, 5])
x = np.array([5e-6, 5e-7, 5e-8, 5e-9, 5e-10])
y = np.array([5e-7, 5e-8, 5e-9, 5e-10, 5e-11])
#plt.plot(t, x)
#plt.plot(t, x-y)
#plt.plot(t, x+y)
plt.errorbar(data.T[0], data.T[1], yerr=data.T[2])
plt.show()
