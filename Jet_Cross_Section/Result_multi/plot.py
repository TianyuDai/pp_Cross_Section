import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt("jet_multiplicity_err")
exp_data = np.loadtxt("exp")
def_data = np.loadtxt("MCplot")

plt.figure()
x = [2, 3, 4, 5, 6]
plt.errorbar(x, data.T[0], yerr=data.T[1], fmt='o', markersize=1, label='Pythia', capsize=2, elinewidth=0.5)
plt.yscale("log")
plt.xlim = (1, 7)
plt.errorbar(exp_data.T[1], exp_data.T[3], yerr=exp_data.T[4], fmt='o', markersize=1, label='ATLAS(2011)', capsize=2, elinewidth=0.5)
plt.errorbar(def_data.T[1], def_data.T[3], yerr=def_data.T[4], fmt='o', markersize=1, label='default', capsize=2, elinewidth=0.5)
plt.ylabel("$sigma$ [pb]")
plt.xlabel("$N_text{jet}$")
plt.legend()
plt.title("Inclusive jet multiplicity ($R=0.4$)")
plt.savefig("JetMulti.pdf")
plt.show()

