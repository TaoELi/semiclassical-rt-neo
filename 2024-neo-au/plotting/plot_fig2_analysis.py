import enum
from threading import local
import numpy as np
import columnplots as clp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.pyplot as plt

AU2Ev = 27.211324570273

# prepare data
data = {}
data_energy = np.loadtxt("../data_proton_squeezing/au_junction/energy.txt")
data_dipole = np.loadtxt("../data_proton_squeezing/au_junction/dipole.txt")
# R(Au-Au)  Escf  Eneo  Eneo(noke)
e_base = np.min(data_energy[:,1:])
data["R(Au-Au)"], data["Escf"], data["Eneo"], data["Eneo(noke)"] = data_energy[:, 0], data_energy[:, 1], data_energy[:, 2], data_energy[:, 3]
data["Escf"] = (data["Escf"]-e_base)*AU2Ev
data["Eneo"] = (data["Eneo"]-e_base)*AU2Ev
data["Eneo(noke)"] = (data["Eneo(noke)"]-e_base)*AU2Ev
data["Eke"] = data["Eneo"] - data["Eneo(noke)"]
data["Ezpe"] = data["Eneo"] - data["Escf"]
# # R(Au-Au)  proton<x>  proton<y>  proton<z>  proton<x^2>  proton<y^2>  proton<z^2>
data["px"], data["py"], data["pz"], data["px2"], data["py2"], data["pz2"] = data_dipole[:, 1], data_dipole[:, 2], data_dipole[:, 3], data_dipole[:, 4], data_dipole[:, 5], data_dipole[:, 6]


axes = clp.initialize(3, 1, width=3.2, height=3.2*0.618*3, LaTeX=False, fontsize=12, sharex=True,
                      labelthem=True, labelsize=14, labelthemPosition=[-0.05, 1.05])

# plot energy relation
xs_e = [data["R(Au-Au)"]*2]*2
ys_e = [data["Escf"], data["Eneo"]]
colors = [clp.black, clp.red]
labels = ["conventional DFT", "NEO-DFT"]

clp.plotone(xs_e, ys_e, axes[0], colors=colors, labels=labels, ylabel="SCF energy [eV]")

# plot kinetic energy relation
xs_e = [data["R(Au-Au)"]*2]*2
ys_e = [data["Eke"], data["Ezpe"], data["Ezpe"]-data["Eke"]]
print("R(Au-Au) = \n", data["R(Au-Au)"]*2)
print("Eke = \n", data["Eke"])
colors = ["0.5", clp.navy_blue]
labels = ["proton KE", "proton ZPE", "diff"]

clp.plotone(xs_e, ys_e, axes[1], colors=colors, labels=labels, ylabel="proton energy [eV]")

xs_p = [data["R(Au-Au)"]*2]*3
px_var = data["px2"]-data["px"]**2
py_var = data["py2"]-data["py"]**2
pz_var = data["pz2"]-data["pz"]**2
print("R(Au-Au) = \n", data["R(Au-Au)"]*2)
print("x_var = \n", px_var)
print("y_var = \n", py_var)
print("z_var = \n", pz_var)
pvar_avg = (px_var+py_var+pz_var)/3.0
ys_dvar = [px_var, py_var, pz_var, pvar_avg]
ys_dvar = [y*100 for y in ys_dvar]
colors = [clp.red, clp.dark_green, clp.brown, clp.black]
labels = [r"$\sigma_x^2$", r"$\sigma_y^2$", r"$\sigma_z^2$", r"avg"]
clp.plotone(xs_p, ys_dvar, axes[2], colors=colors, linestyles=["--", "--", "--", "-"], labels=labels, xlabel=r"$R_{\rm{Au-Au}}$ distance [Å]",
            ylabel=r"proton variance [$\times 10^{-2}$ a.u.]", xlim=[2.75, 4.85])

clp.adjust(tight_layout=True, savefile="fig2_analysis.pdf")
