import enum
from threading import local
import numpy as np
import columnplots as clp
import matplotlib.pyplot as plt

au2fs = 0.02418884254
Angstrom2AU = 1.8897259885789


methods = [ "../proton_transfer/neo_eh_3cen/",
            "../proton_transfer/neo_boeh_3cen/", 
            "../proton_transfer/neo_boeh_3cen_smalldt/",
            "../proton_transfer/eh_cl/", 
            "../proton_transfer/boeh_cl/"]
dts = [0.4, 4.2, 0.4, 0.4, 4.2]

glob_data = {}

for idx, dir in enumerate(methods):
    # open the file for energy dynamics
    filename = "%s/e_tot.txt" %(dir)
    print("working under %s" %filename)
    data = np.loadtxt(filename)
    # create a dictionary to store important info
    local_data = {}
    local_data["t"] = data[:,0] * dts[idx] * au2fs
    local_data["e_tot"] = data[:,1] - data[0,1]
    # do an additional truncation to keep the initial 100 fs
    nidx = np.sum(local_data["t"] < 100)
    local_data["t"] = local_data["t"][:nidx]
    local_data["e_tot"] = local_data["e_tot"][:nidx]
    # calculate distance to OA and OD
    glob_data[dir] = local_data

# Now we do simple plotting
colors = [clp.black, clp.red, clp.yellow]
labels_eh2 = [r"RT-NEO-Ehrenfest ($\Delta t_{\rm q}$= 0.010 fs)", r"BO-RT-NEO-Ehrenfest ($\Delta t_{\rm q}$= 0.102 fs)", r"BO-RT-NEO-Ehrenfest ($\Delta t_{\rm q}$= 0.010 fs)"]
labels_cl = [r"Ehrenfest classical proton ($\Delta t_{\rm q}$= 0.010 fs)", r"BOMD classical proton ($\Delta t_{\rm q}$= 0.102 fs)"]

axes = clp.initialize(2, 1, width=4.3, height=4.3*0.618*2, LaTeX=True, fontsize=12, 
    labelsize=14, labelthem=True, labelthemPosition=[-0.02, 1.02], sharey=True)

xs, ys = [], []
for idx2, dir in enumerate(methods[0:3]):
    coeff = 1 if idx2 == 1 else 1.0
    filename = dir
    xs.append(glob_data[filename]["t"]*coeff)
    ys.append(glob_data[filename]["e_tot"])
clp.plotone(xs, ys, axes[0], colors=colors, labels=labels_eh2, lw=1.5,
        showlegend=True,
        xlim=[0,100], ylim=[-6e-4, 9.4e-4], legendFontSize=9,
        ylabel="relative energy [a.u.]", alphaspacing=0.)

xs, ys = [], []
for dir in methods[3:]:
    filename = dir
    xs.append(glob_data[filename]["t"])
    ys.append(glob_data[filename]["e_tot"])
clp.plotone(xs, ys, axes[1], colors=colors, labels=labels_cl, lw=1.5,
    showlegend=True,
    xlabel="time [fs]", legendFontSize=9,
    xlim=[0,100],
    ylabel="relative energy [a.u.]", alphaspacing=0.)

clp.adjust(tight_layout=True, savefile="figs2_pt_energy.pdf")
