import enum
from threading import local
import numpy as np
import columnplots as clp
import matplotlib.pyplot as plt

au2fs = 0.02418884254
Angstrom2AU = 1.8897259885789


def get_distance(data):
    print("data shape", np.shape(data))
    OD = data[:, 18:21]
    OA = data[:, 24:27]
    H = data[:, -3:]
    # calculate HOD, HOA distance
    R_HOD = np.sqrt(np.sum((H - OD)**2, axis=1))
    R_HOA = np.sqrt(np.sum((H - OA)**2, axis=1))
    return R_HOD, R_HOA


methods = ["../data/4cen_rescf_0_rerun", "../data/pgrad_0_mv_13_rescf_0_dtsm/"]

dts = [0.4, 0.04]

methods_bc = methods[1:]

dt_lst = ["0.08"]

glob_data = {}

for idx, dir in enumerate(methods):
    for dt in dt_lst:
        # open the file for proton dipole moment
        filename = "%s/config.txt" %(dir)
        print("working under %s" %filename)
        data = np.loadtxt(filename)
        # create a dictionary to store important info
        local_data = {}
        local_data["t"] = np.arange(np.shape(data)[0]) * 4.0 * au2fs
        R_HOD, R_HOA = get_distance(data)
        local_data["R_HOD"] = R_HOD
        local_data["R_HOA"] = R_HOA
        # calculate distance to OA and OD
        glob_data[filename] = local_data

for idx, dir in enumerate(methods_bc):
    for dt in dt_lst:
        # open the file for proton dipole moment
        filename = "%s/config_bc.txt" %(dir)
        print("working under %s" %filename)
        data = np.loadtxt(filename)
        # create a dictionary to store important info
        local_data = {}
        local_data["t"] = np.arange(np.shape(data)[0]) * 4.0 * au2fs
        R_HOD, R_HOA = get_distance(data)
        local_data["R_HOD"] = R_HOD
        local_data["R_HOA"] = R_HOA
        # calculate distance to OA and OD
        glob_data[filename] = local_data

# get energy of each method
for idx, dir in enumerate(methods):
    for dt in dt_lst:
        # open the file for proton dipole moment
        filename = "%s/e_tot.txt" %(dir)
        print("working under %s" %filename)
        data = np.loadtxt(filename)
        filename2 = "%s/e_kq.txt" %(dir)
        print("working under %s" %filename2)
        data2 = np.loadtxt(filename2)
        # create a dictionary to store important info
        local_data = {}
        local_data["te"] = data[:,0] * dts[idx] * au2fs
        local_data["e_tot"] = data[:,1] - data[0,1]
        local_data["e_ext"] = local_data["e_tot"] + data2[:,1]
        glob_data[dir] = local_data

# Now we do simple plotting

colors_fig1 = [clp.black, clp.red, "0.5"]
colors_fig2 = [clp.black, clp.red]
linestyles_dashed = ["-", "-", "-", "--", "--", "--"]

labels_fig1 = [r"FPB", r"TPB $\gamma$-thermostat", r"TPB basis center"]

filenames_fig1 = [methods[0] + "/config.txt", methods[1] + "/config.txt", methods[1] + "/config_bc.txt"]
coeffs = [1.0, 0.1, 0.1]

figure, axes = clp.initialize(3, 1, width=4.3, height=4.3*0.618*3, LaTeX=True, fontsize=13,
    labelsize=14, labelthem=True, labelthemPosition=[-0.02, 1.02],
    sharex=True, 
    sharey=False, return_fig_args=True)


xs, y1s, y2s = [], [], []
for idx, filename in enumerate(filenames_fig1):
    coeff = coeffs[idx]
    xs.append(glob_data[filename]["t"]*coeff)
    y1s.append(glob_data[filename]["R_HOD"])
    y2s.append(glob_data[filename]["R_HOA"])
clp.plotone(xs, y1s, axes[0], colors=colors_fig1, labels=labels_fig1, lw=1.5,
        showlegend=True, 
        xlim=[0,27], ylim=[0.9, 1.79],  legendFontSize=8,
        ylabel="H-O [Angstrom]", alphaspacing=0.1)
clp.plotone(xs, y2s, axes[0], colors=colors_fig1, labels=labels_fig1, lw=1.5,
        linestyles=linestyles_dashed, #xlabel="time [fs]",
        showlegend=False, alphaspacing=0.1)

#axes[0].annotate("basis function center", xy=(19.34, 1.63), xytext=(5.5,1.68), arrowprops=dict(facecolor="0.5", edgecolor="0.5", shrink=0.05), color="0.5", fontsize=11)

# plot energy
x3s, y3s, y4s  = [], [], []
for dir in methods:
    x3s.append(glob_data[dir]["te"][::10])
    y3s.append(glob_data[dir]["e_tot"][::10]*1e3)
    y4s.append(glob_data[dir]["e_ext"][::10]*1e3)

clp.plotone(x3s, y3s, axes[1], colors=colors_fig2, labels=labels_fig1, lw=1.5, ylim=[-8.0, 1],
        showlegend=False, legendFontSize=10, alphaspacing=0.05,  ylabel=r"$\Delta E_{\rm tot}$ [$\times 10^{-3}$a.u.]")

clp.plotone(x3s, y4s, axes[2], colors=colors_fig2, labels=labels_fig1, lw=1.5, ylim=[-2.0, 12.0],
        showlegend=False, alphaspacing=0.05, xlabel="time [fs]", ylabel=r"$\Delta E_{\rm ext}$ [$\times 10^{-3}$a.u.]")

clp.adjust(tight_layout=True, savefile="traj_ohba_sub.pdf")
