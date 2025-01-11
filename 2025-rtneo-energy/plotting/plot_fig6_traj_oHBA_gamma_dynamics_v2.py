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

methods = ["../data/pgrad_0_mv_13_rescf_0_dtsm/"]

figure, ax = clp.initialize(1, 1, width=4.3, height=4.3*0.618, LaTeX=True, fontsize=12,
    labelsize=14, return_fig_args=True)


y = np.loadtxt(methods[0] + "/tau_dynamics.txt")
x = np.linspace(0,  np.size(y) * 0.04 * au2fs, np.size(y))
clp.plotone([x], [y-1], ax, colors=[clp.red], lw=0.8, ylim=[-0.48, 0.48], xlim=[0, 27],
                xlabel="time [fs]", ylabel=r"rescaling factor $\gamma$", showlegend=False)

clp.adjust(tight_layout=True, savefile="traj_ohba_gamma.pdf")
