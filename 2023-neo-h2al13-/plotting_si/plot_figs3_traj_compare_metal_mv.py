import enum
from threading import local
import numpy as np
import columnplots as clp
import matplotlib.pyplot as plt

'''
Pade Transform
'''
au2fs = 0.02418884254

'''
Compare the performance of two methods under different dt
'''

methods = ["%s" %x for x in ["../run_pbe_6-31gd/eq_rt_pulse_2e-2_mv_partial/",
                             "../run_pbe_6-31gd/eq_rt_pulse_2e-2_mv_all/"
                             ]]
dt = 0.1

glob_data = {}

for dir_name in methods:
    # open the file for proton dipole moment
    filename = "%s/xn1.txt" %(dir_name)
    print("working under %s" %filename)
    data1 = np.loadtxt(filename)
    filename2 = "%s/xn2.txt" %(dir_name)
    print("working under %s" %filename2)
    data2 = np.loadtxt(filename2)
    # Now let us calculate the bond length
    data = np.zeros((np.shape(data1)[0], 2))
    data[:,0] = data1[:,0]
    x1, y1, z1 = data1[:,1], data1[:,2], data1[:,3]
    x2, y2, z2 = data2[:,1], data2[:,2], data2[:,3]
    data[:,1] = np.sqrt((x1 - x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    # create a dictionary to store important info
    local_data = {}
    local_data["t"] = data[:,0] * dt * au2fs
    local_data["r"] = data[:,1] * 0.529177
    glob_data[dir_name] = local_data

# Now we do simple plotting

colors = [clp.black, clp.sky_blue]
linestyles = ["-", "--"]
labels = ["RT-Ehrenfest H$_2$ fixed Al atoms", "RT-Ehrenfest H$_2$ moving Al atoms"]

ax = clp.initialize(1, 1, width=4.3, height=3., LaTeX=False, fontsize=12)

xs, ys = [], []
freq_lst, sp_lst = [], []

for dir_name in methods:
    xs.append(glob_data[dir_name]["t"])
    ys.append(glob_data[dir_name]["r"])

# plot bond length dynamics
clp.plotone(xs, ys, ax, colors=colors, labels=labels, linestyles=linestyles,
            showlegend=True, xlim=[0, 20], xlabel="time [fs]", legendFontSize=9,
            ylabel=r"H-H interatomic distance [Å]", ylim=[0.71, 2.6],
            lw=1.5)


clp.adjust(tight_layout=False, savefile="figs3_bond_length_dynamics_compare_metal_mv.pdf")
