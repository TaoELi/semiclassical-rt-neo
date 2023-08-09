import numpy as np
import columnplots as clp


def func_exp(x, A, k):
    return A*np.exp(k * x)


def get_fitted_rate(x, y):
    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(func_exp, x, y)
    return popt[1]

au2fs = 0.02418884254
au2angstrom = 0.529177
dt = 0.1

nstep_10fs = int(10.0 / au2fs / dt)

dirs = ["../rerun_largebasis/eq_rt_pulse_5e-3", "../rerun_largebasis/eq_neo_pulse_5e-3", "../rerun_largebasis/eq_rt_pulse_5e-3_D", "../rerun_largebasis/eq_neo_pulse_5e-3_D",
        "../rerun_largebasis/eq_rt_pulse_1.0e-2", "../rerun_largebasis/eq_neo_pulse_1.0e-2", "../rerun_largebasis/eq_rt_pulse_1.0e-2_D", "../rerun_largebasis/eq_neo_pulse_1.0e-2_D",
        "../rerun_largebasis/eq_rt_pulse_week", "../rerun_largebasis/eq_neo_pulse_week", "../rerun_largebasis/eq_rt_pulse_2e-2_D", "../rerun_largebasis/eq_neo_pulse_2e-2_D"
        ]

Al13_coords = np.array([ [-0.0076262410,     -0.0594344247,     -0.0612487404],
             [0.1101047247,      0.0390858125,      2.6016674667],
             [2.2970077678,      0.2548228141,     -1.3686533275],
             [1.7400738130,      1.6579651251,      0.9918966950],
             [-0.3410875708,     -2.3674204126,      1.2364042776],
             [1.0144132679,     -2.2350572960,     -1.2197570802],
             [-2.3116819862,     -0.3722848018,      1.2465533498],
             [-0.1230961130,     -0.1567459002,     -2.7242976040],
             [-1.7499048525,     -1.7826984255,     -1.1136573362],
             [-1.0274624357,      2.1166186879,      1.0958454890],
             [0.3222485310,      2.2498105967,     -1.3562881495],
             [2.1616752527,     -1.1179926394,      1.0752485190],
             [-2.1796729702,      0.9897352428,     -1.1998904537]])

Al13_xc, Al13_yc, Al13_zc = np.mean(Al13_coords, axis=0)

glob_data = {}

for dir in dirs:
    # open the file for proton dipole moment
    filename = "%s/xn1.txt" %(dir)
    print("working under %s" %filename)
    data1 = np.loadtxt(filename)
    filename2 = "%s/xn2.txt" %(dir)
    print("working under %s" %filename2)
    data2 = np.loadtxt(filename2)
    # Now let us calculate the bond length
    data = np.zeros((np.shape(data1)[0], 4))
    data[:,0] = data1[:,0]
    x1, y1, z1 = data1[:,1] * au2angstrom, data1[:,2] * au2angstrom, data1[:,3] * au2angstrom
    x2, y2, z2 = data2[:,1] * au2angstrom, data2[:,2] * au2angstrom, data2[:,3] * au2angstrom
    data[:,1] = np.sqrt((x1 - x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    data[:,2] = np.sqrt((x1 - Al13_xc)**2 + (y1 - Al13_yc)**2 + (z1 - Al13_zc)**2)
    data[:,3] = np.sqrt((x2 - Al13_xc)**2 + (y2 - Al13_yc)**2 + (z2 - Al13_zc)**2)
    # create a dictionary to store important info
    local_data = {}
    local_data["t"] = data[:,0] * dt * au2fs
    local_data["r"] = data[:,1]
    local_data["separation1"] = data[:,2]
    local_data["separation2"] = data[:,3]
    # now we calculate the bond length dynamics rates
    k = get_fitted_rate(local_data["t"][:nstep_10fs], local_data["r"][:nstep_10fs] - local_data["r"][0])
    local_data["k"] = k
    # print("R displacement is ", local_data["r"][nstep_10fs] - local_data["r"][0])
    glob_data[dir] = local_data

# Now we do simple plotting

colors = [clp.black, clp.red, clp.black, clp.red]
linestyles = ["-", "-", "--", "--"]
labels = ["RT-Ehrenfest H$_2$", "RT-NEO-TDDFT H$_2$", "RT-Ehrenfest D$_2$", "RT-NEO-TDDFT D$_2$"]

axes = clp.initialize(3, 1, width=3.25, height=3.25*0.618*3, LaTeX=False, fontsize=10,
                      labelthem=True, labelthemPosition=[0.1, 0.95], labelsize=11)

xs, ys, separations_1, separations_2 = [], [], [], []

for dir in dirs:
    xs.append(glob_data[dir]["t"])
    ys.append(glob_data[dir]["r"])
    separations_1.append(glob_data[dir]["separation1"])
    separations_2.append(glob_data[dir]["separation2"])

# plot bond length dynamics
clp.plotone(xs[0:4], ys[0:4], axes[0], colors=colors, labels=labels, linestyles=linestyles,
            showlegend=True, xlim=[0, 20], legendloc="center right", legendFontSize=9,
            ylim=[0.71, 2.0],
            lw=1.5)

clp.plotone(xs[4:8], ys[4:8], axes[1], colors=colors, labels=labels, linestyles=linestyles,
            showlegend=False, xlim=[0, 20], # xlabel="time [fs]",
            ylabel=r"H-H or D-D interatomic distance [Å]",
            ylim=[0.71, 2.0],
            lw=1.5)

clp.plotone(xs[8:12], ys[8:12], axes[2], colors=colors, labels=labels, linestyles=linestyles,
            showlegend=False, xlim=[0, 20], xlabel="time [fs]",
            # ylabel="H$_{2}$ bond length [Angstrom]",
            ylim=[0.71, 2.0],
            lw=1.5)

E0_lst = ["5e-3", "1e-2", "2e-2"]
E0_lst = [51.42 * float(E0) for E0 in E0_lst]
for i, E0 in enumerate(E0_lst):
    axes[i].text(0.115, 0.87, r"$E_0$ = %.2f V/Å" %(E0), fontsize=10,
                 transform=axes[i].transAxes)

# finally, we insert some images
# clp.add_figure(imag_filename="Al13H2_eq.png", ax=axes[1], zoom=0.13, location=(5.2, 1.2))

clp.adjust(tight_layout=False, savefile="figs14_bond_length_dynamics_pulse_intensity_sc_largebasis.pdf")
