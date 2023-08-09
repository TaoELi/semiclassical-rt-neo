import numpy as np
import columnplots as clp

def get_charge(pathname):
    data = np.loadtxt(pathname + "chg.txt")[:850,:]
    time = np.arange(len(data[:,0])) * 0.1 * 0.02418884254
    charge = np.sum(data, axis=1)
    # charge -= charge[0]
    charge_A = data[:,0]
    charge_B = data[:,1]
    return time, charge, charge_A, charge_B

def get_mo_coeff(pathname, states=[67, 70, 108, 116, 120]):
    states = np.array(states, dtype=np.int64)
    data = np.loadtxt(pathname + "mo_coeff.txt")[:850,:]
    print("shape is", data.shape)
    time = np.arange(len(data[:,0])) * 0.1 * 0.02418884254
    data -= data[0,:]
    data_states = [data[:,n-1] for n in states]
    time_states = [time] * len(states)
    return time_states, data_states

axes = clp.initialize(1, 2, width=8.4, height=4.3*0.618, LaTeX=False, fontsize=12,
                      labelsize=15, labelthem=True, labelthemPosition=[0.1, 0.96])

# Fig A
time_cl, charge_cl, cA_cl, cB_cl = get_charge(pathname="../run_pbe_6-31gd/eq_rt_pulse_2e-2_rerun_fix/")
clp.plotone([time_cl], [charge_cl], axes[0], colors=[clp.black], labels=["RT-TDDFT"],
            xlim=[0,2],
            xlabel="time [fs]", ylabel="CHELPG charge of H$_2$ [a.u.]")

# Fig B
states = [67, 70, 108, 116, 120]
xs, ys = get_mo_coeff("../"
                      "run_pbe_6-31gd/eq_rt_pulse_2e-2_rerun_fix/", states=states)
labels = ["MO {}".format(x) for x in states]
colors = ["k", "k--", clp.red, clp.sky_blue, clp.dark_green]
clp.plotone(xs, ys, axes[1],
            xlabel="time [fs]", ylabel="population difference [a.u.]", xlim=[0,2],
            colors=colors, labels=labels, lw=1.6)

clp.adjust(tight_layout=True, savefile="figs7_mo_coeff_diff.png")