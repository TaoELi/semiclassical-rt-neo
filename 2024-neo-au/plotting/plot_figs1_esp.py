import numpy as np
import columnplots as clp


def get_esp(filename):
    data = np.loadtxt(filename)
    x, esp = data[:,0], data[:,-1]
    return x, esp


dis_lst = [1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]

xs, ys = [], []
for idx, dis in enumerate(dis_lst):
    print("dealing with dis = ", dis)
    filename = "../data_proton_squeezing/au_junction_efield_grid/d_%.1f/esp.txt" % dis
    x, esp = get_esp(filename)
    xs.append(x)
    ys.append(esp)

ax = clp.initialize(1, 1, width=4.4, height=4.*0.618, LaTeX=False, fontsize=12)

labels = [r"%.1f Å" %(2*x) for x in dis_lst]

clp.plotone(xs, ys, ax, ylabel="ESP [a.u.]", xlim=[-2, 2], ylim=[0.25, 2], rainbowColor=True,
            xlabel=r"$x$ [Å]", labels=labels,
            showlegend=True)

ax.text(-0.14, 0.4, r"H", fontsize=12, color="r")

# Put a legend to the right of the current axis
lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), facecolor='inherit', edgecolor='inherit')

clp.adjust(tight_layout=True, savefile="figs1_esp.pdf", includelegend=lgd)
