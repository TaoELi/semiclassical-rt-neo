import numpy as np
import columnplots as clp


neo_mv = {'total_time': 310.23, 'scf_time': 12.0/60.0, 'neo_scf_time': 360.75/60.0}
cl_mv = {'total_time': 236.24, 'scf_time': 12.0/60.0, 'neo_scf_time': 24.8/60.0}
cl_fix = {'total_time': 113.52, 'scf_time': 13.0/60.0, 'neo_scf_time': 25.91/60.0}
neo_fix = {'total_time': 145.42, 'scf_time': 12.0/60.0, 'neo_scf_time': 373.03/60.0}
neo_2c = {'total_time': 228.34, 'scf_time': 48.0/60.0, 'neo_scf_time': 1326.32/60.0}
neo_3c = {'total_time': 332.16, 'scf_time': 57.0/60.0, 'neo_scf_time': 2814.73/60.0}
data = {"neo_mv": neo_mv, "cl_mv": cl_mv, "cl_fix": cl_fix, "neo_fix": neo_fix, "neo_2c": neo_2c, "neo_3c": neo_3c}

for key, value in data.items():
    value['dynamics_time'] = value['total_time'] - value['scf_time'] - value['neo_scf_time']

# calculate time cost per time step
for key, value in data.items():
    value['dynamics_time'] /= 1000.0
    print(key, value['dynamics_time'])

ax = clp.initialize(1, 1, width=3.25, height=3.25*0.618*1.3, fontsize=8, LaTeX=False)

clp.plotone([np.zeros(0)], [np.zeros(0)], ax=ax, xlim=[-0.5, 5.5], showlegend=False,
            ylabel="computational cost [min/step]")
# add some colors to the background
ax.axvspan(-0.5, 1.5, alpha=0.2, color='0.5')
ax.axvspan(1.5, 3.5, alpha=0.2, color='blue')
ax.axvspan(3.5, 5.5, alpha=0.2, color='green')
ax.text(-0.2, 0.29, "Fig. 1", fontsize=8, color='0.5', alpha=0.6)
ax.text(1.8, 0.29, "Fig. 2", fontsize=8, color='blue', alpha=0.6)
ax.text(3.8, 0.29, "Fig. S2", fontsize=8, color='green', alpha=0.6)

ax.bar("RT-TDDFT", data["cl_fix"]["dynamics_time"], color=clp.black, edgecolor='k', linewidth=1, width=0.3,
       fill=False, hatch='//')
ax.bar("RT-NEO FPB/1c", data["neo_fix"]["dynamics_time"], color=clp.red, edgecolor=clp.red, linewidth=1, width=0.3,
       fill=False, hatch='//')
ax.bar("RT-Ehrenfest", data["cl_mv"]["dynamics_time"], color=clp.black, edgecolor='k', linewidth=1, width=0.3)
ax.bar("RT-NEO TPB/1c", data["neo_mv"]["dynamics_time"], color=clp.red, edgecolor='k', linewidth=1, width=0.3)
ax.bar("RT-NEO FPB/2c", data["neo_2c"]["dynamics_time"], color=clp.sky_blue, edgecolor='k', linewidth=1, width=0.3)
ax.bar("RT-NEO FPB/3c", data["neo_3c"]["dynamics_time"], color=clp.dark_green, edgecolor='k', linewidth=1, width=0.3)

import matplotlib.pyplot as plt
plt.xticks(rotation=45, ha='right')

clp.adjust(tight_layout=True, savefile="figs2_time_compare.pdf")