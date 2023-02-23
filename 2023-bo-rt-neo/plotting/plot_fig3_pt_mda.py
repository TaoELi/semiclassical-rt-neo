import enum
from threading import local
import numpy as np
import columnplots as clp
import matplotlib.pyplot as plt

au2fs = 0.02418884254
Angstrom2AU = 1.8897259885789

def get_distance(data):
    OD = data[:, 0:3]
    OA = data[:, 3:6]
    H = data[:, -3:]
    # calculate HOD, HOA distance
    R_HOD = np.sqrt(np.sum((H - OD)**2, axis=1)) 
    R_HOA = np.sqrt(np.sum((H - OA)**2, axis=1)) 
    return R_HOD, R_HOA

methods = [ "../proton_transfer/neo_eh_3cen/",
            "../proton_transfer/neo_boeh_3cen/", 
            "../proton_transfer/eh_cl/", 
            "../proton_transfer/boeh_cl/"]

# coeffs used to rescale the time step
# the prefactor 4.2/4.0 is because all BO simulations output molecular geometry every 4.2 a.u., 
# while the non-BO simulations output molecular geometry every 4.0 a.u.
coeffs = [1.0, 4.2/4.0, 1.0, 4.2/4.0]

glob_data = {}

for idx, dir in enumerate(methods):
    # open the file for molecular geometry
    filename = "%s/config.txt" %(dir)
    print("working under %s" %filename)
    data = np.loadtxt(filename)
    # create a dictionary to store important info
    local_data = {}
    local_data["t"] = np.arange(np.shape(data)[0]) * 4.0 * au2fs * coeffs[idx]
    R_HOD, R_HOA = get_distance(data)
    local_data["R_HOD"] = R_HOD
    local_data["R_HOA"] = R_HOA
    # do an additional truncation to keep the initial 100 fs
    nidx = np.sum(local_data["t"] < 100)
    local_data["t"] = local_data["t"][:nidx]
    local_data["R_HOD"] = local_data["R_HOD"][:nidx]
    local_data["R_HOA"] = local_data["R_HOA"][:nidx]
    # calculate distance to OA and OD
    glob_data[filename] = local_data

# Now we do simple plotting

colors = [clp.black, clp.red, clp.sky_blue]
linestyles = ["-", "-",  "--", "--"]
labels_eh = [r"RT-NEO-Ehrenfest", r"BO-RT-NEO-Ehrenfest"]
labels_cl = [r"Ehrenfest classical proton", r"BOMD classical proton"]

figure, axes = clp.initialize(2, 1, width=4.3, height=4.3*0.618*2.4, LaTeX=True, fontsize=12, 
    labelsize=14, labelthem=True, labelthemPosition=[-0.02, 1.02], sharey=True, return_fig_args=True)


xs, y1s, y2s = [], [], []
for idx2, dir in enumerate(methods[0:2]):
    filename = "%s/config.txt" %(dir)
    xs.append(glob_data[filename]["t"])
    y1s.append(glob_data[filename]["R_HOD"])
    y2s.append(glob_data[filename]["R_HOA"])
clp.plotone(xs, y1s, axes[0], colors=colors, labels=labels_eh, lw=1.5,
    showlegend=True if idx == 0 else False,
    xlim=[0,100], ylim=[0.935, 1.735],
    ylabel="H-O distance [Angstrom]", alphaspacing=0.2)
clp.plotone(xs, y2s, axes[0], colors=colors, linestyles=linestyles,  labels=labels_eh, lw=1.5, 
    showlegend=False, alphaspacing=0.2)

xs, y1s, y2s = [], [], []
for idx2, dir in enumerate(methods[2:]):
    filename = "%s/config.txt" %(dir)
    xs.append(glob_data[filename]["t"])
    y1s.append(glob_data[filename]["R_HOD"])
    y2s.append(glob_data[filename]["R_HOA"])
clp.plotone(xs, y1s, axes[1], colors=colors, labels=labels_cl, lw=1.5,
    showlegend=True if idx == 0 else False,
    xlabel="time [fs]",
    xlim=[0,100],
    ylabel="H-O distance [Angstrom]", alphaspacing=0.2)
clp.plotone(xs, y2s, axes[1], colors=colors, linestyles=linestyles, labels=labels_cl, lw=1.5, 
    showlegend=False, alphaspacing=0.2)

# finally, we insert the image
import matplotlib.image as mpimg
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
import matplotlib.patches as mpatches

data_img = mpimg.imread('mda_3cen_0.png')
imagebox = OffsetImage(data_img, zoom=0.1)
ab = AnnotationBbox(imagebox, (10, 1.24), frameon=False, annotation_clip=False)
ab.set(zorder=-1)
axes[0].add_artist(ab)

clp.adjust(tight_layout=False, savefile="fig3_pt_eh.pdf")
