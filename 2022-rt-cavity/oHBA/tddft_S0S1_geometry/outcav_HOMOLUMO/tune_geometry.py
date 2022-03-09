import numpy as np

file_rgs = "oHBA_rgs_raw.xyz"
file_grd = "oHBA_ground_raw.xyz"

data_rgs = np.loadtxt(file_rgs)
data_grd = np.loadtxt(file_grd)

Ntot = 20
for i in range(Ntot+1):
    ratio = i / Ntot
    print("Create a geometry with %.2f ground and %.2f rgs" %(ratio, 1 - ratio))
    data_new = ratio * data_grd + (1-ratio) * data_rgs
    filename = "oHBA_%d.txt" %i
    np.savetxt(filename, data_new, fmt="%.6f")
