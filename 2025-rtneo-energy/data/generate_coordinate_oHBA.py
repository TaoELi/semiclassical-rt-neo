import numpy as np
import sys

AU2Angstrom = 1.0 / 1.8897259886

pathname = sys.argv[1]

data_x = np.loadtxt(pathname + "/x.txt")
data_y = np.loadtxt(pathname + "/y.txt")
data_z = np.loadtxt(pathname + "/z.txt")
data_xn = np.loadtxt(pathname + "/xn.txt")

nframes = data_x.shape[0]
natoms = data_x.shape[1]-1
print(nframes, natoms)
atom_string = ["C"]*6 + ["O", "C", "O"] + ["H"]*9

with open(pathname +"/config.xyz", "w") as f:
    for i in range(nframes):
        print("generating frame No.", i)
        f.write("%d\n\n" %natoms)
        xs = data_x[i,1:] * AU2Angstrom
        ys = data_y[i,1:] * AU2Angstrom
        zs = data_z[i,1:] * AU2Angstrom
        # replace the last proton from basis center to expectation value
        data_xn_idx = data_xn[:,0]
        idx_xn = np.where(data_xn_idx == data_x[i,0])[0][0]
        print("idx is", idx_xn)
        r_proton = data_xn[idx_xn, 1:]
        xs[-1] = r_proton[0] * AU2Angstrom
        ys[-1] = r_proton[1] * AU2Angstrom
        zs[-1] = r_proton[2] * AU2Angstrom
        for idx in range(natoms):
            f.write("%s %.6f %.6f %.6f\n" %(atom_string[idx], xs[idx], ys[idx], zs[idx]))

with open(pathname +"/config.txt", "w") as f:
    for i in range(nframes):
        print("generating frame No.", i)
        xs = data_x[i,1:] * AU2Angstrom
        ys = data_y[i,1:] * AU2Angstrom
        zs = data_z[i,1:] * AU2Angstrom
        # replace the last proton from basis center to expectation value
        data_xn_idx = data_xn[:,0]
        idx_xn = np.where(data_xn_idx == data_x[i,0])[0][0]
        print("idx is", idx_xn)
        r_proton = data_xn[idx_xn, 1:]
        xs[-1] = r_proton[0] * AU2Angstrom
        ys[-1] = r_proton[1] * AU2Angstrom
        zs[-1] = r_proton[2] * AU2Angstrom
        for idx in range(natoms):
            f.write("%.6f %.6f %.6f " %(xs[idx], ys[idx], zs[idx]))
        f.write("\n")
