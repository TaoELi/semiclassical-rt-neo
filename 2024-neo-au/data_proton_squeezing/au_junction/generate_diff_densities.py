import numpy as np
import sys

def load_density(filename, skiprows=19):
    # load header
    with open(filename) as f:
        header = f.readlines()[0:skiprows]
    # load data
    data = np.loadtxt(filename, skiprows=skiprows)
    return header, data


def save_density(filename, header, data):
    np.savetxt(filename, data,  fmt='%.6e')
    with open(filename, 'r') as f:
        content = f.read()
    with open(filename, 'w') as f:
        f.writelines(header)
        f.write(content)


def get_modified_density(filename_init, filename_final, filename_new):
    header_init, data_init = load_density(filename_init)
    header_final, data_final = load_density(filename_final)
    data_diff = (data_init - data_final) * 1e3
    # replace atomic number to AU
    header_init = [x.replace("  19", "  79") for x in header_init]
    save_density(filename=filename_new, header=header_init, data=data_diff)


def get_diff_densities(dir, n_files=1):
    # get difference density for electronic components
    for i in range(n_files):
        print("Processing file {} for protons".format(i))
        nstep = 400 * i
        filename_init = dir + "/den_p_{}_0.cube".format(nstep)
        filename_final = "d_1.675/den_p_{}_0.cube".format(nstep)
        filename_new = dir + "/diff_den_p_{}_0.cube".format(nstep)
        get_modified_density(filename_init, filename_final, filename_new)

if __name__ == "__main__":
    dir = sys.argv[-1]
    print("Working under path {}".format(dir))
    get_diff_densities(dir)
