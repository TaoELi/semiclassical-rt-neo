import enum
from threading import local
from matplotlib import transforms
import numpy as np
import columnplots as clp
import matplotlib.pyplot as plt

'''
Pade Transform
'''
au2fs = 0.02418884254
au2eV = 27.211399
eV2cminv = 8065.540106923572

sigma = 1e6
w_step = 2e-6
linewidth = 4.0 / sigma * au2eV
linewidth = 1.7e-3
e_cutoff_ev = 3.0
e_cutoff_au = e_cutoff_ev / au2eV
e_start_ev = 0.0
e_start_au = e_start_ev / au2eV

fsinv2eV = 4.135668 # 1fs-1 to 4.13 eV
fsinv2cminv = fsinv2eV * eV2cminv

def pade(time,signal,sigma=sigma,max_len=None,w_min=e_start_au,w_max=e_cutoff_au,w_step=w_step,read_freq=None):
    """ Routine to take the Fourier transform of a time signal using the method
          of Pade approximants.
        Inputs:
          time:      (list or Numpy NDArray) signal sampling times
          signal:    (list or Numpy NDArray)
        Optional Inputs:
          sigma:     (float) signal damp factor, yields peaks with
                       FWHM of 2/sigma
          max_len:   (int) maximum number of points to use in Fourier transform
          w_min:     (float) lower returned frequency bound
          w_max:     (float) upper returned frequency bound
          w_step:    (float) returned frequency bin width
        Returns:
          fsignal:   (complex NDArray) transformed signal
          frequency: (NDArray) transformed signal frequencies
        From: Bruner, Adam, Daniel LaMaster, and Kenneth Lopata. "Accelerated
          broadband spectra using transition signal decomposition and Pade
          approximants." Journal of chemical theory and computation 12.8
          (2016): 3741-3750.
    """

    # center signal about zero
    signal = np.asarray(signal) - signal[0]

    stepsize = time[1] - time[0]

    # Damp the signal with an exponential decay.
    damp = np.exp(-(stepsize*np.arange(len(signal)))/float(sigma))
    signal *= damp

    M = len(signal)
    N = int(np.floor(M / 2))

    # Check signal length, and truncate if too long
    if max_len:
        if M > max_len:
            N = int(np.floor(max_len / 2))

    # G and d are (N-1) x (N-1)
    # d[k] = -signal[N+k] for k in range(1,N)
    d = -signal[N+1:2*N]

    try:
        from scipy.linalg import toeplitz, solve_toeplitz
        # Instead, form G = (c,r) as toeplitz
        #c = signal[N:2*N-1]
        #r = np.hstack((signal[1],signal[N-1:1:-1]))
        b = solve_toeplitz((signal[N:2*N-1],\
            np.hstack((signal[1],signal[N-1:1:-1]))),d,check_finite=False)
    except (ImportError,np.linalg.linalg.LinAlgError) as e:
        # OLD CODE: sometimes more stable
        # G[k,m] = signal[N - m + k] for m,k in range(1,N)
        G = signal[N + np.arange(1,N)[:,None] - np.arange(1,N)]
        b = np.linalg.solve(G,d)

    # Now make b Nx1 where b0 = 1
    b = np.hstack((1,b))

    # b[m]*signal[k-m] for k in range(0,N), for m in range(k)
    a = np.dot(np.tril(toeplitz(signal[0:N])),b)
    p = np.poly1d(np.flip(a))
    q = np.poly1d(np.flip(b))

    if read_freq is None:
        # choose frequencies to evaluate over
        frequency = np.arange(w_min,w_max,w_step)
    else:
        frequency = read_freq

    W = np.exp(-1j*frequency*stepsize)

    fsignal = p(W)/q(W)

    return fsignal, frequency


'''
Compare the performance of three methods under different dt
'''

methods = ["../benchmark_incav/%s" %x for x in ["rt_delta_proton", "rtBO_delta_proton"]]
dt_lst = ["0.20", "0.40", "1.0", "8.0"]

glob_data = {}

for dir in methods:
    for dt in dt_lst:
        # open the file for proton dipole moment
        filename = "%s/dt_%s/mu_n.txt" %(dir, dt)
        print("working under %s" %filename)
        data = np.loadtxt(filename)
        nskip = int(4.0 / float(dt))
        if nskip < 1:
            nskip = 1
        t, mux, muy, muz = data[::nskip,1], data[::nskip,2], data[::nskip,3], data[::nskip,4]
        # create a dictionary to store important info
        local_data = {}
        local_data["t"] = t * au2fs
        local_data["mux"] = mux
        local_data["muy"] = muy
        local_data["muz"] = muz
        # calculate the spectroscopy from Pade approximation
        spx, freq = pade(t, mux)
        spy, freq = pade(t, muy)
        spz, freq = pade(t, muz)
        sptot = np.abs(spx) + np.abs(spy) + np.abs(spz)
        local_data["freq"] = freq * au2eV * eV2cminv
        local_data["sp"] = sptot / np.max(sptot)
        glob_data[filename] = local_data

# Now we do simple plotting

colors = [clp.black, clp.red, clp.dark_green, clp.sky_blue]
labels = ["%.3f fs" %(float(dt)*au2fs) for dt in dt_lst]

axes = clp.initialize(2, 2, width=9.0, height=4.6, LaTeX=True, fontsize=12, 
    labelthem=True, labelthemPosition=[-0.05, 1.1], labelsize=13)


for idx, dir in enumerate(methods):
    xs, ys = [], []
    freq_lst, sp_lst = [], []
    for dt in dt_lst:
        filename = "%s/dt_%s/mu_n.txt" %(dir, dt)
        xs.append(glob_data[filename]["t"])
        mun = glob_data[filename]["mux"]
        ys.append(mun - mun[0])
        freq_lst.append(glob_data[filename]["freq"])
        sp_lst.append(glob_data[filename]["sp"])
    # add vertical lines for comparision
    axes[idx, 1].axvline(x=0.2215 * eV2cminv, c="0.5", lw=4.0, ls="-",  alpha=0.5)
    axes[idx, 1].axvline(x=0.4569 * eV2cminv, c="0.5", lw=4.0, ls="-",  alpha=0.5)
    
    ys = [y * 1e3 + j*10 for j, y in enumerate(ys)]
    sp_lst = [sp + j*0.5 for j, sp in enumerate(sp_lst)]

    clp.plotone(xs, ys, axes[idx, 0], colors=colors, labels=labels, 
                ylim=[-5.5, 43], 
                showlegend=False,
                lw=1.0, 
                xlim=[-10, 1e3]
                )
    
     # plot dipole spectrum
    clp.plotone(freq_lst, sp_lst, axes[idx, 1], colors=colors, labels=labels, 
                xlim=[0.4569 * eV2cminv-100, 0.4569 * eV2cminv+100], ylim=[-0.001, 2.7],
                showlegend=True if idx == 0 else False,
                lw=1.0)

axes[0,1].legend(#loc='center left', #bbox_to_anchor=(1, 0.45), 
                prop={'size': 10}, edgecolor="k")

axes[1, 0].set_xlabel("time [fs]")
axes[1, 1].set_xlabel("frequency [cm$^{-1}$]")
axes[0, 0].set_ylabel(r"$\mu_x^{\rm n}$ [$\times 10^{-3}$ a.u.]")
axes[1, 0].set_ylabel(r"$\mu_x^{\rm n}$ [$\times 10^{-3}$ a.u.]")
axes[0, 1].set_ylabel("$P_{\\rm n}(\omega)$ [arb. units]")
axes[1, 1].set_ylabel("$P_{\\rm n}(\omega)$ [arb. units]")
axes[0, 0].text(0.04, 0.85, "semiclassical RT-NEO for polaritons", transform=axes[0,0].transAxes, fontsize=12)
axes[1, 0].text(0.04, 0.85, "semiclassical BO-RT-NEO for polaritons", transform=axes[1,0].transAxes, fontsize=12)

# finally, we insert the image
import matplotlib.image as mpimg
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox

data_img = mpimg.imread('HCN_demo_incav.png')
imagebox = OffsetImage(data_img, zoom=0.08)
ab = AnnotationBbox(imagebox, (730, 12), frameon=False)
ab.set(zorder=-1)
axes[0,0].add_artist(ab)

clp.adjust(tight_layout=False, savefile="fig2_mun_incav.pdf")
