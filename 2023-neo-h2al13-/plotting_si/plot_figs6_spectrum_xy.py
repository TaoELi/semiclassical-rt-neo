import enum
from threading import local
import numpy as np
import columnplots as clp
import matplotlib.pyplot as plt

'''
Pade Transform
'''
au2fs = 0.02418884254
au2eV = 27.211399
eV2cminv = 8065.540106923572

sigma = 1e4
w_step = 1e-4
linewidth = 4.0 / sigma * au2eV
linewidth = 1.7e-3
e_cutoff_ev = 20.0
e_cutoff_au = e_cutoff_ev / au2eV
e_start_ev = 2.0
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


methods = ["%s" %x for x in ["../run_pbe_6-31gd/eq_rt_delta_long", "../run_pbe_6-31gd/eq_neo_delta_long"]]

dt = 0.1

glob_data = {}

nmax = 20582
nskip = 2

for dir in methods:
    # open the file for proton dipole moment
    filename = "%s/xe.txt" %(dir)
    print("working under %s" %filename)
    data = np.loadtxt(filename)
    # create a dictionary to store important info
    local_data = {}
    print("size is", np.size(data[:,0]))
    local_data["t"] = data[:nmax:nskip,0] * dt * au2fs
    local_data["mux"] = data[:nmax:nskip,1]
    local_data["muy"] = data[:nmax:nskip,2]
    local_data["muz"] = data[:nmax:nskip,3]
    # calculate the spectroscopy from Pade approximation
    t, mux, muy, muz = data[:nmax:nskip,0] * dt, data[:nmax:nskip,1], data[:nmax:nskip,2], data[:nmax:nskip,3]
    spx, freq = pade(t, mux)
    spy, freq = pade(t, muy)
    spz, freq = pade(t, muz)
    sptot = (np.abs(spx) + np.abs(spy) + np.abs(spz))
    local_data["freq"] = freq * au2eV
    local_data["sp"] = sptot / np.max(sptot)
    glob_data[dir] = local_data

# Now we do simple plotting

colors = [clp.black, clp.red]
linestyles = ["-", "-"]
labels = ["conventional RT-TDDFT", "RT-NEO-TDDFT"]

fig, axes = clp.initialize(2, 2, width=7., height=3.5*0.618*2, LaTeX=False, fontsize=10,
                           labelthem=True, labelthemPosition=[0.10, 0.95], labelsize=11,
                           return_fig_args=True)

xs, ys = [], []
freq_lst, sp_lst = [], []

signal_x, signal_y, signal_z = [], [], []

for idx, dir in enumerate(methods):
    xs.append(glob_data[dir]["t"])
    muex = glob_data[dir]["mux"]
    muey = glob_data[dir]["muy"]
    muez = glob_data[dir]["muz"]
    ys.append(muex - muex[0])
    signal_x.append(muex - muex[0])
    signal_y.append(muey - muey[0])
    signal_z.append(muez - muez[0])
    freq_lst.append(glob_data[dir]["freq"])
    sp_lst.append(glob_data[dir]["sp"])

# plot dipole dynamics
clp.plotone(xs, signal_x, axes[0,0], colors=colors, labels=labels, ylim=[-0.015, 0.015], 
            xlim=[0, 50], xlabel="time [fs]", ylabel="$\mu_x^{\\rm e}$ [a.u.]",
            showlegend=True, alphaspacing=0.2, linestyles=linestyles, legendFontSize=8.,
            lws=[1., 0.8])
clp.plotone(xs, signal_y, axes[0,1], colors=colors, labels=labels, ylim=[-0.015, 0.015], 
            xlim=[0, 50], xlabel="time [fs]", ylabel="$\mu_y^{\\rm e}$ [a.u.]",
            showlegend=False, alphaspacing=0.2, linestyles=linestyles,
            lws=[1., 0.8])
clp.plotone(xs, signal_z, axes[1,0], colors=colors, labels=labels, ylim=[-0.015, 0.015], 
            xlim=[0, 50], xlabel="time [fs]", ylabel="$\mu_z^{\\rm e}$ [a.u.]",
            showlegend=False, alphaspacing=0.2, linestyles=linestyles,
            lws=[1., 0.8])
# plot dipole spectrum
clp.plotone(freq_lst, sp_lst, axes[1,1], colors=colors, labels=labels, ylim=[-0.04, 1.05],
            xlim=[6, 7], xlabel="frequency [eV]", ylabel="$P_{\\rm e}(\omega)$ [arb. units]",
            linestyles=linestyles, showlegend=False,
            lws=[1., 0.8],
            alphaspacing=0.2,)

# finally, we insert some images
#clp.add_figure(imag_filename="Al13H2_eq.png", ax=axes[1], zoom=0.15, location=(13.5, 0.4))

#clp.subplots_adjust()

clp.adjust(savefile="figs6_spectrum_full.pdf", tight_layout=True)
