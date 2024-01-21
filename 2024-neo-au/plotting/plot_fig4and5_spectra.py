import enum
from threading import local
import numpy as np
import columnplots as clp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.pyplot as plt

'''
Pade Transform
'''
au2fs = 0.02418884254
au2eV = 27.211399
eV2cminv = 8065.540106923572

sigma = 1e4
w_step = 2e-5
linewidth = 4.0 / sigma * au2eV
linewidth = 1.7e-3
e_cutoff_ev = 2.5
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
Compare the performance of two methods under different dt
'''

N_lst = [16, 20, 24, 26, 28, 30, 32, 34]
#N_lst = [20, 28, 34]
N_lst2 = N_lst * 2
methods = ["../data_ir_energy_transfer/N_%d_d1_0.0_excite_all/mu_n.txt" %N for N in N_lst] + ["../data_ir_energy_transfer/N_%d_d1_0.0_excite_all/mu_e.txt" %N for N in N_lst]

glob_data = {}

for filename in methods:
    # open the file for proton dipole moment
    print("working under %s" %filename)
    data = np.loadtxt(filename)
    dt = data[1,0] - data[0,0]
    nskip = int(4.0 / float(dt))
    if nskip < 1:
        nskip = 1
    t, mux, muy, muz = data[::nskip,0], data[::nskip,1], data[::nskip,2], data[::nskip,3]
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
    sptot = freq * (np.abs(spx) + np.abs(spy) + np.abs(spz))
    local_data["freq"] = freq * au2eV #* eV2cminv
    local_data["sp"] = sptot / np.max(sptot)
    glob_data[filename] = local_data

# Now we do simple plotting

colors = None
labels = [r"$N_{\rm Au}={%d}$" %N for N in N_lst]

axes = clp.initialize(2, 2, width=11.0, height=4.6, LaTeX=False, fontsize=12,
    labelthem=True, labelthemPosition=[-0.05, 1.1], labelsize=13)


xs_e, ys_e = [], []
freq_lst_e, sp_lst_e = [], []
xs_n, ys_n = [], []
freq_lst_n, sp_lst_n = [], []
for idx, filename in enumerate(methods):
    if "mu_n" in filename:
        xs_n.append(glob_data[filename]["t"])
        mun = glob_data[filename]["muz"]
        ys_n.append(mun - mun[0])
        freq_lst_n.append(glob_data[filename]["freq"])
        sp_lst_n.append(glob_data[filename]["sp"])
    elif "mu_e" in filename:
        xs_e.append(glob_data[filename]["t"])
        mun = glob_data[filename]["muz"]
        ys_e.append((mun - mun[0])) # / N_lst2[idx])
        freq_lst_e.append(glob_data[filename]["freq"])
        sp_lst_e.append(glob_data[filename]["sp"])
    # plot dipole dynamics
    #axes[idx, -1].vlines(x=0.2215 * eV2cminv, color="0.5", lw=2, ls="--", ymax=5, ymin=0, alpha=0.5)
    #axes[idx, -1].vlines(x=0.4569 * eV2cminv, color="0.5", lw=2, ls="--", ymax=5, ymin=0, alpha=0.5)

ys_e = [y + j * np.max(ys_e[0])*0 for j, y in enumerate(ys_e)]
sp_lst_e = [sp + j for j, sp in enumerate(sp_lst_e)]
ys_n = [y + j * np.max(ys_n[0])*0 for j, y in enumerate(ys_n)]
sp_lst_n = [sp + j for j, sp in enumerate(sp_lst_n)]

clp.plotone(xs_e, ys_e, axes[0, 0], colors=colors, labels=labels,
            showlegend=False, xlim=[0,30],
            lw=1.0)
# plot dipole spectrum
clp.plotone(freq_lst_e, sp_lst_e, axes[0, 1], colors=colors, labels=labels,
            showlegend=False, xlim=[0.3, 0.9],
            lw=1.0)
ys_n = [y*1e4 for y in ys_n]
clp.plotone(xs_n, ys_n, axes[1, 0], colors=colors, labels=labels,
            showlegend=False, xlim=[0,30],
            lw=1.0, sharewhichx=axes[0,0])
# plot dipole spectrum
clp.plotone(freq_lst_n, sp_lst_n, axes[1, 1], colors=colors, labels=labels,
            showlegend=False, xlim=[0.3, 0.9],
            lw=1.0)

axes[1, 0].set_xlabel("time [fs]")
axes[1, -1].set_xlabel("frequency [eV]")
# here we use mu_x instead of mu_z because in the paper the Au chain is placed along the x direction, while in our calculation it is placed along the z-direction
axes[0, 0].set_ylabel("$\mu_x^{\\rm e}$ [a.u.]")
axes[1, 0].set_ylabel(r"$\mu_x^{\rm n}$ [$\times 10^{-4}$ a.u.]")
axes[0, -1].set_ylabel("$P_{\\rm e}(\omega)$ [arb. units]")
axes[1, -1].set_ylabel("$P_{\\rm n}(\omega)$ [arb. units]")

# reset line colors
colormap = plt.cm.hot
colors = [colormap(i) for i in np.linspace(0, 0.6,len(freq_lst_e))]
for i,j in enumerate(axes[0, 0].lines):
    j.set_color(colors[i])
for i,j in enumerate(axes[0, 1].lines):
    j.set_color(colors[i])
for i,j in enumerate(axes[1, 0].lines):
    j.set_color(colors[i])
for i,j in enumerate(axes[1, 1].lines):
    j.set_color(colors[i])

# axes[0,1].legend(loc='center left', bbox_to_anchor=(1.015, 0.45), prop={'size': 10}, edgecolor="k")

axes[1, 1].axvline(x=0.6033, color=clp.navy_blue, linestyle="--", lw=1., alpha=0.6)
axes[1, 1].axvline(x=0.3358, color="0.6", linestyle="--", lw=1., alpha=0.6)
# assign arrows to indicate the plasmon frequency
pla_lst = np.array([0.8708, 0.7494, 0.6574, 0.6210, 0.5894, 0.5573, 0.5323, 0.5078])
for idx, pla_freq in enumerate(pla_lst):
    axes[0, 1].annotate("", xy=(pla_freq, float(idx) - 0.4), xytext=(pla_freq, 0.4 + float(idx)),
                        arrowprops=dict(arrowstyle="<-", color=colors[idx], lw=1.))
    axes[1, 1].annotate("", xy=(pla_freq, float(idx) - 0.4), xytext=(pla_freq, 0.4 + float(idx)),
                        arrowprops=dict(arrowstyle="<-", color=colors[idx], lw=1.))
    # add text to indicate the Au chain length
    axes[0, 1].text(0.31, sp_lst_e[idx][0]+0.2, labels[idx], color=colors[idx], fontsize=11)

clp.adjust(tight_layout=False, savefile="fig4_mu_linear_response.pdf")



# also plot the maximal protonic value
LPn_lst = np.array([0.6032, 0.6030, 0.6019, 0.5986, 0.5900, 0.5595, 0.5317, 0.5078])
UPn_lst = np.array([0.8713, 0.7500, 0.6601, 0.6248, 0.6416, 0.6025, 0.6019, 0.6019])
pla_lst = np.array([0.8708, 0.7494, 0.6574, 0.6210, 0.5894, 0.5573, 0.5323, 0.5078])
exc_lst = np.array([0.6032] * len(N_lst))

dt = xs_n[0][1] - xs_n[0][0]
nmax = int(27.0 // dt)
print("nmax = ", nmax)
osc_n_lst = np.array([np.max(y[:nmax]) for y in ys_n])

xs = [np.flip(pla_lst)]
ys = [osc_n_lst]

ax = clp.initialize(1, 1, width=4.3, height=4.3*0.618, LaTeX=False, fontsize=12)
clp.plotone(xs, ys, ax, colors=["r-o"], showlegend=False, xlabel="electronic transition [eV]",
    ylabel="Max[$\mu_x^{\\rm p}(t<27 \\rm{\ fs})$] [a.u.]", alpha=0.5)
ax.axvline(x=0.6032, color=clp.navy_blue, linestyle="--", lw=1.5, alpha=0.5)
ax.text(0.6232, 4e-4, "proton frequency", color=clp.navy_blue, alpha=0.5, fontsize=12)
clp.adjust(tight_layout=False, savefile="fig5_maximal.tiff")


