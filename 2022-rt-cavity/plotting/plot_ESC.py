'''
This script contains the following two figures

(a) Plot HCN molecule under electronic strong coupling
1. Outcav muz(t) 2. Outcav electronic spectrum
3. Incav muz(t) 4. Incav electronic spectrum

(b) Rabi splitting dependence and cavity loss dependence

(c), (d) Same plots for VSC
'''

import numpy as np
import sys
import columnplots as clp
import matplotlib.cm as cm
import matplotlib.pyplot as plt

au2fs = 0.02418884254
au2eV = 27.211399
eV2cminv = 8065.540106923572

sigma = 1e5
w_step = 1e-5
linewidth = 4.0 / sigma * au2eV
linewidth = 1.7e-3
print("linewidth is", linewidth)
e_cutoff_ev = 20.0
e_cutoff_au = e_cutoff_ev / au2eV
e_start_ev = 10.0
e_start_au = e_start_ev / au2eV

fsinv2eV = 4.135668 # 1fs-1 to 4.13 eV

nskip = 10

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

def lorentz(x, w0, gamma, I):
    return I / np.pi * (0.5*gamma) / ((x-w0)**2 + (0.5*gamma)**2)

def gaussian(x, w0, gamma, I):
    return I / (2.0 * np.pi * gamma**2)**(0.5) * np.exp(-(x-w0)**2/2.0/gamma**2)

def get_LR_spectrum(filename):
    data = np.loadtxt(filename)
    energy, e_osc, p_osc, ep_osc = data[:,0], data[:,1], data[:,2], data[:,3]
    freq_ev = np.linspace(0, e_cutoff_ev, int(e_cutoff_ev / w_step)+1)
    n = np.size(freq_ev)
    e_sp, p_sp, ep_sp = np.zeros(n), np.zeros(n), np.zeros(n)
    for idx in range(len(energy)):
        w0 = energy[idx]
        I_e = e_osc[idx]
        I_p = p_osc[idx]
        I_ep = ep_osc[idx]
        e_sp += lorentz(freq_ev, w0, linewidth, I_e)
        p_sp += lorentz(freq_ev, w0, linewidth, I_p)
        ep_sp += lorentz(freq_ev, w0, linewidth, I_ep)
    return freq_ev, e_sp, p_sp, ep_sp

def get_dipole(filename):
    data = np.loadtxt(filename)
    nst = int(np.shape(data)[0] * 0.)
    nend = int(np.shape(data)[0] * 1.0)
    t, mux, muy, muz = data[nst:nend:nskip,0], data[nst:nend:nskip,1], data[nst:nend:nskip,2], data[nst:nend:nskip,3]
    t_fs = t * au2fs
    # we truncate t_fs end to 25 fs
    '''
    if t_fs[-1] > 100:
        print("exceeding 25 fs for %s and will truncate" %filename)
        df_fs = t_fs[1] - t_fs[0]
        nend = int(25.0 / df_fs)
        t = t[0:nend]
        mux = mux[0:nend]
        muy = muy[0:nend]
        muz = muz[0:nend]
        t_fs = t_fs[0:nend]
    else:
        print("Less than 25 fs")
    '''
    spx, freq = pade(t, mux)
    spy, freq = pade(t, muy)
    spz, freq = pade(t, muz)

    #mu_norm = np.sqrt(mux**2 + muy**2 + muz**2)
    mu_norm = muz
    freq_ev = freq * au2eV
    sp_tot = np.abs(spx) + np.abs(spy) + np.abs(spz)
    #sp_tot = freq * np.abs(spz)
    df = freq_ev[1] - freq_ev[0]
    # output only selectively number of data
    nmax = int(e_cutoff_ev // df)
    freq_ev = freq_ev[0:nmax]
    sp_tot = sp_tot[0:nmax]
    return t_fs, mu_norm, freq_ev, sp_tot

# get linear response data
freq_lr, e_sp_lr, p_sp_lr, ep_sp_lr = get_LR_spectrum(filename="../HCN/tddft_outcav/excitations_LR.txt")

def plot_dynamics():
    # get real-time data
    path_outcav = "../HCN/tddft_outcav/"
    t_outcav, mu_outcav, freq_outcav, sp_outcav = get_dipole(filename=path_outcav + "/mu_e.txt")
    path_incav = "../HCN/tddft_ESC/Freq_13.3344764488_CavLifetime_1e8_Coupling_4e-3/"
    t_incav, mu_incav, freq_incav, sp_incav = get_dipole(filename=path_incav + "/mu_e.txt")

    fig, axes = clp.initialize(2, 2, width=8, height=4, return_fig_args=True,
        fontsize=12, LaTeX=True, labelthem=True, labelthemPosition=[0.1, 0.963], labelsize=13)

    labels = ["RT-NEO-TDDFT", "LR-NEO-TDDFT"]
    colors = [clp.navy_blue, "k--"]

    x1s, y1s = [t_outcav], [mu_outcav]
    clp.plotone(x1s, y1s, axes[0, 0], colors=colors, labels=labels, lw=0.5,
             ylabel=r"$\mu_z^{\rm e}(t)$ [a.u.]", showlegend=False, ylim=[-6e-3, 6e-3], xlim=[0, 50],  yscientific=True)

    x2s, y2s = [freq_outcav, freq_lr], [sp_outcav / np.max(sp_outcav), e_sp_lr / np.max(e_sp_lr)]
    clp.plotone(x2s[0:1], y2s[0:1], axes[0, 1], colors=colors, labels=labels, lw=1,
             ylabel=r"$P_{\rm e}(\omega)$ [arb. units]", xlim=[12, 18], alphaspacing=0.2)
    clp.plotone(x2s[1:], y2s[1:], axes[0, 1], colors=colors[1:], labels=labels[1:], lw=0.5,
             ylabel=r"$P_{\rm e}(\omega)$ [arb. units]", xlim=[12, 18], alphaspacing=0.2)

    labels = ["sc RT-NEO-TDDFT"]
    colors = [clp.red]

    x3s, y3s = [t_incav], [mu_incav]
    clp.plotone(x3s, y3s, axes[1, 0], colors=colors, labels=labels, lw=0.5,
            xlabel="time [fs]", ylabel=r"$\mu_z^{\rm e}(t)$ [a.u.]", showlegend=False, ylim=[-3e-3, 3e-3], xlim=[0, 50], yscientific=True)

    x4s, y4s = [freq_incav], [sp_incav / np.max(sp_incav)]
    clp.plotone(x4s, y4s, axes[1, 1], colors=colors, labels=labels, lw=1,
            xlabel="energy [eV]", ylabel=r"$P_{\rm e}(\omega)$ [arb. units]", xlim=[12, 18], showlegend=True)

    axes[1,1].axvline(x=13.3344764488, c=clp.navy_blue, lw=1.0, ls="--")

    axes[0,0].text(0.75, 0.85, "cavity off", transform=axes[0,0].transAxes, fontsize=14, c=clp.navy_blue, fontweight="heavy")
    axes[1,0].text(0.75, 0.85, "cavity on", transform=axes[1,0].transAxes, fontsize=14, c=clp.red, fontweight="heavy")

    # finally, we insert some images
    import matplotlib.image as mpimg
    from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox

    data_img = mpimg.imread('hcn_illustration.png')
    imagebox = OffsetImage(data_img, zoom=0.1)
    ab = AnnotationBbox(imagebox, (16.5, 0.34), frameon=False)
    ab.set(zorder=-1)
    axes[1,1].add_artist(ab)

    clp.adjust(tight_layout=True, savefile="ESC_dynamics.pdf")

def prepare_Rabi_versus_coupling():
    coupling_lst = ["1e-3", "2e-3", "3e-3", "4e-3", "5e-3", "6e-3"]
    freq_lst, sp_lst = [], []
    for coupling in coupling_lst:
        path_incav = "../HCN/tddft_ESC/Freq_13.3344764488_CavLifetime_1e2_Coupling_%s/" %coupling
        t_incav, mu_incav, freq_incav, sp_incav = get_dipole(filename=path_incav + "/mu_e.txt")
        # truncate to a small range and normalize the spectrum
        freq_start = 12 - e_start_ev
        freq_end = 14.5 - e_start_ev
        df = freq_incav[1] - freq_incav[0]
        n_start = int(freq_start/df)
        n_end = int(freq_end/df)
        print("dealing with %s" %path_incav)
        freq_incav = freq_incav[n_start:n_end]
        sp_incav = sp_incav[n_start:n_end]
        sp_incav /= np.max(sp_incav)
        freq_lst.append(freq_incav)
        sp_lst.append(sp_incav)
    # label
    labels = ["$\\varepsilon = %d \\times 10^{-3}$ a.u." %(n+1) for n in range(len(coupling_lst))]
    return freq_lst, sp_lst, labels

def prepare_Rabi_versus_loss():
    loss_lst = ["50", "1e2", "5e2", "1e3"]
    loss_lst = loss_lst[::-1]
    freq_lst, sp_lst = [], []
    for loss in loss_lst:
        path_incav = "../HCN/tddft_ESC/Freq_13.3344764488_CavLifetime_%s_Coupling_4e-3/" %loss
        t_incav, mu_incav, freq_incav, sp_incav = get_dipole(filename=path_incav + "/mu_e.txt")
        # truncate to a small range and normalize the spectrum
        freq_start = 12 - e_start_ev
        freq_end = 14.5 - e_start_ev
        df = freq_incav[1] - freq_incav[0]
        n_start = int(freq_start/df)
        n_end = int(freq_end/df)
        print("dealing with %s" %path_incav)
        freq_incav = freq_incav[n_start:n_end]
        sp_incav = sp_incav[n_start:n_end]
        sp_incav /= np.max(sp_incav)
        freq_lst.append(freq_incav)
        sp_lst.append(sp_incav)
    # label
    loss_lst_fs = [float(x) * au2fs for x in loss_lst]
    #labels = ["$\\tau_c = $ %.1f fs" %x for x in loss_lst_fs]
    labels = ["$\\gamma_c = $ %.2f eV" %(1.0/x*fsinv2eV) for x in loss_lst_fs]
    return freq_lst, sp_lst, labels

def plot_dependence():
    fig, axes = clp.initialize(1, 2, width=8.6, height=4.3*0.618, return_fig_args=True,
        fontsize=12, LaTeX=True, labelthem=True, labelthemPosition=[0.15, 0.963])

    # plot Rabi versus coupling strength
    freq_lst, sp_lst, labels = prepare_Rabi_versus_coupling()

    clp.plotone(freq_lst, sp_lst, axes[0], labels=labels, xlim=[12, 14.5], ylim=[0,1.05], lw=0.8,
         xlabel="energy [eV]", ylabel=r"$P_{\rm e}(\omega)$ [arb. units]")

    # set the colormap
    colormap = plt.cm.hot
    colors = [colormap(i) for i in np.linspace(0, 0.6,len(freq_lst))]
    for i,j in enumerate(axes[0].lines):
        j.set_color(colors[i])

    axes[0].legend(loc='center left', bbox_to_anchor=(1, 0.45), prop={'size': 10}, edgecolor="k")

    axes[0].axvline(x=13.3344764488, c=clp.navy_blue, lw=1.0, ls="--")


    # plot Rabi versus cavity loss
    freq_lst, sp_lst, labels = prepare_Rabi_versus_loss()

    clp.plotone(freq_lst, sp_lst, axes[1], labels=labels, xlim=[12, 14.5], ylim=[0,1.05], lw=0.8,
        xlabel="energy [eV]")#, ylabel=r"$P_{\rm e}(\omega)$ [arb. units]")

    # set the colormap
    colormap = plt.cm.cool
    colors = [colormap(i) for i in np.linspace(0, 0.8,len(freq_lst))]
    for i,j in enumerate(axes[1].lines):
        j.set_color(colors[i])

    axes[1].legend(loc='center left', bbox_to_anchor=(1, 0.45), prop={'size': 10}, edgecolor="k")

    axes[1].axvline(x=13.3344764488, c=clp.navy_blue, lw=1.0, ls="--")

    gamma_c = 1.0/(1e2 * au2fs) * fsinv2eV
    axes[0].text(1.17, 0.82, r"$\gamma_{\rm c} = %.2f$ eV" %gamma_c, transform=axes[0].transAxes, fontsize=12)
    axes[1].text(1.04, 0.72, "$\\varepsilon = 4\\times 10^{-3}$ a.u.", transform=axes[1].transAxes, fontsize=12)

    clp.adjust(tight_layout=True, savefile="ESC_dependence.pdf")

if __name__ == '__main__':
    plot_dynamics()
    plot_dependence()
