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


ax = clp.initialize(1, 1, width=3.25, LaTeX=False, fontsize=10, labelsize=11)


# add pulse spectrum
t = np.linspace(0, 1000, 1000)
t0 = 0
tau = 800.0 
omega = 6.0 / au2eV
Ex = np.exp(-(t/tau)**2) * np.sin(omega * t)
sp_pulse, freq_pulse = pade(t, Ex)
sp_pulse = np.abs(sp_pulse)
freq_pulse *= au2eV
sp_pulse /= np.max(sp_pulse)
clp.plotone([freq_pulse], [sp_pulse], ax, colors=[clp.yellow], showlegend=False, lw=0.1, alpha=0.5,
            ylim=[-0.04, 1.05],
            xlim=[2, 18], xlabel="frequency [eV]", ylabel="$P(\omega)$ [arb. units]",
            alphaspacing=0.3)
plt.fill_between(freq_pulse, sp_pulse, step="pre", alpha=0.4, color=clp.yellow)

clp.adjust(tight_layout=False, savefile="figs1_pulse.pdf")
