import numpy as np
import columnplots as clp
from numpy import linalg as LA
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt


au2fs = 0.02418884254
au2eV = 27.211399
eV2cminv = 8065.540106923572

sigma = 10000.0
w_step = 0.00001
linewidth = 4.0 / sigma * au2eV
e_cutoff_ev = 10.0
e_cutoff_au = e_cutoff_ev / au2eV
Angstrom2AU = 1.8897259885789

def get_tddft_energy(filename):
    data = np.loadtxt(filename)
    e_scf = data[0] * au2eV
    e_exc = data[1:3] + e_scf
    e_tot = np.zeros(np.size(e_exc)+1)
    e_tot[0] = e_scf
    e_tot[1:] = e_exc
    return e_tot

def get_tddft_transition_dipole(filename):
    data = np.loadtxt(filename)
    dipole = data[:,1]
    return dipole

def gather_tddft_data():
    OA = np.array([2.682115, -0.219090])
    H0 = np.array([1.844187, 1.135044])
    OA_H0 = OA - H0
    Rs = []
    e0, e1 = [], []
    dipoles = []
    for i in range(30):
        ratio = i*0.01
        filename = "../oHBA/PES_scan/oHBA_rgs_%.2f.out" %ratio
        print("dealing with %s" %filename)
        e_tot = get_tddft_energy(filename+".energy")
        dipole = get_tddft_transition_dipole(filename+".transdip")
        dipoles.append(dipole)
        # calculate the OH distance
        dR = ratio * OA_H0
        R = np.sqrt(np.sum(dR**2))
        Rs.append(R)

        e0.append(e_tot[0])
        e1.append(e_tot[1])
    Rs = np.array(Rs)
    e0 = np.array(e0)
    rs = [Rs] * 2
    es = [e0 - e0[0], e1- e0[0]]
    dipoles = np.array([d[0] for d in dipoles])
    return rs, es, dipoles

def calculate_tddft_polariton(gs, es, ds, coupling=0.0):
    LPs, UPs = [], []
    for i in range(len(gs)):
        g, e, d = gs[i], es[i], ds[i]
        H = [[g, d*coupling], [d*coupling, e]]
        H = np.array(H)
        w, v = LA.eig(H)
        LP = np.min(w)
        UP = np.max(w)
        LPs.append(LP)
        UPs.append(UP)
    LPs = np.array(LPs)
    UPs = np.array(UPs)
    return LPs, UPs

def pade(time,signal,sigma=sigma,max_len=None,w_min=0.0,w_max=e_cutoff_au,w_step=w_step,read_freq=None):
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

def get_dipole_traj(filename, OD, OA, H0):
    data = np.loadtxt(filename)
    t, mux, muy, muz = data[:,0], data[:,1], data[:,2], data[:,3]
    t_fs = t * au2fs
    mu_norm = np.sqrt(mux**2 + muy**2 + muz**2)
    # calculate the distance to the neighboring two Oxygens
    H_disp = data[:,1:4]
    H = H_disp + H0
    R_HOD = np.sqrt(np.sum((H - OD)**2, axis=1)) / Angstrom2AU
    R_HOA = np.sqrt(np.sum((H - OA)**2, axis=1)) / Angstrom2AU
    return t_fs, R_HOD, R_HOA

def get_PT_traj(pathname, OD_Ang, OA_Ang, H0_au):
    OD = OD_Ang * Angstrom2AU
    OA = OA_Ang * Angstrom2AU
    H0 = H0_au
    t_fs, R_HOD, R_HOA = get_dipole_traj(filename=pathname+"/mu_n.txt", OD=OD, OA=OA, H0=H0)
    return t_fs, R_HOD, R_HOA

def prepare_PT_Figb():
    pathname_outcav = "../oHBA/tddft_S0S1_geometry/outcav_pulse/Geom_10_Amp_2e-2/"
    pathname_incav = "../oHBA/tddft_S0S1_geometry/incav_pulse/Geom_10_Amp_2e-2/"
    OD_rgs_Ang = np.array([1.014139, 1.694260,   0.000000])
    OA_rgs_Ang = np.array([2.668602, -0.278478,   0.000000])
    OD_OA_dis = np.sum((OD_rgs_Ang - OA_rgs_Ang)**2)**0.5
    print("OD-OA separation is %.6f" %OD_OA_dis)
    H0_au = np.array([3.6168e+00,   2.2602e+00,   8.6128e-16])
    t_fs_out, R_HOD_out, R_HOA_out = get_PT_traj(pathname_outcav, OD_rgs_Ang, OA_rgs_Ang, H0_au)
    t_fs_in, R_HOD_in, R_HOA_in = get_PT_traj(pathname_incav, OD_rgs_Ang, OA_rgs_Ang, H0_au)
    xs = [t_fs_out, t_fs_out, t_fs_in, t_fs_in]
    ys = [R_HOD_out, R_HOA_out, R_HOD_in, R_HOA_in]
    colors = [clp.navy_blue, clp.navy_blue, clp.red, clp.red]
    labels = [r"H--O$_{\rm D}$ cavity off", r"H--O$_{\rm A}$ cavity off", r"H--O$_{\rm D}$ cavity on", r"H--O$_{\rm A}$ cavity on"]
    return xs, ys, colors, labels

def prepare_PT_Figc():
    pathname_outcav = "../oHBA/tddft_S0S1_geometry/outcav_pulse/Geom_10_Amp_4e-2/"
    pathname_incav = "../oHBA/tddft_S0S1_geometry/incav_pulse/Geom_10_Amp_4e-2/"
    OD_rgs_Ang = np.array([1.014139, 1.694260,   0.000000])
    OA_rgs_Ang = np.array([2.668602, -0.278478,   0.000000])
    OD_OA_dis = np.sum((OD_rgs_Ang - OA_rgs_Ang)**2)**0.5
    print("OD-OA separation is %.6f" %OD_OA_dis)
    H0_au = np.array([3.6168e+00,   2.2602e+00,   8.6128e-16])
    t_fs_out, R_HOD_out, R_HOA_out = get_PT_traj(pathname_outcav, OD_rgs_Ang, OA_rgs_Ang, H0_au)
    t_fs_in, R_HOD_in, R_HOA_in = get_PT_traj(pathname_incav, OD_rgs_Ang, OA_rgs_Ang, H0_au)
    xs = [t_fs_out, t_fs_out, t_fs_in, t_fs_in]
    ys = [R_HOD_out, R_HOA_out, R_HOD_in, R_HOA_in]
    colors = [clp.navy_blue, clp.navy_blue, clp.red, clp.red]
    labels = [r"H--O$_{\rm D}$ cavity off", r"H--O$_{\rm A}$ cavity off", r"H--O$_{\rm D}$ cavity on", r"H--O$_{\rm A}$ cavity on"]
    return xs, ys, colors, labels

def prepare_PT_Figd():
    pathname_outcav = "../oHBA/tddft_S0S1_geometry/outcav_pulse/Geom_10_Amp_8e-2/"
    pathname_incav = "../oHBA/tddft_S0S1_geometry/incav_pulse/Geom_10_Amp_8e-2/"
    OD_rgs_Ang = np.array([1.014139, 1.694260,   0.000000])
    OA_rgs_Ang = np.array([2.668602, -0.278478,   0.000000])
    OD_OA_dis = np.sum((OD_rgs_Ang - OA_rgs_Ang)**2)**0.5
    print("OD-OA separation is %.6f" %OD_OA_dis)
    H0_au = np.array([3.6168e+00,   2.2602e+00,   8.6128e-16])
    t_fs_out, R_HOD_out, R_HOA_out = get_PT_traj(pathname_outcav, OD_rgs_Ang, OA_rgs_Ang, H0_au)
    t_fs_in, R_HOD_in, R_HOA_in = get_PT_traj(pathname_incav, OD_rgs_Ang, OA_rgs_Ang, H0_au)
    xs = [t_fs_out, t_fs_out, t_fs_in, t_fs_in]
    ys = [R_HOD_out, R_HOA_out, R_HOD_in, R_HOA_in]
    colors = [clp.navy_blue, clp.navy_blue, clp.red, clp.red]
    labels = [r"H--O$_{\rm D}$ cavity off", r"H--O$_{\rm A}$ cavity off", r"H--O$_{\rm D}$ cavity on", r"H--O$_{\rm A}$ cavity on"]
    return xs, ys, colors, labels

def get_dipole(filename):
    data = np.loadtxt(filename)
    nskip = 2
    nst = int(np.shape(data)[0] * 0.1)
    nend = int(np.shape(data)[0] * 1)
    t, mux, muy, muz = data[nst:nend:nskip,0], data[nst:nend:nskip,1], data[nst:nend:nskip,2], data[nst:nend:nskip,3]
    spx, freq = pade(t, mux)
    spy, freq = pade(t, muy)
    spz, freq = pade(t, muz)

    t_fs = t * au2fs
    #mu_norm = np.sqrt(mux**2 + muy**2 + muz**2)
    mu_norm = muz
    freq_ev = freq * au2eV
    sp_tot = freq * (np.abs(spx) + np.abs(spy) + np.abs(spz))
    #sp_tot = freq * np.abs(spz)
    df = freq_ev[1] - freq_ev[0]
    # output only selectively number of data
    nmax = int(e_cutoff_ev // df)
    freq_ev = freq_ev[0:nmax]
    sp_tot = sp_tot[0:nmax]
    return t_fs, mu_norm, freq_ev, sp_tot

def get_polariton_spectrum(pathname, w_start=2.5, w_end=4.5):
    print("dealing with %s" %pathname)
    t_e, mu_e, freq_e, sp_e = get_dipole(filename=pathname+"/mu_e.txt")
    df = freq_e[1] - freq_e[0]
    n_start = int(w_start // df)
    n_end = int(w_end // df)
    sp_truncated = sp_e[n_start:n_end]
    freq_truncated = freq_e[n_start:n_end]
    sp_truncated /= np.sum(sp_truncated) * df
    return sp_truncated, freq_truncated

def get_Efield_profile(t_au, t0=0, sigma=400.0, omega_eV=3.2959, amp_efield=1e-3):
    omega = omega_eV / au2eV
    pulse = amp_efield * np.exp(-(t_au-t0)**2/sigma**2) * np.cos(omega*t_au)
    return pulse

def get_pulse_spectrum(sigma=400.0, omega_eV=3.111):
    t_au = np.linspace(0, 500, 501)
    pulse = get_Efield_profile(t_au, t0=0, sigma=sigma, omega_eV=omega_eV)
    sp, freq = pade(t_au, pulse, sigma=2e3)
    freq_eV = freq * au2eV
    sp = np.abs(sp) * freq
    sp /= np.max(sp)
    return sp * 8, freq_eV

def prepare_spectrum(idx):
    pathname_outcav = "../oHBA/tddft_S0S1_geometry/outcav_LR/Geom_%d_LR/" %idx
    pathname_incav = "../oHBA/tddft_S0S1_geometry/incav_LR/Geom_%d_LR/" %idx
    sp_out, freq_out = get_polariton_spectrum(pathname_outcav)
    sp_in, freq_in = get_polariton_spectrum(pathname_incav)
    xs = [freq_out, freq_in]
    ys = [sp_out + 10, sp_in]
    colors = [clp.navy_blue, clp.red]
    labels = ["cavity off", "cavity on"]
    return xs, ys, colors, labels

# we also plot the pulse spectrum as a function of time
sp_pulse_out, freq_pulse = get_pulse_spectrum(sigma=400.0, omega_eV=3.611)
sp_pulse_in, freq_pulse = get_pulse_spectrum(sigma=400.0, omega_eV=3.355)
xs_pulse = [freq_pulse, freq_pulse]
ys_pulse = [sp_pulse_out + 10, sp_pulse_in]

#ax = clp.initialize()
#clp.plotone(xs_pulse, ys_pulse, ax)
#clp.adjust()

# Here, I will plot a figure to collect all results in the proton transfer simulation

fig, gs = clp.initialize_gridSpec(col=4, row=4, width=10, height=6,
        LaTeX=True, fontsize=12)

# Fig. a. Polaritonic energy surface
ax_a1 = fig.add_subplot(gs[0:2, 0])
ax_a2 = fig.add_subplot(gs[2:, 0])
ax_as = [ax_a1, ax_a2]

rs, es, dipoles = gather_tddft_data()
rs2 = [rs[0]]
es2 = [es[0] + es[1][0]]
LPs, UPs = calculate_tddft_polariton(es[0] + es[1][0], es[1], dipoles, coupling=0.3)

for ax_a in ax_as:
    clp.plotone(rs, es, ax_a, colors=['k']*len(rs), lw=1, showlegend=False, xlim=[0, 0.45])
    clp.plotone(rs2, es2, ax_a, colors=[clp.brown], lw=1, showlegend=False, xlim=[0, 0.45])
    clp.plotone(rs2*2, [LPs, UPs], ax_a, colors=[clp.red, clp.red], lw=1,
    #xlabel="proton coordinate [Angstrom]", ylabel="potential [eV]",
    showlegend=False, xlim=[0, 0.45])
clp.broken_y(ax=ax_a1, ax2=ax_a2, d=0.015, ymin_0=0, ymax_0=1.3, ymin_1=2.9, ymax_1=4.2)

ax_a2.set_xlabel("proton displacement $\Delta R$ [Angstrom]")
fig.text(-0.014, 0.5, "potential energy [eV]", va='center', rotation='vertical', fontsize=12)
# add annotation to Fig. a
ax_a2.text(0.2, 0.4, "S$_0$", color="k", fontsize=12)
ax_a1.text(0.2, 3.1, "LP", color=clp.red, fontsize=12)
ax_a1.text(0.2, 3.9, "UP", color=clp.red, fontsize=12)
ax_a1.text(0.216, 3.6, "S$_0$ + 1 photon", color=clp.brown, fontsize=12)
ax_a1.text(0.2, 3.4, "S$_1$", color="k", fontsize=12)

ax_a1.text(-0.05, 1.05, "(a)", fontsize=14, transform=ax_a1.transAxes)

# Fig. b inset: plot polariton spectrum
axins_b = fig.add_subplot(gs[0, 1])
xs, ys, colors, labels = prepare_spectrum(idx=10)
clp.plotone(xs, ys, axins_b, colors=colors, labels=labels, showlegend=True,
            lw=1, xlabel="energy [eV]", ylabel="spectrum [arb. units]",
            xlim=[2.5, 4.5], ylim=[0, 43], legendEdgeColor='k')

clp.plotone(xs_pulse, ys_pulse, axins_b, colors=[clp.yellow, clp.yellow], showlegend=False, lw=1)
plt.fill_between(xs_pulse[0], ys_pulse[0], np.ones(np.shape(xs_pulse[0]))*10, step="pre", alpha=0.4, color=clp.yellow)
plt.fill_between(xs_pulse[1], ys_pulse[1], step="pre", alpha=0.4, color=clp.yellow)

axins_b.text(-0.05, 1.05, "(b)", fontsize=14, transform=axins_b.transAxes)
axins_b.legend(prop={'size': 10}, edgecolor="k")
# Fig. b. The corresponding proton transfer inside versus outside the cavity: Weak PT
ax_b = fig.add_subplot(gs[1:, 1])
xs, ys, colors, labels = prepare_PT_Figb()
clp.plotone(xs, ys, ax_b, colors=colors, labels=labels, lw=1, ylim=[1.0, 1.8],
            xlabel="time [fs]", ylabel="distance [Angstrom]", xlim=[0, 10],
            linestyles=["-", "--", "-", "--"], legendEdgeColor='k')
ax_b.text(-0.05, 1.05, "(c)", fontsize=14, transform=ax_b.transAxes)
ax_b.legend(prop={'size': 10}, edgecolor="k")

# Fig. c inset: plot polariton spectrum
axins_c = fig.add_subplot(gs[0, 2])
xs, ys, colors, labels = prepare_spectrum(idx=10)
clp.plotone(xs, ys, axins_c, colors=colors, labels=labels, showlegend=False,
            lw=1, xlabel="energy [eV]", xlim=[2.5, 4.5], ylim=[0, 43])
axins_c.text(-0.05, 1.05, "(d)", fontsize=14, transform=axins_c.transAxes)

clp.plotone(xs_pulse, ys_pulse, axins_c, colors=[clp.yellow, clp.yellow], showlegend=False, lw=1)
plt.fill_between(xs_pulse[0], ys_pulse[0], np.ones(np.shape(xs_pulse[0]))*10, step="pre", alpha=0.4, color=clp.yellow)
plt.fill_between(xs_pulse[1], ys_pulse[1], step="pre", alpha=0.4, color=clp.yellow)

# Fig. c. Proton transfer under a more relaxed geometry: Large PT
ax_c = fig.add_subplot(gs[1:, 2])
xs, ys, colors, labels = prepare_PT_Figc()
clp.plotone(xs, ys, ax_c, colors=colors, labels=labels, lw=1,
            xlabel="time [fs]",  xlim=[0, 10], ylim=[1.0, 1.8],
            showlegend=False,
            linestyles=["-", "--", "-", "--"])
ax_c.text(-0.05, 1.05, "(e)", fontsize=14, transform=ax_c.transAxes)

# Fig. d inset: plot polariton spectrum
axins_d = fig.add_subplot(gs[0, 3])
xs, ys, colors, labels = prepare_spectrum(idx=10)
clp.plotone(xs, ys, axins_d, colors=colors, labels=labels, showlegend=False,
            lw=1, xlabel="energy [eV]", xlim=[2.5, 4.5], ylim=[0, 43])
axins_d.text(-0.05, 1.05, "(f)", fontsize=14, transform=axins_d.transAxes)

clp.plotone(xs_pulse, ys_pulse, axins_d, colors=[clp.yellow, clp.yellow], showlegend=False, lw=1)
plt.fill_between(xs_pulse[0], ys_pulse[0], np.ones(np.shape(xs_pulse[0]))*10, step="pre", alpha=0.4, color=clp.yellow)
plt.fill_between(xs_pulse[1], ys_pulse[1], step="pre", alpha=0.4, color=clp.yellow)

# Fig. d. Proton transfer under optimized geometry: No PT
ax_d = fig.add_subplot(gs[1:, 3])
xs, ys, colors, labels = prepare_PT_Figd()
clp.plotone(xs, ys, ax_d, colors=colors, labels=labels, lw=1,
            xlabel="time [fs]", xlim=[0, 10], ylim=[1.0, 1.8],
            showlegend=False,
            linestyles=["-", "--", "-", "--"])
ax_d.text(-0.05, 1.05, "(g)", fontsize=14, transform=ax_d.transAxes)

#axins_b.axvline(x=3.295, c=clp.navy_blue, lw=1.0, ls="--")
#axins_c.axvline(x=3.611, c=clp.navy_blue, lw=1.0, ls="--")
#axins_d.axvline(x=3.890, c=clp.navy_blue, lw=1.0, ls="--")

# finally, we insert some images
import matplotlib.image as mpimg
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox

axis_img = mpimg.imread('axis_xy.png')
imagebox = OffsetImage(axis_img, zoom=0.15)
ab_axis = AnnotationBbox(imagebox, (0.067, 0.4), frameon=False)
ab_axis.set(zorder=-1)
ax_a2.add_artist(ab_axis)

data_img = mpimg.imread('oHBA_rgs_0_noghost.ps')
data_img = data_img[230:791-250,100:612-100,:]
imagebox = OffsetImage(data_img, zoom=0.25)
ab = AnnotationBbox(imagebox, (0.2, 0.8), frameon=False)
ab.set(zorder=-1)
ax_a2.add_artist(ab)
ax_a2.text(0.269, 1.004, r"O$_{\rm D}$", color="k", fontsize=10)
ax_a2.text(0.309, 0.914, "H", color="k", fontsize=10)
ax_a2.text(0.346, 0.802, r"O$_{\rm A}$", color="k", fontsize=10)

'''
data_img = mpimg.imread('oHBA_rgs_0.ps')
data_img = data_img[230:791-250,100:612-100,:]
imagebox = OffsetImage(data_img, zoom=0.28)
ab = AnnotationBbox(imagebox, (1.3, 1.3), frameon=False)
ab.set(zorder=-1)
ax_b.add_artist(ab)

data_img = mpimg.imread('oHBA_rgs_10.ps')
data_img = data_img[230:791-250,100:612-100,:]
imagebox = OffsetImage(data_img, zoom=0.28)
ab = AnnotationBbox(imagebox, (2.2, 1.32), frameon=False)
ab.set(zorder=-1)
ax_c.add_artist(ab)

data_img = mpimg.imread('oHBA_rgs_20.ps')
data_img = data_img[230:791-250,100:612-100,:]
imagebox = OffsetImage(data_img, zoom=0.28)
ab = AnnotationBbox(imagebox, (2.2, 1.35), frameon=False)
ab.set(zorder=-1)
ax_d.add_artist(ab)
'''

# finally draw an arrow
import matplotlib.patches
ax0tr = ax_b.transData # Axis 0 -> Display
ax1tr = ax_d.transData # Axis 1 -> Display
figtr = fig.transFigure.inverted() # Display -> Figure
# 2. Transform arrow start point from axis 0 to figure coordinates
ptB = figtr.transform(ax0tr.transform((3., 0.8)))
# 3. Transform arrow end point from axis 1 to figure coordinates
ptE = figtr.transform(ax1tr.transform((10., 0.8)))
# 4. Create the patch
arrow = matplotlib.patches.FancyArrowPatch(
    ptB, ptE, transform=fig.transFigure,  # Place arrow in figure coord system
    fc = clp.sky_blue, connectionstyle="arc3,rad=0.0", arrowstyle='simple',
    ec="k",
    mutation_scale = 40.
)

fig.patches.append(arrow)

fig.text(0.5, -0.08, r"increasing pulse amplitude", fontsize=14, color=clp.sky_blue)

# add additional patch to connect H and O_D
ax0tr = ax_a2.transData # Axis 0 -> Display
ax1tr = ax_a2.transData # Axis 1 -> Display
figtr = fig.transFigure.inverted() # Display -> Figure
# 2. Transform arrow start point from axis 0 to figure coordinates
ptE = figtr.transform(ax0tr.transform((0.135, 0.8)))
# 3. Transform arrow end point from axis 1 to figure coordinates
ptB = figtr.transform(ax1tr.transform((0.106, 0.89)))
# 4. Create the patch
arrow = matplotlib.patches.FancyArrowPatch(
    ptB, ptE, transform=fig.transFigure,  # Place arrow in figure coord system
    fc = "b", connectionstyle="arc3,rad=0.0", arrowstyle='simple', alpha = 0.4,
    mutation_scale = 10.
)

fig.patches.append(arrow)

clp.adjust(savefile="PT_pulse.pdf")
