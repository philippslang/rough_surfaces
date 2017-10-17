import numpy as np
import scipy.integrate as scin
import scipy.optimize as scop
import warnings


def self_affine_psd(q, pref, hurst, onedim=False):
    """Ideal self-affine power spectrum, dependent only on prefactor and Hurst exponent."""
    exp = -2 * (hurst + 1)
    if onedim:
        exp = -1 - 2 * hurst
    return pref * q**exp


def self_affine_psd_fit(q, C, onedim=False):
    """Returns best-fit pre-factor and Hurst exponent for given roughness spectrum."""
    # we have to fit here in log space, else the high freq components are irrelevant
    # to least squares, magnitude-wise
    def self_affine_psd_log(q, pref, hurst):
        return np.log10(self_affine_psd(q, pref, hurst, onedim))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # ignoring unphysical vals during fit
        popt, pcov = scop.curve_fit(self_affine_psd_log, q, np.log10(C))
    return popt


def radially_integrated_psd(q, C):
    """Assumes no roll-off/cut-off in PSD."""
    if q[0] > q[1]:
        q = q[::-1]
        C = C[::-1]
    return 2 * np.pi * scin.simps(np.multiply(q, C), q)


def radially_averaged_psd(surface, window=True):
    """
    Computes the power spectral density/roughness spectrum from roughness readings using angular averaging.
    This radial averaging makes this algorithm appropriate for isotropic surfaces.

    Parameters
    ----------
    h: brown.surface.Surface (L)
        A square 2D array containing the surface height profile, equidistant measurements assumed in both directions.
    window: bool
        Hanning window applied to data before Discrete Fourier Transform to adhere to inherent periodicity. Defaults
        to True, not needed for periodic surfaces.

    Returns
    -------
    frequency, power: array_like
        Two 1D numpy arrays containing the roughness frequency (L^-1) vs. power (L^4).

    Notes
    -----
    The routine uses numpy's FFT package and follows the algorithm in Persson (2005), 10.1088/0953-8984/17/1/R01,
    and power normalization as in Elson & Bennet (1995), 10.1364/AO.34.000201.
    Todo
    ----
    Normalize.
    """
    h = np.copy(surface)
    N = h.shape[0]
    N_center = int(N / 2)
    dxy = surface.dxy
    L = N * dxy
    q_L = 2 * np.pi / L

    # windowing
    if window:
        h = h * np.outer(np.hanning(N), np.hanning(N))

    h_a_q = (dxy**2 / (2 * np.pi)**2) * np.fft.fft2(h)    # Equ. D.1
    C_q = ((2 * np.pi)**2 / L**2) * np.abs(h_a_q)**2      # Equ. C.6

    # first quadrant of FFT, rest follows in periodic averaging
    m = np.array(range(1, N_center + 1))         # frequency multiplier
    mxy = np.array(range(0, N_center))         # grid indices
    q = np.array(m) * q_L             # frequencies
    xv, yv = np.meshgrid(mxy, mxy)  # grid index meshes
    mxyabs = np.sqrt(xv**2 + yv**2)   # grid index vector lengths

    # averaging
    N = N - 1
    C = np.zeros_like(q)
    for i, mi in enumerate(m):
        mind, maxd = mi - 0.5, mi + 0.5
        idcs = np.where((mind < mxyabs) & (mxyabs < maxd))
        C_sum = 0.
        for j in range(len(idcs[0])):
            mx, my = idcs[0][j], idcs[1][j]
            C_sum += C_q[mx, my] + C_q[N - mx, my] + C_q[mx, N - my] + C_q[N - mx, N - my]
        C[i] = C_sum / (2 * np.pi * mi)

    # Parseval's theorem based power normalization
    meansqu_target = np.var(h)
    meansqu_is = radially_integrated_psd(q, C)
    C *= meansqu_target / meansqu_is
    return q, C


def axis_averaged_psd(h, dx, window=True, axis=0):
    """
    Average over 1D power spectra in x-direction (rows) if axis == 0, y-direction (columns) otherwise.
    Notes
    -----
    The routine uses numpy's FFT package and follows the algorithm in Elson & Bennet (1995), 10.1364/AO.34.000201.
    """
    h = np.copy(h)
    N = h.shape[0]
    L = N * dx
    if axis == 1:
        h = h.T

    # psds
    def power(r):
        if window:
            r *= np.hanning(N)
        return dx * np.abs(np.fft.fft(r))**2 / N  # Equ. 8a
    powers = np.array([power(line) for line in h])

    # freqs
    q = np.fft.fftfreq(N, dx) * 2 * np.pi
    q = q[1:int(N / 2)]

    # avg
    C = np.average(powers, axis=0)[1:int(N / 2)]
    return q, C
