import sys
import numpy as np
from . import surface as rs


def make_self_affine(saparams, power_of_two, seed=None):
    """
    Periodic self-affine surface based on an ideal radially averaged PSD. Roll-off only.

    Parameters
    ----------
    saparams:
        params.SelfAffineParameters
    power_of_two:
        System size such that number of points generated along each axis is N = 2^power_of_two
    seed: int
        Randomn number seed. If None, system time.

    Returns
    -------
    surface.Surface

    Notes
    -----
    See Fractal Imaging, Chapter 2
    """
    np.random.seed(seed)
    lambda_L_over_lambda_0 = 1 if saparams.lambda_L_over_lambda_0 is None else saparams.lambda_L_over_lambda_0
    lambda_L_over_lambda_1 = sys.maxsize if saparams.lambda_L_over_lambda_1 is None else saparams.lambda_L_over_lambda_1
    N, L = 2**power_of_two, saparams.dimensions[0]
    power = -(saparams.hurst + 1.0)
    f_L = 1.0 / N  # releative frequency
    f_0, f_1 = f_L * lambda_L_over_lambda_0, f_L * lambda_L_over_lambda_1
    A = np.zeros((N, N), dtype=complex)
    center = int(N / 2)
    rand_norm_1, rand_norm_2 = np.random.randn(center + 1, center + 1), np.random.randn(center + 1, center + 1)
    rand_unif_1, rand_unif_2 = np.random.rand(center + 1, center + 1), np.random.rand(center + 1, center + 1)
    for i in range(0, int(N / 2 + 1)):
        for j in range(0, int(N / 2 + 1)):
            phase = 2.0 * np.pi * rand_unif_1[i, j]
            rad = 0.0
            f = np.sqrt((float(i) / N)**2 + saparams.anisotropy * (float(j) / N)**2)
            if i != 0 or j != 0:
                f = f if f > f_0 else f_0          # hi pass for f_0 power
                rad = rand_norm_1[i, j] * f**power
            if f > f_1:                            # lo pass for f_1
                rad, phase = 0.0, 0.0
            A[i, j] = rad * np.cos(phase) + rad * np.sin(phase) * 1j
            i0 = 0 if i == 0 else N - i
            j0 = 0 if j == 0 else N - j
            A[i0, j0] = rad * np.cos(phase) - rad * np.sin(phase) * 1j
    A[center, 0] = A[center, 0].real + 0j
    A[0, center] = A[0, center].real + 0j
    A[center, center] = A[center, center].real + 0j
    for i in range(1, center):
        for j in range(1, center):
            phase = 2. * np.pi * rand_unif_2[i, j]
            f = np.sqrt((float(i) / N)**2 + saparams.anisotropy * (float(j) / N)**2)
            f = f if f > f_0 else f_0              # hi pass for f_0 power
            rad = rand_norm_2[i, j] * f**power
            if f > f_1:                            # lo pass for f_1
                rad, phase = 0.0, 0.0
            A[i, N - j] = rad * np.cos(phase) + rad * np.sin(phase) * 1j
            A[N - i, j] = rad * np.cos(phase) - rad * np.sin(phase) * 1j
    H = np.real(np.fft.ifft2((A)))
    s = rs.Surface(H, L / float(N))
    s = rs.scale_to_rms(s, saparams.hrms)
    s = rs.shift_to_zero_mean(s)
    return s


def sphere(N, edge_length, radius, scaling=1.0):
    """
    Creates a 2D Hertzian contact domain, a half sphere embedded in a plane.

    >>> r = 1.0
    >>> s = sphere(100, 1., r)
    >>> np.isclose(np.max(s.h), r, rtol=1.0E-4)
    True
    """
    x = np.linspace(-edge_length / 2.0, edge_length / 2.0, N)
    XX, YY = np.meshgrid(x, x)
    r = np.sqrt(XX * XX + YY * YY)
    h = np.zeros(r.shape)
    h[r <= radius] = np.sqrt(radius**2 - r[r <= radius]**2) * scaling
    dxy = edge_length / float(N)
    s = rs.Surface(h, dxy)
    return s

if __name__ == '__main__':
    import doctest
    doctest.testmod()
