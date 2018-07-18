import numpy as np
import numpy.fft as dft
from . import surface as bs


def homogeneous_composite_modulus(E, nu):
    """Returns elastic composite modulus for homogeneous system, i.e. both sides of fracture of same material."""
    return 1.0 / (2.0 * (1.0 - nu**2) / E)


class Results:
    """Elastic contact results, pressure and displacement field."""

    def __init__(self):
        self.p = None
        self.u = None

    @property
    def displacement(self):
        return self.u

    @displacement.setter
    def displacement(self, value):
        self.u = value

    def contact_area(self, dxy):
        """Absolute area of contact in [L2], assuming uniform grid spacing."""
        return len(self.p[self.p > 0.0]) * dxy**2

    def average_aperture(self, rigid_surface):
        """
        The average mechanical aperture as between the elastic body
        and the rigid surface.
        """
        # the stored displacements are only relative
        shifted_displacement = self.displacement - (np.max(self.displacement) - np.max(rigid_surface))
        aperture = shifted_displacement - rigid_surface
        return np.mean(aperture)


class ProgressBar:
    def __init__(self):
        self.fill = '█'
        self.length = 25
        self.decimals = 1
        self.min_abs_residual = None
        self.initial_fraction = None
        self.iteration = 1

    def __enter__(self):
        self.update()
        return self

    def __exit__(self, *args):
        print()

    def err_lim(self, min_abs_residual):
        self.min_abs_residual = min_abs_residual

    def update(self, abs_residual=None):
        if not abs_residual:
            return
        current_fraction = np.log10(abs_residual / self.min_abs_residual)
        if abs_residual and not self.initial_fraction:
            self.initial_fraction = current_fraction
        fraction = 1.0 - current_fraction / self.initial_fraction
        num_filled = abs(int(self.length * fraction))
        bar = self.fill * num_filled + '-' * (self.length - num_filled)
        print('\rConvergence |{0}| {1:3d}%, res={3:.2e}(max={4:.2e}), it={2:<4d}'.format(bar,
              abs(int(100.0 * fraction)), self.iteration, abs_residual, self.min_abs_residual), end='\r')
        self.iteration += 1


def _recon_FFT(x, y):
    """
    Internal method. Aides the convolution with the log of radius function.

    Parameters
    ----------
    x: array_like
        A dense 2D array filled with x coords.
    y: array_like
        A dense 2D array filled with y coords.

    Returns
    -------
    reco: array_like
        The convolution with the log of radius function.

    Notes
    -----
    Used in stiffness matrix computation of CG/FFT based contact algorithms.
    """
    r = np.sqrt(x * x + y * y)
    reco = y * np.log(x + r) + x * np.log(y + r)
    return reco


def _stiffness_FFT(N, d, E, nu):
    """
    Internal method. Computes the convolved cyclic influence (stiffness) matrix for a
    2D (3D in space) half-space elastic contact problem for FFT algorithms.

    Parameters
    ----------
    N:
        System size, Nx = Ny.
    d:
        Grid-spacing, dx = dy.
    E:
        Young's modulus, E1 = E2.
    nu:
        Poisson's ratio, nu1 = nu2

    Returns
    -------
    B: ndarray
        2D numpy cyclic influence matrix of size (2N,2N), means
        computational domain four times data domain.
    """
    Ee = homogeneous_composite_modulus(E, nu)
    C = 1.0 / (np.pi * Ee)
    p2 = d / 2.0
    L = d * N
    xyb = np.arange(0.0, L - p2, d)
    yyb, xxb = np.meshgrid(xyb, xyb)
    xxm = xxb - p2
    xxp = xxb + p2
    yym = yyb - p2
    yyp = yyb + p2
    A = _recon_FFT(xxm, yym) + _recon_FFT(xxp, yyp) - _recon_FFT(xxm, yyp) - _recon_FFT(xxp, yym)
    A = C * A
    B = np.zeros((2 * N, 2 * N))
    B[:N, :N] = A[:]
    B[N, :] = B[N - 1, :]
    B[N + 1:2 * N, :] = B[N - 1:0:-1, :]
    B[:, N] = B[:, N - 1]
    B[:, N + 1:2 * N] = B[:, N - 1:0:-1]
    return B


def contact_FFT(s, nominal_stress, E, nu, it_max=1000, err_lim=1.0E-10, initial_penetration_fraction=0.1, verbose=0):
    """
    Solves the elastic half-space contact problem (elastic plane-rigid composite profile) for a regular,
    rectangular grid of equi-distant spacing using a FFT/Conjugate Gradient algorithm.

    Parameters
    ----------
    s:
        An instance of surface.Surface concept, constituting the rough rigid
        composite profile, any [L].
    nominal_stress:
        Nominal normal stress acting on elastic surface [Pa].
    E:
        Young's modulus. Assumes E1 = E2 [Pa].
    nu:
        Poisson's ratio. Assumes nu1 = nu2 [-].
    it_max: otpional
        Upper limit of iterations.
    err_lim:
        The relative error limit beyond which convergence is achieved.
    initial_penetration_fraction: optional, float
        Initial guess for mean penetration of elastic plane w.r.t. maximum height of surface.
    verbose: bool, optional
        If true, prints iteration and residual. Defaults to false.

    Returns
    -------
    A Results instance with contact normal stress and displacement solution.

    Notes
    -----
    See Stanley and Kato (1999, Tribology);
    Liu et al. (2000, Wear); Sainsot and Lubrecht (2011, Eng. Tribology) for a descriptions of the algorithm.
    Here, we assume that both contacting surfaces have identical mechanical properties (Young's, Poisson's).
    Further, the initial guess is improved to be more problem specific. Also, the error limit is relative
    and depends on the first iteration error and goes err_lim below that.

    Todo
    ----
    Faster FFT that retains initial settings.
    """
    # getting discretization props
    N, dxy = s.shape[0], s.dxy
    # retrieving kernel (influence coefficients)
    A = _stiffness_FFT(N, dxy, E, nu)
    # and transformint from spatial to frequency domain
    fA = dft.fft2(A)

    # discretizing P and assigning initial guess as uniformly
    # distributed from target load (charge..[N])
    charge = nominal_stress * bs.nominal_area(s)
    # imporved initial guess,
    h_max = np.max(s)
    P = np.zeros((N, N))
    ic = np.where(s > h_max - h_max * initial_penetration_fraction)
    contact_area_initial_guess = ic[0].size * dxy**2
    P[ic] = charge / contact_area_initial_guess

    # convergence auxilliaries
    errlim = err_lim  # rel limit set first to initial absolute error
    pk = np.zeros((N, N))
    err, errs = 1., []
    it, gold = 0, 1.
    s *= -1.
    if verbose > 1:
        print('CG/FFT elastic contact algorithm on %d X %d grid, sigma_0 = %8.2e, E = %8.2e, nu = %5.2f' % (N, N, nominal_stress, E, nu))

    # CG loop
    with ProgressBar() as progress_bar:
        while err > errlim:
            if it > it_max:
                if verbose:
                        print('stopped due to maximum iteration constraint of %d' % (it_max))
                break
            it = it + 1
            sy = np.where(P > 0.0)  # contact area
            sn = np.where(P <= 0.0)  # free area
            # compute displacement in frequency domain and transform back to spatial domain
            dd = dft.ifft2(fA * dft.fft2(P, s=(2 * N, 2 * N)))
            dd = dd.real
            u = dd[0:N, 0:N]
            rk = u + s
            do = np.mean(rk[sy])
            rk = rk - do
            # norm
            G = np.sum(rk[sy] * rk[sy])
            # slopes
            pk[sy] = rk[sy] + (G / gold) * pk[sy]
            pk[sn] = 0.0
            gold = G
            # qk
            dd = dft.ifft2(fA * dft.fft2(pk, s=(2 * N, 2 * N)))
            dd = dd.real
            qk = dd[0:N, 0:N]
            # rb is the adjustment in approach
            rb = np.mean(qk[sy])
            qk = qk - rb
            # coeffs
            dp = np.sum(rk[sy] * pk[sy]) / np.sum(qk[sy] * pk[sy])
            P[sy] = P[sy] - dp * pk[sy]
            P[P < 0.0] = 0.0
            sol = np.where((P == 0.0) & (rk < 0.0))
            P[sol] = P[sol] - dp * rk[sol]
            W = np.sum(P) * dxy**2
            P = (charge / W) * P
            err = np.sqrt(gold * dxy**2)
            errs.append(err)
            if it == 1:
                errlim = err * errlim
                progress_bar.err_lim(errlim)
            if verbose:
                progress_bar.update(err)

    # return
    s *= -1.
    contact = Results()
    contact.p, contact.u = P, u
    return contact


def stiffness(nominal_stresses, rigid_surface, E, nu, contact=contact_FFT, **kwargs):
    """
    Computes stiffness given stresses, kwargs go to contact function.
    """
    first_call = [True]

    def average_aperture(nominal_stress):
        verbosity = 1
        if first_call[0]:
            verbosity = 2
            first_call[0] = False
        result = contact(rigid_surface, nominal_stress, E, nu, verbose=verbosity, **kwargs)
        return result.average_aperture(rigid_surface)

    apertures = [average_aperture(nominal_stress) for nominal_stress in nominal_stresses]

    # kappa = dp / d<a_m>, Eq. 2.8, negate aperture delta
    # here instead of abs to highlight errors
    stiffness = np.diff(nominal_stresses) / (-1.0 * np.diff(apertures))
    return stiffness


if __name__ == '__main__':
    import doctest
    doctest.testmod()
