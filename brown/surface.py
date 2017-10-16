import numpy as np


class Surface(np.ndarray):
    """
    One- or two-dimensional surface height representation.

    The assumption upon which this framework is based is a uniform lattice size in both directions.
    This is tightly integrated here. 'Surface' is the fundamental class that most modules build
    upon. It usually represents the model or computational domain, as it may discretize either,
    individual and composite surfaces, i.e., rough surfaces and aperture fields.

    Standard initialization is from two-dimensional ndarray and lattice size:
    >>> import numpy as np
    >>> N, dxy = 100, 0.1
    >>> h = np.zeros((N,N))
    >>> s = Surface(h, dxy)
    >>> length(s)  # egde length in x-direction
    10.0
    >>> length(s, 1) # egde length in y-direction
    10.0

    Surfaces can also be one-dimensional, e.g., represent traces or cross-sections:
    >>> import numpy as np
    >>> N, dxy = 100, 0.1
    >>> h = np.zeros((N))
    >>> s = Surface(h, dxy)
    >>> length(s) # length
    10.0
    >>> length(s, 1) # there is no second axis for one-dimensional surfaces
    Traceback (most recent call last):
        ...
    IndexError: tuple index out of range
    """

    def __new__(cls, input_array, dxy):
        obj = np.asarray(input_array).view(cls)
        obj.dxy = float(dxy)
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            self.dxy = getattr(obj, 'dxy', None)


def rms(surface):
    """"Returns root-mean-square roughness [L]."""
    return np.sqrt(np.mean(surface**2))


def length(surface, axis=0):
    """"Returns length [L] of surface in x- or y-direction, for axis=0 and 1, respectively."""
    return surface.shape[axis] * surface.dxy


def nominal_area(surface):
    """"Returns length() [L] for 1D, area [L^2] for 2D."""
    a = 1.0
    for i in range(len(surface.shape)):
        a *= length(surface)
    return a


def shift_to_zero_mean(surface):
    """"Returns shifted surface such that <h> = 0."""
    return Surface(surface - np.mean(surface), surface.dxy)


def mean_aperture(surface):
    """"Composite surface assumption: mean of difference field to highest point."""
    return np.mean(np.abs(np.subtract(surface, np.max(surface))))


def pore_volume(surface):
    """"Composite surface assumption: mean aperture times area (2D-->[L^3]) or length (1D-->[L^2])."""
    return mean_aperture(surface) * nominal_area(surface)


def scale_to_rms(surface, rms_target):
    """
    Scales height to fit target property, which must be name of scalar returning method.
    """
    rms_current = rms(surface)
    return Surface(surface * (rms_target / rms_current), surface.dxy)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
