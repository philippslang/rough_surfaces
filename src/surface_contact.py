import numpy as np


def homogeneous_composite_modulus(E, nu):
    """Returns elastic composite modulus for homogeneous system, i.e. both sides of fracture of same material."""
    return 1./( 2.*((1.-nu**2)/(E)))


class Results:
    """Elastic contact results, pressure and displacement field."""
    
    def __init__(self):
        self.p = None
        self.u = None
        
    def contact_area(self, dxy):
        """Absolute area of contact in [L2], assuming uniform grid spacing."""
        return len(self.p[self.p>0.])*dxy**2
    
    
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
    r = np.sqrt(x*x+y*y)
    reco = y*np.log(x+r) + x*np.log(y+r)
    return reco


if __name__ == '__main__':
    import doctest
    doctest.testmod()