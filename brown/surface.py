import numpy as np


class Surface:
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
    >>> s
    Surface((100x100), dxy=0.1)
    >>> s.length()  # egde length in x-direction
    10.0
    >>> s.length(1) # egde length in y-direction
    10.0
    
    Surfaces can also be one-dimensional, e.g., represent traces or cross-sections:
    >>> import numpy as np
    >>> N, dxy = 100, 0.1
    >>> h = np.zeros((N))
    >>> s = Surface(h, dxy)
    >>> s
    Surface((100), dxy=0.1)
    >>> s.length() # length
    10.0
    >>> s.length(1) # there is no second axis for one-dimensional surfaces
    Traceback (most recent call last):
        ...
    IndexError: tuple index out of range
    """

    def __init__(self, input_array, dxy):
        """"Input array [L] is copied to ndarray, and dxy [L] is cell size."""
        self.h = np.array(input_array)
        self.dxy = float(dxy)

    def __repr__(self):
        cls_name = type(self).__name__
        msg_shape = '{!r}x{!r}' if len(self.h.shape) == 2 else '{!r}'
        msg = '{}((' + msg_shape + '), dxy={!r})'
        return msg.format(cls_name, *self.h.shape, self.dxy)
    
    def __str__(self):
        return '\n'.join([self.__repr__(), str(self.h)]) 
    
    def __eq__(self, other):
        """"Hard comparison (no tol), first cheap lattice size, then element wise."""
        if not self.dxy == other.dxy:
            return False
        return np.equal(self.h, other.h).all()
    
    def __getitem__(self, pos):
        """"Surface height at index"""
        return self.h[pos]
    
    def __setitem__(self, pos, val):
        """"Surface height at index"""
        self.h[pos] = val
    
    def save(self, file):    
        np.savez(file, h=self.h, dxy=self.dxy)
        
    @classmethod
    def load(cls, file):
        d = np.load(file)
        return cls(d['h'], d['dxy'])
    
    @property
    def shape(self):
        return self.h.shape
    
    def rms(self):
        """"Returns root-mean-square roughness [L]."""
        return np.sqrt(np.mean(self.h**2))
    
    def length(self, axis=0):
        """"Returns length [L] of surface in x- or y-direction, for axis=0 and 1, respectively."""
        return self.shape[axis]*self.dxy
    
    def nominal_area(self):
        """"Returns length() [L] for 1D, area [L^2] for 2D."""        
        a = 1.0
        for i in range(len(self.shape)):
            a *= self.length(i)
        return a
    
    def shift_to_zero_mean(self):
        """"Shifts surface such that <h> = 0."""
        self.h = self.h-np.mean(self.h)
        
    def mean_aperture(self):
        """"Composite surface assumption: mean of difference field to highest point."""
        return np.mean(np.abs(np.subtract(self.h, np.max(self.h))))
    
    def pore_volume(self):
        """"Composite surface assumption: mean aperture times area (2D-->[L^3]) or length (1D-->[L^2])."""
        return self.mean_aperture()*self.nominal_area()

    def scale_to_property(self, prop, target_value):
        """
        Scales height to fit target property, which must be name of scalar returning method.
        
        >>> import numpy as np
        >>> N, dxy = 100, 0.1
        >>> h = np.ones((N))
        >>> s = Surface(h, dxy)
        >>> s.rms()
        1.0
        >>> s.scale_to_property('rms', 0.5)
        >>> s.rms()
        0.5
        """
        is_value = getattr(self, prop)()        
        self.h *= target_value/is_value


if __name__ == '__main__':
    import doctest
    doctest.testmod()