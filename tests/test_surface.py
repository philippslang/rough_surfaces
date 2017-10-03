import unittest
import numpy as np
import brown.surface as sr

class TestSurface(unittest.TestCase):

    def setUp(self):
        self.N, self.dxy, self.val_init = 100, 1.0, 1.0
        self.L = self.dxy*self.N        
        
    def _test_common(self, s):
        assert(self.val_init == s.rms())
        s.shift_to_zero_mean()
        assert(s.rms() == 0.)
        fname = '/tmp/save_surface.npz'
        s.save(fname)
        sclone = sr.Surface.load(fname)
        assert(type(s) == type(sclone))
        assert(s == sclone)
        with self.assertRaises(AttributeError):
            s.shape = 10

    def test_interface_2D(self):
        a = np.full((self.N,self.N), self.val_init)
        s = sr.Surface(a, self.dxy)
        assert(s[0,0] == self.val_init)
        assert(s[0,self.N-1] == self.val_init)
        with self.assertRaises(IndexError):
            s[0,self.N]
        new_val = 2.0
        s[0,self.N-1] = new_val
        assert(s[0,self.N-1] == new_val)
        s[0,self.N-1] = self.val_init
        assert(s.dxy == self.dxy)
        assert(s.length() == self.L)
        assert(s.length(1) == self.L)
        assert(s.nominal_area() == self.L**2) 
        self._test_common(s)
        
    def test_interface_1D(self):
        a = np.full((self.N), self.val_init)
        s = sr.Surface(a, self.dxy)        
        assert(s.dxy == self.dxy)
        assert(s.length() == self.L)
        assert(s.nominal_area() == self.L)
        assert(s[0] == self.val_init)
        new_val = 2.0
        s[0] = new_val
        assert(s[0] == new_val)
        s[0] = self.val_init
        with self.assertRaises(IndexError):
            s.length(1)
        self._test_common(s)


if __name__ == '__main__':
    unittest.main()

