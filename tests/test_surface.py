import unittest
import numpy as np
import brown.surface as bs


class TestSurface(unittest.TestCase):

    def setUp(self):
        self.N, self.dxy, self.val_init = 100, 1.0, 1.0
        self.L = self.dxy * self.N

    def _test_common(self, s):
        assert(self.val_init == bs.rms(s))
        s = bs.shift_to_zero_mean(s)
        assert(bs.rms(s) == 0.)
        """
        fname = '/tmp/save_surface.npz'
        bs.save(fname)
        sclone = sr.Surface.load(fname)
        assert(type(s) == type(sclone))
        assert(s == sclone)
        with self.assertRaises(AttributeError):
            s.shape = 10
        """

    def test_interface_2D(self):
        a = np.full((self.N, self.N), self.val_init)
        s = bs.Surface(a, self.dxy)
        assert(s[0, 0] == self.val_init)
        assert(s[0, self.N - 1] == self.val_init)
        with self.assertRaises(IndexError):
            s[0, self.N]
        new_val = 2.0
        s[0, self.N - 1] = new_val
        assert(s[0, self.N - 1] == new_val)
        s[0, self.N - 1] = self.val_init
        assert(s.dxy == self.dxy)
        assert(bs.length(s) == self.L)
        assert(bs.length(s, 1) == self.L)
        assert(bs.nominal_area(s) == self.L**2)
        self._test_common(s)

    def test_interface_1D(self):
        a = np.full((self.N), self.val_init)
        s = bs.Surface(a, self.dxy)
        assert(s.dxy == self.dxy)
        assert(bs.length(s) == self.L)
        assert(bs.nominal_area(s) == self.L)
        assert(s[0] == self.val_init)
        new_val = 2.0
        s[0] = new_val
        assert(s[0] == new_val)
        s[0] = self.val_init
        with self.assertRaises(IndexError):
            bs.length(s, 1)
        self._test_common(s)
