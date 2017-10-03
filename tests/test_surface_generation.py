import unittest
import numpy as np
import brown.generate as sg


class TestSurfaceGeneration(unittest.TestCase):

    def setUp(self):
        self.N, self.dxy, self.val_init = 100, 1.0, 1.0
        self.L = self.dxy*self.N                

    def test_self_affine(self):
        pass
     

if __name__ == '__main__':
    unittest.main()

