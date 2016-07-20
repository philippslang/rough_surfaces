import unittest
import numpy as np
import surface as sr

class TestSurfaceContact(unittest.TestCase):

    def setUp(self):
        self.N, self.dxy, self.val_init = 100, 1.0, 1.0
        self.L = self.dxy*self.N        

    def test_hertz(self):
        assert(1 == 1)
        


if __name__ == '__main__':
    unittest.main()

