import unittest
import numpy as np
import surfflow as sf
import surfcont as sc
import surfunits as su

verbose = False

class test_surfgen(unittest.TestCase):

    def test_base(self):
        su.length_unit_is_meter()

        N = 100
        a_field = np.zeros((N,N))
        aperture = 1.0E-5
        a_field[:,:] = aperture
        tolerance_places = 15
        a_analytic = aperture
        results = sf.hydraulic_aperture(a_field, 7.5E-5)
        if verbose:
            print 'a  analytic: %4.2E' % a_analytic, 'a numeric: %4.2E' % np.mean([results[0].a,results[1].a])
        self.assertAlmostEqual(results[0].a, a_analytic, places=tolerance_places)
        self.assertAlmostEqual(results[1].a, a_analytic, places=tolerance_places)

        tolerance_places = 7
        aperture_right = aperture
        aperture_left = aperture/2.
        a_field[:,:N/2] = aperture_left
        T_analytic_y = np.mean([aperture_left**3/12., aperture_right**3/12.])
        a_analytic_y = np.power(T_analytic_y*12., 1./3.)
        T_analytic_x = 2./(1./(aperture_left**3/12.)+1./(aperture_right**3/12.))  
        a_analytic_x = np.power(12.*T_analytic_x, 1./3.)
        results = sf.hydraulic_aperture(a_field, 7.5E-5)
        if verbose:
            print 'ax analytic: %4.2E' % a_analytic_x, 'ax numeric: %4.2E' % results[0].a, 'ay analytic: %4.2E' % a_analytic_y, 'ay numeric: %4.2E' % results[1].a
        self.assertAlmostEqual(results[0].a, a_analytic_x, places=tolerance_places)
        self.assertAlmostEqual(results[1].a, a_analytic_y, places=tolerance_places)



if __name__ == '__main__':
    verbose = True
    unittest.TextTestRunner().run(unittest.TestLoader().loadTestsFromTestCase(test_surfgen))
    raw_input('done...')