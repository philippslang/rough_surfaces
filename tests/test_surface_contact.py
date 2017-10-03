import unittest
import numpy as np
import brown.surface_generation as sg
import brown.surface_contact as sc


verbose = 0


def _rel_err(num, ref):
    return abs(num-ref)/abs(ref)


class TestSurfaceContact(unittest.TestCase):

    def setUp(self):
        self.F = [1.0E5,1.0E6]
        self.E = 1.0E+9
        self.nu = 0.3
        self.radius = 0.5 
        self.edge_length = 2.4*self.radius
        self.N = 300     

    def test_hertz(self):
        s = sg.sphere(self.N, self.edge_length, self.radius, 2.)        
        composite_radius = self.radius/2.
        composite_modulus = sc.homogeneous_composite_modulus(self.E, self.nu)
        for force in self.F:
            # Hertz solution
            contact_radius = (3*force*composite_radius/(4*composite_modulus))**(1/3)
            approach = contact_radius**2/composite_radius
            p_max = 3*force/(2*np.pi*contact_radius**2)
            # numerical solution
            nominal_stress = force/self.edge_length**2
            contact = sc.contact_FFT(s, nominal_stress, self.E, self.nu, verbose=verbose)
            contact_radius_num = np.sqrt(contact.contact_area(s.dxy)/np.pi)
            approach_num = np.max(contact.u)
            p_max_num = np.max(contact.p)
            # check
            if verbose:
                print('Relative error radius = {:4.1f}%'.format(100*_rel_err(contact_radius_num, contact_radius)))
                print('Relative error approach = {:4.1f}%'.format(100*_rel_err(approach_num, approach)))
                print('Relative error pressure = {:4.1f}%'.format(100*_rel_err(p_max_num, p_max)))
            assert(np.isclose(contact_radius_num, contact_radius, rtol=0.01))
            assert(np.isclose(approach_num, approach, rtol=0.01))
            assert(np.isclose(p_max_num, p_max, rtol=0.01))


if __name__ == '__main__':
    verbose = 1
    unittest.main()

