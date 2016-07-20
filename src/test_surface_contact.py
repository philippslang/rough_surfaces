import unittest
import numpy as np
import surface_generation as sg
import surface_contact as sc


class TestSurfaceContact(unittest.TestCase):

    def setUp(self):
        self.F = [1.0E5]
        self.E = 1.0E+9
        self.nu = 0.3
        self.radius = 0.5 
        self.edge_length = 2.4*self.radius
        self.N = 256        

    def test_hertz(self):
        s = sg.sphere(self.N, self.edge_length, self.radius, 2.)        
        composite_radius = 1./(2.*(1./self.radius))
        composite_modulus = sc.homogeneous_composite_modulus(self.E, self.nu)
        for force in self.F:
            # Hertz solution
            contact_radius = (3*force*composite_radius/(4*composite_modulus))**(1/3)
            approach = contact_radius**2/composite_radius
            p_max = 3*force/(2*np.pi*contact_radius**2)
            # numerical solution
            contact_radius_num = np.sqrt(contact.contact_area(s.dxy)/np.pi)
            approach_num = np.max(contact.u)
            p_max_num = np.max(contact.p)
            # check
            np.isclose(contact_radius, contact_radius_num)
            np.isclose(approach, approach_num)
            np.isclose(p_max, p_max_num)


if __name__ == '__main__':
    unittest.main()

