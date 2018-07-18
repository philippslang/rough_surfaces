import pytest
import collections
import numpy as np
import rough_surfaces.surface as rs
import rough_surfaces.generate as rg
import rough_surfaces.contact as rc


Parameters = collections.namedtuple('Parameters', 'force E nu radius edge_length N dxy')


@pytest.fixture(params=[1.0E5, 1.0E6])
def contact_numeric(request):
    force = request.param
    radius = 0.5
    edge_length = 2.4 * radius
    N = 300
    dxy = edge_length / N
    params = Parameters(force, E=1.0E+9, nu=0.3, radius=radius, edge_length=edge_length, N=N, dxy=dxy)
    surface = rg.sphere(params.N, params.edge_length, params.radius, 2.0)
    nominal_stress = params.force / params.edge_length**2
    contact = rc.contact_FFT(surface, nominal_stress, params.E, params.nu)
    return contact, params


def _rel_err(num, ref):  # also dont do this TODO use np.allclose
    return abs(num - ref) / abs(ref)


def test_hertz(contact_numeric):
    # numeric solution
    contact, params = contact_numeric
    # Hertz solution
    composite_radius = params.radius / 2.0
    composite_modulus = rc.homogeneous_composite_modulus(params.E, params.nu)
    contact_radius = (3.0 * params.force * composite_radius / (4 * composite_modulus))**(1 / 3)
    approach = contact_radius**2 / composite_radius
    p_max = 3.0 * params.force / (2 * np.pi * contact_radius**2)
    # compare
    contact_radius_num = np.sqrt(contact.contact_area(params.dxy) / np.pi)
    approach_num = np.max(contact.u)
    p_max_num = np.max(contact.p)
    assert np.isclose(contact_radius_num, contact_radius, rtol=0.01)
    assert np.isclose(approach_num, approach, rtol=0.01)
    assert np.isclose(p_max_num, p_max, rtol=0.01)


def test_aperture():
    dim = 10000
    displacement = np.ones((dim, dim))
    rigid_surface = np.zeros_like(displacement)
    aperture = 0.5
    rigid_surface[0, 0] = aperture
    contact_results = rc.Results()
    contact_results.displacement = displacement
    aperture_to_test = contact_results.average_aperture(rigid_surface)
    assert np.isclose(aperture, aperture_to_test, rtol=(1.0 / dim**2))


def test_stiffness():
    dim = 10000
    rigid_surface = rs.Surface(np.zeros((dim, dim)), 1.0)
    max_height = 11.0
    rigid_surface[0, 0] = max_height
    nominal_stresses = np.arange(1.0, 10.0)

    def contact_mock(rigid_surface, nominal_stress, E, nu, verbose):
        result = rc.Results()
        result.displacement = np.zeros_like(rigid_surface)
        result.displacement.fill(nominal_stress)
        result.displacement[0, 0] = max_height
        return result

    stiffness_to_test = rc.stiffness(nominal_stresses, rigid_surface, None, None, contact_mock)
    assert np.allclose(stiffness_to_test, -1.0 * np.ones_like(stiffness_to_test))
