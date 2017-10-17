import pytest
import numpy as np
import brown.analyse as ba
import brown.params as bp
import brown.generate as bg
import brown.surface as bs


@pytest.fixture
def rough_surface():
    # TODO more hurst parameters
    # TODO discretization parameter abstraction
    N_power_of_two = 9
    surface_params = bp.SelfAffineParameters()
    surface = bg.self_affine(surface_params, N_power_of_two)
    return surface, surface_params, N_power_of_two


def test_self_affine(rough_surface):
    surface, surface_params, N_power_of_two = rough_surface
    N = 2**N_power_of_two
    assert surface.shape[0] == surface.shape[1] == N
    assert np.isclose(bs.rms(surface), bp.SelfAffineParameters().hrms)
    # TODO this should be able to take the surface directly, make it numpy array compatible
    surface_spectrum = ba.radially_averaged_psd(surface)
    surface_invariants = ba.self_affine_psd_fit(*surface_spectrum)
    surface_hurst = surface_invariants[1]
    assert np.isclose(surface_params.hurst, surface_hurst, rtol=0.1)
