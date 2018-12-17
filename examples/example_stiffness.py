import numpy as np
import rough_surfaces.params as rp
import rough_surfaces.generate as rg
import rough_surfaces.analyse as ra
import rough_surfaces.contact as rc
import matplotlib.pyplot as plt
import rough_surfaces.plot as rplt


if 1:
    N_power_of_two = 9
    surface_params = rp.SelfAffineParameters()
    surface_params.hrms = 0.005
    surface = rg.make_self_affine(surface_params, N_power_of_two, seed=0)
    E, nu = 1.0E+9, 0.3
    dxy = 1.0E-3
    nominal_stress = np.logspace(6, 8, 15)
    stiffness = rc.stiffness(nominal_stress, surface, E, nu, err_lim=1.0E-8)
    np.savetxt('stiffness', stiffness)
else:
    stiffness = np.loadtxt('stiffness')
    
if 1:
    fig, ax = plt.subplots()
    stress_at_stiffness = nominal_stress[:-1] + np.diff(nominal_stress)
    ax.semilogy(stress_at_stiffness, stiffness, 'x')
    ax.set_xlabel('Pressure (Pa)')
    ax.set_ylabel('$\kappa$ (Pa m$^{-1}$)')
    plt.show()
