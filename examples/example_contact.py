import numpy as np
import brown.params as bp
import brown.generate as bg
import brown.analyse as ba
import brown.contact as bc
import matplotlib.pyplot as plt
from matplotlib import rcParams
import brown.plot as bplt
# TODO dont do this
rcParams['savefig.dpi'] = 300
rcParams['legend.loc'] = 'upper right'
rcParams['image.cmap'] = 'hot'

N_power_of_two = 9
surface_params = bp.self_affine_default_parameters()
surface = bg.self_affine(surface_params, N_power_of_two, seed=0)
E, nu = 1.0E+9, 0.3
composite_modulus = bc.homogeneous_composite_modulus(E, nu)
dxy = 1.0E-3
nominal_stress = 1.0E7
contact = bc.contact_FFT(surface, nominal_stress, E, nu, verbose=1, err_lim=1.0E-8)


if 0:
    fig, ax = plt.subplots()
    bplt.traces(ax, surface.h, [contact.u], 128)
    unit_den = '(m)'
    ax.set_xlabel('x ' + unit_den)
    ax.set_ylabel('y ' + unit_den)
    plt.legend()
    plt.show()

if 1:
    fig, ax = plt.subplots()
    N = surface.shape[0]
    L = surface.length()
    x = np.linspace(-L/2., L/2., N)
    XX, YY = np.meshgrid(x, x)
    pressure_plot = ax.pcolor(XX, YY, contact.p) 
    ax.axis('equal')
    unit_den = '(m)'
    ax.set_xlabel('x ' + unit_den)
    ax.set_ylabel('y ' + unit_den)
    cbar = plt.colorbar(pressure_plot)
    cbar.set_label('Pressure (Pa)', rotation=270)
    plt.show()

# TODO cross section