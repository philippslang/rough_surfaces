import numpy as np
import rough_surfaces.params as rp
import rough_surfaces.generate as rg
import rough_surfaces.analyse as ra
import rough_surfaces.surface as rs
import rough_surfaces.contact as rc
import matplotlib.pyplot as plt
from matplotlib import rcParams
import rough_surfaces.plot as rplt
# TODO dont do this
rcParams['savefig.dpi'] = 300
rcParams['legend.loc'] = 'upper right'
rcParams['image.cmap'] = 'hot'

N_power_of_two = 9
surface_params = rp.SelfAffineParameters()
surface = rg.make_self_affine(surface_params, N_power_of_two, seed=0)
E, nu = 1.0E+9, 0.3
dxy = 1.0E-3
nominal_stress = 1.0E7
contact = rc.contact_FFT(surface, nominal_stress, E, nu, verbose=2, err_lim=1.0E-8)


if 1:
    fig, ax = plt.subplots()
    rplt.traces(ax, surface, [contact.u], 128)
    unit_den = '(m)'
    ax.set_xlabel('x ' + unit_den)
    ax.set_ylabel('y ' + unit_den)
    plt.legend()
    plt.show()

if 1:
    fig, ax = plt.subplots()
    N = surface.shape[0]
    L = rs.length(surface)
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
