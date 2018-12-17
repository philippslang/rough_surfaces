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

N_power_of_two, dxy = 9, 1

fig, axes = plt.subplots(1, 2, figsize=(12, 6))

surface_params = rp.SelfAffineParameters()
surface_params.lambda_L_over_lambda_0 = 2
seed = 43

surface_isotropic = rg.make_self_affine(surface_params, N_power_of_two, seed=seed)
ax = axes[0]
ax.set_title('Isotropic') 
rplt.roughness(ax, surface_isotropic, dxy, 'L')

surface_params.anisotropy = 10
surface_anisotropic = rg.make_self_affine(surface_params, N_power_of_two, seed=seed)
ax = axes[1]
ax.set_title('Anisotropic') 
rplt.roughness(ax, surface_anisotropic, dxy, 'L')

plt.show()
