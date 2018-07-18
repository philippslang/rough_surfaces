import numpy as np
import rough_surfaces.surface as rs
import rough_surfaces.analyse as ra
import rough_surfaces.plot as rp
import matplotlib.pyplot as plt


fnames = ['../doc/isotropic.txt', '../doc/anisotropic.txt']
dxys = [5./2**9] * len(fnames) # sample spacing, here the same for both decks
window = True
length_unit = 'm'

for dxy, fname in zip(dxys, fnames):
    h = np.loadtxt(fname)        
    surface = rs.Surface(h, dxy)

    fig, axes = plt.subplots(2, 2, figsize=(14,9))
    xtext, ytext = 0.65, 0.9

    ax = axes[0,0]
    ax.set_title('Roughness Profile') 
    rp.roughness(ax, h, dxy, length_unit)

    ax = axes[0,1]
    ax.set_title('Radially Averaged PSD') 
    f, C = ra.radially_averaged_psd(surface, window)
    pref, hurst = ra.self_affine_psd_fit(f, C)
    rp.roughness_spectrum(ax, f, C, length_unit)
    rp.roughness_spectrum(ax, f, ra.self_affine_psd(f, pref, hurst), length_unit)
    ax.text(xtext, ytext, ''.join(['$\mathit{H}$',' = {:.2f}'.format(hurst)]), transform=ax.transAxes, va='top')

    ax = axes[1,0]
    ax.set_title('Axially Averaged PSDs') 
    fx, Cx = ra.axis_averaged_psd(h, dxy, window, 0)
    fy, Cy = ra.axis_averaged_psd(h, dxy, window, 1)
    prefy, hursty = ra.self_affine_psd_fit(fy, Cy, onedim=True)
    prefx, hurstx = ra.self_affine_psd_fit(fx, Cx, onedim=True)
    rp.roughness_spectrum(ax, fx, Cx, length_unit, onedim=True)
    rp.roughness_spectrum(ax, fx, ra.self_affine_psd(fx, prefx, hurstx, True), length_unit, onedim=True)
    
    rp.roughness_spectrum(ax, fy, Cy, length_unit, onedim=True)
    rp.roughness_spectrum(ax, fy, ra.self_affine_psd(fy, prefy, hursty, True), length_unit, onedim=True)
    ax.text(xtext, ytext, ''.join(['$\mathit{H_x}$',' = {:.2f}'.format(hurstx), '\n', 
                                '$\mathit{H_y}$',' = {:.2f}'.format(hursty), '\n', 
                                '$\mathit{H_x/H_y}$',' = {:.2f}'.format(hurstx/hursty)]), va='top', transform=ax.transAxes)
    
    ax = axes[1,1]
    ax.set_title('Roughness Distribution') 
    rp.roughness_histogram(ax, h, length_unit)
    ax.text(xtext, ytext, ''.join(['$\mathit{\sigma}$',' = {:.1e} '.format(np.std(h)), length_unit]), transform=ax.transAxes, va='top') # biased estimator

    fname = fname.replace('txt', 'png')
    plt.tight_layout()
    plt.savefig(fname)
    plt.show()