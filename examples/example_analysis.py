import numpy as np
import surfplot as sp
import surfanalysis as sa
import matplotlib.pyplot as plt


if __name__ == '__main__':
    fnames = ['isotropic.txt', 'anisotropic.txt']
    dxys = [5./2**9 for i in range(len(fnames))] # sample spacing
    window = True
    length_unit = 'm'

    for dxy, fname in zip(dxys, fnames):
        h = np.loadtxt(fname)        

        fig, axes = plt.subplots(2, 2, figsize=(14,9))
        xtext, ytext = 0.65, 0.9

        ax = axes[0,0]
        ax.set_title('Roughness Profile') 
        sp.roughness(ax, h, dxy, length_unit)

        ax = axes[0,1]
        ax.set_title('Radially Averaged PSD') 
        f, C = sa.radially_averaged_psd(h, dxy, window)
        pref, hurst = sa.self_affine_psd_fit(f, C)
        sp.roughness_spectrum(ax, f, C, length_unit)
        sp.roughness_spectrum(ax, f, sa.self_affine_psd(f, pref, hurst), length_unit)
        ax.text(xtext, ytext, ''.join(['$\mathit{H}$',' = {:.2f}'.format(hurst)]), transform=ax.transAxes, va='top')

        ax = axes[1,0]
        ax.set_title('Axially Averaged PSDs') 
        fx, Cx = sa.axis_averaged_psd(h, dxy, window, 0)
        fy, Cy = sa.axis_averaged_psd(h, dxy, window, 1)
        prefy, hursty = sa.self_affine_psd_fit(fy, Cy, onedim=True)
        prefx, hurstx = sa.self_affine_psd_fit(fx, Cx, onedim=True)
        sp.roughness_spectrum(ax, fx, Cx, length_unit, onedim=True)
        sp.roughness_spectrum(ax, fx, sa.self_affine_psd(fx, prefx, hurstx, True), length_unit, onedim=True)
        
        sp.roughness_spectrum(ax, fy, Cy, length_unit, onedim=True)
        sp.roughness_spectrum(ax, fy, sa.self_affine_psd(fy, prefy, hursty, True), length_unit, onedim=True)
        ax.text(xtext, ytext, ''.join(['$\mathit{H_x}$',' = {:.2f}'.format(hurstx), '\n', 
                                   '$\mathit{H_y}$',' = {:.2f}'.format(hursty), '\n', 
                                   '$\mathit{H_x/H_y}$',' = {:.2f}'.format(hurstx/hursty)]), va='top', transform=ax.transAxes)
      
        ax = axes[1,1]
        ax.set_title('Roughness Distribution') 
        sp.roughness_histogram(ax, h, length_unit)
        ax.text(xtext, ytext, ''.join(['$\mathit{\sigma}$',' = {:.1e} '.format(np.std(h)), length_unit]), transform=ax.transAxes, va='top') # biased estimator

        fname = fname.replace('.txt', '.png')
        plt.tight_layout()
        plt.savefig(fname)
        plt.show()