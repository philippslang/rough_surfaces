import numpy as np
from matplotlib import rcParams
import scipy.stats as scst
rcParams['font.family'] = 'sans-serif'
rcParams['mathtext.default'] = 'regular'
rcParams['font.sans-serif'] = ['Arial']
rcParams['font.size'] = 14
rcParams['legend.fontsize'] = 10
rcParams['axes.color_cycle'] = ['black', 'red', 'blue']
rcParams['savefig.dpi'] = 300
rcParams['legend.loc'] = 'upper right'
rcParams['image.cmap'] = 'hot'


def roughness_spectrum(ax, q, C, lenght_unit, onedim=False):
    ax.loglog(q, C)
    ax.set_xlabel('q (' + lenght_unit +  '$^{-1}$' + ')')
    ax.set_ylabel('C (' + lenght_unit +  '$^{4}$' + ')')
    if onedim:
        ax.set_ylabel('C (' + lenght_unit +  '$^{3}$' + ')')


def roughness(ax, h, dxy, length_unit):
    N = h.shape[0]
    L = dxy*N
    x = np.linspace(-L/2., L/2., N)
    XX,YY = np.meshgrid(x,x)
    ax.pcolor(XX, YY, h)
    ax.axis('equal')
    unit_den = ''.join(['(', length_unit, ')'])
    ax.set_xlabel('x ' + unit_den)
    ax.set_ylabel('y ' + unit_den)


def roughness_histogram(ax, h, length_unit):
    bins = 30
    ax.hist(h.flatten(), bins, normed=1, color='gray', ec='white')
    hsigma = np.std(h)
    hspace = np.linspace(h.min(), h.max(), 100)
    ax.plot(hspace, scst.norm.pdf(hspace, np.mean(h), hsigma), lw=3, alpha=0.5)
    unit_den = ''.join(['(', length_unit, ')'])
    ax.set_xlabel('Surface Height ' + unit_den)
    ax.set_ylabel('Relative Probability')


def slope_histogram(ax, h):
    bins = 30
    g = np.gradient(h)
    ax.hist(np.ravel(g), bins, normed=1, color='gray', ec='white')
    #hsigma = np.std(h)
    #hspace = np.linspace(h.min(), h.max(), 100)
    #ax.plot(hspace, scst.norm.pdf(hspace, np.mean(h), hsigma), lw=3, alpha=0.5)
    ax.set_xlabel('Surface Slope (-)')
    ax.set_ylabel('Relative Probability')