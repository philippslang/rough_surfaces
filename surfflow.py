# SURFACE FLOW SIMULATION MODULE

import numpy as np
import pickle as pc
from fipy import *
import warnings
import matplotlib.pyplot as plt
import surfplot as splt
import surfunits as su
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import fipy 

def water_viscosity(T=273.15):
    """Returns [Pa s] from [K]."""
    return 2.414E-5*10**(247.8/(T-140.))


def flow_rate(T, press_diff, W, L, ah):
    """Returns transport volumetric flow rate [m3 s-1]. T [K]; press_diff [Pa]; W, L, ah [m]."""
    mu = water_viscosity(T)
    return ah**3 * W * press_diff / (12.0 * mu * L) # [m3 s-1] 


def differential_pressure(T, Q, W, L, ah):
    """Returns pressure difference [Pa] over length [L]. T [K]; W, L, ah [m]."""
    mu = water_viscosity(T)
    return Q*L*mu*12./(W*ah**3) # [m3 s-1]


def interstitial_velocity(Q, W, am):
    """Returns transport velocity [m s-1]. Q [m3 s-1]; W, am [m]."""
    return Q/(W*am) # [m s-1] 


class Results:
    """ 
    A container and analysis class for 2D flow solutions.
    
    Attributes
    ----------
    dxy:
        Isotropic grid spacing.
    p: array_like
        Fluid pressure, 2D numpy array.
    a: float
        Effective hydraulic aperture [L].
    v_x, v_y: ndarray
        Fluid velocity components, 2D numpy array.

    Notes
    -----
    The term 'pressure' herein refers to fluid pressure in cells.
    """
    
    def __init__(self, dxy, p, a=None, v_x=None, v_y=None):
        """ 
        Initializes result containers.

        Parameters
        ----------
        dxy:
            Isotropic grid spacing.
        p: array_like
            Fluid pressure, 2D numpy array.
        a: float
            Effective hydraulic aperture [L].
        v_x, v_y: ndarray
            Fluid velocity components, 2D numpy array.
        """
        self.dxy = dxy
        self.p = p
        self.a = a
        self.v_x = v_x
        self.v_y = v_y

    def __repr__(self):
        return "a_h = %f"%self.a


def save(flow, fname):
    """Saves Results instance to file."""
    f = open(fname, 'wb')
    pc.dump(flow, f)
    f.close()


def load(fname):
    """Loads Results instance from file."""
    f = open(fname, 'rb')
    results = pc.load(f)
    f.close()
    return results


def hydraulic_aperture(a_field, dxy, w=1.0E-9, verbose=False, y_only=False, indent=False):
    """
    Solves the 2D Laplace problem on the aperture profile of the fracture surface, returning the x,y axis effective permeability runs (local cubic law approximation).

    Parameters
    ----------
    a_field: 2D numpy array [L]
        Cell-wise aperture values over a rectangular, uniform grid. Will be converted to [m] internally.
    dxy: float
        Cell size in x and y direction.
    w: float, optional [m]
        Grain contact interface thickness, defaults to 1.0E-8 [m] and represents minimum value imposed on a_field.

    Returns
    -------
    results_x, results_y:
        Tuple of results instances along x- and y-axis, respectively.
    """
    dp = 2.
    nxy = a_field.shape[0]
    L = su.to_meter(dxy) * nxy
    mesh = Grid2D(dx=dxy, dy=dxy, nx=nxy, ny=nxy)
    a_field = su.to_meter(a_field)
    a_field[a_field < w] = w
    k_field = np.divide( np.power(a_field,2), 12. )
    k = CellVariable(name="permeability", mesh=mesh) 
    k.setValue(np.ravel(k_field))
    t_field = np.divide( np.power(a_field,3), 12. )
    t = CellVariable(name="transmissivity", mesh=mesh) 
    t.setValue(np.ravel(t_field))
    a = CellVariable(name="aperture", mesh=mesh) 
    a.setValue(np.ravel(a_field))
    xcc, ycc = mesh.cellCenters
    maxiter = 10000

    p_in = dp/2.
    p_out = -dp/2.

    eq = DiffusionTerm(coeff=t)
    p_y = CellVariable(name="pressure flow y", mesh=mesh)
    p_y.constrain(p_in, mesh.facesBottom)
    p_y.constrain(p_out, mesh.facesTop)
    if verbose:
        if indent:
            print '\tsolving flow problem in y...'
        else:
            print 'solving flow problem in y...'
    eq.solve(var=p_y, solver=LinearPCGSolver(iterations=maxiter))
    q_y = -t*p_y.grad*dxy
    q_yyv = q_y.value[1]
    q_yyv_y0 = q_yyv[ycc<(np.min(ycc)+dxy/4.)]
    Q_y_y0 = np.sum(q_yyv_y0)
    J_y = Q_y_y0/L
    T_y = J_y / (dp/L)
    a_y = np.power(T_y*12., 1./3.)
    a_y = su.from_meter(a_y)
    if y_only:
        return Results(dxy, np.reshape(p_y.value, a_field.shape), a_y, np.reshape(q_y.value[0], a_field.shape), np.reshape(q_y.value[1], a_field.shape))

    eq = DiffusionTerm(coeff=t)
    p_x = CellVariable(name="pressure flow x", mesh=mesh)
    p_x.constrain(p_in, mesh.facesLeft)
    p_x.constrain(p_out, mesh.facesRight)
    if verbose:
        if indent:
            print '\tsolving flow problem in x...'
        else:
            print 'solving flow problem in x...'
    eq.solve(var=p_x, solver=LinearPCGSolver(iterations=maxiter))
    q_x = -t*p_x.grad*dxy
    q_xxv = q_x.value[0]
    q_xxv_x0 = q_xxv[xcc<(np.min(xcc)+dxy/4.)]
    Q_x_x0 = np.sum(q_xxv_x0)
    J_x = Q_x_x0/L
    T_x = J_x / (dp/L)
    a_x = np.power(T_x*12., 1./3.)
    a_x = su.from_meter(a_x)    

    results_x = Results(dxy, np.reshape(p_x.value, a_field.shape), a_x, np.reshape(q_x.value[0], a_field.shape), np.reshape(q_x.value[1], a_field.shape))
    results_y = Results(dxy, np.reshape(p_y.value, a_field.shape), a_y, np.reshape(q_y.value[0], a_field.shape), np.reshape(q_y.value[1], a_field.shape))
    return results_x, results_y


def isotropic_hydraulic_aperture(a_field, dxy, w=1.0E-8):
    results_x, results_y = hydraulic_aperture(a_field, dxy, w)
    return np.mean([results_x.a,results_y.a])


def plot_pressure_field(p, dxy, show=True, save_as='', vmin=None, vmax=None):
    """
    Two-dimensional plot of pressure distribution field and colorbar. Assumes regular grid and square domain.

    Parameters
    ----------
    p: ndarray
        A 2D numpy array containing pressure values [Pa].
    dxy: float [L]
        Cell length of  domain in real units (as current setting in surfunits).
    show: boolean, optional
        If True, shows on screen, supressed otherwise (but can still be saved to file).
    save_as: string, opional
        If non-empty, saves image to file using the content. Must not inlcude extension, as this is retrieved from the plot module.
    vmin, vmax: optional
        Minimum and maximum values for data range (colorplot/colorbar range).
    """
    plt.clf()
    nxy = len(p[0])
    edge_length = dxy*nxy
    x = np.linspace( -edge_length/2., edge_length/2., nxy )
    XX,YY = np.meshgrid(x,x)
    M = np.zeros(p.shape, dtype='bool')
    M[p==0.] = True
    p = np.ma.masked_array(p, mask=M)
    if vmin == None and vmax == None:
        plt.pcolor(XX, YY, p)
    else:
        plt.pcolor( XX, YY, p,vmin=vmin,vmax=vmax)
    cb = plt.colorbar()
    cb.set_label( 'Pressure [Pa]' )
    plt.axis('equal')
    x_label = splt.length_label('Lateral Extension x')
    y_label = splt.length_label('Lateral Extension y')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    majorLocator = MultipleLocator(200000)
    plt.gca().xaxis.set_major_locator(majorLocator)
    plt.gca().yaxis.set_major_locator(majorLocator)
    plt.tight_layout()
    if len(save_as) != 0:
        fname = splt.file_name(save_as)
        plt.savefig(fname)
    if show:
        plt.show()
    plt.close()


def plot_velocity_field( v_x, v_y, dxy, show = True, save_as = '', magn_threshold = None ):
    """
    Two-dimensional plot of velocity vector field. Assumes regular grid and square domain.

    Parameters
    ----------
    v_x, v_y: ndarray
        A 2D numpy array containing velocity component values.
    dxy: float [L]
        Cell length of  domain in real units (as current setting in surfunits).
    show: boolean, optional
        If True, shows on screen, supressed otherwise (but can still be saved to file).
    save_as: string, opional
        If non-empty, saves image to file using the content. Must not inlcude extension, as this is retrieved from the plot module.
    magn_threshold: optional, float
        Minimum magnitude of velocity to be displayed in form of fraction (0..1).
    """
    plt.clf()
    nxy = len(v_x[0])
    edge_length = dxy*nxy
    x = np.linspace( -edge_length/2., edge_length/2., nxy )
    XX,YY = np.meshgrid(x,x)    
    if magn_threshold != None:
        M = np.zeros(v_x.shape, dtype='bool')
        magn = np.sqrt(v_x**2+v_y**2)
        magn_max = np.max(magn)
        magn_threshold = magn_max*magn_threshold
        M[magn<magn_threshold] = True
        v_x = np.ma.masked_array( v_x, mask=M )
        v_y = np.ma.masked_array( v_y, mask=M )
    plt.quiver( XX, YY, v_x, v_y ) 
    plt.axis('equal')
    x_label = splt.length_label('Lateral Extension x')
    y_label = splt.length_label('Lateral Extension y')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.tight_layout()
    if len(save_as) != 0:
        fname = splt.file_name(save_as)
        plt.savefig(fname)
    if show:
        plt.show()
    plt.close()


def plot_contact_velocity_field(p, v_x, v_y, dxy, show=True, save_as='', magn_threshold=None, a=None, nc=4, glyph_scale=None, 
                                every=1, pmin=None, pmax=None, peff=True):
    """
    Two-dimensional plot of contact zones and velocity field. Assumes regular grid and square domain.

    Parameters
    ----------
    p: ndarray
        A 2D numpy array containing contact pressure component values.
    v_x, v_y: ndarray
        A 2D numpy array containing velocity component values.
    dxy: float [L]
        Cell length of  domain in real units (as current setting in surfunits).
    show: boolean, optional
        If True, shows on screen, supressed otherwise (but can still be saved to file).
    save_as: string, opional
        If non-empty, saves image to file using the content. Must not inlcude extension, as this is retrieved from the plot module.
    magn_threshold: optional, float
        Minimum magnitude of velocity to be displayed in form of fraction (0..1).
    a: optional, ndarrayy
        The 2D numpy array of aperture values. If provided, overlays aperture contour plot.
    nc: int, optional
        Number of aperture contour lines.
    glyph_scale: float, opional
        Forward to scale in quiver plot.
    every: integer
        Plots every other arrow in glyph, defaults to 1 (all glyphs).
    """
    plt.clf()
    nxy = len(v_x[0])
    edge_length = dxy*nxy
    x = np.linspace(-edge_length/2., edge_length/2., nxy)
    XX,YY = np.meshgrid(x,x)    
    if magn_threshold != None:
        M = np.zeros(v_x.shape, dtype='bool')
        magn = np.sqrt(v_x**2+v_y**2)
        magn_max = np.max(magn)
        magn_threshold = magn_max*magn_threshold
        M[magn<magn_threshold] = True
        v_x = np.ma.masked_array(v_x, mask=M)
        v_y = np.ma.masked_array(v_y, mask=M)
    else:
        M = np.ones(v_x.shape, dtype='bool')
        M[p==0.] = False
        v_x = np.ma.masked_array(v_x, mask=M)
        v_y = np.ma.masked_array(v_y, mask=M)
    plt.quiver(XX[::every, ::every], YY[::every, ::every], v_x[::every, ::every], v_y[::every, ::every], zorder=0, scale=glyph_scale)
    M = np.zeros(p.shape, dtype='bool')
    M[p == 0.] = True
    p = np.ma.masked_array(p, mask=M)
    if pmin == None and pmax == None:
        plt.pcolor(XX, YY, p)
    else:
        plt.pcolor(XX, YY, p, vmin=pmin, vmax=pmax)
    cb = plt.colorbar()
    if peff:
        cb.set_label(r'$\mathit{\sigma^\prime}$ (Pa)')
    else:
        cb.set_label(r'$\mathit{\sigma}$ (Pa)')
    if a != None:
        ac = plt.contour(XX, YY, a, nc)
    plt.axis('equal')
    x_label = splt.length_label(r'$\mathit{x}$')
    y_label = splt.length_label(r'$\mathit{y}$')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.tight_layout()
    if len(save_as) != 0:
        fname = splt.file_name(save_as)
        plt.savefig(fname)
    if show:
        plt.show()
    plt.close()


def plot_contact_velocity_field_dualmesh(p, pdxy, v_x, v_y, vdxy, show=True, save_as='', skinny=False, every=1, figsize=None, glyph_scale=1, cbaror='horizontal', hook=None):
    """
    Two-dimensional plot of contact zones and velocity field. Assumes regular grid and square domain, works with different
    meshes for contact and flow.

    Parameters
    ----------
    p: ndarray
        A 2D numpy array containing contact pressure component values.
    v_x, v_y: ndarray
        A 2D numpy array containing velocity component values.
    dxy: float [L]
        Cell length of  domain in real units (as current setting in surfunits).
    show: boolean, optional
        If True, shows on screen, supressed otherwise (but can still be saved to file).
    save_as: string, opional
        If non-empty, saves image to file using the content. Must not inlcude extension, as this is retrieved from the plot module.
    """
    # flow
    fig = plt.figure(figsize=figsize)
    nxy = len(v_x[0])
    edge_length = vdxy*nxy
    x = np.linspace(-edge_length/2., edge_length/2., nxy)
    XX,YY = np.meshgrid(x,x) 
    if every == 1:   
        plt.quiver(XX, YY, v_x, v_y, zorder=0)
    else:
        plt.quiver(XX[::every,::every], YY[::every,::every], v_x[::every,::every], v_y[::every,::every], zorder=0, scale=glyph_scale)
    # contact
    nxy = p.shape[0]
    edge_length = pdxy*nxy
    x = np.linspace(-edge_length/2., edge_length/2., nxy)
    XX,YY = np.meshgrid(x,x)   
    M = np.zeros(p.shape, dtype='bool')
    M[p == 0.] = True
    p = np.ma.masked_array(p, mask=M)
    plt.pcolor(XX, YY, p, zorder=1)
    plt.axis('equal')
    if not skinny:
        cb = plt.colorbar(orientation=cbaror)
        cb.set_label(r'$\mathit{\sigma}$ (Pa)')
        plt.xlabel(splt.length_label('$\mathit{x}$'))
        plt.ylabel(splt.length_label('$\mathit{y}$'))
        plt.tight_layout()
    else:
        plt.axis('off') 
    if hook:
        hook(plt)  
    if len(save_as) != 0:
        fname = splt.file_name(save_as)
        if not skinny:
            plt.savefig(fname)
        else:
            plt.savefig(fname, bbox_inches='tight', pad_inches=0, dpi=150)
    if show:
        plt.show()
    plt.close()