# PLOT PARAMETERS MODULE

import itertools
import surfunits as su
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['mathtext.default'] = 'regular'
rcParams['font.sans-serif'] = ['Arial']
rcParams['axes.color_cycle'] = ['b', 'r', 'g', 'c', 'm', 'y', 'k']
rcParams['figure.figsize'] = 8,6
rcParams['font.size'] = 18
rcParams['legend.fontsize'] = 14
rcParams['lines.linewidth'] = 2
rcParams['savefig.dpi'] = 300
rcParams['legend.loc'] = 'upper right'
rcParams['image.cmap'] = 'hot'
# TODO
# http://matplotlib.org/users/style_sheets.html replace by style sheet
# 1.) all plots should setup plt contents and then send to a central print/plot function here hat takes care of formatting, saving etc...


# global plot parameters
__extension = '.png'


def legend_font_size():
    """We want the legend font size smaller."""
    return rcParams['font.size']-4


def plot_format_is_png():
    """Sets current figure output format to PNG."""
    global __extension
    __extension = '.png'


def plot_format_is_pdf():
    """Sets current figure output format to PDF."""
    global __extension
    __extension = '.pdf'

def plot_format_is_tif():
    """Sets current figure output format to TIFF."""
    global __extension
    __extension = '.tif'


def plot_format_is_pgf():
    """Sets current figure output format to LaTeX PGF Figure."""
    global __extension
    __extension = '.pgf'


def plot_format_is_eps():
    """Sets current figure output format to post-script."""
    global __extension
    __extension = '.pgf'


def file_name(fname):
    """Returns fname with current extension for output format."""
    global __extension
    return str(fname+__extension)


def length_label(label):
    """Returns label with unit appended (plus whitespace in between). Empty unit if settings not recognized."""
    if su.length_unit() == 0:
        return str( label + r' ($\mu$m)' )
    if su.length_unit() == 1:
        return str( label + r' (mm)' )
    if su.length_unit() == 2:
        return str( label + r' (cm)' )
    if su.length_unit() == 3:
        return str( label + r' (m)' )
    return str( label + ' ()' )

def stiffness_label(label):
    """Returns label with unit appended (plus whitespace in between). Empty unit if settings not recognized."""
    if su.length_unit() == 0:
        return str( label + r' (Pa $\mu$m$^{-1}$)' )
    if su.length_unit() == 1:
        return str( label + r' (Pa mm$^{-1}$)' )
    if su.length_unit() == 2:
        return str( label + r' (Pa cm$^{-1}$)' )
    if su.length_unit() == 3:
        return str( label + r' (Pa m$^{-1}$)' )
    return str( label + ' ()' )


def area_label(label):
    """Returns label with unit appended (plus whitespace in between). Empty unit if settings not recognized."""
    if su.length_unit() == 0:
        return str( label + r' ($\mu$m$^2$)' )
    if su.length_unit() == 1:
        return str( label + r' (mm$^2$)' )
    if su.length_unit() == 2:
        return str( label + r' (cm$^2$)' )
    if su.length_unit() == 3:
        return str( label + r' (m$^2$)' )
    return str( label + ' []' )


def frequency_label(label):
    """Returns label with unit appended (plus whitespace in between). Empty unit if settings not recognized."""
    if su.length_unit() == 0:
        return str( label + r' ($\mu$m$^{-1}$)' )
    if su.length_unit() == 1:
        return str( label + r' (mm$^{-1}$)' )
    if su.length_unit() == 2:
        return str( label + r' (cm$^{-1}$)' )
    if su.length_unit() == 3:
        return str( label + r' (m$^{-1}$)' )
    return str( label + ' ()' )


def power_label2D(label):
    """Returns label with unit appended (plus whitespace in between). Empty unit if settings not recognized."""
    if su.length_unit() == 0:
        return str( label + r' ($\mu$m$^{4}$)' )
    if su.length_unit() == 1:
        return str( label + r' (mm$^{4}$)' )
    if su.length_unit() == 2:
        return str( label + r' (cm$^{4}$)' )
    if su.length_unit() == 3:
        return str( label + r' (m$^{4}$)' )
    return str( label + ' ()' )


def power_label1D(label):
    """Returns label with unit appended (plus whitespace in between). Empty unit if settings not recognized."""
    if su.length_unit() == 0:
        return str( label + r' ($\mu$m$^{3}$)' )
    if su.length_unit() == 1:
        return str( label + r' (mm$^{3}$)' )
    if su.length_unit() == 2:
        return str( label + r' (cm$^{3}$)' )
    if su.length_unit() == 3:
        return str( label + r' (m$^{3}$)' )
    return str( label + ' ()' )


def scientific(x, pos):
    """The two args are the value and tick position"""
    return '%.1e' % x


def ccycle():
    return itertools.cycle(['#0000FF', '#FF0000', '#00CC00', '#000000', '#FF9900', '#787878', '#990099'])

def mecycle():
    return itertools.cycle(['x', '+', '1', '.', '2'])

def mfcycle():
    return itertools.cycle(['^', 's', 'o', 'D', 'v'])