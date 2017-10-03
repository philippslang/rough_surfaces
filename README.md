# Brown

A python3 module for the analysis, elastic contact and fluid flow simulation of rock fractures

Platform | Build status
---------|-------------:
Linux, Python 3.6 | [![Build Status](https://travis-ci.org/plang85/rough_surfaces.svg?branch=master)](https://travis-ci.org/plang85/rough_surfaces)

## Install/Contribute

If you're willing to help out, first spread the word ;), then clone the repo and cd into it using a virtual environment
```
git clone https://github.com/plang85/rough_surfaces.git
python3 -m venv rsenv
. rsenv/bin/activate
```
and verify that the tests are passing
```
cd rough_surfaces
pip install -e .[test] 
pytest
```
do some damage and just PR against master. 

To just use the module, do
```
cd rough_surfaces
pip install .
```
instead.

## Examples

### Generation

We can generate an isotropic, self-affine surface like this
```
import brown.params as bp
import brown.generate as bg

surface_params = bp.self_affine_default_parameters()
N_power_of_two = 9
surface = bg.self_affine(surface_params, N_power_of_two)
```
where `surface` is a two-dimensional numpy array duck.

### Analysis

See `example_analysis.py` for a more complete overview and the provided plotting functions. In short, given a two-dimensional array `h` that represents discrete surface height uniformly spaced by `dxy`, the radially averaged power-spectrum can be obtained like this (this API is in flux)
```
import brown.analyse as ba

spectrum = ba.radially_averaged_psd(h, dxy)
invariants = ba.self_affine_psd_fit(*surface_spectrum)
print('Hurst = {0:.2f}'.format(invariants[1]))
```

An isotropic surface is characterized by a near-ideal straight line radially averaged PSD and same-slope curves for the angularly averaged spectra.
<p align="left">
  <img src="https://raw.githubusercontent.com/plang85/rough_surfaces/master/doc/isotropic.png" height="400">
  <br/>
</p>
An anisotropic surface is characterized by a less linear radially averaged PSD and different-slope curves for the angularly averaged spectra.
<p align="left">
  <img src="https://raw.githubusercontent.com/plang85/rough_surfaces/master/doc/anisotropic.png" height="400">
  <br/>
</p>

### Contact (elastic frictionless)

We can solve the elastic frictionless contact between two rough surfaces by solving the equivalent problem of a rigid composite surface against an elastic, flat body of composite properties
```
import brown.contact as bc

dxy = 1.0E-3
nominal_stress = 1.0E7
contact = bc.contact_FFT(surface, nominal_stress, E, nu, verbose=1)
```
For a more detailed snippet see `example_analysis.py`.

<p align="left">
  <img src="https://raw.githubusercontent.com/plang85/rough_surfaces/master/doc/contact.png" height="400">
  <br/>
</p>

### Flow (steady-state laminar flow fluid pressure)

Coming soon...