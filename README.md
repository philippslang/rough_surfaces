## Analysis, elastic contact and fluid flow simulation for rock fractures

:construction: this is under heavy construction right now :construction:

[![Build Status](https://travis-ci.org/plang85/rough_surfaces.svg?branch=master)](https://travis-ci.org/plang85/rough_surfaces)


If you'd be willing to help out, first spread the word ;), then clone the repo and cd into it using a virtual environment

```
git clone https://github.com/plang85/rough_surfaces.git
python3 -m venv rsenv
. rsenv/bin/activate
```

Then verify that the tests are passing

```
cd rough_surfaces
pip install -e .[test] 
pytest
```

and do some damage. 

## Examples

See `example_analysis.py` for a more complete overview and the provided plotting functions.

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
