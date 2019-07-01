# MDInteract

[![Build Status](https://travis-ci.org/MD-Studio/MDInteract.svg?branch=master)](https://travis-ci.org/MD-Studio/MDInteract)
[![PyPI version](https://badge.fury.io/py/mdinteract.svg)](https://badge.fury.io/py/mdinteract)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/MD-Studio/MDInteract/master?filepath=examples)

Analyze biomolecular structures using the popular [Pandas](https://pandas.pydata.org) data analysis toolkit.

Using MDInteract you can:

- Directly analyze atomic structure data in popular file formats with full support for Molecular Dynamics 
  trajectory format with the help of [MDTraj](http://mdtraj.org/1.9.3/).
- Query and analyze structures data using the powerful Pandas query syntax and data analysis methods.
- Analyze specific interaction types such as charged interactions, hydrogen bonding, water bridging, halogen 
  interactions and hydrophobic effects among others.
- Calculate many geometric parameters (distance, angles, RMSD, vectors etc.)

### Try before you install

Get a feel of the basic functionality of MDInteract by running interactive demonstration Jupyter notebooks right
in your browser. Just click the "launch binder" button, no registration required!


### Installation

Install the stable release directly from the Python Package Index (PyPI)

    pip install mdinteract

Or, clone or download the cutting edge MDInteract repository and install using pip:

    pip install MDInteract/
