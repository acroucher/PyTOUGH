![Unit tests](https://github.com/acroucher/PyTOUGH/workflows/Unit%20tests/badge.svg) [![Documentation Status](https://readthedocs.org/projects/PyTOUGH/badge/?version=latest)](https://PyTOUGH.readthedocs.io/en/latest/?badge=latest) [![PyPI version](https://badge.fury.io/py/PyTOUGH.svg)](https://badge.fury.io/py/PyTOUGH) ![PyPI - Downloads](https://img.shields.io/pypi/dm/PyTOUGH)


# What is PyTOUGH?

PyTOUGH (Python TOUGH) is a Python library for simplifying, extending and automating the use of the [TOUGH2](http://esd.lbl.gov/research/projects/tough/) subsurface fluid and heat flow simulator. Using PyTOUGH, it is possible to automate the creation and editing of TOUGH2 model grids and data files, and the analysis and display of model simulation results, using Python scripts.

# Installing PyTOUGH:

From version 1.6.0, PyTOUGH can be installed via the [`pip`](https://pip.pypa.io) Python package installer:

```
pip install PyTOUGH
```
You can also install a particular version of PyTOUGH, e.g. to install version 1.6.0:

```
pip install PyTOUGH==1.6.0
```
To uninstall PyTOUGH:
```
pip uninstall PyTOUGH
```
To install the `testing` branch, to get the most recent changes being tested for the next stable release:
```
pip install git+https://github.com/acroucher/PyTOUGH.git@testing
```

# More information:

For more detailed information on PyTOUGH, consult the [user guide](https://pytough.readthedocs.io) (html, or you can download PDF or Epub versions) and the PyTOUGH [wiki](https://github.com/acroucher/PyTOUGH/wiki/), which has links to published articles on PyTOUGH.

# What's new in PyTOUGH?

The latest stable version is 1.6.2, which has:

* a new `mulgrid` naming convention 3, with 3 characters for columns and 2 characters for layers

* the `mulgrid` `column_track()` method has been re-written with a different algorithm, enabling it (and the `slice_plot()` method) to handle non-contiguous slices, e.g. if the slice passes out of the mesh and back in
