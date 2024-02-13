![Unit tests](https://github.com/acroucher/PyTOUGH/workflows/Unit%20tests/badge.svg)

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

For more detailed information on PyTOUGH, consult the [user guide](https://github.com/acroucher/PyTOUGH/blob/master/doc/PyTOUGH-guide.pdf) (PDF format) and the PyTOUGH [wiki](https://github.com/acroucher/PyTOUGH/wiki/), which has links to published articles on PyTOUGH.

# What's new in PyTOUGH?

The latest stable version is 1.5.7, which has:

* support for exporting source networks (e.g. reinjection) to Waiwera JSON input

* support for exporting unlimited number of time steps to Waiwera JSON input

* a bugfix for parsing TOUGHREACT-OMP listing files

