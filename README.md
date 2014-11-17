# What is PyTOUGH?

PyTOUGH (Python TOUGH) is a Python library for simplifying, extending and automating the use of the [TOUGH2](http://esd.lbl.gov/research/projects/tough/) subsurface fluid and heat flow simulator. Using PyTOUGH, it is possible to automate the creation and editing of TOUGH2 model grids and data files, and the analysis and display of model simulation results, using Python scripts.

# Installing PyTOUGH:

First, make sure you have [Python](http://www.python.org) (note: use Python 2.x- Python 3.x is not supported at present) and the [Numerical Python](http://numpy.scipy.org/) library installed on your machine.  (For some features you will need other libraries such as [Scientific Python](http://www.scipy.org/) or [Matplotlib](http://matplotlib.sourceforge.net/)- consult the user guide for details.)

Click the [Download ZIP](https://github.com/acroucher/PyTOUGH/archive/master.zip) button at the right of the PyTOUGH web page, and save the .zip file to your computer.  Unzip this to any directory on your computer.  This will create a directory containing a file called `setup.py`.  At the command line type `python setup.py install`.

# More information:

For more detailed information on PyTOUGH, consult the user guide (PDF format, in the 'doc' directory of your PyTOUGH install) and the PyTOUGH [wiki](https://github.com/acroucher/PyTOUGH/wiki/), which has links to published articles on PyTOUGH.

# What's new in PyTOUGH?

The latest stable version is 1.4.0, which has:

* support for EOS modules with more than four primary variables in t2data and t2incons

* support for TOUGH2 NSEQ, NADD parameters in initial conditions in t2data and t2incons

* new mulgrid decompose_columns() method, for subdividing columns with more than four sides into triangles and quadrilaterals

* new mulgrid fit_columns() method for fitting arbitrary (x,y,z) datasets to columns (not just surface elevations)

* fitting (x,y,z) datasets to mulgrids containing columns with more than four sides

* new t2grids rectgeo() method, for reverse-engineering rectangular mulgrid geometry from a TOUGH2 data file (e.g. for working with models created by PetraSim)

* block mappings, for using TOUGH2 meshes with naming conventions not supported by mulgrids

* new region() method in t2thermo and IAPWS97 modules- to determine thermodynamic region for temperature and pressure

* checking thermodynamic function input variables to see if they are within operating range

* extracting t2listing histories with datetimes

* exporting t2listing tables to Pandas DataFrames

* exporting mulgrid geometries to ExodusII mesh file format

as well as various bug fixes and other minor enhancements.

# Where's the user guide?

Since PyTOUGH version 1.3.6, the PyTOUGH user guide (PyTOUGH-guide.pdf) is now included in the 'doc' directory of your PyTOUGH install.  Previously this was available separately from the 'Downloads' section on the website, but GitHub decided to phase out this 'Downloads' section.

