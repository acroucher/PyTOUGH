# What is PyTOUGH?

PyTOUGH (Python TOUGH) is a Python library for simplifying, extending and automating the use of the [TOUGH2](http://esd.lbl.gov/research/projects/tough/) subsurface fluid and heat flow simulator. Using PyTOUGH, it is possible to automate the creation and editing of TOUGH2 model grids and data files, and the analysis and display of model simulation results, using Python scripts.

# Installing PyTOUGH:

First, make sure you have [Python](http://www.python.org) (note: use Python 2.x- Python 3.x is not supported at present) and the [Numerical Python](http://numpy.scipy.org/) library installed on your machine.  (For some features you will need other libraries such as [Scientific Python](http://www.scipy.org/) or [Matplotlib](http://matplotlib.sourceforge.net/)- consult the user guide for details.)

Click the [Download ZIP](https://github.com/acroucher/PyTOUGH/archive/master.zip) button at the right of the PyTOUGH web page, and save the .zip file to your computer.  Unzip this to any directory on your computer.  This will create a directory containing a file called `setup.py`.  At the command line type `python setup.py install`.

# More information:

For more detailed information on PyTOUGH, consult the user guide (PDF format, in the 'doc' directory of your PyTOUGH install) and the PyTOUGH [wiki](https://github.com/acroucher/PyTOUGH/wiki/), which has links to published articles on PyTOUGH.

# What's new in PyTOUGH?

The latest stable version is 1.4.2, which has:

* A fix for a bug in the t2grid minc() method, caused by basing the code on a version of GMINC which had the same bug- block order was reversed in MINC connections

* MINC processing ability via the new t2grid minc() method, including partial MINC (applied to only part of grid), flexible block and rocktype naming options, and initial conditions handling

* ability to extract results over a mulgrid geometry for a particular MINC level via new minc_array() method

* ability to plot block names in mulgrid layer_plot(), and support for blockmaps

* better support for mulgrid ExodusII export

* remaining print statements in t2listing.py turned into Python exceptions

as well as various bug fixes and other minor enhancements.

# Where's the user guide?

Since PyTOUGH version 1.3.6, the PyTOUGH user guide (PyTOUGH-guide.pdf) is now included in the 'doc' directory of your PyTOUGH install.  Previously this was available separately from the 'Downloads' section on the website, but GitHub decided to phase out this 'Downloads' section.

