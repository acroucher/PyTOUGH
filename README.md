![Unit tests](https://github.com/acroucher/PyTOUGH/workflows/Unit%20tests/badge.svg)

# What is PyTOUGH?

PyTOUGH (Python TOUGH) is a Python library for simplifying, extending and automating the use of the [TOUGH2](http://esd.lbl.gov/research/projects/tough/) subsurface fluid and heat flow simulator. Using PyTOUGH, it is possible to automate the creation and editing of TOUGH2 model grids and data files, and the analysis and display of model simulation results, using Python scripts.

# Installing PyTOUGH:

First, make sure you have [Python](http://www.python.org) and the [Numerical Python](http://numpy.scipy.org/) library installed on your machine.  (For some features you will need other libraries such as [Scientific Python](http://www.scipy.org/) or [Matplotlib](http://matplotlib.sourceforge.net/)- consult the user guide for details.)

Click the _Clone or download_ button at the right of the PyTOUGH web page and then click [Download ZIP](https://github.com/acroucher/PyTOUGH/archive/master.zip) , and save the .zip file to your computer.  Unzip this to any directory on your computer.  This will create a directory containing a file called `setup.py`.  At the command line type `python setup.py install`.

(Alternatively, if you are confident using the Git version control system, you can clone the PyTOUGH repository instead of downloading a .zip file.)

# More information:

For more detailed information on PyTOUGH, consult the user guide (PDF format, in the 'doc' directory of your PyTOUGH install) and the PyTOUGH [wiki](https://github.com/acroucher/PyTOUGH/wiki/), which has links to published articles on PyTOUGH.

# What's new in PyTOUGH?

The latest stable version is 1.5.6, which has:

* bugfixes and support for GMSH 4.1 format in mulgrid from_gmsh() import
* ability to export mulgrid geometries to [Layermesh](https://github.com/acroucher/layermesh) objects
* support for tracers in t2data json() method, for exporting models to [Waiwera](https://waiwera.github.io/)
* griddata function is now imported from SciPy rather than Matplotlib

# Where's the user guide?

Since PyTOUGH version 1.3.6, the PyTOUGH user guide (PyTOUGH-guide.pdf) is now included in the 'doc' directory of your PyTOUGH install.  Previously this was available separately from the 'Downloads' section on the website, but GitHub decided to phase out this 'Downloads' section.

