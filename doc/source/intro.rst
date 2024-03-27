:tocdepth: 3

Introduction
============

What is PyTOUGH?
----------------

.. index:: PyTOUGH

PyTOUGH (**Py**\ thon **TOUGH**) is a set of Python software routines
for making it easier to use the TOUGH2 geothermal reservoir simulator.
Using PyTOUGH, it is possible to automate the creation and editing of
TOUGH2 model grids and data files, and the analysis and display of model
simulation results.

What are TOUGH2 and AUTOUGH2?
-----------------------------
.. index:: TOUGH2
.. index:: TOUGH2; AUTOUGH2

`TOUGH2 <https://tough.lbl.gov/>`_ is a general-purpose simulator for
modelling subsurface fluid and heat flow, often used for simulating
geothermal reservoirs.

AUTOUGH2 is the `University of Auckland
<https://www.geothermal.auckland.ac.nz/>`_ version of TOUGH2. The main
differences between AUTOUGH2 and TOUGH2 are:

-  **EOS handling**: AUTOUGH2 includes all different equations of state
   (EOSes) in a single executable program, whereas TOUGH2 uses different
   executables for each EOS. As a result, the main input data file for
   an AUTOUGH2 simulation also includes extra data blocks to specify
   which EOS is to be used.

-  **Generator types**: AUTOUGH2 includes a variety of extra generator
   types developed for geothermal reservoir simulation (e.g. makeup and
   reinjection wells).

.. index:: TOUGH2; TOUGH2-MP
.. index:: TOUGH2; TOUGH+
.. index:: TOUGH3

TOUGH2_MP is a multi-processor version of TOUGH2. TOUGH+ is a
redeveloped version of TOUGH2, with a more modular code structure
implemented in Fortran-95. TOUGH3 is another parallelized
redevelopment of TOUGH2.

TOUGH2 data files
~~~~~~~~~~~~~~~~~

.. index:: TOUGH2 data files
.. index:: TOUGH2
.. index:: TOUGH2; AUTOUGH2

TOUGH2 takes its main input from a **data file**, which contains
information about the model grid, simulation parameters, time
stepping, sources of heat and mass etc. The data file formats for
TOUGH2 and AUTOUGH2 are almost identical, with minor
differences. TOUGH2_MP can read TOUGH2 data files, but also supports
some extensions (e.g. for 8-character instead of 5-character block
names) to this format. PyTOUGH does not currently support the
TOUGH2_MP extensions. TOUGH+ and TOUGH3 data files can also have some
extensions, which PyTOUGH does not support as yet.

Because TOUGH2 uses a finite volume formulation, the only model grid
data it needs are the volumes of the grid blocks and the distances and
areas associated with the connections between blocks. Hence, the TOUGH2
data file need not contain any information about the specific locations
of the blocks in space, and it contains no information about the
locations of the vertices or edges of the blocks. This makes it easy to
use TOUGH2 to simulate one-, two- or three-dimensional models, all with
the same format of data file. However, this lack of reference to any
coordinate system also makes it more difficult to generate model grids,
and to visualise simulation results in space.

MULgraph geometry files
~~~~~~~~~~~~~~~~~~~~~~~

.. index:: MULgraph geometry; files

For this reason, a separate **geometry file** can be used to create
grids for TOUGH2 simulations and visualise simulation results. The
geometry file contains information about the locations of the grid
block vertices. The geometry file can be used to visualise results
using the `TIM <https://tim.readthedocs.io/>`_ graphical interface for
TOUGH2 and AUTOUGH, developed at the University of Auckland. (This
file format was originally designed for use with TIM's predecessor,
MULgraph).

The MULgraph geometry file assumes the grid has a layered structure,
with blocks arranged in layers and columns, and the same arrangement of
columns on each layer. (At the top of the model grid, blocks in some
columns may be missing, to allow the grid to follow the surface
topography.)

If you do not have a MULgraph geometry file for your model, it is easy
to create one for a rectangular grid. In fact, PyTOUGH is able to
:ref:`reverse-engineer<sec:t2grid:rectgeo>` a MULgraph geometry from a
TOUGH2 data file containing a rectangular grid.

A specification of the MULgraph geometry file format can be found
:ref:`here<geometry_file_format>`.

TOUGH2 listing files
~~~~~~~~~~~~~~~~~~~~

.. index:: TOUGH2 listing files

The output from TOUGH2 is written to a **listing file**, which is a text
file containing tables of results for each time step (or only selected
time steps, if preferred). At each time step there is an 'element
table', containing results for block properties (e.g. pressure,
temperature etc.). There may also be a 'connection table', containing
results for flows between blocks, and a 'generation table', containing
results (e.g. flow rates) at the generators in the model (e.g. wells).

The formats of the listing files produced by TOUGH2, AUTOUGH2,
TOUGH2_MP, TOUGH+ and TOUGH3 are all slightly different, and also vary
depending on the EOS used. However, PyTOUGH attempts to detect and
read all of these formats.

What is Python?
---------------

.. index:: Python
.. index:: Python; 3.x

Python is a general-purpose programming language. It is free and
open-source, and runs on many different computer operating systems
(Linux, Windows, Mac OS X and others). Python can be downloaded from
the Python `website <http://www.python.org>`_, which also contains
detailed reference material about the Python language. If you are
using Linux you probably already have Python, as it is included in
most Linux distributions.

PyTOUGH should run on any version of Python 2.x newer than 2.4 (though
version 2.6 or newer is recommended). PyTOUGH version 1.5 or later
also runs on Python 3.x.

.. index:: Python; tutorials

If you are unfamiliar with Python (even if you have used another
programming language before), it is highly recommended that you do one
of the many Python tutorials available online, e.g.

-  http://docs.python.org/tutorial/

-  http://wiki.python.org/moin/BeginnersGuide

Python basics
~~~~~~~~~~~~~

Objects
^^^^^^^

.. index:: Python; objects

Python is what is known as an **object-oriented** language, which means
that it is possible to create special customised data types, or
'classes', to encapsulate all the properties and behaviour of the things
(objects) we are dealing with in a program. This is a very useful way of
simplifying complex programs. (In fact, in Python, everything is treated
as an object, even simple things like integers and strings.)

For example, in a TOUGH2 model grid we have collections of grid blocks,
and we need to store the names of these blocks and their volumes and
rock types. In a non-object-oriented language, these could be stored in
three separate arrays: a string array for the names, a real (or 'float')
array for the volumes and another string array for the rock types. In an
object-oriented language like Python, we can define a new data type (or
'class') for blocks, which holds the name, volume and rock type of the
block. If we declare an object called ``blk`` of this block class, we
can access or edit its volume by referring to ``blk.volume``. In this
way, we can store our blocks in one single array of block objects. When
we add or delete blocks from our grid, we can just add or delete block
objects from the array, rather than having to keep track of three
separate arrays.

In general, an object not only has **properties** (like
``blk.volume``) but also **methods**, which are functions the object
can carry out. For example, if we wanted to rotate a MULgraph geometry
file by 30°, we could do this in PyTOUGH by declaring a MULgraph
geometry file object called ``geo``, and calling its ``rotate``
method: ``geo.rotate(30)``. The methods of an object are accessed in
the same way that its properties are accessed: by adding a dot (.)
after the object's name and then adding the name of the property or
method. Any arguments of the method (e.g. the angle in the ``rotate``
function above) are added in parentheses afterwards.

Lists, dictionaries, tuples and sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: Python; lists
.. index:: Python; dictionaries
.. index:: Python; tuples
.. index:: Python; sets

Most programming languages have simple data types built in, e.g. float,
double precision or integer numbers, strings, and arrays of these.
Python has some other data types which are very useful and are used a
lot.

The first of these is the **list**. A list can contain any ordered
collection of objects, of any type, or even of different types, and is
delimited by square brackets. So for example we can declare a list
``things = [1, 'two',3.0]`` containing an integer, a string and a float.
We can access the list's elements in much the same way as we access the
elements of an array, for example ``things[1]`` would return the value
``'two'`` (note that in Python, as in most other languages besides
Fortran, the indices of arrays and lists start at 0, not 1). Additional
elements can be added to a list at any time, without having to
re-declare the size of the list: for example, ``things.append('IV')``
would add an extra element to the end of the list, giving it the value
``[1, 'two', 3.0, 'IV']``. It is also possible to remove elements from a
list, e.g. ``things.remove(3.0)``, which would give our list the value
``[1, 'two', 'IV']``.

Another useful Python data type is the **dictionary**. Dictionaries are
mainly used to store collections of objects (again, of any type or of
different types) that are referenced by name rather than by index (as in
an array or list). A dictionary is delimited by curly brackets. So for
example we can declare a dictionary
``phone = {'Eric':8155, 'Fred':2350, 'Wilma':4667}`` and then find
Fred's phone number from ``phone['Fred']``, which would return ``2350``.
For TOUGH2 models, blocks, generators, rock types and other objects are
often referred to by name rather than index, so dictionaries are an
appropriate way to store them.

A third Python data type, similar to a list, is the **tuple**. A tuple
is essentially a list that cannot be changed, and is often used just for
grouping objects together. A tuple is delimited by parentheses. For
example, ``things = (1, 'two', 3.0)`` declares a tuple with three
elements. We can still refer to the elements of a tuple using e.g.
``a[1]``, but we cannot assign new values to these elements or add or
remove elements from the tuple once it has been declared.

Python also has a **set** data type, which represents a mathematical set
- an unordered collection of objects. One of the useful aspects of sets
is that they cannot contain duplicate items. As a result, for example,
duplicate items can be removed from a list ``x`` simply by converting it
to a set, and then back to a list: ``x = list(set(x))``.

How to run Python
~~~~~~~~~~~~~~~~~

.. index:: Python; running

Python can be run either interactively or via scripts.

.. _python_interactive:

Running Python interactively
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simplest way to run Python interactively is just by typing
``python`` (or possibly ``python3``) at the command line. (On Windows
the directory that Python was installed into may have to be added to
your ``PATH`` environment variable first.) The command line then becomes
an interactive Python environment in which you can type Python commands
at the Python command prompt ``>>>``, e.g.:

::

   bob@superbox:~$ python3
   Python 3.10.12 (main, Nov 20 2023, 15:14:05) [GCC 11.4.0] on linux
   Type "help", "copyright", "credits" or "license" for more information.
   >>> things = [1, 'two', 3.0]
   >>> print(things[1])
   two
   >>> exit()
   bob@superbox:~$

In the interactive Python environment, you can view help on the
properties and methods of any Python object by typing
``help(objectname)``, where ``objectname`` is the name of an object that
has been declared. This will list the properties and methods of the
object and a description of each one.

You can exit the interactive Python environment by typing ``exit()`` or
``Ctrl-Z`` on Windows, or ``Ctrl-D`` on Linux.

Python scripts
^^^^^^^^^^^^^^

.. index:: Python; scripts

The real power of Python, however, lies in using it to write **scripts**
to automate repetitive or complex tasks. You can just type Python
commands into a text file, save it with the file extension ``.py``, and
execute it by typing ``python filename.py``, where ``filename.py`` is
the name of the file. (Once again, on Windows the directory that Python
was installed into may have to be added to your ``PATH`` environment
variable first.)

You can also debug a Python script using the 'pdb' command-line
debugger. Typing ``python -m pdb filename.py`` will start debugging the
script *filename.py*.

It is also possible to run a Python script from within the interactive
Python environment. From the Python environment command line, typing
``execfile('filename.py')`` will execute the script ``filename.py``.

.. _pylibraries:

Python libraries
~~~~~~~~~~~~~~~~

.. index:: Python; libraries

Python comes with a large number of features already built in, but for
specialised tasks, additional **libraries** of Python software can be
imported into Python as you need them. PyTOUGH itself is a set of such
libraries, and it in turn makes use of some other third-party Python
libraries. The most important of these are as follows:

Numerical Python (“NumPy”)
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: Python; numpy

`NumPy <https://numpy.org/>`_ adds a special ``numpy.array`` class for
fast multi-dimensional arrays, which PyTOUGH makes heavy use of, and a
whole range of other features, e.g. linear algebra routines, Fourier
transforms and statistics.

Scientific Python (“SciPy”)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: Python; scipy

`SciPy <http://www.scipy.org/>`_ is a library of advanced mathematical
functions (e.g. interpolation, calculus, optimisation), needed for some
PyTOUGH functionality.

Matplotlib
^^^^^^^^^^

.. index:: Python; matplotlib

`Matplotlib <https://matplotlib.org/>`_ is a library of graphical
plotting routines, which can be used for 2-D PyTOUGH visualization
tools like layer and slice plots.

Other libraries
^^^^^^^^^^^^^^^

.. index:: Visualization Tool Kit (VTK)
.. index:: Python; meshio

Some parts of PyTOUGH use other Python libraries. You do not need to
install these libraries unless you are using the parts of PyTOUGH that
depend on them. If you try to use parts of PyTOUGH that need these
libraries, and you don't have them installed, it will tell you so.

Examples:

- **VTK**, a Python interface to the `Visualization Tool Kit <http://www.vtk.org/>`_,
  a library for 3D visualisation of data via
  VTK itself, or software such as `ParaView <https://www.paraview.org/>`_,
  `Mayavi <http://docs.enthought.com/mayavi/mayavi/>`_ etc.

- **meshio**, a `library <https://pypi.org/project/meshio/>`_ for 3D
  mesh handling – used for exporting PyTOUGH grids to other formats

Importing libraries
^^^^^^^^^^^^^^^^^^^

.. index:: Python; importing

To use any Python library, you just need to **import** it first. For
example, once you have installed Numerical Python, you can make it
available (in the interactive Python environment or in a Python script)
by typing the command ``import numpy``, or alternatively
``from numpy import *``. This imports all classes and commands from
Numerical Python and makes them available for use. (You can also import
only parts of a library rather than the whole thing, e.g.
``from numpy import linalg`` just imports the linear algebra routines
from Numerical Python.)

When you import a library, you can also change its name. For example,
PyTOUGH imports Numerical Python using the command
``import numpy as np``, which renames ``numpy`` as the abbreviated
``np``. This means it can, for example, access the Numerical Python
``numpy.array`` data type as ``np.array``. It also means you have access
to Numerical Python as ``np`` in your own scripts and in the interactive
Python environment, without having to import it yourself.

.. _installing:

Installing PyTOUGH
------------------

.. index:: PyTOUGH; installing

From version 1.6.0, the easiest way to install PyTOUGH is via the
``pip`` Python package installer:

::

   pip install PyTOUGH

or

::

   python -m pip install PyTOUGH

either of which will install the latest version of PyTOUGH, together
with its main dependency libraries (``numpy``, ``scipy`` and
``matplotlib``) if these are not already detected on your system.

You can also install a particular version of PyTOUGH, e.g. to install
version 1.6.0:

::

   pip install PyTOUGH==1.6.0

or upgrade your existing version of PyTOUGH:

::

   pip install --upgrade PyTOUGH

There are various ways of configuring the installation of packages with
``pip``, which may be suitable for your particular system – consult the
``pip`` `documentation <https://pip.pypa.io>`_ for details.

After installing, you should be able to import the PyTOUGH libraries
into the Python interactive environment or your Python scripts, from any
directory on your computer. For example, you can import the MULgraph
geometry library using ``from mulgrids import *`` (see :ref:`mulgrids`).

To uninstall PyTOUGH:

::

   pip uninstall PyTOUGH

Installing the testing branch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: PyTOUGH; testing branch

The PyTOUGH code exists in two main “branches”: the ``master`` branch,
which contains the latest stable release, and the ``testing`` branch,
which includes the most recent changes being tested for inclusion in the
next stable release.

If you need these most recent changes and can't wait for the next stable
release, it is possible to install the ``testing`` branch of PyTOUGH
using e.g.:

::

   pip install git+https://github.com/acroucher/PyTOUGH.git@testing

.. _unittests:

Testing PyTOUGH
---------------

.. index:: PyTOUGH; unit tests

PyTOUGH includes a suite of “unit tests” which can be used to verify
that it is working correctly. These are located in the ``tests/``
directory of the PyTOUGH repository, which includes a number of Python
scripts for testing individual PyTOUGH modules.

First you will the PyTOUGH repository on your machine. This is available
`here <https://github.com/acroucher/PyTOUGH>`_. Click the ``Code`` button
which gives various options for downloading the repository, via e.g. zip
file or Git clone.

The unit test modules in the ``tests/`` directory may be run
individually, the same way as any other Python script would be run. If
the tests in the script all pass, the last message printed out to the
console will read ``OK``. If not, details will be output regarding which
tests did not pass.

It is also possible to run the unit tests for all modules by running the
following command in the ``tests/`` directory:

::

   python -m unittest discover

or with the ``-v`` (verbose) flag to output more detail on which tests
are being run:

::

   python -m unittest discover -v

Licensing
---------

.. index:: PyTOUGH; license

PyTOUGH is free software, distributed under the GNU Lesser General
Public License (LGPL). More information is available
`here <http://www.gnu.org/licenses/>`_.
