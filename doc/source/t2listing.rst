:tocdepth: 3

.. _listingfiles:

TOUGH2 listing files
====================

.. index:: TOUGH2 listing files

.. _introduction-5:

Introduction
------------

The ``t2listing`` library in PyTOUGH contains classes and routines for
reading TOUGH2 listing files. It can be imported using the command:

::

      from t2listing import *

Listing files produced by AUTOUGH2, TOUGH2, TOUGH2_MP, TOUGH+ and TOUGH3
have different formats but are all supported. The main listing files
produced by TOUGHREACT are also supported. (There is also a separate
:ref:`toughreact_tecplot <toughreact_tecplot>` class for handling
the additional Tecplot output files produced by TOUGHREACT.)

``t2listing`` objects
---------------------

The ``t2listing`` library defines a ``t2listing`` class, used for
representing TOUGH2 listing files.

**Example:**

::

   lst = t2listing()

creates an empty ``t2listing`` object called ``lst``.

::

   lst = t2listing(filename)

creates a ``t2listing`` object called ``lst`` and reads its contents
from file ``filename``.

.. _t2listing_properties:

Properties
~~~~~~~~~~

The main properties of a ``t2listing`` object are listed in the
:ref:`table <tb:t2listing_properties>` below.

Element, connection and generation tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 listing files; tables

There are three main 'table' properties, corresponding to the
**element**, **connection** and **generation** tables in the listing
file. These are all of type :ref:`listingtable <listingtableobjects>` and
provide access to the simulation results. Not all of these tables will
necessarily be present - this depends on the settings in the data file
which produced the results. For TOUGH2 results, a fourth **primary**
table may also be present, containing primary variables and their
changes, if the KDATA parameter is set to 3. TOUGH+ results can also
contain additional element tables containing other calculated
quantities; these are named **element1**, **element2** etc. A list of
names of all available tables is given by the  property.

For example, for a ``t2listing`` object ``lst``,
``lst.element['AR210']['Temperature']`` gives the temperature at block
'AR210', at the current time. Blocks can also be identified by index
rather than name, so that ``lst.element[120]['Pressure']`` gives the
pressure at the block with (zero-based) index 120.

These tables can also be accessed to give all results for a given block,
or for a given column in the table. For example,
``lst.element['AR210']`` returns a dictionary containing all results at
block 'AR210', referred to by the name of each table column.
``lst.element['Temperature']`` returns an ``np.array`` containing the
temperatures at all blocks in the model. (Hence,
``lst.element['Pressure'][120]`` gives the same result as
``lst.element[120]['Pressure']``.)

The connection and generation tables work very similarly to the element
table, except that connections are referred to by tuples of block names
(rather than single block names), and generators are referred to by
tuples of block names and generator names. So for example, the mass flow
rate between blocks 'AB300' and 'AC300' might be given by
``lst.connection['AB300', 'AC300']['Mass flow']``.

The names of the columns for each table are read directly from the
listing file, and will depend on the TOUGH2 equation of state (EOS)
being used.

Skipping tables
^^^^^^^^^^^^^^^

The default behaviour is for a ``t2listing`` object to read all tables
present in the listing file. However, it is possible to skip the reading
of specified tables if required. This can be useful for speeding up
reading of large listing files where not all tables are required. For
example, sometimes the connection data are not required, but for large
models the connection table is often much bigger than the others, so
skipping it can make reading significantly faster. Data in skipped
tables are not available either via their corresponding properties or
via the :ref:`history() <sec:t2listing:history>` method.

To skip tables, specify their table names (``element``, ``connection``
etc.) in the optional ``skip_tables`` parameter when creating the
``t2listing`` object. (By default, this parameter is an empty list.) For
example, to read a listing file with name 'output.listing' into the
object ``lst`` and skip reading the connection and generation tables:

::

   lst = t2listing('output.listing', skip_tables = ['connection', 'generation'])

File encoding
^^^^^^^^^^^^^

It is possible to specify the file encoding for the listing file using
the optional ``encoding`` parameter when creating the ``t2listing``
object. The default for this parameter is "latin-1" encoding which
should be fine for reading in most listing files. If you encounter
exotic characters in your listing files which are not read correctly
using the default encoding you may want to try other encodings.

Full and short output
^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 listing files; short output

AUTOUGH2 allows the use of 'short' output, in which a specified
selection of block, connection or generator properties are printed at
time steps between normal full output. A ``t2listing`` object will read
short output results, if they are present, when producing time histories
using the :ref:`history() <sec:t2listing:history>` method.
However it is not possible to navigate to short output results or access
them via the ``t2listing`` table properties above.

TOUGH2, TOUGH2_MP, TOUGHREACT, TOUGH+ and TOUGH3 do not support short
output.

Navigating in time using ``time``, ``index`` and ``step``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 listing files; navigation

The ``time`` property returns the time (in seconds) corresponding to the
current set of results. It is also possible to set the ``time`` property
to navigate to a specific set of full results. For example,
``lst.time=1.e9`` navigates to the set of full results with time closest
to :math:`10^9`\ s.

The ``index`` property gives the index of the current set of results,
and can take any value between 0 and ``num_fulltimes``-1. The value of
``index`` can also be set to change to a different set of results in the
listing file (e.g. ``lst.index=12``). It can be incremented and
decremented like any other Python integer variable, e.g.
``lst.index+=1`` or ``lst.index-=2`` to go to the next set of results,
or the second to last set respectively.

The ``step`` property gives the time step number for the current set of
results. This is the number of time steps carried out in the simulation
up to the current set of results (recall that results are not
necessarily written to the listing file at every time step). Again, its
value can be set to navigate through the results, e.g. ``lst.step=100``
navigates to the set of full results with time step number nearest to
100.

The ``times`` property returns an ``np.array`` of all times at which
results (including short output) are given in the listing file. It has
length equal to ``num_times``. The ``fulltimes`` property returns an
``np.array`` of times at which full results are given (not including
short output), and has length equal to ``num_fulltimes``.

A ``t2listing`` object also has :ref:`methods <t2listingmethods>` (as
well as properties) for navigating through time.

Listing diagnostics
^^^^^^^^^^^^^^^^^^^

``t2listing`` objects have two properties that provide diagnostics on
the results of the TOUGH2 run.

The ``convergence`` property is a dictionary of the maximum absolute
differences in the element table between the second to last and last
sets of results in the listing file. This can be used to check
convergence of steady-state simulations. For example:

::

   lst.convergence['Temperature']

gives the largest absolute temperature change between the second to last
and last sets of results.

The ``reductions`` property is a list of tuples of time step indices at
which the time step size was reduced during the simulation, and the
block name at which the maximum residual occurred prior to each
reduction. This gives an indication of problematic times and blocks
which caused time step reductions.

.. container::
   :name: tb:t2listing_properties

   .. table:: Properties of a ``t2listing`` object

      +-------------------+-------------------------------------------+-----------------------+
      | **Property**      | **Type**                                  | **Description**       |
      +===================+===========================================+=======================+
      | ``connection``    | :ref:`listingtable <listingtableobjects>` | connection table for  |
      |                   |                                           | current set of        |
      |                   |                                           | results               |
      |                   |                                           |                       |
      +-------------------+-------------------------------------------+-----------------------+
      | ``convergence``   | dictionary                                | maximum differences   |
      |                   |                                           | in element table      |
      |                   |                                           | between second to     |
      |                   |                                           | last and last sets of |
      |                   |                                           | results               |
      +-------------------+-------------------------------------------+-----------------------+
      | ``element``       | :ref:`listingtable <listingtableobjects>` | element table for     |
      |                   |                                           | current set of        |
      |                   |                                           | results               |
      +-------------------+-------------------------------------------+-----------------------+
      | ``element1`` etc. | :ref:`listingtable <listingtableobjects>` | additional element    |
      |                   |                                           | table for current set |
      |                   |                                           | of results (TOUGH+    |
      |                   |                                           | only)                 |
      +-------------------+-------------------------------------------+-----------------------+
      | ``filename``      | string                                    | name of listing file  |
      |                   |                                           | on disk               |
      +-------------------+-------------------------------------------+-----------------------+
      | ``fullsteps``     | ``np.array``                              | array of time step    |
      |                   |                                           | numbers (integer) for |
      |                   |                                           | full results          |
      +-------------------+-------------------------------------------+-----------------------+
      | ``fulltimes``     | ``np.array``                              | array of times        |
      |                   |                                           | (float) for full      |
      |                   |                                           | results               |
      +-------------------+-------------------------------------------+-----------------------+
      | ``generation``    | :ref:`listingtable <listingtableobjects>` | generation table for  |
      |                   |                                           | current set of        |
      |                   |                                           | results               |
      |                   |                                           |                       |
      +-------------------+-------------------------------------------+-----------------------+
      | ``index``         | integer                                   | index of current set  |
      |                   |                                           | of results            |
      +-------------------+-------------------------------------------+-----------------------+
      | ``num_fulltimes`` | integer                                   | number of sets of     |
      |                   |                                           | full results          |
      +-------------------+-------------------------------------------+-----------------------+
      | ``num_times``     | integer                                   | number of sets of all |
      |                   |                                           | results (full and     |
      |                   |                                           | short)                |
      +-------------------+-------------------------------------------+-----------------------+
      | ``primary``       | :ref:`listingtable <listingtableobjects>` | primary variable      |
      |                   |                                           | table for current set |
      |                   |                                           | of results (TOUGH2    |
      |                   |                                           | only)                 |
      +-------------------+-------------------------------------------+-----------------------+
      | ``reductions``    | list                                      | time step indices at  |
      |                   |                                           | which time step was   |
      |                   |                                           | reduced during the    |
      |                   |                                           | simulation            |
      +-------------------+-------------------------------------------+-----------------------+
      | ``short_types``   | list of string                            | types of short output |
      |                   |                                           | present               |
      +-------------------+-------------------------------------------+-----------------------+
      | ``simulator``     | string                                    | detected simulator    |
      |                   |                                           | ('AUTOUGH2', 'TOUGH2' |
      |                   |                                           | etc.)                 |
      +-------------------+-------------------------------------------+-----------------------+
      | ``step``          | integer                                   | time step number of   |
      |                   |                                           | current set of        |
      |                   |                                           | results               |
      +-------------------+-------------------------------------------+-----------------------+
      | ``steps``         | ``np.array``                              | array of time step    |
      |                   |                                           | numbers (integer) for |
      |                   |                                           | all results (full and |
      |                   |                                           | short)                |
      +-------------------+-------------------------------------------+-----------------------+
      | ``table_names``   | list                                      | names of available    |
      |                   |                                           | tables                |
      +-------------------+-------------------------------------------+-----------------------+
      | ``time``          | float                                     | time of current set   |
      |                   |                                           | of results            |
      +-------------------+-------------------------------------------+-----------------------+
      | ``times``         | ``np.array``                              | array of times        |
      |                   |                                           | (float) for all       |
      |                   |                                           | results (full and     |
      |                   |                                           | short)                |
      +-------------------+-------------------------------------------+-----------------------+
      | ``title``         | string                                    | simulation title      |
      +-------------------+-------------------------------------------+-----------------------+

.. _t2listingmethods:

Methods
~~~~~~~

The main methods of a ``t2listing`` object are listed in the
:ref:`table <tb:t2listing_methods>` below.

.. container::
   :name: tb:t2listing_methods

   .. table:: Methods of a ``t2listing`` object

      +------------------------------------------------------------+---------------+-------------------------+
      | **Method**                                                 | **Type**      | **Description**         |
      +============================================================+===============+=========================+
      | :ref:`add_side_recharge <sec:t2listing:add_side_recharge>` | –             | adds side recharge      |
      |                                                            |               | generators to a         |
      |                                                            |               | ``t2data`` object       |
      +------------------------------------------------------------+---------------+-------------------------+
      | :ref:`close <sec:t2listing:close>`                         | –             | closes listing file     |
      |                                                            |               |                         |
      |                                                            |               |                         |
      +------------------------------------------------------------+---------------+-------------------------+
      | :ref:`first<sec:t2listing:first>`                          | –             | navigates to the first  |
      |                                                            |               | set of full results     |
      |                                                            |               |                         |
      +------------------------------------------------------------+---------------+-------------------------+
      | :ref:`get_difference<sec:t2listing:get_difference>`        | dictionary    | maximum differences in  |
      |                                                            |               | element table between   |
      |                                                            |               | two sets of results     |
      +------------------------------------------------------------+---------------+-------------------------+
      | :ref:`history<sec:t2listing:history>`                      | list or tuple | time history for a      |
      |                                                            |               | selection of locations  |
      |                                                            |               | and table columns       |
      +------------------------------------------------------------+---------------+-------------------------+
      | :ref:`last<sec:t2listing:last>`                            | –             | navigates to the last   |
      |                                                            |               | set of full results     |
      |                                                            |               |                         |
      +------------------------------------------------------------+---------------+-------------------------+
      | :ref:`next<sec:t2listing:next>`                            | Boolean       | navigates to the next   |
      |                                                            |               | set of full results     |
      |                                                            |               |                         |
      +------------------------------------------------------------+---------------+-------------------------+
      | :ref:`prev<sec:t2listing:prev>`                            | Boolean       | navigates to the        |
      |                                                            |               | previous set of full    |
      |                                                            |               | results                 |
      +------------------------------------------------------------+---------------+-------------------------+
      | :ref:`write_vtk<sec:t2listing:write_vtk>`                  | –             | writes results to VTK   |
      |                                                            |               | file                    |
      |                                                            |               |                         |
      +------------------------------------------------------------+---------------+-------------------------+

Details of these methods are as follows.

----

.. _sec:t2listing:add_side_recharge:

``add_side_recharge(geo, dat)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adds side recharge generators to a ``t2data`` object ``dat`` for a
production run, calculated according to the final results in the
listing. These generators represent side inflows due to pressure changes
in the blocks on the model's horizontal boundaries. Recharge generators
are given the names of their blocks- any existing generators with the
same names will be overwritten.

**Parameters:**

-  | **geo**: :ref:`mulgrid <mulgrids>`
   | Geometry object associated with the listing.

-  | **dat**: :ref:`t2data <datafiles>`
   | TOUGH2 data object for the side recharge generators to be added to.

----

.. _sec:t2listing:close:

``close()``
^^^^^^^^^^^

Closes the listing file after use.

----

.. _sec:t2listing:first:

``first()``
^^^^^^^^^^^

Navigates to the first set of full results in the listing file.

----

.. _sec:t2listing:get_difference:

``get_difference(indexa=None, indexb=None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns dictionary of maximum differences, and locations of difference,
of all element table properties between two sets of results.

**Parameters:**

-  | **indexa**, **indexb**: integer or ``None``
   | Indices of results between which the maximum differences are to be
     calculated. If both indexa and indexb are provided, the result is
     the difference between these two result indices. If only one index
     is given, the result is the difference between the given index and
     the one before that. If neither are given, the result is the
     difference between the last and penultimate sets of results.

----

.. _sec:t2listing:history:

``history(selection, short=True, start_datetime=None``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 listing files; time histories

Returns a list of time histories (as ``np.arrays``) for specified
locations and table columns in the element, connection or generation
tables. For each selection, a tuple of two ``np.arrays`` is returned,
one each for times and values. Short output (AUTOUGH2 only) can be
omitted from the history results by setting the ``short`` parameter to
``False``. If the ``start_datetime`` parameter is given, times in the
output are given as datetimes rather than seconds from the start.

**Parameters:**

-  | **selection**: list of tuples
   | Selection of listing tables, locations (or indices) and table
     columns to produce histories for. Each tuple contains three
     elements: the listing **table type** ('e', 'c', 'p' or 'g' for
     element, connection, primary or generation table respectively), the
     **block/ connection/ generator name** (or index) and the **table
     column name**. (If only a single tuple is given instead of a list
     of tuples, just the single tuple of times and values for that
     selection is returned.) For history of additional element tables in
     TOUGH+ results, use 'e1', 'e2' etc. instead of 'e'. Note that, as
     for listing tables, connection and generator names (or 'keys') are
     specified as two-element tuples (see
     :ref:`Keys for different listing table types <tb:listing_table_keys>`).
     If the second element of a selection tuple is an integer, it will
     be interpreted as the (zero-based) index of the block, connection
     or generator in the corresponding table.

-  | **short**: Boolean
   | Whether short output (AUTOUGH2 only) is to be included in the
     history results - default is ``True``.

-  | **start_datetime**: datetime or ``None``
   | Datetime of the start of the simulation. If ``None`` (the default),
     output times are given as seconds from the start of the simulation.
     If a Python datetime is given, then output times are given as
     datetimes.

**Examples:**

::

   [(tt,temp), (tq,q), (tg,g)] = lst.history([('e', 'AR210', 'Temperature'),
   ('c', ('AB300','AC300'), 'Mass flow'), ('g', ('BR110','SO  1'), 'Generation rate')])

returns a list of three tuples of ``np.arrays``, ``(tt,temp)``,
``(tq,q)`` and ``(tg,g)``, giving the times and values of temperature at
block 'AR210', mass flow at the connection between blocks 'AB300' and
'AC300', and generation rate in the generator 'SO 1' in block 'BR110'
respectively.

::

     from datetime import datetime
     t0 = datetime(1955, 1, 1)
     t,T = lst.history(('e', 'AR210', 'Temperature'), start_datetime = t0)

returns ``T`` as an ``np.array`` of temperature values, and ``t`` as an
``np.array`` of Python datetimes, starting at 1 January 1955.

----

.. _sec:t2listing:last:

``last()``
^^^^^^^^^^

Navigates to the last set of full results in the listing file.

----

.. _sec:t2listing:next:

``next()``
^^^^^^^^^^

Navigates to the next set of full results in the listing file. Returns
``False`` if already at the last set of results (and ``True``
otherwise).

----

.. _sec:t2listing:prev:

``prev()``
^^^^^^^^^^

Navigates to the previous set of full results in the listing file.
Returns ``False`` if already at the first set of results (and ``True``
otherwise).

----

.. _sec:t2listing:write_vtk:

``write_vtk(geo, filename, grid=None, indices=None, flows=False, wells=False, start_time=0, time_unit='s', flux_matrix=None, blockmap = {}, surface_snap=0.1)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 listing files; VTK

Writes a ``t2listing`` object to a set of VTK files on disk, for
visualisation with VTK, Paraview, Mayavi etc. The results in the listing
object are written as an 'unstructured grid' VTK object with data arrays
defined on cells. The data arrays written correspond to the variables
given in the columns of the element table of the ``t2listing`` object.
(For TOUGH+ results, variables from the additional element tables are
also included.) In addition, data arrays from an associated ``mulgrid``
and (optionally) ``t2grid`` objects can be included.

If ``flows`` is ``True`` (and a ``grid`` is specified and the listing
contains connection data), approximate block-average flux vectors at the
centre of each block are also written, for all variables in the
connection table with names ending in 'flow'.

One \*.vtu file is produced for each time step in the ``t2listing``
object at which full results are present, and a \*.pvd file is also
written. This is usually the file that should actually be opened in
Paraview or other software as it contains time information associated
with each \*.vtu file.

Optionally, only a subset of the time indices present in the
``t2listing`` can be written, according to the ``indices`` parameter. A
start time and time unit for the output can optionally be specified.

**Parameters:**

-  | **geo**: :ref:`mulgrid <mulgrids>`
   | The ``mulgrid`` geometry object associated with the results. For
     flexibility, this geometry need not be fully compatible with the
     results – for example, it may contain only a subset of the blocks
     for which results are present, or the blocks may be in a different
     order. However, if it is not fully compatible, the writing process
     will be slower, and flux vectors will not be written (even if
     ``flows`` is set to ``True``).

-  | **filename**: string
   | Name of the \*.pvd file to be written. Names of the individual
     \*.vtu files for each time step are similar but with a time index
     appended and the file extension changed.

-  | **grid**: :ref:`t2grid <t2grids>`
   | Name of optional ``t2grid`` object associated with the results.

-  | **indices**: list or tuple
   | Optional specification of time indices to include in the output. If
     set to ``None`` (the default), all time indices will be included.

-  | **flows**: Boolean
   | Set to ``True`` if approximate block-centred flux vectors are to be
     calculated and written, for visualising flows. Default is
     ``False``. **Note**: flow vectors can only be calculated if a
     **grid** is specified.

-  | **wells**: Boolean
   | Set to ``True`` if a separate VTK file for the wells in the
     :ref:`mulgrid <mulgrids>` object is to be written. Default is
     ``False``.

-  | **start_time**: float
   | Optional start time of the simulation, i.e. time associated with
     the first set of results. Default is zero.

-  | **time_unit**: string
   | Optional time unit for the output. TOUGH2 results are given at
     times in seconds, but this option allows them to be converted to
     other units. Options are: 's', 'h', 'd' and 'y', for seconds,
     hours, days and years respectively. Default is 's'.

-  | **flux_matrix**: ``scipy.sparse.lil_matrix``
   | Sparse matrix that multiplies a vector of connection values to
     produce a partition vector of 3-D block average flows at the
     (underground) block centres. One of these can be produced using the
     ``t2grid.flux_matrix()`` method, and a corresponding ``mulgrid``
     object. A flux matrix will be calculated internally if not
     supplied.

-  | **blockmap**: dictionary
   | Dictionary mapping the block names in the geometry to the block
     naming system used in the listing.

-  | **surface_snap**: float
   | Tolerance for specifying how close column surface elevations need
     to be before being considered "equal" when constructing surface
     nodes.

----

.. _listingtableobjects:

``listingtable`` objects
------------------------

.. index:: TOUGH2 listing files; tables

A ``listingtable`` object represents a table of results in a TOUGH2
listing file (whether it is an element, connection or generation table).
The column headings of the table are taken directly from the
corresponding table in the listing file. The rows of the table may be
accessed either by (zero-based) index, or by the 'key' for the table
row, which depends on the table type (see :ref:`table <tb:listing_table_keys>`
below).

.. container::
   :name: tb:listing_table_keys

   .. table:: Keys for different listing table types

      ============== ============================
      **Table type** **Key**
      ============== ============================
      ``element``    block name
      ``connection`` (block name 1, block name 2)
      ``generation`` (block name, generator name)
      ============== ============================

Hence, the value in the element table for a given block and column can
be accessed by ``lst.element[blockname][columnname]``, or by
``lst.element[blockindex][columnname]`` (for a ``t2listing`` object
``lst``). Note that for connection and generation tables, the keys are
tuples of two strings. For connection tables, the order of these two
strings (the block names) is not important; if the listing file contains
results for (block1, block2), then results for (block2, block1) can be
accessed via the corresponding ``listingtable`` object (though the
results will have the opposite sign to those in the file, as they will
represent flows in the opposite direction).

The values for an entire row or column of the table can also be
accessed, for example ``lst.element[blockname]`` gives the row in the
table for a specified block, with the values arranged in a dictionary
which can be accessed using the column names of the table (e.g.
``lst.element['AR231']['Temperature']``). This dictionary for each row
also contains an additional ``'key'`` item which returns the key for
that row. Conversely, ``lst.element[columnname]`` gives the column in
the table for a specified column name, with the values returned in an
``np.array`` (one value for each block in the grid, for an element
table).

``listingtable`` properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The properties of a ``listingtable`` object are given in the
:ref:`table <tb:listingtable_properties>` below. The entire list of
key values for a ``listingtable`` may be accessed via the ``row_name``
property, which contains the key value for each row. The column
headings of the table can similarly be accessed via the
``column_name`` list property. The ``num_rows`` and ``num_columns``
properties of a ``listingtable`` object return the numbers of rows and
columns respectively. The ``num_keys`` property just returns the
number of keys used to identify each row - generally 1 for an element
table and 2 for connection and generation tables.

.. container::
   :name: tb:listingtable_properties

   .. table:: Properties of a ``listingtable`` object

      +-----------------+------------------+--------------------------+
      | **Property**    | **Type**         |  **Description**         |
      +=================+==================+==========================+
      | ``column_name`` | list             | column headings          |
      +-----------------+------------------+--------------------------+
      | ``DataFrame``   | Pandas DataFrame | data in DataFrame format |
      +-----------------+------------------+--------------------------+
      | ``num_columns`` | integer          | number of columns        |
      +-----------------+------------------+--------------------------+
      | ``num_keys``    | integer          | number of keys per row   |
      +-----------------+------------------+--------------------------+
      | ``num_rows``    | integer          | number of rows           |
      +-----------------+------------------+--------------------------+
      | ``row_name``    | list             | keys for each row        |
      +-----------------+------------------+--------------------------+

Adding and subtracting
~~~~~~~~~~~~~~~~~~~~~~

It is possible to perform addition and subtraction operations on
``listingtable`` objects. Subtraction can be useful, for example, when
comparing results from different runs. These operations can only be
carried out when the row and column names of the two tables are
identical. The resulting table will have the same row and column names
as the original tables, but will contain the element-wise sums or
differences.

Converting to DataFrames
~~~~~~~~~~~~~~~~~~~~~~~~

A ``listingtable`` object has a ``DataFrame`` property which returns
the entire table in the form of a `Pandas
<http://pandas.pydata.org/>`_ DataFrame object. Pandas is a Python
library for data analysis, which you will need to have installed
before you can use the ``DataFrame`` property. With Pandas you can do
advanced data analysis on your TOUGH2 results. See the Pandas
documentation for more details.

``listingtable`` methods
~~~~~~~~~~~~~~~~~~~~~~~~

``listingtable`` objects have one method as described below.

----

``rows_matching(pattern, index=0, match_any=False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns a list of rows in the table with keys matching the specified
regular expression string, ``pattern``.

For tables with multiple keys, ``pattern`` can be a list or tuple of
regular expressions. If a single string pattern is given for a
multiple-key table, the pattern is matched on the index\ :math:`^{th}`
key (and any value of the other key - unless the ``match_any`` option is
used; see below).

If ``match_any`` is set to ``True``, rows are returned with keys
matching any of the specified patterns (instead of all of them). If this
option is used in conjunction with a single string pattern, the
specified pattern is applied to all keys.

**Parameters:**

-  | **pattern**: string, list or tuple
   | Regular expression string specifying the pattern to match. For
     multiple-key tables, this can be a list or tuple of regular
     expression strings.

-  | **index**: integer
   | Index of the key to which the pattern is to be applied, for
     multiple-key tables and when ``pattern`` is specified as a single
     string.

-  | **match_any**: Boolean
   | If ``False``, return only rows with keys matching *all* of their
     corresponding patterns. If ``True``, return rows with keys matching
     *any* of the specified patterns - and if a single string pattern is
     given, apply this to all keys.

----

``t2historyfile`` objects
-------------------------

.. index:: TOUGH2 data files; FOFT
.. index:: TOUGH2 data files; COFT
.. index:: TOUGH2 data files; GOFT

In addition to the main listing file, TOUGH2 can optionally produce
extra files containing time history data from selected blocks,
connections or generators, named ``FOFT``, ``COFT`` and ``GOFT`` files
respectively. TOUGH+ can optionally name these files
``Elem_Time_Series``, ``Conx_Time_Series`` and ``SS_Time_Series``
instead. (AUTOUGH2 does not produce separate history files, but can
instead produce 'short output' at selected blocks, connections or
generators within the listing file itself.)

The ``t2listing`` module contains a ``t2historyfile`` class for reading
and manipulating these history files. History files produced by TOUGH2,
TOUGH2_MP and TOUGH+ are supported, although they all have different
formats. The same class is used for FOFT, COFT and GOFT files. A history
file of any of these types can be opened using a command such as:

::

   hist = t2historyfile(filename)

where ``filename`` is the name of the file. It may contain wildcards (*)
so that several files matching a pattern are read in to the same object.
This is useful for reading output from TOUGH2_MP, which creates separate
history files for each processor used in the calculation (e.g.
``FOFT_P.000``, ``FOFT_P.001``, etc.). It is assumed that all files
opened are however of the same type (FOFT, COFT or GOFT).

Once a history file has been read in, history results for a particular
key (i.e. block, connection or generator) can be extracted. For
TOUGH2_MP, the keys are the block names for FOFT files, tuples of block
names for COFT files, and tuples of block names and source names for
GOFT files. For example:

::

   foft = t2historyfile('FOFT_P.*')
   blockname = 'fmq20'
   results = foft[blockname]

This will return a dictionary containing an ``np.array`` for each column
in the file, indexed by the column name. For example the temperature
history at this block would be given by:

::

   temp = foft[blockname]['TEMPERATURE']

Results at a particular time can also be found:

::

   time = 3.156e7
   result = foft[blockname, time]

Again, this will return a dictionary with one item for each column, but
in this case each item is just a single floating point number instead of
an array.

For **TOUGH2** rather than TOUGH2_MP, the keys are integer indices of
blocks, connections or generators, rather than names or tuples of names.
Similarly, the column names are just integers. This is because the key
names and column names are not given in TOUGH2 history files. Aside from
these differences, they can be used in the same way as TOUGH2_MP history
files, for example:

::

   foft = t2historyfile('FOFT')
   blkindex = 123
   temp = foft[blkindex][1]

For **TOUGH+** connection and generator history files (``COFT`` and
``GOFT``, or ``Conx_Time_Series`` and ``SS_Time_Series``), multiple
connections and generators can be specified as usual in the TOUGH2 input
data file, but individual results for them are not written to the
history file. Instead, the results for them are summed. As a result,
there are no 'keys' as such for accessing individual results, and the
``t2historyfile`` works a little differently. An array containing the
data in each column can be accessed by specifying the column name, for
example:

::

   ct = t2historyfile('Conx_Time_Series')
   qh = ct['HeatFlow']

The properties of a ``t2historyfile`` object are given in the
:ref:`table <tb:t2historyfile_properties>` below.

.. container::
   :name: tb:t2historyfile_properties

   .. table:: Properties of a ``t2historyfile`` object

      +-----------------+--------------+----------------------------------------------+
      | **Property**    | **Type**     | **Description**                              |
      +=================+==============+==============================================+
      | ``column_name`` | list         | column headings                              |
      +-----------------+--------------+----------------------------------------------+
      | ``key_name``    | list         | names of keys                                |
      +-----------------+--------------+----------------------------------------------+
      | ``num_times``   | integer      | number of times                              |
      +-----------------+--------------+----------------------------------------------+
      | ``num_columns`` | integer      | number of data columns                       |
      +-----------------+--------------+----------------------------------------------+
      | ``num_rows``    | integer      | total number of data (for all keys)          |
      +-----------------+--------------+----------------------------------------------+
      | ``simulator``   | string       | detected simulator ('TOUGH2' or 'TOUGH2_MP') |
      +-----------------+--------------+----------------------------------------------+
      | ``times``       | ``np.array`` | times at which results are given             |
      +-----------------+--------------+----------------------------------------------+
      | ``type``        | string       | history type ('FOFT', 'COFT' or 'GOFT')      |
      +-----------------+--------------+----------------------------------------------+

.. _toughreact_tecplot:

``toughreact_tecplot`` objects
------------------------------

.. index:: TOUGH2; TOUGHREACT

The ``t2listing`` library also defines a ``toughreact_tecplot`` class,
used for representing the additional Tecplot output files produced by
TOUGHREACT.

**Example:**

::

   tp = toughreact_tecplot(filename, blocks)

creates a ``toughreact_tecplot`` object called ``tp`` and reads its
contents from file ``filename``. The ``blocks`` object passed in as a
second parameter specifies the block names (see
:ref:`Specifying block names <toughreact_tecplot_blocknames>`).

Differences from ``t2listing`` objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A ``toughreact_tecplot`` object is similar to a
:ref:`t2listing <listingfiles>` object in many respects. Apart from the need to specify
the block names on creation (see
:ref:`Specifying block names <toughreact_tecplot_blocknames>`), the other main difference
is that unlike a ``t2listing`` object, which usually contains several
``listingtable`` objects, a ``toughreact_tecplot`` object contains only
one: the ``element`` table. Because of this, when using the ``history``
method, tables need not be specified.

These Tecplot files do not contain any information about time step
numbers, so ``t2listing`` properties like ``step`` and ``steps`` are not
present in a ``toughreact_tecplot`` object. There is also no ``title``
property, as this is not present in the Tecplot file.

There is also no 'short' output in a Tecplot file, so a
``toughreact_tecplot`` object does not have properties like
``fulltimes``, as this would just be the same as the ``times`` property.
There are also no diagnostic methods like ``convergence`` or
``reductions``.

.. _toughreact_tecplot_blocknames:

Specifying block names
~~~~~~~~~~~~~~~~~~~~~~

In the Tecplot file, results are not associated with block names, though
they appear in the same order as in the TOUGH2 data file used to
generate the results. To make results accessible by block name, a second
parameter containing the block names must be specified when a
``toughreact_tecplot`` object is created. This parameter is not
optional. It can be either:

-  a list of strings specifying the block names

-  a :ref:`mulgrid <mulgrids>` geometry object

-  a :ref:`t2grid <t2grids>` object

.. _properties-4:

Properties
~~~~~~~~~~

The main properties of a ``toughreact_tecplot`` object are listed in
the :ref:`table <tb:toughreact_tecplot_properties>` below. For more details, see
the corresponding properties of the :ref:`t2listing <t2listing_properties>` class.

.. container::
   :name: tb:toughreact_tecplot_properties

   .. table:: Properties of a ``toughreact_tecplot`` object

      +---------------+--------------------------------------------+-------------------------+
      | **Property**  | **Type**                                   | **Description**         |
      +===============+============================================+=========================+
      | ``element``   | :ref:`listingtable <listingtableobjects>`  | element table for       |
      |               |                                            | current set of results  |
      |               |                                            |                         |
      +---------------+--------------------------------------------+-------------------------+
      | ``filename``  | string                                     | name of listing file on |
      |               |                                            | disk                    |
      +---------------+--------------------------------------------+-------------------------+
      | ``index``     | integer                                    | index of current set of |
      |               |                                            | results                 |
      +---------------+--------------------------------------------+-------------------------+
      | ``num_times`` | integer                                    | number of sets of       |
      |               |                                            | results                 |
      +---------------+--------------------------------------------+-------------------------+
      | ``time``      | float                                      | time of current set of  |
      |               |                                            | results                 |
      +---------------+--------------------------------------------+-------------------------+
      | ``times``     | ``np.array``                               | array of times (float)  |
      |               |                                            | for all results         |
      +---------------+--------------------------------------------+-------------------------+

.. _methods-2:

Methods
~~~~~~~

The methods of a ``toughreact_tecplot`` object are listed in the
:ref:`table <tb:toughreact_tecplot_methods>` below.

.. container::
   :name: tb:toughreact_tecplot_methods

   .. table:: Methods of a ``toughreact_tecplot`` object

      +------------------------------------------------------+---------------+-------------------------+
      | **Method**                                           | **Type**      | **Description**         |
      +======================================================+===============+=========================+
      | :ref:`close <sec:toughreact_tecplot:close>`          | –             | closes file             |
      |                                                      |               |                         |
      |                                                      |               |                         |
      +------------------------------------------------------+---------------+-------------------------+
      | :ref:`first <sec:toughreact_tecplot:first>`          | –             | navigates to the first  |
      |                                                      |               | set of full results     |
      |                                                      |               |                         |
      +------------------------------------------------------+---------------+-------------------------+
      | :ref:`history <sec:toughreact_tecplot:history>`      | list or tuple | time history for a      |
      |                                                      |               | selection of locations  |
      |                                                      |               | and table columns       |
      +------------------------------------------------------+---------------+-------------------------+
      | :ref:`last <sec:toughreact_tecplot:last>`            | –             | navigates to the last   |
      |                                                      |               | set of full results     |
      |                                                      |               |                         |
      +------------------------------------------------------+---------------+-------------------------+
      | :ref:`next <sec:toughreact_tecplot:next>`            | Boolean       | navigates to the next   |
      |                                                      |               | set of full results     |
      |                                                      |               |                         |
      +------------------------------------------------------+---------------+-------------------------+
      | :ref:`prev<sec:toughreact_tecplot:prev>`             | Boolean       | navigates to the        |
      |                                                      |               | previous set of full    |
      |                                                      |               | results                 |
      +------------------------------------------------------+---------------+-------------------------+
      | :ref:`write_vtk <sec:toughreact_tecplot:write_vtk>`  | –             | writes results to VTK   |
      |                                                      |               | file                    |
      |                                                      |               |                         |
      |                                                      |               |                         |
      +------------------------------------------------------+---------------+-------------------------+

Details of these methods are as follows.

----

.. _sec:toughreact_tecplot:close:

``close()``
^^^^^^^^^^^

Closes the file after use.

----

.. _sec:toughreact_tecplot:first:

``first()``
^^^^^^^^^^^

Navigates to the first set of results in the Tecplot file.

----

.. _sec:toughreact_tecplot:history:

``history(selection)``
^^^^^^^^^^^^^^^^^^^^^^

Returns a list of time histories (as ``np.arrays``) for specified
locations and table columns in the element table. For each selection, a
tuple of two ``np.arrays`` is returned, one each for times and values.

**Parameters:**

-  | **selection**: list of tuples
   | Selection of locations (or indices) and table columns to produce
     histories for. Each tuple contains two elements: **block name** and
     **table column name**. (If only a single tuple is given instead of
     a list of tuples, just the single tuple of times and values for
     that selection is returned.)

----

.. _sec:toughreact_tecplot:last:

``last()``
^^^^^^^^^^

Navigates to the last set of results in the Tecplot file.

----

.. _sec:toughreact_tecplot:next:

``next()``
^^^^^^^^^^

Navigates to the next set of results in the Tecplot file. Returns
``False`` if already at the last set of results (and ``True``
otherwise).

----

.. _sec:toughreact_tecplot:prev:

``prev()``
^^^^^^^^^^

Navigates to the previous set of results in the Tecplot file. Returns
``False`` if already at the first set of results (and ``True``
otherwise).

----

.. _sec:toughreact_tecplot:write_vtk:

``write_vtk(geo, filename, grid=None, indices=None, start_time=0, time_unit='s', blockmap = {}, surface_snap=0.1)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Writes a ``toughreact_tecplot`` object to a set of VTK files on disk,
for visualisation with VTK, Paraview, Mayavi etc. The results in the
element table of the Tecplot file object are written as an 'unstructured
grid' VTK object with data arrays defined on cells. The data arrays
written correspond to the variables given in the columns of the element
table of the ``toughreact_tecplot`` object. In addition, data arrays
from an associated ``mulgrid`` and (optionally) ``t2grid`` objects can
be included.

One \*.vtu file is produced for each time step in the
``toughreact_tecplot`` object, and a \*.pvd file is also written. This
is usually the file that should actually be opened in Paraview or other
software as it contains time information associated with each \*.vtu
file.

Optionally, only a subset of the time indices present in the
``toughreact_tecplot`` can be written, according to the ``indices``
parameter. A start time and time unit for the output can optionally be
specified.

**Parameters:**

-  | **geo**: :ref:`mulgrid <mulgrids>`
   | The ``mulgrid`` geometry object associated with the results. For
     flexibility, this geometry need not be fully compatible with the
     results – for example, it may contain only a subset of the blocks
     for which results are present, or the blocks may be in a different
     order. However, if it is not fully compatible, the writing process
     will be slower.

-  | **filename**: string
   | Name of the \*.pvd file to be written. Names of the individual
     \*.vtu files for each time step are similar but with a time index
     appended and the file extension changed.

-  | **grid**: :ref:`t2grid <t2grids>`
   | Name of optional ``t2grid`` object associated with the results.

-  | **indices**: list or tuple
   | Optional specification of time indices to include in the output. If
     set to ``None`` (the default), all time indices will be included.

-  | **start_time**: float
   | Optional start time of the simulation, i.e. time associated with
     the first set of results. Default is zero.

-  | **time_unit**: string
   | Optional time unit for the output. TOUGHREACT Tecplot results are
     given at times in years, but this option allows them to be
     converted to other units. Options are: 's', 'h', 'd' and 'y', for
     seconds, hours, days and years respectively. Default is 's'.

-  | **blockmap**: dictionary
   | Dictionary mapping the block names in the geometry to the block
     naming system used in the Tecplot output.

-  | **surface_snap**: float
   | Tolerance for specifying how close column surface elevations need
     to be before being considered "equal" when constructing surface
     nodes.

----

Examples
--------

Slice plot of drawdown
~~~~~~~~~~~~~~~~~~~~~~

This script shows a vertical slice plot along the model's *x*-axis of
the difference in pressure (i.e. drawdown) between the start and end of
a simulation.

::

   from mulgrids import *
   from t2listing import *
   from copy import copy

   geo = mulgrid('gmodel.dat')
   results = t2listing('model.listing')

   results.first()
   p0 = copy(results.element['Pressure'])
   results.last()
   p1 = results.element['Pressure']

   geo.slice_plot('x', (p1-p0)/1.e5, 'Pressure\ difference', 'bar')

(Note: the ``copy`` command is needed, otherwise the arrays ``p0`` and
``p1`` would both contain the final values of pressure after the
``results.last()`` command.)

Pressure-temperature diagram
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This script plots model results from a specified block on a
pressure-temperature diagram.

::

   from t2listing import *
   import matplotlib.pyplot as plt

   lst = t2listing('model.listing')
   blk = ' n 60'
   [(tp,p), (tt,t)] = lst.history([('e', blk, 'Pressure'), ('e', blk, 'Temperature')])

   plt.plot(t, p/1.e5, 'o-')
   plt.xlabel('T ($\degree$C)')
   plt.ylabel('P (bar)')
   plt.show()

.. _comparison_example:

Comparing results of two models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This script reads grids and results for two different models, a coarse
model and a fine one, and produces a comparison plot of the time history
of temperature for both models at a given point.

::

   from mulgrids import *
   from t2listing import *
   import matplotlib.pyplot as plt

   geoc, geof = mulgrid('gcoarse.dat'), mulgrid('gfine.dat')
   coarse, fine = t2listing('coarse.listing'), t2listing('fine.listing')

   p = [47.e3, 0.0, -7000.0]
   blkc = geoc.block_name_containing_point(p)
   blkf = geof.block_name_containing_point(p)

   tc, tempc = coarse.history(('e', blkc, 'Temperature'))
   tf, temp = fine.history(('e', blkf, 'Temperature'))

   plt.plot(tc, tempc, 'o-', label = 'coarse model')
   plt.plot(tf, tempf, 's-', label = 'fine model')
   plt.xlabel('time (s)')
   plt.ylabel('Temperature ($\degree$C)')
   plt.legend()

   plt.show()

