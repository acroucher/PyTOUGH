:tocdepth: 3

.. _t2grids:

TOUGH2 grids
============

.. index:: TOUGH2 grids

.. _introduction-2:

Introduction
------------

The ``t2grids`` library in PyTOUGH contains classes and routines for
manipulating TOUGH2 grids. It can be imported using the command:

::

       from t2grids import *

``t2grid`` objects
------------------

The ``t2grids`` library defines a ``t2grid`` class, used for
representing TOUGH2 grids. This gives access via Python to the grid's
rock types, blocks, connections and other parameters.

Normally a TOUGH2 grid is not created directly, but is either read from
a TOUGH2 data file, or constructed from a :ref:`mulgrid <mulgrids>`
geometry object using the :ref:`fromgeo() <sec:t2grid:fromgeo>` method.

Printing a ``t2grid`` object (e.g. ``print(grid)``) displays a summary
of information about the grid: how many rock types, blocks and
connections it contains.

.. _properties-1:

Properties
~~~~~~~~~~

The main properties of a ``t2grid`` object are listed in the 
:ref:`table <tb:t2grid_properties>` below. Essentially a ``t2grid`` object
contains collections of blocks, rock types and connections, each
accessible either by name or by index. For example, block 'AB 20' in a
``t2grid`` called ``grid`` is given by ``grid.block['AB 20']``.

Connections are slightly different from blocks or rock types, in that
they are not named individually. However, they can be accessed by the
names of the blocks connected by the connection. For example, the
connection between blocks 'aa 10' and 'ab 10' in a ``t2grid`` called
``grid`` is given by ``grid.connection['aa 10','ab 10']``.

The ``rocktype_frequencies`` property gives information about how
frequently each rock type is used (i.e. how many blocks use that rock
type). It returns a list of tuples, the first element of each tuple
being the frequency of use, and the second element being a list of rock
type names with that frequency. The list is given in order of increasing
frequency.

The ``rocktype_indices`` property gives an ``np.array`` containing the
index of the rocktype for each block in the grid. This can be used to
give a plot of rock types, in conjunction with the ``mulgrid`` methods
``layer_plot`` or ``slice_plot``.

.. container::
   :name: tb:t2grid_properties

   .. table:: Properties of a ``t2grid`` object

      +---------------------------+----------------+-------------------------+
      | **Property**              | **Type**       | **Description**         |
      +===========================+================+=========================+
      | ``atmosphere_blocks``     | list           | atmosphere blocks       |
      +---------------------------+----------------+-------------------------+
      | ``blocklist``             | list           | blocks (by index)       |
      +---------------------------+----------------+-------------------------+
      | ``block``                 | dictionary     | blocks (by name)        |
      +---------------------------+----------------+-------------------------+
      | ``block_centres_defined`` | Boolean        | whether block centres   |
      |                           |                | have been calculated    |
      +---------------------------+----------------+-------------------------+
      | ``connectionlist``        | list           | connections (by index)  |
      +---------------------------+----------------+-------------------------+
      | ``connection``            | dictionary     | connections (by tuples  |
      |                           |                | of block names)         |
      +---------------------------+----------------+-------------------------+
      | ``num_atmosphere_blocks`` | integer        | number of atmosphere    |
      |                           |                | blocks                  |
      +---------------------------+----------------+-------------------------+
      | ``num_blocks``            | integer        | number of blocks        |
      +---------------------------+----------------+-------------------------+
      | ``num_connections``       | integer        | number of connections   |
      +---------------------------+----------------+-------------------------+
      | ``num_rocktypes``         | integer        | number of rock types    |
      +---------------------------+----------------+-------------------------+
      | ``num_underground_blocks``| integer        | number of               |
      |                           |                | non-atmosphere blocks   |
      +---------------------------+----------------+-------------------------+
      | ``rocktypelist``          | list           | rock types (by index)   |
      +---------------------------+----------------+-------------------------+
      | ``rocktype``              | dictionary     | rock types (by name)    |
      +---------------------------+----------------+-------------------------+
      | ``rocktype_frequencies``  | list of tuples | frequencies of rock     |
      |                           |                | types                   |
      +---------------------------+----------------+-------------------------+
      | ``rocktype_indices``      | ``np.array``   | index of rock type for  |
      |                           |                | each block              |
      +---------------------------+----------------+-------------------------+

.. _t2gridmethods:

Methods
~~~~~~~

The main methods of a ``t2grid`` object are listed in the following
:ref:`table <tb:t2grid_methods>`. Details of these methods are given below.

.. container::
   :name: tb:t2grid_methods

   .. table:: Methods of a ``t2grid`` object

      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | **Method**                                                               | **Type**                    | **Description**      |
      +==========================================================================+=============================+======================+
      | :ref:`+ <sec:t2grid:plus>`                                               | :ref:`t2grid <t2grids>`     | adds two grids       |
      |                                                                          |                             | together             |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`add_block <sec:t2grid:add_block>`                                  | –                           | adds a block to the  |
      |                                                                          |                             | grid                 |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`add_connection <sec:t2grid:add_connection>`                        | –                           | adds a connection to |
      |                                                                          |                             | the grid             |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`add_rocktype <sec:t2grid:add_rocktype>`                            | –                           | adds a rock type to  |
      |                                                                          |                             | the grid             |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`blockmap <sec:t2grid:blockmap>`                                    | dictionary                  | returns block name   |
      |                                                                          |                             | mapping from a       |
      |                                                                          |                             | geometry             |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`block_index <sec:t2grid:block_index>`                              | integer                     | returns index of a   |
      |                                                                          |                             | block with a         |
      |                                                                          |                             | specified name       |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`calculate_block_centres <sec:t2grid:calculate_block_centres>`      | –                           | calculates           |
      |                                                                          |                             | geometrical centre   |
      |                                                                          |                             | of all blocks in the |
      |                                                                          |                             | grid                 |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`check <sec:t2grid:check>`                                          | Boolean                     | checks grid for      |
      |                                                                          |                             | errors and           |
      |                                                                          |                             | optionally fixes     |
      |                                                                          |                             | them                 |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`clean_rocktypes <sec:t2grid:clean_rocktypes>`                      | –                           | deletes any unused   |
      |                                                                          |                             | rock types from the  |
      |                                                                          |                             | grid                 |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`connection_index <sec:t2grid:connection_index>`                    | integer                     | returns index of a   |
      |                                                                          |                             | connection with a    |
      |                                                                          |                             | specified pair of    |
      |                                                                          |                             | names                |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`copy_connection_directions <sec:t2grid:copy_connection_directions>`| –                           | copies connection    |
      |                                                                          |                             | permeability         |
      |                                                                          |                             | directions from      |
      |                                                                          |                             | another grid         |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`delete_block <sec:t2grid:delete_block>`                            | –                           | deletes a block from |
      |                                                                          |                             | the grid             |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`delete_connection <sec:t2grid:delete_connection>`                  | –                           | deletes a connection |
      |                                                                          |                             | from the grid        |
      |                                                                          |                             |                      |
      |                                                                          |                             |                      |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`delete_rocktype <sec:t2grid:delete_rocktype>`                      | –                           | deletes a rock type  |
      |                                                                          |                             | from the grid        |
      |                                                                          |                             |                      |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`demote_block <sec:t2grid:demote_block>`                            | –                           | shifts a block (or   |
      |                                                                          |                             | blocks) to the end   |
      |                                                                          |                             | of the blocklist     |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`embed <sec:t2grid:embed>`                                          | :ref:`t2grid <t2grids>`     | embeds a subgrid     |
      |                                                                          |                             | inside one block of  |
      |                                                                          |                             | another              |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`empty <sec:t2grid:empty>`                                          | –                           | empties contents of  |
      |                                                                          |                             | grid                 |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`flux_matrix <sec:t2grid:flux_matrix>`                              | ``scipy.sparse.lil_matrix`` | constructs a sparse  |
      |                                                                          |                             | matrix for           |
      |                                                                          |                             | calculating          |
      |                                                                          |                             | block-average flows  |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`fromgeo <sec:t2grid:fromgeo>`                                      | :ref:`t2grid <t2grids>`     | constructs a TOUGH2  |
      |                                                                          |                             | grid from a          |
      |                                                                          |                             | ``mulgrid`` object   |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`incons <sec:t2grid:incons>`                                        | :ref:`t2incon <incons>`     | constructs initial   |
      |                                                                          |                             | conditions for the   |
      |                                                                          |                             | grid                 |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`minc <sec:t2grid:MINC>`                                            | list                        | creates MINC blocks  |
      |                                                                          |                             | and connections      |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`radial <sec:t2grid:radial>`                                        | :ref:`t2grid <t2grids>`     | constructs a radial  |
      |                                                                          |                             | TOUGH2 grid          |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`rectgeo <sec:t2grid:rectgeo>`                                      | (:ref:`mulgrid <mulgrids>`, | constructs a         |
      |                                                                          | ``dict``)                   | ``mulgrid`` object   |
      |                                                                          |                             | from a rectangular   |
      |                                                                          |                             | TOUGH2 grid          |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`rename_blocks <sec:t2grid:rename_blocks>`                          | –                           | renames blocks the   |
      |                                                                          |                             | grid                 |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`rename_rocktype <sec:t2grid:rename_rocktype>`                      | –                           | renames a rock type  |
      |                                                                          |                             | in the grid          |
      |                                                                          |                             |                      |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`reorder <sec:t2grid:reorder>`                                      | –                           | reorders blocks and  |
      |                                                                          |                             | connections in the   |
      |                                                                          |                             | grid                 |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`rocktype_frequency <sec:t2grid:rocktype_frequency>`                | integer                     | frequency of use of  |
      |                                                                          |                             | a particular rock    |
      |                                                                          |                             | type                 |
      |                                                                          |                             |                      |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`sort_rocktypes <sec:t2grid:sort_rocktypes>`                        | –                           | sorts rock type list |
      |                                                                          |                             | into alphabetical    |
      |                                                                          |                             | order by name        |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+
      | :ref:`write_vtk <sec:t2grid:write_vtk>`                                  | –                           | writes grid to VTK   |
      |                                                                          |                             | file                 |
      |                                                                          |                             |                      |
      +--------------------------------------------------------------------------+-----------------------------+----------------------+

.. _sec:t2grid:plus:

``+``
^^^^^

Adds two grids ``a`` and ``b`` together (i.e. amalgamates them) to form
a new grid ``a+b``. If any rock types, blocks or connections exist in
both grids ``a`` and ``b``, the value from ``b`` is used, so there are
no duplicates. (Technically this is really an 'operator' rather than a
method.)

**Parameters:**

-  | **a, b**: :ref:`t2grid <t2grids>`
   | The two grids to be added together.

----

.. _sec:t2grid:add_block:

``add_block(block)``
^^^^^^^^^^^^^^^^^^^^

Adds a block to the grid. If another block with the same name already
exists, it is replaced.

**Parameters:**

-  | **block**: :ref:`t2block <t2blockobjects>`
   | Block to be added to the grid.

----

.. _sec:t2grid:add_connection:

``add_connection(connection)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adds a connection to the grid. If another connection with the same
column names already exists, it is replaced.

**Parameters:**

-  | **connection**: :ref:`t2connection <t2connectionobjects>`
   | Connection to be added to the grid.

----

.. _sec:t2grid:add_rocktype:

``add_rocktype(rock)``
^^^^^^^^^^^^^^^^^^^^^^

Adds a rock type to the grid. If another rock type with the same name
already exists, it is replaced.

**Parameters:**

-  | **rock**: :ref:`rocktype <rocktypeobjects>`
   | Rock type to be added to the grid.

----

.. _sec:t2grid:block_index:

``block_index(blockname)``
^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns the block index (in the ``blocklist`` list) of a specified block
name.

**Parameters:**

-  | **blockname**: string
   | Name of the block.

----

.. _sec:t2grid:blockmap:

``blockmap(geo, index = None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 grids; block mappings

Returns a mapping from the block name list of the specified geometry
object to the block names in the grid.

**Parameters:**

-  | **geo**: :ref:`mulgrid <mulgrids>`
   | Geometry object.

-  | **index**: list (or ``None``)
   | Specifies a list of integer indices defining which blocks in the
     grid to map to. If ``None``, all blocks are mapped to.

----

.. _sec:t2grid:calculate_block_centres:

``calculate_block_centres(geo)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculates geometrical centres of all blocks in the grid, based on the
specified geometry object ``geo``.

**Parameters:**

-  | **geo**: :ref:`mulgrid <mulgrids>`
   | Geometry object associated with the grid.

----

.. _sec:t2grid:check:

``check(fix=False,silent=False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 grids; checking

Checks a grid for errors and optionally fixes them. Errors checked for
are: blocks not connected to any other blocks, and blocks with isolated
rocktypes (not shared with any neighbouring blocks). Returns ``True`` if
no errors were found, and ``False`` otherwise. If ``fix`` is ``True``,
any identified problems will be fixed. If ``silent`` is ``True``, there
is no printout (only really useful if ``fix`` is ``True``).

Blocks not connected to any others are fixed by deleting them.
Isolated-rocktype blocks are fixed by assigning them the most popular
rocktype of their neighbours. Blocks with large volumes
(:math:`> 10^{20}` m\ :math:`^3`) are never considered isolated (because
they often have a special rocktype, such as an atmosphere one, that
their neighbours will never share).

**Parameters:**

-  | **fix**: Boolean
   | Whether to fix any problems identified.

-  | **silent**: Boolean
   | Whether to print out feedback or not.

----

.. _sec:t2grid:clean_rocktypes:

``clean_rocktypes()``
^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 grids; cleaning rocktypes

Deletes any rock types from the grid which are not assigned to any
block.

----

.. _sec:t2grid:connection_index:

``connection_index(blocknames)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns the connection index (in the ``connectionlist`` list) of the
connection between a specified pair of block names.

**Parameters:**

-  | **blocknames**: tuple
   | A pair of block names, each of type string.

----

.. _sec:t2grid:copy_connection_directions:

``copy_connection_directions(geo,grid)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Copies the connection permeability directions for horizontal connections
from another grid. It is assumed that both grids have the same column
structure, but may have different layer structures.

**Parameters:**

-  | **geo**: :ref:`mulgrid <mulgrids>`
   | Geometry object associated with the source grid.

-  | **grid**: :ref:`t2grid <t2grids>`
   | The source grid from which the connection permeability directions
     are to be copied.

----

.. _sec:t2grid:delete_block:

``delete_block(blockname)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Deletes a block from the grid. This also deletes any connections
involving the specified block.

**Parameters:**

-  | **blockname**: string
   | Name of the block to be deleted from the grid.

----

.. _sec:t2grid:delete_connection:

``delete_connection(connectionname)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Deletes a connection from the grid.

**Parameters:**

-  | **connectionname**: tuple (of string)
   | Pair of block names identifying the connection to be deleted from
     the grid.

----

.. _sec:t2grid:delete_rocktype:

``delete_rocktype(rocktypename)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Deletes a rock type from the grid.

**Parameters:**

-  | **rocktypename**: string
   | Name of the rock type to be deleted from the grid.

----

.. _sec:t2grid:demote_block:

``demote_block(blockname)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Shifts a block (or blocks) to the end of the blocklist. This can be
useful for making blocks inactive - by setting their volumes to zero or
negative, and then shifting them to the end of the list (to avoid all
blocks below them also being treated as inactive).

**Parameters:**

-  | **blockname**: string or list of strings
   | Name(s) of the block(s) to be shifted to the end of the blocklist.

----

.. _sec:t2grid:embed:

``embed(subgrid, connection)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 grids; embedding

Returns a grid with a subgrid embedded inside one of its blocks. The
connection specifies how the two grids are to be connected: the blocks
to be connected and the connection distances, area etc. between them.

**Parameters:**

-  | **subgrid**: :ref:`t2grid <t2grids>`
   | Subgrid to be embedded.

-  | **connection**: :ref:`t2connection <t2connectionobjects>` 
   | Connection specifying how the subgrid is to be embedded, including
     the connection distances and area. The first block should be the
     host block, the second the connecting block in the subgrid.

----

.. _sec:t2grid:empty:

``empty()``
^^^^^^^^^^^

.. index:: TOUGH2 grids; emptying

Empties the grid of all its blocks, rock types and connections.

----

.. _sec:t2grid:flux_matrix:

``flux_matrix(geo, blockmap = {})``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 grids; flux matrices

Takes the grid and a corresponding :ref:`mulgrid <mulgrids>` object,
and constructs a sparse matrix (of type ``scipy.sparse.lil_matrix``)
which can be used to convert connection flow values on the grid to
block-average fluxes (flows per unit area). Specifically, if an array of
connection flow values (one for each connection in the grid) is
multiplied by this sparse matrix, the result is a partitioned array
containing the 3-component block-average flux for each of the
(non-atmosphere) blocks.

The method for constructing the matrix is as follows. For each block, a
distribution of flux is fitted to agree as closely as possible with the
connection flow values. This distribution is either constant or linear,
depending on how many connections the block has (linear for blocks with
at least 6 connections). Fitting the connection values results in a
small linear system to solve, which may be under- or over-determined,
depending on the number of connections and the type of flux
distribution. A pseudo-inverse matrix is calculated which will find the
least-squares solution of this system. The total matrix is formed by
assembling these matrices for each of the blocks into a global matrix.

**Parameters:**

-  | **geo**: :ref:`mulgrid <mulgrids>`
   | The ``mulgrid`` geometry object.

-  | **blockmap**: dictionary
   | Dictionary mapping the block names in the geometry to the block
     naming system used in the grid.

----

.. _sec:t2grid:fromgeo:

``fromgeo(geo)``
^^^^^^^^^^^^^^^^

.. index:: TOUGH2 grids; from MULgraph geometry

Returns a grid constructed from a ``mulgrid`` geometry object. (Any
previous contents of the grid are first emptied.)

**Parameters:**

-  | **geo**: :ref:`mulgrid <mulgrids>`
   | The ``mulgrid`` geometry object.

----

.. _sec:t2grid:incons:

``incons(values=(101.3e3,20.))``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 grids; initial conditions

Returns a :ref:`t2incon <incons>` initial conditions object for the
grid, using the supplied values. Initial conditions can be specified for
only one block, in which case they will be applied to all blocks, or for
each block, in an array.

**Parameters:**

-  | **values**: ``tuple`` or ``np.array``
   | Initial conditions values, either a ``tuple`` of values for one
     block, or an ``np.array`` with each row containing a set of values
     for one block.

----

.. _sec:t2grid:MINC:

``minc(volume_fractions, spacing=50., num_fracture_planes=1, blocks=None, matrix_blockname=None, minc_rockname=None, proximity=None, atmos_volume=1.e25, incon=None, fracture_connection_distance=0.)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 grids; MINC

Creates "Multiple Interacting Continua" (MINC) blocks and connections
in the grid, for simulating fracture flow with matrix blocks attached
to each fracture block. This has capability similar to that of the
`GMINC <https://www.osti.gov/biblio/6065621>`_ program , or of the
MINC part of TOUGH2's :ref:`MESHMAKER <sec:t2data:meshmaker>` section
(except that matrix-matrix flow is not supported).

This function returns a rank-2 integer ``np.array`` with one row for
each MINC level, containing the indices of the blocks for that level.
For example, the first row is a list of all fracture block indices, the
second is a list of all MINC level 1 block indices, etc. This can be
useful for identifying all blocks in a given MINC level, for plotting or
other post-processing.

For example, if the output index array from this method is
``minc_level``, and ``T`` is an array of temperatures computed over the
entire MINC grid (e.g. extracted from the element table of a listing
file), then the temperatures in MINC level ``m`` are given by:

::

   T[minc_level[m]]

Note that plotting MINC results over a :ref:`mulgrid <mulgrids>`
geometry can be made easier (particularly for grids that have MINC
applied over only part of the domain) by using the
:ref:`minc_array() <sec:mulgrid:minc_array>` method to create
the solution vector to plot.

If the ``incon`` parameter is specified as a :ref:`t2incon <incons>`
object (from the original grid), then this method will also return a new
``t2incon`` object for the MINC grid, with values copied from the
original.

Fracture blocks retain the same block name as their original porous
medium blocks. The naming of matrix blocks can be controlled using the
``matrix_blockname`` parameter.

**Parameters:**

-  | **volume_fractions**: list (or ``np.array``)
   | List or array of volume fractions. The first entry corresponds to
     the fractures, with subsequent entries specifying the volume
     fractions for each MINC level. The length of this list or array is
     therefore equal to one plus the number of matrix blocks to be used.
     Entries for all MINC levels must be present, but they need not sum
     to 1- if they do not, they will be scaled so that the sum is 1.
     (This means, for example, that entries may be specified as
     percentage values.)

-  | **spacing**: float or list (or ``np.array``)
   | Fracture spacing parameters. If a float value is specified, this is
     applied to all sets of fracture planes (see below). If a list or
     array is specified, each entry is applied to its corresponding set
     of fracture planes.

-  | **num_fracture_planes**: integer
   | Number of sets of fracture planes (1, 2 or 3).

-  | **blocks**: list (or ``None``)
   | List of blocks or block names, specifying which blocks are to have
     MINC applied. If this parameter is ``None``, all blocks are
     processed (except inactive blocks).

-  | **matrix_blockname**: function (or ``None``)
   | Function returning the name of a MINC matrix block (string), given
     the original block name (string) and MINC level (integer > 0). If
     ``None``, a default function will be used, which simply replaces
     the first character of the original block name with the MINC level.

-  | **minc_rockname**: function (or ``None``)
   | Function returning the MINC rocktype name, given the original
     rocktype name and MINC level (:math:`\geq 0`). If ``None``, a
     default function will be used, which leaves fracture blocks with
     their original rocktype (the properties of which can subsequently
     be edited), and for matrix blocks, simply replaces the first
     character of the original rocktype name with 'X'.

-  | **proximity**: function (or ``None``)
   | Proximity function, returning the total matrix volume within a
     given distance (float) from the fracture faces. If ``None``, a
     default function will be used, corresponding to the
     ``num_fracture_planes`` parameter.

-  | **atmos_volume**: float
   | Maximum block volume for blocks to be considered part of the
     geometrical grid. Blocks with volume greater than this will be
     assumed to be boundary condition blocks and no MINC processing will
     be applied to them.

-  | **incon**: :ref:`t2incon <incons>` (or ``None``)
   | Initial conditions object for the original grid, before MINC
     processing. If not ``None``, then the method returns (as well as
     the block index array) a new ``t2incon`` object for the MINC grid,
     with values for each block copied from the original (for all MINC
     levels).

-  | **fracture_connection_distance**: float
   | Connection distance between fracture and matrix blocks. Default is
     zero, as in MESHMAKER, but in some situations a finite value (e.g.
     :math:`10^{-10}` m) can work better.

----

.. _sec:t2grid:radial:

``radial(rblocks, zblocks, convention=0, atmos_type=2, origin=[0,0], justify='r', case=None, dimension=2, blockmap={}, chars=ascii_lowercase, spaces=True)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 grids; radial

Returns a radial TOUGH2 grid with the specified radial and vertical
block sizes. Grid column and layer naming convention, atmosphere type
and origin can be specified. The optional ``justify`` and ``case``
parameters control the formatting of the character part of the block
names.

The ``dimension`` parameter sets the flow dimension for `"generalized
radial flow" <https://doi.org/10.1029/WR024i010p01796>`_, which can
represent flow in fractured rocks and modifies the block volumes and
areas. The default ``dimension`` = 2 corresponds to standard radial
flow.

**Parameters:**

-  | **rblocks**, **zblocks**: list (or ``np.array``)
   | Lists (or arrays) of block sizes in the *r* and *z* directions.

-  | **convention**: integer
   | Naming convention for grid columns and layers - same as the
     :ref:`naming convention <geometry_format_conventions>` for a
     :ref:`mulgrid <mulgrids>` object.

-  | **atmos_type**: integer
   | Type of atmosphere - also the same as the
     :ref:`atmosphere type <geometry_format_conventions>` for a
     :ref:`mulgrid <mulgrids>` object.

-  | **origin**: list (or ``np.array``)
   | Origin of the grid (of length 2 or 3). The first entry is the
     radial origin, i.e. the starting radius of the grid. The last entry
     is the vertical origin, i.e. the vertical position of the top of
     the grid. If of length 3, the middle entry is ignored.

-  | **justify**: string
   | Specify 'r' for the character part of the block names (first three
     characters) to be right-justified, 'l' for left-justified.

-  | **case**: string
   | Specify 'l' for the character part of the block names (first three
     characters) to be lower case, 'u' for upper case. Alternatively,
     use the more flexible ``chars`` parameter (see below).

-  | **dimension**: float
   | Dimension for 'generalized radial flow', which can take any
     (possibly non-integer) value between 1 and 3. Dimension 1
     corresponds to flow in a linear 'pipe', dimension 2 corresponds to
     standard radial flow in a disc-shaped reservoir and dimension 3
     corresponds to flow in a spherically symmetric reservoir.

-  | **blockmap**: dictionary
   | Dictionary mapping the block names in the geometry to the block
     naming system used in the grid.

-  | **chars**: string
   | Specify a string of characters to be used to form the character
     part of block names. For example, to use both lowercase and
     uppercase characters, set ``chars`` to
     ``ascii_lowercase + ascii_uppercase``, or to use uppercase letters
     only, specify ``ascii_uppercase``.

-  | **spaces**: Boolean
   | Specify ``False`` to disallow spaces in character part of block
     names. In this case, the first element of the ``chars`` parameter
     functions like a 'zero' and replaces spaces.

Visualization of radial :math:`r-z` model grids and results can be done
in PyTOUGH by creating a 'dummy' vertical slice rectangular geometry,
using the ``mulgrid`` :ref:`rectangular() <sec:mulgrid:rectangular>`
method, using its :math:`x` direction for radius (and
having only one block in the :math:`y` direction - which is not used).
The :ref:`slice_plot() <sec:mulgrid:slice_plot>` method can
then be used to plot results.

----

.. _sec:t2grid:rectgeo:

``rectgeo(origin_block=None, atmos_volume=1.e25, remove_inactive=False, convention=0, atmos_type=2, justify='r', chars=ascii_lowercase, spaces=True, layer_snap=0.1, block_order=None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; rectangular

Creates a :ref:`mulgrid <mulgrids>` geometry object from a
rectangular TOUGH2 grid. It also returns a dictionary defining the
mapping from the geometry block names to the grid block names. This
block mapping can be used when the block naming convention used by the
original TOUGH2 grid is not compatible with the layer/column based
:ref:`naming conventions <geometry_format_conventions>` assumed by a
``mulgrid`` geometry.

The method works within the following assumptions:

-  the grid is in fact rectangular (results will not be predictable
   otherwise)

-  block centre coordinates are present for all blocks in the grid

-  the bottom layer of blocks is complete (no missing blocks)

The method should work on rectangular TOUGH2 grids that have been
translated and/or horizontally rotated with respect to the coordinate
axes. Grids with incomplete upper layers (e.g. representing topography)
should also be OK.

**Parameters:**

-  | **origin_block**: string, :ref:`t2block <t2blockobjects>`
     or ``None``
   | The block on the bottom layer of the geometry, at the origin of the
     axes defined by permeability directions 1 and 2. If ``None``, it
     will be detected. Specify it manually if the algorithm does not
     detect it correctly.

-  | **atmos_volume**: float
   | Block volume below which blocks are considered part of the
     geometrical grid. Blocks with volume greater than or equal to this
     value will be assumed to be boundary condition blocks and will not
     be represented geometrically.

-  | **remove_inactive**: Boolean
   | Set ``True`` to remove inactive blocks from the geometry. TOUGH2
     treats all blocks with zero or negative volume, and all subsequent
     blocks in the block list, to be inactive. If this option is used,
     the inactive blocks will be used to detect the surface elevations
     of the columns in the geometry. Otherwise, inactive blocks will be
     retained in the geometry.

-  | **convention**: integer
   | :ref:`Naming convention <geometry_format_conventions>` for grid
     columns and layers in the output geometry.

-  | **atmos_type**: integer
   | :ref:`Atmosphere type <geometry_format_conventions>` for the
     output geometry.

-  | **justify**: string
   | Specify 'r' for the character part of the block names (first three
     characters) to be right-justified, 'l' for left-justified.

-  | **chars**: string
   | Specify a string of characters to be used to form the character
     part of block names. For example, to use both lowercase and
     uppercase characters, set ``chars`` to
     ``ascii_lowercase + ascii_uppercase``, or to use uppercase letters
     only, specify ``ascii_uppercase``.

-  | **spaces**: Boolean
   | Specify ``False`` to disallow spaces in character part of block
     names. In this case, the first element of the ``chars`` parameter
     functions like a 'zero' and replaces spaces.

-  | **layer_snap**: float
   | Smallest desired surface block thickness. Set to a positive value
     to eliminate surface blocks in the geometry with very small
     thicknesses (resulting from column surface elevations that are very
     close to the bottom of a layer). Default value is 0.1 m. Note that
     it is not recommended to use a value of zero, as spurious
     small-thickness surface blocks can arise from rounding errors in
     reading the data file. If this still occurs, try increasing the
     snap value until they disappear.

-  | **block_order**: string or ``None``
   | Specify ``None`` or 'layer_column' for default block ordering by
     layer and column, starting from the atmosphere. Specify 'dmplex' to
     order blocks by geometrical type (8-node hexahedrons first followed
     by 6-node wedges) as in PETSc DMPlex meshes.

----

.. _sec:t2grid:rename_blocks:

``rename_blocks(blockmap = {}, fix_blocknames = True)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Renames blocks in the grid according to the specified block mapping
dictionary. Any block whose name is a key of the block mapping
dictionary is renamed with the corresponding dictionary value. Related
properties such as connections are also renamed.

**Parameters:**

-  | **blockmap**: dictionary
   | Block mapping dictionary, mapping strings to strings.

-  | **fix_blocknames**: Boolean
   | Set ``True`` (the default) to 'fix' block names in the dictionary,
     using the :ref:`fix_blockname() <sec:mulgrid:fix_blockname>` function.

----

.. _sec:t2grid:rename_rocktype:

``rename_rocktype(rockname, newrockname)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Renames a rock type in the grid. An exception is raised if the specified
rocktype name does not exist, or if the new target rocktype name has
already been used.

**Parameters:**

-  | **rockname**: string
   | Name of the rock type to be renamed.

-  | **newrockname**: string
   | New name for the rock type.

----

.. _sec:t2grid:reorder:

``reorder(block_names, connection_names=None, geo=None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 grids; reordering

Reorders the blocks (and optionally connections) in the grid.

**Parameters:**

-  | **block_names**: list of string (or ``None``)
   | List specifying the names of the blocks, in their desired order.
     Each block name must exist in the grid, otherwise an error will be
     raised. If this parameter is ``None`` (the default), blocks are not
     reordered (unless a geometry is specified instead).

-  | **connection_names**: list of string (or ``None``)
   | List specifying the names of the connections, in their desired
     order. Each item in the list should be a tuple of block names. The
     ordering of the block names in any tuple may be reversed with
     respect to the original connection naming. However an error will be
     raised if any tuple of block names in the list does not exist in
     the grid (in either its forward or reverse form). If this parameter
     is ``None`` (the default), connections are not reordered (unless a
     geometry is specified instead).

-  | **geo**: :ref:`mulgrid <mulgrids>` geometry (or ``None``)
   | Geometry object to use for the reordering. If this is specified,
     the geometry's block and connection name lists are used (and the
     previous parameters are ignored). After reordering, the grid's
     blocks and connections will have the same ordering as if the grid
     had been created using the :ref:`fromgeo() <sec:t2grid:fromgeo>`
     method.

----

.. _sec:t2grid:rocktype_frequency:

``rocktype_frequency(rockname)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns the frequency of use of the rock type with the specified name,
i.e. how many blocks are assigned that rock type.

**Parameters:**

-  | **rockname**: string
   | Name of the specified rock type.

----

.. _sec:t2grid:sort_rocktypes:

``sort_rocktypes()``
^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 grids; sorting rocktypes

Sorts the rocktype list into alphabetical order by name.

----

.. _sec:t2grid:write_vtk:

``write_vtk(geo, filename, wells=False, blockmap = {}, surface_snap=0.1)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 grids; VTK

Writes a ``t2grid`` object to a VTK file on disk, for visualisation with
VTK, Paraview, Mayavi etc. The grid is written as an 'unstructured grid'
VTK object with data arrays defined on cells. The data arrays written,
in addition to the defaults arrays for the associated ``mulgrid``
object, are: rock type index, porosity and permeability for each block.
A separate VTK file for the wells in the grid can optionally be written.

**Parameters:**

-  | **geo**: :ref:`mulgrid <mulgrids>`
   | The ``mulgrid`` geometry object associated with the grid. This is
     required as the ``t2grid`` object does not contain any spatial
     information, e.g. locations of block vertices.

-  | **filename**: string
   | Name of the VTK file to be written. This is also required.

-  | **wells**: Boolean
   | Set to ``True`` if the wells from the ``mulgrid`` object are to be
     written to a separate VTK file.

-  | **blockmap**: dictionary
   | Dictionary mapping the block names in the geometry to the block
     naming system used in the grid.

-  | **surface_snap**: float
   | Tolerance for specifying how close column surface elevations need
     to be before being considered "equal" when constructing surface
     nodes.

----

Other objects (``rocktype``, ``t2block`` and ``t2connection``)
--------------------------------------------------------------

A ``t2grid`` object contains lists of other types of objects:
``rocktype``, ``t2block`` and ``t2connection``. These classes are
described below.

.. _rocktypeobjects:

``rocktype`` objects
~~~~~~~~~~~~~~~~~~~~

.. index:: TOUGH2 grids; rocktypes
.. index:: rocktypes

A ``rocktype`` object represents a TOUGH2 rock type. The properties of a
``rocktype`` object, and their default values, are given in the
:ref:`table <tb:rocktype_properties>` below.

.. container::
   :name: tb:rocktype_properties

   .. table:: Properties of a ``rocktype`` object

      +--------------------------+--------------+----------------+-----------------------------------------------+
      | **Property**             | **Type**     | **Description**| **Default**                                   |
      |                          |              |                |                                               |
      +==========================+==============+================+===============================================+
      | ``capillarity``          | dictionary   | capillarity    | –                                             |
      |                          |              | function       |                                               |
      +--------------------------+--------------+----------------+-----------------------------------------------+
      | ``compressibility``      | float        | compressibility|                0 m\ :sup:`2`/N                |
      |                          |              |                |                                               |
      |                          |              |                |                                               |
      |                          |              |                |                                               |
      +--------------------------+--------------+----------------+-----------------------------------------------+
      | ``conductivity``         | float        | heat           | 1.5 W/m/K                                     |
      |                          |              | conductivity   |                                               |
      |                          |              |                |                                               |
      |                          |              |                |                                               |
      +--------------------------+--------------+----------------+-----------------------------------------------+
      | ``density``              | float        | rock grain     | 2600 kg/m\ :sup:`3`                           |
      |                          |              | density        |                                               |
      |                          |              |                |                                               |
      +--------------------------+--------------+----------------+-----------------------------------------------+
      | ``dry_conductivity``     | float        | dry heat       | wet heat                                      |
      |                          |              | conductivity   | conductivity                                  |
      +--------------------------+--------------+----------------+-----------------------------------------------+
      | ``expansivity``          | float        | expansivity    | 0 K\ :sup:`-1`                                |
      |                          |              |                |                                               |
      |                          |              |                |                                               |
      +--------------------------+--------------+----------------+-----------------------------------------------+
      | ``klinkenberg``          | float        | Klinkenberg    | 0 Pa\ :sup:`-1`                               |
      |                          |              | parameter      |                                               |
      |                          |              |                |                                               |
      +--------------------------+--------------+----------------+-----------------------------------------------+
      | ``nad``                  | integer      | number of      | 0                                             |
      |                          |              | extra data     |                                               |
      |                          |              | lines          |                                               |
      +--------------------------+--------------+----------------+-----------------------------------------------+
      | ``name``                 | string       | rock type name | 'dfalt'                                       |
      +--------------------------+--------------+----------------+-----------------------------------------------+
      | ``permeability``         | ``np.array`` | permeability   | ``np.array``\ ([10\ :sup:`-15`]*3) m\ :sup:`2`|
      |                          |              |                |                                               |
      |                          |              |                |                                               |
      |                          |              |                |                                               |
      +--------------------------+--------------+----------------+-----------------------------------------------+
      | ``porosity``             | float        | porosity       | 0.1                                           |
      +--------------------------+--------------+----------------+-----------------------------------------------+
      | ``relative_permeability``| dictionary   | relative       | –                                             |
      |                          |              | permeability   |                                               |
      |                          |              | function       |                                               |
      +--------------------------+--------------+----------------+-----------------------------------------------+
      | ``specific_heat``        | float        | rock grain     | 900 J/kg/K                                    |
      |                          |              | specific heat  |                                               |
      |                          |              |                |                                               |
      |                          |              |                |                                               |
      +--------------------------+--------------+----------------+-----------------------------------------------+
      | ``tortuosity``           | float        | tortuosity     | 0                                             |
      |                          |              | factor         |                                               |
      +--------------------------+--------------+----------------+-----------------------------------------------+
      | ``xkd3``                 | float        | used by EOS7R  | 0 m\ :sup:`3`/kg                              |
      |                          |              |                |                                               |
      |                          |              |                |                                               |
      |                          |              |                |                                               |
      +--------------------------+--------------+----------------+-----------------------------------------------+
      | ``xkd4``                 | float        | used by EOS7R  | 0 m\ :sup:`3`/kg                              |
      |                          |              |                |                                               |
      |                          |              |                |                                               |
      |                          |              |                |                                               |
      +--------------------------+--------------+----------------+-----------------------------------------------+

The main familiar properties of a rock type are referred to in a natural
way, e.g. the porosity of a rock type ``r`` is given by ``r.porosity``.
The permeability property is a 3-element ``np.array``, giving the
permeability in each of the three principal axes of the grid, so e.g.
the vertical permeability of a rock type ``r`` would normally be given
by ``r.permeability[2]`` (recall that array indices in Python are
zero-based, so that the third element has index 2).

Some rock type properties are optional, and only need be specified when
the property ``nad`` is greater than zero. An example is the relative
permeability and capillarity functions that can be specified for a rock
type when ``nad`` :math:`\ge` 2. The way these functions are specified
is described in :ref:`TOUGH2 data files <datafiles>`.

**Example:**

::

   r = rocktype(name = 'ignim', permeability = [10.e-15, 10.e-15, 2.e-15], specific_heat = 850)

declares a rocktype object called ``r`` with name 'ignim', permeability
of 10 mD in the first and second directions and 2 mD in the vertical
direction, and specific heat 850 J.kg\ :math:`^{-1}`.K\ :math:`^{-1}`.

(Note that when declaring rock types, the permeability can for
convenience be specified as a list, which will be converted internally
to an ``np.array``.)


.. _t2blockobjects:

``t2block`` objects
~~~~~~~~~~~~~~~~~~~

.. index:: TOUGH2 grids; blocks
.. index:: blocks

A ``t2block`` object represents a block in a TOUGH2 grid. The properties
of a ``t2block`` object are given in the
:ref:`table <tb:t2block_properties>` below. These reflect the specifications of a
TOUGH2 block as given in a TOUGH2 data file, with the exception of the
``atmosphere``, ``centre``, ``connection_name``, ``neighbour_name`` and
``num_connections`` properties.

.. container::
   :name: tb:t2block_properties

   .. table:: Properties of a ``t2block`` object

      +---------------------+--------------+-----------------------+
      | **Property**        | **Type**     | **Description**       |
      +=====================+==============+=======================+
      | ``ahtx``            | float        | interface area for    |
      |                     |              | heat exchange (TOUGH2 |
      |                     |              | only)                 |
      +---------------------+--------------+-----------------------+
      | ``atmosphere``      | Boolean      | whether block is an   |
      |                     |              | atmosphere block or   |
      |                     |              | not                   |
      +---------------------+--------------+-----------------------+
      | ``centre``          | ``np.array`` | block centre          |
      |                     |              | (optional)            |
      +---------------------+--------------+-----------------------+
      | ``connection_name`` | set          | names of connections  |
      |                     |              | involving the block   |
      +---------------------+--------------+-----------------------+
      | ``nadd``            | integer      | increment between     |
      |                     |              | block numbers in      |
      |                     |              | sequence              |
      +---------------------+--------------+-----------------------+
      | ``name``            | string       | block name            |
      +---------------------+--------------+-----------------------+
      | ``neighbour_name``  | set          | names of neighbouring |
      |                     |              | (connected) blocks    |
      +---------------------+--------------+-----------------------+
      | ``nseq``            | integer      | number of additional  |
      |                     |              | blocks in sequence    |
      +---------------------+--------------+-----------------------+
      | ``num_connections`` | integer      | number of connections |
      |                     |              | containing the block  |
      +---------------------+--------------+-----------------------+
      | ``pmx``             | float        | permeability modifier |
      |                     |              | (TOUGH2 only)         |
      +---------------------+--------------+-----------------------+
      | ``rocktype``        | ``rocktype`` | rock type             |
      +---------------------+--------------+-----------------------+
      | ``volume``          | float        | block volume          |
      +---------------------+--------------+-----------------------+

The ``atmosphere`` property determines whether the block is to be
treated as an atmosphere block. The ``centre`` property can optionally
be used to specify the coordinates of the centre of a block. Block
centres are automatically calculated when a :ref:`t2grid <t2grids>`
object is constructed from a :ref:`mulgrid <mulgrids>` object using
the :ref:`fromgeo <sec:t2grid:fromgeo>` method). The
``connection_name`` property is a set containing the names (as tuples of
strings) of all connections involving the block.

A ``t2block`` object has no methods.

.. _t2connectionobjects:

``t2connection`` objects
~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: TOUGH2 grids; connections
.. index:: connections

A ``t2connection`` object represents a connections between two TOUGH2
blocks. The properties of a ``t2connnection`` object are given in the
:ref:`table <tb:t2connection_properties>` below. These correspond to the
properties of a connection specified in a TOUGH2 data file. Note that
the ``block`` property returns :ref:`t2block <t2blockobjects>`
objects, not just the names of the blocks in the connection. Hence, for
example, the volume of the first block in a connection object ``con`` is
given simply by ``con.block[0].volume``.

A ``t2connection`` object has no methods.

.. container::
   :name: tb:t2connection_properties

   .. table:: Properties of a ``t2connection`` object

      +--------------+----------+----------------------------------------------+       
      | **Property** | **Type** | **Description**                              | 
      +==============+==========+==============================================+
      | ``area``     | float    | connection area                              |
      +--------------+----------+----------------------------------------------+       
      | ``block``    | ``list`` | two-element list of blocks                   |  
      +--------------+----------+----------------------------------------------+       
      | ``dircos``   | float    | gravity direction cosine                     |  
      +--------------+----------+----------------------------------------------+       
      | ``direction``| integer  | permeability direction (1, 2, or 3)          |  
      +--------------+----------+----------------------------------------------+       
      | ``distance`` | ``list`` | two-element list of connection distances     |  
      +--------------+----------+----------------------------------------------+       
      | ``nad1,nad2``| integer  | increments in sequence numbering             |  
      +--------------+----------+----------------------------------------------+       
      | ``nseq``     | integer  | number of additional connections in sequence | 
      +--------------+----------+----------------------------------------------+       
      | ``sigma``    | float    | radiant emittance factor (TOUGH2 only)       | 
      +--------------+----------+----------------------------------------------+       

Example
-------

The following piece of Python script creates a rectangular 2-D slice
TOUGH2 grid with two rock types, and assigns these rock types to blocks
in the grid according to their position along the slice.

::

   from t2grids import *

   geo = mulgrid().rectangular([500]*20, [1000], [100]*20, atmos_type = 0, convention = 2)
   geo.write('2Dgrd.dat')
   grid = t2grid().fromgeo(geo)

   grid.add_rocktype(rocktype('greyw', permeability = [1.e-15]*2 + [0.1e-15]))
   grid.add_rocktype(rocktype('fill ', permeability = [15.e-15]*2 + [5.e-15]))

   for blk in grid.blocklist[1:]:
       if 200 <= blk.centre[0] <= 400: blk.rocktype = grid.rocktype['fill ']
       else: blk.rocktype = grid.rocktype['greyw']

The first line just imports the required PyTOUGH library. (It is not
necessary to import the ``mulgrids`` library explicitly, because it is
used and therefore imported by the ``t2grids`` library.)

The second block of code creates a rectangular MULgraph geometry object
with 20 columns (each 500 m wide) along the slice and 20 layers (each
100 m thick), writes this to a geometry file on disk, and creates a
TOUGH2 grid from it.

Then the two rock types are created, ``'greyw'`` and ``'fill '``. (Note
that rock types are expected by TOUGH2 to have names 5 characters long,
so it is necessary to add spaces to shorter names.)

The final part assigns the rock types to the blocks in the grid. The
loop starts from 1 instead of 0, so that the atmosphere block is
skipped. In this example, the blocks in the grid are assigned the
``'fill '`` rock type if they are between 200 m and 400 m along the
slice. Blocks outside this region are assigned the ``'greyw'`` rock
type.
