:tocdepth: 3

.. _mulgrids:

MULgraph geometry files
=======================

.. _introduction-1:

.. index:: MULgraph geometry

Introduction
------------

The ``mulgrids`` library in PyTOUGH contains classes and routines for
creating, editing and saving MULgraph geometry files. It can be imported
using the command:

::

     from mulgrids import *

``mulgrid`` objects
-------------------

The ``mulgrids`` library defines a ``mulgrid`` class, used for
representing MULgraph geometry files.

**Example:**

::

   geo = mulgrid()

creates an empty ``mulgrid`` object called ``geo``.

::

   geo = mulgrid('geom.dat')

creates a ``mulgrid`` object called ``geo`` and reads its contents from
a file named ``'geom.dat'``.

Printing a ``mulgrid`` object (e.g. ``print(geo)``) displays a summary
of information about the grid: how many nodes, columns, layers, blocks
and wells it contains, as well as its naming convention and atmosphere
type.

A specification of the MULgraph geometry file format can be found
:ref:`here<geometry_file_format>`.

Properties
~~~~~~~~~~

The main properties of a ``mulgrid`` object are listed in the
:ref:`table <tb:mulgrid_properties>` below. Some of these properties
are 'header' information, corresponding to the data at the start of a
MULgraph geometry file (``type``, ``convention``, ``atmosphere_type``,
``atmosphere_volume``, ``atmosphere_connection`` and ``unit_type``).

The most important properties of a ``mulgrid`` object are ``node``,
``column``, ``connection``, ``layer`` and ``well``, which are
dictionaries of the grid nodes, columns, connections, layers and wells,
accessed by name. For example, grid layer 'AA' of a ``mulgrid`` object
``geo`` can be accessed by ``geo.layer['AA']``. (The ``nodelist``,
``columnlist``, ``connectionlist``, ``layerlist`` and ``welllist``
properties offer access to the nodes, columns, connections, layers and
wells by index, which is sometimes useful e.g. for looping over all
columns in the grid.)

Connections are slightly different from nodes, columns etc. in that they
are not named individually. However, they can be accessed by the names
of the columns connected by the connection. For example, the connection
between columns '10' and '11' in a ``mulgrid`` called ``geo`` is given
by ``geo.connection[' 10',' 11']``.

The elements of these lists and dictionaries are of type ``node``,
``column``, ``connection``, ``layer`` and ``well`` respectively. These
are additional object classes to represent nodes, columns, connections,
layers and wells, defined in the ``mulgrids`` library (see
:ref:`Other objects <other_mulgrid_objects>`).

.. container::
   :name: tb:mulgrid_properties

   .. table:: Properties of a ``mulgrid`` object

      +--------------------------------+--------------------+-----------------------+
      | **Property**                   | **Type**           | **Description**       |
      +================================+====================+=======================+
      | ``area``                       | float              | total horizontal area |
      |                                |                    | covered by the grid   |
      +--------------------------------+--------------------+-----------------------+
      | ``atmosphere_connection``      | float              | connection distance   |
      |                                |                    | to atmosphere blocks  |
      +--------------------------------+--------------------+-----------------------+
      | ``atmosphere_type``            | integer            | type of atmosphere    |
      +--------------------------------+--------------------+-----------------------+
      | ``atmosphere_volume``          | float              | volume of atmosphere  |
      |                                |                    | blocks                |
      +--------------------------------+--------------------+-----------------------+
      | ``bad_columns``                | set                | columns that do not   |
      |                                |                    | contain their own     |
      |                                |                    | centres               |
      +--------------------------------+--------------------+-----------------------+
      | ``bad_layers``                 | set                | layers that do not    |
      |                                |                    | contain their own     |
      |                                |                    | centres               |
      +--------------------------------+--------------------+-----------------------+
      | ``block_connection_name_index``| dictionary         | indices of block      |
      |                                |                    | connections (by name) |
      +--------------------------------+--------------------+-----------------------+
      | ``block_connection_name_list`` | list               | names of block        |
      |                                |                    | connections (by       |
      |                                |                    | index)                |
      +--------------------------------+--------------------+-----------------------+
      | ``block_name_index``           | dictionary         | indices of blocks (by |
      |                                |                    | name)                 |
      +--------------------------------+--------------------+-----------------------+
      | ``block_name_list``            | list               | names of blocks (by   |
      |                                |                    | index)                |
      +--------------------------------+--------------------+-----------------------+
      | ``block_order``                | string             | block ordering scheme |
      +--------------------------------+--------------------+-----------------------+
      | ``boundary_columns``           | set                | set of columns on the |
      |                                |                    | outer boundary of the |
      |                                |                    | grid                  |
      +--------------------------------+--------------------+-----------------------+
      | ``boundary_nodes``             | list               | ordered list of nodes |
      |                                |                    | on the outer boundary |
      |                                |                    | of the grid           |
      +--------------------------------+--------------------+-----------------------+
      | ``boundary_polygon``           | list               | list of points        |
      |                                |                    | representing grid     |
      |                                |                    | boundary (extra       |
      |                                |                    | colinear points       |
      |                                |                    | removed)              |
      +--------------------------------+--------------------+-----------------------+
      | ``bounds``                     | list               | [bottom left, top     |
      |                                |                    | right] horizontal     |
      |                                |                    | bounds of grid        |
      +--------------------------------+--------------------+-----------------------+
      | ``centre``                     | ``np.array``       | position of           |
      |                                |                    | horizontal centre of  |
      |                                |                    | the grid              |
      +--------------------------------+--------------------+-----------------------+
      | ``columnlist``                 | list               | columns (by index,    |
      |                                |                    | e.g.                  |
      |                                |                    | ``columnlist[23]``)   |
      +--------------------------------+--------------------+-----------------------+
      | ``column_angle_ratio``         | ``np.array``       | angle ratio for each  |
      |                                |                    | column                |
      +--------------------------------+--------------------+-----------------------+
      | ``column_side_ratio``          | ``np.array``       | side ratio for each   |
      |                                |                    | column                |
      +--------------------------------+--------------------+-----------------------+
      | ``column``                     | dictionary         | columns (by name,     |
      |                                |                    | e.g.                  |
      |                                |                    | ``column['AA']``)     |
      +--------------------------------+--------------------+-----------------------+
      | ``connectionlist``             | list               | connections between   |
      |                                |                    | columns (by index)    |
      +--------------------------------+--------------------+-----------------------+
      | ``connection_angle_cosine``    | ``np.array``       | angle cosines for all |
      |                                |                    | connections           |
      +--------------------------------+--------------------+-----------------------+
      | ``convention``                 | integer            | naming convention for |
      |                                |                    | columns and layers    |
      +--------------------------------+--------------------+-----------------------+
      | ``default_surface``            | Boolean            | ``True`` if all       |
      |                                |                    | columns have default  |
      |                                |                    | surface elevation     |
      +--------------------------------+--------------------+-----------------------+
      | ``extra_connections``          | set                | connections defined   |
      |                                |                    | between columns that  |
      |                                |                    | are not against each  |
      |                                |                    | other                 |
      +--------------------------------+--------------------+-----------------------+
      | ``filename``                   | string             | file name on disk     |
      +--------------------------------+--------------------+-----------------------+
      | ``gdcx``, ``gdcy``             | float              | cosines of angles x-  |
      |                                |                    | and y-axes make with  |
      |                                |                    | gravity vector        |
      +--------------------------------+--------------------+-----------------------+
      | ``node_kdtree``                | ``cKDTree``        | tree structure for    |
      |                                |                    | fast searching for    |
      |                                |                    | nodes                 |
      +--------------------------------+--------------------+-----------------------+
      | ``layerlist``                  | list               | layers (by index)     |
      +--------------------------------+--------------------+-----------------------+
      | ``layermesh``                  | ``layermesh`` mesh | Layermesh library     |
      |                                |                    | mesh object           |
      +--------------------------------+--------------------+-----------------------+
      | ``layer``                      | dictionary         | layers (by name)      |
      +--------------------------------+--------------------+-----------------------+
      | ``min_surface_block_thickness``| (float, string)    | thickness of thinnest |
      |                                |                    | surface block (and    |
      |                                |                    | associated column     |
      |                                |                    | name)                 |
      +--------------------------------+--------------------+-----------------------+
      | ``missing_connections``        | set                | missing connections   |
      |                                |                    | between columns       |
      +--------------------------------+--------------------+-----------------------+
      | ``nodelist``                   | list               | nodes (by index)      |
      +--------------------------------+--------------------+-----------------------+
      | ``node``                       | dictionary         | nodes (by name)       |
      +--------------------------------+--------------------+-----------------------+
      | ``num_atmosphere_blocks``      | integer            | number of atmosphere  |
      |                                |                    | blocks                |
      +--------------------------------+--------------------+-----------------------+
      | ``num_blocks``                 | integer            | total number of       |
      |                                |                    | blocks in the grid    |
      +--------------------------------+--------------------+-----------------------+
      | ``num_block_connections``      | integer            | total number of block |
      |                                |                    | connections in the    |
      |                                |                    | grid                  |
      +--------------------------------+--------------------+-----------------------+
      | ``num_columns``                | integer            | number of columns     |
      +--------------------------------+--------------------+-----------------------+
      | ``num_connections``            | integer            | number of connections |
      |                                |                    | between columns       |
      +--------------------------------+--------------------+-----------------------+
      | ``num_layers``                 | integer            | number of layers      |
      +--------------------------------+--------------------+-----------------------+
      | ``num_nodes``                  | integer            | number of nodes       |
      +--------------------------------+--------------------+-----------------------+
      | ``num_underground_blocks``     | integer            | number of             |
      |                                |                    | non-atmosphere blocks |
      +--------------------------------+--------------------+-----------------------+
      | ``num_wells``                  | integer            | number of wells       |
      +--------------------------------+--------------------+-----------------------+
      | ``orphans``                    | set                | orphaned nodes (nodes |
      |                                |                    | not belonging to any  |
      |                                |                    | column)               |
      +--------------------------------+--------------------+-----------------------+
      | ``permeability_angle``         | float              | rotation angle        |
      |                                |                    | (degrees              |
      |                                |                    | anticlockwise) of     |
      |                                |                    | first horizontal      |
      |                                |                    | permeability          |
      |                                |                    | direction             |
      +--------------------------------+--------------------+-----------------------+
      | ``read_function``              | dictionary         | dictionary of         |
      |                                |                    | functions used to     |
      |                                |                    | read data from file   |
      +--------------------------------+--------------------+-----------------------+
      | ``type``                       | string             | type of geometry      |
      |                                |                    | (currently only       |
      |                                |                    | 'GENER' supported)    |
      +--------------------------------+--------------------+-----------------------+
      | ``unit_type``                  | string             | distance unit (blank  |
      |                                |                    | for metres, 'FEET'    |
      |                                |                    | for ft)               |
      +--------------------------------+--------------------+-----------------------+
      | ``welllist``                   | list               | wells (by index)      |
      +--------------------------------+--------------------+-----------------------+
      | ``well``                       | dictionary         | wells (by name)       |
      +--------------------------------+--------------------+-----------------------+

Grid diagnostics
^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; diagnostics

A ``mulgrid`` object has some properties (and methods) for evaluating
its integrity. The property ``column_angle_ratio`` returns an
``np.array`` of the 'angle ratio' for each column (the ratio of
largest to smallest interior angles - see
:ref:`column objects <columnobjects>`), a measure of skewness.
The ``column_side_ratio`` returns an ``np.array`` of the 'side ratio'
for each column (the ratio of largest to smallest side length), a
measure of elongation. These array properties can be plotted using the
``layer_plot`` method (see :ref:`mulgrid methods <mulgridmethods>`) for a
graphical overview of grid quality.

There is also a ``connection_angle_cosine`` property, which returns an
``np.array`` of the angle cosine for each connection (the cosine of the
angle between a line joining the nodes in the connection and a line
joining the centres of the blocks in the connection). In general it is
desirable for these lines to be as close to perpendicular as possible,
making the cosines close to zero.

The ``bad_columns``, ``bad_layers``, ``missing_connections``,
``extra_connections`` and ``orphans`` properties return actual problems
with the grid which should be fixed. A summary of all these problems is
given by the :ref:`check() <sec:mulgrid:check>` method).

Blocks at the ground surface that have very small vertical thickness can
sometimes cause problems. The ``min_surface_block_thickness`` property
gives a tuple containing the minimum surface block thickness and the
name of the column in which it occurs. Thin surface blocks of this type
can be eliminated using the ``snap_columns_to_layers()`` method.

.. _mulgridreadfunctions:

Functions for reading data from file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; file format
.. index:: MULgraph geometry; reading

A ``mulgrid`` object has a ``read_function`` property which controls how
data are read from file. This property is a dictionary with six keys:
'd', 'f', 'e', 'g', 's' and 'x', denoting respectively integer, float,
exponential, general, string and blank. Each item in the dictionary is a
function which converts a string from the file on disk into the
appropriate value. For example, ``read_function['f']`` converts a string
to a floating point value. By default, the built-in Python ``float``
function is used for this (although it is modified slightly so that it
returns ``None`` if the input string is blank). There is a dictionary of
default reading functions included in PyTOUGH, called
``default_read_function``.

However, the user can specify other functions if needed. In particular,
files produced from Fortran programs sometimes have formatting that is
not readable by the default functions, if some more exotic Fortran
formatting options have been used. For example, a 'd' can also be used
to represent an exponent (like 'e'), or spaces can be included within a
number, or the exponent identifier (e.g. 'e') can be omitted. PyTOUGH
includes a second set of reading functions, called
``fortran_read_function``, for handling Fortran formatting. These are
slightly slower than the default reading functions.

The reading functions for a ``mulgrid`` object can be specified when the
object is being created, e.g.:

::

   geo = mulgrid('geom.dat', read_function = fortran_read_function)

.. _sec:mulgrid:blockordering:

Block ordering schemes
^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; block ordering

By default, the blocks in a TOUGH2 grid created from a ``mulgrid``
geometry are ordered by layer, from the atmosphere down to the bottom of
the model, with the blocks within each layer ordered by column
(following the ordering of the ``columnlist`` property, which is the
same as the column order specified in the geometry file).

It is also possible to sort the blocks according to their geometrical
type (8-node hexahedrons and 6-node wedges, corresponding to 4-node or
3-node columns respectively). This is useful for exporting the model to
Waiwera, which uses the PETSc DMPlex mesh representation, which sorts
cells by cell type in this way.

This can be done by setting the ``block_order`` property of the
geometry. This can be set when the ``mulgrid`` object is created or read
from file, as an optional parameter, e.g.:

::

   geo = mulgrid('geom.dat', block_order = 'dmplex')

It can also be specified after creation. The ``block_order`` property is
a string which can take the value **'layer_column'** for layer/column
block ordering, or **'dmplex'** if the blocks are to be sorted by
geometrical type. It can also take the value ``None`` which gives the
default layer/column ordering.

::

   geo.block_order = 'layer_column'

The block ordering scheme can be stored in the MULgraph geometry file,
via an integer flag in the header (see
:ref:`MULgraph geometry file format <geometry_file_format>`). This flag
is an extension to the original MULgraph geometry file format. If a
``mulgrid`` object is created by reading a file in which this flag is
not present, its ``block_order`` property will be ``None``, in which
case the default layer/column ordering will be used. When a geometry
file is read in, and a block ordering is specified via the
``block_order`` parameter, this will override any block ordering
stored in the file.

Tilted geometries
^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; tilting

Non-horizontal (i.e. tilted) geometries can be constructed by setting
the ``mulgrid`` properties ``gdcx`` and ``gdcy`` non-zero. These
properties represent the cosines of the angles the x- and y-axes make
with the gravity vector. By default they are both zero, giving a
horizontal grid. A geometry with ``gdcx`` = 1 can be used to construct a
2-D vertical slice grid with a non-layered structure. When a ``t2grid``
object is created from a tilted geometry, e.g. using the ``t2grid``
:ref:`fromgeo() <sec:t2grid:fromgeo>` method, only the
gravity cosines of the connections are affected (the ``dircos`` property
of each connection).

Rotating permeability directions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; permeability directions

It is possible to rotate the permeability principal directions of a
``mulgrid`` object with respect to the coordinate axes- for example, to
align permeabilities with a dominant fault direction- by specifying the
``permeability_angle`` property. When a ``t2grid`` object is created,
e.g. using the ``t2grid`` :ref:`fromgeo() <sec:t2grid:fromgeo>`
method, this can change the ``direction`` property of each
connection.

Conversion to and from Layermesh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; Layermesh conversion
.. index:: Layermesh

A ``mulgrid`` geometry may be converted to a `Layermesh
<https://github.com/acroucher/layermesh>`_ mesh simply by accessing
its ``layermesh`` property. Layermesh is a dedicated library for
general layer/column meshes. Its mesh objects have capabilities
similar to those of a ``mulgrid`` object, but it has advantages such
as higher efficiency and a simpler interface. The Layermesh library
must be installed before this property can be used.

Example:

::

     geo = mulgrid('gmymesh.dat')
     m = geo.layermesh # m is a Layermesh mesh object

Conversely, a Layermesh object can be imported into a ``mulgrid`` object
using the :ref:`from_layermesh() <sec:mulgrid:from_layermesh>` method.
                           
.. _mulgridmethods:

Methods
~~~~~~~

The main methods of a ``mulgrid`` object are listed in the following
:ref:`table <tb:mulgrid_methods>`. Details of these methods are
given below.

.. container::
   :name: tb:mulgrid_methods

   .. table:: Methods of a ``mulgrid`` object

      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | **Method**                                                                        | **Type**                     | **Description**      |
      +===================================================================================+==============================+======================+
      | :ref:`add_column <sec:mulgrid:add_column>`                                        | –                            | adds a column to the |
      |                                                                                   |                              | grid                 |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`add_connection <sec:mulgrid:add_connection>`                                | –                            | adds a connection to |
      |                                                                                   |                              | the grid             |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`add_layer <sec:mulgrid:add_layer>`                                          | –                            | adds a layer to the  |
      |                                                                                   |                              | grid                 |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`add_node <sec:mulgrid:add_node>`                                            | –                            | adds a node to the   |
      |                                                                                   |                              | grid                 |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`add_well <sec:mulgrid:add_well>`                                            | –                            | adds a well to the   |
      |                                                                                   |                              | grid                 |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`block_centre <sec:mulgrid:block_centre>`                                    | ``np.array``                 | block centre         |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`block_contains_point <sec:mulgrid:block_contains_point>`                    | Boolean                      | whether a block      |
      |                                                                                   |                              | contains a 3D point  |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`block_mapping <sec:mulgrid:block_mapping>`                                  | dictionary                   | mapping from the     |
      |                                                                                   |                              | blocks of another    |
      |                                                                                   |                              | ``mulgrid`` object   |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`block_name <sec:mulgrid:block_name>`                                        | string                       | name of block at     |
      |                                                                                   |                              | given layer and      |
      |                                                                                   |                              | column               |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`block_name_containing_point <sec:mulgrid:block_name_containing_point>`      | string                       | name of block        |
      |                                                                                   |                              | containing specified |
      |                                                                                   |                              | point                |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`block_surface <sec:mulgrid:block_surface>`                                  | float                        | block top elevation  |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`block_volume <sec:mulgrid:block_volume>`                                    | float                        | block volume         |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`check <sec:mulgrid:check>`                                                  | Boolean                      | checks grid for      |
      |                                                                                   |                              | errors (and          |
      |                                                                                   |                              | optionally fixes     |
      |                                                                                   |                              | them)                |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`column_boundary_nodes <sec:mulgrid:column_boundary_nodes>`                  | list                         | nodes around the     |
      |                                                                                   |                              | outer boundary of a  |
      |                                                                                   |                              | group of columns     |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`column_bounds <sec:mulgrid:column_bounds>`                                  | list                         | bounding rectangle   |
      |                                                                                   |                              | around a list of     |
      |                                                                                   |                              | columns              |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`column_containing_point <sec:mulgrid:column_containing_point>`              | column                       | column containing    |
      |                                                                                   |                              | specified horizontal |
      |                                                                                   |                              | point                |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`column_mapping <sec:mulgrid:column_mapping>`                                | dictionary                   | mapping from the     |
      |                                                                                   |                              | columns of another   |
      |                                                                                   |                              | object               |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`column_name <sec:mulgrid:column_name>`                                      | string                       | column name of a     |
      |                                                                                   |                              | block name           |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`column_neighbour_groups <sec:mulgrid:column_neighbour_groups>`              | list                         | groups connected     |
      |                                                                                   |                              | columns              |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`column_quadtree <sec:mulgrid:column_quadtree>`                              | quadtree                     | quadtree structure   |
      |                                                                                   |                              | for searching        |
      |                                                                                   |                              | columns              |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`column_surface_layer <sec:mulgrid:column_surface_layer>`                    | :ref:`layer <layerobjects>`  | surface layer for a  |
      |                                                                                   |                              | specified column     |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`column_values <sec:mulgrid:column_values>`                                  | tuple                        | values of a variable |
      |                                                                                   |                              | down a column        |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`columns_in_polygon <sec:mulgrid:columns_in_polygon>`                        | list                         | columns inside a     |
      |                                                                                   |                              | specified polygon    |
      |                                                                                   |                              | (or rectangle)       |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`connects <sec:mulgrid:connects>`                                            | Boolean                      | whether the grid has |
      |                                                                                   |                              | a connection between |
      |                                                                                   |                              | two specified        |
      |                                                                                   |                              | columns              |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`copy_layers_from <sec:mulgrid:copy_layers_from>`                            | –                            | copies layer         |
      |                                                                                   |                              | structure from       |
      |                                                                                   |                              | another geometry     |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`copy_wells_from <sec:mulgrid:copy_wells_from>`                              | –                            | copies wells from    |
      |                                                                                   |                              | another geometry     |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`decompose_columns <sec:mulgrid:decompose_columns>`                          | –                            | decomposes columns   |
      |                                                                                   |                              | into triangles and   |
      |                                                                                   |                              | quadrilaterals       |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`delete_column <sec:mulgrid:delete_column>`                                  | –                            | deletes a column     |
      |                                                                                   |                              | from the grid        |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`delete_connection <sec:mulgrid:delete_connection>`                          | –                            | deletes a connection |
      |                                                                                   |                              | from the grid        |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`delete_layer <sec:mulgrid:delete_layer>`                                    | –                            | deletes a layer from |
      |                                                                                   |                              | the grid             |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`delete_node <sec:mulgrid:delete_node>`                                      | –                            | deletes a node from  |
      |                                                                                   |                              | the grid             |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`delete_orphans <sec:mulgrid:delete_orphans>`                                | –                            | deletes any orphaned |
      |                                                                                   |                              | nodes from the grid  |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`delete_orphan_wells <sec:mulgrid:delete_orphan_wells>`                      | –                            | deletes any orphaned |
      |                                                                                   |                              | wells from the grid  |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`delete_well <sec:mulgrid:delete_well>`                                      | –                            | deletes a well from  |
      |                                                                                   |                              | the grid             |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`empty <sec:mulgrid:empty>`                                                  | –                            | empties contents of  |
      |                                                                                   |                              | grid                 |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`export_surfer <sec:mulgrid:export_surfer>`                                  | –                            | exports to various   |
      |                                                                                   |                              | files on disk for    |
      |                                                                                   |                              | visualization in     |
      |                                                                                   |                              | Surfer               |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`fit_columns <sec:mulgrid:fit_columns>`                                      | ``np.array`` or              | fits scattered data  |
      |                                                                                   | dictionary                   | to column centres    |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`fit_surface <sec:mulgrid:fit_surface>`                                      | –                            | fits column surface  |
      |                                                                                   |                              | elevations from data |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`from_amesh <sec:mulgrid:from_amesh>`                                        | (:ref:`mulgrid <mulgrids>`,  | creates Voronoi      |
      |                                                                                   | dict)                        | geometry from AMESH  |
      |                                                                                   |                              | grid                 |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`from_gmsh <sec:mulgrid:from_gmsh>`                                          | :ref:`mulgrid <mulgrids>`    | creates geometry     |
      |                                                                                   |                              | from a ``gmsh`` grid |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`from_layermesh <sec:mulgrid:from_layermesh>`                                | :ref:`mulgrid <mulgrids>`    | creates geometry     |
      |                                                                                   |                              | from a ``Layermesh`` |
      |                                                                                   |                              | grid                 |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`layer_containing_elevation <sec:mulgrid:layer_containing_elevation>`        | layer                        | layer containing     |
      |                                                                                   |                              | specified vertical   |
      |                                                                                   |                              | elevation            |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`layer_mapping <sec:mulgrid:layer_mapping>`                                  | dictionary                   | mapping from the     |
      |                                                                                   |                              | layers of another    |
      |                                                                                   |                              | object               |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`layer_name <sec:mulgrid:layer_name>`                                        | string                       | layer name of a      |
      |                                                                                   |                              | block name           |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`layer_plot <sec:mulgrid:layer_plot>`                                        | –                            | plots a variable     |
      |                                                                                   |                              | over a layer of the  |
      |                                                                                   |                              | grid                 |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`line_plot <sec:mulgrid:line_plot>`                                          | –                            | plots a variable     |
      |                                                                                   |                              | along an arbitrary   |
      |                                                                                   |                              | line through the     |
      |                                                                                   |                              | grid                 |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`line_values <sec:mulgrid:line_values>`                                      | tuple                        | values of a variable |
      |                                                                                   |                              | along an arbitrary   |
      |                                                                                   |                              | line through the     |
      |                                                                                   |                              | grid                 |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`meshio_grid <sec:mulgrid:meshio_grid>`                                      | tuple                        | mesh in ``meshio``   |
      |                                                                                   |                              | format               |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`minc_array <sec:mulgrid:minc_array>`                                        | array                        | values for a         |
      |                                                                                   |                              | particular level in  |
      |                                                                                   |                              | a MINC grid          |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`nodes_in_columns <sec:mulgrid:nodes_in_columns>`                            | list                         | nodes in a specified |
      |                                                                                   |                              | list of columns      |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`nodes_in_polygon <sec:mulgrid:nodes_in_polygon>`                            | list                         | nodes inside a       |
      |                                                                                   |                              | specified polygon    |
      |                                                                                   |                              | (or rectangle)       |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`node_nearest_to <sec:mulgrid:node_nearest_to>`                              | :ref:`node <nodeobjects>`    | node nearest to a    |
      |                                                                                   |                              | specified point      |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`optimize <sec:mulgrid:optimize>`                                            | –                            | adjusts node         |
      |                                                                                   |                              | positions to         |
      |                                                                                   |                              | optimize grid        |
      |                                                                                   |                              | quality              |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`polyline_values <sec:mulgrid:polyline_values>`                              | tuple                        | values of a variable |
      |                                                                                   |                              | along an arbitrary   |
      |                                                                                   |                              | polyline through the |
      |                                                                                   |                              | grid                 |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`read <sec:mulgrid:read>`                                                    | :ref:`mulgrid <mulgrids>`    | reads geometry file  |
      |                                                                                   |                              | from disk            |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`rectangular <sec:mulgrid:rectangular>`                                      | :ref:`mulgrid <mulgrids>`    | creates rectangular  |
      |                                                                                   |                              | grid                 |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`reduce <sec:mulgrid:reduce>`                                                | –                            | reduces a grid to    |
      |                                                                                   |                              | contain only         |
      |                                                                                   |                              | specified columns    |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`refine <sec:mulgrid:refine>`                                                | –                            | refines specified    |
      |                                                                                   |                              | columns in the grid  |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`refine_layers <sec:mulgrid:refine_layers>`                                  | –                            | refines specified    |
      |                                                                                   |                              | layers in the grid   |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`rename_column <sec:mulgrid:rename_column>`                                  | Boolean                      | renames a column     |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`rename_layer <sec:mulgrid:rename_layer>`                                    | Boolean                      | renames a layer      |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`rotate <sec:mulgrid:rotate>`                                                | –                            | rotates a grid in    |
      |                                                                                   |                              | the horizontal plane |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`slice_plot <sec:mulgrid:slice_plot>`                                        | –                            | plots a variable     |
      |                                                                                   |                              | over a vertical      |
      |                                                                                   |                              | slice through the    |
      |                                                                                   |                              | grid                 |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`snap_columns_to_layers <sec:mulgrid:snap_columns_to_layers>`                | –                            | snaps column         |
      |                                                                                   |                              | surfaces to layer    |
      |                                                                                   |                              | bottoms              |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`snap_columns_to_nearest_layers <sec:mulgrid:snap_columns_to_nearest_layers>`| –                            | snaps column         |
      |                                                                                   |                              | surfaces to nearest  |
      |                                                                                   |                              | layer elevations     |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`split_column <sec:mulgrid:split_column>`                                    | Boolean                      | splits a             |
      |                                                                                   |                              | quadrilateral column |
      |                                                                                   |                              | into two triangles   |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`translate <sec:mulgrid:translate>`                                          | –                            | moves a grid by      |
      |                                                                                   |                              | simple translation   |
      |                                                                                   |                              | in 3D                |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`well_values <sec:mulgrid:well_values>`                                      | tuple                        | values of a variable |
      |                                                                                   |                              | down a well          |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`write <sec:mulgrid:write>`                                                  | –                            | writes to geometry   |
      |                                                                                   |                              | file on disk         |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`write_bna <sec:mulgrid:write_bna>`                                          | –                            | writes to Atlas BNA  |
      |                                                                                   |                              | file on disk         |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`write_exodusii <sec:mulgrid:write_exodusii>`                                | –                            | writes to ExodusII   |
      |                                                                                   |                              | file on disk         |
      |                                                                                   |                              |                      |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`write_mesh <sec:mulgrid:write_mesh>`                                        | –                            | writes to mesh file  |
      |                                                                                   |                              | (various formats) on |
      |                                                                                   |                              | disk                 |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+
      | :ref:`write_vtk <sec:mulgrid:write_vtk>`                                          | –                            | writes to VTK file   |
      |                                                                                   |                              | on disk              |
      |                                                                                   |                              |                      |
      +-----------------------------------------------------------------------------------+------------------------------+----------------------+

.. _sec:mulgrid:add_column:

``add_column(col)``
^^^^^^^^^^^^^^^^^^^

Adds a :ref:`column <columnobjects>` object ``col`` to the grid.
If a column with the same name already exists, no new column is added.

----

.. _sec:mulgrid:add_connection:

``add_connection(con)``
^^^^^^^^^^^^^^^^^^^^^^^

Adds a :ref:`connection <connectionobjects>` object ``con`` to
the grid. If a connection with the same name already exists, no new
connection is added.

----

.. _sec:mulgrid:add_layer:

``add_layer(lay)``
^^^^^^^^^^^^^^^^^^

Adds a :ref:`layer <layerobjects>` object ``lay`` to the grid.
If a layer with the same name already exists, no new layer is added.

----

.. _sec:mulgrid:add_node:

``add_node(n)``
^^^^^^^^^^^^^^^

Adds a :ref:`node <nodeobjects>` object ``n`` to the grid. If a
node with the same name already exists, no new node is added.

----

.. _sec:mulgrid:add_well:

``add_well(w)``
^^^^^^^^^^^^^^^

Adds a :ref:`well <wellobjects>` object ``w`` to the grid. If a
well with the same name already exists, no new well is added.

----

.. _sec:mulgrid:block_contains_point:

``block_contains_point(blockname, pos)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; searching

Returns ``True`` if the grid block with the given name contains the 3D
point ``pos``.

**Parameters:**

-  | **blockname**: string
   | The name of the block.

-  | **pos**: ``np.array``
   | 3-element array representing the 3D point.

----
   
.. _sec:mulgrid:block_centre:

``block_centre(lay, col)``
^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns the centre of the block corresponding to the given layer and
column.

The horizontal centre is given by the column centre. The vertical centre
is given by the layer centre, except for surface blocks with column
surface lower than the layer top, in which case it is the midpoint
between the column surface and the layer bottom. (For surface blocks
with column surface higher than the layer top, the vertical centre is
still the layer centre, to give a uniform pressure reference.)

**Parameters:**

-  | **lay**: :ref:`layer <layerobjects>` or string
   | The specified layer or layer name.

-  | **col**: :ref:`column <columnobjects>` or string
   | The specified column or column name.

----
   
.. _sec:mulgrid:block_mapping:

``block_mapping(geo, column_mapping=False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; block mappings

Returns a dictionary mapping each block name in the ``mulgrid`` object
``geo`` to the name of the nearest block in the object's own geometry.
Can optionally also return the associated column mapping.

**Parameters:**

-  | **geo**: :ref:`mulgrid <mulgrids>` 
   | The ``mulgrid`` object to create a block mapping from.

-  | **column_mapping**: Boolean
   | If ``True``, the column mapping will also be returned (i.e. the
     function will return a tuple containing the block mapping and the
     column mapping). Default value is ``False``.

----
     
.. _sec:mulgrid:block_name:

``block_name(layer_name, column_name, blockmap = {})``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Gives the name of the block corresponding to the specified layer and
column names, according to the naming convention of the grid.

An optional block name mapping can be applied.

**Parameters:**

-  | **layer_name**, **column_name**: string
   | Name of layer and column (the widths of these strings are
     determined by the grid's naming convention).

-  | **blockmap**: dictionary
   | Dictionary mapping the block names in the geometry to another block
     naming system. This dictionary need not contain entries for all
     blocks in the geometry- those not included in the mapping will not
     be altered.

----
     
.. _sec:mulgrid:block_name_containing_point:

``block_name_containing_point(pos, qtree=None, blockmap={})``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; searching
.. index:: MULgraph geometry; quadtree

Gives the name of the block containing a specified 3-D position in the
grid (returns ``None`` if the point lies outside the grid).

**Parameters:**

-  | **pos**: ``np.array``
   | Position of point in 3-D

-  | **qtree**: ``quadtree``
   | Quadtree object for fast searching of grid columns (can be
     constructed using the :ref:`column_quadtree() <sec:mulgrid:column_quadtree>`
      method).

-  | **blockmap**: dictionary
   | Dictionary mapping the block names in the geometry to another block
     naming system.

----
     
.. _sec:mulgrid:block_surface:

``block_surface(lay, col)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns the elevation of the top surface of the block corresponding to
the given layer and column.

**Parameters:**

-  | **lay**: :ref:`layer <layerobjects>`
   | The specified layer.

-  | **col**: :ref:`column <columnobjects>`
   | The specified column.

----
   
.. _sec:mulgrid:block_volume:

``block_volume(lay, col)``
^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns the volume of the block corresponding to the given layer and
column.

**Parameters:**

-  | **lay**: :ref:`layer <layerobjects>`
   | The specified layer.

-  | **col**: :ref:`column <columnobjects>`
   | The specified column.

----
   
.. _sec:mulgrid:check:

``check(fix=False,silent=False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; checking

Checks a grid for errors and optionally fixes them. Errors checked for
are: missing connections, extra connections, orphaned nodes, and columns
and layers that do not contain their own centres. Returns ``True`` if no
errors were found, and ``False`` otherwise. If ``fix`` is ``True``, any
identified problems will be fixed. If ``silent`` is ``True``, there is
no printout (only really useful if ``fix`` is ``True``).

**Parameters:**

-  | **fix**: Boolean
   | Whether to fix any problems identified.

-  | **silent**: Boolean
   | Whether to print out feedback or not.

----
   
.. _sec:mulgrid:column_boundary_nodes:

``column_boundary_nodes(columns)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; searching

Returns the nodes around the outer boundary of a list of columns. The
list is ordered, in a counter-clockwise direction.

**Parameters:**

-  | **columns**: list
   | The list of columns for which the boundary is required.

----
   
.. _sec:mulgrid:column_bounds:

``column_bounds(columns)``
^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns a bounding rectangle around a list of columns.

**Parameters:**

-  | **columns**: list
   | The list of columns for which the bounds are required.

----
   
.. _sec:mulgrid:column_containing_point:

``column_containing_point(pos, columns=None, guess=None, bounds=None, qtree=None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; searching
.. index:: MULgraph geometry; quadtree

Returns the grid column containing the specified horizontal point. If
``columns`` is specified, only columns in the given list will be
searched. An initial ``guess`` column can optionally be specified. If
``bounds`` is specified, points outside the given polygon will always
return ``None``. A quadtree structure can also be specified to speed up
searching.

**Parameters:**

-  | **pos**: ``np.array``
   | Horizontal position (*x*, *y*)

-  | **columns**: list of :ref:`column <columnobjects>` (or
     ``None``)
   | List of columns to search. If ``None``, the entire grid will be
     searched.

-  | **guess**: :ref:`column <columnobjects>` (or ``None``)
   | Guess of required column. If specified, this column will be tested
     first, followed (if necessary) by its neighbours; only if none of
     these contain the point will the remaining columns be searched.
     This can speed up the process if data follow a sequential pattern
     in space, e.g. a grid or lines.

-  | **bounds**: list of ``np.array`` (or ``None``)
   | Polygon or rectangle representing e.g. the boundary of the grid:
     points outside this polygon will always return ``None``. If the
     polygon has only two points, it will be interpreted as a rectangle
     [bottom left, top right].

-  | **qtree**: ``quadtree``
   | A quadtree object for searching the columns of the grid. If many
     points are to be located, this option can speed up the search. The
     quadtree can be constructed before searching using the
     :ref:`column_quadtree() <sec:mulgrid:column_quadtree>`
     method.

----
     
.. _sec:mulgrid:column_mapping:

``column_mapping(geo)``
^^^^^^^^^^^^^^^^^^^^^^^

Returns a dictionary mapping each column name in the ``mulgrid`` object
``geo`` to the name of the nearest column in the object's own geometry.
If the SciPy library is available, a KDTree structure is used to speed
searching.

**Parameters:**

-  | **geo**: :ref:`mulgrid <mulgrids>`
   | The ``mulgrid`` object to create a column mapping from.

----
   
.. _sec:mulgrid:column_name:

``column_name(block_name)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Gives the name of the column corresponding to the specified block name,
according to the naming convention of the grid.

**Parameters:**

-  | **block_name**: string
   | Block name.

----
   
.. _sec:mulgrid:column_neighbour_groups:

``column_neighbour_groups(columns)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From the given list or set of columns, finds sets of columns that are
connected together, and returns a list of them.

**Parameters:**

-  | **columns**: list or set
   | List or set of columns to group.

----
   
.. _sec:mulgrid:column_quadtree:

``column_quadtree(columns=None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; quadtree

Returns a quadtree structure for fast searching of grid columns, to find
which column a given point lies in. This can then be passed into various
other ``mulgrid`` methods that do such searching, e.g.
:ref:`block_name_containing_point() <sec:mulgrid:block_name_containing_point>` or
:ref:`well_values() <sec:mulgrid:well_values>`, to speed them
up (useful for large grids).

The quadtree is an instance of a ``quadtree`` class, defined in the
``mulgrids`` module.

**Parameters:**

-  | **columns**: list (or ``None``)
   | A list of columns in the grid, specifying the search area. This
     parameter can be used to further speed searching if it is only
     necessary to search columns in a defined area. If ``None``, the
     search area is the whole grid (all columns).

----
     
.. _sec:mulgrid:column_surface_layer:

``column_surface_layer(col)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns the layer containing the surface elevation of a specified
column.

**Parameters:**

-  | **col**: :ref:`column <columnobjects>`
   | The column for which the surface layer is to be found.

----
   
.. _sec:mulgrid:column_values:

``column_values(col, variable, depth = False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns values of a specified variable down a specified column. The
variable can be a list or ``np.array`` containing a value for every
block in the grid.

The routine returns a tuple of two arrays (``d``,\ ``v``), the first
(``d``) containing the elevation (or depth from surface if the ``depth``
parameter is set to ``True``), and the second (``v``) containing the
value of the variable at each block in the column.

**Parameters:**

-  | **col**: :ref:`column <columnobjects>` or string
   | The column for which values are to be found.

-  | **variable**: list (or ``np.array``)
   | Values of variable, of length equal to the number of blocks in the
     grid.

-  | **depth**: Boolean
   | Set to ``True`` to give depths from surface, instead of elevations,
     as the first returned array.

----
     
.. _sec:mulgrid:columns_in_polygon:

``columns_in_polygon(polygon)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; searching

Returns a list of all columns with centres inside the specified polygon
or rectangle.

**Parameters:**

-  | **polygon**: list (of ``np.array``)
   | List of points defining the polygon (each point is a two-element
     ``np.array``). If the list has only two points, it will be
     interpreted as a rectangle [bottom left, top right].

----
     
.. _sec:mulgrid:connects:

``connects(column1, column2)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns ``True`` if the geometry contains a connection connecting the
two specified columns.

**Parameters:**

-  | **column1, column2**: :ref:`column <columnobjects>` 
   | Two columns in the geometry.

----
   
.. _sec:mulgrid:copy_layers_from:

``copy_layers_from(geo)``
^^^^^^^^^^^^^^^^^^^^^^^^^

Copies the layer structure from the geometry ``geo`` (deleting any
existing layers first).

**Parameters:**

-  | **geo**: :ref:`mulgrid <mulgrids>`
   | The geometry to copy layers from.

----
   
.. _sec:mulgrid:copy_wells_from:

``copy_wells_from(geo)``
^^^^^^^^^^^^^^^^^^^^^^^^

Copies the wells from the geometry ``geo`` (deleting any existing wells
first).

**Parameters:**

-  | **geo**: :ref:`mulgrid <mulgrids>`
   | The geometry to copy wells from.

----
   
.. _sec:mulgrid:decompose_columns:

``decompose_columns(columns = [], mapping = False, chars = ascii_lowercase)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Decomposes columns with more than four sides into triangular and
quadrilateral columns. This can be useful when carrying out calculations
on the geometry that rely on finite element methods (e.g. the
``fit_columns()`` method uses it).

In general, columns are decomposed by adding a node at the column
centroid and forming triangles around it. However, there are special
cases for columns with lower numbers of sides (less than 9) and
'straight' nodes, i.e. nodes on a straight line between their
neighbouring nodes in the column). These make use of simpler
decompositions.

**Parameters:**

-  | **columns**: list
   | List of columns to be decomposed. If the list is empty (the
     default), all columns are decomposed.

-  | **mapping**: Boolean
   | If ``True``, return a dictionary mapping each original column name
     to a list of decomposed columns that replace it.

-  | **chars**: string
   | Specifies a string of characters to use when forming new node and
     column names. Default is lowercase letters.

----
     
.. _sec:mulgrid:delete_column:

``delete_column(colname)``
^^^^^^^^^^^^^^^^^^^^^^^^^^

Deletes the column with the specified name from the grid.

**Parameters:**

-  | **colname**: string
   | Name of the column to be deleted.

----

.. _sec:mulgrid:delete_connection:

``delete_connection(colnames)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Deletes the connection between the specified columns from the grid.

**Parameters:**

-  | **colnames**: tuple of string
   | Tuple of two column names.

----

.. _sec:mulgrid:delete_layer:

``delete_layer(layername)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Deletes the layer with the specified name from the grid.

**Parameters:**

-  | **layername**: string
   | Name of the layer to be deleted.

----

.. _sec:mulgrid:delete_node:

``delete_node(nodename)``
^^^^^^^^^^^^^^^^^^^^^^^^^

Deletes the node with the specified name from the grid.

**Parameters:**

-  | **nodename**: string
   | Name of the node to be deleted.

----

.. _sec:mulgrid:delete_orphans:

``delete_orphans()``
^^^^^^^^^^^^^^^^^^^^

Deletes any orphaned nodes (those not belonging to any column) from the
grid.

----

.. _sec:mulgrid:delete_orphan_wells:

``delete_orphan_wells()``
^^^^^^^^^^^^^^^^^^^^^^^^^

Deletes any orphaned wells (those with wellheads outside the grid).

----

.. _sec:mulgrid:delete_well:

``delete_well(wellname)``
^^^^^^^^^^^^^^^^^^^^^^^^^

Deletes the well with the specified name from the grid.

**Parameters:**

-  | **layername**: string
   | Name of the layer to be deleted.

----

.. _sec:mulgrid:empty:

``empty()``
^^^^^^^^^^^

.. index:: MULgraph geometry; emptying

Empties the grid of all its nodes, columns, layers, wells and
connections. Other properties are unaffected.

----

.. _sec:mulgrid:export_surfer:

``export_surfer(filename='', aspect=8.0, left=0.0)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; exporting

Exports the grid to files on disk useful for visualization in Surfer.
Six files are written out:

-  an Atlas BNA file (``filename.bna``) representing the grid columns

-  a CSV file (``filename_column_names.csv``) containing the column
   names

-  a Golden Software blanking file (``filename_layers.bln``) file
   representing the grid layers

-  a CSV file (``filename_layer_bottom_elevations.csv``) containing the
   bottom elevations of the layers

-  a CSV file (``filename_layer_centres.csv``) containing the elevations
   of the centres of the layers

-  a CSV file (``filename_layer_names.csv``) containing the names of the
   layers

**Parameters:**

-  | **filename**: string
   | Base name for the exported files. If it is not specified, the
     ``filename`` property of the ``mulgrid`` object itself is used
     (unless this is also blank, in which case a default name is used),
     with its extension removed.

-  | **aspect**: float
   | Aspect ratio for the layer plot, so that the width is the total
     height of the grid divided by ``aspect`` (default 8.0).

-  | **left**: float
   | Coordinate value of the left hand side of the layer plot (default
     zero).

----

.. _sec:mulgrid:fit_columns:

``fit_columns(data, alpha=0.1, beta=0.1, columns=[], min_columns=[], grid_boundary=False, silent=False, output_dict=False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; fitting data

Fits scattered data to column centres, using bilinear least-squares
finite element fitting with Sobolev smoothing. Smoothing is useful when
data density is low in some areas of the grid, in which case
least-squares fitting without smoothing can fail (e.g. if there are any
columns which do not contain any data points).

By default, this method returns an ``np.array`` with length given by the
number of columns to be fitted. Each value in the array represents the
fitted data value at the centre of the corresponding column. If the
``output_dict`` parameter is set to ``True``, a dictionary is returned,
with fitted values indexed by column names.

**Parameters:**

-  | **data**: ``np.array``
   | Two-dimensional array of data to fit. Each row of the array should
     contain the x,y co-ordinates for each data point, followed by the
     corresponding data value. Such an array can be conveniently read
     from a text file using the ``np.loadtxt()`` method.

-  | **alpha**: float
   | Smoothing parameter for first derivatives - increasing its value
     results in solutions with lower gradients (but may result in
     extrema being smoothed out).

-  | **beta**: float
   | Smoothing parameter for second derivatives - increasing its value
     results in solutions with lower curvature.

-  | **columns**: list of string or :ref:`column <columnobjects>`     
   | Columns, or names of columns to be fitted. If empty (the default),
     then all columns will be fitted.

-  | **min_columns**: list of string or :ref:`column <columnobjects>`     
   | Columns, or names of columns for which fitted data will be
     determined from the minimum of the fitted nodal values (fitted
     values at all other columns are determined from the average of the
     fitted nodal values).

-  | **grid_boundary**: Boolean
   | If ``True``, test each data point first to see if it lies inside
     the boundary polygon of the grid. This can speed up the fitting
     process if there are many data points outside the grid, and the
     grid has a simple boundary (e.g. a rectangle). In general if there
     are many data points outside the grid, it is best to clip the data
     set before fitting, particularly if it is to be used more than
     once.

-  | **silent**: Boolean
   | Set to ``True`` to suppress printing fitting progress.

-  | **output_dict**: Boolean
   | Set ``True`` to return results as a dictionary of fitted values
     indexed by column names, instead of an array.

----

.. _sec:mulgrid:fit_surface:

``fit_surface(data, alpha=0.1, beta=0.1, columns=[], min_columns=[], grid_boundary=False, layer_snap=0.0, silent=False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; fitting surface

Fits column surface elevations from data, using bilinear least-squares
finite element fitting with Sobolev smoothing (using the
:ref:`fit_columns() <sec:mulgrid:fit_columns>` method).
Smoothing is useful when data density is low in some areas of the grid,
in which case least-squares fitting without smoothing can fail (e.g. if
there are any columns which do not contain any data points). Use the
``layer_snap`` parameter to eliminate surface blocks with very small
thickness.

**Parameters:**

-  | **data**: ``np.array``
   | Two-dimensional array of data to fit. Each row of the array should
     contain the x,y,z values for each data point. Such an array can be
     conveniently read from a text file using the ``np.loadtxt()``
     method.

-  | **alpha**: float
   | Smoothing parameter for first derivatives - increasing its value
     results in solutions with lower gradients (but may result in
     extrema being smoothed out).

-  | **beta**: float
   | Smoothing parameter for second derivatives - increasing its value
     results in solutions with lower curvature.

-  | **columns**: list of string or :ref:`column <columnobjects>`     
   | Columns, or names of columns to be fitted. If empty (the default),
     then all columns will be fitted.

-  | **min_columns**: list of string or :ref:`column <columnobjects>`     
   | Columns, or names of columns for which elevations will be
     determined from the minimum of the fitted nodal elevations
     (elevations at all other columns are determined from the average of
     the fitted nodal elevations).

-  | **grid_boundary**: Boolean
   | If ``True``, test each data point first to see if it lies inside
     the boundary polygon of the grid. This can speed up the fitting
     process if there are many data points outside the grid, and the
     grid has a simple boundary (e.g. a rectangle). In general if there
     are many data points outside the grid, it is best to clip the data
     set before fitting, particularly if it is to be used more than
     once.

-  | **layer_snap**: float
   | Smallest desired surface block thickness. Set to a positive value
     to prevent columns being assigned surface elevations that are very
     close to the bottom of a layer (resulting in very thin surface
     blocks). Default value is zero (i.e. no layer snapping).

-  | **silent**: Boolean
   | Set to ``True`` to suppress printing fitting progress.

----

.. _sec:mulgrid:from_amesh:

``from_amesh(input_filename='in', segment_filename='segmt', convention=0, node_tolerance=None, justify='r', chars=ascii_lowercase, spaces=True, block_order=None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; Voronoi

Returns a ``mulgrid`` object (and a block mapping dictionary) from a
Voronoi mesh previously created by the `AMESH
<https://tough.lbl.gov/licensing-download/free-software-download/>`_
utility, or by other software that uses AMESH (e.g.  WinGridder or
Steinar).

The block naming convention for the output ``mulgrid`` object can be
specified via the ``convention`` parameter. Note that in general this
may not be the same as the block naming convention of the original mesh
created by AMESH. In fact, AMESH can create meshes with block naming
conventions that do not correspond to any of the
:ref:`MULgraph conventions <tb:mulgrid_conventions>`. This is why the
``from_amesh()`` method also returns a block mapping dictionary, which
maps block names in the ``mulgrid`` geometry to the block names in the
original AMESH grid.

The optional ``justify`` and ``case`` parameters control the formatting
of the character part of the block names. Additionally, the characters
used to form node/column or layer names can be specified using the
``chars`` parameter. (This can be useful for example for grids with
large numbers of nodes and/or columns, for which lowercase letters alone
may not be enough.)

The ``from_amesh()`` method assumes the original AMESH grid has layers
of constant thickness (i.e. all blocks in each layer of the AMESH input
file have the same specified thickness). Grids with layers of
non-constant thickness cannot be represented by a ``mulgrid`` object and
will cause an exception to be raised.

**Parameters:**

-  | **input_filename**: string
   | Filename for AMESH input file. Default is 'in'.

-  | **segment_filename**: string
   | Filename for AMESH output segment file. Default is 'segmt'.

-  | **convention**: integer
   | Naming convention for grid columns and layers.

-  | **node_tolerance**: float or ``None``
   | Horizontal tolerance for identifying distinct nodes in the segment
     file. If a node is read in with horizontal distance from an
     existing node less than the tolerance, then the two nodes are
     assumed to be identical. If ``None`` (the default), then the
     tolerance is set to 90% of the smallest segment length. If errors
     are encountered in identifying nodes belonging to the grid columns,
     it may be worth adjusting this parameter.

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

-  | **block_order**: string or ``None``
   | Specify ``None`` or 'layer_column' for default block ordering by
     layer and column, starting from the atmosphere. Specify 'dmplex' to
     order blocks by geometrical type (8-node hexahedrons first followed
     by 6-node wedges) as in PETSc DMPlex meshes.

----

.. _sec:mulgrid:from_gmsh:

``from_gmsh(filename, layers, convention=0, atmosphere_type=2, top_elevation=0, justify = 'r', chars = ascii_lowercase, spaces=True, block_order=None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; from GMSH

Imports a 2-D `Gmsh <http://geuz.org/gmsh/>`_ mesh into a geometry
object. The horizontal structure of the geometry object is created
from the Gmsh mesh, while the layer structure is specified via the
``layers`` parameter, a list of layer thicknesses. The elevation of
the top surface can also be specified, as well as the naming
convention and atmosphere type.

**Parameters:**

-  | **filename**: string
   | Name of the Gmsh mesh file.

-  | **layers**: list
   | List of floats containing the desired layer thicknesses.

-  | **convention**: integer
   | Naming convention for grid columns and layers.

-  | **atmosphere_type**: integer
   | Type of atmosphere.

-  | **top_elevation**: float
   | Elevation of the top surface of the model (default is zero).

-  | **justify**: string
   | Specify 'r' for the character part of the block names (first three
     characters) to be right-justified, 'l' for left-justified.

-  | **chars**: string
   | Specifies a string of characters to use when forming the character
     part of block names. Default is lowercase letters.

-  | **spaces**: Boolean
   | Specify ``False`` to disallow spaces in character part of block
     names. In this case, the first element of the ``chars`` parameter
     functions like a 'zero' and replaces spaces.

-  | **block_order**: string or ``None``
   | Specify ``None`` or 'layer_column' for default block ordering by
     layer and column, starting from the atmosphere. Specify 'dmplex' to
     order blocks by geometrical type (8-node hexahedrons first followed
     by 6-node wedges) as in PETSc DMPlex meshes.

----

.. _sec:mulgrid:from_layermesh:

``from_layermesh(mesh, convention=0, atmosphere_type=2, justify='r',  chars=ascii_lowercase, spaces=True, block_order=None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; from Layermesh

Imports a `Layermesh <https://github.com/acroucher/layermesh>`_ object
into a geometry object.

**Parameters:**

-  | **mesh**: ``layermesh``
   | Layermesh object to import.

-  | **convention**: integer
   | Naming convention for grid columns and layers.

-  | **atmosphere_type**: integer
   | Type of atmosphere.

-  | **justify**: string
   | Specify 'r' for the character part of the block names (first three
     characters) to be right-justified, 'l' for left-justified.

-  | **chars**: string
   | Specifies a string of characters to use when forming the character
     part of block names. Default is lowercase letters.

-  | **spaces**: Boolean
   | Specify ``False`` to disallow spaces in character part of block
     names. In this case, the first element of the ``chars`` parameter
     functions like a 'zero' and replaces spaces.

-  | **block_order**: string or ``None``
   | Specify ``None`` or 'layer_column' for default block ordering by
     layer and column, starting from the atmosphere. Specify 'dmplex' to
     order blocks by geometrical type (8-node hexahedrons first followed
     by 6-node wedges) as in PETSc DMPlex meshes.

----

.. _sec:mulgrid:layer_containing_elevation:

``layer_containing_elevation(elevation)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; searching

Returns the grid layer containing the specified vertical elevation.

**Parameters:**

-  | **elevation**: float
   | Vertical elevation.

----

.. _sec:mulgrid:layer_mapping:

``layer_mapping(geo)``
^^^^^^^^^^^^^^^^^^^^^^

Returns a dictionary mapping each layer name in the ``mulgrid`` object
``geo`` to the name of the nearest layer in the object's own geometry.
(Note: this mapping takes no account of the grid surface, which may
alter which layer is nearest in a given column.)

**Parameters:**

-  | **geo**: :ref:`mulgrid <mulgrids>`
   | The ``mulgrid`` object to create a layer mapping from.

----

.. _sec:mulgrid:layer_name:

``layer_name(block_name)``
^^^^^^^^^^^^^^^^^^^^^^^^^^

Gives the name of the layer corresponding to the specified block name,
according to the naming convention of the grid.

**Parameters:**

-  | **block_name**: string
   | Block name.

----

.. _sec:mulgrid:layer_plot:

``layer_plot(layer, variable=None, variable_name=None, unit=None, column_names=None, node_names=None, column_centres=None, nodes=None, colourmap=None, linewidth=0.2, linecolour='black', aspect='equal', plt=None, subplot=111, title=None, xlabel='x (m)', ylabel='y (m)', contours=False, contour_label_format='%3.0f', contour_grid_divisions=(100,100), connections=None, colourbar_limits=None, plot_limits=None, wells=None, well_names=True, hide_wells_outside=True, wellcolour='blue', welllinewidth=1.0, wellname_bottom=True, rocktypes=None, allrocks=False, rockgroup=None, flow=None, grid=None, flux_matrix=None, flow_variable_name=None, flow_unit=None, flow_scale=None, flow_scale_pos=(0.5, 0.02), flow_arrow_width=None, connection_flows=False, blockmap = {}, block_names=None``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; plotting

Plots a variable over a layer of the grid, using the ``matplotlib``
plotting library. The required layer can be specified by name or as an
elevation (in which case the routine will find the corresponding layer).
Specifying the layer as ``None`` gives a plot over the ground surface of
the geometry (i.e. the surface layer for each column).

The variable can be a list or ``np.array`` containing a value for every
block (or column) in the grid, in the order given by the
``block_name_list`` property of the geometry. If no variable is
specified, only the grid in the layer is plotted, without shading. If
the variable contains a value for each column in the grid, these values
are extended down each column to fill the entire grid.

The name and units of the variable can optionally be specified, and the
names of the columns and nodes can also optionally be displayed on the
plot, as well as the column centres (represented by crosses). The colour
map and limits of the variable shading, the line width of the grid
columns and the aspect ratio of the plot can also be set, as can the
title and x- and y-axis labels, and the plot limits.

When a variable is plotted over the grid, contours at specified levels
can also be drawn, and optionally labelled with their values.

Well tracks can also optionally be plotted. Each well is drawn as a line
following the well track, with the well name at the bottom (or
optionally the top) of the well. For surface plots (``layer`` =
``None``), wells are drawn with solid lines; otherwise, wells are drawn
with dotted lines except where they pass through the specified layer,
where they are drawn with solid lines.

Rock types can be shown on the layer plot by specifying a
:ref:`t2grid <t2grids>` object as the ``rocktypes`` parameter. It is
possible to group similar rock types (e.g. those in the same geological
formation but with slightly different permeabilities) to simplify the
plot if there are a lot of rock types.

Flows can be shown on the layer by specifying an array of connection
flow values (e.g mass flow) as the ``flow`` parameter. Flows will then
be drawn on the slice by arrows at the block centres, each representing
the average flux (flow per unit area) over the block, projected onto the
layer. (For example, connection values of mass flow in kg/s will be
represented as block-average mass fluxes in kg/:math:`m^2`/s.)
Alternatively, flows through the connection faces can be plotted by
setting the ``connection_flows`` parameter to ``True``.

**Parameters:**

-  | **layer**: :ref:`layer <layerobjects>`, string, integer,
     float or ``None``
   | Layer or name (string) of layer to plot, or elevation (float or
     integer). Specifying ``None`` gives a surface plot.

-  | **variable**: list (or ``np.array``)
   | Variable to be plotted, of length equal to the number of blocks or
     columns in the grid (or ``None`` just to plot the grid).

-  | **variable_name**: string
   | Name of the variable (as it will appear on the scale of the plot).

-  | **unit**: string
   | Units of the variable (as it will appear on the scale of the plot).

-  | **column_names**: Boolean or list
   | Set to ``True`` if column names are to be indicated on the plot, or
     to a list of names of columns to be named.

-  | **node_names**: Boolean or list
   | Set to ``True`` if node names are to be indicated on the plot, or
     to a list of names of nodes to be named.

-  | **column_centres**: Boolean or list
   | Set to ``True`` if column centres are to be indicated on the plot
     (as crosses), or to a list of names of columns to be indicated.

-  | **nodes**: Boolean or list
   | Set to ``True`` if nodes are to be indicated on the plot (as
     crosses), or to a list of names of nodes to be indicated.

-  | **colourmap**: string
   | Name of ``matplotlib`` colour map to use for shading the variable.

-  | **linewidth**: float
   | Line width to use for drawing the grid.

-  | **linecolour**: string
   | Line colour to use for drawing the grid.

-  | **aspect**: string
   | Aspect ratio to use for drawing the grid (default is 'equal' (i.e.
     1:1).

-  | **plt**: ``matplotlib.pyplot`` instance
   | An instance of the ``matplotlib.pyplot`` library, imported in the
     calling script using e.g. ``import matplotlib.pyplot as plt``.

-  | **subplot**: integer
   | Subplot number for multi-plots, e.g. set 223 to draw the third plot
     in a 2-by-2 multiplot (default is 111).

-  | **title**: string
   | Plot title. If set to ``None`` (the default value), a title will be
     constructed from the other plot parameters. Set to for no title.

-  | **xlabel**: string
   | x axis label (default is 'x (m)').

-  | **ylabel**: string
   | y axis label (default is 'y (m)').

-  | **contours**: Boolean, list or ``np.array``
   | Set to ``True`` or to a list or array of contour values to draw
     contours on the plot (default ``False``).

-  | **contour_label_format**: string
   | Format string for contour labels (default '%3.0f').

-  | **contour_grid_divisions**: tuple (of integer)
   | Number of divisions in the x- and y-directions in the regular grid
     superimposed on the model grid, and used to produce the contours
     (default (100,100)).

-  | **connections**: float (or ``None``)
   | Set non-zero to plot connections in the grid, shaded by absolute
     value of the connection angle cosine. The value specifies the lower
     cut-off value, above which connections will be plotted. Connections
     are shaded in greyscale from white (0.0) to black (1.0). This can
     be used to check orthogonality of grid connections, as less
     orthogonal connections (with larger angle cosine) will show up
     darker on the plot. If set to ``None``, no connections will be
     plotted.

-  | **colourbar_limits**: tuple, list, ``np.array`` (or ``None``)
   | Specify a two-element tuple, list or ``np.array`` to set the limits
     of the colour scale. Default (``None``) will auto-scale.

-  | **plot_limits**: tuple or list (or ``None``)
   | Specify a two-element tuple (or list) of plot axis ranges, each
     itself being a tuple (or list) of minimum and maximum values, i.e.
     ((xmin,xmax),(ymin,ymax)). Default is ``False`` which will
     auto-scale.

-  | **wells**: Boolean or list (or ``None``)
   | Specify ``True`` to plot all well tracks, ``False`` or ``None`` not
     to plot them, or a list of wells or well names to specify only
     particular wells.

-  | **well_names**: Boolean or list (or ``None``)
   | Specify ``True`` to label each well with its name , ``False`` or
     ``None`` not to label them, or a list of wells or well names to
     label only particular wells.

-  | **hide_wells_outside**: Boolean
   | Set to ``True`` if wells that do not intersect the specified layer
     are to be hidden.

-  | **wellcolour**: string
   | Colour to use for drawing the wells.

-  | **welllinewidth**: float
   | Line width for drawing the wells.

-  | **wellname_bottom**: Boolean
   | Set to ``False`` to label wells at the wellhead rather than the
     bottom.

-  | **rocktypes**: :ref:`t2grid <t2grids>` (or ``None``)
   | To plot rock types, specify a ``t2grid`` object containing rock
     types for the grid. If ``None``, no rock types will be plotted.

-  | **allrocks**: Boolean
   | If ``False`` (the default), only rock types present on the
     specified layer will be shown in the colour bar; others will be
     omitted. If ``True``, all rocks present in the model grid will be
     shown on the colour bar, regardless of whether they appear in the
     specified layer.

-  | **rockgroup**: tuple, list, string (or ``None``)
   | To group similar rock types into one colour, specify a tuple or
     list of integers, representing the significant characters of the
     rock type names. For example, to group rock types having the same
     first two characters, specify (0,1). Alternatively, specify a
     5-character string mask containing asterisks in positions that are
     not significant, and any other characters in the significant
     positions (e.g. '++**\*').

-  | **flow**: ``np.array`` (or ``None``)
   | To plot flows, specify an array of connection flow values (one
     floating point value for each connection in the grid). These may
     for example be extracted from the columns of the connection table
     in a :ref:`t2listing <listingfiles>` object.

-  | **grid**: :ref:`t2grid <t2grids>` (or ``None``)
   | Specify a ``t2grid`` object associated with the grid, to be used to
     calculate the 'flux matrix' which converts the connection flow
     values to block-average fluxes. If this is not specified (and
     neither is the ``flux_matrix`` parameter), then a ``t2grid`` object
     will be created internally.

-  | **flux_matrix**: ``scipy.sparse.lil_matrix`` (or ``None``)
   | A sparse matrix used to convert the connection flow values to
     block-average fluxes. Such a matrix can be created using the
     :ref:`flux_matrix() <sec:t2grid:flux_matrix>` method of
     a ``t2grid`` object and an appropriate ``mulgrid`` object. If no
     flux matrix is specified, one will be created internally. This can
     be time-consuming for large grids, so for multiple flow plots it is
     faster to pre-calculate a flux matrix in your script and pass it
     via this parameter. If this parameter is specified, there is no
     need also to specify the ``grid`` parameter.

-  | **flow_variable_name**: string (or ``None``)
   | Name of the flow variable (as it will appear on the scale of the
     plot).

-  | **flow_unit**: string (or ``None``)
   | Units of the flow variable (as it will appear on the scale of the
     plot, divided by area).

-  | **flow_scale**: string (or ``None``)
   | Length of flow scale arrow. If not specified, this will be
     calculated.

-  | **flow_scale_pos**: tuple
   | Position of the flow scale on the plot, in units of dimensionless
     plot size. The default (0.5, 0.02) draws the flow scale in the
     horizontal centre of the plot, slightly above the bottom axis. If
     you want the flow scale below the bottom axis (so it doesn't get
     mixed up with the actual flow arrows), specify this parameter with
     a small negative second component, e.g. (0.8, -0.1).

-  | **flow_arrow_width**: float (or ``None``)
   | Width of the flow arrows, in units of dimensionless plot width. If
     not specified, this will be calculated internally.

-  | **connection_flows**: Boolean
   | Set to ``True`` to plot flows through connection faces, rather than
     block-averaged fluxes. In this case, usually the ``grid`` parameter
     should also be specified (but not ``flux_matrix``), otherwise a
     grid will be calculated internally.

-  | **blockmap**: dictionary
   | Dictionary mapping the block names in the geometry to another block
     naming system. This has an effect only on the block names displayed
     on the plot via the ``block_names`` parameter, and on the rock
     types displayed. Note that if a mapping is used, then the
     ``block_names`` list should contain mapped block names.

-  | **block_names**: Boolean or list
   | Set to ``True`` if block names are to be indicated on the plot, or
     to a list of names of blocks to be named.

**Example:**

::

   geo.layer_plot(-500., t, 'Temperature', '$\degree$C', contours = np.arange(100,200,25))

plots the variable ``t`` at elevation -500 m over the grid, with the
values as Temperature (°C), and with contours drawn from 100°C to
200°C with a contour interval of 25°C.

----

.. _sec:mulgrid:line_plot:

``line_plot(start=None, end=None, variable, variable_name=None, unit=None, divisions=100, plt=None, subplot=111, title='', xlabel='distance (m)', coordinate=False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; plotting

Plots a variable along a line through the grid, using the ``matplotlib``
plotting library. The line is specified by its start and end points in
3D. The variable can be a list or ``np.array`` containing a value for
every block (or column) in the grid. If the variable contains a value
for each column in the grid, these values are extended down each column
to fill the entire grid. The name and units of the variable can
optionally be specified, as well as the number of divisions the line is
divided into (default 100), the plot title and the axis labels.

**Parameters:**

-  | **start**, **end**: list, tuple or ``np.array``
   | Start and end point of the line, each of length 3 (``None`` to plot
     across the bounds of the grid).

-  | **variable**: list (or ``np.array``)
   | Variable to be plotted, of length equal to the number of blocks (or
     columns) in the grid.

-  | **variable_name**: string
   | Name of the variable (as it will appear on the scale of the plot).

-  | **unit**: string
   | Units of the variable (as it will appear on the scale of the plot).

-  | **divisions**: integer
   | Number of divisions to divide the line into (default 100).

-  | **plt**: ``matplotlib.pyplot`` instance
   | An instance of the ``matplotlib.pyplot`` library, imported in the
     calling script using e.g. ``import matplotlib.pyplot as plt``.

-  | **subplot**: integer
   | Subplot number for multi-plots, e.g. set 223 to draw the third plot
     in a 2-by-2 multiplot (default is 111).

-  | **title**: string
   | Plot title. If set to ``None`` (the default value), a title will be
     constructed from the other plot parameters. Set to for no title.

-  | **xlabel**: string
   | x axis label (default is 'distance (m)').

-  | **coordinate**: integer or Boolean
   | If ``False``, plot against distance along the line, otherwise plot
     against specified coordinate (0,1 or 2) values.

**Example:**

::

   geo.line_plot([0.,0.,500.], [1000.,0.,500.], t, 'Temperature', '$\degree$C')

plots the variable ``t`` along a line from (0,0,500) to (1000,0,500)
through the grid, with the values as Temperature (°C).

----

.. _sec:mulgrid:line_values:

``line_values(start, end, variable, divisions=100, coordinate=False, qtree=None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns values of a specified variable along an arbitrary line through
the grid. The start and end points of the line (``start`` and ``end``)
are 3-element lists, tuples or ``np.arrays`` specifying points in 3D.
The variable can be a list or ``np.array`` containing a value for every
block in the grid. The number of divisions along the line (default 100)
can be optionally specified.

The routine returns a tuple of two arrays (*l*,\ *v*), the first (*l*)
containing the distance from the start (or the appropriate coordinate
(0,1, or 2) if ``coordinate`` is specified) for each point along the
line, and the second (*v*) containing the value of the variable at that
point. The value of the variable at any point is the (block average)
value at the block containing the point.

**Parameters:**

-  | **start**, **end**: list, tuple or ``np.array`` (of length 3)
   | Start and end points of the line in 3D.

-  | **variable**: list (or ``np.array``)
   | Variable to be plotted, of length equal to the number of blocks in
     the grid.

-  | **divisions**: integer
   | Number of segments the line is divided up into (default 100).

-  | **coordinate**: integer or Boolean
   | If ``False``, return distance along the line in first array,
     otherwise return specified coordinate (0,1 or 2) values.

-  | **qtree**: ``quadtree``
   | Quadtree object for fast searching of grid columns (can be
     constructed using the ``column_quadtree()`` method).

----

.. _sec:mulgrid:meshio_grid:

``meshio_grid(surface_snap = 0.1, dimension = 3, slice = None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; meshio grid

Returns mesh corresponding to the geometry, in the format used by the
``meshio`` library (https://pypi.python.org/pypi/meshio). This consists
of a two-element tuple: firstly, an ``np.array`` of nodal coordinates,
and secondly a dictionary of element definitions, indexed by number of
nodes in the elements.

The primary use of this is as an interchange format for input/output of
meshes in different formats. Note that exporting the geometry directly
to a mesh file can also be done using the
:ref:`write_mesh() <sec:mulgrid:write_mesh>` method (which is
just a wrapper for this one).

**Parameters:**

-  | **surface_snap**: float
   | Tolerance for eliminating elements with very small vertical
     thickness at the top of the mesh.

-  | **dimension**: integer
   | Dimension of the mesh: when set to 3, return the full 3-D mesh.
     When set to 2, return a 2-D mesh, corresponding either to the
     horizontal mesh only (the default), or a vertical slice mesh if the
     ``slice`` parameter is used.

-  | **slice**: list, string, float or ``None``
   | Horizontal line defining the slice for vertical 2-D meshes. This
     can be a list of two horizontal (*x*,\ *y*) points (``np.arrays``)
     defining the endpoints of the slice line, or string 'x' or 'y' to
     specify the *x*- or *y*-axis, or northing (float) through grid
     centre. If set to ``None`` (the default) then the horizontal 2-D
     mesh is returned.

----

.. _sec:mulgrid:minc_array:

``minc_array(vals, minc_indices, level=0, outside=0.0)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; MINC arrays

Returns an array for all blocks in the geometry, with values taken from
the input ``vals`` array, for the specified MINC level. Indexing of MINC
blocks is specified by the ``minc_indices`` array (returned by the
``t2grid`` :ref:`minc() <sec:t2grid:MINC>` method).

**Parameters:**

-  | **vals**: ``np.array``
   | Array of values over the entire MINC grid, with values for all MINC
     levels, obtained e.g. from a column of the element table of a
     :ref:`t2listing <listingfiles>` object.

-  | **minc_indices**: ``np.array`` (of integer)
   | Rank-2 array containing integer indices for each MINC level,
     obtained from the output of the ``t2grid``
     :ref:`minc() <sec:t2grid:MINC>` method.

-  | **level**: integer
   | MINC level, 0 being the fracture level and higher levels being the
     matrix levels.

-  | **outside**: Boolean or float
   | Determines how blocks outside the MINC part of the grid are
     handled. If ``True``, include porous medium values outside the MINC
     part of the grid. If a float value is given, assign that value
     instead. If ``False``, the value zero will be assigned.

----

.. _sec:mulgrid:nodes_in_columns:

``nodes_in_columns(columns)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; searching

Returns a list of all nodes in a specified list of columns.

**Parameters:**

-  | **columns**: list (of :ref:`column <columnobjects>`)
   | List of columns in which to find nodes.

----

.. _sec:mulgrid:nodes_in_polygon:

``nodes_in_polygon(polygon)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; searching

Returns a list of all nodes inside the specified polygon or rectangle.

**Parameters:**

-  | **polygon**: list (of ``np.array``)
   | List of points defining the polygon (each point is a two-element
     ``np.array``). If the list has only two points, it will be
     interpreted as a rectangle [bottom left, top right].

----

.. _sec:mulgrid:node_nearest_to:

``node_nearest_to(point, kdtree=None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; searching

Returns the node nearest to a specified point. An optional kd-tree
structure can be specified to speed searching - useful if searching for
many points.

**Parameters:**

-  | **point**: ``np.array``, list or tuple
   | Array or list of length 2, specifying the required point in 2-D.

-  | **kdtree**: ``cKDTree``
   | kd-tree structure for searching for nodes. Such a tree can be
     constructed using the ``node_kdtree`` property of a ``mulgrid``
     object. You will need the ``scipy`` library installed before you
     can use this property.

----

.. _sec:mulgrid:optimize:

``optimize(nodenames=None, connection_angle_weight=1.0, column_aspect_weight=0.0, column_skewness_weight=0.0, pest=False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; optimizing

Adjusts positions of the specified nodes to optimize grid quality. If no
nodes are specified, all node positions are optimized. Grid quality can
be defined as a combination of connection angle cosine, column aspect
ratio and column skewness. Increasing the weight for any of these
increases its importance in the evaluation of grid quality.

Note that an error will result if the connection angle weight and either
of the other two weights is set to zero - in this case there are not
enough constraints to fit the parameters.

If the ``pest`` parameter is set to ``True``, the `PEST
<http://www.pesthomepage.org/>`_ parameter estimation software is used
to carry out the optimzation (this obviously requires that PEST is
installed on your machine). Otherwise, the ``leastsq`` routine in the
``scipy`` Python library is used. PEST seems to be more robust in some
cases, and often gives better results when nodes on the boundary of
the grid are included in the optimization.  However, when ``leastsq``
does work satisfactorily, it is generally faster (mainly because PEST
has to read the geometry from disk and write it out again each time
the grid quality is evaluated during the optimization). If PEST is
used, a variety of intermediate files (named ``pestmesh.*``) will be
written to the working directory, including the PEST run record file
(``pestmesh.rec``) which contains a detailed record of the
optimization process.

**Parameters:**

-  | **nodenames**: list of string
   | List of names of nodes to optimize. If not specified, all nodes in
     the grid are optimized.

-  | **connection_angle_weight**: float
   | Weighting to be given to connection angle cosines. A higher value
     will place greater priority on making connections perpendicular to
     the column sides.

-  | **column_aspect_weight**: float
   | Weighting to be given to column aspect ratios. A higher value will
     place greater priority on making column side ratios closer to 1.0.

-  | **column_skewness_weight**: float
   | Weighting to be given to column skewness. A higher value will place
     greater priority on making column angle ratios closer to 1.0.

-  | **pest**: Boolean
   | Set ``True`` to use the PEST parameter estimation software to
     perform the optimization.

----

.. _sec:mulgrid:polyline_values:

``polyline_values(polyline, variable, divisions=100, coordinate=False, qtree=None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns values of a specified variable along an arbitrary polyline
through the grid, defined as a list of 3-element lists or ``np.arrays``
specifying points in 3D. The variable can be a list or ``np.array``
containing a value for every block in the grid. The number of divisions
along the line (default 100) can be optionally specified.

The routine returns a tuple of two arrays (``l``,\ ``v``), the first
(``l``) containing the distance from the start (or the appropriate
coordinate (0, 1, or 2) if ``coordinate`` is specified) for each point
along the polyline, and the second (``v``) containing the value of the
variable at that point. The value of the variable at any point is the
(block average) value at the block containing the point.

**Parameters:**

-  | **polyline**: list of 3-element lists or ``np.arrays``
   | Polyline points in 3D.

-  | **variable**: list (or ``np.array``)
   | Variable to be plotted, of length equal to the number of blocks in
     the grid.

-  | **divisions**: integer
   | Number of segments the line is divided up into (default 100).

-  | **coordinate**: integer or Boolean
   | If ``False``, return distance along the line in first array,
     otherwise return specified coordinate (0, 1 or 2) values.

-  | **qtree**: ``quadtree``
   | Quadtree object for fast searching of grid columns (can be
     constructed using the :ref:`column_quadtree() <sec:mulgrid:column_quadtree>`
     method).

----

.. _sec:mulgrid:read:

``read(filename)``
^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; reading

Reads a ``mulgrid`` object from a MULgraph geometry file on disk.

**Parameters:**

-  | **filename**: string
   | Name of the MULgraph geometry file to be read.

**Example:**

::

   geo = mulgrid().read(filename)

creates a ``mulgrid`` object and reads its contents from file
``filename``. This can be done more simply just by passing the filename
into the ``mulgrid`` creation command:

::

   geo = mulgrid(filename)

----

.. _sec:mulgrid:rectangular:

``rectangular(xblocks, yblocks, zblocks, convention=0, atmos_type=2, origin=[0,0,0], justify='r', case=None, chars=ascii_lowercase, spaces=True, block_order=None``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; rectangular

Gives a ``mulgrid`` geometry object a rectangular grid structure. The
grid sizes in the *x*, *y* and *z* directions can be non-uniform, and
the grid column and layer naming convention, atmosphere type and origin
can be specified. The optional ``justify`` and ``case`` parameters
control the formatting of the character part of the block names.
Additionally, the characters used to form node/column or layer names can
be specified using the ``chars`` parameter. (This can be useful for
example for grids with large numbers of nodes and/or columns, for which
lowercase letters alone may not be enough.)

Note that it is also possible to reverse-engineer a rectangular geometry
from an existing TOUGH2 data file or ``t2grid`` object, using the
:ref:`rectgeo() <sec:t2grid:rectgeo>` method.

**Parameters:**

-  | **xblocks**, **yblocks**, **zblocks**: list, tuple or ``np.array``
   | Lists (or arrays) of block sizes (float) in the *x*, *y* and *z*
     directions.

-  | **convention**: integer
   | Naming convention for grid columns and layers.

-  | **atmos_type**: integer
   | Type of atmosphere.

-  | **origin**: list (or ``np.array``)
   | Origin of the grid (of length 3).

-  | **justify**: string
   | Specify 'r' for the character part of the block names (first three
     characters) to be right-justified, 'l' for left-justified.

-  | **case**: string
   | Specify 'l' for the character part of the block names (first three
     characters) to be lower case, 'u' for upper case. Now deprecated -
     using the ``chars`` parameter is more flexible.

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

-  | **block_order**: string or ``None``
   | Specify ``None`` or 'layer_column' for default block ordering by
     layer and column, starting from the atmosphere. Specify 'dmplex' to
     order blocks by geometrical type (8-node hexahedrons first followed
     by 6-node wedges) as in PETSc DMPlex meshes.

**Example:**

::

   geo = mulgrid().rectangular([1000]*10, [500]*20, [100]*5+[200]*10, origin=[0,0,2500])

creates a ``mulgrid`` object called ``geo``, and fills it with a
rectangular grid of 10 blocks of size 1000 m in the *x*-direction, 20
blocks of size 500 m in the *x*-direction, 5 layers at the top of
thickness 100 m and 10 layers underneath of thickness 200 m, and with
origin (0,0,2500) m. The grid will have the default naming convention
(0) and atmosphere type (2).

----

.. _sec:mulgrid:reduce:

``reduce(columns)``
^^^^^^^^^^^^^^^^^^^

Reduces a grid so that it contains only the specified list of columns
(or columns with specified names).

**Parameters:**

-  | **columns**: list
   | List of required columns or column names.

----

.. _sec:mulgrid:refine:

``refine(columns=[], bisect=False, bisect_edge_columns=[], chars = ascii_lowercase, spaces=True)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; refining

Refines the specified columns in the grid. Appropriate transition
columns are created around the refined region. If no columns are
specified, all columns are refined. All columns in the region to be
refined (and in the transition region) must be either triangular or
quadrilateral. Each column in split into four, unless the ``bisect``
parameter is ``True``, in which case each column in split into two. If
``bisect`` is 'x' or 'y', columns are split in the closest direction to
the axis specified; or if ``bisect`` is ``True``, between its longest
sides.

The ``bisect_edge_columns`` parameter can be used to give more desirable
column shapes in the transition region, if the original columns
occupying the transition region have large aspect ratios. By default,
these will become even worse when they are triangulated to form the
transition columns, if they are connected to the refinement region by
their shorter sides. Including them in ``bisect_edge_columns`` means
they will be bisected (parallel to the edge of the refinement region)
before the refinement is carried out, which should improve the aspect
ratios of the transition columns.

**Note**: TOUGH2 implicitly assumes that the connections in its finite
volume grids are orthogonal, i.e. the line joining the centres of two
connected blocks should be perpendicular to the connecting face. The
triangular transition columns generated by the ``refine()`` method
generally give rise to connections that are not orthogonal. However,
they can be modified and made as orthogonal as possible using the
:ref:`optimize() <sec:mulgrid:optimize>` method.

**Parameters:**

-  | **columns**: list
   | List of columns or column names to be refined.

-  | **bisect**: Boolean or string
   | Set to ``True`` if columns are to be split into two, between their
     longest sides, instead of four (the default). Set to 'x' or 'y' to
     split columns along the specified axis.

-  | **bisect_edge_columns**: list
   | List of columns or column names in the transition region (just
     outside the refinement area) to be bisected prior to the
     refinement, to improve the aspect ratios of the transition columns.

-  | **chars**: string
   | Specifies a string of characters to use when forming the character
     part of block names. Default is lowercase letters.

-  | **spaces**: Boolean
   | Specify ``False`` to disallow spaces in character part of block
     names. In this case, the first element of the ``chars`` parameter
     functions like a 'zero' and replaces spaces.

----

.. _sec:mulgrid:refine_layers:

``refine_layers(layers=[], factor=2, chars = ascii_lowercase, spaces=True)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; refining

Refines the specified layers in the grid. If no layers are specified,
all layers are refined. Each layer is refined by the specified integer
factor.

**Parameters:**

-  | **layers**: list
   | List of layers or layer names to be refined.

-  | **factor**: integer
   | Refinement factor: default is 2, which bisects each layer.

-  | **chars**: string
   | Specifies a string of characters to use when forming the character
     part of block names. Default is lowercase letters.

-  | **spaces**: Boolean
   | Specify ``False`` to disallow spaces in character part of block
     names. In this case, the first element of the ``chars`` parameter
     functions like a 'zero' and replaces spaces.

----

.. _sec:mulgrid:rename_column:

``rename_column(oldcolname, newcolname)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Renames a grid column. Returns ``True`` if the column was found and
renamed, or ``False`` if the specified column does not exist. Multiple
columns can be renamed at once by specifying lists of old and new column
names - this is faster than calling the method multiple times, and the
block and connection name lists are updated only once.

**Parameters:**

-  | **oldcolname**: string or list of strings
   | Name(s) of the column(s) to rename.

-  | **newcolname**: string or list of strings
   | New name(s) of the column(s).

----

.. _sec:mulgrid:rename_layer:

``rename_layer(oldlayername, newlayername)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Renames a grid layer. Returns ``True`` if the layer was found and
renamed, or ``False`` if the specified layer does not exist. Multiple
layers can be renamed at once by specifying lists of old and new layer
names - this is faster than calling the method multiple times, and the
block and connection name lists are updated only once.

**Parameters:**

-  | **oldlayername**: string or list of strings
   | Name(s) of the layer(s) to rename.

-  | **newlayername**: string or list of strings
   | New name(s) of the layer(s).

----

.. _sec:mulgrid:rotate:

``rotate(angle, centre=None, wells=False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; rotating

Rotates a grid by a specified angle (in degrees) clockwise in the
horizontal plane. Any wells in the grid are also rotated. The centre of
rotation can be optionally specified. If it is not specified, the centre
of the grid is used as the centre of rotation. If the ``wells``
parameter is ``True``, any wells in the grid are also rotated.

**Parameters:**

-  | **angle**: float
   | Angle (in degrees) to rotate the grid, positive for clockwise,
     negative for anti-clockwise.

-  | **centre**: list, tuple or ``np.array``
   | Centre of rotation in the horizontal *x*,\ *y* plane (of length 2).

-  | **wells**: Boolean
   | Set ``True`` to rotate wells.

**Example:**

::

   geo.rotate(30)

rotates the grid ``geo`` clockwise by 30° about its centre in the
horizontal plane.

----

.. _sec:mulgrid:slice_plot:

``slice_plot(line=None, variable=None, variable_name=None, unit=None, block_names=None, colourmap=None, linewidth=0.2, linecolour='black', aspect='auto', plt=None, subplot=111, title=None, xlabel='', ylabel='elevation (m)', contours=False, contour_label_format='%3.0f', contour_grid_divisions=(100,100), colourbar_limits=None, plot_limits=None, column_axis=False, layer_axis=False, wells=None, well_names=True, hide_wells_outside=False, wellcolour='blue', welllinewidth=1.0, wellname_bottom=False, rocktypes=None, allrocks=False, rockgroup=None, flow=None, grid=None, flux_matrix=None, flow_variable_name=None, flow_unit=None, flow_scale=None, flow_scale_pos=(0.5, 0.02), flow_arrow_width=None, connection_flows=False, blockmap = {})``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; plotting

Plots a variable over a vertical slice through the grid, using the
``matplotlib`` plotting library. The required slice is specified by a
horizontal line through the grid, defined as either a two-element list
of (*x*,\ *y*) points (``np.arrays``), or as a string 'x' or 'y' which
defines the *x*- or *y*-axes respectively, or as a northing (in degrees)
through the centre of the grid. If no line is specified, the line is
taken to be across the bounds of the grid. For slice plots along the x-
or y-axis, the horizontal coordinate represents the x- or y-coordinate;
for other slice directions it represents distance along the slice line.

The variable can be a list or ``np.array`` containing a value for every
block (or column) in the grid, in the order given by the
``block_name_list`` property of the geometry. If no variable is
specified, only the grid is plotted, without shading. If the variable
contains a value for each column in the grid, these values are extended
down each column to fill the entire grid.

The name and units of the variable can optionally be specified, and the
name of each block can also optionally be displayed on the plot. The
colour map and limits of the variable shading, the line width of the
grid columns and the aspect ratio of the plot can also be set, as can
the plot title and x- and z-axis labels, and the plot limits.

When a variable is plotted over the grid, contours at specified levels
can also be drawn, and optionally labelled with their values.

Well tracks can also optionally be plotted. Each well is drawn as a line
following the well track, with the well name at the top (or optionally
the bottom) of the well. If ``hide_wells_outside`` is specified as a
floating point number, wells that do not pass within the specified
distance from the slice line are not shown. Well tracks are shown as
solid lines over sections within the specified distance from the slice
line, and dotted lines otherwise.

Rock types can be shown on the slice plot by specifying a ``t2grid``
object as the ``rocktypes`` parameter. It is possible to group similar
rock types (e.g. those in the same geological formation but with
slightly different permeabilities) to simplify the plot if there are a
lot of rock types.

Flows can be shown on the slice by specifying an array of connection
flow values (e.g mass flow) as the ``flow`` parameter. Flows will then
be drawn on the slice by arrows at the block centres, each representing
the average flux (flow per unit area) over the block, projected onto the
slice. (For example, connection values of mass flow in kg/s will be
represented as block-average mass fluxes in kg/:math:`m^2`/s.)
Alternatively, flows through the connection faces can be plotted by
setting the ``connection_flows`` parameter to ``True``.

**Parameters:**

-  | **line**: list, string or float
   | List of two horizontal (*x*,\ *y*) points (``np.arrays``) defining
     the endpoints of the line, or string 'x' or 'y' to specify the *x*-
     or *y*-axis, or northing (float) through grid centre.

-  | **variable**: list (or ``np.array``)
   | Variable to be plotted, of length equal to the number of blocks (or
     columns) in the grid (or ``None`` just to plot the grid).

-  | **variable_name**: string
   | Name of the variable (as it will appear on the scale of the plot).

-  | **unit**: string
   | Units of the variable (as it will appear on the scale of the plot).

-  | **block_names**: Boolean or list
   | Set to ``True`` if block names are to be indicated on the plot, or
     to a list of names of blocks to be named.

-  | **colourmap**: string
   | Name of ``matplotlib`` colour map to use for shading the variable.

-  | **linewidth**: float
   | Line width to use for drawing the grid.

-  | **linecolour**: string
   | Line colour to use for drawing the grid.

-  | **aspect**: string
   | Aspect ratio to use for drawing the grid (default is 'auto').

-  | **plt**: ``matplotlib.pyplot`` instance
   | An instance of the ``matplotlib.pyplot`` library, imported in the
     calling script using e.g. ``import matplotlib.pyplot as plt``.

-  | **subplot**: integer
   | Subplot number for multi-plots, e.g. set 223 to draw the third plot
     in a 2-by-2 multiplot (default is 111).

-  | **title**: string
   | Plot title. If set to ``None`` (the default value), a title will be
     constructed from the other plot parameters. Set to for no title.

-  | **xlabel**: string
   | x axis label. If set to ``None`` (the default value), a label will
     be constructed according to the slice orientation- either 'x (m)',
     'y (m)' or 'distance (m)' as appropriate.

-  | **ylabel**: string
   | y axis label (default is 'elevation (m)').

-  | **contours**: Boolean, list or ``np.array``
   | Set to ``True`` or to a list or array of contour values to draw
     contours on the plot (default ``False``).

-  | **contour_label_format**: string
   | Format string for contour labels (default '%3.0f').

-  | **contour_grid_divisions**: tuple (of integer)
   | Number of divisions in the x- and z-directions in the regular grid
     superimposed on the slice, and used to produce the contours
     (default (100,100)).

-  | **colourbar_limits**: tuple, list, ``np.array`` (or ``None``)
   | Specify a two-element tuple, list or ``np.array`` to set the limits
     of the colour scale. Default (``None``) will auto-scale.

-  | **plot_limits**: tuple or list (or ``None``)
   | Specify a two-element tuple (or list) of plot axis ranges, each
     itself being a tuple (or list) of minimum and maximum values, i.e.
     ((xmin,xmax),(zmin,zmax)). Default is ``False`` which will
     auto-scale.

-  | **column_axis**: Boolean
   | If ``True``, show column names instead of coordinates on the
     horizontal axis.

-  | **layer_axis**: Boolean
   | If ``True``, show layer names instead of coordinates on the
     vertical axis.

-  | **wells**: Boolean or list (or ``None``)
   | Specify ``True`` to plot all well tracks, ``False`` or ``None`` not
     to plot them, or a list of wells or well names to specify only
     particular wells.

-  | **well_names**: Boolean or list (or ``None``)
   | Specify ``True`` to label each well with its name , ``False`` or
     ``None`` not to label them, or a list of wells or well names to
     label only particular wells.

-  | **hide_wells_outside**: ``False`` or float
   | Specify distance as a floating point number to hide wells further
     from the slice line than the specified distance.

-  | **wellcolour**: string
   | Colour to use for drawing the wells.

-  | **welllinewidth**: float
   | Line width for drawing the wells.

-  | **wellname_bottom**: Boolean
   | Set to ``True`` to label wells at the bottom rather than the
     wellhead.

-  | **rocktypes**: :ref:`t2grid <t2grids>` (or ``None``)
   | To plot rock types, specify a ``t2grid`` object containing rock
     types for the grid. If ``None``, no rock types will be plotted.

-  | **allrocks**: Boolean
   | If ``False`` (the default), only rock types present on the
     specified slice will be shown in the colour bar; others will be
     omitted. If ``True``, all rocks present in the model grid will be
     shown on the colour bar, regardless of whether they appear in the
     specified slice.

-  | **rockgroup**: tuple, list, string (or ``None``)
   | To group similar rock types into one colour, specify a tuple or
     list of integers, representing the significant characters of the
     rock type names. For example, to group rock types having the same
     first two characters, specify (0,1). Alternatively, specify a
     5-character string mask containing asterisks in positions that are
     not significant, and any other characters in the significant
     positions (e.g. '++**\*').

-  | **flow**: ``np.array`` (or ``None``)
   | To plot flows, specify an array of connection flow values (one
     floating point value for each connection in the grid). These may
     for example be extracted from the columns of the connection table
     in a ``t2listing`` object.

-  | **grid**: :ref:`t2grid <t2grids>` (or ``None``)
   | Specify a ``t2grid`` object associated with the grid, to be used to
     calculate the 'flux matrix' which converts the connection flow
     values to block-average fluxes. If this is not specified (and
     neither is the ``flux_matrix`` parameter), then a ``t2grid`` object
     will be created internally.

-  | **flux_matrix**: ``scipy.sparse.lil_matrix`` (or ``None``)
   | A sparse matrix used to convert the connection flow values to
     block-average fluxes. Such a matrix can be created using the
     ``flux_matrix()`` method of a ``t2grid`` object and an appropriate
     ``mulgrid`` object. If no flux matrix is specified, one will be
     created internally. This can be time-consuming for large grids, so
     for multiple flow plots it is faster to pre-calculate a flux matrix
     in your script and pass it via this parameter. If this parameter is
     specified, there is no need also to specify the ``grid`` parameter.

-  | **flow_variable_name**: string (or ``None``)
   | Name of the flow variable (as it will appear on the scale of the
     plot).

-  | **flow_unit**: string (or ``None``)
   | Units of the flow variable (as it will appear on the scale of the
     plot, divided by area).

-  | **flow_scale**: string (or ``None``)
   | Length of flow scale arrow. If not specified, this will be
     calculated.

-  | **flow_scale_pos**: tuple
   | Position of the flow scale on the plot, in units of dimensionless
     plot size. The default (0.5, 0.02) draws the flow scale in the
     horizontal centre of the plot, slightly above the bottom axis. If
     you want the flow scale below the bottom axis (so it doesn't get
     mixed up with the actual flow arrows), specify this parameter with
     a small negative second component, e.g. (0.8, -0.1).

-  | **flow_arrow_width**: float (or ``None``)
   | Width of the flow arrows, in units of dimensionless plot width. If
     not specified, this will be calculated internally.

-  | **connection_flows**: Boolean
   | Set to ``True`` to plot flows through connection faces, rather than
     block-averaged fluxes. In this case, usually the ``grid`` parameter
     should also be specified (but not ``flux_matrix``), otherwise a
     grid will be calculated internally.

-  | **blockmap**: dictionary
   | Dictionary mapping the block names in the geometry to another block
     naming system. This has an effect only on the block names displayed
     on the plot via the ``block_names`` parameter, and on the rock
     types displayed. Note that if a mapping is used, then the
     ``block_names`` list should contain mapped block names.

**Example:**

::

   geo.slice_plot(45., t, 'Temperature', '$\degree$C', contours = [100,200])

plots the variable ``t`` through a SW–NE vertical slice (heading 45°)
through the grid, with the values as Temperature (°C) and contours
drawn at 100°C and 200°C.

::

   from matplotlib import cm
   cmap = cm.get_cmap('jet', 10)
   geo.slice_plot(45., t, 'Temperature', '$\degree$C',
     colourbar_limits = (0., 250.), colourmap = cmap)

plots the variable ``t`` again, but with a specified discrete colour
scale with 10 divisions from zero to 250°C.

----

.. _sec:mulgrid:snap_columns_to_layers:

``snap_columns_to_layers(min_thickness=1.0, columns=[])``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; snapping

Snaps column surfaces to the bottom of their layers, if the surface
block thickness is smaller than a given value. This can be carried out
over an optional subset of columns in the grid, otherwise over all
columns.

**Parameters:**

-  | **min_thickness**: float
   | Minimum surface block thickness. Blocks with thickness less than
     this value will be eliminated by 'snapping' the column surface
     elevation to the bottom of the surface layer. Values of
     ``min_thickness`` less than or equal to zero will have no effect.

-  | **columns**: list (of :ref:`column <columnobjects>` or
     string)
   | List of columns to process. If empty (the default), process all
     columns.

----

.. _sec:mulgrid:snap_columns_to_nearest_layers:

``snap_columns_to_nearest_layers(columns=[])``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; snapping

Snaps column surfaces to the nearest layer elevation (top or bottom).
This can be carried out over an optional subset of columns in the grid,
otherwise over all columns.

**Parameters:**

-  | **columns**: list (of :ref:`column <columnobjects>` or
     string)
   | List of columns to process. If empty (the default), process all
     columns.

----

.. _sec:mulgrid:split_column:

``split_column(colname, nodename, chars = ascii_lowercase)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Splits a quadrilateral column with specified name into two triangular
columns. The direction of the split is determined by specifying the name
of one of the splitting nodes. The method returns ``True`` if the split
was carried out successfully.

**Parameters:**

-  | **colname**: string
   | Name of the quadrilateral column to be split. If the column is not
     quadrilateral, the method returns ``False`` and nothing is done to
     the column.

-  | **nodename**: string
   | Name of one of the splitting nodes. The column is split across this
     node and the one on the opposite side of the column. If the
     specified node is not in the column, the method returns ``False``
     and nothing is done to the column.

-  | **chars**: string
   | Specifies a string of characters to use when forming the character
     part of block names. Default is lowercase letters.

----

.. _sec:mulgrid:translate:

``translate(shift, wells=False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; translating

Translates a grid by a specified shift in the *x*, *y* and *z*
directions. If the ``wells`` parameter is ``True``, any wells in the
grid are also translated.

**Parameters:**

-  | **shift**: list, tuple or ``np.array``
   | Distance to shift the grid in the *x*, *y* and *z* directions (of
     length 3).

-  | **wells**: Boolean
   |  Set ``True`` to translate wells.

**Example:**

::

   geo.translate([10.e3, 0.0, -1000.0])

translates the grid ``geo`` by 10 km in the *x* direction and down 1 km
in the *z* direction.

----

.. _sec:mulgrid:well_values:

``well_values(well_name, variable, divisions=1, elevation=False, deviations=False, qtree=None, extend=False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns values of a specified variable down a specified well. The
variable can be a list or ``np.array`` containing a value for every
block in the grid. The number of divisions between layer centres or
along each well deviation (default 1) can be optionally specified (this
can be increased to capture detail along a deviation that passes through
several blocks). If ``deviations`` is ``True``, values will be returned
at the nodes of the well track, instead of at grid layer centres. If
``extend`` is ``True``, the well trace is artificially extended to the
bottom of the model.

The routine returns a tuple of two arrays (``d``,\ ``v``), the first
(``d``) containing the measured depth down the well (or elevation if the
``elevation`` parameter is set to ``True``), and the second (``v``)
containing the value of the variable at each point. The value of the
variable at any point is the (block average) value at the block
containing the point.

**Parameters:**

-  | **well_name**: string
   | Name of the well.

-  | **variable**: list (or ``np.array``)
   | Variable to be plotted, of length equal to the number of blocks in
     the grid.

-  | **divisions**: integer
   | Number of divisions each well deviation is divided up into (default
     1).

-  | **elevation**: Boolean
   | Set to ``True`` if elevation rather than measured depth is to be
     returned.

-  | **deviations**: Boolean
   | Set to ``True`` to return values at deviation nodes, rather than
     intersections of layer centres with the well track.

-  | **qtree**: ``quadtree``
   | Quadtree object for fast searching of grid columns (can be
     constructed using the :ref:`column_quadtree() <sec:mulgrid:column_quadtree>`
     method).

-  | **extend**: Boolean
   | Set ``True`` to artificially extend the well trace to the bottom of
     the model.

----

.. _sec:mulgrid:write:

``write(filename='')``
^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; writing

Writes a ``mulgrid`` object to a MULgraph geometry file on disk.

**Parameters:**

-  | **filename**: string
   | Name of the MULgraph geometry file to be written. If no file name
     is specified, the object's own ``filename`` property is used.

----

.. _sec:mulgrid:write_bna:

``write_bna(filename='')``
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; exporting
           
Writes a geometry object to an Atlas BNA file on disk, for visualisation
with Surfer or GIS tools.

**Parameters:**

-  | **filename**: string
   | Name of the BNA file to be written. If no file name is specified,
     the object's own ``filename`` property is used, with the extension
     changed to \*.bna. If the object's ``filename`` property is not
     set, the default name 'geometry.bna' is used.

----

.. _sec:mulgrid:write_exodusii:

``write_exodusii(filename='', arrays=None, blockmap={})``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; exporting

Writes a ``mulgrid`` object to an ExodusII file on disk, for
visualisation or export to other software.

This method uses the VTK-Python library, so you will need that installed
on your machine before you can use it. An alternative is to use the
:ref:`write_mesh <sec:mulgrid:write_mesh>` method instead,
which can also write meshes to ExodusII format (as well as others), and
does not need the VTK-Python library (though you will need the
``meshio`` library).

**Parameters:**

-  | **filename**: string
   | Name of the ExodusII file to be written. If no file name is
     specified, the object's own ``filename`` property is used, with the
     extension changed to \*.exo. If the object's ``filename`` property
     is not set, the default name 'geometry.exo' is used.

-  | **arrays**: dictionary or ``None``
   | Data arrays to be included in the ExodusII file. If set to
     ``None``, default arrays (block name, layer index, column index,
     column area, column elevation, block number and volume) are
     included.

-  | **blockmap**: dictionary
   | Dictionary mapping the block names in the geometry to another block
     naming system.

----

.. _sec:mulgrid:write_mesh:

``write_mesh(filename, surface_snap = 0.1, dimension = 3, slice = None, file_format = None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; exporting

Writes a ``mulgrid`` object to a mesh file on disk, with the specific
format determined by the file extension of the specified
filename. This method uses the `meshio
<https://pypi.python.org/pypi/meshio>`_ library, which must be
installed on your machine, and supports various mesh output formats
including Dolfin XML, ExodusII, MSH, VTK, XDMF and others. The
``meshio`` library may be installed from PyPI (using e.g.  ``pip
install meshio``).

Note that many of these formats do not support columns with more than
four sides.

**Parameters:**

-  | **filename**: string
   | Name of the mesh file to be written.

-  | **surface_snap**: float
   | Tolerance for eliminating elements with very small vertical
     thickness at the top of the mesh (3-D meshes only).

-  | **dimension**: integer
   | Dimension of the mesh: when set to 3 (the default), write the full
     3-D mesh. When set to 2, write a 2-D mesh, corresponding either to
     the horizontal mesh only (the default), or a vertical slice mesh if
     the ``slice`` parameter is used.

-  | **slice**: list, string, float or ``None``
   | Horizontal line defining the slice for vertical 2-D meshes. This
     can be a list of two horizontal (*x*,\ *y*) points (``np.arrays``)
     defining the endpoints of the slice line, or string 'x' or 'y' to
     specify the *x*- or *y*-axis, or northing (float) through grid
     centre. If set to ``None`` (the default) then the horizontal 2-D
     mesh is written.

-  | **file_format**: string or ``None``
   | File format for mesh output. If ``None``, the file format will be
     decided from the filename extension (e.g. if the filename is
     'mesh.exo' then the mesh will be written in ExodusII format). See
     the ``meshio`` documentation for details.

----

.. _sec:mulgrid:write_vtk:

``write_vtk(filename='', arrays=None, wells=False, blockmap={}, surface_snap=0.1)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; exporting

Writes a ``mulgrid`` object to a VTK file on disk, for visualisation
with VTK, Paraview, Mayavi etc. The grid is written as an 'unstructured
grid' VTK object with optional data arrays defined on cells. A separate
VTK file for the wells in the grid can optionally be written.

**Parameters:**

-  | **filename**: string
   | Name of the VTK file to be written. If no file name is specified,
     the object's own ``filename`` property is used, with the extension
     changed to \*.vtu. If the object's ``filename`` property is not
     set, the default name 'geometry.vtu' is used.

-  | **arrays**: dictionary or ``None``
   | Data arrays to be included in the VTK file. If set to ``None``,
     default arrays (block name, layer index, column index, column area,
     column elevation, block number and volume) are included.

-  | **wells**: Boolean
   | If set to ``True``, a separate VTK file is written representing the
     wells in the grid.

-  | **blockmap**: dictionary
   | Dictionary mapping the block names in the geometry to another block
     naming system.

-  | **surface_snap**: float
   | Tolerance for specifying how close column surface elevations need
     to be before being considered "equal" when constructing surface
     nodes.

----

.. _other_mulgrid_objects:

Other objects (``node``, ``column``, ``layer``, ``connection`` and ``well``)
----------------------------------------------------------------------------

A ``mulgrid`` object contains lists of other types of objects:
:ref:`node <nodeobjects>`, :ref:`column <columnobjects>`,
:ref:`layer <layerobjects>`,
:ref:`connection <connectionobjects>` and
:ref:`well <wellobjects>` objects. These classes are described below.

.. _nodeobjects:

``node`` objects
~~~~~~~~~~~~~~~~

.. index:: MULgraph geometry; nodes
.. index:: nodes

A ``node`` object represents a node (i.e. vertex) in a ``mulgrid``
object. A ``node`` object has three properties: ``name``, which is a
string property containing the name of the node, ``pos`` which is an
``np.array`` with three elements, containing the node's position in 3D,
and ``column`` which is a set of the columns the node belongs to. A
``node`` object does not have any methods.

A ``node`` object ``n`` can be created for example using the command
``n = node(name,pos)`` where ``name`` is the node name and pos is an
``np.array`` (or list, or tuple) representing the node's position.

.. _columnobjects:

``column`` objects
~~~~~~~~~~~~~~~~~~

.. index:: MULgraph geometry; columns
.. index:: columns

A ``column`` object represents a column in a ``mulgrid`` object. The
properties of a ``column`` object are listed in the
:ref:`table <tb:column_properties>` below.

.. container::
   :name: tb:column_properties

   .. table:: Properties of a ``column`` object

      +--------------------+--------------+--------------------------------+
      | **Property**       | **Type**     | **Description**                |
      +====================+==============+================================+
      | ``angle_ratio``    | float        | ratio of largest to smallest   |
      |                    |              | interior angles                |
      +--------------------+--------------+--------------------------------+
      | ``area``           | float        | horizontal area of the column  |
      +--------------------+--------------+--------------------------------+
      | ``centre``         | ``np.array`` | horizontal centre of the       |
      |                    |              | column                         |
      +--------------------+--------------+--------------------------------+
      | ``centroid``       | ``np.array`` | average position of the        |
      |                    |              | column's vertices              |
      +--------------------+--------------+--------------------------------+
      | ``connection``     | set          | connections the column is in   |
      +--------------------+--------------+--------------------------------+
      | ``name``           | string       | name of the column             |
      +--------------------+--------------+--------------------------------+
      | ``neighbour``      | set          | set of neighbouring columns    |
      +--------------------+--------------+--------------------------------+
      | ``neighbourlist``  | list         | ordered list of neighbouring   |
      |                    |              | columns                        |
      +--------------------+--------------+--------------------------------+
      | ``node``           | list         | list of nodes (vertices)       |
      |                    |              | belonging to the column        |
      +--------------------+--------------+--------------------------------+
      | ``num_neighbours`` | integer      | number of neighbouring columns |
      +--------------------+--------------+--------------------------------+
      | ``num_nodes``      | integer      | number of nodes belonging to   |
      |                    |              | the column                     |
      +--------------------+--------------+--------------------------------+
      | ``num_layers``     | integer      | number of layers in the column |
      |                    |              | below the ground surface       |
      +--------------------+--------------+--------------------------------+
      | ``side_ratio``     | float        | ratio of largest to smallest   |
      |                    |              | side length                    |
      +--------------------+--------------+--------------------------------+
      | ``surface``        | float        | surface elevation of the       |
      |                    |              | column (``None`` if not        |
      |                    |              | specified)                     |
      +--------------------+--------------+--------------------------------+

The main properties defining a column are its ``name`` and ``node``
properties. The ``name`` is specified according to the naming
convention of the ``mulgrid`` object that the column belongs to. The
``node`` property is a list of ``node`` objects (not node names) that
belong to the column. A ``column``\ 's ``neighbour`` property is a set
of other ``columns`` connected to that column via a :ref:`connection
<connectionobjects>`), and its property is a set of connections the
column is part of. The ``neighbourlist`` property is a list of
neighbouring columns, with each item corresponding to a column edge
(``None`` if the edge is on a grid boundary). A ``column``\ 's
``centroid`` property returns the average of the positions of its
vertices - which is what the ``centre`` property is set to, unless
otherwise specified.

A ``column`` object has two properties measuring 'grid quality'. The
``angle_ratio`` property returns the ratio of largest to smallest
interior angles in the column. The ``side_ratio`` property returns the
ratio of largest to smallest side lengths (a generalisation of 'aspect
ratio' to columns with any number of sides). Values as close as possible
to 1.0 for both these measures are desirable (their values are both
exactly 1.0 for any regular polygon, e.g. an equilateral triangle or
square). Columns with large angle ratios will be highly skewed, while
those with large side ratios will be typically highly elongated in one
direction.

A ``column`` object ``col`` can be created for example using the
command:

::

   col = column(name, nodes, centre, surface)

where ``name`` is the column name and ``nodes`` is a list of
:ref:`node <nodeobjects>` objects defining the column. The
``centre`` and ``surface`` parameters are optional.

The methods of a ``column`` object are listed in the :ref:`table
<tb:column_methods>` below.

.. container::
   :name: tb:column_methods

   .. table:: Methods of a ``column`` object

      +------------------------------------------------------+---------------------------+--------------------------+
      | **Method**                                           | **Type**                  | **Description**          |
      +======================================================+===========================+==========================+
      | :ref:`contains_point() <sec:column:contains_point>`  | Boolean                   | if column contains point |
      |                                                      |                           |                          |
      |                                                      |                           |                          |
      +------------------------------------------------------+---------------------------+--------------------------+
      | :ref:`in_polygon() <sec:column:in_polygon>`          |          Boolean          | if column centre is      |
      |                                                      |                           | within a given polygon   |
      |                                                      |                           |                          |
      +------------------------------------------------------+---------------------------+--------------------------+
      |     :ref:`is_against() <sec:column:is_against>`      | Boolean                   | if two columns are       |
      |                                                      |                           | adjacent                 |
      |                                                      |                           |                          |
      +------------------------------------------------------+---------------------------+--------------------------+

----

.. _sec:column:contains_point:

``contains_point(pos)``
^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; searching

Returns ``True`` if a 2D point lies inside the column, and ``False``
otherwise.

**Parameters:**

-  | **pos**: ``np.array``
   | Horizontal position of the point.

----

.. _sec:column:in_polygon:

``in_polygon(polygon)``
^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; searching

Returns ``true`` if the column centre is inside the specified polygon or
rectangle.

**Parameters:**

-  | **polygon**: list (of ``np.array``)
   | List of points defining the polygon (each point is a two-element
     ``np.array``). If the list has only two points, it will be
     interpreted as a rectangle [bottom left, top right].

----

.. _sec:column:is_against:

``is_against(othercolumn)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns ``true`` if the column is 'against' ``othercolumn`` – that is,
if it shares more than one node with it.

**Parameters:**

-  | **othercolumn**: ``column``)
   | Any other column in the geometry.

----

.. _layerobjects:

``layer`` objects
~~~~~~~~~~~~~~~~~

.. index:: MULgraph geometry; layers
.. index:: layers

A ``layer`` object represents a layer in a ``mulgrid`` object. The
properties of a ``layer`` object are given in the 
:ref:`table <tb:layer_properties>` below.

.. container::
   :name: tb:layer_properties

   .. table:: Properties of a ``layer`` object

      +-------------+----------+--------------------------------------+
      |**Property** | **Type** | **Description**                      |
      +=============+==========+======================================+
      |``bottom``   | float    | elevation of the bottom of the layer |
      +-------------+----------+--------------------------------------+
      |``centre``   | float    | elevation of the centre of the layer |
      +-------------+----------+--------------------------------------+
      |``thickness``| float    | layer thickness (top - bottom)       |
      +-------------+----------+--------------------------------------+
      |``top``      | float    | elevation of the top of the layer    |
      +-------------+----------+--------------------------------------+
      |``name``     | string   | name of the layer                    |
      +-------------+----------+--------------------------------------+

A ``layer`` object ``lay`` can be created for example using the command:

::

   lay = layer(name, bottom, centre, top)

where ``name`` is the layer name and ``bottom``, ``centre`` and ``top``
specify the vertical position of the layer.

The methods of a ``layer`` object are given in the
:ref:`table <tb:layer_methods>` below.

.. container::
   :name: tb:layer_methods

   .. table:: Methods of a ``layer`` object

      +------------------------------------------------------------+---------------------------+--------------------------+
      |**Method**                                                  |**Type**                   |**Description**           |
      +============================================================+===========================+==========================+
      |:ref:`contains_elevation() <sec:layer:contains_elevation>`  |Boolean                    |if layer contains         |
      |                                                            |                           |elevation                 |
      |                                                            |                           |                          |
      +------------------------------------------------------------+---------------------------+--------------------------+
      |:ref:`translate() <sec:layer:translate>`                    |–                          |translate layer up or down|
      |                                                            |                           |                          |
      |                                                            |                           |                          |
      +------------------------------------------------------------+---------------------------+--------------------------+

----

.. _sec:layer:contains_elevation:

``contains_elevation(z)``
^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: MULgraph geometry; searching

Returns ``True`` if a point at a given elevation lies inside the layer,
and ``False`` otherwise.

**Parameters:**

-  | **z**: float
   | Elevation of the point.

----

.. _sec:layer:translate:

``translate(shift)``
^^^^^^^^^^^^^^^^^^^^

Translates a layer up or down by a specified distance.

**Parameters:**

-  | **shift**: float
   | Distance to shift the layer (positive for up, negative for down).

----

.. _connectionobjects:

``connection`` objects
~~~~~~~~~~~~~~~~~~~~~~

.. index:: MULgraph geometry; connections
.. index:: connections

A ``connection`` object represents a connection between ``columns`` in a
``mulgrid`` object. It has three properties: ``column``, which contains
a two-element list of the ``column`` objects making up the connection,
``node``, which contains a two-element list of the ``nodes`` on the face
joining the two columns in the connection, and ``angle_cosine``, which
gives the cosine of the angle between a line joining the nodes in the
connection and a line joining the centres of the two columns. This is
used as a measure of grid quality, these two lines should ideally be as
close to perpendicular as possible, making the cosine of the angle zero.
A ``connection`` has no methods.

A ``connection`` object ``con`` can be created for example using the
command

::

   con = connection(cols)

where ``cols`` is a two-element list of the ``column`` objects in the
connection.

.. _wellobjects:

``well`` objects
~~~~~~~~~~~~~~~~

.. index:: MULgraph geometry; wells
.. index:: wells

A ``well`` object represents a well in a ``mulgrid`` object. The
properties of a ``well`` object are given in the
:ref:`table <tb:well_properties>` below.

.. container::
   :name: tb:well_properties

   .. table:: Properties of a ``well`` object

      +--------------------+--------------+--------------------------------------------+
      | **Property**       | **Type**     | **Description**                            |
      +====================+==============+============================================+
      | ``bottom``         | ``np.array`` | well bottom position                       |
      +--------------------+--------------+--------------------------------------------+
      | ``deviated``       | Boolean      | whether well is deviated                   |
      +--------------------+--------------+--------------------------------------------+
      | ``head``           | ``np.array`` | well head position                         |
      +--------------------+--------------+--------------------------------------------+
      | ``name``           | string       | well name                                  |
      +--------------------+--------------+--------------------------------------------+
      | ``num_deviations`` | integer      | number of deviations                       |
      +--------------------+--------------+--------------------------------------------+
      | ``num_pos``        | integer      | number of well track nodes                 |
      +--------------------+--------------+--------------------------------------------+
      | ``pos``            | list         | positions (3-D arrays) of well track nodes |
      +--------------------+--------------+--------------------------------------------+
      | ``pos_depth``      | ``np.array`` | downhole depths along well track           |
      +--------------------+--------------+--------------------------------------------+

The well track can be deviated, and is defined as a list ``pos`` of (at
least two) 3D positions (``np.arrays``). The ``num_deviations`` property
returns the number of deviations in the track (one less than the
``num_pos`` property, which is the number of nodes in the ``pos`` list).
The ``deviated`` property returns ``True`` if there is more than one
deviation. The ``pos_depth`` property returns an array of the downhole
depths at each node along the well track.

A ``well`` object ``w`` can be created simply with the command
``w = well(name,pos)``, where ``name`` is the well name and ``pos`` is a
list of 3-element ``np.arrays`` (or lists, or tuples) representing the
well trace (starting from the wellhead).

The methods of a ``well`` object are listed in the
:ref:`table <tb:well_methods>` and described below.

.. container::
   :name: tb:well_methods

   .. table:: Methods of a ``well`` object

      +---------------------------------------------------+---------------------------+--------------------------+
      | **Method**                                        | **Type**                  | **Description**          |
      +===================================================+===========================+==========================+
      | :ref:`depth_elevation <sec:well:depth_elevation>` | float                     | elevation for a given    |
      |                                                   |                           | downhole depth           |
      |                                                   |                           |                          |
      +---------------------------------------------------+---------------------------+--------------------------+
      | :ref:`depth_pos <sec:well:depth_pos>`             |       ``np.array``        | position on well track   |
      |                                                   |                           | for a given downhole     |
      |                                                   |                           | depth                    |
      +---------------------------------------------------+---------------------------+--------------------------+
      | :ref:`elevation_depth <sec:well:elevation_depth>` | float                     | downhole depth for a     |
      |                                                   |                           | given elevation          |
      |                                                   |                           |                          |
      +---------------------------------------------------+---------------------------+--------------------------+
      | :ref:`elevation_pos <sec:well:elevation_pos>`     | ``np.array``              | position on well track   |
      |                                                   |                           | for a given elevation    |
      |                                                   |                           |                          |
      +---------------------------------------------------+---------------------------+--------------------------+
      | :ref:`pos_coordinate <sec:well:pos_coordinate>`   | ``np.array``              | array of coordinates for |
      |                                                   |                           | a given index            |
      |                                                   |                           |                          |
      +---------------------------------------------------+---------------------------+--------------------------+

.. _sec:well:depth_elevation:

``depth_elevation(depth)``
^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns the elevation corresponding to the specified downhole ``depth``
(or ``None`` if ``depth`` is above the wellhead or below the bottom).

**Parameters:**

-  | **depth**: float
   | Downhole depth.

----

.. _sec:well:depth_pos:

``depth_pos(depth)``
^^^^^^^^^^^^^^^^^^^^

Returns the 3D position of the point in the well with specified downhole
``depth`` (or ``None`` if ``depth`` is above the wellhead or below the
bottom). The position is interpolated between the deviation locations.

**Parameters:**

-  | **depth**: float
   | Downhole depth of the required point.

----

.. _sec:well:elevation_depth:

``elevation_depth(elevation)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns the downhole depth corresponding to the specified ``elevation``
(or ``None`` if ``elevation`` is above the wellhead or below the
bottom).

**Parameters:**

-  | **elevation**: float
   | Elevation.

----

.. _sec:well:elevation_pos:

``elevation_pos(elevation, extend=False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns the 3D position of the point in the well with specified
``elevation`` (or ``None`` if ``elevation`` is above the wellhead or
below the bottom). The position is interpolated between the deviation
locations. If ``extend`` is ``True``, return extrapolated positions for
elevations below the bottom of the well.

**Parameters:**

-  | **elevation**: float
   | Elevation of the required point.

-  | **extend**: Boolean
   | If ``True``, extrapolated positions will be returned for elevations
     below the bottom of the well (otherwise ``None`` will be returned).

----

.. _sec:well:pos_coordinate:

``pos_coordinate(index)``
^^^^^^^^^^^^^^^^^^^^^^^^^

Returns an ``np.array`` of the well track node coordinates for the given
index (0, 1 or 2). For example, ``pos_coordinate(2)`` returns an array
containing the elevations of all well track nodes.

**Parameters:**

-  | **index**: integer
   | Index required (0, 1 or 2).

----

Other functions: block name conversions
---------------------------------------

The ``mulgrids`` library contains two other functions connected with
working with geometry files and TOUGH2 grids:

.. _sec:mulgrid:fix_blockname:

``fix_blockname(name)``
~~~~~~~~~~~~~~~~~~~~~~~

.. index:: MULgraph geometry; block names

TOUGH2 always assumes that the last two characters of a block name
represent a two-digit number. However, if that number is less than 10,
the fourth character is not padded with zeros, so for example 'AA101'
becomes 'AA1 1' when processed by TOUGH2.

The ``fix_blockname`` function corrects this by padding the fourth
character of a block name with a zero if necessary. This is only done if
the third character is also a digit, e.g. when naming convention 2 is
used (two characters for layer followed by three digits for column).

**Parameters:**

-  | **name**: string
   | Block name.

----

``unfix_blockname(name)``
~~~~~~~~~~~~~~~~~~~~~~~~~

.. index:: MULgraph geometry; block names

This function reverses the effect of ``fix_blockname()``.

**Parameters:**

-  | **name**: string
   | Block name.

----

.. _sec:mulgrid:blockmappings:

Block mappings: handling other block naming conventions
-------------------------------------------------------

.. index:: MULgraph geometry; block mappings

The MULgraph geometry format names blocks according to one of its
:ref:`naming conventions <geometry_format_conventions>`.
All of these conventions use part of the block name to indicate the layer
and part of it to indicate the column.

However, in PyTOUGH it is possible to make a ``mulgrid`` object handle
other block naming conventions by means of a **block mapping**. This is
simply a dictionary that maps the block names in a ``mulgrid`` to block
names in a ``t2grid`` object. The block names in the ``t2grid`` can
follow an arbitrary convention, not based on layers and columns. For
example, blocks in TOUGH2 grids created by PetraSim may be simply
numbered.

A block mapping dictionary can be passed in as an optional parameter to
many PyTOUGH methods that involve both a MULgraph geometry and TOUGH2
grid, for example the ``mulgrid`` :ref:`block_name() <sec:mulgrid:block_name>`,
:ref:`slice_plot() <sec:mulgrid:slice_plot>` and
:ref:`write_vtk() <sec:mulgrid:write_vtk>` methods, and the ``write_vtk()``
methods of the :ref:`t2grid <sec:t2grid:write_vtk>` and
:ref:`t2listing <sec:t2listing:write_vtk>` classes.

When the :ref:`rectgeo() <sec:t2grid:rectgeo>` method is used
to create a ``mulgrid`` object from a ``t2grid``, a block mapping is
also created, and may be used in the PyTOUGH methods that can accept a
block mapping.

A block mapping need not contain entries for all blocks. If for example
a model follows the naming convention of a MULgraph geometry in most
blocks, and only a few are different, then only entries for the
different block names need be present in the mapping dictionary.

Block mappings can be saved to and loaded from disk (like any other
Python object) using the ``pickle`` library. This is part of the
standard Python library collection. For example a block mapping called
``blockmap`` can be saved to a file called ``'blockmap.pkl'`` as
follows:

::

     import pickle
     pickle.dump(blockmap, open('blockmap.pkl', 'w'))

It can be loaded back in again like this:

::

     blockmap = pickle.load(open('blockmap.pkl'))

