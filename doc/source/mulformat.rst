:tocdepth: 3

.. _geometry_file_format:

MULgraph geometry file format
=============================

.. _introduction-8:

Introduction
------------

This section gives a format specification of the MULgraph geometry
file. These files can be used to give a geometrical description of a
TOUGH2 model grid, useful for creating grids and visualizing simulation
results.

MULgraph geometry files were originally developed for use with
MULgraph, a graphical interface for TOUGH2 and AUTOUGH2 developed at
the University of Auckland in the 1990s, and subsequently adopted by
the `TIM <https://tim.readthedocs.io/>`_ graphical interface. However,
MULgraph geometry files can be used independently of MULgraph or
TIM. PyTOUGH is able to represent the contents of a MULgraph geometry
file in a Python script via the :ref:`mulgrid <mulgrids>` class.

Grid structure
--------------

Layers and columns
~~~~~~~~~~~~~~~~~~

MULgraph geometry files implicitly assume a layered structure, with
blocks arranged in layers and columns, and the same arrangement of
columns in each layer. The only exception to this is at the top surface
of the model, where layers are allowed to be incomplete (i.e. not
contain all columns) in order to represent topography.

The layers are always of constant vertical thickness. However, the
blocks in the top layer are allowed to vary in height, again to
represent variations in ground surface elevation.

Atmosphere blocks
~~~~~~~~~~~~~~~~~

The blocks in the top layer may optionally be connected to the
atmosphere- either a single atmosphere block connected to all columns,
or a separate atmosphere block over each column (see
:ref:`Naming conventions and atmosphere types <geometry_format_conventions>`).

.. _tilted-geometries-1:

Tilted geometries
~~~~~~~~~~~~~~~~~

It is possible to tilt the geometry coordinate axes with respect to the
vertical, to represent non-horizontal geometries. When a TOUGH2 grid is
created from such a tilted geometry, only the gravity cosines of the
grid connections are affected.

.. _rotating-permeability-directions-1:

Rotating permeability directions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is also possible to rotate the permeability principal directions with
respect to the coordinate axes- for example, to align permeabilities
with a dominant fault direction. When a TOUGH2 grid is created, this can
change the permeability index associated with each connection.

Geometry types
--------------

The original MULgraph file specification allowed for three types of
geometry: 'general', 'rectangular' and 'radial'. Only the 'general'
geometry type is supported by PyTOUGH. It is intended for representing
general grids with arbitrary, possibly unstructured horizontal column
arrangements.

The 'rectangular' type was a special type for grids with rectangular
horizontal column structures. These can also be represented using the
'general' geometry type. Since PyTOUGH contains
:ref:`methods <sec:mulgrid:rectangular>` for constructing
rectangular grids within the 'general' geometry type, there is usually
no longer any significant benefit from using the 'rectangular' type.

The 'radial' type was intended for grids with radial horizontal column
structure. PyTOUGH also contains :ref:`methods <sec:t2grid:radial>`
for creating radial TOUGH2 grids. Simulation results from radial
models can also be visualized using a simple one- or two-dimensional
rectangular 'general' geometry to represent the grid structure in the
radial direction.

.. _geometry_format_conventions:

Naming conventions and atmosphere types
---------------------------------------

The grid block naming convention and atmosphere type used in a
MULgraph geometry file are both integers and can be given values in
the range 0 -- 3 and 0 -- 2 respectively.  The meanings of these
values are shown in the tables below.

Note that the grid nodes (vertices) are also named according to the
column part of the block naming convention. If naming nodes, columns or
layers manually, while the names can in principle be arbitrary (within
the naming convention), it is safest to right-justify them.

The MULgraph block naming conventions all use part of the block name to
indicate the layer, and part of it to indicate the column. In PyTOUGH,
it is also possible to use MULgraph geometry files in conjunction with
TOUGH2 grids that follow other naming conventions, by means of a
:ref:`block mapping <sec:mulgrid:blockmappings>` dictionary.

Block naming convention 3 was not supported by the original MULgraph
geometry file format, and produces block names which do not conform to
the TOUGH2 block naming requirements (having numbers in the last two
characters). It can be used to produce grids for other simulators such
as Waiwera which do not have these requirements. An alternative tool
for creating such grids is the `Layermesh
<https://github.com/acroucher/layermesh>`_ library.

.. container::
   :name: tb:mulgrid_conventions

   .. table:: MULgraph geometry file naming conventions

      +------------+-------------------------------------------------------+
      | Convention | Meaning                                               |
      +============+=======================================================+
      | 0          |3 characters for column followed by 2 digits for layer |
      +------------+-------------------------------------------------------+
      | 1          |3 characters for layer followed by 2 digits for column |
      +------------+-------------------------------------------------------+
      | 2          |2 characters for layer followed by 3 digits for column |
      +------------+-------------------------------------------------------+
      | 3          |3 characters for column followed by 2 characters for   |
      |            |layer                                                  |
      +------------+-------------------------------------------------------+

.. container::
   :name: tb:mulgrid_atmosphere_types

   .. table:: MULgraph geometry file atmosphere types

      +------+---------------------------------------+
      | Type | Meaning                               |
      +======+=======================================+
      | 0    | A single atmosphere block             |
      +------+---------------------------------------+
      | 1    | One atmosphere block over each column |
      +------+---------------------------------------+
      | 2    | No atmosphere blocks                  |
      +------+---------------------------------------+

File format
-----------

MUlgraph geometry files are simple formatted ASCII text files with a
header line at the top, followed by a number of sections. Each section
begins with a keyword and ends with a blank line. Each line has
**fixed** format, so the different values have to be specified in the
right text columns.

If you use PyTOUGH scripts to create and manipulate your grid
geometries, you don't need to know anything about the format of a
MULgraph geometry file, because PyTOUGH will handle reading and writing
them for you. If, however, for some reason you do need to know how these
files are structured, the format specification for a 'general' type
geometry file is given below.

Header
~~~~~~

This is a single line containing a number of global parameters of the
geometry. Its format is given in the
:ref:`table <tb:mulgraph_format_header>` below.

Note that the block ordering parameter is an extension to the original
MULgraph file format.

.. container::
   :name: tb:mulgraph_format_header

   .. table:: MULgraph geometry file header line format

      +---------------+-----------+------------+-------------+-------------------------------------------+
      | **Name**      | **Type**  | **Length** | **Columns** | **Description**                           |
      |               |           |            |             |                                           |
      +===============+===========+============+=============+===========================================+
      | **Geometry    | character | 5          | 1–5         | 'GENER' for general geometry type;        |
      | type**        |           |            |             | 'RECTA' or 'RADIA' for other types        |
      |               |           |            |             | (but these are not supported by           |
      |               |           |            |             | PyTOUGH)                                  |
      +---------------+-----------+------------+-------------+-------------------------------------------+
      | **Naming      | integer   | 1          | 6           | Block naming                              |
      | convention**  |           |            |             | :ref:`convention<tb:mulgrid_conventions>` |
      |               |           |            |             |                                           |
      |               |           |            |             |                                           |
      |               |           |            |             |                                           |
      +---------------+-----------+------------+-------------+-------------------------------------------+
      | **Atmosphere  | integer   | 1          | 7           | :ref:`Type <tb:mulgrid_atmosphere_types>` |
      | type**        |           |            |             | of atmosphere                             |
      |               |           |            |             |                                           |
      |               |           |            |             |                                           |
      |               |           |            |             |                                           |
      |               |           |            |             |                                           |
      +---------------+-----------+------------+-------------+-------------------------------------------+
      | **Atmosphere  | float     | 10         | 8–17        | Volume of each atmosphere block           |
      | volume**      |           |            |             | (default 10\ :sup:`20` m\ :sup:`3`)       |
      +---------------+-----------+------------+-------------+-------------------------------------------+
      | **Atmosphere  | float     | 10         | 18–27       | Connection distance for each              |
      | connection    |           |            |             | atmosphere block (default                 |
      | distance**    |           |            |             | 10\ :sup:`-6` m)                          |
      +---------------+-----------+------------+-------------+-------------------------------------------+
      | **Length      | character | 5          | 28–32       | Default is metres (blank); for            |
      | unit**        |           |            |             | feet specify 'FEET'                       |
      +---------------+-----------+------------+-------------+-------------------------------------------+
      | **x-direction | float     | 10         | 33–42       | Cosine of angle between x-axis and        |
      | cosine**      |           |            |             | gravity vector (default zero); set        |
      |               |           |            |             | positive for tilt in the x-direction      |
      +---------------+-----------+------------+-------------+-------------------------------------------+
      | **y-direction | float     | 10         | 43–52       | Cosine of angle between                   |
      | cosine**      |           |            |             | y-axis and gravity vector (default        |
      |               |           |            |             | zero); set positive for tilt in the       |
      |               |           |            |             | y-direction                               |
      +---------------+-----------+------------+-------------+-------------------------------------------+
      | **Connection  | integer   | 1          | 53          | Method of calculating connection          |
      | type**        |           |            |             | parameters (default zero)- not            |
      |               |           |            |             | supported by PyTOUGH                      |
      +---------------+-----------+------------+-------------+-------------------------------------------+
      | **Permeability| float     | 10         | 54–63       | Horizontal angle (degrees                 |
      | angle**       |           |            |             | anti-clockwise) between first             |
      |               |           |            |             | permeability direction and x-axis         |
      +---------------+-----------+------------+-------------+-------------------------------------------+
      | **Block       | integer   | 2          | 64–65       | Block ordering scheme: 0 for original     |
      | ordering**    |           |            |             | MULgraph layer/column ordering; 1 for     |
      |               |           |            |             | PETSc DMPlex ordering (sorted by          |
      |               |           |            |             | block type)                               |
      +---------------+-----------+------------+-------------+-------------------------------------------+

Vertices
~~~~~~~~

This section defines the horizontal locations of the grid vertices
(nodes), at the corners of the columns. The first line just contains the
keyword 'VERTI'. Each subsequent line defines the position of a vertex,
and has the format given in the
:ref:`table <tb:mulgraph_format_vertices>` below. The vertices section is
terminated by a blank line.

.. container::
   :name: tb:mulgraph_format_vertices

   .. table:: MULgraph geometry file vertices format

      +--------------+-----------+------------+-------------+-------------------------------------------------+
      | **Name**     | **Type**  | **Length** | **Columns** |                 **Description**                 |
      |              |           |            |             |                                                 |
      +==============+===========+============+=============+=================================================+
      | **Vertex     | character | 3          | 1–3         | Name of the vertex (honouring the column naming |
      | name**       |           |            |             | :ref:`convention <tb:mulgrid_conventions>`      |
      |              |           |            |             |                                                 |
      |              |           |            |             |                                                 |
      |              |           |            |             |                                                 |
      |              |           |            |             |                                                 |
      |              |           |            |             |                                                 |
      |              |           |            |             |                                                 |
      |              |           |            |             |                                                 |
      |              |           |            |             |                                                 |
      +--------------+-----------+------------+-------------+-------------------------------------------------+
      | **x**        | float     | 10         | 4–13        | x-coordinate                                    |
      |              |           |            |             | of the                                          |
      |              |           |            |             | vertex                                          |
      +--------------+-----------+------------+-------------+-------------------------------------------------+
      | **y**        | float     | 10         | 14–23       | y-coordinate                                    |
      |              |           |            |             | of the                                          |
      |              |           |            |             | vertex                                          |
      +--------------+-----------+------------+-------------+-------------------------------------------------+

Grid
~~~~

This section specifies the vertices making up each column. The first
line just contains the keyword 'GRID'.

For each grid column, there is then a sub-header line with information
about the column, followed by a line for each vertex making up the
column. The lines for the sub-header and each vertex have the formats
given in the tables below. There are no blank lines between the
definitions of the grid columns, but there is a blank line at the end
of the section.

.. container::
   :name: tb:mulgraph_format_column_header

   .. table:: MULgraph geometry file column header format

      +--------------+-----------+------------+-------------+---------------------------------------------------+
      | **Name**     | **Type**  | **Length** | **Columns** |                  **Description**                  |
      |              |           |            |             |                                                   |
      +==============+===========+============+=============+===================================================+
      | **Column     | character | 3          | 1–3         | Name of the column (honouring the column naming   |
      | name**       |           |            |             | :ref:`convention <tb:mulgrid_conventions>`)       |
      |              |           |            |             |                                                   |
      |              |           |            |             |                                                   |
      |              |           |            |             |                                                   |
      |              |           |            |             |                                                   |
      |              |           |            |             |                                                   |
      |              |           |            |             |                                                   |
      |              |           |            |             |                                                   |
      |              |           |            |             |                                                   |
      |              |           |            |             |                                                   |
      +--------------+-----------+------------+-------------+---------------------------------------------------+
      | **Centre     | integer   | 1          | 4–5         | Set non-zero                                      |
      | specified**  |           |            |             | to specify                                        |
      |              |           |            |             | the column                                        |
      |              |           |            |             | centre                                            |
      |              |           |            |             | location, or                                      |
      |              |           |            |             | zero                                              |
      |              |           |            |             | (default) to                                      |
      |              |           |            |             | calculate it                                      |
      |              |           |            |             | as the                                            |
      |              |           |            |             | centroid of                                       |
      |              |           |            |             | the column                                        |
      +--------------+-----------+------------+-------------+---------------------------------------------------+
      | **Number of  | integer   | 2          | 6–7         | Number of                                         |
      | vertices**   |           |            |             | vertices in                                       |
      |              |           |            |             | the column                                        |
      +--------------+-----------+------------+-------------+---------------------------------------------------+
      | **Column     | float     | 10         | 8–17        | x-coordinate                                      |
      | centre x**   |           |            |             | of column                                         |
      |              |           |            |             | centre                                            |
      +--------------+-----------+------------+-------------+---------------------------------------------------+
      | **Column     | float     | 10         | 18–27       | y-coordinate                                      |
      | centre y**   |           |            |             | of column                                         |
      |              |           |            |             | centre                                            |
      +--------------+-----------+------------+-------------+---------------------------------------------------+

.. container::
   :name: tb:mulgraph_format_column_vertex

   .. table:: MULgraph geometry file column vertex format

      +--------------+-----------+------------+-------------+----------------+
      | **Name**     | **Type**  | **Length** | **Columns** | **Description**|
      |              |           |            |             |                |
      +==============+===========+============+=============+================+
      | **Vertex     | character | 3          | 1–3         | Name of the    |
      | name**       |           |            |             | vertex, as     |
      |              |           |            |             | specified in   |
      |              |           |            |             | the vertices   |
      |              |           |            |             | section        |
      +--------------+-----------+------------+-------------+----------------+

Connections
~~~~~~~~~~~

This section defines the horizontal connections between columns. The
first line just contains the keyword 'CONNE'.

Each subsequent line defines a connection between two columns, and has
the format given in the :ref:`table <tb:mulgraph_format_connection>` below.
There is a blank line at the end of the section.

.. container::
   :name: tb:mulgraph_format_connection

   .. table:: MULgraph geometry file connection format

      +--------------+-----------+------------+-------------+----------------+
      | **Name**     | **Type**  | **Length** | **Columns** | **Description**|
      |              |           |            |             |                |
      +==============+===========+============+=============+================+
      | **First      | character | 3          | 1–3         | Name of the    |
      | column       |           |            |             | first column   |
      | name**       |           |            |             |                |
      +--------------+-----------+------------+-------------+----------------+
      | **Second     | character | 3          | 4–6         | Name of the    |
      | column       |           |            |             | second         |
      | name**       |           |            |             | column         |
      +--------------+-----------+------------+-------------+----------------+

Layers
~~~~~~

This section defines the grid layers. The first line just contains the
keyword 'LAYER'.

Each subsequent line defines a layer, with format given in the
:ref:`table <tb:mulgraph_format_layer>` below. There are no blank lines between
layers, but there is a blank line at the end of the section.

.. container::
   :name: tb:mulgraph_format_layer

   .. table:: MULgraph geometry file layer format

      +--------------+-----------+------------+-------------+----------------------------------------------------+
      | **Name**     | **Type**  | **Length** | **Columns** |                  **Description**                   |
      |              |           |            |             |                                                    |
      +==============+===========+============+=============+====================================================+
      | **Layer      | character | 3          | 1–3         | Name of the layer (honouring the layer             |
      | name**       |           |            |             | naming :ref:`convention <tb:mulgrid_conventions>`) |
      |              |           |            |             |                                                    |
      |              |           |            |             |                                                    |
      |              |           |            |             |                                                    |
      |              |           |            |             |                                                    |
      |              |           |            |             |                                                    |
      |              |           |            |             |                                                    |
      |              |           |            |             |                                                    |
      |              |           |            |             |                                                    |
      |              |           |            |             |                                                    |
      |              |           |            |             |                                                    |
      |              |           |            |             |                                                    |
      |              |           |            |             |                                                    |
      |              |           |            |             |                                                    |
      +--------------+-----------+------------+-------------+----------------------------------------------------+
      | **Bottom     | float     | 10         | 4–13        | Elevation of                                       |
      | elevation**  |           |            |             | the bottom                                         |
      |              |           |            |             | of the layer                                       |
      +--------------+-----------+------------+-------------+----------------------------------------------------+
      | **Centre     | float     | 10         | 14–23       | Elevation of                                       |
      | elevation**  |           |            |             | the centre                                         |
      |              |           |            |             | of the layer                                       |
      +--------------+-----------+------------+-------------+----------------------------------------------------+

Surface elevation
~~~~~~~~~~~~~~~~~

This section is optional, and can be used to define the surface
elevation at any or all columns in the grid, to represent topography.
The first line just contains the keyword 'SURFA'.

Each subsequent line defines the surface elevation at a column, with
format given in the :ref:`table <tb:mulgraph_format_surface>` below. There is a
blank line at the end of the section.

.. container::
   :name: tb:mulgraph_format_surface

   .. table:: MULgraph geometry file surface elevation format

      +--------------+-----------+------------+-------------+----------------+
      | **Name**     | **Type**  | **Length** | **Columns** | **Description**|
      |              |           |            |             |                |
      +==============+===========+============+=============+================+
      | **Column     | character | 3          | 1–3         | Name of the    |
      | name**       |           |            |             | column         |
      +--------------+-----------+------------+-------------+----------------+
      | **Surface    | float     | 10         | 4–13        | Surface        |
      | elevation**  |           |            |             | elevation of   |
      |              |           |            |             | the column     |
      +--------------+-----------+------------+-------------+----------------+

Wells
~~~~~

This section is optional, and can be used to define the positions of
wells (including their tracks) within the geometry. Deviated wells are
supported. The first line of the section just contains the keyword
'WELLS'.

Each subsequent line defines the location of one point on a well track,
with format given in the :ref:`table <tb:mulgraph_format_wells>` below. At
least two points are required to define each well (one for the wellhead
and one for the bottom), with more than two points needed to define a
deviated well. There is a blank line at the end of the section.

.. container::
   :name: tb:mulgraph_format_wells

   .. table:: MULgraph geometry file well format

      +--------------+-----------+------------+-------------+----------------+
      | **Name**     | **Type**  | **Length** | **Columns** |**Description** |
      |              |           |            |             |                |
      +==============+===========+============+=============+================+
      | **Well       | character | 5          | 1–5         | Name of the    |
      | name**       |           |            |             | well           |
      +--------------+-----------+------------+-------------+----------------+
      | **x**        | float     | 10         | 6–15        | x-coordinate   |
      |              |           |            |             | of the well    |
      |              |           |            |             | location       |
      +--------------+-----------+------------+-------------+----------------+
      | **y**        | float     | 10         | 16–25       | y-coordinate   |
      |              |           |            |             | of the well    |
      |              |           |            |             | location       |
      +--------------+-----------+------------+-------------+----------------+
      | **z**        | float     | 10         | 26–35       | z-coordinate   |
      |              |           |            |             | of the well    |
      |              |           |            |             | location       |
      +--------------+-----------+------------+-------------+----------------+
