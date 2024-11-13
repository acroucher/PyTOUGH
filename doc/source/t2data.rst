:tocdepth: 3

.. _datafiles:

TOUGH2 data files
=================

.. index:: TOUGH2 data files

.. _introduction-3:

Introduction
------------

The ``t2data`` library in PyTOUGH contains classes and routines for
creating, editing and saving TOUGH2 or AUTOUGH2 data files. It can be
imported using the command:

::

      from t2data import *

``t2data`` objects
------------------

The ``t2data`` library defines a ``t2data`` class, used for representing
TOUGH2 data files.

**Example:**

::

   dat = t2data()

creates an empty ``t2data`` object called ``dat``.

::

   dat = t2data(filename)

creates a ``t2data`` object called ``dat`` and reads its contents from
file ``filename``. (It is also possible to read the mesh part of the
``t2data`` object from separate files - see below.)

Because a ``t2data`` object contains a large number of different
parameters, it is usually easier to load one from an existing TOUGH2
data file and edit it, rather than creating a new one from scratch.

.. _properties-2:

Properties
~~~~~~~~~~

The main properties of a ``t2data`` object are listed in the
:ref:`table <tb:t2data_properties>` below. In general, each of
these properties corresponds to an input block in a TOUGH2 data file.
Most of these input blocks contain a number of different parameters, so
that the ``t2data`` property corresponding to each input block is
usually in the form of a dictionary, containing a number of keys
representing sub-properties.

For example, the maximum number of time steps for the simulation is
controlled by ``max_timesteps`` key in the ``parameter`` property, which
for a ``t2data`` object called ``dat`` would be accessed by
``dat.parameter['max_timesteps']``.

.. container::
   :name: tb:t2data_properties

   .. table:: Properties of a ``t2data`` object

      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | **Property**                                                    | **Type**               | **Description**  | **Input block**|
      |                                                                 |                        |                  |                |
      +=================================================================+========================+==================+================+
      | :ref:`capillarity <sec:t2data:capillarity>`                     | dictionary             | capillarity      | RELP           |
      |                                                                 |                        | function         |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`diffusion <sec:t2data:diffusion>`                         | list                   | diffusion        | DIFFU          |
      |                                                                 |                        | coefficients     |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`echo_extra_precision <sec:t2data:echo_extra_precision>`   | Boolean                | echoing extra    | –              |
      |                                                                 |                        | precision        |                |
      |                                                                 |                        | sections to      |                |
      |                                                                 |                        | main data file   |                |
      |                                                                 |                        | (AUTOUGH2        |                |
      |                                                                 |                        | only)            |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`end_keyword <sec:t2data:end_keyword>`                     | string                 | keyword to end   | ENDCY or ENDFI |
      |                                                                 |                        | file             |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`extra_precision <sec:t2data:extra_precision>`             | list                   | data sections    | –              |
      |                                                                 |                        | read from        |                |
      |                                                                 |                        | extra            |                |
      |                                                                 |                        | precision        |                |
      |                                                                 |                        | auxiliary file   |                |
      |                                                                 |                        | (AUTOUGH2        |                |
      |                                                                 |                        | only)            |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`filename <sec:t2data:filename>`                           | string                 | file name on     | –              |
      |                                                                 |                        | disk             |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`generator <sec:t2data:generator>`                         | dictionary             | generators (by   | GENER          |
      |                                                                 |                        | block name and   |                |
      |                                                                 |                        | generator        |                |
      |                                                                 |                        | name)            |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`generatorlist <sec:t2data:generatorlist>`                 | list                   | generators (by   | GENER          |
      |                                                                 |                        | index)           |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`grid <sec:t2data:grid>`                                   | :ref:`t2grid <t2grids>`| model grid       | ELEME, CONNE   |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`history_block <sec:t2data:history_block>`                 | list                   | history blocks   | FOFT           |
      |                                                                 |                        | (TOUGH2 only)    |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`history_connection <sec:t2data:history_connection>`       | list                   | history          | COFT           |
      |                                                                 |                        | connections      |                |
      |                                                                 |                        | (TOUGH2 only)    |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`history_generator <sec:t2data:history_generator>`         | list                   | history          | GOFT           |
      |                                                                 |                        | generators       |                |
      |                                                                 |                        | (TOUGH2 only)    |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`incon <sec:t2data:incon>`                                 | dictionary             | initial          | INCON          |
      |                                                                 |                        | conditions       |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`indom <sec:t2data:indom>`                                 | dictionary             | rocktype-specific| INDOM          |
      |                                                                 |                        | initial          |                |
      |                                                                 |                        | conditions       |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`lineq <sec:t2data:lineq>`                                 | dictionary             | linear           | LINEQ          |
      |                                                                 |                        | equation         |                |
      |                                                                 |                        | solver options   |                |
      |                                                                 |                        | (AUTOUGH2        |                |
      |                                                                 |                        | only)            |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`meshfilename <sec:t2data:meshfilename>`                   | string or              | file name(s)     | –              |
      |                                                                 | tuple                  | on disk          |                |
      |                                                                 |                        | containing       |                |
      |                                                                 |                        | mesh data        |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`meshmaker <sec:t2data:meshmaker>`                         | list                   | mesh             | MESHM          |
      |                                                                 |                        | generation       |                |
      |                                                                 |                        | options          |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`more_option <sec:t2data:more_options>`                    | array of               | additional       | MOMOP          |
      |                                                                 | integer                | parameter        |                |
      |                                                                 |                        | options          |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`multi <sec:t2data:multi>`                                 | dictionary             | EOS              | MULTI          |
      |                                                                 |                        | configuration    |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`noversion <sec:t2data:noversion>`                         | Boolean                | suppressing      | NOVER          |
      |                                                                 |                        | printing of      |                |
      |                                                                 |                        | version          |                |
      |                                                                 |                        | summary          |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`num_generators <sec:t2data:num_generators>`               | integer                | number of        | –              |
      |                                                                 |                        | generators       |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`output_times <sec:t2data:output_times>`                   | dictionary             | times to write   | TIMES          |
      |                                                                 |                        | output           |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`parameter <sec:t2data:parameter>`                         | dictionary             | run-time         | PARAM          |
      |                                                                 |                        | parameters       |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`relative_permeability <sec:t2data:relative_permeability>` | dictionary             | relative         | RELP           |
      |                                                                 |                        | permeability     |                |
      |                                                                 |                        | function         |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`selection <sec:t2data:selection>`                         | dictionary             | selection        | SELEC          |
      |                                                                 |                        | parameters       |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`short_output <sec:t2data:short_output>`                   | dictionary             | short output     | SHORT          |
      |                                                                 |                        | (AUTOUGH2        |                |
      |                                                                 |                        | only)            |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`simulator <sec:t2data:simulator>`                         | string                 | simulator name   | SIMUL          |
      |                                                                 |                        | (AUTOUGH2        |                |
      |                                                                 |                        | only)            |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`solver <sec:t2data:solver>`                               | dictionary             | linear           | SOLVR          |
      |                                                                 |                        | equation         |                |
      |                                                                 |                        | solver options   |                |
      |                                                                 |                        | (TOUGH2 only)    |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`start <sec:t2data:start>`                                 | Boolean                | run              | START          |
      |                                                                 |                        | initialisation   |                |
      |                                                                 |                        | option           |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`title <sec:t2data:title>`                                 | string                 | simulation       | TITLE          |
      |                                                                 |                        | title            |                |
      |                                                                 |                        |                  |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+
      | :ref:`type <sec:t2data:type>`                                   | string                 | simulator type   | –              |
      |                                                                 |                        | (AUTOUGH2 or     |                |
      |                                                                 |                        | TOUGH2)          |                |
      |                                                                 |                        |                  |                |
      +-----------------------------------------------------------------+------------------------+------------------+----------------+

The details of the ``t2data`` properties are as follows.

----

.. _sec:t2data:capillarity:

``capillarity`` property
^^^^^^^^^^^^^^^^^^^^^^^^

A dictionary property specifying the capillarity function used,
corresponding to the second line of the **RPCAP** input block in the
TOUGH2 data file. The individual keys of this property are given in
the :ref:`table <tb:capillarity>` below.

.. container::
   :name: tb:capillarity

   .. table:: ``capillarity`` property keys

      +----------------+----------------+----------------+----------------------+
      | **Key**        | **Type**       | **Description**| **TOUGH2 parameter** |
      |                |                |                |                      |
      +================+================+================+======================+
      | ``parameters`` | array (7) of   | function       | CP                   |
      |                | float          | parameters     |                      |
      +----------------+----------------+----------------+----------------------+
      | ``type``       | integer        | type of        | ICP                  |
      |                |                | capillarity    |                      |
      |                |                | function       |                      |
      +----------------+----------------+----------------+----------------------+

----

.. _sec:t2data:diffusion:

``diffusion`` property
^^^^^^^^^^^^^^^^^^^^^^

A list property specifying diffusion coefficients for each mass
component simulated, corresponding to the **DIFFU** input block in the
TOUGH2 data file. The list has length ``multi['num_components']`` (i.e.
NK in TOUGH2 terminology), and each element is a list of the diffusion
coefficients for each component (with length ``multi['num_phases']``, or
NPH).

----

.. _sec:t2data:echo_extra_precision:

``echo_extra_precision`` property
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A Boolean property (AUTOUGH2 only) governing whether data written to an
auxiliary extra-precision file is also echoed to the main data file. If
``True``, all extra-precision data sections are echoed to the main file.

----

.. _sec:t2data:end_keyword:

``end_keyword`` property
^^^^^^^^^^^^^^^^^^^^^^^^

A string property containing the keyword used in the data file to end
the file. Normally this is 'ENDCY', but 'ENDFI' can also be used.

----

.. _sec:t2data:extra_precision:

``extra_precision`` property
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A list property determining which data sections will be written to an
auxiliary extra-precision file (AUTOUGH2 only). Recent versions of
AUTOUGH2 support an additional data file containing some data written
with extra precision. Possible extra-precision data sections are ROCKS,
ELEME, CONNE, RPCAP and GENER. Typical usage of this extra-precision
data is for automatic model calibration using PEST or similar software,
where calculation of derivatives of model outputs with respect to model
parameters requires higher precision than is possible with the standard
TOUGH2 data file format.

The ``extra_precision`` parameter may be a list containing names of
sections to be written in extra precision (e.g. ['RPCAP', 'GENER']), or
set to ``False`` to disable extra precision (equivalent to []), or to
``True`` to specify that all possible sections should be written in
extra precision.

The :ref:`read() <sec:t2data:read>` method of a ``t2data``
object determines whether extra precision data are available by
searching for an additional file with the same base name as the data
file itself, but with a '.pdat' or '.PDAT' extension (depending on the
case of the main data file name). If no such file exists, then no extra
precision data will be read.

----

.. _sec:t2data:filename:

``filename`` property
^^^^^^^^^^^^^^^^^^^^^

A string property containing the name of the TOUGH2 data file on disk.
(This does not correspond to any parameter in the TOUGH2 data file.)

----

.. _sec:t2data:generator:

``generator`` property
^^^^^^^^^^^^^^^^^^^^^^

A dictionary property containing the generators for the simulation,
accessed by tuples of block name and generator name. Each generator is
an object of type :ref:`t2generator <t2generatorobjects>`.

----

.. _sec:t2data:generatorlist:

``generatorlist`` property
^^^^^^^^^^^^^^^^^^^^^^^^^^

A list property containing the generators for the simulation, accessed
by index.

----

.. _sec:t2data:grid:

``grid`` property
^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 data files; grid

A :ref:`t2grid <t2grids>` object representing the simulation grid,
corresponding to the **ELEME** and **CONNE** input blocks in a TOUGH2
data file.

----

.. _sec:t2data:history_block:

``history_block`` property
^^^^^^^^^^^^^^^^^^^^^^^^^^

A list property containing blocks for which time history output is
required, corresponding to the **FOFT** input block in a TOUGH2 data
file. If the ``t2data`` object contains grid data, the items in this
list are :ref:`t2block <t2blockobjects>` objects; otherwise,
they are block names (i.e. strings).

----

.. _sec:t2data:history_connection:

``history_connection`` property
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A list property containing connections for which time history output is
required, corresponding to the **COFT** input block in a TOUGH2 data
file. If the ``t2data`` object contains grid data, the items in this
list are :ref:`t2connection <t2connectionobjects>` objects;
otherwise, they are tuples of block names (i.e. tuples of strings).

----

.. _sec:t2data:history_generator:

``history_generator`` property
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A list property containing blocks in which generators are defined and
for which time history output is required, corresponding to the **GOFT**
input block in a TOUGH2 data file. If the ``t2data`` object contains
grid data, the items in this list are :ref:`t2block <t2blockobjects>`
objects; otherwise, they are block names (i.e. strings).

----

.. _sec:t2data:incon:

``incon`` property
^^^^^^^^^^^^^^^^^^

A dictionary property representing the initial conditions for the
simulation, accessed by block name, corresponding to the **INCON** input
block in a TOUGH2 data file. The value of each element of the dictionary
is a list consisting of the porosity of the block, followed by a list of
the specified initial primary thermodynamic variables in the block. If
the TOUGH2 NSEQ and NADD values are used, these are stored after the
thermodynamic variables. If they are not used, they can either be set to
``None`` or simply omitted.

For example, to specify porosity 0.1 and initial conditions (101.3E3,
20.0) in block ``'AB105'`` of a ``t2data`` object called ``dat``, set
``dat.incon['AB105'] = [0.1, [101.3e3, 20.0]]``.

To specify these same conditions but with NSEQ = 10 and NADD = 2, set
``dat.incon['AB105'] = [0.1, [101.3e3, 20.0], 10, 2]``.

Porosity can be specified as ``None`` if default porosity (from the
rocktype) is to be used.

----

.. _sec:t2data:indom:

``indom`` property
^^^^^^^^^^^^^^^^^^

A dictionary property representing the initial conditions for the
simulation, accessed by rocktype name, corresponding to the **INDOM**
input block in a TOUGH2 data file. The value of each element of the
dictionary is a list consisting of the specified initial primary
thermodynamic variables for the rocktype.

----

.. _sec:t2data:lineq:

``lineq`` property
^^^^^^^^^^^^^^^^^^

A dictionary property representing linear equation solver options,
corresponding to the **LINEQ** input block in an AUTOUGH2 data file. The
individual keys of this property are given in the :ref:`table <tb:lineq>`
below.

.. container::
   :name: tb:lineq

   .. table:: ``lineq`` property keys

      +-------------------+----------+------------------+------------------+
      | **Key**           | **Type** | **Description**  | **AUTOUGH2       |
      |                   |          |                  | parameter**      |
      +===================+==========+==================+==================+
      | ``epsilon``       | float    | solver tolerance | EPN              |
      +-------------------+----------+------------------+------------------+
      | ``gauss``         | integer  | Gauss            | IGAUSS           |
      |                   |          | elimination      |                  |
      |                   |          | parameter        |                  |
      +-------------------+----------+------------------+------------------+
      | ``max_iterations``| integer  | max. number of   | MAXIT            |
      |                   |          | iterations       |                  |
      +-------------------+----------+------------------+------------------+
      | ``num_orthog``    | integer  | number of        | NORTH            |
      |                   |          | or               |                  |
      |                   |          | thogonalisations |                  |
      +-------------------+----------+------------------+------------------+
      | ``type``          | integer  | type of solver   | ISOLVR           |
      |                   |          | (1 or 2)         |                  |
      +-------------------+----------+------------------+------------------+

----

.. _sec:t2data:meshfilename:

``meshfilename`` property
^^^^^^^^^^^^^^^^^^^^^^^^^

A string property (or tuple of strings) containing the name(s) of files
on disk containing the mesh data. (This does not correspond to any
parameter in the TOUGH2 data file.) Its default value is an empty string
which means mesh data will be read from the main data file.

If ``meshfilename`` is a single (non-empty) string, this is interpreted
as the name of a formatted text file containing 'ELEME' and 'CONNE'
sections specifying the mesh (e.g. the 'MESH' file created by TOUGH2 or
TOUGH2_MP).

If ``meshfilename`` is a tuple of two strings, these are interpreted as
the names of two binary files containing the mesh data, e.g. the 'MESHA'
and 'MESHB' files created by TOUGH2_MP.

----

.. _sec:t2data:meshmaker:

``meshmaker`` property
^^^^^^^^^^^^^^^^^^^^^^

A list property representing mesh generation options, corresponding to
the **MESHM** input block in a TOUGH2 data file. For more detail on the
use of **MESHM** data, consult the TOUGH2 users' guide.

The **MESHM** data may contain multiple sections (e.g. creation of a
rectilinear XYZ grid followed by MINC processing), so the ``meshmaker``
property is structured as a list of two-element tuples, each containing
the type of section (``rz2d``, ``xyz`` or ``minc``) followed by the
section data itself.

The form of the section data varies depending on the section type. For
the ``rz2d`` type it is also structured as a list, as these types may
contain variable numbers of sub-sections. (For example, data for the
``rz2d`` type may contain multiple ``logar`` sub-sections for different
logarithmic radial parts of the mesh.) Each sub-section is again a
two-element tuple, consisting of the sub-section type (a string)
followed by a dictionary containing the data for the sub-section.

Data for the ``xyz`` type are also structured as a list, with the first
element containing the stand-alone ``deg`` parameter (a float), followed
by the other sub-sections, corresponding to the **NX**, **NY** and
**NZ** sub-sections in the TOUGH2 data file. The ``minc`` type does not
have sub-sections so MINC data are not structured as a list but simply a
dictionary.

Possible sub-section types for ``rz2d`` data are ``radii``, ``equid``,
``logar`` and ``layer``, corresponding to their (uppercase) keyword
counterparts in the TOUGH2 data file. Data keys for these types are
given in the :ref:`rz2d data keys <tb:rz2d>` table. Data keys for the
``xyz`` and ``minc`` data are given in :ref:`xyz data keys <tb:xyz>`
and :ref:`minc data keys <tb:minc>` tables.

**Example**: The easiest way to understand how the ``meshmaker``
property works is to read some example input data into a ``t2data``
object and examine the result. The **MESHM** data for the standard
TOUGH2 test problem 'rhbc' ('Production from a geothermal reservoir with
hypersaline brine') is represented as a ``t2data`` ``meshmaker``
property as follows:

::

   [('rz2d',[
    ('radii', {'radii': [5.0]}),
    ('equid', {'dr': 2.0, 'nequ': 1}),
    ('logar', {'rlog': 100.0, 'nlog': 50}),
    ('logar', {'rlog': 1000.0, 'nlog': 20}),
    ('equid', {'dr': 0.0, 'nequ': 1}),
    ('layer', {'layer': [500.0]})
    ])
   ]

.. container::
   :name: tb:rz2d

   .. table:: ``rz2d`` data keys

      +---------------+-----------+----------+-----------------------------+----------------------+
      |**Sub-section**| **Key**   | **Type** | **Description**             | **TOUGH2 parameter** |
      |               |           |          |                             |                      |
      +===============+===========+==========+=============================+======================+
      |**radii**      | ``radii`` | list     | specified mesh radii        | RC                   |
      +---------------+-----------+----------+-----------------------------+----------------------+
      |**equid**      | ``dr``    | float    | radial increment            | DR                   |
      |               +-----------+----------+-----------------------------+----------------------+
      |               | ``nequ``  | integer  | number of equidistant radii | NEQU                 |
      +---------------+-----------+----------+-----------------------------+----------------------+
      |**logar**      | ``dr``    | float    | reference radial increment  | DR                   |
      |               +-----------+----------+-----------------------------+----------------------+
      |               | ``nlog``  | integer  | number of logarithmic radii | NLOG                 |
      |               +-----------+----------+-----------------------------+----------------------+
      |               | ``rlog``  | float    | largest radius              | RLOG                 |
      +---------------+-----------+----------+-----------------------------+----------------------+
      |**layer**      | ``layer`` | list     | layer thicknesses           | H                    |
      +---------------+-----------+----------+-----------------------------+----------------------+

.. container::
   :name: tb:xyz

   .. table:: ``xyz`` data keys

      +-----------+----------+-------------------------------------+----------------------+
      | **Key**   | **Type** | **Description**                     | **TOUGH2 parameter** |
      +===========+==========+=====================================+======================+
      | ``deg``   | float    | angle between y-axis and horizontal | DEG                  |
      +-----------+----------+-------------------------------------+----------------------+
      | ``del``   | float    | constant grid increment             | DEL                  |
      +-----------+----------+-------------------------------------+----------------------+
      | ``deli``  | list     | variable grid increments            | DEL                  |
      +-----------+----------+-------------------------------------+----------------------+
      | ``no``    | integer  | number of grid increments           | DR                   |
      +-----------+----------+-------------------------------------+----------------------+
      | ``ntype`` | string   | axis direction ('NX', 'NY' or 'NZ') | NTYPE                |
      +-----------+----------+-------------------------------------+----------------------+

.. container::
   :name: tb:minc

   .. table:: ``minc`` data keys

      +------------------+----------+------------------+---------------------+
      | **Key**          | **Type** | **Description**  | **TOUGH2 parameter**|
      |                  |          |                  |                     |
      +==================+==========+==================+=====================+
      | ``dual``         | string   | treatment of     | DUAL                |
      |                  |          | global           |                     |
      |                  |          | matrix-matrix    |                     |
      |                  |          | flow             |                     |
      +------------------+----------+------------------+---------------------+
      | ``num_continua`` | integer  | number of        | J                   |
      |                  |          | interacting      |                     |
      |                  |          | continua         |                     |
      +------------------+----------+------------------+---------------------+
      | ``spacing``      | list     | fracture         | PAR                 |
      |                  |          | spacings         |                     |
      +------------------+----------+------------------+---------------------+
      | ``type``         | string   | proximity        | TYPE                |
      |                  |          | function type    |                     |
      +------------------+----------+------------------+---------------------+
      | ``vol``          | list     | volume fractions | VOL                 |
      +------------------+----------+------------------+---------------------+
      | ``where``        | string   | direction of     | WHERE               |
      |                  |          | volume fraction  |                     |
      |                  |          | specification    |                     |
      +------------------+----------+------------------+---------------------+

----

.. _sec:t2data:more_options:

``more_option`` property
^^^^^^^^^^^^^^^^^^^^^^^^

An array property containing additional integer parameter options,
corresponding to the **MOMOP** input block in a TOUGH2 data file (it is
not recognised by AUTOUGH2). Introduced by iTOUGH2, this is an extension
of the ``parameter.option`` property. It is of length 21 and is
populated with zeros by default. Like the ``parameter.option`` property,
values are accessed using 1-based (not zero-based) indices.

----

.. _sec:t2data:multi:

``multi`` property
^^^^^^^^^^^^^^^^^^

A dictionary property selecting the equation of state (EOS) module used
and setting associated parameters, corresponding to the **MULTI** input
block in a TOUGH2 or AUTOUGH2 data file. The individual keys of this
property are given in the :ref:`table <tb:multi>` below.

.. container::
   :name: tb:multi

   .. table:: ``multi`` property keys

      +-----------------------------+----------+------------------+---------------------+
      | **Key**                     | **Type** | **Description**  | **TOUGH2 parameter**|
      |                             |          |                  |                     |
      +=============================+==========+==================+=====================+
      | ``eos``                     | string   | EOS name         | NAMEOS              |
      |                             |          | (AUTOUGH2 only)  |                     |
      +-----------------------------+----------+------------------+---------------------+
      | ``num_components``          | integer  | number of        | NK                  |
      |                             |          | components       |                     |
      +-----------------------------+----------+------------------+---------------------+
      | ``num_equations``           | integer  | number of        | NEQ                 |
      |                             |          | equations        |                     |
      +-----------------------------+----------+------------------+---------------------+
      | ``num_inc``                 | integer  | number of mass   | NKIN                |
      |                             |          | components in    |                     |
      |                             |          | INCON data       |                     |
      |                             |          | (TOUGH2 only)    |                     |
      +-----------------------------+----------+------------------+---------------------+
      | ``num_phases``              | integer  | number of phases | NPH                 |
      +-----------------------------+----------+------------------+---------------------+
      | ``num_secondary_parameters``| integer  | number of        | NB                  |
      |                             |          | secondary        |                     |
      |                             |          | parameters       |                     |
      +-----------------------------+----------+------------------+---------------------+

----

.. _sec:t2data:noversion:

``noversion`` property
^^^^^^^^^^^^^^^^^^^^^^

A Boolean property specifying whether to suppress printing of version
and date information, corresponding to the **NOVER** input block in a
TOUGH2 data file.

----

.. _sec:t2data:num_generators:

``num_generators`` property
^^^^^^^^^^^^^^^^^^^^^^^^^^^

A read-only integer property returning the number of generators.

----

.. _sec:t2data:output_times:

``output_times`` property
^^^^^^^^^^^^^^^^^^^^^^^^^

A dictionary property specifying the times at which model output is
required, corresponding to the **TIMES** input block in a TOUGH2 data
file. The individual keys of this property are given in the
:ref:`table <tb:outputtimes>` below.

.. container::
   :name: tb:outputtimes

   .. table:: ``output_times`` property keys

      +------------------------+---------------+----------------+----------------+
      | **Key**                | **Type**      | **Description**| **TOUGH2       |
      |                        |               |                | parameter**    |
      +========================+===============+================+================+
      | ``max_timestep``       | float         | maximum time   | DELAF          |
      |                        |               | step           |                |
      +------------------------+---------------+----------------+----------------+
      | ``num_times_specified``| integer       | number of      | ITI            |
      |                        |               | times          |                |
      |                        |               | specified      |                |
      +------------------------+---------------+----------------+----------------+
      | ``num_times``          | integer       | total number   | ITE            |
      |                        |               | of times       |                |
      +------------------------+---------------+----------------+----------------+
      | ``time``               | list of float | times at which | TIS            |
      |                        |               | output is      |                |
      |                        |               | required       |                |
      +------------------------+---------------+----------------+----------------+
      | ``time_increment``     | float         | time increment | TINTER         |
      |                        |               | after          |                |
      |                        |               | specified      |                |
      |                        |               | times          |                |
      +------------------------+---------------+----------------+----------------+

----

.. _sec:t2data:parameter:

``parameter`` property
^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 data files; simulation parameters

A dictionary property specifying run-time parameters, corresponding to
the **PARAM** input block in a TOUGH2 data file. The individual keys of
this property are given in the :ref:`table <tb:parameter>` below.

The ``option`` parameter (MOP array in TOUGH2) is an array of 24
integers, and has a 1-based index so that its indices are the same as
those in the TOUGH2 documentation. (In fact it is really zero-based,
like all other Python arrays, but has an extra unused
zero\ :sup:`th` element).

.. container::
   :name: tb:parameter

   .. table:: ``parameter`` property keys

      +-------------------------+----------------+----------------+----------------+
      | **Key**                 | **Type**       | **Description**| **TOUGH2       |
      |                         |                |                | parameter**    |
      +=========================+================+================+================+
      | ``absolute_error``      | float          | absolute       | RE2            |
      |                         |                | convergence    |                |
      |                         |                | tolerance      |                |
      +-------------------------+----------------+----------------+----------------+
      | ``be``                  | float          | enhanced       | BE             |
      |                         |                | vapour         |                |
      |                         |                | diffusion      |                |
      +-------------------------+----------------+----------------+----------------+
      | ``const_timestep``      | float          | time step      | DELTEN         |
      |                         |                | length         |                |
      +-------------------------+----------------+----------------+----------------+
      | ``default_incons``      | list of float  | default        | DEP            |
      |                         |                | initial        |                |
      |                         |                | conditions     |                |
      +-------------------------+----------------+----------------+----------------+
      | ``derivative_increment``| float          | numerical      | DFAC           |
      |                         |                | derivative     |                |
      |                         |                | increment      |                |
      |                         |                | factor         |                |
      +-------------------------+----------------+----------------+----------------+
      | ``diff0``               | float          | diffusive      | DIFF0          |
      |                         |                | vapour flux    |                |
      |                         |                | (AUTOUGH2      |                |
      |                         |                | only)          |                |
      +-------------------------+----------------+----------------+----------------+
      | ``gravity``             | float          | gravitational  | GF             |
      |                         |                | acceleration   |                |
      +-------------------------+----------------+----------------+----------------+
      | ``max_duration``        | integer        | maximum        | MSEC           |
      |                         |                | simulation     |                |
      |                         |                | duration       |                |
      |                         |                | (machine       |                |
      |                         |                | seconds)       |                |
      +-------------------------+----------------+----------------+----------------+
      | ``max_iterations``      | integer        | maximum number | NOITE          |
      |                         |                | of iterations  |                |
      |                         |                | per time step  |                |
      +-------------------------+----------------+----------------+----------------+
      | ``max_timesteps``       | integer        | maximum number | MCYC           |
      |                         |                | of time steps  |                |
      +-------------------------+----------------+----------------+----------------+
      | ``max_timestep``        | float          | maximum time   | DELTMX         |
      |                         |                | step size      |                |
      +-------------------------+----------------+----------------+----------------+
      | ``newton_weight``       | float          | Newton-Raphson | WNR            |
      |                         |                | weighting      |                |
      |                         |                | factor         |                |
      +-------------------------+----------------+----------------+----------------+
      | ``option``              | array(24) of   | simulation     | MOP            |
      |                         | integer        | options        |                |
      +-------------------------+----------------+----------------+----------------+
      | ``pivot``               | float          | pivoting       | U              |
      |                         |                | parameter for  |                |
      |                         |                | linear solver  |                |
      +-------------------------+----------------+----------------+----------------+
      | ``print_block``         | string         | block name for | ELST           |
      |                         |                | short printout |                |
      +-------------------------+----------------+----------------+----------------+
      | ``print_interval``      | integer        | time step      | MCYPR          |
      |                         |                | interval for   |                |
      |                         |                | printing       |                |
      +-------------------------+----------------+----------------+----------------+
      | ``print_level``         | integer        | amount of      | KDATA          |
      |                         |                | printout       |                |
      +-------------------------+----------------+----------------+----------------+
      | ``relative_error``      | float          | relative       | RE1            |
      |                         |                | convergence    |                |
      |                         |                | tolerance      |                |
      +-------------------------+----------------+----------------+----------------+
      | ``scale``               | float          | grid scale     | SCALE          |
      |                         |                | factor         |                |
      +-------------------------+----------------+----------------+----------------+
      | ``texp``                | float          | binary         | TEXP           |
      |                         |                | diffusion      |                |
      |                         |                | temperature    |                |
      |                         |                | parameter      |                |
      +-------------------------+----------------+----------------+----------------+
      | ``timestep_reduction``  | float          | time step      | REDLT          |
      |                         |                | reduction      |                |
      |                         |                | factor         |                |
      +-------------------------+----------------+----------------+----------------+
      | ``timestep``            | list of float  | specified time | DLT            |
      |                         |                | step sizes     |                |
      +-------------------------+----------------+----------------+----------------+
      | ``tstart``              | float          | start time     | TSTART         |
      |                         |                | (seconds)      |                |
      +-------------------------+----------------+----------------+----------------+
      | ``tstop``               | float          | stop time      | TIMAX          |
      +-------------------------+----------------+----------------+----------------+
      | ``upstream_weight``     | float          | upstream       | WUP            |
      |                         |                | weighting      |                |
      |                         |                | factor         |                |
      +-------------------------+----------------+----------------+----------------+

----

.. _sec:t2data:relative_permeability:

``relative_permeability`` property
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A dictionary property specifying the relative permeability function
used, corresponding to the first line of the **RPCAP** input block in
the TOUGH2 data file. The individual keys of this property are given in
the :ref:`table <tb:relativepermeability>` below.

.. container::
   :name: tb:relativepermeability

   .. table:: ``relative_permeability`` property keys

      +----------------+----------------+----------------+----------------+
      | **Key**        | **Type**       | **Description**| **TOUGH2       |
      |                |                |                | parameter**    |
      +================+================+================+================+
      | ``parameters`` | array (7) of   | function       | RP             |
      |                | float          | parameters     |                |
      +----------------+----------------+----------------+----------------+
      | ``type``       | integer        | type of        | IRP            |
      |                |                | relative       |                |
      |                |                | permeability   |                |
      |                |                | function       |                |
      +----------------+----------------+----------------+----------------+

----

.. _sec:t2data:selection:

``selection`` property
^^^^^^^^^^^^^^^^^^^^^^

A dictionary property representing selection parameters for the
simulation (only used by some EOS modules, e.g. EOS7, EOS7R, EWASG),
corresponding to the **SELEC** block in the TOUGH2 data file.

The dictionary contains two keys: 'integer' and 'float', the first of
which accesses a list of the integer selection parameters (the first
line of the **SELEC** block), while the second accesses a list of the
float selection parameters (the remaining lines of the **SELEC** block).

----

.. _sec:t2data:short_output:

``short_output`` property
^^^^^^^^^^^^^^^^^^^^^^^^^

A dictionary property representing blocks, connections and generators
for which short output is required, corresponding to the **SHORT** input
block in an AUTOUGH2 data file.

The dictionary contains four keys: 'frequency', 'block', 'connection'
and 'generator'. The last three of these access lists of blocks,
connections and generators respectively for short output. (Note that
each of these lists contains :ref:`t2block <t2blockobjects>`,
:ref:`t2connection <t2connectionobjects>` or
:ref:`t2generator <t2generatorobjects>` objects, rather than
names.) The 'frequency' key accesses the time step frequency (an
integer) for which short output is required.

----

.. _sec:t2data:simulator:

``simulator`` property
^^^^^^^^^^^^^^^^^^^^^^

A string property specifying the type of simulator, corresponding to the
**SIMUL** input block in an AUTOUGH2 data file.

----

.. _sec:t2data:solver:

``solver`` property
^^^^^^^^^^^^^^^^^^^

A dictionary property representing linear equation solver options,
corresponding to the **SOLVR** input block in a TOUGH2 data file. The
individual keys of this property are given in the
:ref:`table <tb:solver>` below.

.. container::
   :name: tb:solver

   .. table::  ``solver`` property keys

      +----------------------------+----------+------------------+------------------+
      | **Key**                    | **Type** | **Description**  | **TOUGH2         |
      |                            |          |                  | parameter**      |
      +============================+==========+==================+==================+
      | ``closure``                | float    | convergence      | CLOSUR           |
      |                            |          | criterion        |                  |
      +----------------------------+----------+------------------+------------------+
      | ``relative_max_iterations``| float    | relative max.    | RITMAX           |
      |                            |          | number of        |                  |
      |                            |          | iterations       |                  |
      +----------------------------+----------+------------------+------------------+
      | ``type``                   | integer  | solver type      | MATSLV           |
      +----------------------------+----------+------------------+------------------+
      | ``o_precond``              | string   | O                | OPROCS           |
      |                            |          | -preconditioning |                  |
      |                            |          | type             |                  |
      +----------------------------+----------+------------------+------------------+
      | ``z_precond``              | string   | Z                | ZPROCS           |
      |                            |          | -preconditioning |                  |
      |                            |          | type             |                  |
      +----------------------------+----------+------------------+------------------+

----

.. _sec:t2data:start:

``start`` property
^^^^^^^^^^^^^^^^^^

A Boolean property specifying whether the flexible start option is used,
corresponding to the **START** input block in a TOUGH2 data file.

----

.. _sec:t2data:title:

``title`` property
^^^^^^^^^^^^^^^^^^

A string property containing the simulation title, corresponding to the
**TITLE** input block in a TOUGH2 data file.

----

.. _sec:t2data:type:

``type`` property
^^^^^^^^^^^^^^^^^

A string property specifying the simulator type ('AUTOUGH2' or
'TOUGH2'). Changing the value of this property will cause one of the
:ref:`convert_to_TOUGH2() <sec:t2data:convert_to_TOUGH2>` or
:ref:`convert_to_AUTOUGH2() <sec:t2data:convert_to_AUTOUGH2>`
methods to be executed, with default method parameters. Hence, changing
the ``type`` property to 'AUTOUGH2' causes the EOS to be set to the
default 'EW'. It is also not possible to specify TOUGH2_MP options when
setting ``type``. For more control over how the conversion is carried
out, use the conversion methods directly instead of setting ``type``.

Functions for reading data from file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to specify customized functions to control how data are
read from a TOUGH2 data file. This is done using the optional
``read_function`` parameter when a ``t2data`` object is created- in
exactly the same way it is done for a ``mulgrid`` object. For more
details, see the corresponding
:ref:`documentation <mulgridreadfunctions>` for ``mulgrid`` objects.
By default, the read functions for ``t2data`` objects are given by the
``default_read_function`` dictionary.

Methods
~~~~~~~

The main methods of a ``t2data`` object are listed in the 
:ref:`table <tb:t2data_methods>` below.

.. container::
   :name: tb:t2data_methods

   .. table:: Methods of a ``t2data`` object

      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      | **Method**                                                         | **Type**                            | **Description**      |
      +====================================================================+=====================================+======================+
      | :ref:`add_generator <sec:t2data:add_generator>`                    | –                                   | adds a generator     |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      | :ref:`clear_generators <sec:t2data:clear_generators>`              | –                                   | deletes all          |
      |                                                                    |                                     | generators           |
      |                                                                    |                                     |                      |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      |    :ref:`convert_to_AUTOUGH2 <sec:t2data:convert_to_AUTOUGH2>`     | –                                   | converts from TOUGH2 |
      |                                                                    |                                     | input to AUTOUGH2    |
      |                                                                    |                                     |                      |
      |                                                                    |                                     |                      |
      |                                                                    |                                     |                      |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      | :ref:`convert_to_TOUGH2 <sec:t2data:convert_to_TOUGH2>`            | –                                   | converts from        |
      |                                                                    |                                     | AUTOUGH2 input to    |
      |                                                                    |                                     | TOUGH2               |
      |                                                                    |                                     |                      |
      |                                                                    |                                     |                      |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      | :ref:`delete_generator <sec:t2data:delete_generator>`              | –                                   | deletes a generator  |
      |                                                                    |                                     |                      |
      |                                                                    |                                     |                      |
      |                                                                    |                                     |                      |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      | :ref:`delete_orphan_generators <sec:t2data:del_orphan_geners>`     | –                                   | deletes orphaned     |
      |                                                                    |                                     | generators           |
      |                                                                    |                                     |                      |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      | :ref:`effective_incons <sec:t2data:effective_incons>`              | list or                             | effective initial    |
      |                                                                    | :ref:`t2incon <incons>`             | conditions           |
      |                                                                    |                                     |                      |
      |                                                                    |                                     |                      |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      | :ref:`generator_index <sec:t2data:generator_index>`                | integer                             | returns index of     |
      |                                                                    |                                     | generator with       |
      |                                                                    |                                     | specified name and   |
      |                                                                    |                                     | block name           |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      | :ref:`json <sec:t2data:json>`                                      | dictionary                          | Waiwera JSON input   |
      |                                                                    |                                     |                      |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      | :ref:`read <sec:t2data:read>`                                      | :ref:`t2data <datafiles>`           | reads data file from |
      |                                                                    |                                     | disk                 |
      |                                                                    |                                     |                      |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      | :ref:`rename_blocks <sec:t2data:rename_blocks>`                    | –                                   | renames blocks       |
      |                                                                    |                                     |                      |
      |                                                                    |                                     |                      |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      | :ref:`run <sec:t2data:run>`                                        | –                                   | runs a TOUGH2        |
      |                                                                    |                                     | simulation           |
      |                                                                    |                                     |                      |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      |    :ref:`specific_generation <sec:t2data:specific_generation>`     | ``np.array``                        | generation per unit  |
      |                                                                    |                                     | volume in each block |
      |                                                                    |                                     |                      |
      |                                                                    |                                     |                      |
      |                                                                    |                                     |                      |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      | :ref:`total_generation <sec:t2data:total_generation>`              | ``np.array``                        | total generation in  |
      |                                                                    |                                     | each block           |
      |                                                                    |                                     |                      |
      |                                                                    |                                     |                      |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      | :ref:`transfer_from <sec:t2data:transfer_from>`                    | –                                   | transfers data from  |
      |                                                                    |                                     | another              |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+
      | :ref:`write <sec:t2data:write>`                                    | –                                   | writes to data file  |
      |                                                                    |                                     | on disk              |
      +--------------------------------------------------------------------+-------------------------------------+----------------------+

Details of these methods are as follows.

----

.. _sec:t2data:add_generator:

``add_generator(generator)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adds a generator to the data file object.

**Parameters:**

-  | **generator**: :ref:`t2generator <t2generatorobjects>`
   | Generator to be added to the data file object.

----

.. _sec:t2data:convert_to_AUTOUGH2:

``convert_to_AUTOUGH2(warn=True, MP=False, simulator='AUTOUGH2.2', eos='EW')``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 data files; converting

Converts a TOUGH2 (or TOUGH2_MP) data file for use with AUTOUGH2.
Various parameter options are altered to try to make the AUTOUGH2
simulation give similar results to the original TOUGH2 simulation. If
necessary, the ``filename`` property is changed to end in '.dat' (or
'.DAT', depending on the case of the base file name), as required by
AUTOUGH2.

The simulator and EOS name can also be specified, as AUTOUGH2 data files
contain this information in the SIMUL and MULTI sections.

**Parameters:**

-  | **warn**: Boolean
   | If ``True``, warnings will be printed regarding TOUGH2 options used
     in the original data file which are not supported in AUTOUGH2.

-  | **MP**: Boolean
   | if ``True``, treats the original ``t2data`` object as a TOUGH2_MP
     data file, which uses some of the parameters differently (e.g.
     MOP(20)).

-  | **simulator**: string
   | Simulator name, used for the leading part of the AUTOUGH2 SIMUL
     data section. Possible values are 'MULKOM', 'TOUGH2', 'TOUGH2.2',
     'AUTOUGH2' and 'AUTOUGH2.2'.

-  | **eos**: string
   | EOS name, used for the trailing part of the AUTOUGH2 SIMUL data
     section (e.g. 'EW', 'EWC', 'EWA', 'EWAV' etc.)

----

.. _sec:t2data:convert_to_TOUGH2:

``convert_to_TOUGH2(warn=True, MP=False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 data files; converting

Converts an AUTOUGH2 data file for use with TOUGH2 (or compatible
simulators such as TOUGH2_MP). Various parameter options are altered to
try to make the TOUGH2 simulation give similar results to the original
AUTOUGH2 simulation. This particularly affects AUTOUGH2 options related
to backward compatibility with MULKOM. In particular, if these are used
then the heat conductivities in the ROCKS block have to be altered to
give the same results. Data blocks specific to AUTOUGH2 (e.g. SIMULATOR,
LINEQ, and SHORT) are removed, and AUTOUGH2-specific generator types are
converted to their TOUGH2 equivalents if possible, or otherwise deleted.

**Parameters:**

-  | **warn**: Boolean
   | If ``True``, warnings will be printed regarding AUTOUGH2 options
     used in the original data file which are not supported in TOUGH2.

-  | **MP**: Boolean
   | if ``True``, converts to a TOUGH2_MP data file, which treats some
     of the parameters differently (e.g. MOP(20)). The ``filename``
     property is also changed to INFILE, as required by TOUGH2_MP.

----

.. _sec:t2data:clear_generators:

``clear_generators()``
^^^^^^^^^^^^^^^^^^^^^^

Deletes all generators from the data file object.

----

.. _sec:t2data:delete_generator:

``delete_generator(blocksourcenames)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Deletes the generator with the specified block and generator (source)
name, if it exists.

**Parameters:**

-  | **blocksourcenames**: tuple
   | Tuple of block name and generator name (both strings) of the
     generator to be deleted.

----

.. _sec:t2data:del_orphan_geners:

``delete_orphan_generators()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Deletes all generators with block names that are not in the grid.

----

.. _sec:t2data:effective_incons:

``effective_incons(incons = None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns effective initial conditions, based on on the specified initial
conditions in combination with any initial conditions specified in the
``t2data`` object itself – whether as default initial conditions
specified via the :ref:`parameter <sec:t2data:parameter>` 
property, or via the :ref:`incon <sec:t2data:incon>` 
property, or the :ref:`indom <sec:t2data:indom>` property (or
any combination of these).

Any ``indom`` specifications override the defaults in the ``parameter``
property. Values in the ``incon`` property override both the defaults
and values in ``indom``. Finally, values passed into this method via the
``incons`` parameter override any other specifications. Note that any of
these may contain incomplete specifications (i.e. values are not
specified for all blocks in the grid).

If only default homogeneous initial conditions are in effect, then a
list of the primary variables is returned. Otherwise, a :ref:`t2incon <incons>`
object is returned with initial conditions values for every
block.

**Parameters:**

-  | **incons**: ``t2incon`` or ``None``
   | Initial conditions object, usually representing the contents of a
     separate initial conditions file.

----

.. _sec:t2data:generator_index:

``generator_index(blocksourcenames)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns the index (in the ``generatorlist`` list) of the generator with
the specified block and generator name.

**Parameters:**

-  | **blocksourcenames**: tuple
   | Tuple of block name and generator name (both strings) of the
     generator.

----

.. _sec:t2data:json:

``json(geo, mesh_filename, atmos_volume = 1.e25, incons = None, eos = None, bdy_incons = None, mesh_coords = 'xyz')``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 data files; JSON

Returns a JSON dictionary representing the contents of the ``t2data``
object (and associated mesh geometry), suitable for input to the
`Waiwera <http://waiwera.github.io>`_ simulator.

Sources in the Waiwera JSON dictionary are given names based on the
corresponding TOUGH2 generator names. If the TOUGH2 model has no
duplicate generator names, these are used directly for the source
names. If there are duplicate generator names, the block names are
prepended to the generator names to form the source names. If there
are duplicate generator names within the same block, the source names
will have "_1", "_2" etc. appended to them as needed to make them
unique.

**Parameters:**

-  | **geo**: :ref:`mulgrid <mulgrids>`
   | Geometry object. Note that geometric meshes with column surface
     elevations that do not correspond to layer elevations are not
     supported in Waiwera. For meshes of this type, the column surface
     elevations can be "snapped" to layer elevations using the
     :ref:`snap_columns_to_nearest_layers() <sec:mulgrid:snap_columns_to_nearest_layers>`
      method. In that case the
     ``t2grid`` in the ``t2data`` object must be updated so it
     corresponds to the snapped mesh geometry, and other parts of the
     data file updated to reference the new mesh (e.g. using the
     :ref:`transfer_from() <sec:t2data:transfer_from>` 
     method). The geometry's :ref:`block_order <sec:mulgrid:blockordering>`
      property should be set to 'dmplex', particularly if
     it contains mixtures of 3- and 4-sided columns.

-  | **mesh_filename**: string
   | The filename of the mesh file (e.g. ExodusII or GMSH mesh) for the
     Waiwera simulation.

-  | **atmos_volume**: float
   | Maximum block volume for blocks to be considered part of the
     geometric grid. Blocks with volume greater than this value (or
     zero) will be treated as boundary condition (e.g. atmosphere)
     blocks rather than part of the simulation mesh.

-  | **incons**: :ref:`t2incon <incons>`, string, or ``None``
   | Initial conditions for the Waiwera model. If specified as a string,
     this should be the filename of the Waiwera HDF5 output file for
     restarting the simulation from the output of a previous run. If
     ``None`` is specified, then default initial conditions will be
     applied from the ``parameter``
     :ref:`property <sec:t2data:parameter>`.

-  | **eos**: string, integer or ``None``
   | Equation of state used for the simulation. For AUTOUGH2
     simulations, this can generally be set to ``None``, and the EOS
     will be read from the ``t2data`` ``simulator`` or ``multi``
     properties. Otherwise, it can be specified as an integer
     corresponding to the EOS number (1 being pure water, 2 being water
     / CO\ :math:`_2` etc.) or as a string corresponding to the AUTOUGH2
     EOS names (EOS1 being 'EW', EOS2 being 'EWC' etc.). Note that for
     integer values, only EOS modules 1, 2 and 4 are supported. For
     AUTOUGH2 EOS names, these correspond to 'W', 'EW', 'EWC' and
     'EWAV'. The AUTOUGH2 passive tracer EOS modules 'EWT' and 'ETD' are
     also supported (the latter supporting only constant diffusivity,
     i.e. all elements of the ``diffusion`` property must be negative
     and equal).

-  | **bdy_incons**: :ref:`t2incon <incons>`, or ``None``
   | TOUGH2 initial conditions from which boundary conditions are to be
     derived. In many cases this parameter is not needed, because
     boundary conditions are taken from the ``incons`` parameter: if the
     ``incons`` parameter is specified as a ``t2incon`` object, then the
     ``bdy_incons`` parameter can be set to ``None``. If, however,
     ``incons`` is a string or ``None``, then it will not contain
     boundary condition data, in which case boundary conditions can be
     specified by passing a ``t2incon`` object as the ``bdy_incons``
     parameter; otherwise, if this is set to ``None`` then default
     boundary conditions will be applied from the default initial
     conditions in the ``t2data`` ``parameter`` property. Faces on which
     to apply boundary conditions are identified by the presence of
     connections to blocks with either zero or large volume (above the
     volume specified by the ``atmos_volume`` parameter). Note that for
     side boundary conditions (with horizontal connections), the
     boundary blocks must have centres defined, otherwise it is not
     possible to calculate the appropriate normal vector for the
     boundary condition.

-  | **mesh_coords**: string
   | String representing the coordinate system to be used in the Waiwera
     model. 3-D Cartesian meshes are identified as 'xyz'. 2-D Cartesian
     meshes may be identified as either 'xy', 'xz', or 'yz' (depending
     on orientation), while 2-D radial meshes are identified as 'rz'.

----

.. _sec:t2data:read:

``read(filename, meshfilename='')``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 data files; reading

Reads a ``t2data`` object from a TOUGH2 data file on disk. The mesh data
may optionally be read from auxiliary files, if it is not present in the
main data file. (Note that if the main data file does contain mesh
information (the 'ELEME' and 'CONNE' sections), any auxiliary mesh files
will not be read.)

**Parameters:**

-  | **filename**: string
   | Name of the TOUGH2 data file to be read.

-  | **meshfilename**: string or tuple
   | Name of separate mesh file(s) to read, containing element and
     connection data. If empty, then mesh data will be read from the
     main data file. If a non-empty string is given, this is interpreted
     as the name of a formatted text file containing 'ELEME' and 'CONNE'
     data sections (as in the 'MESH' files created by TOUGH2 and
     TOUGH2_MP). If a tuple of two filenames is given, these are
     interpreted as the names of the two binary MESHA and MESHB files
     used by TOUGH2_MP.

Note that it is possible to create a ``t2data`` object and read its
contents in from disk files in one step, e.g.:
``dat = t2data(filename,meshfilename)``.

----

.. _sec:t2data:rename_blocks:

``rename_blocks(blockmap={}, invert=False, fix_blocknames = True)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Renames blocks in the model according to the specified block mapping
dictionary. Any block whose name is a key of the block mapping
dictionary is renamed with the corresponding dictionary value. The
blocks in the :ref:`t2grid <t2grids>` object are renamed using its
own :ref:`rename_blocks() <sec:t2grid:rename_blocks>` method.
Other ``t2data`` properties such as generators, initial conditions and
history specifications are similarly renamed.

**Parameters:**

-  | **blockmap**: dictionary
   | Block mapping dictionary, mapping strings to strings.

-  | **invert**: Boolean
   | Set ``True`` to invert the block mapping dictionary, i.e. to map
     its values to its keys. This can be used, for example, to rename
     the blocks to correspond to a geometry created using the
     :ref:`t2grid <t2grids>` :ref:`rectgeo() <sec:t2grid:rectgeo>`
      method, via the block mapping dictionary also created
     by that method.

-  | **fix_blocknames**: Boolean
   | Set ``True`` (the default) to 'fix' block names in the dictionary,
     using the :ref:`fix_blockname() <sec:mulgrid:fix_blockname>` function.

----

.. _sec:t2data:run:

``run(save_filename='', incon_filename='', simulator='AUTOUGH2_2', silent=False, output_filename='')``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 data files; running

Runs an AUTOUGH2 or TOUGH2 (but not TOUGH2_MP) simulation using the data
file corresponding to a ``t2data`` object. The contents of the
``t2data`` object must first have been written to disk using the
``write`` function. If the file names for the save file or initial
conditions file are not specified, they are constructed by changing the
file extension of the data file name. The name of the TOUGH2 executable
can be specified.

For running TOUGH2 (rather than AUTOUGH2), the name of the TOUGH2
executable must be specified via the ``simulator`` parameter. However,
the ``save_filename`` and ``incon_filename`` parameters do not need to
be specified. Initial conditions will be read from the file INCON and
final results written to SAVE. The listing file name will be the same as
the data file name, but with the extension changed to \*.listing, unless
the ``output_filename`` is specified.

Running TOUGH2_MP is generally done via MPI rather than directly, and
the exact syntax for doing so may vary with different implementations of
MPI (OpenMPI, MPICH2 etc.) It is also necessary to specify the number of
processors to use. However it is still possible to run TOUGH2_MP from a
Python script using a system call, e.g.:

::

   from os import system
   system("mpirun -np 16 t2eos1_mp")

**Parameters:**

-  | **save_filename**: string
   | Name of the save file to be written to disk during the simulation
     (AUTOUGH2 only). Default is 'base.save' where the AUTOUGH2 data
     file name is 'base.dat'.

-  | **incon_filename**: string
   | Name of the initial conditions file for the simulation (AUTOUGH2
     only). Default is 'base.incon' where the AUTOUGH2 data file name is
     'base.dat'.

-  | **simulator**: string
   | Name of the AUTOUGH2 or TOUGH2 executable. Default is 'AUTOUGH2_2'.

-  | **silent**: Boolean
   | Set to ``True`` to suppress output to the display while running
     (default is ``False``).

-  | **output_filename**: string
   | Name of the output listing file for the simulation (TOUGH2 only).
     Default is 'base.listing' where the base name of the TOUGH2 data
     file (without file extension) is 'base'.

----

.. _sec:t2data:specific_generation:

``specific_generation(type='MASS', name='')``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns an ``np.array`` containing the total specific generation rate in
each block (i.e. generation rate per unit volume) for the specified
generator type and name.

**Parameters:**

-  | **type**: string
   | Generation type ('HEAT', 'MASS' etc.) – default is 'MASS'.

-  | **name**: string
   | Regular expression to match generator names (e.g. 'SP...' (or
     '^SP') will match all generators with names beginning with 'SP'.)

----

.. _sec:t2data:transfer_from:

``transfer_from(source, sourcegeo, geo, top_generator=[], bottom_generator=[], sourceinconfilename='', inconfilename='', rename_generators=False, preserve_generation_totals=False)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 data files; transferring

Transfers data from another ``t2data`` object, and its associated
``mulgrid`` object. Parameters, rock types and rock type assignments,
and optionally initial conditions files are transferred. In general the
data for a given block in the geometry is found by identifying the
nearest block in the source geometry and transferring data from that
block. There are, however, exceptions, such as for generators that need
to remain on the surface or bottom of the model. The ``top_generator``
and ``bottom_generator`` lists specify the 'layer' part of the generator
name for generators that should remain on the top or bottom of the
model, respectively.

For generator types in which the ``gx`` and ``rate`` properties
represent generation rates (as opposed to other types for which these
properties are used to represent other things, e.g. productivity index
for wells on deliverability), the values of ``gx`` and ``rate`` are
scaled to account for the different volume of the block the generator
has been mapped into. If ``preserve_generation_totals`` is ``True``, and
a generator with generation rate :math:`G` is mapped into :math:`n`
blocks with volumes :math:`V_1, V_2,\ldots, V_n`, then the generation
rate for the new generator in block :math:`i` will be
:math:`G V_i/\sum_{k=1}^{n}{V_k}`. This should preserve the total
generation rate over the model. (For generator types matching the
``bottom_generator`` or ``top_generator`` specifications, the column
area instead of the block volume is used to determine the appropriate
scaling.) Note that of the columns a top or bottom generator is mapped
into, only those with centres inside the source geometry are included in
the scaling calculations. The generator types for which this scaling is
carried out are: 'AIR', 'COM1', 'COM2', 'COM3', 'COM4', 'COM5', 'HEAT',
'MASS', 'NACL', 'TRAC' and 'VOL'.

If both ``sourceinconfilename`` and ``inconfilename`` are specified, a
new initial conditions file with filename ``inconfilename`` is written
to disk, with initial conditions transferred from the file
``sourceinconfilename``.

**Parameters:**

-  | **source**: :ref:`t2data <datafiles>` 
   | The ``t2data`` object to transfer data from.

-  | **sourcegeo**: :ref:`mulgrid <mulgrids>` 
   | The ``mulgrid`` object corresponding to ``source``.

-  | **geo**: :ref:`mulgrid <mulgrids>` 
   | The ``mulgrid`` object corresponding to the destination ``t2data``
     object.

-  | **top_generator**: list
   | A list of generator 'layer' identifier strings for generators that
     need to be kept at the top of the model (e.g. rain generators).

-  | **bottom_generator**: list
   | A list of generator 'layer' identifier strings for generators that
     need to be kept at the bottom of the model (e.g. basement heat and
     mass inputs).

-  | **sourceinconfilename**: string
   | Name of the (optional) initial conditions file to transfer initial
     conditions data from (corresponding to ``source``).

-  | **inconfilename**: string
   | Name of the (optional) initial conditions file to write,
     corresponding to the destination ``t2data`` object.

-  | **rename_generators**: Boolean
   | If ``False``, generators other than those at the top and bottom of
     the model retain their original names. Otherwise, they will be
     renamed according to their column names in the new grid.

-  | **preserve_generation_totals**: Boolean
   | If ``False`` (the default), the transfer of generators will attempt
     to preserve the distribution of specific generation of the original
     model; otherwise, it will attempt to preserve the total generation
     over the model.

----

.. _sec:t2data:total_generation:

``total_generation(type='MASS', name='')``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns an ``np.array`` containing the total generation rate in each
block for the specified generator type and name.

**Parameters:**

-  | **type**: string
   | Generation type ('HEAT', 'MASS' etc.) – default is 'MASS'.

-  | **name**: string
   | Regular expression to match generator names (e.g. 'SP...' (or
     '^SP') will match all generators with names beginning with 'SP'.)

----

.. _sec:t2data:write:

``write(filename='', meshfilename='', extra_precision=None, echo_extra_precision=None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: TOUGH2 data files; writing

Writes a ``t2data`` object to a TOUGH2 data file on disk. If the
``meshfilename`` parameter is used, mesh information can be written to
auxiliary mesh files.

**Parameters:**

-  | **filename**: string
   | Name of the TOUGH2 data file to be written. If no file name is
     specified, the object's own ``filename`` property is used.

-  | **meshfilename**: string or tuple
   | Name of auxiliary mesh file(s) to be written. If this is empty (the
     default), the object's own ``meshfilename`` property is used.
     Otherwise, if a single (non-empty) string is given, this in
     interpreted as the name of a file to write formatted mesh
     information to (as in the 'MESH' files produced by TOUGH2 and
     TOUGH2_MP). If a tuple of two strings is given, this in interpreted
     as the names of two binary files (as in the 'MESHA' and 'MESHB'
     files produced by TOUGH2_MP).

-  | **extra_precision**: list or Boolean
   | Controls whether to write extra precision data to auxiliary file
     (AUTOUGH2 only). If set to ``True``, then all possible sections
     will be written to the extra precision file. Currently the possible
     extra-precision sections are the ROCKS, ELEME, CONNE, RPCAP and
     GENER sections. If set to ``False`` or [], then no extra-precision
     data will be written. If set to a list of section names (e.g.
     ['RPCAP', 'GENER']), then only those sections will be written in
     extra precision. If set to ``None`` (the default), then the value
     of the data object's ``extra_precision`` property is used.
     Otherwise, the value of this property is overwritten by the value
     specified here.

-  | **echo_extra_precision**: Boolean or None
   | Controls whether to echo all extra-precision data sections to the
     main data file (AUTOUGH2 only). If ``None``, the value of the data
     object's ``echo_extra_precision`` property is used. Otherwise, the
     value of this property is overwritten by the value specified here.

----

.. _t2generatorobjects:

``t2generator`` objects
-----------------------

.. index:: TOUGH2 data files; generators
.. index:: generators

A ``t2generator`` object represents a generator in a TOUGH2 simulation
(i.e. an item in the generation table). The properties of a
``t2generator`` object are given in the
:ref:`table <tb:t2generator_properties>` below. These correspond closely to the
parameters specified in the TOUGH2 **GENER** input block. A
``t2generator`` object has no methods.

.. container::
   :name: tb:t2generator_properties

   .. table:: Properties of a ``t2generator`` object

      +--------------+---------------+--------------------------------+-----------------+
      | **Property** | **Type**      |        **Description**         | **TOUGH2        |
      |              |               |                                | parameter**     |
      +==============+===============+================================+=================+
      | ``block``    | string        | name of block                  | EL, NE          |
      |              |               | containing the                 |                 |
      |              |               | generator                      |                 |
      +--------------+---------------+--------------------------------+-----------------+
      | ``enthalpy`` | list of float | generation enthalpies          | F3              |
      |              |               | (\|ltab\|>1, itab<>'')         |                 |
      |              |               |                                |                 |
      |              |               |                                |                 |
      |              |               |                                |                 |
      |              |               |                                |                 |
      |              |               |                                |                 |
      +--------------+---------------+--------------------------------+-----------------+
      | ``ex``       | float         | enthalpy for                   | EX              |
      |              |               | injection                      |                 |
      +--------------+---------------+--------------------------------+-----------------+
      | ``gx``       | float         | generation rate                | GX              |
      |              |               | (or                            |                 |
      |              |               | productivity                   |                 |
      |              |               | index for                      |                 |
      |              |               | deliverability)                |                 |
      +--------------+---------------+--------------------------------+-----------------+
      | ``hg``       | float         | layer thickness                | HG              |
      |              |               | for                            |                 |
      |              |               | deliverability                 |                 |
      +--------------+---------------+--------------------------------+-----------------+
      | ``fg``       | float         | separator                      | FG              |
      |              |               | pressure/                      |                 |
      |              |               | injectivity                    |                 |
      |              |               | etc.                           |                 |
      +--------------+---------------+--------------------------------+-----------------+
      | ``itab``     | string        | blank unless                   | ITAB            |
      |              |               | table of                       |                 |
      |              |               | specific                       |                 |
      |              |               | enthalpies                     |                 |
      |              |               | specified                      |                 |
      +--------------+---------------+--------------------------------+-----------------+
      | ``ltab``     | integer       | number of                      | LTAB            |
      |              |               | generation                     |                 |
      |              |               | times (or open                 |                 |
      |              |               | layers for                     |                 |
      |              |               | deliverability)                |                 |
      +--------------+---------------+--------------------------------+-----------------+
      | ``nadd``     | integer       | successive                     | NADD            |
      |              |               | block increment                |                 |
      +--------------+---------------+--------------------------------+-----------------+
      | ``nads``     | integer       | successive                     | NADS            |
      |              |               | generator                      |                 |
      |              |               | increment                      |                 |
      +--------------+---------------+--------------------------------+-----------------+
      | ``name``     | string        | generator name                 | SL, NS          |
      +--------------+---------------+--------------------------------+-----------------+
      | ``nseq``     | integer       | number of                      | NSEQ            |
      |              |               | additional                     |                 |
      |              |               | generators                     |                 |
      +--------------+---------------+--------------------------------+-----------------+
      | ``rate``     | list of float |generation rates (\|ltab\|>1)   | F2              |
      |              |               |                                |                 |
      +--------------+---------------+--------------------------------+-----------------+
      | ``time``     | list of float |generation times (\|ltab\|>1)   | F1              |
      |              |               |                                |                 |
      +--------------+---------------+--------------------------------+-----------------+
      | ``type``     | string        |generator type (default 'MASS') | TYPE            |
      |              |               |                                |                 |
      +--------------+---------------+--------------------------------+-----------------+

.. _example-1:

Example
-------

The following piece of Python script opens a MULgraph geometry file and
TOUGH2 data file, changes some TOUGH2 run-time parameters and assigns
heat generators to the blocks in the bottom layer inside a defined area,
with the specified total heat divided uniformly amongst the generators.

::

   geo = mulgrid('gmodel.dat')
   dat = t2data('model.dat')

   dat.parameter['max_timesteps'] = 300
   dat.parameter['print_interval'] = dat.parameter['max_timesteps']/10
   dat.parameter['option'][16] = 5 # time step control

   dat.clear_generators()
   totalheat = 10.e6
   layer = geo.layerlist[-1]  # bottom layer
   cols = [col for col in geo.columnlist if 10.e3 <= col.centre[0] <= 20.e3]
   totalarea = sum([col.area for col in cols])
   q = totalheat / totalarea

   for col in cols:
       blockname = geo.block_name(layer.name, col.name)
       gen = t2generator(name = ' q'+col.name, block = blockname, type = 'HEAT', gx = q*col.area)
       dat.add_generator(gen)

   dat.write()
