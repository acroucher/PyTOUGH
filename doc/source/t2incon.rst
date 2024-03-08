:tocdepth: 3

.. _incons:

TOUGH2 initial conditions
=========================

.. index:: TOUGH2 initial conditions

.. _introduction-4:

Introduction
------------

The ``t2incons`` library in PyTOUGH contains classes and routines for
reading, editing and writing TOUGH2 initial conditions and files. It can
be imported using the command:

::

      from t2incons import *

The initial conditions files used by TOUGH2 and AUTOUGH2 have the same
format. PyTOUGH also supports TOUGHREACT initial conditions files, which
have a slightly different format – permeabilities are included for each
block, and timing information at the bottom of the file is formatted
differently.

``t2incon`` objects
-------------------

The ``t2incons`` library defines a ``t2incon`` class, used for
representing TOUGH2 initial conditions.

**Example:**

::

   inc = t2incon()

creates an empty ``t2incon`` object called ``inc``.

::

   inc = t2incon(filename)

creates a ``t2incon`` object called ``inc`` and reads its contents from
file ``filename``.

.. _properties-3:

Properties
~~~~~~~~~~

The main properties of a ``t2incon`` object are listed in the
:ref:`table <tb:t2incon_properties>` below. Once a set of initial conditions is
loaded into a ``t2incon`` object, conditions for individual blocks can
be accessed by block name or index. For example, for a ``t2incon``
object ``inc``, the initial conditions in block 'blockname' are given
simply by ``inc[blockname]``. This returns a :ref:`t2blockincon <t2blockincons>`
object. Similarly, ``inc[i]`` returns the initial conditions at the block with
(zero-based) index ``i``.

Each column in the initial conditions file can be accessed by adding an
integer (zero-based) index after the ``t2blockincon`` object, so for
example:

::

   t = inc['aa 20'][1]

assigns the variable ``t`` the value of the second primary thermodynamic
variable (index 1) in block ``'AA 20'``. Initial conditions can be
edited in a similar way, for example:

::

   inc['aa 20'][0] = p

assigns the value of ``p`` to the first primary variable (usually
pressure) in block ``'AA 20'``. For convenience, initial conditions for
a given block can also be specified as a simple list or tuple of values,
for example:

::

   inc['ab 25'] = (101.3e5,25.0)

sets the initial conditions at block ``'ab 25'`` to the specified
values. This will work even if no initial conditions have been
previously specified for the given block.

An ``np.array`` of the values of the variables at all blocks can be
found from the ``variable`` property. For example:

::

   inc.variable[:,2]

returns an ``np.array`` of the third variable (index 2) in each block.
The ``variable`` property can also be set to a given array. Note,
however, that the whole array must be set, not just part of it. For
example, adding an offset ``P0`` to all pressures (variable 0) in the
initial conditions could be done by:

::

   v = inc.variable
   v[:,0] += P0
   inc.variable = v

The ``porosity`` property may be set to assign values of porosity to all
blocks. The assigned value may be an ``np.array`` with a value for each
block, or a scalar float (in which case the same value is assigned to
all blocks), or ``None`` which assigns the value in each block to
``None``.

Similarly, for TOUGHREACT initial conditions files, the ``permeability``
property can be used to read or assign permeabilities for all blocks.
When assigning this property, the value can be an ``np.array`` of shape
(``num_blocks``, 3), (i.e. a row for each block), or a single
``np.array`` with 3 elements, to be applied to all blocks, a single
scalar float (to assign isotropic permeabilities to all blocks) or
``None`` which assigns ``None`` to all block permeabilities.

The ``timing`` property of a ``t2incon`` object contains the optional
timing information at the end of the file. This is a dictionary property
with keys ``'kcyc'``, ``'iter'``, ``'nm'``, ``'tstart'`` and
``'sumtim'``, corresponding to the values stored on this line.

The ``simulator`` string property is 'TOUGH2' by default, and is set to
'TOUGHREACT' if permeabilities are detected while reading from file.
Setting this property back to 'TOUGH2' will cause the file to be written
out in TOUGH2 format (no permeabilities, and different format for timing
information) if the ``write()`` method is executed.

.. container::
   :name: tb:t2incon_properties

   .. table:: Properties of a ``t2incon`` object

      +-------------------+--------------+---------------------------------+
      | **Property**      | **Type**     | **Description**                 |
      +===================+==============+=================================+
      | ``blocklist``     | list         | ordered list of block names in  |
      |                   |              | the initial conditions file     |
      +-------------------+--------------+---------------------------------+
      | ``num_blocks``    | integer      | number of blocks at which       |
      |                   |              | conditions are specified        |
      +-------------------+--------------+---------------------------------+
      | ``num_variables`` | integer      | number of thermodynamic         |
      |                   |              | variables specified at each     |
      |                   |              | block                           |
      +-------------------+--------------+---------------------------------+
      | ``permeability``  | ``np.array`` | array of permeability values    |
      |                   |              | specified at each block         |
      |                   |              | (TOUGHREACT only)               |
      +-------------------+--------------+---------------------------------+
      | ``porosity``      | ``np.array`` | array of porosity values        |
      |                   |              | specified at each block         |
      +-------------------+--------------+---------------------------------+
      | ``simulator``     | string       | simulator type ('TOUGH2' or     |
      |                   |              | 'TOUGHREACT')                   |
      +-------------------+--------------+---------------------------------+
      | ``timing``        | dictionary   | additional timing information   |
      |                   |              | for restarting                  |
      +-------------------+--------------+---------------------------------+
      | ``variable``      | ``np.array`` | two-dimensional array of        |
      |                   |              | thermodynamic variable values   |
      |                   |              | at each block                   |
      +-------------------+--------------+---------------------------------+

.. _functions-for-reading-data-from-file-1:

Functions for reading data from file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to specify customized functions to control how data are
read from a TOUGH2 initial conditions file. This is done using the
optional ``read_function`` parameter when a ``t2incon`` object is
created- in exactly the same way it is done for a ``mulgrid`` object.
For more details, see the corresponding
:ref:`documentation <mulgridreadfunctions>` for ``mulgrid``
objects. By default, the read functions for ``t2incon`` objects are given
by the ``fortran_read_function`` dictionary.

Specifying the number of primary variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Most common TOUGH2 EOS modules have no more than four primary variables,
in which case the variables for a given block all fit on one line in the
initial conditions file. However, some EOS modules (e.g. EOS7c and
EOS7r) have more than four primary variables. For these, the variables
for a given block are specified over multiple lines in the initial
conditions file.

In this case, it is not possible for PyTOUGH to reliably detect the
number of primary variables, as it does when there are no more than four
variables. Instead, the number of primary variables must be specified
when the ``t2incon`` object is created (or its
:ref:`read() <sec:t2incon:read>`  method is executed). This can
be done by setting the optional integer ``num_variables`` parameter,
which defaults to ``None`` (meaning PyTOUGH will detect the number of
variables). For example:

::

   from t2incons import *
   inc = t2incon('model.incon', num_variables = 6)

opens initial conditions for an EOS using six primary variables.

For writing initial conditions files with more than four primary
variables, no extra parameters need be set, as the data stored in the
``t2incon`` object determines the number of primary variables, and they
will be written out over multiple lines as required.

Checking block names
^^^^^^^^^^^^^^^^^^^^

By default, when a ``t2incon`` object is read from file, the block
names are checked to make sure they are valid TOUGH2 block names (3
characters plus 2 digits). However these checks can be skipped by
setting the optional ``check_blocknames`` parameter to ``False``. For example:

::

   from t2incons import *
   inc = t2incon('model.incon', check_blocknames = False)

.. _methods-1:

Methods
~~~~~~~

The main methods of a ``t2incon`` object are listed in the
:ref:`table <tb:t2incon_methods>` below.

.. container::
   :name: tb:t2incon_methods

   .. table:: Methods of a ``t2incon`` object

      +--------------------------------------------------+----------+----------------------------+
      | **Method**                                       | **Type** | **Description**            |
      +==================================================+==========+============================+
      | :ref:`add_incon <sec:t2incon:add_incon>`         | –        | adds a set of initial      |
      |                                                  |          | conditions for one block   |
      |                                                  |          |                            |
      +--------------------------------------------------+----------+----------------------------+
      | :ref:`delete_incon <sec:t2incon:delete_incon>`   | –        | deletes the initial        |
      |                                                  |          | conditions for one block   |
      |                                                  |          |                            |
      +--------------------------------------------------+----------+----------------------------+
      | :ref:`empty <sec:t2incon:empty>`                 | –        | deletes all initial        |
      |                                                  |          | conditions from the object |
      |                                                  |          |                            |
      |                                                  |          |                            |
      +--------------------------------------------------+----------+----------------------------+
      | :ref:`insert_incon <sec:t2incon:insert_incon>`   | –        | inserts initial conditions |
      |                                                  |          | for one block at a         |
      |                                                  |          | specified index            |
      +--------------------------------------------------+----------+----------------------------+
      | :ref:`read <sec:t2incon:read>`                   | –        | reads initial conditions   |
      |                                                  |          | from file                  |
      |                                                  |          |                            |
      +--------------------------------------------------+----------+----------------------------+
      | :ref:`transfer_from <sec:t2incon:transfer_from>` | –        | transfers initial          |
      |                                                  |          | conditions from one grid   |
      |                                                  |          | to another                 |
      +--------------------------------------------------+----------+----------------------------+
      | :ref:`write <sec:t2incon:write>`                 | –        | writes initial conditions  |
      |                                                  |          | to file                    |
      |                                                  |          |                            |
      +--------------------------------------------------+----------+----------------------------+

Details of these methods are as follows.

----

.. _sec:t2incon:add_incon:

``add_incon(incon)``
^^^^^^^^^^^^^^^^^^^^

Adds a set of initial conditions for a single block.

**Parameters:**

-  | **incon**: :ref:`t2blockincon <t2blockincons>` 
   | Initial conditions for the block.

----

.. _sec:t2incon:delete_incon:

``delete_incon(blockname)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Deletes a set of initial conditions for a single block.

**Parameters:**

-  | **blockname**: string
   | Name of the block at which initial conditions are to be deleted.

----

.. _sec:t2incon:empty:

``empty()``
^^^^^^^^^^^

Deletes initial conditions for all blocks.

----

.. _sec:t2incon:insert_incon:

``insert_incon(index,incon)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Inserts a set of initial conditions for a single block at the specified
index.

**Parameters:**

-  | **index**: integer
   | Index (zero-based) at which to insert the initial conditions.

-  | **incon**: :ref:`t2blockincon <t2blockincons>`
   | Initial conditions for the block.

----

.. _sec:t2incon:read:

``read(filename, num_variables = None, check_blocknames = True)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reads initial conditions from file.

**Parameters:**

-  | **filename**: string
   | Name of the initial conditions file to be read.

-  | **num_variables**: integer or ``None``
   | If reading initial conditions files with more than four primary
     variables, set to the number of primary variables. Otherwise, the
     default ``None`` value can be used, in which case the number of
     primary variables will be detected automatically.

-  | **check_blocknames**: Boolean
   | Whether to check if block names in the file are valid TOUGH2 block
     names (3 characters followed by 2 digits).

----

.. _sec:t2incon:transfer_from:

``transfer_from(sourceinc, sourcegeo, geo, mapping={}, colmapping={})``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Transfers initial conditions from another ``t2incon`` object
``sourceinc``, using the two corresponding ``mulgrid`` geometry objects
``sourcegeo`` and ``geo``, and optionally the block and column mappings
between the two grids (which are created if not specified).

**Parameters:**

-  | **sourceinc**: :ref:`t2incon <incons>`
   | Source initial conditions object.

-  | **sourcegeo**: :ref:`mulgrid <mulgrids>`
   | Geometry object corresponding to the source initial conditions.

-  | **geo**: :ref:`mulgrid <mulgrids>`
   | Geometry object for the grid to be transferred to.

-  | **mapping**: dictionary
   | Dictionary mapping block names from ``geo`` to ``sourcegeo``.

-  | **colmapping**: dictionary
   | Dictionary mapping column names from ``geo`` to ``sourcegeo``.

----

.. _sec:t2incon:write:

``write(filename, reset=True)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Writes initial conditions to file.

**Parameters:**

-  | **filename**: string
   | Name of the initial conditions file to be written.

-  | **reset**: Boolean
   | Set to ``False`` if timing information is not to be reset - e.g. if
     restarting a transient simulation.

----

.. _t2blockincons:

``t2blockincon`` objects
------------------------

A ``t2blockincon`` object represents the initial conditions for a
particular block. The properties of a ``t2blockincon`` object are given
in the :ref:`table <tb:t2blockincon_properties>` below. The ``permeability``
property is used only by TOUGHREACT. If no values are specified for
``porosity``, ``permeability``, ``nseq`` or ``nadd``, their values are
``None``. A ``t2blockincon`` object has no methods.

The ``variable`` property of a ``t2blockincon`` can be more easily
accessed simply by adding the required (zero-based) variable index after
the object. For example, for a ``t2blockincon`` object ``b``, the value
of the second variable is given simply by ``b[1]``.

To create a new ``t2blockincon`` object, simply invoke the class name
with values of the desired properties, e.g.:

::

     binc = t2blockincon(block = 'abc10', porosity = 0.1, variable = [101.3e3, 28.])

.. container::
   :name: tb:t2blockincon_properties

   .. table:: Properties of a ``t2blockincon`` object

      +------------------+------------------------+------------------------+
      | **Property**     | **Type**               | **Description**        |
      +==================+========================+========================+
      | ``block``        | string                 | block name             |
      +------------------+------------------------+------------------------+
      | ``nadd``         | integer or ``None``    | optional block index   |
      |                  |                        | increment between      |
      |                  |                        | additional blocks with |
      |                  |                        | the same initial       |
      |                  |                        | conditions             |
      +------------------+------------------------+------------------------+
      | ``nseq``         | integer or ``None``    | optional number of     |
      |                  |                        | additional blocks with |
      |                  |                        | the same initial       |
      |                  |                        | conditions             |
      +------------------+------------------------+------------------------+
      | ``permeability`` | ``np.array`` or        | optional permeability  |
      |                  | ``None``               | for the block          |
      |                  |                        | (TOUGHREACT only)      |
      +------------------+------------------------+------------------------+
      | ``porosity``     | float or ``None``      | optional porosity for  |
      |                  |                        | the block              |
      +------------------+------------------------+------------------------+
      | ``variable``     | list                   | list of thermodynamic  |
      |                  |                        | variable values for    |
      |                  |                        | the block              |
      +------------------+------------------------+------------------------+

Reading save files and converting to initial conditions
-------------------------------------------------------

TOUGH2 writes a save file (SAVE, or \*.save for AUTOUGH2) at the end of
the simulation, which has a format almost the same as that of an initial
conditions file and can be used to start a subsequent run. A save file
generally has some extra timing information at the end which can be used
to restart a simulation at a particular time. However, in many cases,
e.g when running natural state simulations, we want to restart at the
original start time and this timing information must be discarded.

PyTOUGH will read a save file into a ``t2incon`` object. This can then
be written to file, providing a simple way to convert save files into
incon files. By default, the timing information is discarded when
writing (it can be retained by setting the ``reset`` parameter of the
``write`` method to ``False``). For example:

::

   t2incon('model1.save').write('model2.incon')

will read the save file ``'model1.save'``, convert it to initial
conditions, and write it to the initial conditions file
``'model2.incon'``.

.. _example-2:

Example
-------

The following piece of Python script reads in a save file and prints out
a table of block names and temperatures for the first 10 blocks. It then
adds an extra variable to each initial condition and gives it a constant
value (giving a new column in the initial conditions file), and finally
writes out the edited initial conditions to a new file.

Adding a new variable to each initial condition can be useful when e.g.
changing from one TOUGH2 equation of state (EOS) module to another, as
different EOS modules may have different numbers of primary
thermodynamic variables.

::

   from t2incons import *
   inc = t2incon('model1.save')
   for blk in inc[0:10]:
       print('Block %5s: temperature = %5.1f' % (blk.block,blk[1]))
   patm = 101.3e3
   for blk in inc: blk.variable.append(patm)
   inc.write('model2.incon')

