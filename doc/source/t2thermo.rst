:tocdepth: 3

.. _t2thermo:

TOUGH2 thermodynamics
=====================

.. index:: thermodynamics; IFC-67

.. _introduction-6:

Introduction
------------

The ``t2thermo`` library in PyTOUGH contains a Python implementation
of the thermodynamic routines used in TOUGH2. These can be used to
calculate the thermodynamic properties of water and steam under a
range of conditions. They are based on a subset of the IFC-67
thermodynamic formulation.

The ``t2thermo`` library can be imported using the command:

::

      from t2thermo import *

The functions available through the ``t2thermo`` library are listed in
the :ref:`table <tb:t2thermo_functions>` below.

.. container::
   :name: tb:t2thermo_functions

   .. table:: ``t2thermo`` functions

      +------------------------------------------------------------------------+----------+----------------------------+
      | **Function**                                                           | **Type** | **Description**            |
      +========================================================================+==========+============================+
      | :ref:`cowat<sec:t2thermo:cowat>`                                       | tuple    | density and internal       |
      |                                                                        |          | energy of liquid water     |
      |                                                                        |          |                            |
      +------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`sat<sec:t2thermo:sat>`                                           | float    | saturation pressure as a   |
      |                                                                        |          | function of temperature    |
      |                                                                        |          |                            |
      +------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`region<sec:t2thermo:region>`                                     | integer  | thermodynamic region       |
      |                                                                        |          |                            |
      |                                                                        |          |                            |
      +------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`separated_steam_fraction<sec:t2thermo:separated_steam_fraction>` | float    | separated steam fraction   |
      |                                                                        |          | for given enthalpy and     |
      |                                                                        |          | separator pressure         |
      |                                                                        |          |                            |
      +------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`supst<sec:t2thermo:supst>`                                       | tuple    | density and internal       |
      |                                                                        |          | energy of dry steam        |
      |                                                                        |          |                            |
      +------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`tsat<sec:t2thermo:tsat>`                                         | float    | saturation temperature as  |
      |                                                                        |          | a function of pressure     |
      |                                                                        |          |                            |
      +------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`visw<sec:t2thermo:visw>`                                         | float    | dynamic viscosity of water |
      |                                                                        |          |                            |
      |                                                                        |          |                            |
      +------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`viss<sec:t2thermo:viss>`                                         | float    | dynamic viscosity of steam |
      |                                                                        |          |                            |
      |                                                                        |          |                            |
      +------------------------------------------------------------------------+----------+----------------------------+

Thermodynamic functions
-----------------------

The thermodynamic routines used in TOUGH2 provide functions for liquid
water and dry steam. These functions calculate secondary parameters from
the primary thermodynamic variables. Their names follow the subroutine
names used in the TOUGH2 code.

----

.. _sec:t2thermo:cowat:

Liquid water: ``cowat(t, p, bounds = False)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``cowat`` function returns a two-element tuple (``d``,\ ``u``) of
density (kg/m\ :math:`^3`) and internal energy (J/kg) of liquid water as
a function of temperature ``t`` (°C) and pressure ``p``
(Pa).

**Parameters:**

-  | **t**: float
   | Temperature (°C)

-  | **p**: float
   | Pressure (Pa)

-  | **bounds**: Boolean
   | If ``True``, return ``None`` if the input temperature and pressure
     are outside the operating range of the routine (as defined by
     thermodynamic region 1 of the IFC-67 specification).

----

.. _sec:t2thermo:supst:

Dry steam: ``supst(t, p, bounds = False)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``supst`` function returns a two-element tuple (``d``,\ ``u``) of
density (kg/m\ :math:`^3`) and internal energy (J/kg) of dry steam as a
function of temperature ``t`` (°C) and pressure ``p``
(Pa).

**Parameters:**

-  | **t**: float
   | Temperature (°C)

-  | **p**: float
   | Pressure (Pa)

-  | **bounds**: Boolean
   | If ``True``, return ``None`` if the input temperature and pressure
     are outside the operating range of the routine (as defined by
     thermodynamic region 2 of the IFC-67 specification).

----

Viscosity
---------

.. _sec:t2thermo:visw:

Liquid water: ``visw(t,p,ps)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``visw`` function returns the dynamic viscosity (Pa.s) of liquid
water as a function of temperature ``t`` (°C), pressure (Pa) and
saturation pressure (Pa).

**Parameters:**

-  | **t**: float
   | Temperature (°C)

-  | **p**: float
   | Pressure (Pa)

-  | **ps**: float
   | Saturation pressure (Pa), calculated for example using the ``sat``
     function.

----

.. _sec:t2thermo:viss:

Dry steam: ``viss(t,d)``
~~~~~~~~~~~~~~~~~~~~~~~~

The ``viss`` function returns the dynamic viscosity (Pa.s) of dry steam
as a function of temperature ``t`` (°C) and density
``d`` (kg/m\ :math:`^3`).

**Parameters:**

-  | **t**: float
   | Temperature (°C)

-  | **d**: float
   | Density (kg/m\ :math:`^3`)

----

Saturation line: ``sat(t)`` and ``tsat(p)``
-------------------------------------------

.. _sec:t2thermo:sat:

``sat(t, bounds = False)``
~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``sat`` function returns the saturation pressure (Pa) at a given
temperature ``t`` (°C), for temperatures below the
critical temperature.

**Parameters:**

-  | **t**: float
   | Temperature (°C)

-  | **bounds**: Boolean
   | If ``True``, return ``None`` if the input temperature is outside
     the operating range of the routine (i.e. less than
     0.01 °C or greater than the critical temperature,
     374.15 °C ).

----

.. _sec:t2thermo:tsat:

``tsat(p, bounds = False)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``tsat`` function returns the saturation temperature
(°C) at a given pressure ``p`` (Pa), for pressures below
the critical pressure.

Note that the IFC-67 formulation did not include an explicit formula for
calculating saturation temperature as a function of pressure, so here
(as in TOUGH2) this is calculated using an iterative root-finding
process on the ``sat`` function. The root-finding function is from the
``scipy`` library, so this library must be installed before the ``tsat``
function will work.

**Parameters:**

-  | **p**: float
   | Pressure (Pa)

-  | **bounds**: Boolean
   | If ``True``, return ``None`` if the input pressure is outside the
     operating range of the routine (i.e. less than ``sat(0.01)`` or
     greater than the critical pressure, 22.12 MPa).

----

Other functions
---------------

Separated steam fraction
~~~~~~~~~~~~~~~~~~~~~~~~

.. _sec:t2thermo:separated_steam_fraction:

``separated_steam_fraction(h, separator_pressure, separator_pressure2 = None)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns the separated steam fraction for a given enthalpy ``h`` and
separator pressure. A second separator pressure may be specified in the
case of two-stage flash.

**Parameters:**

-  | **h**: float
   | Enthalpy (J/kg)

-  | **separator_pressure**: float
   | Steam separator pressure (Pa)

-  | **separator_pressure2**: float (or ``None``)
   | Second separator pressure (Pa) for two-stage flash – set to
     ``None`` (the default) for single-stage.

----

Determining thermodynamic region
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _sec:t2thermo:region:

``region(t, p)``
^^^^^^^^^^^^^^^^

Returns the thermodynamic region (integer, or ``None``) corresponding to
the given temperature (°C) and pressure (Pa), as defined
by the IFC-67 specification. The regions are:

#. liquid water

#. dry steam

#. supercritical

#. near-critical

If the input temperature and/or pressure are outside the operating range
of the IFC-67 formulation, the routine will return ``None``.

**Parameters:**

-  | **t**: float
   | Temperature (°C)

-  | **Pressure**: float
   | Pressure (Pa)

----

.. _example-3:

Example
-------

The following script reads in a geometry file and writes an initial
conditions file with approximate hydrostatic conditions corresponding to
a specified vertical temperature gradient. In this case, the model has a
simple flat surface, so that each column has the same number of layers.
The ``cowat`` function is used to calculate the fluid density at each
layer, and hence the approximate vertical pressure distribution.

::

   from mulgrids import *
   from t2thermo import *

   geo = mulgrid('gmodel.dat')

   patm, tatm = 101.325e3, 15.0
   ptblk = np.zeros((geo.num_blocks, 2))
   ptblk[:,0] = patm; ptblk[:,1] = tatm

   g = 9.8
   p, t = patm, tatm
   thick = 0.0
   tgradient = 30 # deg C/km
   for lay in geo.layerlist[1:]:
       d = cowat(t, p)[0]
       thisthick = lay.top - lay.bottom
       h = 0.5 * (thick + thisthick)
       p += d * g * h
       t += tgradient / 1.e3 * h
       thick = thisthick
       for col in geo.columnlist:
           blkname = geo.block_name(lay.name, col.name)
           iblk = geo.block_name_index[blkname]
           ptblk[iblk] = [p, t]
   inc = dat.grid.incons(ptblk)
   inc.write('model.incon')

