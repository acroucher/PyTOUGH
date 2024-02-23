:tocdepth: 3

.. _iapws97:

IAPWS-97 thermodynamics
=======================

.. index:: thermodynamics; IAPWS-97

.. _introduction-7:

Introduction
------------

The ``IAPWS97`` library in PyTOUGH contains a Python implementation of
the main functions of the International Association for the Properties
of Water and Steam (`IAPWS <http://www.iapws.org/>`_) 1997
thermodynamic formulation. These can be used to calculate the
thermodynamic properties of water, steam and supercritical water. The
IAPWS-97 supersedes the :ref:`IFC-67 <t2thermo>` formulation used in
TOUGH2, being generally faster and more accurate, as well as having a
simpler representation of the thermodynamic region around the critical
point.

The operating range of the IAPWS-97 formulation is shown in the
pressure-temperature plot below. It covers temperatures up to 800°C
and pressures up to 100 MPa, and is divided into four thermodynamic
regions:

#. liquid water

#. dry steam

#. supercritical fluid

#. two-phase

The two-phase region (4) follows the saturation line on the
pressure-temperature plot (the boundary between liquid water and dry
steam), up to the critical point :math:`C` (:math:`T` = 373.946 °C,
:math:`P` = 22.064 MPa), where the distinction between liquid water
and steam disappears. Region 3 covers supercritcal fluid (above the
critical point) and also near-critical fluid, just below the critical
point. The boundary between regions 1 and 3 (liquid water and
supercritical) is aribitrarily set at :math:`T` = 350 °C. The boundary
between regions 2 and 3 (dry steam and supercritical) is described by
the ``b23p`` and ``b23t`` :ref:`functions <region23_boundary>`.

.. image:: iapws_regions.*
   :alt: IAPWS-97 thermodynamics operating range
   :width: 500
   :name: fg:iapws97_range

The ``IAPWS97`` library can be imported using the command:

::

      from IAPWS97 import *

The functions available through the ``IAPWS97`` library are listed in
the :ref:`table <tb:iapws97_functions>` below.

.. container::
   :name: tb:iapws97_functions

   .. table:: ``IAPWS97`` functions

      +----------------------------------------------------------------------------+----------+----------------------------+
      | **Function**                                                               | **Type** | **Description**            |
      +============================================================================+==========+============================+
      | :ref:`b23p <sec:iapws97:b23p>`                                             | float    | pressure on boundary       |
      |                                                                            |          | between steam and          |
      |                                                                            |          | supercritical regions, as  |
      |                                                                            |          | a function of temperature  |
      +----------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`b23t <sec:iapws97:b23t>`                                             | float    | temperature on boundary    |
      |                                                                            |          | between steam and          |
      |                                                                            |          | supercritical regions, as  |
      |                                                                            |          | a function of pressure     |
      +----------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`cowat <sec:iapws97:cowat>`                                           | tuple    | density and internal       |
      |                                                                            |          | energy of liquid water     |
      |                                                                            |          |                            |
      +----------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`density_temperature_plot <sec:iapws97:density_temperature_plot>`     | –        | draws region boundaries on |
      |                                                                            |          | a density-temperature plot |
      |                                                                            |          |                            |
      |                                                                            |          |                            |
      +----------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`pressure_temperature_plot<sec:iapws97:pressure_temperature_plot>`    | –        | draws region boundaries on |
      |                                                                            |          | a pressure-temperature     |
      |                                                                            |          | plot                       |
      |                                                                            |          |                            |
      +----------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`region<sec:iapws97:region>`                                          | integer  | thermodynamic region       |
      |                                                                            |          |                            |
      |                                                                            |          |                            |
      +----------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`sat<sec:iapws97:sat>`                                                | float    | saturation pressure as a   |
      |                                                                            |          | function of temperature    |
      |                                                                            |          |                            |
      +----------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`super<sec:iapws97:super>`                                            | tuple    | pressure and internal      |
      |                                                                            |          | energy of supercritical    |
      |                                                                            |          | fluid                      |
      +----------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`supst<sec:iapws97:supst>`                                            | tuple    | density and internal       |
      |                                                                            |          | energy of dry steam        |
      |                                                                            |          |                            |
      +----------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`tsat<sec:iapws97:tsat>`                                              | float    | saturation temperature as  |
      |                                                                            |          | a function of pressure     |
      |                                                                            |          |                            |
      +----------------------------------------------------------------------------+----------+----------------------------+
      | :ref:`visc<sec:iapws97:visc>`                                              | float    | dynamic viscosity of       |
      |                                                                            |          | water, steam or            |
      |                                                                            |          | supercritical fluid        |
      +----------------------------------------------------------------------------+----------+----------------------------+

.. _thermodynamic-functions-1:

Thermodynamic functions
-----------------------

The IAPWS-97 formulation provides thermodynamic functions for liquid
water, dry steam and supercritical fluid. These functions calculate
secondary parameters from the primary thermodynamic variables.

----

.. _sec:iapws97:cowat:

Liquid water: ``cowat(t,p)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``cowat`` function returns a two-element tuple (``d``,\ ``u``) of
density (kg/m\ :math:`^3`) and internal energy (J/kg) of liquid water as
a function of temperature ``t`` (°C) and pressure ``p``
(Pa).

**Parameters:**

-  | **t**: float
   | Temperature (°C)

-  | **p**: float
   | Pressure (Pa)

----

.. _sec:iapws97:supst:

Dry steam: ``supst(t,p)``
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``supst`` function returns a two-element tuple (``d``,\ ``u``) of
density (kg/m\ :math:`^3`) and internal energy (J/kg) of dry steam as a
function of temperature ``t`` (°C) and pressure ``p``
(Pa).

**Parameters:**

-  | **t**: float
   | Temperature (°C)

-  | **p**: float
   | Pressure (Pa)

----

.. _sec:iapws97:super:

Supercritical fluid: ``super(d,t)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``super`` function returns a two-element tuple (``p``,\ ``u``) of
pressure (Pa) and internal energy (J/kg) of supercritical fluid as a
function of density ``d`` (kg/m\ :math:`^3`) and temperature ``t``
(°C).

**Parameters:**

-  | **d**: float
   | Density (kg/m\ :math:`^3`)

-  | **t**: float
   | Temperature (°C)

----

.. _sec:iapws97:visc:

Viscosity: ``visc(d,t)``
------------------------

The ``visc`` function returns the dynamic viscosity (Pa.s) of liquid
water, dry steam or supercritical fluid as a function of density ``d``
(kg/m\ :math:`^3`) and temperature ``t`` (°C). This function is based
on the supplementary "IAPWS Formulation 2008 for the Viscosity of
Ordinary Water Substance", without the critical enhancement of
viscosity near the critical point.

**Parameters:**

-  | **d**: float
   | Density (kg/m\ :math:`^3`)

-  | **t**: float
   | Temperature (°C)

----

Region boundaries
-----------------

These functions describe the boundaries between the four thermodynamic
:ref:`regions <fg:iapws97_range>` of the IAPWS-97 formulation. There
is no equation for the boundary between regions 1 and 3 as this is
simply the line :math:`T` = 350 °C.

.. _saturation-line-satt-and-tsatp-1:

Saturation line: ``sat(t)`` and ``tsat(p)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _sec:iapws97:sat:

``sat(t)``
^^^^^^^^^^

The ``sat`` function returns the saturation pressure (Pa) at a given
temperature ``t`` (°C), for temperatures below the
critical temperature.

**Parameters:**

-  | **t**: float
   | Temperature (°C)

----

.. _sec:iapws97:tsat:

``tsat(p)``
^^^^^^^^^^^

The ``tsat`` function returns the saturation temperature
(°C) at a given pressure ``p`` (Pa), for pressures below
the critical pressure.

**Parameters:**

-  | **p**: float
   | Pressure (Pa)

----

.. _region23_boundary:

Steam/supercritical boundary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _sec:iapws97:b23p:

``b23p(t)``
^^^^^^^^^^^

The ``b23p`` function returns the pressure (Pa) on the boundary of the
dry steam and supercritical regions (regions 2 and 3) at a given
temperature ``t`` (°C).

**Parameters:**

-  | **t**: float
   | Temperature (°C)

----

.. _sec:iapws97:b23t:

``b23t(p)``
^^^^^^^^^^^

The ``b23t`` function returns the temperature (°C) on
the boundary of the dry steam and supercritical regions (regions 2 and
3) at a given pressure ``p`` (Pa).

**Parameters:**

-  | **p**: float
   | Pressure (Pa)

----

.. _determining-thermodynamic-region-1:

Determining thermodynamic region
--------------------------------

.. _sec:iapws97:region:

``region(t, p)``
~~~~~~~~~~~~~~~~

Returns the thermodynamic region (integer, or ``None``) corresponding to
the given temperature (°C) and pressure (Pa), as defined
by the IAPWS-97 specification. The regions are:

#. liquid water

#. dry steam

#. supercritical

If the input temperature and/or pressure are outside the operating range
of the IAPWS-97 formulation, the routine will return ``None``.

**Parameters:**

-  | **t**: float
   | Temperature (°C)

-  | **Pressure**: float
   | Pressure (Pa)

----

Plotting functions
------------------

The ``IAPWS97`` library contains two functions used for including the
IAPWS-97 thermodynamic region boundaries on plots.

----

.. _sec:iapws97:pressure_temperature_plot:

``pressure_temperature_plot(plt)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Draws the IAPWS-97 thermodynamic region boundaries on a
pressure-temperature diagram.

**Parameters:**

-  | **plt**: ``matplotlib.pyplot`` instance
   | An instance of the ``matplotlib.pyplot`` library, imported in the
     calling script using e.g. ``import matplotlib.pyplot as plt``.

----

.. _sec:iapws97:density_temperature_plot:

``density_temperature_plot(plt)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Draws the IAPWS-97 thermodynamic region boundaries on a
density-temperature diagram. (This function requires the Scientific
Python (``scipy``) library to be installed.)

**Parameters:**

-  | **plt**: ``matplotlib.pyplot`` instance
   | An instance of the ``matplotlib.pyplot`` library, imported in the
     calling script using e.g. ``import matplotlib.pyplot as plt``.
