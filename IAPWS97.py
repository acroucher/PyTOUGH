"""IAPWS-97 industrial thermodynamic formulation, as described by:

Wagner, W., Cooper, J.R., Dittman, A., Kijima, J., Kretzschmar, H.-J., Kruse, A., Mares, R.,
Oguchi, K., Sato, H., Stocker, I., Sifner, O., Takaishi, Y., Tanishita, I., Trubenbach, J. and
Willkommen, Th. (2000).  The IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
Water and Steam.  Trans. ASME 122, 150-182.

Viscosity function visc() is described by:
IAPWS, 2008.  Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance.

Copyright 2011 University of Auckland.

This file is part of PyTOUGH.

PyTOUGH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PyTOUGH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PyTOUGH.  If not, see <http://www.gnu.org/licenses/>."""

try:
    import numpy as np
    from numpy import float64
except ImportError: # try importing Numeric on old installs
    import Numeric as np
    from Numeric import Float64 as float64
from math import sqrt,exp

# General constants:
rconst = 0.461526e3       # Gas constant
tc_k = 273.15             # Conversion from Celsius to Kelvin
tcriticalk = 647.096      # Critical temperature (Kelvin)
tcritical = tcriticalk - tc_k
dcritical = 322.0         # Critical density (kg/m3)
pcritical = 22.064e6      # Critical pressure (Pa)

# -- Region 1 constants: --------------------------------------------------

# scaling values for pressure and temperature:
pstar1, tstar1 = 16.53e6, 1386.0

# coefficients n:
nr1 = np.array([
        0.14632971213167, -0.84548187169114, -0.37563603672040e1,
        0.33855169168385e1, -0.95791963387872, 0.15772038513228,
        -0.16616417199501e-1,0.81214629983568e-3, 0.28319080123804e-3,
        -0.60706301565874e-3, -0.18990068218419e-1,-0.32529748770505e-1,
        -0.21841717175414e-1, -0.52838357969930e-4, -0.47184321073267e-3,
        -0.30001780793026e-3, 0.47661393906987e-4, -0.44141845330846e-5,
        -0.72694996297594e-15,-0.31679644845054e-4, -0.28270797985312e-5,
        -0.85205128120103e-9, -0.22425281908000e-5,-0.65171222895601e-6,
        -0.14341729937924e-12, -0.40516996860117e-6, -0.12734301741641e-8,
        -0.17424871230634e-9, -0.68762131295531e-18, 0.14478307828521e-19,
        0.26335781662795e-22,-0.11947622640071e-22, 0.18228094581404e-23,
        -0.93537087292458e-25], float64)

# powers i:
ir1 = np.array([
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
    2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8,
    21, 23, 29, 30, 31, 32])

# powers j:
jr1 = np.array([
    -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1,
    3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8,
    -11, -6, -29, -31, -38, -39, -40, -41])

pc1 = ((2, (1, 1)), (3, (2, 1)), (4, (2, 2)), (5, (3, 2)), (7, (4, 3)), (8, (4, 4)),
       (20, (8, 8, 4)), (21, (20, 1)), (22, (20, 2)), (23, (20, 3)), (28, (23, 5)),
       (29, (22, 7)), (30, (22, 8)), (31, (23, 8)), (32, (28, 4)))

tc1 = ((2, (1, 1)), (3, (2, 1)), (4, (2, 2)), (5, (3, 2)), (6, (3, 3)), (9, (5, 4)),
       (10, (5, 5)), (16, (10, 6)), (17, (16, 1)),
       (-2, (-1, -1)), (-3, (-2, -1)), (-4, (-2, -2)), (-5, (-3, -2)),
       (-6, (-3, -3)), (-7, (-4, -3)), (-8, (-4, -4)), (-9, (-5, -4)),
       (-10, (-5, -5)), (-11, (-6, -5)), (-12, (-6, -6)), (-29, (-10, -10, -9)),
       (-30, (-29, -1)), (-31, (-29, -2)), (-32, (-29, -3)),
       (-38, (-32, -6)), (-39, (-32, -7)), (-40, (-32, -8)),
       (-41, (-32, -9)), (-42, (-39, -3)))

# -- Region 2 constants: --------------------------------------------------

# scaling values for pressure and temperature:
pstar2, tstar2 = 1.0e6, 540.

# coefficients n0:
n0r2 = np.array([
    -0.96927686500217e1, 0.10086655968018e2,-0.56087911283020e-2,
    0.71452738081455e-1, -0.40710498223928, 0.14240819171444e1,
    -0.43839511319450e1,  -0.28408632460772, 0.21268463753307e-1], float64)

# coefficients n:
nr2 = np.array([
    -0.17731742473213e-2, -0.17834862292358e-1,-0.45996013696365e-1, -0.57581259083432e-1,
    -0.50325278727930e-1, -0.33032641670203e-4,-0.18948987516315e-3, -0.39392777243355e-2,
    -0.43797295650573e-1, -0.26674547914087e-4, 0.20481737692309e-7,  0.43870667284435e-6,
    -0.32277677238570e-4, -0.15033924542148e-2, -0.40668253562649e-1, -0.78847309559367e-9,
    0.12790717852285e-7,  0.48225372718507e-6, 0.22922076337661e-5, -0.16714766451061e-10,
    -0.21171472321355e-2, -0.23895741934104e2, -0.59059564324270e-17,-0.12621808899101e-5,
    -0.38946842435739e-1,  0.11256211360459e-10, -0.82311340897998e1,   0.19809712802088e-7,
    0.10406965210174e-18,-0.10234747095929e-12, -0.10018179379511e-8, -0.80882908646985e-10,
    0.10693031879409,  -0.33662250574171, 0.89185845355421e-24, 0.30629316876232e-12,
    -0.42002467698208e-5, -0.59056029685639e-25, 0.37826947613457e-5, -0.12768608934681e-14,
    0.73087610595061e-28, 0.55414715350778e-16, -0.94369707241210e-6], float64)

# powers j0:
j0r2 = np.array([0, 1, -5, -4, -3, -2, -1, 2, 3])

# powers i
ir2 = np.array([
    1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
    4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10,
    10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24])

# powers j:
jr2 = np.array([
    0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6,
    35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36,
    13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53,
    39, 26, 40, 58])

tc2 = ((2, (1, 1)), (-2, (-1, -1)), (-3, (-2, -1)),
       (-4, (-2, -2)), (-5, (-3, -2)), (-6, (-3, -3)))

pc2 = ((2, (1, 1)), (3, (2, 1)), (4, (2, 2)), (5, (3, 2)),
       (6, (3, 3)), (7, (4, 3)), (8, (4, 4)), (9, (5, 4)), 
       (10, (5, 5)), (15, (8, 7)), (16, (8, 8)), (17, (9, 8)),
       (18, (9, 9)), (19, (10, 9)), (20, (10, 10)), 
       (21, (15, 6)), (22, (15, 7)), (23, (15, 8)), (24, (15, 9)))

tsc2 = ((2, (1, 1)), (3, (2, 1)), (4, (2, 2)), (5, (3, 2)),
        (6, (3, 3)), (7, (4, 3)), (8, (4, 4)), (9, (5, 4)), 
        (10, (5, 5)), (11, (6, 5)), (12, (6, 6)), (13, (7, 6)),
        (14, (7, 7)), (15, (8, 7)), (16, (8, 8)), (19, (10, 9)), 
        (20, (10, 10)), (21, (11, 10)), (24, (12, 12)), (25, (13, 12)),
        (26, (13, 13)), (28, (14, 14)), (29, (15, 14)), 
        (34, (19, 15)), (35, (19, 16)), (36, (20, 16)), (38, (19, 19)),
        (39, (20, 19)), (40, (20, 20)), (47, (26, 21)), 
        (48, (24, 24)), (49, (25, 24)), (50, (25, 25)), (52, (26, 26)),
        (53, (28, 25)), (56, (28, 28)), (57, (29, 28)), (58, (29, 29)))

# -- Region 3 constants: --------------------------------------------------

# scaling values for density and temperature:
dstar3, tstar3 = dcritical, tcriticalk

# coefficients n
nr3 = np.array([
    0.10658070028513e1, -0.15732845290239e2,0.20944396974307e2,-0.76867707878716e1,
    0.26185947787954e1, -0.28080781148620e1,0.12053369696517e1,-0.84566812812502e-2,
    -0.12654315477714e1, -0.11524407806681e1, 0.88521043984318,-0.64207765181607,
    0.38493460186671, -0.85214708824206, 0.48972281541877e1,-0.30502617256965e1,
    0.39420536879154e-1, 0.12558408424308, -0.27999329698710,0.13899799569460e1,
    -0.20189915023570e1, -0.82147637173963e-2, -0.47596035734923,0.43984074473500e-1,
    -0.44476435428739, 0.90572070719733, 0.70522450087967,0.10770512626332,
    -0.32913623258954, -0.50871062041158, -0.22175400873096e-1,0.94260751665092e-1,
    0.16436278447961, -0.13503372241348e-1, -0.14834345352472e-1,0.57922953628084e-3,
    0.32308904703711e-2, 0.80964802996215e-4, -0.16557679795037e-3,-0.44923899061815e-4],
               float64)

# powers i:
ir3 = np.array([
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
    3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9,
    9 , 10, 10, 11])

# powers j:
jr3 = np.array([
    0, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22,
    26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2,
    26, 2, 26, 0, 1, 26])

tc3 = ((2, (1, 1)), (3, (2, 1)), (4, (2, 2)), (5, (3, 2)), (6, (3, 3)),
       (7, (4, 3)), (9, (5, 4)), (10, (5, 5)), (11, (6, 5)),
       (12, (6, 6)), (14, (7, 7)), (15, (9, 6)), (16, (9, 7)),
       (17, (10, 7)), (21, (11, 10)), (22, (11, 11)),
       (23, (12, 11)), (25, (14, 11)), (26, (14, 12)))

dc3 = ((2, (1, 1)), (3, (2, 1)), (4, (2, 2)), (5, (3, 2)), (6, (3, 3)),
       (7, (4, 3)), (8, (4, 4)), (9, (5, 4)), (10, (5, 5)), (11, (6, 5)))

# -- Region 4 constants: --------------------------------------------------

pstar4 = 1.0e6  # scaling value for pressure

# coefficients n
nr4 = np.array([
    0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
    0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
    -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
    0.65017534844798e3], float64)

#  -- Boundary between regions 2 & 3: --------------------------------------

nr23 = np.array([
    0.34805185628969e3, -0.11671859879975e1, 0.10192970039326e-2,
    0.57254459862746e3, 0.13918839778870e2],float64)

# -- Constants for dynamic viscosity calculation: -------------------------

mustar = 1.00e-6

h0v = np.array([1.67752, 2.20462, 0.6366564, -0.241605], float64)

h1v = np.array([
    5.20094e-1, 8.50895e-2, -1.08374, -2.89555e-1, 2.22531e-1,
    9.99115e-1, 1.88797, 1.26613, 1.20573e-1, -2.81378e-1,
    -9.06851e-1, -7.72479e-1, -4.89837e-1, -2.57040e-1, 1.61913e-1,
    2.57399e-1, -3.25372e-2, 6.98452e-2, 8.72102e-3, -4.35673e-3,
    -5.93264e-4], float64)

ivs = np.array([0, 1, 2, 3, 0, 1, 2, 3, 5, 0, 1, 2, 3, 4, 0, 1, 0, 3, 4, 3, 5])
jvs = np.array([0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6])

ticv = ((2, (1, 1)), (3, (2, 1)))
tscv = ((2, (1, 1)), (3, (2, 1)), (4, (2, 2)), (5, (3, 2)))
dscv = ((2, (1, 1)), (3, (2, 1)), (4, (2, 2)), (5, (3, 2)), (6, (3, 3)))

#------------------------------------------------------------------------

def power_array(value, combination):
    """Calculates array of powers of 'value', using the specified list of combinations to use
    to calculate them."""

    ppowers = [item[0] for item in combination if item[0] > 0]
    # for simplicity, always calculate -1st power (usually need it anyway):
    npowers = [item[0] for item in combination if item[0] < 0] + [-1]
    nneg, npos = -min(npowers), max(ppowers)
    # store negative powers at end of array, so can use negative indexing:
    p = np.zeros(1 + npos + nneg, float64)
    p[0], p[1], p[-1] = 1.0, value, 1.0 / value
    for c in combination:
        p[c[0]] = p[c[1][0]]
        for mult in c[1][1:]:
            p[c[0]] *= p[mult]
    return p

#------------------------------------------------------------------------

def cowat(t, p):
    """Density d and internal energy u of liquid water as a function of
    temperature t (deg C) and pressure p (Pa)."""

    if t <= 350.0 and p <= 100.e6:

        tk = t + tc_k
        pi = p / pstar1
        tau = tstar1 / tk
        pspow = power_array(7.1 - pi, pc1)
        tspow = power_array(tau - 1.222, tc1)

        gampi = -sum([n * i * pspow[i - 1] * tspow[j] for
                      (i, j, n) in zip(ir1, jr1, nr1)])
        gamt = sum([n * pspow[i] * j * tspow[j - 1] for
                    (i, j, n) in zip(ir1, jr1, nr1)])

        rt = rconst * tk
        d = pstar1 / (rt * gampi)
        u = rt * (tau * gamt - pi * gampi)
        return (d, u)

    else:
        return None

#------------------------------------------------------------------------

def supst(t, p):
    """Density d and internal energy u of dry steam as a function of
    temperature t (deg C) and pressure p (Pa)"""

    if t <= 1000.0 and p <= 100.e6:

        tk = t + tc_k
        pi = p / pstar2
        tau = tstar2 / tk
        taupow = power_array(tau, tc2)
        pipow = power_array(pi, pc2)
        tspow = power_array(tau - 0.5, tsc2)

        gampi0 = 1. / pi
        gamt0 = sum([n * j * taupow[j - 1] for
                     (j, n) in zip(j0r2, n0r2)])
        gampir = sum([n * i * pipow[i - 1] * tspow[j] for
                      (i, j, n) in zip(ir2, jr2, nr2)])
        gamtr = sum([n * pipow[i] * j * tspow[j - 1] for
                     (i, j, n) in zip(ir2, jr2, nr2)])

        gampi = gampi0 + gampir
        rt = rconst * tk

        d = pstar2 / (rt * gampi)
        u = rt * (tau * (gamt0 + gamtr) - pi * gampi)
        return (d, u)

    else:
        return None

#------------------------------------------------------------------------

def super(d, t):
    """Pressure p and internal energy u of supercritical water/steam
    as a function of density d and temperature t (deg C)."""

    tk = t + tc_k
    tau = tstar3 / tk
    delta = d / dstar3
    taupow = power_array(tau, tc3)
    delpow = power_array(delta, dc3)

    phidelta = nr3[0] * delpow[-1] + sum([n * i * delpow[i - 1] * taupow[j] for
                                          (i, j, n) in zip(ir3, jr3, nr3)])
    phitau = sum([n * delpow[i] * j * taupow[j - 1] for
                  (i, j, n) in zip(ir3, jr3, nr3)])

    rt = rconst * tk
    p = d * rt * delta * phidelta
    u = rt * tau * phitau

    return (p, u)

#------------------------------------------------------------------------

def sat(t):
    """Saturation pressure as a function of temperature."""

    if 0. <= t <= tcritical:

        tk = t + tc_k
        theta = tk + nr4[8] / (tk - nr4[9])
        theta2 = theta * theta
        a = theta2 + nr4[0] * theta + nr4[1]
        b = nr4[2] * theta2 + nr4[3] * theta + nr4[4]
        c = nr4[5] * theta2 + nr4[6] * theta + nr4[7]
        x = 2. * c / (-b + sqrt(b * b - 4. * a * c))
        x = x * x
        p = pstar4 * x * x
        return p

    else: return None

#------------------------------------------------------------------------

def tsat(p):
    """Saturation temperature (deg C) as a function of pressure.  Returns
    false if called outside its operating range (611.213 Pa <= p <=
    critical pressure)."""

    if 611.213 <= p <= pcritical:

        beta2 = sqrt(p / pstar4)
        beta = sqrt(beta2)
        e = beta2 + nr4[2] * beta + nr4[5]
        f = nr4[0] * beta2 + nr4[3] * beta + nr4[6]
        g = nr4[1] * beta2 + nr4[4] * beta + nr4[7]
        d = 2.0 * g / (-f - sqrt(f * f - 4. * e * g))
        x = nr4[9] + d
        t = 0.5 * (nr4[9] + d - sqrt(x * x - 4. * (nr4[8] + nr4[9] * d))) - tc_k
        return t

    else: return None

#------------------------------------------------------------------------

def visc(d, t):
    """Calculates dynamic viscosity of water or steam, given the density
    d and temperature t, using the IAPWS industrial formulation 2008.
    Critical enhancement of viscosity near the critical point is not included."""

    tk = t + tc_k
    tau = tk / tcriticalk
    delta = d / dcritical

    tauipow = power_array(1. / tau, ticv)
    tspow = power_array(tauipow[1] - 1., tscv)
    dspow = power_array(delta - 1., dscv)

    # Viscosity in dilute-gas limit:
    s0 = np.dot(h0v, tauipow[0:4])
    mu0 = 100. * sqrt(tau) / s0

    # Contribution due to finite density:
    s1 = sum([tspow[i] * h * dspow[j] for (i, j, h) in zip(ivs, jvs, h1v)])
    mu1 = exp(delta * s1)
    visc = mustar * mu0 * mu1
    return visc

#-----------------------------------------------------------------------

def b23p(t):
    """Returns the pressure on the boundary between regions 2 and 3,
    given a temperature t (deg C)."""

    tk = t + tc_k
    return 1.e6 * (nr23[0] + tk * (nr23[1] + tk * nr23[2]))

#-----------------------------------------------------------------------

def b23t(p):
    """Returns the temperature on the boundary between regions 2 and 3,
    given a pressure p (Pa)."""

    return nr23[3] + sqrt((p / 1.e6 - nr23[4]) / nr23[2]) - tc_k

#------------------------------------------------------------------------

def region(t, p):
    """Returns thermodynamic region corresponding to the given temperature and pressure,
    or None if out of bounds."""

    if (0.01 <= t <= 800.) and (0. <= p <= 100.e6):
        if t <= 350.:
            return 1 if p > sat(t) else 2
        elif t <= 590.:
            return 3 if p > b23p(t) else 2
        else: return 2
    else: return None

#------------------------------------------------------------------------

def pressure_temperature_plot(plt, subplot = 111):
    """Plots IAPWS-97 pressure-temperature region boundaries on plot plt"""

    plt.subplot(subplot)
    plt.xlabel(r'T ($^o$C)')
    plt.ylabel('P (MPa)')
    plt.grid(True)
    MPa = 1.e-6

    # Saturation curve:
    t = np.linspace(0., tcritical, 100)
    p = np.array([sat(tx) for tx in t])
    plt.plot(t, p * MPa, color = 'k', marker = '', linestyle = '-')

    # Critical point:
    plt.plot(np.array([tcritical]), np.array([pcritical]) * MPa,
             color = 'k', marker = 'o', linestyle = '')

    # Region 3 boundaries:
    t0 = 350.
    p0 = sat(t0)
    p = np.linspace(p0, 100.e6, 100)
    t = np.array([b23t(px) for px in p])
    plt.plot(t, p * MPa, color = 'k', marker = '', linestyle = '--')
    plt.plot(np.array([t0, t0]), np.array([p0, 100.e6]) * MPa,
             color = 'k', marker = '', linestyle = '--')

    plt.axis([0, 800, 0, 100])

def density_temperature_plot(plt, subplot = 111):
    """Plots IAPWS-97 density-temperature region boundaries on plot plt"""
    plt.subplot(subplot)

    from scipy.optimize import fsolve
    eps = 1.e-3

    plt.xlabel(r'T ($^o$C)')
    plt.ylabel('$\\rho (kg/m^3)$')
    plt.grid(True)

# Region 4 boundary:
    t = np.linspace(0, 350, 100)
    p = np.array([sat(tx) for tx in t])
    ds = np.array([supst(tx, px)[0] for tx, px in zip(t, p)])
    dw = np.array([cowat(tx, px)[0] for tx, px in zip(t, p)])
    plt.plot(t, ds, color = 'k', marker = '', linestyle = '-')
    plt.plot(t, dw, color = 'k', marker = '', linestyle = '-')
    t = np.linspace(350, tcritical, 50)
    p = np.array([sat(tx) for tx in t])
    def f(dx): return super(dx, tx)[0] - px
    ds = np.array([fsolve(f, 120. + 25. / 200 * (tx - 350.)) for tx, px in zip(t, p)])
    dw = np.array([fsolve(f, 560. - 25. / 200 * (tx - 350.)) for tx, px in zip(t, p)])
    plt.plot(t, ds, color = 'k', marker = '', linestyle = '-')
    plt.plot(t, dw, color = 'k', marker = '', linestyle = '-')

# Critical point:
    plt.plot(np.array([tcritical]), np.array([dcritical]),
             color = 'k', marker = 'o', linestyle = '')

# Region 1/3 boundary:
    t = np.array([350] * 2)
    d = np.array([cowat(350, sat(350))[0], cowat(350, 100.e6)[0]])
    plt.plot(t, d, color = 'k', marker = '', linestyle = '--')

# Region 2/3 boundary:
    tmax = b23t(100.e6) - eps
    t = np.linspace(350, tmax, 100)
    p = np.array([b23p(tx) for tx in t])
    d = np.array([supst(tx, px)[0] for tx, px in zip(t, p)])
    plt.plot(t, d, color = 'k', marker = '', linestyle = '--')

# 100 MPa boundary:
    t = np.linspace(0, 350, 100)
    p = 100.e6
    d = np.array([cowat(tx, p)[0] for tx in t])
    plt.plot(t, d, color = 'k', marker = '', linestyle = '--')
    t = np.linspace(350, 800, 100)
    def g(dx): return super(dx, tx)[0] - p
    d = np.array([fsolve(g, 1100 - tx) for tx in t])
    plt.plot(t, d, color = 'k', marker = '', linestyle = '--')

    plt.axis([0, 800, 0, 1100])
