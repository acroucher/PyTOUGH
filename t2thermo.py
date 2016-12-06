"""IFC-67 thermodynamic formulation, as used in AUTOUGH2.  There are the 'fast' versions
of the thermodynamic routines originally developed by Mike O'Sullivan.

Copyright 2011 University of Auckland.

This file is part of PyTOUGH.

PyTOUGH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PyTOUGH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PyTOUGH.  If not, see <http://www.gnu.org/licenses/>."""

# data arrays:
cowat_a = [
    0., 6.824687741e3, -5.422063673e2, -2.096666205e4, 3.941286787e4,
    -13.466555478e4, 29.707143084e4, -4.375647096e5, 42.954208335e4,
    -27.067012452e4, 9.926972482e4, -16.138168904e3, 7.982692717e0,
    -2.616571843e-2, 1.522411790e-3, 2.284279054e-2, 2.421647003e2,
    1.269716088e-10, 2.074838328e-7, 2.174020350e-8, 1.105710498e-9,
    1.293441934e1, 1.308119072e-5, 6.047626338e-14]
cowat_sa = [
    0., 8.438375405e-1, 5.362162162e-4, 1.720000000e0, 7.342278489e-2,
    4.975858870e-2, 6.537154300e-1, 1.150e-6, 1.51080e-5,
    1.41880e-1, 7.002753165e0, 2.995284926e-4, 2.040e-1]
supst_b_index = [
    0, 1, 2, 3, 4, 5, 11, 12, 21, 22, 23, 31, 32, 41, 42, 51, 52,
    53, 61, 62, 71, 72, 81, 82, 90, 91, 92, 93, 94, 95, 96]
supst_b_data = [
    16.83599274, 28.56067796, -54.38923329, 0.4330662834, -0.6547711697, 8.565182058e-2,
    6.670375918e-2, 1.388983801, 8.390104328e-2, 2.614670893e-2,
    -3.373439453e-2, 4.520918904e-1, 1.069036614e-1, -5.975336707e-1,
    -8.847535804e-2, 5.958051609e-1, -5.159303373e-1, 2.075021122e-1, 1.190610271e-1,
    -9.867174132e-2, 1.683998803e-1, -5.809438001e-2, 6.552390126e-3,
    5.710218649e-4, 1.936587558e2, -1.388522425e3, 4.126607219e3, -6.508211677e3,
    5.745984054e3, -2.693088365e3, 5.235718623e2]
supst_b = dict(zip(supst_b_index, supst_b_data))
supst_sb = {
    0: 7.633333333e-1, 61: 4.006073948e-1, 71: 8.636081627e-2,
    81: -8.532322921e-1, 82: 3.460208861e-1}

Pc1 = 22120000.
Tc1 = 647.3
L0, L1, L2 = 1.574373327e1, -3.417061978e1, 1.931380707e1
tc_k = 273.15
Tc1_C = Tc1 - tc_k

from math import sqrt, exp

def cowat(t, p, bounds = False):
    """Density d and internal energy u of liquid water as a function of
    temperature t (deg C) and pressure p (Pa).
    If bounds is True, return None values if pressure and/or temperature
    are outside the operating bounds of the routine."""
    if bounds:
        if (0.01 <= t <= 350.) and (p <= 1.e8): ok = (p >= sat(t))
        else: ok = False
    else: ok = True
    if ok:
        TKR = (t + 273.15) / 647.3
        TKR2 = TKR * TKR
        TKR3 = TKR * TKR2
        TKR4 = TKR2 * TKR2
        TKR5 = TKR2 * TKR3
        TKR6 = TKR4 * TKR2
        TKR7 = TKR4 * TKR3
        TKR8 = TKR4 * TKR4
        TKR10 = TKR4 * TKR6
        TKR11 = TKR * TKR10
        TKR19 = TKR8 * TKR11
        TKR18 = TKR8 * TKR10
        TKR20 = TKR10 * TKR10
        PNMR = p / 2.212e7
        PNMR2 = PNMR * PNMR
        PNMR3 = PNMR * PNMR2
        PNMR4 = PNMR * PNMR3
        Y = 1.0 - cowat_sa[1] * TKR2 - cowat_sa[2] / TKR6
        ZP = cowat_sa[3] * Y * Y - 2.0 * cowat_sa[4] * TKR + 2.0 * cowat_sa[5] * PNMR
        if ZP >= 0.0:
            Z = Y + sqrt(ZP)
            CZ = Z ** (5. / 17.)
            PAR1 = cowat_a[12] * cowat_sa[5] / CZ
            CC1 = cowat_sa[6] - TKR
            CC2 = CC1 * CC1
            CC4 = CC2 * CC2
            CC8 = CC4 * CC4
            CC10 = CC2 * CC8
            AA1 = cowat_sa[7] + TKR19
            PAR2 = cowat_a[13] + cowat_a[14] * TKR + cowat_a[15] * TKR2 + \
                   cowat_a[16] * CC10 + cowat_a[17] / AA1
            PAR3 = (cowat_a[18] + 2. * cowat_a[19] * PNMR + 3. * \
                    cowat_a[20] * PNMR2) / (cowat_sa[8] + TKR11)
            DD1 = cowat_sa[10] + PNMR
            DD2 = DD1 * DD1
            DD4 = DD2 * DD2
            PAR4 = cowat_a[21] * TKR18 * (cowat_sa[9] + TKR2) * (-3.0 / DD4 + cowat_sa[11])
            PAR5 = 3.0 * cowat_a[22] * (cowat_sa[12] - TKR) * \
                   PNMR2 + 4.0 * cowat_a[23] / TKR20 * PNMR3
            VMKR = PAR1 + PAR2 - PAR3 - PAR4 + PAR5
            V = VMKR * 3.17e-3
            D = 1.0 / V
            YD = -2.0 * cowat_sa[1] * TKR + 6.0 * cowat_sa[2] / TKR7
            SNUM = cowat_a[10] + cowat_a[11] * TKR
            SNUM = SNUM * TKR + cowat_a[9]
            SNUM = SNUM * TKR + cowat_a[8]
            SNUM = SNUM * TKR + cowat_a[7]
            SNUM = SNUM * TKR + cowat_a[6]
            SNUM = SNUM * TKR + cowat_a[5]
            SNUM = SNUM * TKR + cowat_a[4]
            SNUM = SNUM * TKR2 - cowat_a[2]
            PRT1 = cowat_a[12] * (Z * (17.0 * (Z / 29.0 - Y / 12.0) + \
                                       5.0 * TKR * YD / 12.0) + cowat_sa[4] * \
                                  TKR - (cowat_sa[3] - 1.0) * TKR * Y * YD) / CZ
            PRT2 = PNMR * (cowat_a[13] - cowat_a[15] * TKR2 + cowat_a[16] * \
                           (9.0 * TKR + cowat_sa[6]) * CC8 * CC1 + cowat_a[17] * \
                           (19.0 * TKR19 + AA1) / (AA1 * AA1))
            BB1 = cowat_sa[8] + TKR11
            BB2 = BB1 * BB1
            PRT3 = (11.0 * TKR11 + BB1) / BB2 * (cowat_a[18] * PNMR + cowat_a[19] * \
                                                 PNMR2 + cowat_a[20] * PNMR3)
            EE1 = cowat_sa[10] + PNMR
            EE3 = EE1 * EE1 * EE1
            PRT4 = cowat_a[21] * TKR18 * (17.0 * cowat_sa[9] + 19.0 * TKR2) * \
                   (1.0 / EE3 + cowat_sa[11] * PNMR)
            PRT5 = cowat_a[22] * cowat_sa[12] * PNMR3 + 21.0 * cowat_a[23] / TKR20 * PNMR4
            ENTR = cowat_a[1] * TKR - SNUM + PRT1 + PRT2 - PRT3 + PRT4 + PRT5
            H = ENTR * 70120.4
            U = H - p * V
            return D, U
        else: return None, None
    else: return None, None

def supst(t, p, bounds = False):
    """Density d and internal energy u of dry steam as a function of
    temperature t (deg C) and pressure p (Pa).
    If bounds is True, return None values if pressure and / or temperature
    are outside the operating bounds of the routine."""
    if bounds:
        if (0.01 <= t <= 800.) and (0 <= p):
            if t <= Tc1_C: ok = (p <= sat(t))
            elif t <= 590.: ok = (p <= b23p(t))
            else: ok = (p <= 1.e8)
        else: ok = False
    else: ok = True
    if ok:
        THETA = (t + 273.15) / 647.3
        BETA = p / 2.212e7
        RI1 = 4.260321148
        X = exp(supst_sb[0] * (1.0 - THETA))
        X2 = X * X
        X3 = X2 * X
        X4 = X3 * X
        X5 = X4 * X
        X6 = X5 * X
        X8 = X6 * X2
        X10 = X6 * X4
        X11 = X10 * X
        X14 = X10 * X4
        X18 = X14 * X4
        X19 = X18 * X
        X24 = X18 * X6
        X27 = X24 * X3
        THETA2 = THETA * THETA
        THETA3 = THETA2 * THETA
        THETA4 = THETA3 * THETA
        BETA2 = BETA * BETA
        BETA3 = BETA2 * BETA
        BETA4 = BETA3 * BETA
        BETA5 = BETA4 * BETA
        BETA6 = BETA5 * BETA
        BETA7 = BETA6 * BETA
        BETAL = 15.74373327 - 34.17061978 * THETA + 19.31380707 * THETA2
        DBETAL =  -34.17061978 + 38.62761414 * THETA
        R = BETA / BETAL
        R2 = R * R
        R4 = R2 * R2
        R6 = R4 * R2
        R10 = R6 * R4
        CHI2 = RI1 * THETA / BETA
        SC = (supst_b[11] * X10 + supst_b[12]) * X3
        CHI2 = CHI2 - SC
        SC = supst_b[21] * X18 + supst_b[22] * X2 + supst_b[23] * X
        CHI2 = CHI2 - 2 * BETA * SC
        SC = (supst_b[31] * X8 + supst_b[32]) * X10
        CHI2 = CHI2 - 3 * BETA2 * SC
        SC = (supst_b[41] * X11 + supst_b[42]) * X14
        CHI2 = CHI2 - 4 * BETA3 * SC
        SC = (supst_b[51] * X8 + supst_b[52] * X4 + supst_b[53]) * X24
        CHI2 = CHI2 - 5 * BETA4 * SC
        SD1 = 1.0 / BETA4 + supst_sb[61] * X14
        SD2 = 1.0 / BETA5 + supst_sb[71] * X19
        SD3 = 1.0 / BETA6 + (supst_sb[81] * X27 + supst_sb[82]) * X27
        SD12 = SD1 * SD1
        SD22 = SD2 * SD2
        SD32 = SD3 * SD3
        SN = (supst_b[61] * X + supst_b[62]) * X11
        CHI2 = CHI2 - SN / SD12 * 4 / BETA5
        SN = (supst_b[71] * X6 + supst_b[72]) * X18
        CHI2 = CHI2 - SN / SD22 * 5 / BETA6
        SN = (supst_b[81] * X10 + supst_b[82]) * X14
        CHI2 = CHI2 - SN / SD32 * 6 / BETA7
        SC = supst_b[96]
        SC = SC * X + supst_b[95]
        SC = SC * X + supst_b[94]
        SC = SC * X + supst_b[93]
        SC = SC * X + supst_b[92]
        SC = SC * X + supst_b[91]
        SC = SC * X + supst_b[90]
        CHI2 = CHI2 + 11.0 * R10 * SC
        V = CHI2 * 0.00317
        D = 1.0 / V
        OS1 = supst_sb[0] * THETA
        EPS2 = supst_b[0] * THETA - (-supst_b[1] + supst_b[3] * THETA2 + 2.0 * \
                                     supst_b[4] * THETA3 + 3.0 * supst_b[5] * THETA4)
        SC = (supst_b[11] * (1.0 + 13.0 * OS1) * X10 + supst_b[12] * (1.0 + 3.0 * OS1)) * X3
        EPS2 = EPS2 - BETA * SC
        SC = supst_b[21] * (1.0 + 18.0 * OS1) * X18 + \
             supst_b[22] * (1.0 + 2.0 * OS1) * X2 + supst_b[23] * (1.0 + OS1) * X
        EPS2 = EPS2 - BETA2 * SC
        SC = (supst_b[31] * (1.0 + 18.0 * OS1) * X8 + supst_b[32] * (1.0 + 10.0 * OS1)) * X10
        EPS2 = EPS2 - BETA3 * SC
        SC = (supst_b[41] * (1.0 + 25.0 * OS1) * X11 + supst_b[42] * (1.0 + 14.0 * OS1)) * X14
        EPS2 = EPS2 - BETA4 * SC
        SC = (supst_b[51] * (1.0 + 32.0 * OS1) * X8 + \
              supst_b[52] * (1.0 + 28.0 * OS1) * X4 + supst_b[53] * (1.0 + 24.0 * OS1)) * X24
        EPS2 = EPS2 - BETA5 * SC
        SN6 = 14.0 * supst_sb[61] * X14
        SN7 = 19.0 * supst_sb[71] * X19
        SN8 = (54.0 * supst_sb[81] * X27 + 27.0 * supst_sb[82]) * X27
        OS5 = 1.0 + 11.0 * OS1 - OS1 * SN6 / SD1
        SC = (supst_b[61] * X * (OS1 + OS5) + supst_b[62] * OS5) * (X11 / SD1)
        EPS2 = EPS2 - SC
        OS6 = 1.0 + 24.0 * OS1 - OS1 * SN7 / SD2
        SC = (supst_b[71] * X6 * OS6 + supst_b[72] * (OS6 - 6.0 * OS1)) * (X18 / SD2)
        EPS2 = EPS2 - SC
        OS7 = 1.0 + 24.0 * OS1 - OS1 * SN8 / SD3
        SC = (supst_b[81] * X10 * OS7 + supst_b[82] * (OS7 - 10.0 *  OS1)) * (X14 / SD3)
        EPS2 = EPS2 - SC
        OS2 = 1.0 + THETA * 10.0 * DBETAL / BETAL
        SC = (OS2 + 6.0 * OS1) * supst_b[96]
        SC = SC * X + (OS2 + 5.0 * OS1) * supst_b[95]
        SC = SC * X + (OS2 + 4.0 * OS1) * supst_b[94]
        SC = SC * X + (OS2 + 3.0 * OS1) * supst_b[93]
        SC = SC * X + (OS2 + 2.0 * OS1) * supst_b[92]
        SC = SC * X + (OS2 + OS1) * supst_b[91]
        SC = SC * X + OS2 * supst_b[90]
        EPS2 = EPS2 + BETA * R10 * SC
        H = EPS2 * 70120.4
        U = H - p * V
        return D, U
    else: return None, None

def sat(t, bounds = False):
    """Saturation pressure (Pa) as a function of temperature (deg C).
    If bounds is True, returns None if the temperature is out of 
    range."""
    a = [0., -7.691234564, -2.608023696e1, -1.681706546e2, 6.423285504e1,
       -1.189646225e2, 4.167117320, 2.097506760e1, 1.0e9, 6.0]
    if bounds:
        ok = (0.01 <= t <= Tc1_C)
    else: ok = True
    if ok:
        if (0.01 <= t <= 500.0): # arbitrary upper limit in TOUGH2 implementation
            TC = (t + 273.15) / 647.3
            X1 = 1.0 - TC
            X2 = X1 * X1
            SC = a[5] * X1 + a[4]
            SC = SC * X1 + a[3]
            SC = SC * X1 + a[2]
            SC = SC * X1 + a[1]
            SC = SC * X1
            PC = exp(SC / (TC * (1.0 + a[6] * X1 + a[7] * X2)) - \
                     X1 / (a[8] * X2 + a[9]))
            return PC * 2.212e7
        else: return None
    else: return None

def tsat(p, bounds = False):
    """Saturation temperature (deg C) as a function of pressure (Pa).
    If bounds is True, return None if the pressure is out of range."""
    if bounds:
        ok = (sat(0.01) <= p <= Pc1)
    else: ok = True
    if ok:
        from scipy.optimize import fsolve
        def f(t): return sat(t) - p
        from math import log
        t0 = max(4606.0 / (24.02 - log(p)) - 273.15, 5.0) # starting estimate
        t = fsolve(f, t0)
        # need to check this as some versions of SciPy return an array from fsolve:
        import collections
        if isinstance(t, collections.Iterable): return t[0]
        else: return t
    else: return None

def visw(t, p, ps):
    """Viscosity of liquid water as a function of temperature (deg C) and
    pressure and saturation pressure (Pa)."""
    EX = 247.8 / (t + 133.15)
    PHI = 1.0467 * (t - 31.85)
    AM = 1.0 + PHI * (p - ps) * 1.0e-11
    return 1.0e-7 * AM * 241.4 * 10.0 ** EX
    
def viss(t, d):
    """Viscosity of vapour as a function of temperature (deg C) and
    density (kg / m3)."""
    V1 = 0.407 * t + 80.4
    if t <= 350.0: return 1.0e-7 * (V1 - d * (1858.0 - 5.9 * t) * 1.0e-3)
    else: return 1.0e-7 * (V1 + d * (0.353 + d * (676.5e-6 + d * 102.1e-9)))
    
def separated_steam_fraction(h,  separator_pressure,  separator_pressure2 = None):
    """Return separated steam fraction from given enthalpy h and separator
    pressure.  Specify a second separator pressure for two-stage
    flash.
    """
    def enth(t, p, f):
        d, u = f(t, p)
        return u + p / d
    def hlhs(p):
        ts = tsat(p)
        return enth(ts, p, cowat), enth(ts, p, supst)
    if separator_pressure2 is None:
        # if single stage
        hl1, hs1 = hlhs(separator_pressure)
        hl = 1.0 / (hs1 - hl1)
        hs = -hl1 / (hs1 - hl1)
    else:
        # if two stages
        hl1, hs1 = hlhs(separator_pressure)
        hl2, hs2 = hlhs(separator_pressure2)
        hl = (hs2 - hl1) / ((hs1 - hl1) * (hs2 - hl2))
        hs = (hs1 * (hl1 - hl2) - hl1 * (hs2 - hl2)) / ((hs1 - hl1) * (hs2 - hl2))
    # steam fraction is hl*h+hs
    frac = hl * h + hs
    return max(min(frac, 1.0), 0.0)

def b23p(t):
    """Returns the pressure on the boundary between regions 2 and 3,
    given a temperature t (deg C)."""
    theta = (t + tc_k) / Tc1
    return Pc1 * (L0 + theta * (L1 + L2 * theta))

def region(t, p):
    """Returns thermodynamic region corresponding to the given temperature and pressure,
    or None if out of bounds."""
    if (0.01 <= t <= 800.) and (0. <= p <= 100.e6):
        if t <= 350.:
            return 2 if p < sat(t) else 1
        elif t <= Tc1_C:
            if p <= b23p(t): return 2
            else: return 3 if p < sat(t) else 4
        elif t <= 590.:
            return 2 if p < b23p(t) else 3
        else: return 2
    else: return None
