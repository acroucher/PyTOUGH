import unittest
from IAPWS97 import *

class IAPWS97TestCase(unittest.TestCase):

    def test_cowat(self):
        """Tests cowat() method on values from IAPWS97 paper."""
        cases = (((300.,3.e6), (0.100215168e-2,0.112324818e6)),
                 ((300.,80.e6), (0.971180894e-3,0.106448356e6)),
                 ((500.,3.e6), (0.120241800e-2,0.971934985e6)))
        for ((tk,p), (v,u)) in cases:
            dc,uc = cowat(tk-tc_k, p)
            self.assertAlmostEqual(1./dc, v, places = 10)
            self.assertAlmostEqual(uc, u, places = 3)

    def test_supst(self):
        """Tests supst() method on values from IAPWS97 paper."""
        cases = (((300.,0.0035e6), (0.394913866e2,0.241169160e7)),
                 ((700.,0.0035e6), (0.923015898e2,0.301262819e7)),
                 ((700.,30.e6), (0.542946619e-2,0.246861076e7)))
        for ((tk,p), (v,u)) in cases:
            dc,uc = supst(tk-tc_k, p)
            self.assertAlmostEqual(1./dc, v, places = 7)
            self.assertAlmostEqual(uc, u, places = 1)

    def test_super(self):
        """Tests super() method on values from IAPWS97 paper."""
        cases = (((650.,500.), (0.255837018e8,0.181226279e7)),
                 ((650.,200.), (0.222930643e8,0.226365868e7)),
                 ((750.,500.), (0.783095639e8,0.210206932e7)))
        for ((tk,d), (p,u)) in cases:
            pc,uc = super(d, tk-tc_k)
            self.assertAlmostEqual(pc, p, places = 1)
            self.assertAlmostEqual(uc, u, places = 2)

    def test_sat(self):
        """Tests sat() method on values from IAPWS97 paper."""
        cases = ((300.,0.353658941e4),
                 (500.,0.263889776e7),
                 (600.,0.123443146e8))
        for (tk, ps) in cases:
            psc = sat(tk-tc_k)
            self.assertAlmostEqual(psc, ps, places = 1)

    def test_visc(self):
        """Tests visc() method on values from IAPWS97 paper."""
        cases = (((298.15,998.), 889.735100e-6),
                 ((298.15,1200.), 1437.649467e-6),
                 ((373.15,1000.), 307.883622e-6),
                 ((873.15,600.), 77.430195e-6),
                 ((1173.15,100.), 47.640433e-6))
        for (tk,d), v in cases:
            vc = visc(d, tk-tc_k)
            self.assertAlmostEqual(vc, v, places = 1)

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(IAPWS97TestCase)
    unittest.TextTestRunner(verbosity = 1).run(suite)
