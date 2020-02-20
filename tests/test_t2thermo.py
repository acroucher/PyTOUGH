import unittest
from t2thermo import *

class t2thermoTestCase(unittest.TestCase):

    def test_cowat(self):
        """Tests cowat()."""
        cases = (((300., 3.e6), (997.95721560998174, 112247.43313085975)),
                 ((300., 80.e6), (1029.7256888266911, 106310.47344628950)),
                 ((500., 3.e6), (831.84196191567298, 971985.91117384087)))
        for ((tk,p), (d,u)) in cases:
            dc, uc = cowat(tk - tc_k, p)
            self.assertAlmostEqual(dc, d, places = 3)
            self.assertAlmostEqual(uc, u, places = 3)

    def test_supst(self):
        """Tests supst()."""
        cases = (((300., 0.0035e6), (2.5316826343790743e-2, 2412405.0932077002)),
                 ((700., 0.0035e6), (1.0834441421293962e-2, 3012229.4965919587)),
                 ((700., 30.e6), (183.90041953968711, 2474981.3799304822)))
        for ((tk, p), (d, u)) in cases:
            dc, uc = supst(tk-tc_k, p)
            self.assertAlmostEqual(dc, d, places = 6)
            self.assertAlmostEqual(uc, u, places = 1)

    def test_sat(self):
        """Tests sat() method."""
        cases = ((300., 0.35323426e4),
                 (500., 0.263961572e7),
                 (600., 0.123493902e8))
        for (tk, ps) in cases:
            ts = tk - tc_k
            psc = sat(ts)
            self.assertAlmostEqual(psc, ps, places = 1)
            tsc = tsat(psc)
            self.assertAlmostEqual(tsc, ts, places = 3)

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(t2thermoTestCase)
    unittest.TextTestRunner(verbosity = 1).run(suite)
