"""Tests PyTOUGH t2incon module using a range of different incon/save files."""

import unittest
import os
from t2incons import *

class t2inconTestCase(unittest.TestCase):

    def check_variable(self):
        v = np.load(os.path.join(self.path, 'variable.npy'))
        self.assertTrue(np.allclose(self.inc.variable, v),
                        'variable table has incorrect entries')

    def check_porosity(self):
        p = np.load(os.path.join(self.path, 'porosity.npy'))
        self.assertTrue(np.allclose(self.inc.porosity, p),
                        'porosity table has incorrect entries')

    def check_no_permeabilities(self):
        self.assertTrue(all([k is None for k in self.inc.permeability]))

    def check_permeabilities(self):
        k = np.load(os.path.join(self.path, 'permeability.npy'))
        self.assertTrue(np.allclose(self.inc.permeability, k),
                        'permeability table has incorrect entries')

    def check_timing(self, expected):
        msg = 'Timing information incorrect'
        self.assertEqual(self.inc.timing['kcyc'], expected['kcyc'], msg)
        self.assertEqual(self.inc.timing['iter'], expected['iter'], msg)
        self.assertEqual(self.inc.timing['nm'], expected['nm'], msg)
        self.assertAlmostEqual(self.inc.timing['tstart'], expected['tstart'], msg=msg)
        self.assertAlmostEqual(self.inc.timing['sumtim'], expected['sumtim'], msg=msg)

#--------------------------------------------------------------------------------

    def test_AUTOUGH2_1(self):
        """AUTOUGH2 case1"""
        self.path = os.path.join('incon', 'AUTOUGH2', '1')
        self.inc = t2incon(os.path.join(self.path,'case1.incon'))
        self.assertTrue(self.inc.simulator == 'TOUGH2', 'incorrect simulator')
        self.check_variable()
        self.check_porosity()
        self.check_no_permeabilities()
        self.assertTrue(self.inc.timing is None)

        filename = os.path.join(self.path, 'test.incon')
        self.inc.write(filename)
        self.inc = t2incon(filename)
        self.assertTrue(self.inc.simulator == 'TOUGH2', 'incorrect simulator')
        self.check_variable()
        self.check_porosity()
        self.check_no_permeabilities()
        self.assertTrue(self.inc.timing is None)
        os.remove(filename)

    def test_AUTOUGH2_2(self):
        """AUTOUGH2 case2"""
        self.path = os.path.join('incon', 'AUTOUGH2', '2')
        self.inc = t2incon(os.path.join(self.path, 'case2.incon'))
        self.assertTrue(self.inc.simulator == 'TOUGH2', 'incorrect simulator')
        self.check_variable()
        self.check_porosity()
        self.check_no_permeabilities()
        self.assertTrue(self.inc.timing is None)

        filename = os.path.join(self.path, 'test.incon')
        self.inc.write(filename)
        self.inc = t2incon(filename)
        self.assertTrue(self.inc.simulator == 'TOUGH2', 'incorrect simulator')
        self.check_variable()
        self.check_porosity()
        self.check_no_permeabilities()
        self.assertTrue(self.inc.timing is None)
        os.remove(filename)

    def test_AUTOUGH2_3(self):
        """AUTOUGH2 case 3"""
        self.path = os.path.join('incon', 'AUTOUGH2', '3')
        self.inc = t2incon(os.path.join(self.path, 'case3.incon'))
        self.assertTrue(self.inc.simulator == 'TOUGH2', 'incorrect simulator')
        self.check_variable()
        self.check_porosity()
        self.check_no_permeabilities()
        timing = {'kcyc': 30, 'iter': 145, 'nm': 34, 'tstart': 0.0, 'sumtim': 0.106496000E+17}
        self.check_timing(timing)

        filename = os.path.join(self.path, 'test.incon')
        self.inc.write(filename, reset = False)
        self.inc = t2incon(filename)
        self.assertTrue(self.inc.simulator == 'TOUGH2', 'incorrect simulator')
        self.check_variable()
        self.check_porosity()
        self.check_no_permeabilities()
        self.check_timing(timing)
        os.remove(filename)

    def test_TOUGH2_1(self):
        """Synthetic TOUGH2 with strange exponents"""
        self.path = os.path.join('incon', 'TOUGH2', '1')
        self.inc = t2incon(os.path.join(self.path, 'case1.incon'))
        self.assertTrue(self.inc.simulator == 'TOUGH2', 'incorrect simulator')
        self.check_variable()
        self.check_no_permeabilities()
        self.assertTrue(self.inc.timing is None)

        filename = os.path.join(self.path, 'test.incon')
        self.inc.write(filename, reset = False)
        self.inc = t2incon(filename)
        self.assertTrue(self.inc.simulator == 'TOUGH2', 'incorrect simulator')
        self.check_variable()
        self.check_no_permeabilities()
        self.assertTrue(self.inc.timing is None)
        os.remove(filename)

    def test_TOUGH2_2(self):
        """TOUGH2 EOS7c"""
        num_vars = 6
        self.path = os.path.join('incon', 'TOUGH2', '2')
        self.inc = t2incon(os.path.join(self.path, 'INCON'), num_variables = num_vars)
        self.assertTrue(self.inc.simulator == 'TOUGH2', 'incorrect simulator')
        self.check_variable()
        self.check_porosity()
        self.check_no_permeabilities()
        self.assertTrue(self.inc.timing is None)

        filename = os.path.join(self.path, 'test.incon')
        self.inc.write(filename, reset = False)
        self.inc = t2incon(filename, num_variables = num_vars)
        self.assertTrue(self.inc.simulator == 'TOUGH2', 'incorrect simulator')
        self.check_variable()
        self.check_porosity()
        self.check_no_permeabilities()
        self.assertTrue(self.inc.timing is None)
        os.remove(filename)

    def test_TOUGH2_3(self):
        """NSEQ, NADD"""
        self.path = os.path.join('incon', 'TOUGH2', '3')
        self.inc = t2incon(os.path.join(self.path, 'test.incon'))
        self.assertTrue(self.inc.simulator == 'TOUGH2', 'incorrect simulator')
        self.assertTrue(self.inc.timing is None)
        self.assertEqual(self.inc['AAA 1'].nseq, None)
        self.assertEqual(self.inc['AAA 1'].nadd, None)
        self.assertEqual(self.inc['AAA 3'].nseq, 3)
        self.assertEqual(self.inc['AAA 3'].nadd, 2)

        filename = os.path.join(self.path + 'out.incon')
        self.inc.write(filename, reset = False)
        self.inc = t2incon(filename)
        self.assertTrue(self.inc.simulator == 'TOUGH2', 'incorrect simulator')
        self.assertTrue(self.inc.timing is None)
        self.assertEqual(self.inc['AAA 1'].nseq, None)
        self.assertEqual(self.inc['AAA 1'].nadd, None)
        self.assertEqual(self.inc['AAA 3'].nseq, 3)
        self.assertEqual(self.inc['AAA 3'].nadd, 2)
        os.remove(filename)

    def test_TOUGHREACT_1(self):
        """TOUGHREACT save file 1"""
        self.path = os.path.join('incon', 'TOUGHREACT', '1')
        self.inc = t2incon(os.path.join(self.path, 'SAVE_1'))
        self.assertTrue(self.inc.simulator == 'TOUGHREACT', 'incorrect simulator')
        self.check_variable()
        self.check_porosity()
        self.check_permeabilities()
        timing = {'kcyc': 11100, 'iter': 40102, 'nm': 1, 'tstart': 0.0, 'sumtim': 0.52710494E+05}
        self.check_timing(timing)

        filename = os.path.join(self.path, 'test.incon')
        self.inc.write(filename, reset = False)
        self.inc = t2incon(filename)
        self.assertTrue(self.inc.simulator == 'TOUGHREACT', 'incorrect simulator')
        self.check_variable()
        self.check_porosity()
        self.check_permeabilities()
        self.check_timing(timing)
        os.remove(filename)

    def test_transfer_from(self):
        """transfer_from() test"""
        import mulgrids
        dx, dy, dz = [100.]*10, [120.]*8, [10.]*5
        geo = mulgrids.mulgrid().rectangular(dx, dy, dz)
        inc = t2incon()
        for blkname in geo.block_name_list:
            layername = geo.layer_name(blkname)
            colname = geo.column_name(blkname)
            pos = geo.block_centre(layername, colname)
            P = 1.e5 + 0.1 * pos[0] * pos[1] - pos[2]
            T = 240. - 1.e-4 * pos[0] * (1000. - pos[1]) + 0.1 * pos[2]
            inc[blkname] = t2blockincon(block = blkname, variable = [P, T])
        newgeo = mulgrids.mulgrid().rectangular(dx, dy, dz)
        newgeo.refine()
        newinc = t2incon()
        newinc.transfer_from(inc, geo, newgeo)
        test_blks = ['  a 1', ' bs 2', ' cb 3', '  j 4']
        for blkname in test_blks:
            val = inc[blkname].variable
            layername = geo.layer_name(blkname)
            colname = geo.column_name(blkname)
            pos = geo.block_centre(layername, colname)
            newblk = newgeo.block_name_containing_point(pos)
            newval = newinc[newblk].variable
            self.assertEqual(newval, val)

#--------------------------------------------------------------------------------
        
if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(t2inconTestCase)
    unittest.TextTestRunner(verbosity = 1).run(suite)
