"""Tests PyTOUGH t2listing module by opening a range of different listing files and checking results against
expected values."""

import unittest
import os
from t2listing import *

class t2listingTestCase(unittest.TestCase):

    def table_test(self, step, tables):
        """Tests entire table data against expected values stored as numpy *.npy files."""
        self.listing.step = step
        base = self.base + '_' + str(self.listing.step) + '_'
        for tablename in tables:
            if tablename in self.listing._table:
                self.assertTrue(np.allclose(self.listing._table[tablename]._data,
                                 np.load(base + tablename + '.npy')),
                                tablename + ' table has incorrect entries')
            else: self.fail(tablename + ' table not found in listing')

    def table_spot_test(self, step, tablename, row, col, value):
        """Tests a listing table entry against an expected value."""
        self.listing.step = step
        self.assertEqual(self.listing._table[tablename][row][col], value)

    def history_test(self, specs, value):
        """Tests a listing history against an expected value."""
        hist = self.listing.history(specs)
        self.assertTrue(all([np.allclose(h, v) for h,v in zip(hist, value)]))

#--------------------------------------------------------------------------------
# AUTOUGH2 tests

    def test_AUTOUGH2_1(self):
        """AUTOUGH2 case 1"""

        self.base = os.path.join('listing', 'AUTOUGH2', '1', 'case1')
        self.listing = t2listing(self.base + '.listing')

        self.assertEqual(self.listing.simulator, 'AUTOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 1)

        self.table_test(27, ['element', 'connection', 'generation'])

        self.table_spot_test(27, 'element', 'AJ190', 'Temperature', 0.16235E+02)
        self.table_spot_test(27, 'connection', ('AF 95','AE 95'),'Enthalpy', 0.658346E+05)
        self.table_spot_test(27, 'generation', ('AI615','RW615'), 'Generation rate', 0.26160E+01)
        self.history_test(('e', 'AK668', 'Temperature'), ((0.1e16), (0.18803E+02)))

        self.listing.close()

    def test_AUTOUGH2_2(self):
        """AUTOUGH2 case 2"""

        self.base = os.path.join('listing', 'AUTOUGH2', '2', 'case2')
        self.listing = t2listing(self.base + '.listing')

        self.assertEqual(self.listing.simulator, 'AUTOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 26)

        self.table_test(26, ['element', 'generation'])

        self.table_spot_test(17, 'element', 'AA 12', 'Temperature', 0.39926E+03)
        self.table_spot_test(26, 'generation', ('AA  1','WEL 1'), 'Enthalpy', 0.22066E+07)

        self.listing.close()

    def test_AUTOUGH2_2_skip(self):
        """Test skip_table parameter using AUTOUGH2 case 2"""

        self.base = os.path.join('listing', 'AUTOUGH2', '2', 'case2')

        self.listing = t2listing(self.base + '.listing', skip_tables = 'generation')
        self.assertEqual(self.listing.simulator, 'AUTOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 26)

        self.table_test(26, ['element'])
        self.table_spot_test(17, 'element', 'AA 12', 'Temperature', 0.39926E+03)

        self.listing.close()

        self.listing = t2listing(self.base + '.listing', skip_tables = 'element')
        self.assertEqual(self.listing.simulator, 'AUTOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 26)

        self.table_test(26, ['generation'])
        self.table_spot_test(26, 'generation', ('AA  1', 'WEL 1'), 'Enthalpy', 0.22066E+07)

        self.listing.close()

    def test_AUTOUGH2_3(self):
        """AUTOUGH2 case 3"""

        self.base = os.path.join('listing', 'AUTOUGH2', '3', 'case3')
        self.listing = t2listing(self.base + '.listing')

        self.assertEqual(self.listing.simulator, 'AUTOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 2)

        self.table_test(40, ['element', 'connection', 'generation'])

        self.table_spot_test(40, 'element', 'EE929', 'Liquid density', 0.85194E+03)
        self.table_spot_test(40, 'connection', ('EE 62', 'DD 62'), 'Liquid mass flow', -0.488999E-01)

        self.listing.close()

    def test_AUTOUGH2_4(self):
        """AUTOUGH2 case 4"""

        self.base = os.path.join('listing', 'AUTOUGH2', '4', 'case4')
        self.listing = t2listing(self.base + '.listing')

        self.assertEqual(self.listing.simulator, 'AUTOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 6)

        self.table_test(68, ['element', 'connection', 'generation'])

        self.table_spot_test(68, 'element', 'AH929', 'Pressure', 0.10178E+06)
        self.table_spot_test(55, 'connection', ('AN672', 'GS672'), 'Heat flow', -0.223477E+07)
        gh = np.load(self.base + '_history.npy')
        self.history_test(('g', ('AT310', 'OO 12'), 'Generation rate'), (gh[:, 0], gh[:, 1]))

        self.listing.close()

    def test_AUTOUGH2_5(self):
        """AUTOUGH2 case 5"""

        self.base = os.path.join('listing', 'AUTOUGH2', '5', 'case5')
        self.listing = t2listing(self.base + '.listing')

        self.assertEqual(self.listing.simulator, 'AUTOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 1)

        self.table_test(7, ['element', 'connection', 'generation'])

        self.table_spot_test(7, 'element', 'gea 7', 'Temperature', 76.508)
        self.table_spot_test(7, 'connection', ('eer 2', 'ATM 0'), 'Liquid velocity', -1.379604E-07)
        self.history_test(('g', ('fav 1', 'D1015'), 'Generation rate'),
                                    ((86400., 259200., 604800., 1296000., 2678400., 5443200., 7862000.),
                                     (-0.00659804, -0.00590907, -0.00566929,
                                      -0.00549556, -0.00534694, -0.00523802, -0.00519849)))

        self.listing.close()

    def test_AUTOUGH2_6(self):
        """AUTOUGH2 test with no simulation title"""

        self.base = os.path.join('listing', 'AUTOUGH2', '6', 'case6')
        self.listing = t2listing(self.base + '.listing')

        self.assertEqual(self.listing.simulator, 'AUTOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 1)
        self.assertEqual(self.listing.num_times, 11)

        self.table_test(5027, ['element', 'connection', 'generation'])


        self.table_spot_test(5027, 'element', 'asr 5', 'CO2 mass fractio', 4.57248E-03)
        self.table_spot_test(5027, 'connection', ('ags 6', 'ahs 6'), 'Heat flow', -3.527687E+01)
        gh = np.load(self.base + '_history.npy')
        self.history_test(('g', ('apn24', 'NPB 9'), 'Generation rate'), (gh[:, 0],  gh[:, 1]))

        self.listing.close()

    def test_AUTOUGH2_7(self):
        """AUTOUGH2 case 7"""

        self.base = os.path.join('listing', 'AUTOUGH2', '7', 'case7')
        self.listing = t2listing(self.base+'.listing')

        self.assertEqual(self.listing.simulator,'AUTOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 1)
        self.assertEqual(self.listing.num_times, 6)

        self.table_test(1551, ['element', 'generation'])

        self.table_spot_test(1551,'element','arp 6','CO2 partial pres', 2.16650E+05)
        self.table_spot_test(1551,'generation',('Auw10','pcp 3'),'Generation rate', -2.11)

        # This test used to fail because it has SHORT items with different order from that in the
        # main tables:
        generKeys = [('Auw18', 'pbp21'), ('Zwc18', 'pbp 8'), ('Asx17', 'pbp 9')]
        historyList = []
        cols = ['Enthalpy']
        for k in generKeys:
            for c in cols: historyList += [('g', k, c)]
        expected = np.load(self.base+'_history.npy')
        results = self.listing.history(historyList)
        for i, (t,v) in enumerate(results):
            self.assertTrue(np.allclose(t, expected[:,0]))
            self.assertTrue(np.allclose(v, expected[:, i + 1]))

        # This test used to fail because the generation table has a missing index:
        gh = np.load(self.base+'_history_pap9.npy')
        self.history_test(('g',('Avx 8','pap 9'),'Generation rate'),(gh[:,0], gh[:,1]))

        self.listing.close()

    def test_AUTOUGH2_8(self):
        """AUTOUGH2 radial first cell history"""

        self.base = os.path.join('listing', 'AUTOUGH2', '8', 'case8')
        self.listing = t2listing(self.base + '.listing')

        self.assertEqual(self.listing.simulator, 'AUTOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 24)
        self.assertEqual(self.listing.num_times, 24)

        t, p = self.listing.history(('e', '  a 1', 'Pressure'))
        pexpected = np.load(self.base + '_history.npy')
        self.assertTrue(np.allclose(p, pexpected))

        self.listing.close()

#--------------------------------------------------------------------------------
# TOUGH2 tests

    def test_TOUGH2_1(self):
        """TOUGH2 EOS1 r1q sample problem"""

        self.base = os.path.join('listing', 'TOUGH2', '1', 'r1q')
        self.listing = t2listing(self.base + '.listing')

        self.assertEqual(self.listing.simulator, 'TOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 1)

        self.table_test(31, ['element', 'generation'])

        self.table_spot_test(31, 'element', 'A1 45', 'T', 0.27473E+03)
        self.table_spot_test(31, 'generation', ('A1  1', 'wel 1'), 'P(WB)', 0.48316E+07)

        self.listing.close()

    def test_TOUGH2_2(self):
        """TOUGH2 EOS1 rfp (5-spot) sample problem"""

        self.base = os.path.join('listing', 'TOUGH2', '2', 'rfp')
        self.listing = t2listing(self.base + '.listing')

        self.assertEqual(self.listing.simulator, 'TOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 3)

        self.table_test(63, ['element', 'generation'])

        self.table_spot_test(63, 'element', ' FF 1', 'P', 0.84657E+07)
        self.history_test(('e', ' IC 1', 'P'),
                                    ((0.15779E+09, 0.78894E+09, 0.11519E+10),
                                     (0.85906E+07, 0.82426E+07, 0.81216E+07)))
        self.history_test(('g',  1,  'GENERATION RATE'),
                                    ((0.15779E+09, 0.78894E+09, 0.11519E+10),
                                     (-3.75, -3.75, -3.75)))

        self.listing.close()

    def test_TOUGH2_3(self):
        """TOUGH2 case 3"""

        self.base = os.path.join('listing', 'TOUGH2', '3', 'OUTFILE')
        self.listing = t2listing(self.base)

        self.assertEqual(self.listing.simulator, 'TOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 1)

        self.table_test(100, ['element', 'connection', 'generation', 'primary'])

        self.table_spot_test(100, 'element', 'CA275', 'DL', 0.99918E+03)
        self.table_spot_test(100, 'connection', ('AC 85', 'AC175'), 'FLOH/FLOF', -0.915165E+05)
        self.table_spot_test(100, 'generation', ('AT347', 'SP347'), 'ENTHALPY', 0.52221E+06)

        self.listing.close()

    def test_TOUGH2_4(self):
        """TOUGH2 case 4"""

        self.base = os.path.join('listing', 'TOUGH2', '4', 'case4')
        self.listing = t2listing(self.base+'.out')

        self.assertEqual(self.listing.simulator, 'TOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 4)

        self.table_test(200, ['element', 'connection', 'generation', 'primary'])

        self.table_spot_test(10, 'element', ' cx 1', 'DW', 998.723)
        self.table_spot_test(11, 'element', ' bi 1', 'P', 0.97856E+06)

        self.listing.close()

    def test_TOUGH2_4_skip(self):
        """Test skip_tables parameter with TOUGH2 case 4"""

        self.base = os.path.join('listing', 'TOUGH2', '4', 'case4')

        self.listing = t2listing(self.base + '.out', skip_tables = 'connection')
        self.assertEqual(self.listing.simulator, 'TOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 4)
        self.table_test(200, ['element', 'generation', 'primary'])
        self.table_spot_test(10, 'element', ' cx 1', 'DW', 998.723)
        self.table_spot_test(11, 'element', ' bi 1', 'P', 0.97856E+06)
        self.listing.close()

        self.listing = t2listing(self.base + '.out',
                                 skip_tables = ['element', 'connection', 'generation'])
        self.assertEqual(self.listing.simulator, 'TOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 4)
        self.table_test(200, ['primary'])

        self.listing.close()

    def test_TOUGH2_5(self):
        """TOUGH2 case 5"""

        self.base = os.path.join('listing', 'TOUGH2', '5', 'case5')
        self.listing = t2listing(self.base + '.out')

        self.assertEqual(self.listing.simulator, 'TOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 4)

        self.table_test(106, ['element', 'generation'])

        self.table_spot_test(41, 'element', ' ao 2', 'T', 45.0)
        self.table_spot_test(106, 'element', '  a 6', 'XCO2aq', 0.13145E-35)
        self.table_spot_test(106, 'generation', (' at25', 'inj 1'), 'GENERATION RATE', 100.)

        self.listing.close()

    def test_TOUGH2_6(self):
        """TOUGH2 case 6"""

        self.base = os.path.join('listing', 'TOUGH2', '6', 'case6')
        self.listing = t2listing(self.base)

        self.assertEqual(self.listing.simulator, 'TOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 1)

        self.table_test(10, ['element', 'connection', 'generation', 'primary'])

        self.table_spot_test(10, 'element', 'BB101', 'XAIR(LIQ)', 0.15822E-04)
        self.table_spot_test(10, 'connection', ('B7101', 'B7201'), 'FLOH', -.212171E-01)
        self.table_spot_test(10, 'primary', 'AM601', 'K(GAS)', 0.11031)
        self.table_spot_test(10, 'generation', ('AQ101', '    0'), 'GENERATION RATE', 0.281E-04)

        self.listing.close()

    def test_TOUGH2_7(self):
        """TOUGH2 case 7 (Petrasim) """

        # This listing comes from PetraSim's EOS3 and contains extra
        # blank lines between the headers and units lines in the
        # connection tables.

        self.base = os.path.join('listing', 'TOUGH2', '7', 'case7')
        self.listing = t2listing(self.base + '.out')

        self.assertEqual(self.listing.simulator, 'TOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 2)

        self.table_test(561, ['element', 'connection', 'generation'])

        self.table_spot_test(561, 'element', '2  50', 'T', 0.10542E+03)
        self.table_spot_test(561, 'connection', ('  224', '  225'), 'FLOF', -0.106635E+01)
        self.table_spot_test(561, 'generation', (' 1253', '   54'), 'ENTHALPY', 0.12762E+07)
        self.history_test(('e', '2  24', 'DL'),
                                    ((1.0, 0.41000E+09),
                                     (0.90533E+03, 0.90519E+03)))

        self.listing.close()

    def test_TOUGH2_8(self):
        """TOUGH2 case 8"""

        self.base = os.path.join('listing', 'TOUGH2', '8', 'OUTFILE')
        self.listing = t2listing(self.base)

        self.assertEqual(self.listing.simulator, 'TOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 5)

        self.table_test(20, ['element', 'connection', 'primary', 'generation'])

        self.table_spot_test(5, 'element', 'A1 30', 'PCO2', 0.472742E-49)
        self.table_spot_test(25, 'connection', ('A1 54', 'A1 55'), 'FLO(GAS)', 0.183886E+01)
        self.table_spot_test(10, 'primary', 'A1 20', 'X2', 0.134142)
        self.history_test(('e', 'A1 10', 'PCO2'),
                          ((0.160000E+07, 0.120000E+08, 0.536000E+08, 0.232800E+09, 0.100000E+10),
                           (0.120001E-35, 0.357205E-36, 0.117162E-38, 0.452751E-44, 0.265833E-52)))

        self.listing.close()

    def test_TOUGH2_9(self):
        """TOUGH2 case 9"""

        # generation tables have only one value per line

        self.base = os.path.join('listing', 'TOUGH2', '9', 'OUTFILE')
        self.listing = t2listing(self.base)

        self.assertEqual(self.listing.simulator, 'TOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 3)

        self.table_test(40, ['element', 'connection', 'primary', 'generation'])

        self.table_spot_test(20, 'element', ' at 1', 'P', 0.57058E+06)
        self.table_spot_test(43, 'connection', (' ce 1', ' co 1'), 'FLOH/FLOF', -0.400631E+07)
        self.table_spot_test(43, 'primary', ' cg 3', 'X1',  0.25510E+07)
        self.table_spot_test(20, 'generation', (' as10', ' as10'), 'GENERATION RATE',  0.1e4)

        self.history_test(('e', ' cu 3', 'P'),
                          ((0.10486E+09, 0.48104E+14, 0.75591E+14),
                           (0.24513E+07, 0.25510E+07, 0.25510E+07)))

        self.listing.close()

    def test_TOUGH2_10(self):
        """TOUGH2 case 10"""

        # from iTOUGH2, with non-exponential formatted pressures

        self.base = os.path.join('listing', 'TOUGH2', '10', 'case10')
        self.listing = t2listing(self.base + '.listing')

        self.assertEqual(self.listing.simulator, 'TOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 1)

        self.table_test(1, ['element', 'connection', 'generation'])

        self.listing.close()

    def test_TOUGH2_11(self):
        """TOUGH2 case 11"""

        # from iTOUGH2, with extra primary table on last time step:
        # this table will not be parsed because it is not present in
        # the first set of results

        self.base = os.path.join('listing', 'TOUGH2', '11', 'case11')
        self.listing = t2listing(self.base + '.listing')

        self.assertEqual(self.listing.simulator, 'TOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 3)

        self.table_test(16974, ['element', 'connection', 'generation'])

        self.assertFalse(hasattr(self.listing, 'primary'))

        self.listing.close()

    def test_TOUGH2_12(self):
        """TOUGH2 case 12"""

        # iTOUGH2 ECO2M, with 'I' column in element table

        self.base = os.path.join('listing', 'TOUGH2', '12', 'case12')
        self.listing = t2listing(self.base + '.out')

        self.assertEqual(self.listing.simulator, 'TOUGH2')
        self.assertEqual(self.listing.num_fulltimes, 2)

        self.table_test(515, ['element'])
        self.history_test(('e','A1080','P'),
                          ((2.592000E+06, 8.640000E+08),
                           (0.197965E+08, 0.215416E+08)))

        self.listing.close()

#--------------------------------------------------------------------------------
# TOUGH2-MP tests

    def test_TOUGH2MP_1(self):
        """TOUGH2-MP case 1"""

        self.base = os.path.join('listing', 'TOUGH2-MP', '1', 'OUTPUT_DATA')
        self.listing = t2listing(self.base)

        self.assertEqual(self.listing.simulator, 'TOUGH2_MP')
        self.assertEqual(self.listing.num_fulltimes, 1)

        self.table_test(144, ['element', 'connection'])

        self.table_spot_test(144, 'element', 'dec12', 'XAIRL', 0.11783E-26)

        self.listing.close()

    def test_TOUGH2MP_2(self):
        """TOUGH2-MP case 2"""

        self.base = os.path.join('listing', 'TOUGH2-MP', '2', 'OUTPUT_DATA')
        self.listing = t2listing(self.base)

        self.assertEqual(self.listing.simulator, 'TOUGH2_MP')
        self.assertEqual(self.listing.num_fulltimes, 1)

        self.table_test(55, ['element', 'connection'])

        self.table_spot_test(55, 'element', 'AZ332', 'T', 0.22122E+03)
        self.table_spot_test(55, 'connection', ('AZ145', 'AZ325'),
                             'FLO(LIQ.)', -0.186357E-01)
        self.history_test(('e', 'AZ153', 'XAIRL'), ((0.79000E+15), (0.59224E-05)))

        self.listing.close()

    def test_TOUGH2MP_3(self):
        """TOUGH2-MP case 3"""

        self.base = os.path.join('listing', 'TOUGH2-MP', '3', 'OUTPUT_DATA')
        self.listing = t2listing(self.base)

        self.assertEqual(self.listing.simulator, 'TOUGH2_MP')
        self.assertEqual(self.listing.num_fulltimes, 1)

        self.table_test(40, ['element', 'connection'])

        self.table_spot_test(40, 'element', 'dca 1', 'P', 0.596969E+06)
        self.table_spot_test(40, 'connection', ('eoi 1', 'eok 1'),
                             'ENTHALPY', 0.564181E+06)

        self.listing.close()

    def test_TOUGH2MP_4(self):
        """TOUGH2-MP case 4"""

        self.base = os.path.join('listing', 'TOUGH2-MP', '4', 'OUTPUT_DATA')
        self.listing = t2listing(self.base)

        self.assertEqual(self.listing.simulator, 'TOUGH2_MP')
        self.assertEqual(self.listing.num_fulltimes, 1)

        self.table_test(257, ['element', 'connection'])

        self.table_spot_test(257, 'element', 'dgs 1', 'T', 44.3)
        self.table_spot_test(257, 'element', 'dgq 1', 'XNACL', 0.15579E-02)
        self.table_spot_test(257, 'connection', ('dkg13', 'dki13'), 'FLO(NaCl)', -.138999E-14)

        self.listing.close()

    def test_TOUGH2MP_5(self):
        """TOUGH2-MP case 5"""

        self.base = os.path.join('listing', 'TOUGH2-MP', '5', 'OUTPUT_DATA')
        self.listing = t2listing(self.base)

        self.assertEqual(self.listing.simulator, 'TOUGH2_MP')
        self.assertEqual(self.listing.num_fulltimes, 2)

        self.table_test(368, ['element', 'connection'])

        self.table_spot_test(534, 'element', 'A1609', 'DG',  776.35)
        self.table_spot_test(534, 'connection', ('A1604', 'A1704'), 'FLOH/FLOF',  0.584530E+06)

        self.history_test(('e', 'A1406', 'SG'), ((0.315576E+08, 0.631152E+08),
                                                 (0.44087E+00, 0.49543E+00)))

        self.listing.close()

    def test_TOUGH2MP_6(self):
        """TOUGH2-MP case 6"""

        self.base = os.path.join('listing', 'TOUGH2-MP', '6', 'OUTPUT_DATA')
        self.listing = t2listing(self.base)

        self.assertEqual(self.listing.simulator, 'TOUGH2_MP')
        self.assertEqual(self.listing.num_fulltimes, 6)

        self.table_test(25, ['element', 'connection'])

        self.table_spot_test(10, 'element', 'A1 16', 'T',  0.27549E+03)
        self.table_spot_test(25, 'connection', ('A1  1', 'A1  2'), 'VEL(LIQ.)',  0.33260E-05)

        self.history_test(('e', 'A1100', 'XHYDG'),
                          ((0.13188E+07, 0.77734E+07, 0.32441E+08,
                            0.12218E+09, 0.47250E+09, 0.10000E+10),
                           (0.12662E-04, 0.12662E-04, 0.12662E-04,
                            0.12661E-04, 0.12379E-04, 0.11075E-04)))

        self.listing.close()

    def test_TOUGH2MP_7(self):
        """TOUGH2-MP case 7"""

        self.base = os.path.join('listing', 'TOUGH2-MP', '7', 'OUTPUT_DATA')
        self.listing = t2listing(self.base)

        self.assertEqual(self.listing.simulator, 'TOUGH2_MP')
        self.assertEqual(self.listing.num_fulltimes, 2)

        self.table_test(10, ['element', 'connection'])

        self.table_spot_test(1, 'element', 'A1303', 'PRES',  0.81483E+05)
        self.table_spot_test(1, 'connection', ('A1103', 'A1203'), 'VEL(fract)', -0.53098E-03)

        self.history_test(('e', 'A1301', 'S(liq)'),
                          ((0.10000E-08, 0.50258E+05), (0.14990, 0.22059)))

        self.listing.close()

    def test_TOUGH2MP_7_skip(self):
        """Test skip_tables parameter with TOUGH2-MP case 7"""

        self.base = os.path.join('listing', 'TOUGH2-MP', '7', 'OUTPUT_DATA')

        self.listing = t2listing(self.base, skip_tables = 'connection')
        self.assertEqual(self.listing.simulator, 'TOUGH2_MP')
        self.assertEqual(self.listing.num_fulltimes, 2)
        self.table_test(10, ['element'])
        self.assertFalse(hasattr(self.listing, 'connection'))
        self.table_spot_test(1, 'element', 'A1303', 'PRES',  0.81483E+05)
        self.history_test(('e', 'A1301', 'S(liq)'),
                          ((0.10000E-08, 0.50258E+05), (0.14990, 0.22059)))
        self.listing.close()

        self.listing = t2listing(self.base, skip_tables = 'element')
        self.assertEqual(self.listing.simulator, 'TOUGH2_MP')
        self.assertEqual(self.listing.num_fulltimes, 2)
        self.table_test(10, ['connection'])
        self.assertFalse(hasattr(self.listing, 'element'))
        self.table_spot_test(1, 'connection', ('A1103', 'A1203'), 'VEL(fract)', -0.53098E-03)
        self.listing.close()

#--------------------------------------------------------------------------------
# TOUGH+ tests

    def test_TOUGHplus_1(self):
        """TOUGH+ case 1"""

        self.base = os.path.join('listing', 'TOUGHplus', '1', 'case1')
        self.listing = t2listing(self.base + '.dat')

        self.assertEqual(self.listing.simulator, 'TOUGH+')
        self.assertEqual(self.listing.num_fulltimes, 3)

        self.table_test(5000, ['element', 'element1', 'element2', 'primary', 'connection'])

        self.table_spot_test(3000, 'element', '  a22', 'Temperature', 18.1)
        self.table_spot_test(5000, 'connection', ('  a10', '  a 9'), 'Heat Flow', 6.1611)

        self.listing.close()

    def test_TOUGHplus_1_skip(self):
        """Test skip_tables parameter with TOUGH+ case 1"""

        self.base = os.path.join('listing', 'TOUGHplus', '1', 'case1')

        self.listing = t2listing(self.base + '.dat', skip_tables = 'connection')
        self.assertEqual(self.listing.simulator, 'TOUGH+')
        self.assertEqual(self.listing.num_fulltimes, 3)
        self.table_test(5000, ['element', 'element1', 'element2', 'primary'])
        self.table_spot_test(3000, 'element', '  a22', 'Temperature', 18.1)
        self.assertFalse(hasattr(self.listing, 'connection'))
        self.listing.close()

        self.listing = t2listing(self.base + '.dat', skip_tables = 'primary')
        self.assertEqual(self.listing.simulator, 'TOUGH+')
        self.assertEqual(self.listing.num_fulltimes, 3)
        self.table_test(5000, ['element', 'element1', 'element2', 'connection'])
        self.table_spot_test(3000,'element','  a22','Temperature',18.1)
        self.table_spot_test(5000,'connection',('  a10','  a 9'),'Heat Flow',6.1611)
        self.assertFalse(hasattr(self.listing, 'primary'))
        self.listing.close()

        skip = ['element', 'element1', 'element2']
        self.listing = t2listing(self.base+'.dat', skip_tables = skip)
        self.assertEqual(self.listing.simulator, 'TOUGH+')
        self.assertEqual(self.listing.num_fulltimes, 3)
        self.table_test(5000, ['primary', 'connection'])
        self.table_spot_test(5000,'connection',('  a10','  a 9'),'Heat Flow',6.1611)
        for name in skip:
            self.assertFalse(hasattr(self.listing, name))
        self.listing.close()

    def test_TOUGHplus_2(self):
        """TOUGH+ case 2"""

        self.base = os.path.join('listing', 'TOUGHplus', '2', 'case2')
        self.listing = t2listing(self.base + '.dat')

        self.assertEqual(self.listing.simulator, 'TOUGH+')
        self.assertEqual(self.listing.num_fulltimes, 2)

        self.table_test(800, ['element', 'element1', 'element2', 'primary', 'connection'])

        self.table_spot_test(979, 'element', ' bf 2', 'S_aqueous', 1.2152E-01)
        self.table_spot_test(979, 'element1', ' cm 2', 'Dens_Gas', 1.9497E+02)

        self.listing.close()

    def test_TOUGHplus_3(self):
        """TOUGH+ case 3"""

        self.base = os.path.join('listing', 'TOUGHplus', '3', '1p_out')
        self.listing = t2listing(self.base + '.dat')

        self.assertEqual(self.listing.simulator, 'TOUGH+')
        self.assertEqual(self.listing.num_fulltimes, 17)

        self.table_test(26, ['element', 'element1', 'element2', 'primary', 'connection'])

        self.table_spot_test(18, 'element', 'ina 0', 'P_EqHydr', 5.3720246E+06)

        self.listing.close()

    def test_TOUGHplus_4(self):
        """TOUGH+ case 4"""

        self.base = os.path.join('listing', 'TOUGHplus', '4', 't3T_out')
        self.listing = t2listing(self.base + '.dat')

        self.assertEqual(self.listing.simulator, 'TOUGH+')
        self.assertEqual(self.listing.num_fulltimes, 2)

        self.table_test(20, ['element', 'element1', 'element2', 'primary', 'connection'])

        self.table_spot_test(39, 'connection', ('A0031', 'A0032'), 'Aqu Veloc', 2.1206E-07)
        self.table_spot_test(39, 'primary', 'A0020', 'X2', 4.852208E-01)
        self.history_test(('c', ('A0009', 'A0010'), 'Heat Flow'),
                          ((2.088930E+03, 4.065000E+05),
                           (9.5523E+04, 5.7209E+04)))
        self.listing.close()

#--------------------------------------------------------------------------------
# TOUGHREACT tests

    def test_TOUGHREACT_1(self):
        """TOUGHREACT case 1"""

        self.base = os.path.join('listing', 'TOUGHREACT', '1', 'case1')
        self.listing = t2listing(self.base + '.out')

        self.assertEqual(self.listing.simulator, 'TOUGHREACT')
        self.assertEqual(self.listing.num_fulltimes, 2)

        self.table_test(186, ['element', 'connection', 'generation'])

        self.table_spot_test(186, 'element', 'A  13', 'PCAP', -.524952E+04)
        self.table_spot_test(5000, 'connection', ('B  28', 'B  29'),
                             'FLO(CO2)', 0.181833E-08)
        self.history_test(('e', 'B  49', 'T'), (( 0.864000E+07, 0.315576E+09),
                                                (0.219377E+03, 0.219377E+03)))

        self.listing.close()

    def test_TOUGHREACT_2(self):
        """TOUGHREACT case 2"""

        self.base = os.path.join('listing', 'TOUGHREACT', '2', 'case2')
        self.listing = t2listing(self.base + '.out')

        self.assertEqual(self.listing.simulator, 'TOUGHREACT')
        self.assertEqual(self.listing.num_fulltimes, 3)

        self.table_test(25, ['element'])

        self.table_spot_test(27, 'element', ' ag 3', 'P', 0.13243E+07)
        self.history_test(('e', '  c 5', 'P'), ((0.100000E+01, 0.315360E+08, 0.630720E+08),
                                                (0.10130E+06, 0.23016E+07, 0.23030E+07)))

        self.listing.close()

#--------------------------------------------------------------------------------
# TOUGH3 tests

    def test_TOUGH3_1(self):
        """TOUGH3 case 1 (r1q)"""

        self.base = os.path.join('listing', 'TOUGH3', '1', 'OUTPUT')
        self.listing = t2listing(self.base)

        self.assertEqual(self.listing.simulator, 'TOUGH3')
        self.assertEqual(self.listing.num_fulltimes, 1)

        self.table_test(31, ['element'])

        self.table_spot_test(1, 'element', 'A1 55', 'SAT_G', 0.2133)

        self.listing.close()

    def test_TOUGH3_2(self):
        """TOUGH3 case 2 (rfp)"""

        self.base = os.path.join('listing', 'TOUGH3', '2', 'OUTPUT')
        self.listing = t2listing(self.base)

        self.assertEqual(self.listing.simulator, 'TOUGH3')
        self.assertEqual(self.listing.num_fulltimes, 2)

        self.table_test(20, ['element'])

        self.table_spot_test(55, 'element', '2CB 1', 'X_CO2_G', 0.6434)

        self.listing.close()

    def test_TOUGH3_3(self):
        """TOUGH3 case 3 (rvf)"""

        self.base = os.path.join('listing', 'TOUGH3', '3', 'OUTPUT')
        self.listing = t2listing(self.base)

        self.assertEqual(self.listing.simulator, 'TOUGH3')
        self.assertEqual(self.listing.num_fulltimes, 1)

        self.table_test(46, ['element', 'connection'])

        self.table_spot_test(46, 'element', 'A1311', 'TEMP', 190.5)
        self.table_spot_test(46, 'connection', ('A1803', 'A1804'), 'FLOW_L', -0.6251)

        self.listing.close()

    def test_TOUGH3_4(self):
        """TOUGH3 case 4"""

        self.base = os.path.join('listing', 'TOUGH3', '4', 'OUTPUT')
        self.listing = t2listing(self.base)

        self.assertEqual(self.listing.simulator, 'TOUGH3')
        self.assertEqual(self.listing.num_fulltimes, 1)

        self.table_test(345, ['element', 'connection', 'generation'])

        self.table_spot_test(345, 'element', 'djg38', 'TEMP', 187.1)
        self.table_spot_test(345, 'connection', ('dgi 6', 'dgi 0'), 'FLOW', 0.3067E+02)
        self.table_spot_test(345, 'generation', ('apo38', 'apo98'), 'GEN', 2.371)

        self.listing.close()

#--------------------------------------------------------------------------------

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(t2listingTestCase)
    unittest.TextTestRunner(verbosity = 1).run(suite)
