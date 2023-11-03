import unittest
import os
from t2data import *

def column_nan_to_num(x):
    """Returns numpy dtype array with NaNs in floating point columns
    converted to zero."""
    if x.dtype.name.startswith('float'):
        return np.nan_to_num(x)
    else:
        for key,dt in x.dtype.fields.items():
            if dt[0].name.startswith('float'): x[key] = np.nan_to_num(x[key])
        return x

default_separator_pressure = 0.55e6
default_2_stage_separator_pressure = [1.45e6, 0.55e6]

class t2data_stats(t2data):
    """Variant of t2data class with extra properties for vital statistics
    of a t2data object- for comparisons between t2data objects that
    should be the same.
    """
    def get_block_names(self):
        return [blk.name for blk in self.grid.blocklist]
    block_names = property(get_block_names)
    def get_block_volumes(self):
        return np.array([blk.volume for blk in self.grid.blocklist])
    block_volumes = property(get_block_volumes)
    def get_rock_names(self):
        return [rt.name for rt in self.grid.rocktypelist]
    rock_names = property(get_rock_names)
    def get_rock_assignments(self):
        return [blk.rocktype.name for blk in self.grid.blocklist]
    rock_assignments = property(get_rock_assignments)
    def get_connection_blocks(self):
        return [[blk.name for blk in con.block] for con in self.grid.connectionlist]
    connection_blocks = property(get_connection_blocks)
    def get_connection_dists(self):
        return np.array([con.distance for con in self.grid.connectionlist])
    connection_dists = property(get_connection_dists)
    def get_connection_areas(self):
        return np.array([con.area for con in self.grid.connectionlist])
    connection_areas = property(get_connection_areas)
    def get_connection_dirns(self):
        return np.array([con.dircos if con.dircos is not None else 0.
                         for con in self.grid.connectionlist])
    connection_dirns = property(get_connection_dirns)
    def get_generator_names(self):
        return [(gen.name,gen.block) for gen in self.generatorlist]
    generator_names = property(get_generator_names)
    def get_generator_rates(self):
        return np.array([gen.gx for gen in self.generatorlist])
    generator_rates = property(get_generator_rates)

    def compare(self, other, testcase, tol = 1.e-7):
        """Compares self with another t2data_stats object and reports results
        to a unittest testcase.  
        """
        testcase.assertEqual(self.block_names, other.block_names)
        testcase.assertTrue(np.allclose(self.block_volumes, other.block_volumes, rtol = tol))
        testcase.assertEqual(self.rock_names, other.rock_names)
        testcase.assertEqual(self.rock_assignments, other.rock_assignments)
        testcase.assertEqual(self.connection_blocks, other.connection_blocks)
        testcase.assertTrue(np.allclose(self.connection_dists, other.connection_dists, rtol = tol))
        testcase.assertTrue(np.allclose(self.connection_areas, other.connection_areas, rtol = tol))
        testcase.assertTrue(np.allclose(self.connection_dirns, other.connection_dirns, rtol = tol))
        testcase.assertEqual(self.generator_names, other.generator_names)
        testcase.assertTrue(np.allclose(self.generator_rates, other.generator_rates, rtol = tol))
        
class t2dataTestCase(unittest.TestCase):

#------------------------------------------------------------------------

    def test_fromscratch(self):
        """construct data file from scratch"""
        filename = 'test_fromscratch.dat'
        geo = mulgrid(os.path.join('mulgrid', 'g6.dat'))
        dat = t2data_stats()
        dat.title = 'Test constructing data file from scratch'
        dat.simulator = 'AUTOUGH2.2EWAV'

        dat.grid = t2grid().fromgeo(geo)
        self.assertEqual(geo.block_name_list, [blk.name for blk in dat.grid.blocklist])

        r=rocktype(name='ARGIL',permeability=[5.e-15,5.e-15,1.e-15])
        r.porosity = .100000001
        r.density = 2500
        r.conductivity = 2.5
        r.specific_heat=1000
        dat.grid.add_rocktype(r)
        for blk in dat.grid.blocklist[0:]: blk.rocktype = r
        dat.grid.clean_rocktypes()

        dat.parameter.update({
                'max_iterations' : 0,
                'print_level' : 3,
                'max_timesteps' : 999,
                'max_duration' : 0,
                'print_interval' : 999,
                'tstart' : 0.0,
                'tstop' : 315360000,
                'const_timestep' : 864000.,
                'gravity' : 9.81})

        options = {1: 1, 10: 2, 11: 2, 12: 2, 16: 5, 20: 1, 23: 1}
        for i,v in options.items(): dat.parameter['option'][i] = v

        dat.parameter['default_incons'] = [1.01325e5, 25, 9.825e4]

        dat.relative_permeability['type'] = 1.0
        dat.relative_permeability['parameters'] = [0.7, 0.0, 1.0, 0.3]
        dat.capillarity['type'] = 1
        dat.capillarity['parameters'] = [0.0, 0.0, 1.0]

        dat.lineq.update({'type' : 2, 'epsilon' : 1e-11, 'max_iterations' : 999,
                          'num_orthog' : 100, 'gauss' : 1})

        dat.multi.update({'num_components' : 2, 'num_equations' : 3,
                          'num_phases' : 2, 'num_secondary_parameters' : 6, 'eos' : 'EWAV'})

        flux = 1.e-3
        layer = geo.layerlist[-1]
        for col in geo.columnlist:
            blkname = geo.block_name(layer.name, col.name)
            gen = t2generator(blkname,blkname,type = 'HEAT', gx = col.area * flux)
            dat.add_generator(gen)

        dat.write(filename)
        dat2 = t2data_stats(filename)
        dat.compare(dat2, self)
        from os import remove
        remove(filename)

#------------------------------------------------------------------------

    def rocktypes_test(self):
        """Tests rocktype table against expected values stored as numpy *.npy file."""
        expected = np.load(self.base + '_rocks.npy')
        rocks = np.array([tuple([rt.name, rt.density, rt.porosity] + list(rt.permeability) +
                             [rt.conductivity, rt.specific_heat])
                          for rt in self.dat.grid.rocktypelist],
                         dtype = [('f0','S5'), ('f1','f'), ('f2','f'), ('f3','f'), ('f4','f'),
                                  ('f5','f'), ('f6','f'), ('f7','f')])
        self.assertTrue((rocks == expected).all(), 'Error in rocktype table')

    def blocks_test(self):
        """Tests block table against expected values stored as numpy *.npy file."""

        expected = np.load(self.base + '_blocks.npy')
        expected['f0'] = np.array([fix_blockname(blk.decode()) for blk in expected['f0']])
        blocks = np.array([(blk.name, blk.rocktype.name, blk.volume)
                           for blk in self.dat.grid.blocklist],
                          dtype = [('f0','S5'), ('f1','S5'), ('f2','f')])
        self.assertTrue((blocks == expected).all(), 'Error in element table')

    def connections_test(self):
        """Tests connection table against expected values stored as numpy *.npy file."""
        expected = np.load(self.base + '_connections.npy')
        for f in ['f0','f1']:
            expected[f] = np.array([fix_blockname(blk.decode()) for blk in expected[f]])
        con = column_nan_to_num(np.array([(con.block[0].name, con.block[1].name, con.direction,
                         con.distance[0], con.distance[1], con.area, con.dircos)
                        for con in self.dat.grid.connectionlist],
                       dtype = [('f0','S5'), ('f1','S5'), ('f2','d'), ('f3','f'),
                                ('f4','f'), ('f5','f'), ('f6','f')]))
        self.assertTrue((con == expected).all(), 'Error in connection table')

    def generators_test(self):
        """Tests generation table against expected values stored as numpy *.npy file."""
        expected = np.load(self.base + '_generators.npy')
        for f in ['f0','f1']:
            expected[f] = np.array([fix_blockname(blk.decode()) for blk in expected[f]])
        gen = column_nan_to_num(np.array([(gen.block, gen.name, gen.type, gen.gx, gen.ex)
                        for gen in self.dat.generatorlist],
                       dtype = [('f0', 'S5'), ('f1', 'S5'), ('f2', 'S4'),
                                ('f3', 'f'), ('f4', 'f')]))
        self.assertTrue((gen == expected).all(), 'Error in generation table')

    def tables_test(self):
        """Tests all tables- rocktypes, blocks, connections and generators."""
        self.rocktypes_test()
        self.blocks_test()
        self.connections_test()
        self.generators_test()

    def write_read_test(self, tol = 1.e-7):
        """Tests if writing the data file out to disk and reading it back in gives
        the same result."""
        filename = self.base + '.compare'

        if self.dat.meshfilename == '':
            compare_meshfilename = ''
            self.dat.write(filename)
            read_dat = t2data_stats(filename)
        else: # separate mesh file(s)
            if isinstance(self.dat.meshfilename, str):
                compare_meshfilename = self.dat.meshfilename + '.compare'
            else:
                compare_meshfilename = tuple([mname + '.compare'
                                              for mname in self.dat.meshfilename])
            self.dat.write(filename, meshfilename = compare_meshfilename)
            read_dat = t2data_stats(filename, meshfilename = compare_meshfilename)

        self.dat.compare(read_dat, self, tol)

        from os import remove
        remove(filename)
        if compare_meshfilename: # clean up any mesh files created
            if isinstance(compare_meshfilename, str): remove(compare_meshfilename)
            else:
                for mname in compare_meshfilename: remove(mname)

    def dict_test(self, testdict, expected):
        """Tests dictionary against supplied expected values."""
        for key in expected:
            if isinstance(testdict[key], np.ndarray):
                self.assertTrue((testdict[key] == expected[key]).all())
            else: self.assertEqual(testdict[key], expected[key])

#------------------------------------------------------------------------

    def test_AUTOUGH2_1(self):
        """AUTOUGH2 case 1"""

        self.base = os.path.join('data', 'AUTOUGH2', '1', 'case1')
        # This data file needs the Fortran read functions because of
        # the formatting in the generation table:
        self.dat = t2data_stats(self.base + '.dat', read_function = fortran_read_function)

        self.assertEqual(self.dat.simulator, 'AUTOUGH2.2EWAV')
        self.dict_test(self.dat.parameter,
                       {'print_level': 3, 'max_timesteps': 5, 'print_interval': 5, 'max_iterations': None,
                        'option': np.array([0,2,0,0,0,0,0,0,0,0,2,2,2,0,0,0,5,0,0,0,1,0,0,1,0]),
                        'tstart': 0.0, 'timestep': [0.1e16], 'print_block': 'SC 26', 'gravity': 9.8065,
                        'default_incons': [.10135e06, .999, .99013e05]})
        self.dict_test(self.dat.multi,
                       {'eos': '', 'num_phases': 2, 'num_equations': 3, 'num_components': 2,
                        'num_secondary_parameters': 6})
        self.dict_test(self.dat.relative_permeability,
                       {'type': 1, 'parameters': [0.7, 0.0, 1.0, 0.3, None, None, None]})
        self.dict_test(self.dat.capillarity,
                       {'type': 1, 'parameters': [0.0, 0.0, 1.0, None, None, None, None]})
        self.dict_test(self.dat.short_output,
                       {'frequency': 1, 'block': [], 'connection': [], 'generator': []})

        self.tables_test()
        self.write_read_test(tol = 1.e-2)

    def test_AUTOUGH2_2(self):
        """AUTOUGH2 case2"""

        self.base = os.path.join('data', 'AUTOUGH2', '2', 'case2')
        # This data file needs the Fortran read functions because of
        # the formatting in the generation table:
        self.dat = t2data_stats(self.base+'.dat', read_function = fortran_read_function)

        self.assertEqual(self.dat.simulator, 'AUTOUGH2  EWAV')
        self.dict_test(self.dat.parameter,
                       {'print_level': 3, 'max_timesteps': 400,
                        'print_interval': 400, 'max_iterations': None,
                        'option': np.array([0,2,0,0,0,0,0,0,0,0,2,2,2,0,0,0,5,0,0,0,1,0,0,1,0]),
                        'tstart': 0.0, 'tstop': 0.15e17, 'timestep': [0.1e9], 'print_block': 'AR 20',
                        'gravity': 9.8065, 'default_incons': [.10135e06, 20.]})
        self.dict_test(self.dat.multi,
                       {'eos': 'EWAV', 'num_phases': 2, 'num_equations': 3, 'num_components': 2,
                        'num_secondary_parameters': 6})
        self.dict_test(self.dat.relative_permeability,
                       {'type': 1, 'parameters': [0.7, 0.0, 1.0, 0.3, None, None, None]})
        self.dict_test(self.dat.capillarity,
                       {'type': 1, 'parameters': [0.0, 0.0, 1.0, None, None, None, None]})

        self.tables_test()
        self.write_read_test(tol = 1.e-2)

    def test_AUTOUGH2_3(self):
        """AUTOUGH2 extra precision input"""

        self.base = os.path.join('data', 'AUTOUGH2', '3', 'a1')
        self.dat = t2data_stats(self.base+'.dat')

        self.assertEqual(self.dat.simulator, 'AUTOUGH2  EWAV')
        self.dict_test(self.dat.parameter,
                       {'print_level': 3, 'max_timesteps': 9999,
                        'print_interval': 500, 'max_iterations': None,
                        'option': np.array([0,2,0,0,0,0,0,0,0,0,2,2,2,0,0,0,5,0,0,0,1,0,0,1,0]),
                        'tstart': 0.0, 'tstop': 0.15e17, 'timestep': [0.1e9], 'print_block': 'AR 20',
                        'gravity': 9.8065, 'default_incons': [.10135e06, 15.]})
        self.dict_test(self.dat.multi,
                       {'eos': 'EWAV', 'num_phases': 2, 'num_equations': 3, 'num_components': 2,
                        'num_secondary_parameters': 6})
        self.dict_test(self.dat.relative_permeability,
                       {'type': 1, 'parameters': [0.7, 0.0, 1.0, 0.3, None, None, None]})
        self.dict_test(self.dat.capillarity,
                       {'type': 1, 'parameters': [0.0, 0.0, 1.0, None, None, None, None]})

        self.tables_test()
        self.write_read_test()

    def test_TOUGH2_1(self):
        """TOUGH2 r1q, reading grid from separate MESH file"""

        self.base = os.path.join('data', 'TOUGH2', '1', 'r1q')
        self.dat = t2data_stats(self.base,
                                meshfilename = os.path.join('data', 'TOUGH2', '1', 'MESH'))

        self.assertEqual(self.dat.simulator, '')
        expected_meshmaker = [{'radii': [0.0, 100.0]}, {'nequ': 20, 'dr': 10.0},
                              {'nlog': 30, 'rlog': 1.e3}, {'nlog': 30, 'rlog': 1.e4},
                              {'dr': 0.0, 'nequ': 1}, {'layer': [500.0]}]
        for i,section in enumerate(expected_meshmaker):
            self.dict_test(self.dat.meshmaker[0][1][i][1], section)
        self.dict_test(self.dat.parameter,
                       {'print_level': 1, 'max_timesteps': 100,
                        'print_interval': 100, 'max_iterations': None,
                        'option': np.array([0,1,0,0,0,0,0,0,0,0,0,0,0,0,2,0,3,0,0,0,0,3,0,0,0]),
                        'tstart': 0.0, 'tstop': 1.e9, 'timestep': [1.e5], 'print_block': None,
                        'gravity': 0.0, 'default_incons': [60.e5, 0.1]})
        self.dict_test(self.dat.multi, {})
        self.dict_test(self.dat.relative_permeability,
                       {'type': 3, 'parameters': [0.3, 0.05, None, None, None, None, None]})
        self.dict_test(self.dat.capillarity,
                       {'type': 1, 'parameters': [None, 1.0, None, None, None, None, None]})

        self.tables_test()
        self.write_read_test()

    def test_TOUGH2_2(self):
        """TOUGH2 EOS7c"""
        self.base = os.path.join('data', 'TOUGH2', '2', 'eos7c.dat')
        self.dat = t2data(self.base)
        self.assertEqual(len(self.dat.grid.rocktypelist), 5)
        self.assertEqual(self.dat.parameter['default_incons'], [1.e5, 0., 0., 0., 10.5, 15.])
        self.assertEqual(self.dat.selection['integer'][0], 6)
        self.assertEqual(self.dat.selection['float'][40:42], [-330300., 222.])

    def test_TOUGH2_MP_1(self):
        """TOUGH2-MP rfp, reading grid from separate MESHA, MESHB files"""

        d = os.path.join('data', 'TOUGH2-MP', '1')
        self.base = os.path.join(d, 'rfp_nomesh')
        self.dat = t2data_stats(self.base,
                                meshfilename = (os.path.join(d, 'MESHA'),
                                                os.path.join(d, 'MESHB')))

        self.assertEqual(self.dat.simulator, '')
        self.dict_test(self.dat.parameter,
                       {'print_level': 1, 'max_timesteps': 99,
                        'print_interval': 99, 'max_iterations': None,
                        'option': np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0]),
                        'tstart': 0.0, 'tstop': 1.151852E9, 'timestep': [1.e5], 'print_block': ' KA 1',
                        'gravity': 0.0, 'default_incons': [300., 0.01, 5.e5]})
        self.dict_test(self.dat.multi, {})
        self.dict_test(self.dat.relative_permeability,
                       {'type': 3, 'parameters': [0.3, 0.05, None, None, None, None, None]})
        self.dict_test(self.dat.capillarity,
                       {'type': 1, 'parameters': [None, 1.0, None, None, None, None, None]})

        self.tables_test()
        self.write_read_test()

    def test_rectangular_grid_geometry(self):
        """Reverse engineered geometry from rectangular grid"""
        def test_grid(atmos_type = 2, flat = True):
            """Test geometry and associated t2data"""
            dx = np.logspace(0., 1.4, 10)
            dy = np.logspace(0.4, 1.8, 12)
            dz = np.logspace(-1, 1.2, 14)
            def sf(pos): return -0.02 * pos[1] - 0.04 * pos[0]
            geo = mulgrid().rectangular(dx, dy, dz, atmos_type = atmos_type)
            if not flat:
                for col in geo.columnlist:
                    surface = sf(col.centre)
                    col.surface = surface
                    geo.set_column_num_layers(col)
                geo.setup_block_name_index()
                geo.setup_block_connection_name_index()
            rotation_angle = 30.
            geo.rotate(rotation_angle, np.zeros(2))
            geo.permeability_angle = -rotation_angle
            geo.translate(np.array([1000., -2000., 0.]))
            dat = t2data_stats()
            dat.grid = t2grid().fromgeo(geo)
            return geo, dat

        for atmos_type in [0,1,2]:
            for flat in [True, False]:
                geo, dat = test_grid(atmos_type, flat)
                geo1, blockmap = dat.grid.rectgeo(atmos_type = atmos_type, layer_snap = 1.e-3)
                dat1 = t2data_stats()
                dat1.grid = t2grid().fromgeo(geo1)
                dat.compare(dat1, self)

        # Pseudo-PetraSim data file (blocks numbered down columns):
        d = os.path.join('grid', 'rectgeo')
        dat = t2data(os.path.join(d, 'data.dat'))
        geo, mapping = dat.grid.rectgeo(atmos_type = 2)
        self.assertEqual(geo.num_nodes, 16)
        self.assertEqual(geo.num_columns, 9)
        self.assertEqual(geo.num_layers, 4)
        self.assertEqual(geo.atmosphere_type, 2)
        node_pos = np.array([n.pos for n in geo.nodelist])
        expected_pos = np.array([[  0.,   0.], [ 10.,   0.], [ 30.,   0.],
                                 [ 60.,   0.], [  0.,  35.], [ 10.,  35.],
                                 [ 30.,  35.], [ 60.,  35.], [  0.,  60.],
                                 [ 10.,  60.], [ 30.,  60.], [ 60.,  60.],
                                 [  0.,  75.], [ 10.,  75.], [ 30.,  75.],
                                 [ 60.,  75.]])
        self.assertTrue(np.allclose(node_pos, expected_pos))
        layer_centres = np.array([lay.centre for lay in geo.layerlist])
        expected_centres = np.array([ 0. , -0.5, -2. , -4.5])
        self.assertTrue(np.allclose(layer_centres, expected_centres))
        surf = np.array([col.surface for col in geo.columnlist])
        expected_surf = np.array([-0.5, 0., 0., 0., 0., 0., 0., 0., 0])
        self.assertTrue(np.allclose(surf, expected_surf))
        layer_blk_nums = [1, 4, 7, 10, 13, 16, 19, 22, 25]
        blk_nums = layer_blk_nums + \
                   [i+1 for i in layer_blk_nums] + \
                   [i+2 for i in layer_blk_nums]
        block_names = ['%5d' % i for i in blk_nums]
        expected_mapping = dict(zip(geo.block_name_list, block_names))
        self.assertEqual(mapping, expected_mapping)

    def test_rename_blocks(self):

        dx, dy, dz = [10.]*5, [12.]*4, [2.]*4
        geo1 = mulgrid().rectangular(dx, dy, dz, convention = 0)
        geo2 = mulgrid().rectangular(dx, dy, dz, convention = 2)
        dat1 = t2data(); dat2 = t2data()
        dat1.grid = t2grid().fromgeo(geo1)
        dat2.grid = t2grid().fromgeo(geo2)
        inc = [1.e5, 20.]
        dat1.incon = {blk.name: [blk.rocktype.porosity, inc]
                      for blk in dat1.grid.blocklist}
        dat2.incon = {blk.name: [blk.rocktype.porosity, inc]
                      for blk in dat2.grid.blocklist}
        gen1 = t2generator(name = 'gen 1', block = dat1.grid.blocklist[-1].name)
        dat1.add_generator(gen1)
        gen2 = t2generator(name = 'gen 1', block = dat2.grid.blocklist[-1].name)
        dat2.add_generator(gen2)
        dat1.history_blocks = [blk.name for blk in dat1.grid.blocklist[:5]]
        dat2.history_blocks = [blk.name for blk in dat2.grid.blocklist[:5]]
        
        blockmap = {b1.name:b2.name for b1,b2 in zip(dat1.grid.blocklist,
                                                     dat2.grid.blocklist)}
        dat1.rename_blocks(blockmap)

        self.assertEqual([blk.name for blk in dat1.grid.blocklist],
                         [blk.name for blk in dat2.grid.blocklist])
        keys1 = set(dat1.grid.block.keys())
        keys2 = set(dat2.grid.block.keys())
        self.assertEqual(keys1, keys2)

        for blk1, blk2 in zip(dat1.grid.blocklist, dat2.grid.blocklist):
            self.assertEqual(blk1.connection_name, blk2.connection_name)

        keys1 = set(dat1.grid.connection.keys())
        keys2 = set(dat2.grid.connection.keys())
        self.assertEqual(keys1, keys2)

        keys1 = set(dat1.incon.keys())
        keys2 = set(dat2.incon.keys())
        self.assertEqual(keys1, keys2)

        for gen1, gen2 in zip(dat1.generatorlist, dat2.generatorlist):
            self.assertEqual(gen1.block, gen2.block)
        keys1 = set(dat1.generator.keys())
        keys2 = set(dat2.generator.keys())
        self.assertEqual(keys1, keys2)

        self.assertEqual(dat1.history_block, dat2.history_block)

    def test_effective_incons(self):

        nx, ny, nz = 2, 2, 1
        geo = mulgrid().rectangular([10.]*nx, [10.]*ny, [2.]*nz)
        dat = t2data()
        dat.grid = t2grid().fromgeo(geo)

        incs = dat.effective_incons()
        self.assertEqual(incs, [])

        vals = [1.e5, 20.]
        dat.parameter['default_incons'] = vals
        incs = dat.effective_incons()
        self.assertEqual(incs, vals)
        
        vals = [1.e5, 20., None]
        dat.parameter['default_incons'] = vals
        incs = dat.effective_incons()
        self.assertEqual(incs, vals[:-1])

        vals = [20.e5, 0.2, 0.]
        dat.parameter['default_incons'] = vals
        incs = dat.effective_incons()
        self.assertEqual(incs, vals)

        vals = [1.e5, 20.]
        dat.parameter['default_incons'] = vals
        rockname = 'foo 1'
        rt = rocktype(name = rockname)
        dat.grid.add_rocktype(rt)
        for blk in dat.grid.blocklist[-2:]:
            blk.rocktype = rt
        dat.indom = {rockname: [2.e5, 30.]}
        incs = dat.effective_incons()
        expected_p = np.array([1.e5] * 2 + [2.e5] * 2)
        expected_t = np.array([20.] * 2 + [30.] * 2)
        self.assertTrue(np.allclose(incs.variable[:,0], expected_p))
        self.assertTrue(np.allclose(incs.variable[:,1], expected_t))

        dat.indom = {}
        blknames = [blk.name for blk in dat.grid.blocklist[-2:]]
        phi = 0.1
        dat.incon = dict(zip(blknames, [[phi, [2.e5, 30.]]]))
        dat.incon[blknames[0]] = [phi, [2.e5, 30.]]
        dat.incon[blknames[1]] = [phi, [3.e5, 40.]]
        incs = dat.effective_incons()
        self.assertTrue(np.allclose(incs.variable[:,0], np.array([1.e5, 1.e5, 2.e5, 3.e5])))
        self.assertTrue(np.allclose(incs.variable[:,1], np.array([20., 20., 30., 40.])))

        dat.incon = {}
        vals = np.array([
            [1.e5, 20.], [2.e5, 30.], [3.e5, 40.], [5.e5, 60.]])
        incons = dat.grid.incons(vals)
        incs = dat.effective_incons(incons)
        self.assertTrue(np.allclose(incons.variable, incs.variable))

        dat.indom = {rockname: [10.e5, 100.]}
        incs = dat.effective_incons(incons)
        self.assertTrue(np.allclose(incons.variable, incs.variable))

        dat.grid.blocklist[-1].rocktype = dat.grid.rocktype['dfalt']
        incons = t2incon()
        incons[dat.grid.blocklist[1].name] = [2.e5, 30.]
        incons[dat.grid.blocklist[3].name] = [5.e5, 60.]
        incs = dat.effective_incons(incons)
        self.assertTrue(np.allclose(incs.variable[:,0], np.array([1.e5, 2.e5, 10.e5, 5.e5])))
        self.assertTrue(np.allclose(incs.variable[:,1], np.array([20., 30., 100., 60.])))
        self.assertEqual(list(incs.variable[:,0]), [1.e5, 2.e5, 10.e5, 5.e5])

    def test_json(self):

        import json
        nx, ny, nz = 4, 5, 6
        geo = mulgrid().rectangular([10.]*nx, [10.]*ny, [2.]*nz)
        dat = t2data()
        dat.grid = t2grid().fromgeo(geo)

        title = 'test json'
        filename_base = 'test'
        mesh_filename = 'test.exo'
        gravity = 9.80665
        dat.title = title
        dat.filename = filename_base + '.dat'
        dat.parameter['gravity'] = gravity
        dat.multi = {'eos': 'EW'}
        dat.diffusion = [[-1e-6, -1e-6], [-1e-6, -1e-6]]

        def basic_test():
            j = dat.json(geo, mesh_filename)
            self.assertEqual(j['title'], title)
            self.assertEqual(j['gravity'], gravity)
            self.assertEqual(j['thermodynamics'], 'ifc67')
            json.dumps(j)

        def mesh_test():

            # default permeability directions
            geo = mulgrid().rectangular([10.]*nx, [10.]*ny, [2.]*nz)
            dat = t2data()
            dat.grid = t2grid().fromgeo(geo)
            j = dat.mesh_json(geo, mesh_filename)
            self.assertEqual(j['mesh'], {'filename': mesh_filename})
            self.assertFalse('permeability_angle' in j['mesh'])
            self.assertFalse('faces' in j['mesh'])
            json.dumps(j)

            # permeability angle
            geo.permeability_angle = 90.
            dat.grid = t2grid().fromgeo(geo)
            j = dat.mesh_json(geo, mesh_filename)
            self.assertEqual(j['mesh']['permeability_angle'], geo.permeability_angle)
            self.assertFalse('faces' in j['mesh'])
            json.dumps(j)

            # permeability directions swapped 1 <-> 2 in bottom layer
            geo.permeability_angle = 0.
            dat.grid = t2grid().fromgeo(geo)
            new_direction = {1:2, 2:1}
            bottom_layername = geo.layerlist[-1].name
            overridden = []
            for con in dat.grid.connectionlist:
                blknames = [blk.name for blk in con.block]
                if all([geo.layer_name(blkname) == bottom_layername for blkname in blknames]):
                    con.direction = new_direction[con.direction]
                    blk_indices = [geo.block_name_index[blkname] for blkname in blknames]
                    cell_indices = set([i - geo.num_atmosphere_blocks for i in blk_indices])
                    overridden.append(cell_indices)
            j = dat.mesh_json(geo, mesh_filename)
            self.assertEqual(len(j['mesh']['faces']), len(overridden))
            self.assertTrue(all([set(f['cells']) in overridden for f in j['mesh']['faces']]))
            json.dumps(j)

            # rotated mesh and permeability angle, default permeability directions
            angle = 135
            geo.rotate(-angle)
            geo.permeability_angle = angle
            dat.grid = t2grid().fromgeo(geo)
            j = dat.mesh_json(geo, mesh_filename)
            self.assertEqual(j['mesh']['permeability_angle'], angle)
            self.assertFalse('faces' in j['mesh'])
            json.dumps(j)

            # rotated mesh, vertical permeability directions = 1 in bottom layer
            angle = 70
            geo.rotate(-angle)
            geo.permeability_angle = angle
            dat.grid = t2grid().fromgeo(geo)
            overridden = []
            bot_laynames = set([lay.name for lay in geo.layerlist[-2:]])
            for con in dat.grid.connectionlist:
                blknames = [blk.name for blk in con.block]
                laynames = set([geo.layer_name(blkname) for blkname in blknames])
                if laynames == bot_laynames:
                    con.direction = 1
                    blk_indices = [geo.block_name_index[blkname] for blkname in blknames]
                    cell_indices = set([i - geo.num_atmosphere_blocks for i in blk_indices])
                    overridden.append(cell_indices)
            j = dat.mesh_json(geo, mesh_filename)
            self.assertEqual(len(j['mesh']['faces']), len(overridden))
            self.assertTrue(all([f['permeability_direction'] == 1 for f in j['mesh']['faces']]))
            self.assertTrue(all([set(f['cells']) in overridden for f in j['mesh']['faces']]))
            json.dumps(j)

        def eos_test():
            eos = None
            eos_data, tracer_data = dat.eos_json(eos)
            self.assertEqual(eos_data['eos'], {'name': 'we'})
            self.assertIsNone(tracer_data)

            eos = 2
            eos_data, tracer_data = dat.eos_json(eos)
            self.assertEqual(eos_data['eos'], {'name': 'wce'})
            self.assertIsNone(tracer_data)

            eos = 'EWAV'
            eos_data, tracer_data = dat.eos_json(eos)
            self.assertEqual(eos_data['eos'], {'name': 'wae'})
            self.assertIsNone(tracer_data)

            eos = 'EWT'
            eos_data, tracer_data = dat.eos_json(eos)
            self.assertEqual(eos_data['eos'], {'name': 'we'})
            self.assertEqual(tracer_data['tracer'],
                             {'name': 'tracer', 'phase': 'liquid'})

            eos = 'EWTD'
            dat.diffusion = [[-1e-6, -1e-6], [-1e-6, -1e-6]]
            eos_data, tracer_data = dat.eos_json(eos)
            self.assertEqual(eos_data['eos'], {'name': 'we'})
            self.assertEqual(tracer_data['tracer'],
                             {'name': 'tracer', 'phase': 'liquid', 'diffusion': 1e-6})

            dat.diffusion = [[1e-5, 1e-6], [1e-6, 1e-5]]
            with self.assertRaises(Exception):
                eos_data, tracer_data = dat.eos_json(eos)

            eos = 3
            with self.assertRaises(Exception):
                dat.eos_json(eos)

        def output_test():

            dat = t2data()
            dat.filename = filename_base + '.dat'
            dat.parameter['max_timesteps'] = 100

            dat.parameter['print_interval'] = 20
            j = dat.output_json()
            self.assertEqual(j['output']['filename'], filename_base + '.h5')
            self.assertEqual(j['output']['frequency'], dat.parameter['print_interval'])
            self.assertFalse(j['output']['initial'])
            self.assertTrue(j['output']['final'])
            json.dumps(j)

            dat.parameter['print_interval'] = 200
            j = dat.output_json()
            self.assertEqual(j['output']['filename'], filename_base + '.h5')
            self.assertEqual(j['output']['frequency'], 0)
            self.assertFalse(j['output']['initial'])
            self.assertTrue(j['output']['final'])
            json.dumps(j)

            dat.parameter['print_interval'] = 20
            dat.parameter['option'][24] = 1
            j = dat.output_json()
            self.assertTrue(j['output']['initial'])
            json.dumps(j)

            dat.parameter['option'][24] = 0
            dat.output_times = {'time': [0., 20., 135., 400.]}
            j = dat.output_json()
            self.assertEqual(j['output']['checkpoint']['time'], [0., 20., 135., 400.])
            self.assertEqual(j['output']['checkpoint']['tolerance'], 0.)
            if 'repeat' in j['output']['checkpoint']:
                self.assertFalse(j['output']['checkpoint']['repeat'])
            json.dumps(j)

            dat.output_times['num_times'] = 6
            dat.output_times['time_increment'] = 500.
            j = dat.output_json()
            self.assertEqual(j['output']['checkpoint']['time'],
                             [0., 20., 135., 400., 900., 1400.])
            json.dumps(j)

            dat.output_times['num_times_specified'] = 3
            j = dat.output_json()
            self.assertEqual(j['output']['checkpoint']['time'],
                             [0., 20., 135., 635., 1135., 1635.])
            json.dumps(j)
 
            dat.output_times = {'num_times': 100, 'num_times_specified': 1,
                                'time': [10.], 'time_increment': 10.}
            j = dat.output_json()
            self.assertEqual(j['output']['checkpoint']['time'], [10.])
            self.assertEqual(j['output']['checkpoint']['repeat'], 100)
            json.dumps(j)

            dat.output_times = {'num_times': 4, 'num_times_specified': -4,
                                'time': [10., 10., 20., 30.]}
            j = dat.output_json()
            self.assertEqual(j['output']['checkpoint']['step'],
                             [10., 10., 20., 30.])
            self.assertFalse('time' in j['output']['checkpoint'])
            json.dumps(j)

            dat.output_times = {'num_times': 4, 'num_times_specified': -2,
                                'time': [10., 10.], 'time_increment': 30.}
            j = dat.output_json()
            self.assertEqual(j['output']['checkpoint']['step'],
                             [10., 10., 30., 30.])
            self.assertFalse('time' in j['output']['checkpoint'])
            json.dumps(j)

            dat.type = 'AUTOUGH2'
            j = dat.output_json()
            self.assertEqual(j['output']['checkpoint']['tolerance'], 0.1)
            json.dumps(j)

        def timestepping_test():
            method = 'beuler'
            stop_time = 1.e16
            max_timesteps = 1000
            amplification_limit = 5
            amplification = 2.
            reduction_tough2 = 0.25
            reduction_limit = 7
            start_timestep = 1.e6
            dat.parameter['tstop'] = stop_time
            dat.parameter['max_timesteps'] = max_timesteps
            dat.parameter['option'][16] = amplification_limit
            dat.parameter['max_iterations'] = reduction_limit
            dat.parameter['const_timestep'] = start_timestep

            j = dat.timestepping_json()
            self.assertEqual(j['time']['step']['method'], method)
            self.assertEqual(j['time']['stop'], stop_time)
            self.assertEqual(j['time']['step']['maximum']['number'], max_timesteps)
            self.assertEqual(j['time']['step']['adapt']['on'], True)
            self.assertEqual(j['time']['step']['adapt']['minimum'], amplification_limit)
            self.assertEqual(j['time']['step']['adapt']['maximum'], reduction_limit)
            self.assertEqual(j['time']['step']['adapt']['reduction'], reduction_tough2)
            self.assertEqual(j['time']['step']['adapt']['amplification'], amplification)
            self.assertEqual(j['time']['step']['size'], start_timestep)
            json.dumps(j)

            dat.parameter['max_timesteps'] = None
            j = dat.timestepping_json()
            self.assertEqual(j['time']['step']['maximum']['number'], 0)
            json.dumps(j)

            dat.parameter['max_timesteps'] = -1
            j = dat.timestepping_json()
            self.assertIsNone(j['time']['step']['maximum']['number'])
            json.dumps(j)

            dat.type = 'AUTOUGH2'
            reduction_autough2 = 0.2
            j = dat.timestepping_json()
            self.assertEqual(j['time']['step']['adapt']['reduction'], reduction_autough2)
            json.dumps(j)

            # fixed timestep sizes:
            dat.parameter['option'][16] = 0
            dat.parameter['const_timestep'] = -1
            timesteps = [2000., 3000., 5000., 7000.]
            dat.parameter['timestep'] = timesteps
            j = dat.timestepping_json()
            self.assertEqual(j['time']['step']['adapt']['on'], False)
            self.assertEqual(j['time']['step']['size'], timesteps)
            json.dumps(j)

        def relative_permeability_test():

            pars = [0.1, 0.9, 0.2, 0.8]
            dat.relative_permeability = {'type': 1, 'parameters': pars}
            j = dat.relative_permeability_json()
            rp = j['relative_permeability']
            self.assertEqual(rp['type'], 'linear')
            self.assertEqual(rp['liquid'], [pars[0], pars[2]])
            self.assertEqual(rp['vapour'], [pars[1], pars[3]])
            json.dumps(j)

            pars = [1.1, 0., 0., 0.]
            dat.relative_permeability = {'type': 2, 'parameters': pars}
            j = dat.relative_permeability_json()
            rp = j['relative_permeability']
            self.assertEqual(rp['type'], 'pickens')
            self.assertEqual(rp['power'], pars[0])
            json.dumps(j)

            pars = [0.3, 0.05, 0., 0.]
            dat.relative_permeability = {'type': 3, 'parameters': pars}
            j = dat.relative_permeability_json()
            rp = j['relative_permeability']
            self.assertEqual(rp['type'], 'corey')
            self.assertEqual(rp['slr'], pars[0])
            self.assertEqual(rp['ssr'], pars[1])
            json.dumps(j)

            dat.relative_permeability['type'] = 4
            j = dat.relative_permeability_json()
            rp = j['relative_permeability']
            self.assertEqual(rp['type'], 'grant')
            self.assertEqual(rp['slr'], pars[0])
            self.assertEqual(rp['ssr'], pars[1])
            json.dumps(j)

            pars = [0.45, 1.e-3, 1., 0.6]
            dat.relative_permeability = {'type': 7, 'parameters': pars}
            j = dat.relative_permeability_json()
            rp = j['relative_permeability']
            self.assertEqual(rp['type'], 'van Genuchten')
            self.assertEqual(rp['lambda'], pars[0])
            self.assertEqual(rp['slr'], pars[1])
            self.assertEqual(rp['sls'], pars[2])
            self.assertEqual(rp['ssr'], pars[3])
            self.assertFalse(rp['sum_unity'])
            json.dumps(j)

            pars[-1] = 0.
            j = dat.relative_permeability_json()
            rp = j['relative_permeability']
            self.assertTrue(rp['sum_unity'])
            json.dumps(j)

            dat.type = 'AUTOUGH2'
            pars = [0.1, 0.05, 0.8, 0.9, 0.02, 0.8, 0.9]
            dat.relative_permeability = {'type': 19, 'parameters': pars}
            j = dat.relative_permeability_json()
            rp = j['relative_permeability']
            self.assertEqual(rp['type'], 'table')
            self.assertEqual(rp['liquid'], [[0,0], [0.1, 0.02], [0.8, 0.9], [1,1]])
            self.assertEqual(rp['vapour'], [[0,0], [0.05, 0.0], [0.9, 0.8], [1,1]])
            json.dumps(j)

            dat.type = 'TOUGH2'
            self.assertRaises(Exception, dat.relative_permeability_json)

            dat.relative_permeability = {'type': 6, 'parameters': [0., 1., 0., 1.]}
            self.assertRaises(Exception, dat.relative_permeability_json)

        def capillary_pressure_test():

            P0 = 0.1e5
            S = [0.1, 0.9]
            dat.capillarity = {'type': 1, 'parameters': [P0, S[0], S[1]]}
            j = dat.capillary_pressure_json()
            cp = j['capillary_pressure']
            self.assertEqual(cp['type'], 'linear')
            self.assertEqual(cp['pressure'], P0)
            self.assertEqual(cp['saturation_limits'], S)
            json.dumps(j)

            pars = [0.5, 1.e-3, 8.e-5, 10.e5, 1.]
            dat.capillarity = {'type': 7, 'parameters': pars}
            j = dat.capillary_pressure_json()
            cp = j['capillary_pressure']
            self.assertEqual(cp['type'], 'van Genuchten')
            self.assertEqual(cp['lambda'], pars[0])
            self.assertEqual(cp['slr'], pars[1])
            self.assertEqual(cp['P0'], 1. / pars[2])
            self.assertEqual(cp['Pmax'], pars[3])
            self.assertEqual(cp['sls'], pars[4])
            json.dumps(j)

            dat.capillarity = {'type': 8, 'parameters': [0., 0., 0., 0.]}
            j = dat.capillary_pressure_json()
            self.assertIsNone(j['capillary_pressure'])

            dat.capillarity = {'type': 2, 'parameters': [0., 1., 0., 1.]}
            self.assertRaises(Exception, dat.capillary_pressure_json)

        def primary_to_region_test():
            self.assertEqual(primary_to_region_we([2.e5, 20.]), 1)
            self.assertEqual(primary_to_region_we([0.5e5, 100.]), 2)
            self.assertEqual(primary_to_region_we([2.e5, 0.5]), 4)
            self.assertEqual(primary_to_region_wge([2.e5, 20, 0.1e5]), 1)
            self.assertEqual(primary_to_region_wge([2.e5, 100, 1.5e5]), 2)
            self.assertEqual(primary_to_region_wge([1.e5, 0.1, 0.5e5]), 4)

        def initial_test():

            nblks = nx * ny * nz

            eos = 'w'

            incons = [50.e5, 20.]
            j = dat.initial_json(geo, incons, eos)
            self.assertEqual(j['initial']['primary'], incons[:1])
            self.assertEqual(j['initial']['region'], 1)
            json.dumps(j)

            primary1 = [2.e5, 15.]
            primary = [primary1 for i in range(nblks)]
            incons = dat.grid.incons(primary)
            j = dat.initial_json(geo, incons, eos)
            self.assertEqual(j['initial']['primary'], primary1[:1])
            self.assertEqual(j['initial']['region'], 1)
            json.dumps(j)

            primary = np.zeros((nblks, 2))
            primary1 = [0.1e5, 200.]
            primary2 = [0.2e5, 200.]
            n2 = nblks // 2
            primary[:n2] = primary1
            primary[n2:] = primary2
            incons = dat.grid.incons(primary)
            j = dat.initial_json(geo, incons, eos)
            self.assertEqual(len(j['initial']['primary']), nblks)
            self.assertEqual(j['initial']['region'], 2)
            self.assertTrue(all([j['initial']['primary'][i] == primary1[:1]
                                 for i in range(n2)]))
            self.assertTrue(all([j['initial']['primary'][i] == primary2[:1]
                                 for i in range(n2, nblks)]))
            json.dumps(j)

            eos = 'we'

            incons = 'model_ns.h5'
            j = dat.initial_json(geo, incons, eos)
            self.assertEqual(j['initial']['filename'], incons)
            json.dumps(j)

            incons = [3.e5, 35.]
            j = dat.initial_json(geo, incons, eos)
            self.assertEqual(j['initial']['primary'], incons)
            self.assertEqual(j['initial']['region'], 1)
            json.dumps(j)

            incons = [1.e5, 130.]
            j = dat.initial_json(geo, incons, eos)
            self.assertEqual(j['initial']['primary'], incons)
            self.assertEqual(j['initial']['region'], 2)
            json.dumps(j)

            incons = [10.e5, 0.6]
            j = dat.initial_json(geo, incons, eos)
            self.assertEqual(j['initial']['primary'], incons)
            self.assertEqual(j['initial']['region'], 4)
            json.dumps(j)
            
            eos = 'we'
            primary = [2.e5, 15.]
            incons = dat.grid.incons(primary)
            j = dat.initial_json(geo, incons, eos)
            self.assertEqual(j['initial']['primary'], primary)
            self.assertEqual(j['initial']['region'], 1)
            json.dumps(j)

            primary = [0.3e5, 110.]
            incons = dat.grid.incons(primary)
            j = dat.initial_json(geo, incons, eos)
            self.assertEqual(j['initial']['primary'], primary)
            self.assertEqual(j['initial']['region'], 2)
            json.dumps(j)

            primary = np.zeros((nblks, 2))
            primary1 = [2.e5, 15.]
            primary2 = [0.3e5, 110.]
            n2 = nblks // 2
            primary[:n2] = primary1
            primary[n2:] = primary2
            incons = dat.grid.incons(primary)
            j = dat.initial_json(geo, incons, eos)
            self.assertEqual(len(j['initial']['primary']), nblks)
            self.assertEqual(len(j['initial']['region']), nblks)
            self.assertTrue(all([j['initial']['primary'][i] == primary1
                                 for i in range(n2)]))
            self.assertTrue(all([j['initial']['region'][i] == 1
                                 for i in range(n2)]))
            self.assertTrue(all([j['initial']['primary'][i] == primary2
                                 for i in range(n2, nblks)]))
            self.assertTrue(all([j['initial']['region'][i] == 2
                                 for i in range(n2, nblks)]))
            json.dumps(j)

            dat.multi['eos'] = 'EWT'
            primary = [2.e5, 15., 1e-6]
            incons = dat.grid.incons(primary)
            j = dat.initial_json(geo, incons, 'we', {'name': 'tracer'})
            self.assertEqual(j['initial']['primary'], primary[:2])
            self.assertEqual(j['initial']['region'], 1)
            self.assertEqual(j['initial']['tracer'], 1e-6)
            json.dumps(j)

        def generators_test():

            dat.parameter['option'][12] = 0

            def generator_json(gen, eos = 'we', tracer = None):
                dat.clear_generators()
                dat.add_generator(gen)
                j = dat.generators_json(geo, eos, tracer)
                self.assertEqual(len(j['source']), 1)
                return j['source'][0]

            # mass production
            q = -10.
            blkname = '  c 3'
            name = 'gen 1'
            gen = t2generator(name = name, block = blkname,
                              type = 'MASS', gx = q)
            g = generator_json(gen)
            self.assertEqual(g['rate'], q)
            self.assertFalse('component' in g)
            i = geo.block_name_index[blkname] - geo.num_atmosphere_blocks
            self.assertEqual(g['cell'], i)
            self.assertEqual(g['name'], name)
            json.dumps(g)

            # mass injection
            q = 5.
            gen = t2generator(name = name, block = blkname,
                              type = 'MASS', gx = q)
            g = generator_json(gen)
            self.assertEqual(g['rate'], q)
            self.assertEqual(g['component'], 1)
            json.dumps(g)

            # COM2
            gen = t2generator(name = name, block = blkname,
                              type = 'COM2', gx = q)
            g = generator_json(gen)
            self.assertEqual(g['rate'], q)
            self.assertEqual(g['component'], 2)
            json.dumps(g)

            # COM2 tracer
            gen = t2generator(name = name, block = blkname,
                              type = 'COM2', gx = q)
            g = generator_json(gen, eos = 'we', tracer = {'name': 'foo'})
            self.assertFalse('rate' in g)
            self.assertEqual(g['tracer'], q)
            self.assertFalse('component' in g)
            json.dumps(g)

            # TRAC
            gen = t2generator(name = name, block = blkname,
                              type = 'TRAC', gx = q)
            g = generator_json(gen, eos = 'we', tracer = {'name': 'foo'})
            self.assertFalse('rate' in g)
            self.assertEqual(g['tracer'], q)
            self.assertFalse('component' in g)
            json.dumps(g)

            # TRAC table
            t = [0., 10., 240., 350.]
            q = [2.e-6, 2.5e-6, 2.9e-6, 3.1e-6]
            gen = t2generator(name = name, block = blkname,
                              type = 'TRAC', time = t, rate = q)
            g = generator_json(gen, eos = 'we', tracer = {'name': 'foo'})
            self.assertEqual(g['tracer'], [list(r) for r in zip(t, q)])
            self.assertFalse('rate' in g)
            self.assertFalse('component' in g)
            self.assertFalse('enthalpy' in g)
            json.dumps(g)

            # heat
            q = 1000.
            gen = t2generator(name = name, block = blkname,
                              type = 'HEAT', gx = q)
            g = generator_json(gen)
            self.assertEqual(g['rate'], q)
            self.assertEqual(g['component'], 2)
            json.dumps(g)

            # heat, 3-eqn EOS
            g = generator_json(gen, 'wce')
            self.assertEqual(g['rate'], q)
            self.assertEqual(g['component'], 3)
            json.dumps(g)

            # mass production table
            t = [0., 10., 240., 350.]
            q = [-2., -2.5, -2.9, -3.1]
            gen = t2generator(name = name, block = blkname,
                              type = 'MASS', time = t, rate = q)
            g = generator_json(gen)
            self.assertEqual(g['rate'], [list(r) for r in zip(t, q)])
            self.assertFalse('component' in g)
            json.dumps(g)

            # mass injection table
            q = [2., 2.5, 2.9, 3.1]
            h = 350.e3
            gen = t2generator(name = name, block = blkname,
                              type = 'MASS', ex = h, time = t, rate = q)
            g = generator_json(gen)
            self.assertEqual(g['rate'], [list(r) for r in zip(t, q)])
            self.assertEqual(g['enthalpy'], h)
            self.assertEqual(g['component'], 1)
            json.dumps(g)

            # mass injection enthalpy table
            h = [300.e3, 330.e3, 350.e3, 335.e3]
            gen = t2generator(name = name, block = blkname,
                              type = 'MASS', ex = h,
                              time = t, rate = q, enthalpy = h)
            g = generator_json(gen)
            self.assertEqual(g['rate'], [list(r) for r in zip(t, q)])
            self.assertEqual(g['enthalpy'], [list(r) for r in zip(t, h)])
            self.assertEqual(g['component'], 1)
            json.dumps(g)

            # deliverability
            PI = 1.e-12
            Pwb = 2.e5
            gen = t2generator(name = name, block = blkname,
                              type = 'DELV', gx = PI, ex = Pwb)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'], {'pressure': Pwb, 'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertEqual(g['separator'], {'pressure': default_separator_pressure})
            json.dumps(g)

            # DELG
            gen = t2generator(name = name, block = blkname,
                              type = 'DELG', gx = PI, ex = Pwb)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'], {'pressure': Pwb, 'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertFalse('limiter' in g)
            self.assertEqual(g['separator'], {'pressure': default_separator_pressure})
            json.dumps(g)

            # DELG with initial rate
            q0 = -3.
            gen = t2generator(name = name, block = blkname,
                              type = 'DELG', ex = Pwb, hg = q0)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'], {'pressure': Pwb})
            self.assertEqual(g['direction'], 'production')
            self.assertEqual(g['rate'], q0)
            self.assertFalse('limiter' in g)
            self.assertEqual(g['separator'], {'pressure': default_separator_pressure})
            json.dumps(g)

            # DELG with steam limiter
            qmax = 10.
            gen = t2generator(name = name, block = blkname,
                              type = 'DELG', gx = PI, ex = Pwb, hg = qmax)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'], {'pressure': Pwb, 'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertEqual(g['limiter'], {'steam': qmax})
            self.assertEqual(g['separator'], {'pressure': default_separator_pressure})
            json.dumps(g)

            # DELG with steam limiter and separator pressure
            qmax = 5.
            psep = 10.e5
            gen = t2generator(name = name, block = blkname,
                              type = 'DELG', gx = PI, ex = Pwb, fg = psep, hg = qmax)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'], {'pressure': Pwb, 'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertEqual(g['limiter'], {'steam': qmax})
            self.assertEqual(g['separator'], {'pressure': psep})
            json.dumps(g)

            # DELG with steam limiter and 2-stage separator
            qmax = 5.
            gen = t2generator(name = name, block = blkname,
                              type = 'DELG', gx = PI, ex = Pwb, fg = -1, hg = qmax)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'], {'pressure': Pwb, 'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertEqual(g['limiter'], {'steam': qmax})
            self.assertEqual(g['separator'], {'pressure': default_2_stage_separator_pressure})
            json.dumps(g)

            # DELG with table of PI vs time
            t = [0., 10., 240., 350., 750.]
            PI = [1.e-11, 8.e-12, 7.e-12, 4.e-12, 3.5e-12]
            gen = t2generator(name = name, block = blkname,
                              type = 'DELG', ex = Pwb, ltab = 1,
                              time = t, rate = PI)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'],
                             {'pressure': Pwb, 'productivity':
                              {'time': [list(r) for r in zip(t, PI)]}})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertFalse('limiter' in g)
            self.assertEqual(g['separator'], {'pressure': default_separator_pressure})
            json.dumps(g)

            # DELG with table of cutoff pressure vs enthalpy
            h = [100.e3, 300.e3, 600.e3, 900.e3]
            Pwb = [5.e5, 4.e5, 2.e5, 1.5e5]
            PI = 2.e-12
            gen = t2generator(name = name, block = blkname,
                              type = 'DELG', ltab = 0, gx = PI,
                              time = h, rate = Pwb)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'],
                             {'pressure': {'enthalpy': [list(r) for r in zip(h, Pwb)]},
                              'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertFalse('limiter' in g)
            self.assertEqual(g['separator'], {'pressure': default_separator_pressure})
            json.dumps(g)

            # DELG with table of cutoff pressure vs enthalpy and steam limiter
            qmax = 5.
            gen = t2generator(name = name, block = blkname,
                              type = 'DELG', ltab = 0, gx = PI,
                              time = h, rate = Pwb, hg = qmax)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'],
                             {'pressure': {'enthalpy': [list(r) for r in zip(h, Pwb)]},
                              'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertEqual(g['limiter'], {'steam': qmax})
            self.assertEqual(g['separator'], {'pressure': default_separator_pressure})
            json.dumps(g)

            # DELS
            PI = 1.e-13
            Pwb = 5.e5
            gen = t2generator(name = name, block = blkname,
                              type = 'DELS', gx = PI, ex = Pwb)
            g = generator_json(gen)
            self.assertEqual(g['production_component'], 2)
            self.assertEqual(g['deliverability'], {'pressure': Pwb, 'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertFalse('limiter' in g)
            json.dumps(g)

            # DELS with steam limiter and separator pressure
            qmax = 5.
            psep = 6.e5
            gen = t2generator(name = name, block = blkname,
                              type = 'DELS', gx = PI, ex = Pwb, fg = psep, hg = qmax)
            g = generator_json(gen)
            self.assertEqual(g['production_component'], 2)
            self.assertEqual(g['deliverability'], {'pressure': Pwb, 'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertEqual(g['limiter'], {'steam': qmax})
            self.assertEqual(g['separator'], {'pressure': psep})
            json.dumps(g)

            # DELT with total flow limiter
            PI = 1.e-13
            Pwb = 5.e5
            qmax = 20.
            gen = t2generator(name = name, block = blkname,
                              type = 'DELT', gx = PI, ex = Pwb, hg = qmax)
            g = generator_json(gen)
            self.assertFalse('production_component' in g)
            self.assertEqual(g['deliverability'], {'pressure': Pwb, 'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertEqual(g['limiter'], {'total': qmax})
            self.assertEqual(g['separator'], {'pressure': default_separator_pressure})
            json.dumps(g)

            # DELT with negative limit specified
            PI = 1.e-13
            Pwb = 5.e5
            qmax = -1.
            gen = t2generator(name = name, block = blkname,
                              type = 'DELT', gx = PI, ex = Pwb, hg = qmax)
            g = generator_json(gen)
            self.assertFalse('production_component' in g)
            self.assertEqual(g['deliverability'], {'pressure': Pwb, 'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertFalse('limiter' in g)
            self.assertEqual(g['separator'], {'pressure': default_separator_pressure})
            json.dumps(g)

            # DELT with table of PI vs time
            t = [0., 10., 240., 350., 750.]
            PI = [1.e-11, 8.e-12, 7.e-12, 4.e-12, 3.5e-12]
            gen = t2generator(name = name, block = blkname,
                              type = 'DELT', ex = Pwb, ltab = 1,
                              time = t, rate = PI)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'],
                             {'pressure': Pwb, 'productivity':
                              {'time': [list(r) for r in zip(t, PI)]}})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertFalse('limiter' in g)
            self.assertEqual(g['separator'], {'pressure': default_separator_pressure})
            json.dumps(g)

            # DELT with table of cutoff pressure vs enthalpy and total flow limiter
            PI = 2.e-12
            h = [100.e3, 300.e3, 600.e3, 900.e3]
            Pwb = [5.e5, 4.e5, 2.e5, 1.5e5]
            qmax = 5.
            gen = t2generator(name = name, block = blkname,
                              type = 'DELT', ltab = 0, gx = PI,
                              time = h, rate = Pwb, hg = qmax)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'],
                             {'pressure': {'enthalpy': [list(r) for r in zip(h, Pwb)]},
                              'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertEqual(g['limiter'], {'total': qmax})
            self.assertEqual(g['separator'], {'pressure': default_separator_pressure})
            json.dumps(g)

            # DELW with liquid flow limiter
            PI = 1.e-13
            Pwb = 5.e5
            qmax = 10.
            gen = t2generator(name = name, block = blkname,
                              type = 'DELW', gx = PI, ex = Pwb, hg = qmax)
            g = generator_json(gen)
            self.assertFalse('production_component' in g)
            self.assertEqual(g['deliverability'], {'pressure': Pwb, 'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertEqual(g['limiter'], {'water': qmax})
            json.dumps(g)

            # DELW with negative limit specified
            PI = 1.e-13
            Pwb = 5.e5
            qmax = -1.
            gen = t2generator(name = name, block = blkname,
                              type = 'DELW', gx = PI, ex = Pwb, hg = qmax)
            g = generator_json(gen)
            self.assertFalse('production_component' in g)
            self.assertEqual(g['deliverability'], {'pressure': Pwb, 'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertFalse('limiter' in g)
            self.assertEqual(g['separator'], {'pressure': default_separator_pressure})
            json.dumps(g)

            # DELW with table of PI vs time
            t = [0., 10., 240., 350., 750.]
            PI = [1.e-11, 8.e-12, 7.e-12, 4.e-12, 3.5e-12]
            gen = t2generator(name = name, block = blkname,
                              type = 'DELW', ex = Pwb, ltab = 1,
                              time = t, rate = PI)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'],
                             {'pressure': Pwb, 'productivity':
                              {'time': [list(r) for r in zip(t, PI)]}})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertFalse('limiter' in g)
            self.assertEqual(g['separator'], {'pressure': default_separator_pressure})
            json.dumps(g)

            # DELW with table of cutoff pressure vs enthalpy and liquid flow limiter
            PI = 2.e-12
            h = [100.e3, 300.e3, 600.e3, 900.e3]
            Pwb = [5.e5, 4.e5, 2.e5, 1.5e5]
            qmax = 5.
            gen = t2generator(name = name, block = blkname,
                              type = 'DELW', ltab = 0, gx = PI,
                              time = h, rate = Pwb, hg = qmax)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'],
                             {'pressure': {'enthalpy': [list(r) for r in zip(h, Pwb)]},
                              'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertEqual(g['limiter'], {'water': qmax})
            self.assertEqual(g['separator'], {'pressure': default_separator_pressure})
            json.dumps(g)

            # RECH with specified mass flow (effectively same as MASS)
            q = 3.5
            h = 400.e3
            gen = t2generator(name = name, block = blkname,
                              type = 'RECH', gx = q, ex = h, hg = 0)
            g = generator_json(gen)
            self.assertEqual(g['rate'], q)
            self.assertEqual(g['enthalpy'], h)
            json.dumps(g)

            # RECH with reference pressure from initial conditions
            coef = 1.e-6
            gen = t2generator(name = name, block = blkname,
                              type = 'RECH', gx = coef, ex = h, hg = -1)
            g = generator_json(gen)
            self.assertFalse('rate' in g)
            self.assertEqual(g['recharge'], {'coefficient': coef, 'pressure': 'initial'})
            self.assertEqual(g['enthalpy'], h)
            self.assertEqual(g['direction'], 'both')
            json.dumps(g)

            # RECH with specified reference pressure
            coef = 1.e-6
            P0 = 2.e5
            gen = t2generator(name = name, block = blkname,
                              type = 'RECH', gx = coef, ex = h, hg = P0)
            g = generator_json(gen)
            self.assertFalse('rate' in g)
            self.assertEqual(g['recharge'], {'coefficient': coef, 'pressure': P0})
            self.assertEqual(g['enthalpy'], h)
            self.assertEqual(g['direction'], 'both')
            json.dumps(g)

            # RECH with specified reference pressure, no inflow
            gen = t2generator(name = name, block = blkname,
                              type = 'RECH', gx = coef, ex = h, hg = P0, fg = -1)
            g = generator_json(gen)
            self.assertFalse('rate' in g)
            self.assertEqual(g['recharge'], {'coefficient': coef, 'pressure': P0})
            self.assertEqual(g['enthalpy'], h)
            self.assertEqual(g['direction'], 'out')
            json.dumps(g)

            # RECH with specified reference pressure, no outflow
            gen = t2generator(name = name, block = blkname,
                              type = 'RECH', gx = coef, ex = h, hg = P0, fg = 1)
            g = generator_json(gen)
            self.assertFalse('rate' in g)
            self.assertEqual(g['recharge'], {'coefficient': coef, 'pressure': P0})
            self.assertEqual(g['enthalpy'], h)
            self.assertEqual(g['direction'], 'in')
            json.dumps(g)

            # MASD
            q = -10.
            PI = 1.e-12
            Pwb = 2.5e5
            Pthreshold = 3.e5
            gen = t2generator(name = name, block = blkname,
                              type = 'MASD', gx = q, ex = PI, fg = Pwb, hg = Pthreshold)
            g = generator_json(gen)
            self.assertEqual(g['rate'], q)
            self.assertEqual(g['deliverability'],
                             {'pressure': Pwb, 'productivity': PI, 'threshold': Pthreshold})
            self.assertEqual(g['direction'], 'production')
            self.assertEqual(g['limiter'], {'total': abs(q)})
            self.assertEqual(g['separator'], {'pressure': Pwb})
            json.dumps(g)

            # MASD with mass table
            t = [0., 10., 240., 350., 750.]
            q = [-10., -12., -14., -13., -16.]
            gx = 16
            PI = 1.e-12
            Pwb = 2.5e5
            Pthreshold = 3.e5
            gen = t2generator(name = name, block = blkname,
                              type = 'MASD', ex = PI, fg = Pwb, hg = Pthreshold,
                              time = t, rate = q, gx = gx)
            g = generator_json(gen)
            self.assertEqual(g['rate'], [list(r) for r in zip(t, q)])
            self.assertEqual(g['deliverability'],
                             {'pressure': Pwb, 'productivity': PI, 'threshold': Pthreshold})
            self.assertEqual(g['direction'], 'production')
            self.assertEqual(g['limiter'], {'total': gx})
            self.assertEqual(g['separator'], {'pressure': Pwb})
            json.dumps(g)

            # DMAK
            PI = 1.e-12
            Pwb = 2.e5
            qmax = 10.
            gen = t2generator(name = name, block = blkname,
                              type = 'DMAK', gx = PI, ex = Pwb, hg = qmax)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'], {'pressure': Pwb, 'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertEqual(g['limiter'], {'steam': qmax})
            self.assertEqual(g['separator'], {'pressure': default_separator_pressure})
            json.dumps(g)

            # DMAK with 2-stage separator
            gen = t2generator(name = name, block = blkname,
                              type = 'DMAK', gx = PI, ex = Pwb, fg = -1, hg = qmax)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'], {'pressure': Pwb, 'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertEqual(g['limiter'], {'steam': qmax})
            self.assertEqual(g['separator'], {'pressure': default_2_stage_separator_pressure})
            json.dumps(g)

            # DMAT
            PI = 1.e-12
            Pwb = 2.e5
            qmax = 20.
            psep = 65.e5
            gen = t2generator(name = name, block = blkname,
                              type = 'DMAT', gx = PI, ex = Pwb, hg = qmax, fg = psep)
            g = generator_json(gen)
            self.assertEqual(g['deliverability'], {'pressure': Pwb, 'productivity': PI})
            self.assertEqual(g['direction'], 'production')
            self.assertFalse('rate' in g)
            self.assertEqual(g['limiter'], {'total': qmax})
            self.assertEqual(g['separator'], {'pressure': psep})
            json.dumps(g)

            # IMAK
            inj = 1.e-4
            P0 = 45.e5
            qmax = 20.
            h = 4.4e5
            gen = t2generator(name = name, block = blkname,
                              type = 'IMAK', gx = qmax, ex = h, hg = P0, fg = inj)
            g = generator_json(gen)
            self.assertFalse('rate' in g)
            self.assertEqual(g['injectivity'], {'coefficient': inj, 'pressure': P0})
            self.assertFalse('enthalpy' in g) # enthalpy defined in reinjector output
            self.assertEqual(g['direction'], 'injection')
            self.assertEqual(g['limiter'], {'total': qmax})
            json.dumps(g)

            # XINJ (with no limiter)
            inj = 1.e-4
            P0 = 45.e5
            h = 4.4e5
            gen = t2generator(name = name, block = blkname,
                              type = 'XINJ', ex = h, hg = P0, fg = inj)
            g = generator_json(gen)
            self.assertFalse('rate' in g)
            self.assertEqual(g['injectivity'], {'coefficient': inj, 'pressure': P0})
            self.assertEqual(g['enthalpy'], h)
            self.assertEqual(g['direction'], 'injection')
            self.assertFalse('limiter' in g)
            json.dumps(g)

        def network_test():

            dat.clear_generators()
            gen = t2generator(name = 'abc 1', block = '  a 1', type = 'DELG')
            dat.add_generator(gen)
            gen = t2generator(name = 'abc 2', block = '  a 2', type = 'DELG')
            dat.add_generator(gen)
            j = dat.generators_json(geo, 'we')
            self.assertFalse('network' in j)

            # reset reinjection
            gen = t2generator(name = 'chk 1', block = 'chk99', type = 'FINJ',
                              gx = 1e5, ex = 85e3, hg = 1., fg = 1.)
            dat.add_generator(gen)
            # add two makeup wells
            gen = t2generator(name = 'foo 1', block = '  a 1', type = 'DMAK')
            dat.add_generator(gen)
            gen = t2generator(name = 'foo 2', block = '  a 2', type = 'DMAK')
            dat.add_generator(gen)
            j = dat.generators_json(geo, 'we')
            self.assertTrue('network' in j)
            self.assertEqual(len(j['network']['group']), 1)
            grp = j['network']['group'][-1]
            self.assertEqual(grp['name'], 'reinjector group 1')
            self.assertEqual(grp['in'], ['abc 1', 'abc 2'])
            self.assertFalse('scaling' in grp)
            self.assertFalse('limiter' in grp)

            # add TMAK
            gen = t2generator(block = '  a 1', type = 'TMAK',
                              gx = 100., hg = -1)
            dat.add_generator(gen)
            j = dat.generators_json(geo, 'we')
            self.assertEqual(len(j['source']), 5)
            self.assertEqual(len(j['network']['group']), 2)
            grp = j['network']['group'][-1]
            self.assertEqual(grp['name'], 'makeup 1')
            self.assertEqual(grp['in'], ['foo 1', 'foo 2'])
            self.assertEqual(grp['scaling'], 'uniform')
            self.assertEqual(grp['limiter'], {'total': 100})

            # add two water reinjection wells
            q1, h1 = 1.1, 85.e3
            gen = t2generator(name = 'inj 1', block = '  b 1', type = 'FINJ',
                              gx = q1, ex = h1, hg = 1.)
            dat.add_generator(gen)
            f2, h2 = 0.3, 90.e3
            gen = t2generator(name = 'inj 2', block = '  b 2', type = 'PINJ',
                              ex = h2, hg = f2)
            dat.add_generator(gen)
            j = dat.generators_json(geo, 'we')
            self.assertEqual(len(j['source']), 7)
            self.assertEqual(len(j['network']['group']), 2)
            self.assertEqual(len(j['network']['reinject']), 2)
            r = j['network']['reinject'][-1]
            self.assertEqual(r['name'], 'reinjector 2')
            self.assertEqual(r['in'], 'makeup 1')
            self.assertEqual(len(r['water']), 2)
            self.assertEqual(len(r['steam']), 0)
            self.assertEqual(r['water'][0], {'out': 'inj 1', 'rate': q1, 'enthalpy': h1})
            self.assertEqual(r['water'][1], {'out': 'inj 2', 'proportion': f2, 'enthalpy': h2})
            self.assertFalse('overflow' in r)

            # add two RINJ reinjection wells
            h3, f3 = 82e3, 0.2
            gen = t2generator(name = 'inj 3', block = '  b 3', type = 'RINJ',
                              ex = h3, hg = f3)
            dat.add_generator(gen)
            h4, f4 = 77.e3, 0.35
            gen = t2generator(name = 'inj 4', block = '  b 4', type = 'RINJ',
                              ex = h4, hg = f4, fg = 1.)
            dat.add_generator(gen)
            j = dat.generators_json(geo, 'we')
            self.assertEqual(len(j['source']), 9)
            self.assertEqual(len(j['network']['group']), 2)
            self.assertEqual(len(j['network']['reinject']), 3)
            r = j['network']['reinject'][1]
            self.assertEqual(r['overflow'], 'reinjector 3')
            r = j['network']['reinject'][2]
            self.assertEqual(r['name'], 'reinjector 3')
            self.assertFalse('in' in r)
            self.assertEqual(len(r['water']), 2)
            self.assertEqual(len(r['steam']), 0)
            self.assertEqual(r['water'][0], {'out': 'inj 3', 'proportion': f3, 'enthalpy': h3})
            self.assertEqual(r['water'][1], {'out': 'inj 4', 'proportion': f4, 'enthalpy': h4})

            # add more production wells and a second TMAK
            gen = t2generator(name = 'foo 3', block = '  a 3', type = 'DELG')
            dat.add_generator(gen)
            gen = t2generator(name = 'foo 4', block = '  a 4', type = 'DMAK')
            dat.add_generator(gen)
            gen = t2generator(name = 'foo 5', block = '  a 5', type = 'DELG')
            dat.add_generator(gen)
            gen = t2generator(name = 'foo 6', block = '  a 6', type = 'DMAT')
            dat.add_generator(gen)
            gen = t2generator(name = 'tmk 2', block = '  a 1', type = 'TMAK',
                              gx = 50., ex = 20, hg = -2)
            dat.add_generator(gen)
            j = dat.generators_json(geo, 'we')
            self.assertEqual(len(j['network']['group']), 3)
            grp = j['network']['group'][2]
            self.assertEqual(grp['name'], 'tmk 2')
            self.assertEqual(grp['in'], ['foo 4', 'foo 6'])
            self.assertEqual(grp['scaling'], 'progressive')
            self.assertEqual(grp['limiter'], {'total': 50, 'steam': 20})

            # add three more reinjection wells including one IMAK
            h5, q5 = 87e3, 1.5
            gen = t2generator(name = 'inj 5', block = '  b 5', type = 'FINJ',
                              gx = q5, ex = h5, hg = 1.)
            dat.add_generator(gen)
            q6, h6, P6, finj6 = 10., 1200.e3, 5e5, 1e-7
            gen = t2generator(name = 'inj 6', block = '  b 6', type = 'IMAK',
                              gx = q6, ex = h6, hg = P6, fg = -finj6)
            dat.add_generator(gen)
            f7, h7 = 0.1, 1400.e3
            gen = t2generator(name = 'inj 7', block = '  c 1', type = 'PINJ',
                              ex = h7, hg = -f7, fg = 1.)
            dat.add_generator(gen)
            j = dat.generators_json(geo, 'we')
            q = j['source'][-2]
            self.assertEqual(q['direction'], 'injection')
            self.assertEqual(q['limiter'], {'total': q6})
            self.assertEqual(q['injectivity'], {'pressure': P6, 'coefficient': finj6})
            self.assertEqual(len(j['network']['group']), 4)
            grp = j['network']['group'][3]
            self.assertEqual(grp['name'], 'reinjector group 4')
            self.assertEqual(grp['in'], ['foo 3', 'foo 5', 'tmk 2'])
            self.assertFalse('scaling' in grp)
            self.assertFalse('limiter' in grp)
            self.assertEqual(len(j['source']), 16)
            self.assertEqual(len(j['network']['reinject']), 4)
            r = j['network']['reinject'][3]
            self.assertEqual(r['name'], 'reinjector 4')
            self.assertEqual(r['in'], 'reinjector group 4')
            self.assertEqual(len(r['water']), 1)
            self.assertEqual(len(r['steam']), 2)
            self.assertEqual(r['water'][0], {'out': 'inj 5', 'rate': q5, 'enthalpy': h5})
            self.assertEqual(r['steam'][0], {'out': 'inj 6', 'enthalpy': h6})
            self.assertEqual(r['steam'][1], {'out': 'inj 7', 'proportion': f7, 'enthalpy': h7})
            self.assertFalse('overflow' in r)

            # two TMAKs contributing to a single reinjector
            dat.clear_generators()
            gen = t2generator(name = 'foo 1', block = '  a 4', type = 'DMAK')
            dat.add_generator(gen)
            gen = t2generator(name = 'foo 2', block = '  a 5', type = 'DMAK')
            dat.add_generator(gen)
            gen = t2generator(name = 'tmk 1', block = '  a 1', type = 'TMAK',
                              gx = 50, ex = 20, hg = -2)
            dat.add_generator(gen)
            gen = t2generator(name = 'foo 3', block = '  a 6', type = 'DMAK')
            dat.add_generator(gen)
            gen = t2generator(name = 'foo 4', block = '  a 3', type = 'DMAK')
            dat.add_generator(gen)
            gen = t2generator(name = 'tmk 2', block = '  a 2', type = 'TMAK',
                              gx = 60, ex = 30, hg = -2)
            dat.add_generator(gen)
            h1, q1 = 87e3, 1.5
            gen = t2generator(name = 'inj 1', block = '  b 1', type = 'FINJ',
                              gx = q1, ex = h1, hg = 1.)
            dat.add_generator(gen)
            f2, h2 = 0.3, 90.e3
            gen = t2generator(name = 'inj 2', block = '  c 1', type = 'PINJ',
                              ex = h2, hg = f2, fg = 1.)
            dat.add_generator(gen)
            j = dat.generators_json(geo, 'we')
            self.assertEqual(len(j['source']), 6)
            self.assertEqual(len(j['network']['group']), 3)
            self.assertEqual(len(j['network']['reinject']), 1)
            grp = j['network']['group'][0]
            self.assertEqual(grp['name'], 'tmk 1')
            self.assertEqual(grp['in'], ['foo 1', 'foo 2'])
            self.assertEqual(grp['scaling'], 'progressive')
            self.assertEqual(grp['limiter'], {'total': 50, 'steam': 20})
            grp = j['network']['group'][1]
            self.assertEqual(grp['name'], 'tmk 2')
            self.assertEqual(grp['in'], ['foo 3', 'foo 4'])
            self.assertEqual(grp['scaling'], 'progressive')
            self.assertEqual(grp['limiter'], {'total': 60, 'steam': 30})
            grp = j['network']['group'][2]
            self.assertEqual(grp['name'], 'reinjector group 1')
            self.assertEqual(grp['in'], ['tmk 1', 'tmk 2'])
            self.assertFalse('scaling' in grp)
            r = j['network']['reinject'][0]
            self.assertEqual(r['name'], 'reinjector 1')
            self.assertEqual(r['in'], 'reinjector group 1')
            self.assertEqual(len(r['water']), 2)
            self.assertFalse('steam' in r)
            self.assertFalse('overflow' in r)

        def boundaries_test():

            nx, ny, nz = 2, 2, 3
            dx, dy, dz = 100., 100., 10.
            geo = mulgrid().rectangular([dx]*nx, [dy]*ny, [dz]*nz, atmos_type = 1)
            dat = t2data()
            dat.grid = t2grid().fromgeo(geo)
            atmos_volume = 1.e25
            mesh_coords = 'xyz'

            # isothermal water, liquid top BCs
            eos = 'w'
            P0, T0 = 1.e5, 15.
            bdy_incons = dat.grid.incons((P0, T0))
            j = dat.boundaries_json(geo, bdy_incons, atmos_volume, eos, mesh_coords)
            self.assertEqual(len(j['boundaries']), 1)
            self.assertEqual(j['boundaries'][0]['primary'], [P0])
            self.assertEqual(j['boundaries'][0]['region'], 1)
            self.assertFalse('tracer' in j['boundaries'][0])
            self.assertEqual(j['boundaries'][0]['faces']['normal'], [0, 0, 1])
            self.assertEqual(j['boundaries'][0]['faces']['cells'], [0, 1, 2, 3])
            json.dumps(j)

            # pure water, liquid top BCs
            eos = 'we'
            P0, T0 = 1.e5, 15.
            bdy_incons = dat.grid.incons((P0, T0))
            j = dat.boundaries_json(geo, bdy_incons, atmos_volume, eos, mesh_coords)
            self.assertEqual(len(j['boundaries']), 1)
            self.assertEqual(j['boundaries'][0]['primary'], [P0, T0])
            self.assertEqual(j['boundaries'][0]['region'], 1)
            self.assertFalse('tracer' in j['boundaries'][0])
            self.assertEqual(j['boundaries'][0]['faces']['normal'], [0, 0, 1])
            self.assertEqual(j['boundaries'][0]['faces']['cells'], [0, 1, 2, 3])
            json.dumps(j)

            # pure water + tracer, liquid top BCs
            dat.multi['eos'] = 'EWT'
            eos = 'we'
            P0, T0, X0 = 1.e5, 15., 1.e-6
            bdy_incons = dat.grid.incons((P0, T0, X0))
            j = dat.boundaries_json(geo, bdy_incons, atmos_volume, eos, mesh_coords,
                                    {'name': 'tracer'})
            self.assertEqual(len(j['boundaries']), 1)
            self.assertEqual(j['boundaries'][0]['primary'], [P0, T0])
            self.assertEqual(j['boundaries'][0]['region'], 1)
            self.assertEqual(j['boundaries'][0]['tracer'], X0)
            self.assertEqual(j['boundaries'][0]['faces']['normal'], [0, 0, 1])
            self.assertEqual(j['boundaries'][0]['faces']['cells'], [0, 1, 2, 3])
            json.dumps(j)
            dat.multi['eos'] = 'EW'

            # pure water, dry steam top BCs
            eos = 'we'
            P0, T0 = 0.8e5, 100.
            bdy_incons = dat.grid.incons((P0, T0))
            j = dat.boundaries_json(geo, bdy_incons, atmos_volume, eos, mesh_coords)
            self.assertEqual(len(j['boundaries']), 1)
            self.assertEqual(j['boundaries'][0]['primary'], [P0, T0])
            self.assertEqual(j['boundaries'][0]['region'], 2)
            self.assertFalse('tracer' in j['boundaries'][0])
            self.assertEqual(j['boundaries'][0]['faces']['normal'], [0, 0, 1])
            self.assertEqual(j['boundaries'][0]['faces']['cells'], [0, 1, 2, 3])
            json.dumps(j)

            # pure water, two-phase top BCs
            eos = 'we'
            P0, Sv0 = 3.e5, 0.2
            bdy_incons = dat.grid.incons((P0, Sv0))
            j = dat.boundaries_json(geo, bdy_incons, atmos_volume, eos, mesh_coords)
            self.assertEqual(len(j['boundaries']), 1)
            self.assertEqual(j['boundaries'][0]['primary'], [P0, Sv0])
            self.assertEqual(j['boundaries'][0]['region'], 4)
            self.assertFalse('tracer' in j['boundaries'][0])
            self.assertEqual(j['boundaries'][0]['faces']['normal'], [0, 0, 1])
            self.assertEqual(j['boundaries'][0]['faces']['cells'], [0, 1, 2, 3])
            json.dumps(j)

            # Air, liquid top BCs
            eos = 'wae'
            P0, T0, Pa0 = 1.e5, 15., 0.1e5
            bdy_incons = dat.grid.incons((P0, T0, Pa0))
            j = dat.boundaries_json(geo, bdy_incons, atmos_volume, eos, mesh_coords)
            self.assertEqual(len(j['boundaries']), 1)
            self.assertEqual(j['boundaries'][0]['primary'], [P0, T0, Pa0])
            self.assertEqual(j['boundaries'][0]['region'], 1)
            self.assertFalse('tracer' in j['boundaries'][0])
            self.assertEqual(j['boundaries'][0]['faces']['normal'], [0, 0, 1])
            self.assertEqual(j['boundaries'][0]['faces']['cells'], [0, 1, 2, 3])
            json.dumps(j)

            # added hotplate BCs on bottom
            for col in geo.columnlist:
                blkname = geo.block_name('99', col.name)
                blk = t2block(blkname, atmos_volume, dat.grid.rocktypelist[0])
                dat.grid.add_block(blk)
                interior_blk_name = geo.block_name(geo.layerlist[-1].name, col.name)
                interior_blk = dat.grid.block[interior_blk_name]
                conblocks = [interior_blk, blk]
                dists = [0.5 * geo.layerlist[1].thickness, 0.]
                con = t2connection(conblocks, dists, col.area, dircos = 1.)
                dat.grid.add_connection(con)

            eos = 'we'
            P0, T0 = 1.e5, 15.
            bdy_incons = dat.grid.incons((P0, T0))
            Pb, Tb = 4.e5, 25.
            for blk in dat.grid.blocklist[-geo.num_columns:]:
                bdy_incons[blk.name] = (Pb, Tb)
            j = dat.boundaries_json(geo, bdy_incons, atmos_volume, eos, mesh_coords)
            self.assertEqual(len(j['boundaries']), 8)
            for ibc, bc in enumerate(j['boundaries'][:4]):
                self.assertEqual(bc,
                                 {'primary': [P0, T0], 'region': 1,
                                  'faces': {'cells': [ibc], 'normal': [0., 0., 1.]}})
            for ibc, bc in enumerate(j['boundaries'][4:]):
                cellindex = geo.num_underground_blocks - geo.num_columns + ibc
                self.assertEqual(bc,
                                 {'primary': [Pb, Tb], 'region': 1,
                                  'faces': {'cells': [cellindex], 'normal': [0., 0., -1.]}})
            json.dumps(j)

            # add side BCs
            cell_indices = [2, 3, 6, 7]
            interior_block_indices = [cell_index + geo.num_atmosphere_blocks
                                      for cell_index in cell_indices]
            interior_blknames = [dat.grid.blocklist[interior_block_index].name for
                                 interior_block_index in interior_block_indices]
            new_blknames = [' bc%2d' %i for i in cell_indices]
            for interior_blkname, new_blkname in zip(interior_blknames, new_blknames):
                centre = dat.grid.block[interior_blkname].centre + np.array([0.5 * dx, 0., 0.])
                blk = t2block(new_blkname, atmos_volume, dat.grid.rocktypelist[0], centre)
                dat.grid.add_block(blk)
                conblocks = [dat.grid.block[blkname] for blkname in [interior_blkname, new_blkname]]
                dists = [0.5 * dx, 0.]
                area = dy * dz
                con = t2connection(conblocks, dists, area, dircos = 0.)
                dat.grid.add_connection(con)
            bdy_incons = dat.grid.incons((P0, T0))
            Ps, Ts = 3.e5, 20.
            for blk in dat.grid.blocklist[-len(cell_indices):]:
                bdy_incons[blk.name] = (Ps, Ts)
            j = dat.boundaries_json(geo, bdy_incons, atmos_volume, eos, mesh_coords)
            self.assertEqual(len(j['boundaries']), 12)
            for ibc, bc in enumerate(j['boundaries'][-4:]):
                cellindex = cell_indices[ibc]
                self.assertEqual(bc,
                                 {'primary': [Ps, Ts], 'region': 1,
                                  'faces': {'cells': [cellindex], 'normal': [1., 0., 0.]}})
            json.dumps(j)

            # grid with no block centres:
            dat.grid = t2grid().fromgeo(geo)
            for blk in dat.grid.blocklist: blk.centre = None
            eos = 'we'
            P0, T0 = 1.e5, 15.
            bdy_incons = dat.grid.incons((P0, T0))
            j = dat.boundaries_json(geo, bdy_incons, atmos_volume, eos, mesh_coords)
            self.assertEqual(len(j['boundaries']), 1)
            self.assertEqual(j['boundaries'][0]['primary'], [P0, T0])
            self.assertEqual(j['boundaries'][0]['region'], 1)
            self.assertEqual(j['boundaries'][0]['faces']['normal'], [0, 0, 1])
            self.assertEqual(j['boundaries'][0]['faces']['cells'], [0, 1, 2, 3])
            json.dumps(j)

            # added hotplate BCs on bottom
            for col in geo.columnlist:
                blkname = geo.block_name('99', col.name)
                blk = t2block(blkname, atmos_volume, dat.grid.rocktypelist[0])
                dat.grid.add_block(blk)
                interior_blk_name = geo.block_name(geo.layerlist[-1].name, col.name)
                interior_blk = dat.grid.block[interior_blk_name]
                conblocks = [interior_blk, blk]
                dists = [0.5 * geo.layerlist[1].thickness, 0.]
                con = t2connection(conblocks, dists, col.area, dircos = 1.)
                dat.grid.add_connection(con)
            bdy_incons = dat.grid.incons((P0, T0))
            Pb, Tb = 4.e5, 25.
            for blk in dat.grid.blocklist[-geo.num_columns:]:
                bdy_incons[blk.name] = (Pb, Tb)
            j = dat.boundaries_json(geo, bdy_incons, atmos_volume, eos, mesh_coords)
            self.assertEqual(len(j['boundaries']), 8)
            for ibc, bc in enumerate(j['boundaries'][:4]):
                self.assertEqual(bc,
                                 {'primary': [P0, T0], 'region': 1,
                                  'faces': {'cells': [ibc], 'normal': [0., 0., 1.]}})
            for ibc, bc in enumerate(j['boundaries'][4:]):
                cellindex = geo.num_underground_blocks - geo.num_columns + ibc
                self.assertEqual(bc,
                                 {'primary': [Pb, Tb], 'region': 1,
                                  'faces': {'cells': [cellindex], 'normal': [0., 0., -1.]}})
            json.dumps(j)

            # 2D xy mesh with BCs at x = 200:
            nx, ny, nz = 2, 2, 1
            dx, dy, dz = 100., 100., 1.
            mesh_coords = 'xy'
            geo = mulgrid().rectangular([dx]*nx, [dy]*ny, [dz]*nz)
            dat = t2data()
            dat.grid = t2grid().fromgeo(geo)
            cell_indices = [2, 3]
            interior_block_indices = [cell_index + geo.num_atmosphere_blocks
                                      for cell_index in cell_indices]
            interior_blknames = [dat.grid.blocklist[interior_block_index].name for
                                 interior_block_index in interior_block_indices]
            new_blknames = [' bc%2d' %i for i in cell_indices]
            for interior_blkname, new_blkname in zip(interior_blknames, new_blknames):
                centre = dat.grid.block[interior_blkname].centre + np.array([0.5 * dx, 0., 0.])
                blk = t2block(new_blkname, atmos_volume, dat.grid.rocktypelist[0], centre)
                dat.grid.add_block(blk)
                conblocks = [dat.grid.block[blkname] for blkname in [interior_blkname, new_blkname]]
                dists = [0.5 * dx, 0.]
                area = dy * dz
                con = t2connection(conblocks, dists, area, dircos = 0.)
                dat.grid.add_connection(con)
            bdy_incons = dat.grid.incons((P0, T0))
            Ps, Ts = 3.e5, 20.
            for blk in dat.grid.blocklist[-len(cell_indices):]:
                bdy_incons[blk.name] = (Ps, Ts)
            j = dat.boundaries_json(geo, bdy_incons, atmos_volume, eos, mesh_coords)
            self.assertEqual(len(j['boundaries']), 1)
            self.assertEqual(j['boundaries'][0]['primary'], [Ps, Ts])
            self.assertEqual(j['boundaries'][0]['region'], 1)
            self.assertFalse('tracer' in j['boundaries'][0])
            self.assertEqual(j['boundaries'][0]['faces']['normal'], [1, 0])
            self.assertEqual(j['boundaries'][0]['faces']['cells'], cell_indices)
            json.dumps(j)

        basic_test()
        mesh_test()
        eos_test()
        output_test()
        timestepping_test()
        relative_permeability_test()
        capillary_pressure_test()
        primary_to_region_test()
        initial_test()
        generators_test()
        network_test()
        boundaries_test()

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(t2dataTestCase)
    unittest.TextTestRunner(verbosity = 1).run(suite)

