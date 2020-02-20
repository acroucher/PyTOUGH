import unittest
import os
from t2data import *

class t2gridTestCase(unittest.TestCase):

#------------------------------------------------------------------------

    def test_block_centres(self):
        """test calculate_block_centres() method, with a block mapping"""
        geo = mulgrid().rectangular([1., 2., 3.], [1.], [1.])
        blknames = ['    1', '   2', '    3']
        blockmap = dict(zip(geo.block_name_list, blknames))
        grd = t2grid().fromgeo(geo, blockmap)
        grd.calculate_block_centres(geo, blockmap)
        xcentres, y, z = [0.5, 2., 4.5], 0.5, -0.5
        for name, x in zip(blknames, xcentres):
            self.assertEqual(grd.block[name].centre[0], x)
            self.assertEqual(grd.block[name].centre[1], y)
            self.assertEqual(grd.block[name].centre[2], z)

    def minc_compare(self, grid1, grid2, num_frac):
        """Compare two MINC grids for a specified number of MINC fracture planes."""
        vol_tol = 1.e-2
        dist_tol = 1.e-2
        area_tol = 1.e-2
        nstr = str(num_frac) + " fracture planes"
        self.assertEqual(grid1.num_blocks, grid2.num_blocks,
                         "MINC num blocks, " + nstr)
        self.assertTrue(all([blk.rocktype.name == grid2.block[blk.name].rocktype.name
                             for blk in grid1.blocklist]))
        v1 = np.array([blk.volume for blk in grid1.blocklist])
        v2 = np.array([grid2.block[blk.name].volume for blk in grid1.blocklist])
        self.assertTrue(np.allclose(v1, v2, rtol = vol_tol),
                        "MINC volumes, level " + nstr)
        self.assertEqual(grid1.num_connections, grid2.num_connections,
                         "MINC num connections, " + nstr)
        for i in range(2):
            name1 = [con.block[i].name for con in grid1.connectionlist]
            name2 = [grid2.connection[tuple([blk.name for blk in con.block])].block[i].name
                           for con in grid1.connectionlist]
            self.assertTrue(all([n1 == n2 for n1, n2 in zip(name1, name2)]))
        for i in range(2):
            d1 = np.array([con.distance[i] for con in grid1.connectionlist])
            d2 = np.array([grid2.connection[tuple([blk.name for blk in con.block])].distance[i]
                           for con in grid1.connectionlist])
            self.assertTrue(np.allclose(d1, d2, rtol = dist_tol),
                            "MINC connections distance " + str(i) + " ," + nstr)
        a1 = np.array([con.area for con in grid1.connectionlist])
        a2 = np.array([grid2.connection[tuple([blk.name for blk in con.block])].area
                           for con in grid1.connectionlist])
        self.assertTrue(np.allclose(a1, a2, rtol = area_tol),
                            "MINC areas, level " + nstr)

    def test_minc(self):
        """test minc() method, 1-3 fracture planes"""
        path = os.path.join('grid', 'minc')
        def minc_name(name, i):
            first = ['a', 'w', 's', 'o', 'k', 'g']
            return first[i] + name[1:]
        dat = t2data(os.path.join(path, 'orig.dat'))
        nf = 1
        vol = [0.1, 0.4, 0.5]
        spacing = 25.
        dat.grid.minc(vol, spacing, nf, matrix_blockname = minc_name)
        dat2 = t2data(os.path.join(path, 'minc1.dat'))
        self.minc_compare(dat.grid, dat2.grid, nf)

        dat = t2data(os.path.join(path, 'orig.dat'))
        nf = 2
        vol = [0.05, 0.1, 0.3, 0.55]
        spacing = 35.
        dat.grid.minc(vol, spacing, nf, matrix_blockname = minc_name)
        dat2 = t2data(os.path.join(path, 'minc2.dat'))
        self.minc_compare(dat.grid, dat2.grid, nf)

        dat = t2data(os.path.join(path, 'orig.dat'))
        nf = 3
        vol = [0.05, 0.1, 0.2, 0.25, 0.4]
        spacing = 50.
        dat.grid.minc(vol, spacing, nf, matrix_blockname = minc_name)
        dat2 = t2data(os.path.join(path, 'minc3.dat'))
        self.minc_compare(dat.grid, dat2.grid, nf)

    def test_partial_minc(self):
        """test minc() over part of grid"""
        d = 100.
        nblks = 10
        geo = mulgrid().rectangular([d]*nblks, [d], [d])
        grid = t2grid().fromgeo(geo)

        P0, T0 = 2.e5, 20.
        inc = grid.incons([P0, T0])
        for i, blk in enumerate(inc):
            inc[blk.block] = [P0, T0 + i]
        nf = 2
        vol = [10, 20, 30, 40]
        spacing = 0.5 * d
        blk_indices = [3, 4, 7, 8]
        blks = [geo.block_name_list[i] for i in blk_indices]
        nmat = len(vol) - 1
        nfblks = len(blks)
        iblk, mincon = grid.minc(vol, spacing, nf, blocks = blks, incon = inc)

        n = geo.num_blocks + nfblks * nmat
        self.assertEqual(grid.num_blocks, n, "Partial MINC num blocks")

        expected_iblk = np.array([[3, 4, 7, 8], [10, 13, 16, 19],
                                 [11, 14, 17, 20], [12, 15, 18, 21]])
        self.assertTrue(np.all(iblk == expected_iblk), "Partial MINC block index array")

        expected_t = list(np.arange(T0, T0 + nblks, 1.))
        for i in blk_indices:
            expected_t += [T0 + i] * nmat
        expected_t = np.array(expected_t)
        tol = 1.e-6
        self.assertTrue(np.allclose(mincon.variable[:,1],expected_t, rtol = tol),
                        "Partial MINC incons")

    def test_rectgeo(self):
        """test rectgeo()"""

        def model(atmos_type = 2, flat = True):
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
            dat = t2data()
            grid = t2grid().fromgeo(geo)
            return geo, grid

        def grid_test(geo, grid, geo1, grid1):

            stol = 1.e-2
            s = np.array([col.surface for col in geo.columnlist])
            s1 = np.array([col.surface for col in geo1.columnlist])
            self.assertEqual(len(s), len(s1))
            self.assertTrue(np.allclose(s, s1, rtol = stol))

            vtol = 1.e-3
            v = np.array([blk.volume for blk in grid.blocklist])
            v1 = np.array([blk.volume for blk in grid1.blocklist])
            self.assertEqual(len(v), len(v1))
            self.assertTrue(np.allclose(v, v1, rtol = vtol))

        snap = 1.e-3

        atmos_type, flat = 2, True
        geo, grid = model(atmos_type = atmos_type, flat = flat)
        geo1, blockmap = grid.rectgeo(remove_inactive = False, atmos_type = atmos_type,
                                      layer_snap = snap)
        grid1 = t2grid().fromgeo(geo1)
        grid_test(geo, grid, geo1, grid1)

        atmos_type, flat = 1, True
        geo, grid = model(atmos_type = atmos_type, flat = flat)
        geo1, blockmap = grid.rectgeo(remove_inactive = False, atmos_type = atmos_type,
                                      layer_snap = snap)
        grid1 = t2grid().fromgeo(geo1)
        grid_test(geo, grid, geo1, grid1)

        atmos_type, flat = 1, False
        geo, grid = model(atmos_type = atmos_type, flat = flat)
        geo1, blockmap = grid.rectgeo(remove_inactive = False, atmos_type = atmos_type,
                                      layer_snap = snap)
        grid1 = t2grid().fromgeo(geo1)
        grid_test(geo, grid, geo1, grid1)

        atmos_type, flat = 0, False
        geo, grid = model(atmos_type = atmos_type, flat = flat)
        geo1, blockmap = grid.rectgeo(remove_inactive = False, atmos_type = atmos_type,
                                      layer_snap = snap)
        grid1 = t2grid().fromgeo(geo1)
        grid_test(geo, grid, geo1, grid1)

    def test_reorder(self):
        """test reorder()"""
        d = 100.
        nblks = [4, 4, 3]
        geo = mulgrid().rectangular([d] * nblks[0], [d] * nblks[1], [d] * nblks[2])
        grid = t2grid().fromgeo(geo)
        indices = list(range(geo.num_blocks))
        indices.reverse()
        block_names = [geo.block_name_list[i] for i in indices]
        indices = list(range(geo.num_block_connections))
        indices.reverse()
        connection_names = [geo.block_connection_name_list[i] for i in indices]
        swap_indices = [5, 9, 21, 29]
        for i in swap_indices: connection_names[i] = connection_names[i][::-1]
        grid.reorder(block_names, connection_names)
        grid_block_names = [blk.name for blk in grid.blocklist]
        self.assertEqual(block_names, grid_block_names)
        grid_dict_block_names = [grid.block[name].name for name in geo.block_name_list]
        self.assertEqual(grid_dict_block_names, geo.block_name_list)
        grid_connection_names = [tuple([blk.name for blk in con.block])
                                 for con in grid.connectionlist]
        self.assertEqual(connection_names, grid_connection_names)
        self.assertTrue(grid.check(fix = False, silent = True))
        grid.reorder(geo = geo)
        grid_block_names = [blk.name for blk in grid.blocklist]
        self.assertEqual(geo.block_name_list, grid_block_names)
        grid_dict_block_names = [grid.block[name].name for name in geo.block_name_list]
        self.assertEqual(grid_dict_block_names, geo.block_name_list)
        grid_connection_names = [tuple([blk.name for blk in con.block])
                                 for con in grid.connectionlist]
        self.assertEqual(geo.block_connection_name_list, grid_connection_names)

    def test_create_rocktypes(self):
        """test rocktype creation"""
        grid = t2grid()
        a = rocktype()
        a.name = 'aaa01'
        grid.add_rocktype(a)
        b = rocktype()
        b.name = 'aaa02'
        grid.add_rocktype(b)
        newperm = 1.e-10
        grid.rocktype['aaa01'].permeability[0] = newperm
        self.assertEqual(newperm, grid.rocktypelist[0].permeability[0])
        self.assertEqual(1.e-15, grid.rocktypelist[1].permeability[0])

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(t2gridTestCase)
    unittest.TextTestRunner(verbosity = 1).run(suite)

