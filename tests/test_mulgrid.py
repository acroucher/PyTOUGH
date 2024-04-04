import unittest
import os
from mulgrids import *

class mulgrid_stats(mulgrid):
    """Variant of mulgrid class with extra properties for vital statistics
    of a mulgrid object- for comparisons between mulgrids that should
    be the same.
    """

    def get_node_names(self): return [n.name for n in self.nodelist]
    node_names = property(get_node_names)
    def get_column_names(self): return [col.name for col in self.columnlist]
    column_names = property(get_column_names)
    def get_node_positions(self): return np.array([n.pos for n in self.nodelist])
    node_positions = property(get_node_positions)
    def get_column_nodes(self): return [[n.name for n in col.node] for col in self.columnlist]
    column_nodes = property(get_column_nodes)
    def get_layer_names(self): return [lay.name for lay in self.layerlist]
    layer_names = property(get_layer_names)
    def get_layer_centres(self): return np.array([lay.centre for lay in self.layerlist])
    layer_centres = property(get_layer_centres)
    def get_well_names(self): return [w.name for w in self.welllist]
    well_names = property(get_well_names)
    def get_well_positions(self): return np.array(sum([w.pos for w in self.welllist],[]))
    well_positions = property(get_well_positions)

    def compare(self, other, testcase):
        """Compares self with another mulgrid_stats object and reports results
        to a unittest testcase."""
        wellpostol = 0.06  # need larger tolerance for well positions, due to roundoff
        msg = ' mismatch for geometry ' + self.filename
        testcase.assertEqual(self.type, other.type, 'type' + msg)
        testcase.assertEqual(self.convention, other.convention, 'convention' + msg)
        testcase.assertEqual(self.atmosphere_type, other.atmosphere_type,
                             'atmosphere type' + msg)
        testcase.assertEqual(self.atmosphere_volume, other.atmosphere_volume,
                             'atmostphere volume' + msg)
        testcase.assertEqual(self.atmosphere_connection, other.atmosphere_connection,
                             'atmosphere connection' + msg)
        testcase.assertEqual(self.unit_type, other.unit_type, 'unit type' + msg)
        testcase.assertEqual(self.gdcx, other.gdcx, 'gdcx' + msg)
        testcase.assertEqual(self.gdcy, other.gdcy, 'gdcy' + msg)
        testcase.assertEqual(self.cntype, other.cntype, 'cntype' + msg)
        testcase.assertEqual(self.permeability_angle, other.permeability_angle,
                             'permeability angle' + msg)

        testcase.assertEqual(self.node_names, other.node_names, 'node names' + msg)
        testcase.assertEqual(self.column_names, other.column_names, 'column names' + msg)
        testcase.assertTrue(np.allclose(self.node_positions, other.node_positions),
                            'node positions' + msg)
        testcase.assertEqual(self.column_nodes, other.column_nodes, 'column nodes' + msg)
        testcase.assertEqual(self.layer_names, other.layer_names, 'layer names' + msg)
        testcase.assertTrue(np.allclose(self.layer_centres, other.layer_centres),
                            'layer centres' + msg)
        testcase.assertEqual(self.well_names, other.well_names, 'well names' + msg)
        testcase.assertTrue(np.max(np.abs(self.well_positions - other.well_positions)) <= wellpostol,
                            'well positions' + msg)

class mulgridTestCase(unittest.TestCase):

    def test_fix_unfix_blockname(self):
        """fix_blockname() and unfix_blockname() functions"""
        blk = 'AA1 6'
        fblk = fix_blockname(blk)
        self.assertNotEqual(fblk[-2], ' ')
        blk2 = unfix_blockname(fblk)
        self.assertEqual(blk, blk2)
        blk = 'A 200'
        self.assertEqual(blk, fix_blockname(blk))
        blk = 'abcde'
        self.assertEqual(blk, fix_blockname(blk))
        self.assertEqual(blk, unfix_blockname(blk))

    def test_valid_blockname(self):
        """valid_blockname() function"""
        self.assertTrue(valid_blockname('AA123'))
        self.assertFalse(valid_blockname('AA12a'))
        self.assertTrue(valid_blockname('A   1'))
        self.assertTrue(valid_blockname('A? 21'))

    def test_fortran_float(self):
        """fortran_float() function"""
        self.assertEqual(fortran_float('   '), 0.0)
        self.assertEqual(fortran_float('1.0'), 1.0)
        self.assertEqual(fortran_float('  -10.3'), -10.3)
        self.assertEqual(fortran_float(' 0.3E+06'), 0.3e6)
        self.assertEqual(fortran_float(' 0.3E+006'), 0.3e6)
        self.assertEqual(fortran_float('1.234-314'), 1.234e-314)
        self.assertEqual(fortran_float('1.234+314'), 1.234e314)
        self.assertEqual(fortran_float('-0.14326-124'), -0.14326e-124)
        self.assertEqual(fortran_float(' 0.3D+06'), 0.3e6)
        self.assertEqual(fortran_float(' 0.3D-06'), 0.3e-6)
        self.assertTrue(np.isnan(fortran_float('1.234x-10')))

    def test_rect_grid_operations(self):
        """rectangular grid operations"""
        dx = np.logspace(1, 1.5, 10)
        dx = np.hstack((dx[::-1], dx))
        dy = np.logspace(0.8, 1.4, 9)
        dy = np.hstack((dy[::-1],dy))
        dz = np.logspace(0.6, 3, 12)
        geo = mulgrid().rectangular(dx, dy, dz)
        # check numbers of columns, layers and blocks:
        num_cols = np.size(dx)*np.size(dy)
        num_layers = np.size(dz)
        num_blocks = num_cols * num_layers
        self.assertEqual(geo.num_columns, num_cols)
        self.assertEqual(geo.num_layers, num_layers + 1)
        self.assertEqual(geo.num_blocks, num_blocks)
        # check area and centre:
        lx = np.sum(dx)
        ly = np.sum(dy)
        a = lx * ly
        c = 0.5*np.array([lx, ly])
        self.assertAlmostEqual(geo.area, a)
        self.assertTrue(np.allclose(geo.centre, c))
        # translate, rotate:
        d = np.array([5321.5, -1245.2, 1000.])
        theta = 37.6
        geo.translate(d)
        c2 = c + d[:2]
        self.assertAlmostEqual(geo.area, a)
        self.assertTrue(np.allclose(geo.centre, c2))
        geo.rotate(theta)
        self.assertAlmostEqual(geo.area, a)
        self.assertTrue(np.allclose(geo.centre, c2))
        geo.rotate(-theta)
        geo.translate(-d)
        self.assertAlmostEqual(geo.area, a)
        self.assertTrue(np.allclose(geo.centre, c))

    def test_read_g2(self):
        """reading g2 geometry"""

        try:
            filename = os.path.join('mulgrid', 'g2.dat')
            geo = mulgrid(filename)
        except: self.fail('Could not open mulgrid: ' + filename)

        self.assertTrue(geo.num_blocks == 28966)
        self.assertTrue(geo.num_columns == 1036)
        self.assertTrue(geo.num_layers == 35)
        self.assertTrue(geo.num_wells == 186)

        self.assertEqual(
                [n.name for n in geo.nodelist[:10]],
                ['mfq', 'aky', 'acy', 'bny', 'mbq', 'aku', 'acu', 'bru', 'mbi', 'akm'])

        self.assertEqual(
                [col.name for col in geo.columnlist[-10:]],
                ['gko', 'goo', 'fbo', 'ffo', 'fjo', 'fno', 'fro', 'fvo', 'eco', 'ego'])

        self.assertTrue(
            (np.array([n.pos for n in geo.nodelist[100:110]]) ==
            np.array([[ 2775720.36,  6283160.18], [ 2775686.41,  6281497.33],
                      [ 2776095.96,  6281052.2 ], [ 2776542.34,  6280565.66],
                      [ 2777019.52,  6280045.87], [ 2777554.27,  6279456.3 ],
                      [ 2778114.41,  6278828.88], [ 2778669.16,  6278196.73],
                      [ 2779219.24,  6277560.52], [ 2779852.28,  6276811.46]])).all())

        self.assertEqual(
            [n.name for n in geo.column['nio'].node],
            ['nkq', 'nko', 'nio', 'niq'])

        self.assertTrue(
            (geo.column['eic'].centre == np.array([2778311.88, 6285540.38])).all())

        self.assertEqual(
            [l.centre for l in geo.layerlist],
             [880.0, 840.0, 775.0, 725.0, 675.0, 625.0, 575.0, 525.0,
              475.0, 425.0, 375.0, 325.0, 275.0, 225.0,
              175.0, 125.0, 75.0, 25.0, -25.0, -75.0, -125.0, -175.0,
              -225.0, -275.0, -350.0, -450.0, -562.5,
              -687.5, -875.0, -1125.0, -1375.0, -1750.0, -2250.0, -2750.0, -3250.0])

        self.assertEqual(
            [geo.column[col].surface for col in \
                 ['fbo', 'ffo', 'fjo', 'fno', 'fro', 'fvo', 'eco', 'ego', 'eko', 'eoo', 'eso']],
            [355.18, 334.37, 333.61, 327.96, 320.19, 385.88, 380.84, 427.08, 617.12, 571.32, 637.08])

        self.assertEqual(geo.welllist[4].name, 'A   5')

        self.assertTrue(
            (geo.well['I 304'].pos_coordinate(2)[:6] == 
             np.array([341.0, 240.0, 189.0, 171.0, 143.0, 121.0])).all())

        self.assertTrue(
            (geo.well['I 310'].pos_coordinate(0)[:5] == 
             np.array([2779465.0, 2779465.0, 2779465.0, 2779467.0, 2779470.0])).all())

    def test_read_write(self):
        """read/write tests"""
        geoms = ['g1.dat', 'g2.dat', 'g3.dat', 'g4.dat']
        from os import remove
        outfilename = 'gread_write_test.dat'
        for geoname in geoms:
            filename = os.path.join('mulgrid', geoname)
            geo = mulgrid_stats(filename)
            geo.write(outfilename)
            geo.filename = filename
            geo2 = mulgrid_stats(outfilename)
            geo.compare(geo2, self)
            remove(outfilename)

    def test_fit_surface(self):
        """fit surface"""
        tol = 1.e-3
        geo = mulgrid(os.path.join('mulgrid', 'g2.dat'))
        p = np.array([2.78153e6, 6.26941e6])
        p2 = p + np.ones(2)*5.e3
        r = [p, p2]
        cols = geo.columns_in_polygon(r)
        d = np.load(os.path.join('mulgrid', 'fit_surface_data.npy'))
        geo.fit_surface(d, alpha = 0.1, beta = 0.1, silent = True)
        s = np.array([col.surface for col in cols])
        r = np.load(os.path.join('mulgrid', 'fit_surface_result.npy'))
        self.assertTrue(np.allclose(s, r, atol = tol))

    def test_refine(self):
        """refine()"""
        tol = 1.e-3
        geo = mulgrid().rectangular([100.]*10, [150.]*8, [10.]*3)
        orig_area = geo.area
        geo.rotate(30.)
        cols = [col for col in geo.columnlist if 400. < col.centre[1] < 600.]
        geo.refine(cols)
        self.assertEqual(geo.area, orig_area)
        self.assertEqual(geo.num_columns, 182)
        self.assertEqual(geo.num_nodes, 169)
        areas = np.sort(np.array([col.area for col in geo.columnlist]))
        a = np.sort(np.load(os.path.join('mulgrid', 'refine_areas.npy')))
        self.assertTrue(np.allclose(areas, a, atol = tol))

    def test_rename_column(self):
        """rename_column()"""
        geo = mulgrid().rectangular([100.]*10, [150.]*8, [10.]*3)
        n = 20
        origcols = [col.name for col in geo.columnlist] 
        oldcols = origcols[:n]
        origblocks  = geo.block_name_list
        newcols = [name.upper() for name in oldcols]
        geo.rename_column(oldcols,newcols)
        expectedcols = ['  ' + l for l in ascii_uppercase[:n]] + origcols[n:]
        self.assertEqual([col.name for col in geo.columnlist], expectedcols)
        expectedblocks = [(blk[:2] + blk[2].upper() + blk[3:]
                           if blk[:2] == '  ' and blk[2] < ascii_lowercase[n]
                           else blk) for blk in origblocks]
        self.assertEqual(geo.block_name_list, expectedblocks)

    def test_column_track(self):
        """column_track()"""
        geo = mulgrid(os.path.join('mulgrid', 'g5.dat'))
        names = ['aag', 'aah', 'abh', 'ach', 'aci', 'adi', 'aei', 'afi', 'afj', 'agj', 'ahj',
                 'aik', 'ajk', 'akk', 'akl', 'all', 'aml', 'anl', 'anm', 'aom', 'apm', 'apn']
        try:
            t = geo.column_track(geo.bounds)
            found_names = [ts[0].name for ts in t]
            err = False
            self.assertEqual(found_names, names)
        except: err = True
        self.assertFalse(err)

        # track through grid with a corner removed:
        geo = mulgrid().rectangular([100]*2, [100]*2, [10])
        geo.translate([1e7, 1e7, 0])
        geo.delete_column(geo.columnlist[0].name)
        geo.rotate(30)
        line = [geo.columnlist[0].centre, geo.columnlist[1].centre]
        t = geo.column_track(line)
        self.assertEqual(len(t), 2)
        if len(t) > 1:
            self.assertEqual(t[0][0].name, geo.columnlist[0].name)
            self.assertEqual(t[1][0].name, geo.columnlist[1].name)

        # M-grid:
        geo = mulgrid().rectangular([100]*5, [100]*3, [10])
        names = [geo.columnlist[i].name for i in [1,3,6,8]]
        for name in names: geo.delete_column(name)
        line = [np.array([0, 50]), np.array([500, 50])]
        t = geo.column_track(line)
        self.assertEqual(len(t), 3)
        if (len(t) == 3):
            self.assertEqual(t[0][0].name, geo.columnlist[0].name)
            self.assertEqual(t[1][0].name, geo.columnlist[1].name)
            self.assertEqual(t[2][0].name, geo.columnlist[2].name)

        line = line[::-1]
        t = geo.column_track(line)
        self.assertEqual(len(t), 3)
        if (len(t) == 3):
            self.assertEqual(t[0][0].name, geo.columnlist[2].name)
            self.assertEqual(t[1][0].name, geo.columnlist[1].name)
            self.assertEqual(t[2][0].name, geo.columnlist[0].name)

        # track within a single column:
        col = geo.columnlist[-1]
        p = col.centre
        d = np.ones(2) * 25
        line = [p - d, p + d]
        t = geo.column_track(line)
        self.assertEqual(len(t), 1)
        if (len(t) == 1):
            self.assertEqual(t[0][0].name, col.name)

        # track outside grid:
        line = [np.array([-10, -10]), np.array([600, -200])]
        t = geo.column_track(line)
        self.assertEqual(len(t), 0)

        # track entirely within grid:
        line = [np.array([380, 270]), np.array([490, 90])]
        t = geo.column_track(line)
        self.assertEqual(len(t), 4)
        self.assertTrue(np.allclose(t[0][1], line[0]))
        self.assertTrue(np.allclose(t[-1][2], line[1]))
        if (len(t) == 4):
            names = [ti[0].name for ti in t]
            col_ind = [9, 10, 5, 2]
            expected_names = [geo.columnlist[i].name for i in col_ind]
            self.assertEqual(names, expected_names)

        # track leaving grid:
        line = [np.array([250, 50]), np.array([300, -100])]
        t = geo.column_track(line)
        self.assertEqual(len(t), 1)
        if (len(t) == 1):
            self.assertEqual(t[0][0].name, geo.columnlist[1].name)
            self.assertTrue(np.allclose(t[0][1], line[0]))

        # track entering grid:
        line = [np.array([270, -50]), np.array([230, 270])]
        t = geo.column_track(line)
        self.assertEqual(len(t), 3)
        if (len(t) == 3):
            names = [ti[0].name for ti in t]
            col_ind = [1, 4, 8]
            expected_names = [geo.columnlist[i].name for i in col_ind]
            self.assertEqual(names, expected_names)
            self.assertTrue(np.allclose(t[-1][2], line[1]))

        # 5x5 grid:
        geo = mulgrid().rectangular([100]*5, [100]*5, [10])
        line = [np.array([100, 100]), np.array([400, 400])]
        t = geo.column_track(line)
        self.assertEqual(len(t), 3)
        if (len(t) == 3):
            self.assertTrue(np.allclose(t[0][1], line[0]))
            self.assertTrue(np.allclose(t[-1][2], line[1]))

    def test_grid3d(self):
        """3D grid"""

        def get_bounding_box(nodes):
            large = 1.e99
            box = []
            for i in range(3): box.append([large, -large])
            for n in nodes:
                for i in range(3):
                    box[i][0] = min(box[i][0], n[i])
                    box[i][1] = max(box[i][1], n[i])
            return np.array(box)

        for atmos_type in range(3):

            geo = mulgrid().rectangular([1000.] * 6, [1000.]*6, [10., 20., 30.],
                                        atmos_type = atmos_type,
                                        block_order = 'dmplex')
            snap = 0.1
            nodes, elts = geo.grid3d(surface_snap = snap)
            self.assertEqual(108, len(elts))
            self.assertEqual(196, len(nodes))
            bounds = get_bounding_box(nodes)
            expected_bounds = np.array([[0, 6000], [0, 6000], [-60, 0]])
            self.assertTrue(np.allclose(bounds, expected_bounds))

            cols = [col.name for col in geo.columnlist if col.centre[0] > 3000]
            geo.refine(cols)
            nodes, elts = geo.grid3d(surface_snap = snap)
            self.assertEqual(306, len(elts))
            self.assertEqual(448, len(nodes))
            self.assertEqual([8]*252 + [6]*54, [len(elt) for elt in elts])
            bounds = get_bounding_box(nodes)
            self.assertTrue(np.allclose(bounds, expected_bounds))

            cols = [col for col in geo.columnlist if col.centre[0] < 1000]
            for col in cols:
                col.surface = 5.
                geo.set_column_num_layers(col)
            geo.setup_block_name_index()
            geo.setup_block_connection_name_index()
            nodes, elts = geo.grid3d(surface_snap = snap)
            self.assertEqual(306, len(elts))
            self.assertEqual(455, len(nodes))
            self.assertEqual([8]*252 + [6]*54, [len(elt) for elt in elts])
            expected_bounds[2] = [-60, 5]
            bounds = get_bounding_box(nodes)
            self.assertTrue(np.allclose(bounds, expected_bounds))

            for col in cols:
                col.surface = -12.
                geo.set_column_num_layers(col)
            geo.setup_block_name_index()
            geo.setup_block_connection_name_index()
            nodes, elts = geo.grid3d(surface_snap = snap)
            self.assertEqual(300, len(elts))
            self.assertEqual(448, len(nodes))
            self.assertEqual([8]*246 + [6]*54, [len(elt) for elt in elts])
            bounds = get_bounding_box(nodes)
            expected_bounds[2] = [-60, 0]
            self.assertTrue(np.allclose(bounds, expected_bounds))

    def test_meshio(self):
        """meshio grid"""
        filename = os.path.join('mulgrid', 'g5')
        geofilename = filename + '.dat'
        geo = mulgrid(geofilename)
        points, cells = geo.meshio_grid()
        self.assertEqual(len(points), 9808)
        self.assertEqual(len(cells['hexahedron']), 7647)

    def test_block_name_containing_point(self):
        geo = mulgrid().rectangular([100.]*10, [150.]*8, [10.]*3)

        p = np.array([350., 500., -5.])
        blkname = geo.block_name_containing_point(p)
        self.assertEqual(' ah 1', blkname)

        # point above mesh:
        p[2] = 5.
        blkname = geo.block_name_containing_point(p)
        self.assertTrue(blkname is None)

        # point below mesh:
        p[2] = -50.
        blkname = geo.block_name_containing_point(p)
        self.assertTrue(blkname is None)

        # point outside horizontal mesh:
        p = np.array([-350., 500., -10.])
        blkname = geo.block_name_containing_point(p)
        self.assertTrue(blkname is None)

    def test_well_values(self):
        geo = mulgrid().rectangular([100.]*10, [150.]*8, [10.]*10)
        for col in geo.columnlist:
            col.surface = -9.
            geo.set_column_num_layers(col)
        geo.setup_block_name_index()
        geo.setup_block_connection_name_index()
        wellname = '  w 1'
        w = well(wellname, [np.array([450., 550., 0.]), np.array([400., 500., -100.])])
        geo.add_well(w)
        t = np.array([20. for name in geo.block_name_list])
        elev, val = geo.well_values(wellname, t, elevation = True)
        
    def test_amesh(self):
        infile = os.path.join('mulgrid', 'in')
        segfile = os.path.join('mulgrid', 'segmt')
        geo, blockmap = mulgrid().from_amesh(infile, segfile)
        self.assertEqual(913, geo.num_nodes)
        self.assertEqual(656, geo.num_columns)
        self.assertEqual(9184, geo.num_blocks)
        self.assertEqual(15, geo.num_layers)

    def test_gmsh(self):
        infile = os.path.join('mulgrid', 'gmsh2_2.msh')
        layers = [1, 2, 3]
        geo = mulgrid().from_gmsh(infile, layers)
        self.assertEqual(117, geo.num_nodes)
        self.assertEqual(96, geo.num_columns)
        self.assertEqual(96 * len(layers), geo.num_blocks)
        self.assertEqual(len(layers) + 1, geo.num_layers)

        infile = os.path.join('mulgrid', 'gmsh4_1.msh')
        layers = [1, 2, 3]
        geo = mulgrid().from_gmsh(infile, layers)
        self.assertEqual(48, geo.num_nodes)
        self.assertEqual(35, geo.num_columns)
        self.assertEqual(35 * len(layers), geo.num_blocks)
        self.assertEqual(len(layers) + 1, geo.num_layers)

    def test_layermesh(self):

        def layermesh_case(geo):
            lm = geo.layermesh
            self.assertEqual(geo.num_nodes, lm.num_nodes)
            self.assertEqual(geo.num_columns, lm.num_columns)
            self.assertEqual(geo.num_layers - 1, lm.num_layers)
            self.assertEqual(geo.num_underground_blocks, lm.num_cells)

            geo2 = geo.from_layermesh(lm, atmosphere_type = geo.atmosphere_type)
            self.assertEqual(geo.num_nodes, geo2.num_nodes)
            self.assertEqual(geo.num_columns, geo2.num_columns)
            self.assertEqual(geo.num_layers, geo2.num_layers)
            self.assertEqual(geo.num_blocks, geo2.num_blocks)

        geo = mulgrid().rectangular([100.]*10, [150.]*8, [10.]*10)
        for i, col in enumerate(geo.columnlist):
            col.surface = -40. + i / 80.
            geo.set_column_num_layers(col)
        geo.snap_columns_to_nearest_layers()
        layermesh_case(geo)

        geo = mulgrid().rectangular([100.]*9, [150.]*10, [10.]*11, atmos_type = 0)
        for i, col in enumerate(geo.columnlist):
            col.surface = -50. + i / 70.
            geo.set_column_num_layers(col)
        geo.snap_columns_to_nearest_layers()
        layermesh_case(geo)

        geo = mulgrid().rectangular([100.]*8, [120.]*7, [10.]*9, atmos_type = 1)
        for i, col in enumerate(geo.columnlist):
            col.surface = -30. + i / 60.
            geo.set_column_num_layers(col)
        geo.snap_columns_to_nearest_layers()
        layermesh_case(geo)

    def test_block_order(self):

        def get_block_nodes(geo):
            return [2 * geo.column[geo.column_name(blkname)].num_nodes
                    for blkname in geo.block_name_list]

        geo = mulgrid().rectangular([100.]*3, [100.], [10.]*2, block_order = 'dmplex')
        geo.refine([geo.columnlist[-1]])
        block_nodes = get_block_nodes(geo)
        self.assertEqual([8]*10 + [6]*6, block_nodes)

        # test re-setting block_order property
        geo = mulgrid().rectangular([100.]*3, [100.], [10.]*2)
        geo.refine([geo.columnlist[-1]])
        geo.block_order = 'dmplex'
        block_nodes = get_block_nodes(geo)
        self.assertEqual([8]*10 + [6]*6, block_nodes)

        # reading block order:
        filename = os.path.join('mulgrid', 'g7.dat')
        geo = mulgrid(filename)
        self.assertIsNone(geo.block_order)
        self.assertIsNone(geo._block_order_int)

        # overriding block order in file:
        for block_order in [None, 'layer_column', 'dmplex']:
            geo = mulgrid(filename, block_order = block_order)
            self.assertEqual(geo.block_order, block_order)

        # writing/reading block order:
        geo = mulgrid(filename)
        filename = os.path.join('mulgrid', 'g7_none.dat')
        geo.write(filename)
        geo = mulgrid(filename)
        self.assertIsNone(geo.block_order)
        self.assertIsNone(geo._block_order_int)
        os.remove(filename)

        geo.block_order = 'dmplex'
        filename = os.path.join('mulgrid', 'g7_dmplex.dat')
        geo.write(filename)
        geo = mulgrid(filename)
        self.assertEqual(geo.block_order, 'dmplex')
        self.assertEqual(geo._block_order_int, 1)

        # override block order in file:
        geo = mulgrid(filename, block_order = 'layer_column')
        self.assertEqual(geo.block_order, 'layer_column')
        self.assertEqual(geo._block_order_int, 0)
        os.remove(filename)

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(mulgridTestCase)
    unittest.TextTestRunner(verbosity = 1).run(suite)
