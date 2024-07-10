"""Class for TOUGH2 data.

Copyright 2011 University of Auckland.

This file is part of PyTOUGH.

PyTOUGH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PyTOUGH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PyTOUGH.  If not, see <http://www.gnu.org/licenses/>."""

from __future__ import print_function
import sys
from fixed_format_file import *
from t2grids import *
from t2incons import *
from math import ceil
import struct
from os.path import splitext

def primary_to_region_we(primary):
    """Returns thermodynamic region deduced from primary variables for EOS we."""
    from t2thermo import region
    if primary[1] < 1.: return 4
    else: return region(primary[1], primary[0])

def primary_to_region_wge(primary):
    """Returns thermodynamic region deduced from primary variables for wge
    (NCG) EOS (wce, wae)."""
    pwater = primary[0] - primary[2]
    return primary_to_region_we([pwater, primary[1]])

primary_to_region_funcs = {'w': primary_to_region_we, 'we': primary_to_region_we,
                           'wce': primary_to_region_wge, 'wae': primary_to_region_wge}
waiwera_eos_num_primary = {'w': 1, 'we': 2, 'wce': 3, 'wae': 3}

def trim_trailing_nones(vals):
    """Trim trailing None values from a list."""
    while vals and vals[-1] is None:
        vals.pop()
    return vals

t2data_format_specification = {
    'title': [['title'], ['80s']],
    'simulator': [['simulator'], ['80s']],
    'rocks1': [['name', 'nad', 'density', 'porosity',
               'k1', 'k2', 'k3', 'conductivity', 'specific_heat'], 
              ['5s', '5d'] + ['10.4e'] * 7],
    'rocks1.1': [['compressibility', 'expansivity', 'dry_conductivity',
                 'tortuosity', 'klinkenberg', 'xkd3', 'xkd4'],
                ['10.4e'] * 7],
    'rocks1.2': [['type', ''] + ['parameter'] * 7, ['5d', '5x'] + ['10.3e'] * 7],
    'rocks1.3': [['type', ''] + ['parameter'] * 7, ['5d', '5x'] + ['10.3e'] * 7],
    'param1_autough2': [['max_iterations', 'print_level', 'max_timesteps',
                        'max_duration', 'print_interval',
                        '_option_str', 'diff0', 'texp', 'be'],
                       ['2d'] * 2 + ['4d'] * 3 + ['24s'] + ['10.3e'] * 3],
    'param1': [['max_iterations', 'print_level', 'max_timesteps',
               'max_duration', 'print_interval', '_option_str', 'texp', 'be'],
              ['2d'] * 2 + ['4d'] * 3 + ['24s'] + ['10.3e'] * 2],
    'param2': [['tstart', 'tstop', 'const_timestep', 'max_timestep',
               'print_block', '', 'gravity', 'timestep_reduction', 'scale'],
              ['10.3e'] * 4 + ['5s', '5x'] + ['10.4e'] * 3],
    'param3': [['relative_error', 'absolute_error', 'pivot', 'upstream_weight',
               'newton_weight', 'derivative_increment'],
              ['10.4e'] * 6],
    '_more_option_str': [['_more_option_str'], ['21s']],
    'timestep': [['timestep'] * 8, ['10.4e'] * 8],
    'multi': [['num_components', 'num_equations', 'num_phases',
               'num_secondary_parameters', 'num_inc'],  ['5d'] * 5],
    'multi_autough2': [['num_components', 'num_equations', 'num_phases',
                        'num_secondary_parameters', 'eos'],  ['5d'] * 4 + ['4s']],
    'lineq': [['type', 'epsilon', 'max_iterations', 'gauss', 'num_orthog'],
             ['2d', '10.4e', '4d', '1d', '4d']],
    'default_incons': [['incon'] * 4, ['20.14e'] * 4],
    'output_times1': [['num_times_specified', 'num_times', 'max_timestep', 'time_increment'],
                     ['5d'] * 2 + ['10.4e'] * 2],
    'output_times2': [['time'] * 8, ['10.4e'] * 8],
    'relative_permeability': [['type', ''] + ['parameter'] * 7, ['5d', '5x'] + ['10.3e'] * 7],
    'capillarity': [['type', ''] + ['parameter'] * 7, ['5d', '5x'] + ['10.3e'] * 7], 
    'blocks': [['name', 'nseq', 'nadd', 'rocktype', 'volume',
               'ahtx', 'pmx', 'x', 'y', 'z'],
              ['5s', '5d', '5d', '5s'] + ['10.4e'] * 3 + ['10.3e'] * 3], 
    'connections': [['block1', 'block2', 'nseq', 'nad1', 'nad2',
                    'direction', 'distance1', 'distance2', 'area', 'dircos', 'sigma'],
                   ['5s'] * 2 + ['5d'] * 4 + ['10.4e'] * 3 + ['10.7f', '10.3e']],
    'generator': [['block', 'name', 'nseq', 'nadd', 'nads', 'ltab',
                  '', 'type', 'itab', 'gx', 'ex', 'hg', 'fg'],
                 ['5s'] * 2 + ['5d'] * 3 + ['5d', '5x', '4s', '1s'] + ['10.3e'] * 4], 
    'generation_times': [['time'] * 4, ['14.7e'] * 4], 
    'generation_rates': [['rate'] * 4, ['14.7e'] * 4], 
    'generation_enthalpy': [['enthalpy'] * 4, ['14.7e'] * 4],
    'short': [['', 'frequency'], ['5x', '2d']],
    'incon1': [['block', 'nseq', 'nadd', 'porosity'], ['5s'] + ['5d'] * 2 + ['15.9e']], 
    'incon2': [['incon'] * 4, ['20.14e'] * 4], 
    'solver': [['type', '', 'z_precond', '', 'o_precond', 'relative_max_iterations', 'closure'],
              ['1d', '2x', '2s', '3x', '2s'] + ['10.4e'] * 2], 
    'indom2': [['indom'] * 4, ['20.13e'] * 4], 
    'diffusion':  [['diff'] * 8, ['10.3e'] * 8],
    'selec1': [['int_selec'] * 16, ['5d'] * 16],
    'selec2': [['float_selec'] * 8, ['10.3e'] * 8], 
    'radii1': [['nrad'], ['5d']],
    'radii2': [['radius'] * 8, ['10.4e'] * 8],
    'equid' : [['nequ', '', 'dr'], ['5d', '5x', '10.4e']],
    'logar' : [['nlog', '', 'rlog', 'dr'], ['5d', '5x'] + ['10.4e'] * 2],
    'layer1': [['nlay'], ['5d']],
    'layer2': [['layer'] * 8, ['10.4e'] * 8],
    'xyz1'  : [['deg'], ['10.4e']],
    'xyz2'  : [['ntype', '', 'no', 'del'], ['2s', '3x', '5d', '10.4e']],
    'xyz3'  : [['deli'] * 8, ['10.4e'] * 8], 
    'minc'  : [['part', 'type', '', 'dual'], ['5s'] * 2 + ['5x', '5s']],
    'part1' : [['num_continua', 'nvol', 'where'] + ['spacing'] * 7,
               ['3d'] * 2 + ['-4s'] + ['10.4e'] * 7],
    'part2' : [['vol'] * 8, ['10.4e'] * 8]
    }

t2data_extra_precision_format_specification = {
    'rocks1': [['name', 'nad', 'density', 'porosity', 'k1', 'k2', 'k3',
               'conductivity', 'specific_heat'], 
              ['5s', '5d'] + ['15.8e'] * 7], 
    'rocks1.1': [['compressibility', 'expansivity', 'dry_conductivity',
                 'tortuosity', 'klinkenberg', 'xkd3', 'xkd4'],  ['15.8e'] * 7],
    'rocks1.2': [['type', ''] + ['parameter'] * 7, ['5d', '5x'] + ['15.8e'] * 7],
    'rocks1.3': [['type', ''] + ['parameter'] * 7, ['5d', '5x'] + ['15.8e'] * 7],
    'blocks': [['name', 'nseq', 'nadd', 'rocktype', 'volume',
               'ahtx', 'pmx', 'x', 'y', 'z'],
              ['5s', '5d', '5d', '5s'] + ['15.8e'] * 3 + ['15.8e'] * 3],
    'connections': [['block1', 'block2', 'nseq', 'nad1', 'nad2',
                    'direction', 'distance1', 'distance2', 'area', 'dircos', 'sigma'],
                   ['5s'] * 2 + ['5d'] * 4 + ['15.8e'] * 3 + ['15.8f', '15.8e']],
    'relative_permeability': [['type', ''] + ['parameter'] * 7, ['5d', '5x'] + ['15.8e'] * 7], 
    'capillarity': [['type', ''] + ['parameter'] * 7, ['5d', '5x'] + ['15.8e'] * 7],
    'generator': [['block', 'name', 'nseq', 'nadd', 'nads', 'ltab',
                  '', 'type', 'itab', 'gx', 'ex', 'hg', 'fg'],
                 ['5s'] * 2 + ['5d'] * 3 + ['5d', '5x', '4s', '1s'] + ['15.8e'] * 4], 
    'generation_times': [['time'] * 4, ['15.8e'] * 4], 
    'generation_rates': [['rate'] * 4, ['15.8e'] * 4], 
    'generation_enthalpy': [['enthalpy'] * 4, ['15.8e'] * 4]}

class t2data_parser(fixed_format_file):
    """Class for parsing TOUGH2 data file."""
    def __init__(self, filename, mode, read_function = default_read_function):
        super(t2data_parser,self).__init__(filename, mode,
                                           t2data_format_specification, read_function)

class t2_extra_precision_data_parser(fixed_format_file):
    """Class for parsing AUTOUGH2 extra-precision auxiliary data file."""
    def __init__(self, filename, mode, read_function = default_read_function):
        super(t2_extra_precision_data_parser,
              self).__init__(filename, mode,
                             t2data_extra_precision_format_specification,
                             read_function)

class fortran_unformatted_file(object):
    """Class for 'unformatted' binary file written by Fortran.  These are
    different from plain binary files in that the byte length of each
    'record' is written at its start and end.
    """

    def __init__(self, filename, mode):
        self.file = open(filename, mode)

    def close(self):
        self.file.close()

    def readrec(self, fmt):
        nb, = struct.unpack('i', self.file.read(4))
        packed = self.file.read(nb)
        self.file.read(4)
        return struct.unpack(fmt, packed)

    def writerec(self, fmt, val):
        nb = struct.calcsize(fmt)
        if isinstance(val, (tuple, list, np.ndarray)):
            packed = struct.pack(fmt,  *val)
        else: packed = struct.pack(fmt, val)
        head = struct.pack('i', nb)
        self.file.write(head)
        self.file.write(packed)
        self.file.write(head)

class t2generator(object):
    """TOUGH2 generator (source or sink)"""

    def __init__(self, name = '     ', block = '     ',
                 nseq = None, nadd = None, nads = None, type = 'MASS', 
                 ltab = 0, itab = '', gx = 0.0, ex = 0.0, hg = 0.0, fg = 0.0,
                 time = None, rate = None, enthalpy = None):
        if time is None: time = []
        if rate is None: rate = []
        if enthalpy is None: enthalpy = []
        self.name = name
        self.block = block
        self.nseq, self.nadd, self.nads = nseq, nadd, nads
        self.type = type
        self.ltab = ltab
        self.itab = itab
        self.gx = gx
        self.ex = ex
        self.hg = hg
        self.fg = fg
        self.time = time
        self.rate = rate
        self.enthalpy = enthalpy
    def __repr__(self): return self.block + ':' + self.name

default_parameters = {
    'max_iterations': None,
    'print_level': None,
    'max_timesteps': None,
    'max_duration': None,
    'print_interval': None, 
    '_option_str': '0' * 24,
    'option': np.zeros(25, int8),
    'diff0': None,
    'texp': None,
    'tstart': 0.0,
    'tstop': None,
    'const_timestep': 0.0,
    'timestep': [],
    'max_timestep': None,
    'print_block': None,
    'gravity': 0.0,
    'timestep_reduction': None,
    'scale': None,
    'relative_error': None,
    'absolute_error': None,
    'pivot': None,
    'upstream_weight': None,
    'newton_weight': None,
    'derivative_increment': None,
    'default_incons': []}

t2data_sections = [
    'SIMUL', 'ROCKS', 'PARAM', 'MOMOP', 'START', 'NOVER', 'RPCAP',
    'LINEQ', 'SOLVR', 'MULTI', 'TIMES', 'SELEC', 'DIFFU',
    'ELEME', 'CONNE', 'MESHM', 'GENER', 'SHORT', 'FOFT',
    'COFT', 'GOFT', 'INCON', 'INDOM']

t2_extra_precision_sections = ['ROCKS', 'ELEME', 'CONNE', 'RPCAP', 'GENER']

class t2data(object):
    """Class for TOUGH2 data."""
    def __init__(self, filename = '', meshfilename = '',
                 read_function = default_read_function):
        from copy import deepcopy
        self.filename = filename
        self.meshfilename = meshfilename
        self.title = ''
        self.simulator = ''
        self.parameter = deepcopy(default_parameters)
        self._more_option_str = '0' * 21,
        self.more_option = np.zeros(22, int8)
        self.multi = {}
        self.start = False
        self.relative_permeability = {}
        self.capillarity = {}
        self.lineq = {}
        self.output_times = {}
        self.grid = t2grid()
        self.generatorlist = []
        self.generator = {}
        self.short_output = {}
        self.incon = {}
        self.solver = {}
        self.history_block = []
        self.history_connection = []
        self.history_generator = []
        self.indom = {}
        self.noversion = False
        self.diffusion = []
        self.selection = {}
        self.meshmaker = []
        self._sections = []
        self.end_keyword = 'ENDCY'
        self._extra_precision, self._echo_extra_precision = [], True
        self.update_read_write_functions()
        self.read_function = read_function
        if self.filename: self.read(filename, meshfilename)

    def get_extra_precision(self): return self._extra_precision
    def set_extra_precision(self, value):
        if value is False: value = []
        elif value is True: value = t2_extra_precision_sections
        elif isinstance(value, str): value = [value]
        # check if removing any extra precision sections:
        for section in set(self._extra_precision) - set(value):
            self.insert_section(section)
        self._extra_precision = value
        self.update_read_write_functions()
    extra_precision = property(get_extra_precision, set_extra_precision)

    def get_echo_extra_precision(self): return self._echo_extra_precision
    def set_echo_extra_precision(self, value):
        if value != self._echo_extra_precision:
            if value is False: # remove previously echoed sections from section list
                for section in self._extra_precision:
                    self.delete_section(section)
            else: # add sections previously not echoed to section list
                for section in self._extra_precision:
                    self.insert_section(section)
            self._echo_extra_precision = value
            self.update_read_write_functions()
    echo_extra_precision = property(get_echo_extra_precision, set_echo_extra_precision)

    def update_read_write_functions(self):
        """Updates functions for reading and writing sections of data file."""

        self.read_fn = dict(zip(
                t2data_sections,
                [self.read_simulator,
                 self.read_rocktypes,
                 self.read_parameters,
                 self.read_more_options,
                 self.read_start,
                 self.read_noversion,
                 self.read_rpcap,
                 self.read_lineq,
                 self.read_solver,
                 self.read_multi,
                 self.read_times,
                 self.read_selection,
                 self.read_diffusion,
                 self.read_blocks,
                 self.read_connections,
                 self.read_meshmaker,
                 self.read_generators,
                 self.read_short_output,
                 self.read_history_blocks,
                 self.read_history_connections,
                 self.read_history_generators,
                 self.read_incons,
                 self.read_indom]))

        self.write_fn = dict(zip(
                t2data_sections,
                [self.write_simulator,
                 self.write_rocktypes,
                 self.write_parameters,
                 self.write_more_options,
                 self.write_start,
                 self.write_noversion,
                 self.write_rpcap,
                 self.write_lineq,
                 self.write_solver,
                 self.write_multi,
                 self.write_times,
                 self.write_selection,
                 self.write_diffusion,
                 self.write_blocks,
                 self.write_connections,
                 self.write_meshmaker,
                 self.write_generators,
                 self.write_short_output,
                 self.write_history_blocks,
                 self.write_history_connections,
                 self.write_history_generators,
                 self.write_incons,
                 self.write_indom]))

        skip_fn = dict(zip(
            t2_extra_precision_sections,
            [self.skip_rocktypes,
             self.skip_blocks,
             self.skip_connections,
             self.skip_rpcap,
             self.skip_generators]))

        if not self.echo_extra_precision:
            for section in self.extra_precision: self.read_fn[section] = skip_fn[section]

    def get_present_sections(self):
        """Returns a list of TOUGH2 section keywords for which there are
        corresponding data in the t2data object."""
        data_present = dict(zip(
            t2data_sections,
            [self.simulator,
             self.grid and self.grid.rocktypelist,
             self.parameter,
             np.any(self.more_option),
             self.start,
             self.noversion,
             self.relative_permeability or self.capillarity,
             self.lineq,
             self.solver,
             self.multi,
             self.output_times,
             self.selection,
             self.diffusion,
             self.grid,
             self.grid,
             self.meshmaker,
             self.generatorlist,
             self.short_output,
             self.history_block,
             self.history_connection,
             self.history_generator,
             self.incon,
             self.indom]))
        return [keyword for keyword in t2data_sections if data_present[keyword]]
    present_sections = property(get_present_sections)

    def insert_section(self, section):
        """Inserts a new section into the internal list of data file
        sections."""
        if section not in self._sections:
            i = self.section_insertion_index(section)
            self._sections.insert(i, section)

    def delete_section(self, section):
        """Deletes a section from the internal list of data file sections."""
        try: self._sections.remove(section)
        except ValueError: pass

    def section_insertion_index(self, section):
        """Determines an appropriate position to insert the specified section
        in the internal list of data file sections.
        """
        try:
            listindex = t2data_sections.index(section)
            if listindex == 0: return 0  # SIMUL section
            else:
                # first look for sections above the one specified,
                # and put new one just after the last found:
                for i in reversed(range(listindex)):
                    try:
                        section_index = self._sections.index(t2data_sections[i])
                        return section_index + 1
                    except ValueError: pass
                # look for sections below the one specified,
                # and put new one just before the first found:
                for i in range(listindex, len(t2data_sections)):
                    try:
                        section_index = self._sections.index(t2data_sections[i])
                        return section_index
                    except ValueError: pass
                return len(self._sections)
        except ValueError: return len(self._sections)

    def update_sections(self):
        """Updates internal section list, based on which properties are present."""
        present = self.present_sections
        missing = [keyword for keyword in present if keyword not in self._sections]
        for keyword in missing: self.insert_section(keyword)
        extra = [keyword for keyword in self._sections if keyword not in present]
        for keyword in extra: self.delete_section(keyword)

    def __repr__(self): return self.title

    def run(self, save_filename = '', incon_filename = '', simulator = 'AUTOUGH2_2',
            silent = False, output_filename = ''):
        """Runs simulation using TOUGH2 or AUTOUGH2.  It's assumed that the
        data object has been written to file using write().  For
        AUTOUGH2, if the filenames for the save file or initial
        conditions file are not specified, they are constructed by
        changing the extensions of the data filename.  Set silent to
        True to suppress screen output. The output_filename applies
        only to TOUGH2, and specifies the name of the main output
        listing file."""
        if self.filename:
            from os.path import splitext, basename
            from os import devnull, remove
            from subprocess import call
            datbase, ext = splitext(self.filename)
            if (self.type == 'AUTOUGH2'):
                if save_filename == '': save_filename = datbase + '.save'
                if incon_filename == '': incon_filename = datbase + '.incon'
                savebase, ext = splitext(save_filename)
                inconbase, ext = splitext(incon_filename)
                runfilename = datbase + '_' + basename(simulator) + '.in'
                open(runfilename, 'w').write('\n'.join([savebase, inconbase, datbase]))
                infile = open(runfilename, 'r')
                cmd = [simulator]
                if silent: outfile = open(devnull, 'w')
                else: outfile = None
                # run AUTOUGH2:
                call(cmd, stdin = infile, stdout = outfile)
                infile.close()
                remove(runfilename)
            else: # run TOUGH2 (need to specify simulator executable name)
                cmd = [simulator]
                infile = open(self.filename, 'r')
                if silent: outfile = None
                else:
                    if output_filename == '': outfilename = datbase + '.listing'
                    else: outfilename = output_filename
                    outfile = open(outfilename, 'w')
                call(cmd, stdin = infile, stdout = outfile)

    def get_type(self):
        """Returns type (TOUGH2 or AUTOUGH2) based on whether the simulator
        has been set."""
        if self.simulator: return 'AUTOUGH2'
        else: return 'TOUGH2'

    def set_type(self, value):
        """Sets type (TOUGH2 or AUTOUGH2), and runs conversion if needed (with
        default options)."""
        if value in ['AUTOUGH2', 'TOUGH2']:
            oldtype = self.type
            if oldtype != value:
                if oldtype == 'AUTOUGH2': self.convert_to_TOUGH2()
                elif oldtype == 'TOUGH2': self.convert_to_AUTOUGH2()
        else: raise Exception('Data file type ' + value + ' is not supported.')

    type = property(get_type, set_type)

    def get_extra_precision_filename(self):
        """Returns name of extra precision data file name, based on the name
        of the data file."""
        from os.path import splitext
        base, ext = splitext(self.filename)
        if base[0].isupper(): pext = 'PDAT'
        else: pext = 'pdat'
        return '.'.join((base, pext))
    extra_precision_filename = property(get_extra_precision_filename)

    def get_num_generators(self):
        return len(self.generatorlist)
    num_generators = property(get_num_generators)

    def generator_index(self, blocksourcenames):
        """Returns index of generator with specified tuple of block and source
        names."""
        if blocksourcenames in self.generator:
            return self.generatorlist.index(self.generator[blocksourcenames])
        else: return None

    def total_generation(self, type = 'MASS', name = ''):
        """Returns array containing total generation in each block of the
        specified generator type and name.  The name parameter
        specifies a regular expression to be matched.
        """
        import re
        tg = np.zeros(self.grid.num_blocks, float64)
        gens = [g for g in self.generatorlist if ((type == g.type) and re.search(name, g.name))]
        for g in gens: tg[self.grid.block_index(g.block)] += g.gx
        return tg

    def specific_generation(self, type = 'MASS', name = ''):
        """Returns array containing total specific generation (i.e. generation
        per unit volume) in each block of the specified generator type
        and name.  The name parameter specifies a regular expression
        to be matched.
        """
        import re
        tg = np.zeros(self.grid.num_blocks, float64)
        gens = [g for g in self.generatorlist if
                ((g.type == type) and (re.search(name, g.name)))]
        for g in gens:
            blkindex = self.grid.block_index(g.block)
            tg[blkindex] += g.gx / self.grid.blocklist[blkindex].volume
        return tg

    def read_title(self, infile):
        """Reads simulation title"""
        infile.read_value_line(self.__dict__, 'title')

    def write_title(self, outfile):
        outfile.write(self.title.strip() + '\n')

    def read_simulator(self, infile):
        """Reads simulator and EOS type.  If the SIMUL section is present,
        check for extra precision data in auxiliary file, and modify
        functions accordingly for reading the main data file.
        """
        infile.read_value_line(self.__dict__, 'simulator')
        if self.type == 'AUTOUGH2': self.read_extra_precision()

    def write_simulator(self, outfile):
        if self.simulator:
            outfile.write('SIMUL\n')
            outfile.write(self.simulator.strip() + '\n')

    def read_rocktypes(self, infile):
        """Reads grid rock types"""
        self.grid.rocktypelist = []
        self.grid.rocktype = {}
        line = padstring(infile.readline())
        while line.strip():
            [name, nad, density, porosity,
             k1, k2, k3, conductivity, specific_heat] = infile.parse_string(line, 'rocks1')
            self.grid.add_rocktype(rocktype(name, nad,
                                            density, porosity,
                                            [k1, k2, k3],
                                            conductivity, specific_heat))
            if nad is None: nad = 0
            if nad >= 1: # additional lines:
                infile.read_value_line(self.grid.rocktype[name].__dict__, 'rocks1.1')
                if nad >=2 :
                    vals = infile.read_values('rocks1.2')
                    self.grid.rocktype[name].relative_permeability['type'] = vals[0]
                    self.grid.rocktype[name].relative_permeability['parameters'] = vals[2: -1]
                    vals = infile.read_values('rocks1.3')
                    self.grid.rocktype[name].capillarity['type'] = vals[0]
                    self.grid.rocktype[name].capillarity['parameters'] = vals[2: -1]
            line = padstring(infile.readline())

    def skip_rocktypes(self, infile):
        """Skips rock type section"""
        while infile.readline().strip(): pass

    def write_rocktypes(self, outfile):
        outfile.write('ROCKS\n')
        for rt in self.grid.rocktypelist:
            vals = [rt.name, rt.nad, rt.density, rt.porosity] + \
                   list(rt.permeability) + [rt.conductivity, rt.specific_heat]
            outfile.write_values(vals, 'rocks1')
            if rt.nad is not None:
                if rt.nad >= 1:
                    outfile.write_value_line(rt.__dict__, 'rocks1.1')
                    if rt.nad >= 2:
                        vals = [rt.relative_permeability['type'], None] + \
                               rt.relative_permeability['parameters']
                        outfile.write_values(vals, 'rocks1.2')
                        vals = [rt.capillarity['type'], None] + rt.capillarity['parameters']
                        outfile.write_values(vals, 'rocks1.2')
        outfile.write('\n')

    def read_parameters(self, infile):
        """Reads simulation parameters"""
        spec = ['param1', 'param1_autough2'][self.type == 'AUTOUGH2']
        infile.read_value_line(self.parameter, spec)
        mops = self.parameter['_option_str'].rstrip().ljust(24).replace(' ', '0')
        self.parameter['option'] = np.array([0] + [int(mop) for mop in mops], int8)
        infile.read_value_line(self.parameter, 'param2')
        if (self.parameter['print_block'] is not None) and \
           (self.parameter['print_block'].strip() == ''):
            self.parameter['print_block'] = None
        self.read_timesteps(infile)
        infile.read_value_line(self.parameter, 'param3')
        for val in infile.read_values('default_incons'):
            self.parameter['default_incons'].append(val)
        self.parameter['default_incons'] = trim_trailing_nones(self.parameter['default_incons'])
        # read any additional lines of default incons:
        more = True
        while more:
            line = padstring(infile.readline())
            if line.strip():
                section = any([line.startswith(keyword) for keyword in t2data_sections])
                if section: more = False
                else:
                    more_incons = infile.parse_string(line, 'default_incons')
                    more_incons = trim_trailing_nones(more_incons)
                    self.parameter['default_incons'] += more_incons
            else: more, line = False, None
        return line

    def write_parameters(self, outfile):
        outfile.write('PARAM\n')
        from copy import copy
        paramw = copy(self.parameter)
        if paramw['print_block'] is not None:
            paramw['print_block'] = unfix_blockname(paramw['print_block'])
        self.parameter['_option_str'] = ''.join([str(m) for m in self.parameter['option'][1:]])
        spec = ['param1', 'param1_autough2'][self.type == 'AUTOUGH2']
        outfile.write_value_line(self.parameter, spec)
        outfile.write_value_line(paramw, 'param2')
        self.write_timesteps(outfile)
        outfile.write_value_line(self.parameter, 'param3')
        num_vars = len(self.parameter['default_incons'])
        if num_vars > 0:
            nlines = int(ceil(num_vars / 4.))
            for i in range(nlines):
                i1, i2 = i * 4, min((i + 1) * 4, num_vars)
                vals = list(self.parameter['default_incons'][i1: i2])
                if len(vals) < 4: vals += [None] * (4 - len(vals))
                outfile.write_values(vals, 'default_incons')
        else: outfile.write('\n')

    def read_more_options(self, infile):
        """Reads additional parameter options"""
        infile.read_value_line(self.__dict__, '_more_option_str')
        momops = self._more_option_str.rstrip().ljust(21).replace(' ', '0')
        self.more_option = np.array([0] + [int(mop) for mop in momops], int8)
        
    def write_more_options(self, outfile):
        """Writes additional parameter options"""
        outfile.write('MOMOP\n')
        self._more_option_str = ''.join([str(m) for m in self.more_option[1:]])
        outfile.write_value_line(self.__dict__, '_more_option_str')
        
    def read_timesteps(self, infile):
        """Reads time step sizes from file"""
        if self.parameter['const_timestep'] >= 0.0:
            self.parameter['timestep'] = [self.parameter['const_timestep']]
        else:
            nlines = -int(self.parameter['const_timestep'])
            self.parameter['timestep'] = []
            for i in range(nlines):
                for val in infile.read_values('timestep'):
                    if val is not None: self.parameter['timestep'].append(val)

    def write_timesteps(self, outfile):
        if self.parameter['const_timestep'] < 0.0:
            nlines = -int(self.parameter['const_timestep'])
            for i in range(nlines):
                i1, i2 = i * 8, min((i + 1) * 8, len(self.parameter['timestep']))
                vals = self.parameter['timestep'][i1: i2]
                if len(vals) < 8: vals += [None] * (8 - len(vals))
                outfile.write_values(vals, 'timestep')

    def read_multi(self, infile):
        """Reads EOS parameters"""
        spec = ['multi', 'multi_autough2'][self.type == 'AUTOUGH2']
        infile.read_value_line(self.multi, spec)
        if 'eos' in self.multi: self.multi['eos'] = self.multi['eos'].strip()

    def write_multi(self, outfile):
        if self.multi != {}:
            outfile.write('MULTI\n')
            spec = ['multi', 'multi_autough2'][self.type == 'AUTOUGH2']
            outfile.write_value_line(self.multi, spec)

    def read_start(self, infile):
        """Sets start parameter"""
        self.start = True

    def write_start(self, outfile):
        if self.start: outfile.write('START\n')

    def read_rpcap(self, infile):
        """Reads relative permeability and capillarity parameters"""
        vals = infile.read_values('relative_permeability')
        self.relative_permeability['type'] = vals[0]
        self.relative_permeability['parameters'] = vals[2:]
        vals = infile.read_values('capillarity')
        self.capillarity['type']= vals[0]
        self.capillarity['parameters'] = vals[2:]

    def skip_rpcap(self, infile):
        """Skips relative permeability and capillarity parameter section."""
        for i in range(2): infile.readline()

    def write_rpcap(self, outfile):
        if self.relative_permeability:
            outfile.write('RPCAP\n')
            vals = [self.relative_permeability['type'], None] + \
                   self.relative_permeability['parameters']
            outfile.write_values(vals, 'relative_permeability')
            vals = [self.capillarity['type'], None] + \
                   self.capillarity['parameters']
            outfile.write_values(vals, 'capillarity')

    def read_lineq(self, infile):
        """Reads linear equation parameters (AUTOUGH2)"""
        infile.read_value_line(self.lineq, 'lineq')

    def write_lineq(self, outfile):
        if self.lineq:
            outfile.write('LINEQ\n')
            outfile.write_value_line(self.lineq, 'lineq')

    def read_solver(self, infile):
        """Reads linear equation parameters (TOUGH2)"""
        infile.read_value_line(self.solver, 'solver')

    def write_solver(self, outfile):
        if self.solver:
            outfile.write('SOLVR\n')
            outfile.write_value_line(self.solver, 'solver')

    def read_blocks(self, infile):
        """Reads grid blocks"""
        self.grid.block, self.grid.blocklist = {}, []
        line = padstring(infile.readline())
        while line.strip():
            [name, nseq, nadd, rockname,
             volume, ahtx, pmx, x, y, z] = infile.parse_string(line, 'blocks')
            name = fix_blockname(name)
            if rockname in self.grid.rocktype:
                rocktype = self.grid.rocktype[rockname]
            elif rockname.strip() == '' and self.grid.num_rocktypes > 0:
                rocktype = self.grid.rocktypelist[0] # default
            else:
                try: # check if rocktype index specified:
                    rockindex = int(rockname) - 1
                    rocktype = self.grid.rocktypelist[rockindex]
                except:
                    raise RuntimeError("Unknown rocktype " + rockname + " in block " + name)
            if (x is not None) and (y is not None) and (z is not None):
                centre = np.array([x, y, z])
            else: centre = None
            if nseq == 0: nseq = None
            if nadd == 0: nadd = None
            self.grid.add_block(t2block(name, volume, rocktype,
                                        centre = centre, ahtx = ahtx,
                                        pmx = pmx, nseq = nseq, nadd = nadd))
            line = padstring(infile.readline())

    def skip_blocks(self, infile):
        """Skips blocks section in file"""
        while infile.readline().strip(): pass

    def write_blocks(self, outfile):
        outfile.write('ELEME\n')
        if self.grid.num_blocks > 0:
            from copy import copy
            for blk in self.grid.blocklist:
                blkw = copy(blk.__dict__)
                blkw['name'] = unfix_blockname(blkw['name'])
                if blk.centre is None: outfile.write_value_line(blkw, 'blocks')
                else:
                    vals = [blkw['name'], blk.nseq, blk.nadd,
                            blk.rocktype.name, blk.volume,
                            blk.ahtx, blk.pmx] + list(blk.centre)
                    outfile.write_values(vals, 'blocks')
        outfile.write('\n')

    def read_connections(self, infile):
        """Reads grid connections"""
        self.grid.connectionlist, self.grid.connection = [], {}
        line = padstring(infile.readline())
        while line.strip() and not line.startswith('+++'):
            [name1, name2, nseq,
             nad1, nad2, isot, d1, d2,
             areax, betax, sigx] = infile.parse_string(line, 'connections')
            name1, name2 = fix_blockname(name1), fix_blockname(name2)
            if nseq == 0: nseq = None
            if nad1 == 0: nad1 = None
            if nad2 == 0: nad2 = None
            self.grid.add_connection(t2connection([self.grid.block[name1],
                                                   self.grid.block[name2]],
                                                  isot, [d1, d2], areax, betax,
                                                  sigx, nseq, nad1, nad2))
            line = padstring(infile.readline())

    def skip_connections(self, infile):
        """Skips connections section in file"""
        while infile.readline().strip(): pass

    def write_connections(self, outfile):
        outfile.write('CONNE\n')
        if self.grid.num_connections > 0:
            for con in self.grid.connectionlist:
                vals = [unfix_blockname(con.block[0].name),
                        unfix_blockname(con.block[1].name),
                        con.nseq, con.nad1, con.nad2, con.direction] + \
                    con.distance + [con.area, con.dircos, con.sigma]
                outfile.write_values(vals, 'connections')
        outfile.write('\n')

    def add_generator(self, generator = None):
        """Adds a generator."""
        if generator is None: generator = t2generator()
        self.generatorlist.append(generator)
        self.generator[(generator.block, generator.name)] = self.generatorlist[-1]

    def delete_generator(self, blocksourcenames):
        i = self.generator_index(blocksourcenames)
        del self.generator[blocksourcenames]
        del self.generatorlist[i]

    def clear_generators(self):
        self.generator.clear()
        del self.generatorlist[:]

    def delete_orphan_generators(self):
        """Deletes any generators specified in blocks which do not exist in
        the grid."""
        delgenindices, delgens = [], set([])
        # Delete items from the generator list and dictionary
        # separately (rather than using self.delete_generator() in
        # case there are multiple generators in the list with the same
        # dictionary key.
        for i, gen in enumerate(self.generatorlist):
            if not (gen.block in self.grid.block):
                delgenindices.append(i)
                delgens.add((gen.block, gen.name))
        for bg in delgens: del self.generator[bg]
        for i in reversed(delgenindices):
            del self.generatorlist[i]

    def read_generator(self, line, infile):
        """Returns generator read from line in file"""
        [block, name, nseq, nadd, nads,
         ltab, empty, gentype, itab,
         gx, ex, hg, fg] = infile.parse_string(line, 'generator')
        block, name = fix_blockname(block), fix_blockname(name)
        time, rate, enthalpy = [], [], []
        if ltab and gentype != 'DELV':
            ntimes = abs(ltab)
            if ntimes > 1:
                nlines = int(ceil(ntimes / 4.))
                for i in range(nlines):
                    for val in infile.read_values('generation_times'):
                        if val is not None: time.append(val)
                for i in range(nlines):
                    for val in infile.read_values('generation_rates'):
                        if val is not None: rate.append(val)
                if itab.strip():
                    for i in range(nlines):
                        for val in infile.read_values('generation_enthalpy'):
                            if val is not None: enthalpy.append(val)
        return t2generator(name = name, block = block,
                           nseq = nseq, nadd = nadd, nads = nads,
                           type = gentype, ltab = ltab, itab = itab,
                           gx = gx, ex = ex, hg = hg, fg = fg,
                           time = time, rate = rate, enthalpy = enthalpy)

    def write_generator(self, gen, outfile):
        from copy import copy
        genw = copy(gen.__dict__)
        genw['name']  = unfix_blockname(genw['name'])
        genw['block'] = unfix_blockname(genw['block'])
        outfile.write_value_line(genw, 'generator')
        if gen.ltab and gen.type != 'DELV': ntimes = abs(gen.ltab)
        else: ntimes = 1
        if ntimes > 1:
            nlines = int(ceil(ntimes / 4.))
            for i in range(nlines):
                i1, i2 = i * 4, min((i + 1) * 4, ntimes)
                vals = list(gen.time[i1: i2])
                if len(vals) < 4: vals += [None] * (4 - len(vals))
                outfile.write_values(vals, 'generation_times')
            for i in range(nlines):
                i1, i2 = i * 4, min((i + 1) * 4, ntimes)
                vals = list(gen.rate[i1: i2])
                if len(vals) < 4: vals += [None] * (4 - len(vals))
                outfile.write_values(vals, 'generation_rates')
            if gen.enthalpy:
                for i in range(nlines):
                    i1, i2 = i * 4, min((i + 1) * 4, ntimes)
                    vals = list(gen.enthalpy[i1: i2])
                    if len(vals) < 4: vals += [None] * (4 - len(vals))
                    outfile.write_values(vals, 'generation_enthalpy')

    def read_generators(self, infile):
        """Reads generators from file"""
        self.generatorlist = []
        line = infile.readline()
        while line.strip():
            self.add_generator(self.read_generator(line, infile))
            line = infile.readline()

    def skip_generators(self, infile):
        """Skips generator section in file"""
        while infile.readline().strip(): pass

    def write_generators(self, outfile):
        if self.generatorlist:
            outfile.write('GENER\n')
            for generator in self.generatorlist:
                self.write_generator(generator, outfile)
            outfile.write('\n')

    def read_times(self, infile):
        """Reads output times from file"""
        infile.read_value_line(self.output_times, 'output_times1')
        self.output_times['time'] = []
        nlines = int(ceil(self.output_times['num_times_specified'] / 8.))
        for i in range(nlines):
            for val in infile.read_values('output_times2'):
                if val is not None: self.output_times['time'].append(val)

    def write_times(self, outfile):
        if self.output_times:
            outfile.write('TIMES\n')
            outfile.write_value_line(self.output_times, 'output_times1')
            nlines = int(ceil(self.output_times['num_times_specified'] / 8.))
            for i in range(nlines):
                i1, i2 = i * 8, min((i + 1) * 8, len(self.output_times['time']))
                vals = self.output_times['time'][i1: i2]
                if len(vals) < 8: vals += [None] * (8 - len(vals))
                outfile.write_values(vals, 'output_times2')

    def read_incons(self, infile):
        """Reads initial conditions from file"""
        line = infile.readline()
        while line.strip():
            [blockname,  nseq,  nadd,
             porosity] = infile.parse_string(line,  'incon1')
            blockname = fix_blockname(blockname)
            variables = infile.read_values('incon2')
            variables = trim_trailing_nones(variables)
            if nseq == 0: nseq = None
            if nadd == 0: nadd = None
            if nseq is None:
                self.incon[blockname] = [porosity, variables]
            else:
                self.incon[blockname] = [porosity, variables, nseq, nadd]
            line = infile.readline()

    def write_incons(self, outfile):
        if self.incon:
            outfile.write('INCON\n')
            for blk in self.grid.blocklist:
                blkname = blk.name
                try:
                    inc = self.incon[blkname]
                except:
                    continue
                if len(inc) >= 4: nseq, nadd = inc[2], inc[3]
                else: nseq, nadd = None, None
                vals = [unfix_blockname(blkname), nseq, nadd, inc[0]]
                outfile.write_values(vals, 'incon1')
                outfile.write_values(inc[1], 'incon2')
            outfile.write('\n')

    def read_short_blocks(self, infile):
        """Reads short output blocks"""
        self.short_output['block'] = []
        badblocks = []
        more = True
        while more:
            line = infile.readline()
            if line.strip():
                if line[0: 5] in ['ELEME', 'CONNE', 'GENER']: more = False
                else:
                    blockname = fix_blockname(line[0: 5])
                    if blockname in self.grid.block:
                        self.short_output['block'].append(self.grid.block[blockname])
                    else: badblocks.append(blockname)
            else: more = False
        if len(badblocks) > 0:
            print('Short output blocks', badblocks, 'do not exist and will be ignored.')
        return line

    def read_short_connections(self, infile):
        """Reads short output connections"""
        self.short_output['connection'] = []
        badcons = []
        more = True
        while more:
            line = infile.readline()
            if line.strip():
                if line[0: 5] in ['ELEME', 'CONNE', 'GENER']:
                    more = False
                else:
                    blknames = (fix_blockname(line[0: 5]), fix_blockname(line[5: 10]))
                    if blknames in self.grid.connection:
                        self.short_output['connection'].append(self.grid.connection[blknames])
                    else: badcons.append(blknames)
            else: more = False
        if len(badcons) > 0:
            print('Short output connections', badcons, 'do not exist and will be ignored.')
        return line

    def read_short_generators(self, infile):
        """Reads short output generators"""
        self.short_output['generator'] = []
        badgens = []
        more = True
        while more:
            line = infile.readline()
            if line.strip():
                if line[0: 5] in ['ELEME', 'CONNE', 'GENER']:
                    more = False
                else:
                    blksourcenames = (fix_blockname(line[0: 5]), fix_blockname(line[5: 10]))
                    if blksourcenames in self.generator:
                        self.short_output['generator'].append(self.generator[blksourcenames])
                    else: badgens.append(blksourcenames)
            else: more = False
        if len(badgens) > 0:
            print('Short output generators', badgens, 'do not exist and will be ignored.')
        return line

    def read_short_output(self, infile, headerline):
        """Reads short output specifications from file.  'headerline' is
        passed in to read the frequency parameter."""
        vals = infile.parse_string(headerline, 'short')
        if len(vals) > 1: self.short_output['frequency'] = vals[1]
        read_fn = {
            'ELEME': self.read_short_blocks,
            'CONNE': self.read_short_connections,
            'GENER': self.read_short_generators}
        more = True
        line = infile.readline()
        while more:
            if line.strip():
                keyword = line[0: 5]
                line = read_fn[keyword](infile)
            else: more = False

    def write_short_output(self, outfile):
        if self.short_output:
            outfile.write('SHORT')
            if 'frequency' in self.short_output:
                if self.short_output['frequency']:
                    outfile.write('%2d' % self.short_output['frequency'])
            outfile.write('\n')
            if 'block' in self.short_output:
                outfile.write('ELEME\n')
                for block in self.short_output['block']:
                    outfile.write(unfix_blockname(block.name) + '\n')
            if 'connection' in self.short_output:
                outfile.write('CONNE\n')
                for con in self.short_output['connection']:
                    outfile.write(unfix_blockname(con.block[0].name) + 
                                  unfix_blockname(con.block[1].name) + '\n')
            if 'generator' in self.short_output:
                outfile.write('GENER\n')
                for gen in self.short_output['generator']:
                    outfile.write(unfix_blockname(gen.block) + \
                                  unfix_blockname(gen.name) + '\n')
            outfile.write('\n')

    def read_history_blocks(self, infile):
        """Reads history blocks (TOUGH2)"""
        self.history_block = []
        badblocks = []
        line = infile.readline()
        if self.grid.num_blocks > 0:
            while line.strip():
                blockname = fix_blockname(line[0: 5])
                if blockname in self.grid.block:
                    self.history_block.append(self.grid.block[blockname])
                else: badblocks.append(blockname)
                line = infile.readline()
            if len(badblocks) > 0:
                print('History blocks', badblocks, 'do not exist and will be ignored.')
        else:
            # no grid- don't check blocks; and store names rather than t2blocks
            while line.strip():
                blockname = fix_blockname(line[0: 5])
                self.history_block.append(blockname)
                line = infile.readline()

    def read_history_connections(self, infile):
        """Reads history connections (TOUGH2)"""
        self.history_connection = []
        badcons = []
        line = infile.readline()
        if self.grid.num_blocks > 0:
            while line.strip():
                blknames = (fix_blockname(line[0: 5]), fix_blockname(line[5: 10]))
                if blknames in self.grid.connection:
                    self.history_connection.append(self.grid.connection[blknames])
                else: badcons.append(blknames)
                line = infile.readline()
            if len(badcons) > 0:
                print('History connections', badcons, 'do not exist and will be ignored.')
        else: # no grid
            while line.strip():
                blknames = (fix_blockname(line[0: 5]), fix_blockname(line[5: 10]))
                self.history_connection.append(blknames)
                line = infile.readline()

    def read_history_generators(self, infile):
        """Reads history generators (TOUGH2)"""
        self.history_generator = []
        badgens = []
        line = infile.readline()
        if self.grid.num_blocks > 0:
            while line.strip():
                blockname = fix_blockname(line[0: 5])
                if blockname in self.grid.block:
                    self.history_generator.append(self.grid.block[blockname])
                else: badgens.append(blockname)
                line = infile.readline()
            if len(badgens) > 0:
                print('History generator blocks', badgens, 'do not exist and will be ignored.')
        else: # no grid
            while line.strip():
                blockname = fix_blockname(line[0: 5])
                self.history_generator.append(blockname)
                line = infile.readline()

    def write_history_blocks(self, outfile):
        if self.history_block:
            outfile.write('FOFT\n')
            for blk in self.history_block:
                if isinstance(blk, str): blkname = blk
                else: blkname = blk.name
                outfile.write(unfix_blockname(blkname) + '\n')
            outfile.write('\n')

    def write_history_connections(self, outfile):
        if self.history_connection:
            outfile.write('COFT\n')
            for con in self.history_connection:
                if isinstance(con, tuple): cname = con
                else: cname = tuple([blk.name for blk in con.block])
                outfile.write(unfix_blockname(cname[0]) + unfix_blockname(cname[1]) + '\n')
            outfile.write('\n')

    def write_history_generators(self, outfile):
        if self.history_generator:
            outfile.write('GOFT\n')
            for blk in self.history_generator:
                if isinstance(blk, str): blkname = blk
                else: blkname = blk.name
                outfile.write(unfix_blockname(blkname) + '\n')
            outfile.write('\n')

    def read_indom(self, infile):
        """Reads rock-specific initial conditions from file"""
        line = infile.readline()
        while line.strip():
            rockname = line[0: 5]
            variables = infile.read_values('indom2')
            variables = trim_trailing_nones(variables)
            self.indom[rockname] = variables
            line = infile.readline()

    def write_indom(self, outfile):
        if self.indom:
            outfile.write('INDOM\n')
            for rockname, inc in self.indom.items():
                outfile.write(rockname + '\n')
                outfile.write_values(inc, 'indom2')
            outfile.write('\n')

    def read_noversion(self, infile):
        """Sets noversion parameter"""
        self.noversion = True

    def write_noversion(self, outfile):
        if self.noversion: outfile.write('NOVER\n')

    def read_diffusion(self, infile):
        """Reads diffusion coefficients from file"""
        if ('num_components' in self.multi) and ('num_phases' in self.multi):
            for comp in range(self.multi['num_components']):
                diffs = infile.read_values('diffusion')[0: self.multi['num_phases']]
                self.diffusion.append(diffs)
        else: print('Unable to read DIFFU block: no MULTI block specified.')

    def write_diffusion(self, outfile):
        if self.diffusion:
            outfile.write('DIFFU\n')
            for comp in self.diffusion: outfile.write_values(comp, 'diffusion')

    def read_selection(self, infile):
        """Reads selection parameters from file"""
        int_selec = infile.read_values('selec1')
        self.selection['integer'] = int_selec
        nlines = int_selec[0]
        float_selec = []
        for i in range(nlines): float_selec += infile.read_values('selec2')
        self.selection['float'] = float_selec

    def write_selection(self, outfile):
        if self.selection:
            outfile.write('SELEC\n')
            outfile.write_values(self.selection['integer'], 'selec1')
            nlines = self.selection['integer'][0]
            for i in range(nlines):
                i1, i2 = i * 8, min((i + 1) * 8, len(self.selection['float']))
                vals = self.selection['float'][i1: i2]
                if len(vals) < 8: vals += [None] * (8 - len(vals))
                outfile.write_values(vals, 'selec2')

    def read_meshmaker(self, infile):
        """Reads meshmaker data"""
        read_fn = {
            'RZ2D': self.read_meshmaker_rz2d,
            'XYZ': self.read_meshmaker_xyz,
            'MINC': self.read_meshmaker_minc}
        more = True
        while more:
            line = infile.readline()
            if line.strip():
                keyword = line[0: 5].strip()
                if keyword in read_fn: read_fn[keyword](infile)
            else: more = False

    def write_meshmaker(self, outfile):
        if self.meshmaker:
            outfile.write('MESHMAKER\n')
            write_fn = {
                'rz2d': self.write_meshmaker_rz2d,
                'xyz': self.write_meshmaker_xyz,
                'minc': self.write_meshmaker_minc}
            for (stype, section) in self.meshmaker:
                write_fn[stype.lower()](section, outfile)
            outfile.write('\n')

    def read_meshmaker_rz2d(self, infile):
        """Reads RZ2D meshmaker data"""
        section = ('rz2d', [])
        more = True
        while more:
            line = infile.readline()
            keyword = line[0: 5].strip()
            subsection = {}
            if keyword == 'RADII':
                nrad = infile.read_values('radii1')[0]
                subsection['radii'] = []
                nlines = int(ceil(nrad / 8.))
                for i in range(nlines):
                    for val in infile.read_values('radii2'):
                        if val is not None: subsection['radii'].append(val)
            elif keyword == 'EQUID': infile.read_value_line(subsection, 'equid')
            elif keyword == 'LOGAR': infile.read_value_line(subsection, 'logar')
            elif keyword == 'LAYER':
                nlayers = infile.read_values('layer1')[0]
                nlines = int(ceil(nlayers / 8.))
                layer = []
                for i in range(nlines): layer += infile.read_values('layer2')
                subsection['layer'] = layer[0: nlayers]
                more = False # LAYER indicates end of RZ2D
            if subsection: section[1].append((keyword.lower(), subsection))
        self.meshmaker.append(section)

    def write_meshmaker_rz2d(self, section, outfile):
        outfile.write('RZ2D\n')
        for stype, subsection in section:
            outfile.write(stype.upper() + '\n')
            if stype == 'radii':
                nrad = len(subsection['radii'])
                outfile.write_values([nrad], 'radii1')
                nlines = int(ceil(nrad / 8.))
                for i in range(nlines):
                    i1, i2 = i * 8, min((i + 1) * 8, nrad)
                    vals = subsection['radii'][i1: i2]
                    if len(vals) < 8: vals += [None] * (8 - len(vals))
                    outfile.write_values(vals, 'radii2')
            elif stype == 'equid': outfile.write_value_line(subsection, 'equid')
            elif stype == 'logar': outfile.write_value_line(subsection, 'logar')
            elif stype == 'layer':
                nlayers = len(subsection['layer'])
                outfile.write_values([nlayers], 'layer1')
                nlines = int(ceil(nlayers / 8.))
                for i in range(nlines):
                    i1, i2 = i * 8, min((i + 1) * 8, nrad)
                    vals = subsection['layer'][i1: i2]
                    if len(vals) < 8: vals += [None] * (8 - len(vals))
                    outfile.write_values(vals, 'layer2')

    def read_meshmaker_xyz(self, infile):
        """Reads XYZ meshmaker data"""
        deg = infile.read_values('xyz1')[0]
        section = ('xyz', [deg])
        more = True
        while more:
            subsection = {}
            line = infile.readline()
            if line.strip():
                [subsection['ntype'], blank,
                 subsection['no'], subsection['del']] = infile.parse_string(line, 'xyz2')
                if subsection['del'] == 0:
                    nlines = int(ceil(subsection['no'] / 8.))
                    deli = []
                    for i in range(nlines):
                        deli += infile.read_values('xyz3')
                    subsection['deli'] = deli[0: subsection['no']]
                section[1].append(subsection)
            else: more = False
        self.meshmaker.append(section)

    def write_meshmaker_xyz(self, section, outfile):
        outfile.write('XYZ\n')
        deg = section[0]
        outfile.write_values([deg], 'xyz1')
        for subsection in section[1:]:
            outfile.write_value_line(subsection, 'xyz2')
            if subsection['del'] == 0:
                nlines = int(ceil(subsection['no'] / 8.))
                for i in range(nlines):
                    i1, i2 = i * 8, min((i + 1) * 8, subsection['no'])
                    vals = subsection['deli'][i1: i2]
                    if len(vals) < 8: vals += [None] * (8 - len(vals))
                    outfile.write_values(vals, 'xyz3')
        outfile.write('\n')

    def read_meshmaker_minc(self, infile):
        """Reads MINC meshmaker data"""
        line = infile.readline().strip()
        keyword = line[0: 5].strip()
        if keyword == 'PART':
            subsection = {}
            [part, subsection['type'],
             dummy, subsection['dual']] = infile.parse_string(line, 'minc')
            vals = infile.read_values('part1')
            subsection['num_continua'], nvol, \
                subsection['where'], subsection['spacing'] = vals[0], vals[1], \
                                                             vals[2], vals[3:]
            nlines = int(ceil(nvol / 8.))
            vol = []
            for i in range(nlines): vol += infile.read_values('part2')
            subsection['vol'] = vol[0: nvol]
            self.meshmaker.append(('minc', subsection))

    def write_meshmaker_minc(self, section, outfile):
        outfile.write('MINC\n')
        outfile.write_values(['PART ', section['type'], '', section['dual']], 'minc')
        nvol = len(section['vol'])
        outfile.write_values([section['num_continua'], nvol,
                              section['where']] + section['spacing'], 'part1')
        nlines = int(ceil(nvol / 8.))
        for i in range(nlines):
            i1, i2 = i * 8, min((i + 1) * 8, nvol)
            vals = section['vol'][i1: i2]
            if len(vals) < 8: vals += [None] * (8 - len(vals))
            outfile.write_values(vals, 'part2')

    def read_meshfile(self, infile):
        """Reads grid from auxiliary ASCII mesh file."""
        mesh_sections = ['ELEME', 'CONNE']
        read_fn = dict(zip(mesh_sections,
                           [self.read_blocks, self.read_connections]))
        more = True
        while more:
            line = infile.readline()
            if line:
                keyword = line[0: 5].strip()
                if keyword in mesh_sections:
                    read_fn[keyword](infile)
                    self._sections.append(keyword)
            else: more = False

    def read_binary_meshfiles(self):
        """Reads grid from auxiliary binary mesh files (e.g. TOUGH2_MP
        'MESHA','MESHB' files)."""
        fa, fb = (fortran_unformatted_file(filename, 'rb') for
                  filename in self.meshfilename)
        nel, = fa.readrec('i')
        ncon, nelb = fb.readrec('2i')
        if nelb < 0:
            # used as a flag to indicate that rocktype names have been
            # replaced by indices
            nelb = -nelb
            rocktype_indices = True
        else: rocktype_indices = False
        if nel == nelb:
            # read MESHA file:
            evol, aht, pmx = (np.array(fa.readrec('%dd' % nel)) for i in range(3))
            gcoord = [np.array(fa.readrec('%dd' % nel)) for i in range(3)]
            gcoord = np.transpose(np.vstack([gc for gc in gcoord]))
            del1, del2, \
                area, beta, sig = (np.array(fa.readrec('%dd' % ncon)) for
                                   i in range(5))
            isox = np.array(fa.readrec('%di' % ncon))
            # read MESHB file:
            elem = [s.decode() for s in fb.readrec('8s' * nel)]
            if rocktype_indices:
                rtype = [self.grid.rocktypelist[i] for
                         i in np.array(fb.readrec('%di' % nel)) - 1]
            else:
                rnames = [s.decode() for s in fb.readrec('5s' * nel)]
                rtype = [self.grid.rocktype[name] for name in rnames]
            nex1, nex2 = (np.array(fb.readrec('%di' % ncon)) - 1 for i in range(2))
            # construct grid:
            self.grid.block, self.grid.blocklist = {}, []
            for i in range(nel):
                name = fix_blockname(elem[i][0: 5])
                self.grid.add_block(t2block(name, evol[i], rtype[i],
                                            centre = gcoord[i, :],
                                            ahtx = aht[i], pmx = pmx[i]))
            del evol, aht, pmx, gcoord, elem
            self.grid.connectionlist, self.grid.connection = [], {}
            for i in range(ncon):
                blk1 = self.grid.blocklist[nex1[i]]
                blk2 = self.grid.blocklist[nex2[i]]
                self.grid.add_connection(t2connection([blk1, blk2],
                                                      isox[i], [del1[i], del2[i]],
                                                      area[i], beta[i], sig[i]))
        else:
            print('Files', self.meshfilename[0], 'and', self.meshfilename[1],
                  'do not contain the same number of blocks (', nel, 'vs.', nelb, ').')
        fa.close(); fb.close()

    def write_binary_meshfiles(self):
        """Writes grid to auxiliary binary mesh files."""
        fa, fb = (fortran_unformatted_file(filename, 'wb') for filename in self.meshfilename)
        nel, ncon = self.grid.num_blocks, self.grid.num_connections
        fa.writerec('i', nel)
        fb.writerec('2i', (ncon, -nel))
        # assemble data into arrays:
        rocknames = [rt.name for rt in self.grid.rocktypelist]
        rockdict = dict([(name, i) for i, name in enumerate(rocknames)])
        block_dt = np.dtype([
            ('index', 'i4'),
            ('name', '|S8'),
            ('name8', '|S8'),
            ('rockindex', 'i4'), 
            ('volume', 'f8'),
            ('ahtx', 'f8'),
            ('pmx', 'f8'),
            ('cx', 'f8'),
            ('cy', 'f8'),
            ('cz', 'f8')])
        blkdata = np.array([(i, blk.name, blk.name.ljust(8).encode(),
                             rockdict[blk.rocktype.name],
                             blk.volume, blk.ahtx, blk.pmx, blk.centre[0],
                             blk.centre[1], blk.centre[2])
                            for i, blk in enumerate(self.grid.blocklist)], dtype = block_dt)
        for var in ['ahtx', 'pmx', 'cx', 'cy', 'cz']:
            # replace nan values with zero:
            blkdata[:][var] = np.nan_to_num(blkdata[:][var])
        blkdict = dict(zip(blkdata[:]['name'], blkdata[:]['index']))
        con_dt = np.dtype([
            ('blk1', '|S8'),
            ('blk2', '|S8'),
            ('blk1index', 'i4'),
            ('blk2index', 'i4'),
            ('d1', 'f8'),
            ('d2', 'f8'),
            ('dirn', 'i4'),
            ('area', 'f8'),
            ('dircos', 'f8'),
            ('sigma', 'f8')])
        condata = np.array([(con.block[0].name.ljust(8).encode(),
                             con.block[1].name.ljust(8).encode(),
                             blkdict[con.block[0].name.encode()],
                             blkdict[con.block[1].name.encode()],
                             con.distance[0], con.distance[1],
                             con.direction, con.area, con.dircos, con.sigma)
                            for con in self.grid.connectionlist], dtype = con_dt)
        condata[:]['sigma'] = np.nan_to_num(condata[:]['sigma'])
        # write MESHA file:
        for var in ['volume', 'ahtx', 'pmx', 'cx', 'cy', 'cz']:
            fa.writerec('%dd' % nel, blkdata[:][var])
        for var in ['d1', 'd2', 'area', 'dircos', 'sigma']:
            fa.writerec('%dd' % ncon, condata[:][var])
        fa.writerec('%di' % ncon, condata[:]['dirn'])
        for var in ['blk1', 'blk2']:
            fa.writerec('8s' * ncon, condata[:][var])
        # write MESHB file:
        fb.writerec('8s' * nel, blkdata[:]['name8'])
        fb.writerec('%di' % nel, blkdata[:]['rockindex'] + 1)
        for var in ['blk1index', 'blk2index']:
            fb.writerec('%di' % ncon, condata[:][var] + 1)
        fa.close()
        fb.close()

    def read_extra_precision(self):
        """Reads extra precision data from auxiliary file, with same base name
        as the data file name, but with file extension .pdat.
        """
        from os.path import exists
        filename = self.extra_precision_filename
        if exists(filename):
            mode = 'r' if sys.version_info > (3,) else 'rU'
            xpfile = t2_extra_precision_data_parser(self.extra_precision_filename, mode,
                                                    read_function = self.read_function)
            read_fn = dict(zip(t2_extra_precision_sections,
                               [self.read_rocktypes, self.read_blocks, self.read_connections,
                                self.read_rpcap, self.read_generators]))
            more = True
            while more:
                line = xpfile.readline()
                if line:
                    keyword = line[0: 5].strip()
                    if keyword in ['ENDCY', 'ENDFI']: more = False
                    elif keyword in read_fn:
                        self.extra_precision.append(keyword)
                        read_fn[keyword](xpfile)
                else: more = False
            xpfile.close()
            if self.extra_precision:
                self.update_read_write_functions()
                self.echo_extra_precision = any([section in self._sections for
                                                 section in self.extra_precision])
            else: self.echo_extra_precision = False

    def write_extra_precision(self, extra_precision = None, echo_extra_precision = None):
        """Writes AUTOUGH2 extra precision data to auxiliary file."""
        if extra_precision is not None:
            self.extra_precision = extra_precision
        if echo_extra_precision is not None:
            self.echo_extra_precision = echo_extra_precision
        if self.extra_precision:
            xpfile = t2_extra_precision_data_parser(self.extra_precision_filename, 'w')
            write_fn = dict(zip(t2_extra_precision_sections,
                                [self.write_rocktypes, self.write_blocks, self.write_connections,
                                 self.write_rpcap, self.write_generators]))
            for section in self.extra_precision:
                write_fn[section](xpfile)
                if section in self._sections and not self.echo_extra_precision:
                    self._sections.remove(section)
            xpfile.close()

    def read(self, filename = '', meshfilename = ''):
        """Reads data from file.  Mesh data can optionally be read from an
        auxiliary file.  Extra precision data will also be read from
        an associated '.pdat' file, if it exists.
        """
        if filename: self.filename = filename
        mode = 'r' if sys.version_info > (3,) else 'rU'
        infile = t2data_parser(self.filename, mode, read_function = self.read_function)
        self.read_title(infile)
        self._sections = []
        more = True
        next_line = None
        while more:
            if next_line: line = next_line
            else: line = infile.readline()
            if line:
                keyword = line[0: 5].strip()
                if keyword in ['ENDCY', 'ENDFI']:
                    more = False
                    self.end_keyword = keyword
                elif keyword in t2data_sections:
                    fn = self.read_fn[keyword]
                    next_line = None
                    if keyword == 'SHORT': fn(infile, line)
                    elif keyword == 'PARAM': next_line = fn(infile)
                    else: fn(infile)
                    self._sections.append(keyword)
            else: more = False
        infile.close()
        if meshfilename and (self.grid.num_blocks == 0):
            self.meshfilename = meshfilename
            if isinstance(meshfilename, str):
                mode = 'r' if sys.version_info > (3,) else 'rU'
                meshfile = t2data_parser(self.meshfilename, mode, read_function = self.read_function)
                self.read_meshfile(meshfile)
                meshfile.close()
            elif isinstance(meshfilename, (list, tuple)):
                if len(meshfilename) == 2: self.read_binary_meshfiles()
            else: print('Mesh filename must be either a string or a two-element tuple or list.')
        return self

    def write(self, filename = '', meshfilename = '',
              extra_precision = None, echo_extra_precision = None):
        """Writes data to file.  Mesh data can optionally be written to an
        auxiliary file.  For AUTOUGH2, if extra_precision is True or a
        list of section names, the corresponding data sections will be
        written to an extra precision file; otherwise, the same
        sections (if any) that were read in as extra precision will
        also be written out as extra precision.  If
        echo_extra_precision is True, the extra precision sections
        will also be written to the main data file.
        """
        if filename: self.filename = filename
        if self.filename =='': self.filename = 't2data.dat'
        self.update_sections()
        mesh_sections = []
        if meshfilename: self.meshfilename = meshfilename
        if self.meshfilename:
            if isinstance(self.meshfilename, str):
                meshfile = t2data_parser(self.meshfilename, 'w')
                self.write_blocks(meshfile)
                self.write_connections(meshfile)
                meshfile.close()
                mesh_sections = ['ELEME', 'CONNE']
            elif isinstance(self.meshfilename, (list, tuple)):
                if len(self.meshfilename) == 2:
                    self.write_binary_meshfiles()
                    mesh_sections = ['ELEME', 'CONNE']
        if self.type == 'AUTOUGH2':
            self.write_extra_precision(extra_precision,  echo_extra_precision)
        outfile = t2data_parser(self.filename, 'w')
        self.write_title(outfile)
        for keyword in self._sections:
            if (keyword not in mesh_sections) and \
                    ((keyword not in self.extra_precision) or
                     (keyword in self.extra_precision and self.echo_extra_precision)):
                    self.write_fn[keyword](outfile)
        outfile.write(self.end_keyword + '\n')
        outfile.close()

    def transfer_rocktypes_from(self, source, mapping):
        """Transfers rock types (definitions and assignments) from another
        t2data object, using the specified block mapping."""
        from copy import deepcopy
        self.grid.rocktypelist = deepcopy(source.grid.rocktypelist)
        self.grid.rocktype = deepcopy(source.grid.rocktype)
        for blk in self.grid.blocklist:
            blk.rocktype = self.grid.rocktype[source.grid.block[mapping[blk.name]].rocktype.name]

    def transfer_generators_from(self, source, sourcegeo, geo,
                                 top_generator = [], bottom_generator = [],
                                 mapping = {}, colmapping = {},
                                 rename = False, preserve_totals = False):
        """Transfers generators from another t2data object, using the
        specified top and bottom generator lists and optional block
        and column mappings.  If the rename parameter is False,
        generators other than those at the top or bottom of the model
        will keep their original names; if True, these generators will
        be renamed according to their new column names.  If
        preserve_totals is True, the transfer will attempt to preserve
        total generation, when a source block is mapped into a set of
        blocks with a different total volume.
        """
        from copy import deepcopy
        tablegens = [' AIR', 'COM1', 'COM2', 'COM3', 'COM4', 'COM5',
                     'HEAT', 'MASS', 'NACL', 'TRAC', ' VOL']
        if (colmapping == {}) or (mapping == {}):
            mapping, colmapping = sourcegeo.block_mapping(geo, True)
        bbox = sourcegeo.bounds
        qt = quadtree(bbox, sourcegeo.columnlist)
        incols = [col for col in geo.columnlist if
                  sourcegeo.column_containing_point(col.centre, qtree = qt) is not None]
        self.clear_generators()
        col_generator = top_generator + bottom_generator
        for sourcegen in source.generatorlist:
            sourcecategory = sourcegeo.layer_name(sourcegen.name)
            sourcecolname = sourcegeo.column_name(sourcegen.block)
            if sourcecategory in col_generator:
                mappedcols = [col for col in incols if
                              colmapping[col.name] == sourcecolname]
                if preserve_totals: area = sum([col.area for col in mappedcols])
                else: area = sourcegeo.column[sourcecolname].area
                for col in mappedcols:
                    gen = deepcopy(sourcegen)
                    area_ratio = col.area / area
                    if gen.ltab: ntimes = abs(gen.ltab)
                    else: ntimes = 1
                    if gen.type in tablegens:
                        if gen.gx: gen.gx *= area_ratio
                        if ntimes > 1: gen.rate = [rate * area_ratio for rate in gen.rate]
                    if geo.convention == sourcegeo.convention:
                        category = sourcecategory
                    else:
                        category = ['%2d' % col_generator.index(sourcecategory),
                                    sourcecategory, sourcecategory][geo.convention]
                    gen.name = geo.block_name(category, col.name)
                    if sourcecategory in top_generator:
                        layername = geo.layerlist[geo.num_layers - col.num_layers].name
                    elif sourcecategory in bottom_generator: layername = geo.layerlist[-1].name
                    gen.block = geo.block_name(layername, col.name)
                    self.add_generator(gen)
            else: # other generators, do block by block:
                sourceblock = source.grid.block[sourcegen.block]
                mappedblocks = [blk for blk in self.grid.blocklist if
                                mapping[blk.name] == sourceblock.name]
                if preserve_totals: vol = sum([blk.volume for blk in mappedblocks])
                else: vol = sourceblock.volume
                for blk in mappedblocks:
                    gen = deepcopy(sourcegen)
                    if gen.ltab: ntimes = abs(gen.ltab)
                    else: ntimes = 1
                    vol_ratio = blk.volume / vol
                    if gen.type in tablegens:
                        if gen.gx: gen.gx *= vol_ratio
                        if ntimes > 1: gen.rate = [rate * vol_ratio for rate in gen.rate]
                    if rename:
                        if geo.convention == sourcegeo.convention: category = sourcecategory
                        else: category = [' 0', sourcecategory, sourcecategory][geo.convention]
                        colname = geo.column_name(blk.name)
                        gen.name = geo.block_name(category, colname)
                    gen.block = blk.name
                    self.add_generator(gen)

    def transfer_from(self, source, sourcegeo, geo,
                      top_generator = [], bottom_generator = [],
                      sourceinconfilename = '', inconfilename = '',
                      rename_generators = False, preserve_generation_totals = False):
        """Copies parameters, rock types and assignments, generators and
        initial conditions from another t2data object (without
        altering the grid structure).  The top_generator and
        bottom_generator lists specify the generators which are to be
        kept at the top or bottom of the model, respectively.  They
        contain the 'layer' part of the generator name for top and
        bottom generators.  If both the inconfilename parameters are
        specified, a new initial conditions file with filename
        'inconfilename' is written to disk, with initial conditions
        transferred from the file 'sourceinconfilename'.  If the
        rename_generators parameter is False, generators (other than
        those at the top and bottom of the model) keep their original
        names- otherwise, they are renamed according to their new
        column names.  If preserve_generation_totals is True,
        generators are transferred in such a way as to preserve total
        generation, even when a source block is mapped into a set of
        blocks with a different total volume.  This can however alter
        the distribution of specific generation.
        """
        mapping, colmapping = sourcegeo.block_mapping(geo, True)
        from copy import copy, deepcopy
        self.grid = t2grid().fromgeo(geo)
        self.simulator = source.simulator
        self.parameter = deepcopy(source.parameter)
        if self.parameter['print_block'] is not None:
            mappedblocks = [blk for blk in self.grid.blocklist if
                            mapping[blk.name] == self.parameter['print_block']]
            if len(mappedblocks) > 0: self.parameter['print_block'] = mappedblocks[0].name
            else: self.parameter['print_block'] = None
        self.multi = copy(source.multi)
        self.start = source.start
        self.noversion = source.noversion
        self.relative_permeability = copy(source.relative_permeability)
        self.capillarity = copy(source.capillarity)
        self.lineq = copy(source.lineq)
        self.solver = copy(source.solver)
        self.diffusion = copy(source.diffusion)
        self.selection = copy(source.selection)
        self.output_times = copy(source.output_times)
        self.transfer_rocktypes_from(source, mapping)
        self.transfer_generators_from(source, sourcegeo, geo,
                                      top_generator, bottom_generator,
                                      mapping, colmapping,
                                      rename_generators, preserve_generation_totals)
        # short output (these can't really be transferred sensibly):
        self.short_output = {}
        self.history_block = {}
        self.history_connection = {}
        self.history_generator = {}
        # incons (within data file):
        self.incon = {}
        for blkname, inc in source.incon.items():
            mappedblocks = [blk for blk in self.grid.blocklist if
                            mapping[blk.name] == blkname]
            for blk in mappedblocks: self.incon[blk.name] = inc
        self.indom = copy(source.indom)
        # incon file:
        if (sourceinconfilename != '') and (inconfilename != ''):
            sourceinc = t2incon(sourceinconfilename)
            inc = t2incon()
            inc.transfer_from(sourceinc, sourcegeo, geo, mapping, colmapping)
            inc.write(inconfilename)

    def convert_mulkom_heat_conductivity(self):
        """Converts MULKOM-style rock heat conductivities to TOUGH2-style.
        MULKOM used a formulation based on wet conductivity and liquid
        conductivity, weighted by porosity (though the liquid value
        appears to have always been zero).  TOUGH2 uses a formulation
        based on wet and dry conductivities, weighted by liquid
        saturation (though by default the wet and dry values are the
        same).  Here we scale the rock heat conductivities to give the
        same effective values as the MULKOM formulation.
        """
        for rt in self.grid.rocktypelist:
            rt.conductivity *= (1. - rt.porosity)

    def convert_AUTOUGH2_parameters_to_TOUGH2(self, warn = True, MP = False):
        """Converts AUTOUGH2 parameters to TOUGH2 parameters, with optional
        warnings about options that aren't supported in TOUGH2.  If MP
        is True, convert to a file suitable for TOUGH2_MP.
        """
        # Modify MULTI:
        if self.multi:
            if 'eos' in self.multi: del self.multi['eos']
            self.multi['num_inc'] = None
        # Convert LINEQ into corresponding MOP(21) option:
        if self.lineq:
            if self.lineq['type'] <= 1: solver_type = 4
            else: solver_type = 5
        else: solver_type = 4
        self.lineq = {}
        self.delete_section('LINEQ')
        # Convert MOPs:
        warnings = []
        if self.parameter['option'][10] == 2:
            self.parameter['option'][10] = 0
            self.convert_mulkom_heat_conductivity()
            warnings.append('MOP(10)=2: MULKOM rock heat conductivities' + \
                            ' (values have been converted to TOUGH2 equivalents)')
        if self.parameter['option'][12] == 2:
            self.parameter['option'][12] = 0
            warnings.append('MOP(12)=2: piecewise linear well table interpolation')
        self.parameter['option'][21] = solver_type
        if self.parameter['option'][22] > 0:
            self.parameter['option'][22] = 0
            warnings.append('MOP(22)>0: USERBC')
        if self.parameter['option'][23] > 0:
            isat2 = self.simulator.startswith('AUTOUGH2') and \
                    (not self.simulator.startswith('AUTOUGH2.2'))
            ismulkom = self.simulator.startswith('MULKOM')
            mulkom_compatibility = self.parameter['option'][23] in [0, 1]
            if (isat2 or ismulkom) and mulkom_compatibility:
                self.convert_mulkom_heat_conductivity()
                warnings.append('MOP(23)>0: MULKOM/TOUGH2 backward compatibility')
            self.parameter['option'][23] = 0
        if self.parameter['option'][24] > 0:
            self.parameter['option'][24] = 0
            warnings.append('MOP(24)>0: initial printout of tables')
        if MP: # these MOPs mean different things in TOUGH2_MP:
            if self.parameter['option'][14] > 0:
                self.parameter['option'][14] = 0
                warnings.append('MOP(14)>0: Pivot failure handling')
            if self.parameter['option'][17] > 0:
                self.parameter['option'][17] = 0
                warnings.append('MOP(17)>0: Jacobian scaling')
            if self.parameter['option'][20] > 0:
                self.parameter['option'][20] = 0
                warnings.append('MOP(20)>0: Disabling vapour pressure lowering')
            self.parameter['option'][21] = 0
        if warn and len(warnings) > 0:
            print('The following options are not supported in TOUGH2:')
            for warning in warnings: print(warning)

    def convert_TOUGH2_parameters_to_AUTOUGH2(self, warn = True, MP = False):
        """Converts TOUGH2 parameters to AUTOUGH2 parameters, with optional
        warnings about options that aren't supported in AUTOUGH2.  If
        MP is True, treat the file as a TOUGH2_MP data file.
        """
        # modify MULTI:
        if self.multi: self.multi['num_inc'] = None
        # set up LINEQ:
        if MP: solver_type = 2
        else:
            if 'type' in self.solver: solver_type = self.solver['type']
            else: solver_type = self.parameter['option'][21]
        self.lineq = {'type': [2, 1, 2, 2, 1, 2, 1][solver_type], 'epsilon': None,
                      'max_iterations': None, 'gauss': None, 'num_orthog': None}
        self.insert_section('LINEQ')
        self.solver = {}
        # Convert MOPs:
        warnings = []
        if self.parameter['option'][12] == 2:
            self.parameter['option'][12] = 0
            warnings.append('MOP(12)=2: rigorous step rate well table interpolation')
        if self.parameter['option'][22] > 0:
            warnings.append('MOP(22)>0: dispersion module T2DM')
        self.parameter['option'][22] = 0
        if self.parameter['option'][23] > 0:
            warnings.append('MOP(23)>0: dispersion module T2DM')
        self.parameter['option'][23] = 0
        if self.parameter['option'][24] > 0:
            warnings.append('MOP(24)>0:  handling of multiphase diffusive fluxes at interfaces')
        self.parameter['option'][24] = 0
        if MP: # these MOPs mean different things in TOUGH2_MP:
            if self.parameter['option'][14] > 0:
                self.parameter['option'][14] = 0
                warnings.append('MOP(14)>0: 8-character element names')
            if self.parameter['option'][17] > 0:
                self.parameter['option'][17] = 0
                warnings.append('MOP(17)>0:  generation of a flow9.dat file for T2R3D')
            if self.parameter['option'][20] > 0:
                self.parameter['option'][20] = 0
                warnings.append('MOP(20)>0: long format for CONNE and GENER indices')
            if self.parameter['option'][21] > 0:
                warnings.append('MOP(21)>0: perform extra Newton iteration after convergence')
        self.parameter['option'][21] = 0 # not used in AUTOUGH2
        if warn and len(warnings) > 0:
            print('The following options are not supported in AUTOUGH2:')
            for warning in warnings: print(warning)

    def convert_AUTOUGH2_generators_to_TOUGH2(self, warn = True):
        """Convert AUTOUGH2 generators to TOUGH2 generators, with optional
        warnings about generator types that aren't supported in TOUGH2
        (and will be deleted).
        """
        allowed = ['HEAT', 'WATE', 'AIR ', 'MASS', 'DELV']
        convert = {'CO2 ':'COM2'}
        delgens = []
        for gen in self.generatorlist:
            if gen.type in convert: gen.type = convert[gen.type]
            elif not ((gen.type in allowed) or gen.type.startswith('COM')):
                delgens.append((gen.block, gen.name))
        if warn and len(delgens) > 0:
            print('The following generators have types not supported' + \
                  ' by TOUGH2 and have been deleted:')
            print(delgens)

    def convert_short_to_history(self):
        """Converts AUTOUGH2 SHORT output to TOUGH2 history (FOFT, COFT, GOFT)."""
        if 'block' in self.short_output:
            self.history_block = self.short_output['block'][:]
        if 'connection' in self.short_output:
            self.history_connection = self.short_output['connection'][:]
        if 'generator' in self.short_output:
            self.history_generator = self.short_output['generator'][:]
        self.short_output = {}

    def convert_history_to_short(self):
        """Converts TOUGH2 history (FOFT, COFT, GOFT) to AUTOUGH2 SHORT
        output.  Items referring to blocks or connections not present
        in the grid are discarded.
        """
        self.short_output = {}
        if self.history_block:
            blks = [blk for blk in self.history_block if isinstance(blk, t2block)]
            if blks: self.short_output['block'] = blks
        if self.history_connection:
            cons = [con for con in self.history_connection if isinstance(con, t2connection)]
            if cons: self.short_output['connection'] = cons
        if self.history_generator:
            gens = [gen for gen in self.history_generator if isinstance(gen, t2generator)]
            if gens: self.short_output['generator'] = gens
        self.history_block = []
        self.history_connection = []
        self.history_generator = []

    def convert_to_TOUGH2(self, warn = True, MP = False):
        """Converts an AUTOUGH2 data file to a TOUGH2 data file.  Various MOP
        parameters are changed to try to make them the TOUGH2
        simulation give similar results to AUTOUGH2 where possible.
        AUTOUGH2-specific input blocks are removed and some generator
        types changed.  If MP is True, the conversion is done to a
        TOUGH2_MP data file, which treats a few of the parameters
        differently.
        """
        if MP: self.filename = 'INFILE'
        self.simulator = ''
        self.delete_section('SIMUL')
        self.convert_AUTOUGH2_parameters_to_TOUGH2(warn, MP)
        self.convert_AUTOUGH2_generators_to_TOUGH2(warn)
        self.convert_short_to_history()

    def convert_to_AUTOUGH2(self, warn = True, MP = False,
                            simulator = 'AUTOUGH2.2', eos = 'EW'):
        """Converts a TOUGH2 data file to an AUTOUGH2 data file."""
        if self.filename:
            if not self.filename.lower().endswith('.dat'):
                if self.filename[0].isupper(): self.filename += '.DAT'
                else: self.filename += '.dat'
        self.simulator = simulator.ljust(10) + eos
        self.insert_section('SIMUL')
        if self.multi: self.multi['eos'] = eos
        self.convert_TOUGH2_parameters_to_AUTOUGH2(warn, MP)
        self.convert_history_to_short()

    def rename_blocks(self, blockmap = {}, invert = False, fix_blocknames = True):
        """Rename blocks in TOUGH2 data file according to specified block
        mapping. The mapping is applied to block names in the grid,
        initial conditions and generators, print block and history
        file specifications. If invert is True, the inverse of the
        specified block mapping is applied.
        """

        if invert: blockmap = {v:k for k,v in blockmap.items()}
        if fix_blocknames: fix_block_mapping(blockmap)

        self.grid.rename_blocks(blockmap, fix_blocknames = False)

        if self.incon:
            for k,v in blockmap.items():
                if k in self.incon:
                    inc = self.incon[k]
                    del self.incon[k]
                    self.incon[v] = inc

        for gen in self.generatorlist:
            if gen.block in blockmap:
                keys = (gen.block, gen.name)
                del self.generator[keys]
                gen.block = blockmap[gen.block]
                newkeys = (gen.block, gen.name)
                self.generator[newkeys] = gen

        if self.parameter['print_block'] in blockmap:
            self.parameter['print_block'] = blockmap[self.parameter['print_block']]

        for i, blk in enumerate(self.history_block):
            if not isinstance(blk, t2block):
                if blk in blockmap:
                    self.history_block[i] = blockmap[blk]

        for i, con in enumerate(self.history_connection):
            if not isinstance(con, t2connection):
                if any([name in blockmap for name in con]):
                    mapped_con = tuple([blockmap[name] if name in blockmap else name
                                          for name in con])
                    self.history_connection[i] = mapped_con

        for i, gen in enumerate(self.history_generator):
            if not isinstance(gen, t2block):
                if gen in blockmap:
                    self.history_generator[i] = blockmap[gen]

    def effective_incons(self, incons = None):
        """Returns effective initial conditions, based on combination of specified incons,
        default initial conditions and INCON and INDOM data blocks in the t2data object.
        If default initial conditions from PARAM are used everywhere, the function returns
        a list with those values. Otherwise, a t2incon object is returned."""

        default_incs = self.parameter['default_incons'][:]
        default_incs = trim_trailing_nones(default_incs)
        effective_incs = default_incs

        if self.indom or self.incon or incons:
            effective_incs = self.grid.incons(default_incs)
            if self.indom:
                for blk in self.grid.blocklist:
                    if blk.rocktype.name in self.indom:
                        effective_incs[blk.name] = self.indom[blk.rocktype.name]
            if self.incon:
                for blkname in self.incon:
                    effective_incs[blkname] = self.incon[blkname][1]
            if isinstance(incons, t2incon):
                for blkinc in incons:
                    effective_incs[blkinc.block] = blkinc.variable

        return effective_incs

    def mesh_json(self, geo, mesh_filename):
        """Converts mesh data to Waiwera JSON dictionary, including overridden
        face permeability directions."""
        jsondata = {}
        jsondata['mesh'] = {'filename': mesh_filename}
        tol = 1.e-6
        if abs(geo.permeability_angle) > tol:
            jsondata['mesh']['permeability_angle'] = geo.permeability_angle
        face_directions = []
        from math import cos, sin, radians
        anglerad = radians(geo.permeability_angle)
        c, s = cos(anglerad), sin(anglerad)
        rotation = np.array([[c, s], [-s, c]])
        for blknames in geo.block_connection_name_list:
            if blknames in self.grid.connection:
                names = blknames
                con = self.grid.connection[blknames]
            else:
                rnames = blknames[::-1]
                if rnames in self.grid.connection:
                    names = rnames
                    con = self.grid.connection[rnames]
                else:
                    raise Exception ('Connection not found: ' + str(blknames))
            blkindices = [geo.block_name_index[blkname] -
                          geo.num_atmosphere_blocks for blkname in names]
            laynames = [geo.layer_name(blkname) for blkname in names]
            if laynames[0] != laynames[1]: # vertical connection
                underground = all([blkindex >= 0 for blkindex in blkindices])
                if underground and con.direction != 3:
                    face_directions.append({
                        "cells": blkindices,
                        "permeability_direction": con.direction})
            else:
                colnames = [geo.column_name(blkname) for blkname in names]
                d = geo.column[colnames[1]].centre - geo.column[colnames[0]].centre
                d2 = np.dot(rotation, d)
                expected_direction = np.argmax(abs(d2)) + 1
                if con.direction != expected_direction:
                    face_directions.append({
                        "cells": blkindices,
                        "permeability_direction": con.direction})
        if face_directions: jsondata['mesh']['faces'] = face_directions
        return jsondata

    def eos_json(self, eos):
        """Converts TOUGH2 EOS data to Waiwera JSON dictionary. Also returns
        a second dictionary with tracer data."""
        jsondata = {}
        supported_eos = {'W': 'w', 'EW': 'we', 'EWC': 'wce', 'EWAV': 'wae',
                         'EWT': 'we', 'EWTD': 'we'}
        aut2eosname = ''
        if eos is None:
            if self.multi:
                if 'eos' in self.multi:
                    if self.multi['eos']: aut2eosname = self.multi['eos'].strip()
            elif self.simulator:
                for eosname in supported_eos.keys():
                    if self.simulator.endswith(eosname):
                        autseosname = eosname
        else:
            if isinstance(eos, int):
                eos_from_index = {1: 'EW', 2: 'EWC', 3: 'EWA', 4: 'EWAV'}
                if eos in eos_from_index: aut2eosname = eos_from_index[eos]
            else: aut2eosname = eos
        if aut2eosname:
            if aut2eosname in supported_eos:
                jsondata['eos'] = {'name': supported_eos[aut2eosname]}
                if jsondata['eos']['name'] == 'w':
                    jsondata['eos']['temperature'] = self.parameter['default_incons'][1]
            else:
                raise Exception ('EOS not supported:' + aut2eosname)
        else:
            raise Exception ('EOS not detected.')
        if aut2eosname in ['EWT', 'EWTD']:
            tracerdata = {'tracer': {'name': 'tracer', 'phase': 'liquid'}}
            if aut2eosname == 'EWTD':
                diffusion = np.array(self.diffusion)
                if np.all(diffusion < 0) and np.allclose(diffusion, diffusion[0][0]):
                    D = -diffusion[0][0]
                    tracerdata['tracer']['diffusion'] = D
                else:
                    raise Exception ('Unhandled diffusion type: %s' % str(self.diffusion))
        else:
            tracerdata = None
        return jsondata, tracerdata

    def timestepping_json(self):
        """Converts TOUGH2 timestepping/ iteration parameters to Waiwera JSON
        dictionary."""
        jsondata = {}
        tstop = self.parameter['tstop']
        if tstop == 0.0: tstop = None
        jsondata['time'] = {'start': self.parameter['tstart'],
                            'stop': tstop}
        maxit = self.parameter['max_iterations']
        if maxit is None or maxit == 0: maxit = 8
        abstol = self.parameter['absolute_error']
        if abstol == 0: abstol = 1.0
        reltol = self.parameter['relative_error']
        if reltol == 0.: reltol = 1.e-5
        jsondata['time']['step'] = \
            {'maximum': {'size': self.parameter['max_timestep']},
             'method': 'beuler',
             'solver': {'nonlinear': {'tolerance': {'function':
                                          {'absolute': abstol, 'relative': reltol}},
                                      'maximum': {'iterations': maxit}}}}
        if self.parameter['max_timesteps'] is None:
            jsondata['time']['step']['maximum']['number'] = 0
        elif self.parameter['max_timesteps'] < 0:
            jsondata['time']['step']['maximum']['number'] = None
        else:
            jsondata['time']['step']['maximum']['number'] = self.parameter['max_timesteps']
        if self.parameter['const_timestep'] < 0. :
            jsondata['time']['step'].update({'size': self.parameter['timestep'],
                                             'adapt': {'on': False}})
        else:
            jsondata['time']['step'].update({'size': self.parameter['const_timestep'],
                                             'adapt': {'on': True}})
        if self.parameter['option'][16] > 0:
            redlt = self.parameter['timestep_reduction']
            if redlt is None or redlt == 0:
                redlt = 5 if self.type == 'AUTOUGH2' else 4
            jsondata['time']['step']['adapt'].update(
                {'method': 'iteration',
                 'reduction': 1. / redlt,
                 'amplification': 2.,
                 'minimum': float(self.parameter['option'][16]), 'maximum': float(maxit)})
        return jsondata

    def rocks_json(self, geo, atmos_volume, mesh_coords):
        """Converts TOUGH2 rocktype definition and assignment data to Waiwera JSON
        dictionary."""
        jsondata = {}
        jsondata['rock'] = {'types': []}
        ir, rock_index = 0, {}
        if mesh_coords == 'xyz': perm_size = 3
        else: perm_size = 2
        for rt in self.grid.rocktypelist:
            rtdata = {'name': rt.name, 'density': rt.density, 'porosity': rt.porosity,
                      'permeability': list(rt.permeability[:perm_size]),
                      'wet_conductivity': rt.conductivity, 'specific_heat': rt.specific_heat}
            dry_cond = rt.dry_conductivity
            if dry_cond is not None and dry_cond > 0.0: rtdata['dry_conductivity'] = dry_cond
            else: rtdata['dry_conductivity'] = rt.conductivity
            rtdata['cells'] = []
            jsondata['rock']['types'].append(rtdata)
            rock_index[rt.name] = ir
            ir += 1
        for blkname in geo.block_name_list:
            blk = self.grid.block[blkname]
            rockname = blk.rocktype.name
            blk_index = geo.block_name_index[blk.name] - geo.num_atmosphere_blocks
            if 0. < blk.volume < atmos_volume:
                jsondata['rock']['types'][rock_index[rockname]]['cells'].append(blk_index)
        return jsondata

    def relative_permeability_json(self):
        """Converts TOUGH2 relative permeability data to Waiwera JSON dictionary."""
        jsondata = {}
        stol = 1.e-9
        if self.relative_permeability:
            rp = {}
            rp_types = {1: 'linear', 2: 'pickens', 3: 'corey', 4: 'grant',
                        5: 'fully mobile', 7: 'van Genuchten'}
            itype = self.relative_permeability['type']
            pars = self.relative_permeability['parameters']
            if itype in rp_types:
                rp['type'] = rp_types[itype]
                if itype == 1:
                    rp['liquid'] = [pars[0], pars[2]]
                    rp['vapour'] = [pars[1], pars[3]]
                elif itype == 2:
                    rp['power'] = pars[0]
                elif itype in [3, 4]:
                    rp['slr'] = pars[0]
                    rp['ssr'] = pars[1]
                elif itype == 7:
                    rp['lambda'] = pars[0]
                    rp['slr'] = pars[1]
                    rp['sls'] = pars[2]
                    if pars[3] > stol:
                        rp['sum_unity'] = False
                        rp['ssr'] = pars[3]
                    else: rp['sum_unity'] = True
                jsondata['relative_permeability'] = rp
            elif self.type == 'AUTOUGH2' and itype == 19:
                # tri-linear: convert to table
                rp['type'] = 'table'
                rp['liquid'] = [[0, 0], [pars[0], pars[4]],
                                [pars[2], pars[6]], [1, 1]]
                rp['vapour'] = [[0, 0], [pars[1], 0],
                                [pars[3], pars[5]], [1,1]]
                jsondata['relative_permeability'] = rp
            else:
                raise Exception ('Unhandled relative permeability type: %d' % itype)
        else: jsondata['relative_permeability'] = {'type': 'fully mobile'}
        return jsondata

    def capillary_pressure_json(self):
        """Converts TOUGH2 capillary pressure data to Waiwera JSON dictionary."""
        jsondata = {}
        stol = 1.e-9
        if self.capillarity:
            cp = {}
            cp_types = {1: 'linear', 7: 'van Genuchten', 8: 'zero'}
            itype = self.capillarity['type']
            pars = self.capillarity['parameters']
            if itype in cp_types:
                cp['type'] = cp_types[itype]
                if itype == 1:
                    cp['pressure'] = pars[0]
                    cp['saturation_limits'] = [pars[1], pars[2]]
                elif itype == 7:
                    cp['lambda'] = pars[0]
                    cp['slr'] = pars[1]
                    cp['P0'] = 1. / pars[2]
                    if pars[3] > stol: cp['Pmax'] = pars[3]
                    cp['sls'] = pars[4]
                elif itype == 8:
                    cp = None
                jsondata['capillary_pressure'] = cp
            else:
                raise Exception ('Unhandled capillary pressure type: %d' % itype)
        else: jsondata['capillary_pressure'] = None
        return jsondata

    def initial_json(self, geo, incons, eos, tracer = None):
        """Converts initial condition specifications to Waiwera JSON dictionary."""
        jsondata = {}

        if isinstance(incons, str):
            jsondata['initial'] = {'filename': incons}
        elif isinstance(incons, list):
            num_primary = waiwera_eos_num_primary[eos]
            if len(incons) >= num_primary:
                jsondata['initial'] = {'primary': incons[:num_primary]}
            if tracer and len(incons) >= num_primary + 1:
                jsondata['initial']['tracer'] = incons[num_primary]
            if incons:
                if eos in primary_to_region_funcs:
                    primary_to_region = primary_to_region_funcs[eos]
                    jsondata['initial']['region'] = primary_to_region(incons)
                else:
                    raise Exception("Finding thermodynamic region from primary variables not yet supported for EOS:" + eos)
        elif isinstance(incons, t2incon):
            num_primary = waiwera_eos_num_primary[eos]
            if eos in primary_to_region_funcs:
                jsondata['initial'] = {'primary': [], 'region': []}
                if tracer: jsondata['initial']['tracer'] = []
                primary_to_region = primary_to_region_funcs[eos]
                for blkname in geo.block_name_list[geo.num_atmosphere_blocks:]:
                    primary = incons[blkname].variable
                    jsondata['initial']['primary'].append(primary[:num_primary])
                    jsondata['initial']['region'].append(primary_to_region(primary))
                    if tracer: jsondata['initial']['tracer'].append(primary[num_primary])
                if np.isclose(jsondata['initial']['primary'],
                              jsondata['initial']['primary'][0], rtol = 1.e-8).all():
                    jsondata['initial']['primary'] = jsondata['initial']['primary'][0]
                if len(set(jsondata['initial']['region'])) == 1:
                    jsondata['initial']['region'] = jsondata['initial']['region'][0]
                if tracer:
                    if np.isclose(jsondata['initial']['tracer'],
                                  jsondata['initial']['tracer'][0], rtol = 1.e-8).all():
                        jsondata['initial']['tracer'] = jsondata['initial']['tracer'][0]
                    else:
                        raise Exception("Inhomogeneous tracer initial conditions not yet supported.")
            else:
                raise Exception("Finding thermodynamic region from primary variables not yet supported for EOS:" + eos)
        return jsondata

    def generators_json(self, geo, eosname, tracer = None):
        """Converts TOUGH2 generator data to Waiwera JSON dictionary, containing
        data for sources and the source network."""

        jsondata = {}
        eos_num_equations = {'w': 1, 'we': 2, 'wce': 3, 'wae': 3}
        num_eqns = eos_num_equations[eosname]
        unsupported_types = {'CO2 ', 'FEED', 'HLOS', 'MAKE', 'POWR',
                             'TOST', 'VOL.', 'WBRE', 'WFLO', 'XIN2'}
        limit_type = {'DELG': 'steam', 'DMAK': 'steam', 'DELS': 'steam',
                      'DELT': 'total', 'DELW': 'water', 'DMAT': 'total'}
        reinjection_contributors = {'DELG', 'DELS', 'DELT', 'DELW', 'DELV',
                                    'DMAK', 'DMAT'}
        used_names = {}
        # prepend block names to generator names if generator names are not unique:
        use_block_names = len(self.generator) < self.num_generators

        def unique_name(gen):
            """Generates a unique generator name not already in used_names."""
            if gen.name == '': return gen.name
            else:
                if use_block_names:
                    name = '%5s%5s' % (gen.block, gen.name)
                else:
                    name = gen.name
                if name in used_names:
                    new_name = '%s_%d' % (name, used_names[name])
                    used_names[name] += 1
                else:
                    new_name = name
                    used_names[name] = 1
                return new_name

        def separator(P):
            if P is None: Psep = 0.55e6
            else:
                if P > 0.: Psep = P
                elif P < 0: Psep = [1.45e6, 0.55e6]
                else: Psep = 0.55e6
            return {'pressure': Psep}

        if self.parameter['option'][12] == 0:
            interp_type, averaging_type = "linear", "endpoint"
        elif self.parameter['option'][12] == 1:
            interp_type, averaging_type = "step", "endpoint"
        else:
            interp_type, averaging_type = "linear", "integrate"

        def generator_json(gen):
            """Converts a single TOUGH2 generator to Waiwera JSON dictionary."""

            mass_component = {'MASS': 1, 'MASD': 1, 'HEAT': num_eqns,
                              'COM1': 1, 'COM2': 2, 'COM3': 3, 'COM4': 4,
                              'COM5': 5, 'WATE': 1, 'AIR ': 2, 'TRAC': 2, 'NACL': 3}

            def specified_injection_generator_json(g, gen):
                """Generators which inject at a specified rate."""
                if tracer and gen.type in ['COM2', 'TRAC']:
                    g['tracer'] = gen.gx
                else:
                    g['rate'] = gen.gx
                    if gen.type == 'MASD': injection = False
                    else:
                        injection = gen.gx > 0. or \
                                    (gen.time and any([r > 0. for r in gen.rate]))
                    if injection:
                        g['component'] = mass_component[gen.type]
                        if gen.type != 'HEAT': g['enthalpy'] = gen.ex
                    else:
                        if gen.type == 'MASS':
                            g['separator'] = separator(gen.hg)
                        elif gen.type == 'MASD':
                            g['deliverability'] = {'productivity': gen.ex,
                                                   'pressure': gen.fg,
                                                   'threshold': gen.hg}
                            g['limiter'] = {'total': abs(gen.gx)}
                            g['separator'] = separator(gen.fg)
                            g['direction'] = 'production'
                return g

            def delv_generator_json(g, gen):
                """DELV generator type."""
                ltab = 0 if gen.ltab is None else gen.ltab
                if ltab > 1:
                    raise Exception('DELV generator with multiple layers not supported.')
                else:
                    g['deliverability'] = {'productivity': gen.gx,
                                           'pressure': gen.ex}
                    if gen.gx >= 0.:
                        g['direction'] = 'production'
                        g['separator'] = separator(gen.fg)
                    else:
                        g['direction'] = 'injection'
                        g['enthalpy'] = gen.fg
                return g

            def geothermal_deliverability_generator_json(g, gen):
                """Geothermal deliverability generator types - DELG etc."""
                g['deliverability'] = {'productivity': gen.gx,
                                       'pressure': gen.ex}
                g['separator'] = separator(gen.fg)
                if gen.hg is not None:
                    if gen.hg > 0.:
                        g['limiter'] = {limit_type[gen.type]: gen.hg}
                    elif gen.hg < 0. and gen.type in ['DELG', 'DMAK', 'DMAT']:
                        g['rate'] = gen.hg # initial rate for computing productivity index
                        del g['deliverability']['productivity']
                if gen.type == 'DELS': g['production_component'] = 2
                g['direction'] = 'production'
                return g

            def recharge_generator_json(g, gen):
                """Recharge generator type."""
                g['enthalpy'] = gen.ex
                if (gen.hg is not None) and gen.hg != 0.:
                    rech = {}
                    g['direction'] = "both"
                    if gen.fg is not None:
                        if gen.fg < 0.: g['direction'] = "out"
                        elif gen.fg > 0.: g['direction'] = "in"
                    if gen.hg > 0.: rech['pressure'] = gen.hg
                    else: rech['pressure'] = 'initial'
                    rech['coefficient'] = gen.gx
                    g['recharge'] = rech
                else:
                    g['rate'] = gen.gx
                return g

            def injectivity_generator_json(g, gen):
                """Generator types which inject against a pressure."""
                if gen.type == 'XINJ': g['enthalpy'] = gen.ex
                g['direction'] = 'injection'
                g['injectivity'] = {'pressure': gen.hg,
                                    'coefficient': abs(gen.fg)}
                if gen.gx > 0:
                    g['limiter'] = {'total': gen.gx}
                return g

            def table_generator_json(g, gen):
                """Generators with tables of values vs. time."""
                g['interpolation'] = interp_type
                g['averaging'] = averaging_type
                data_table = [list(r) for r in zip(gen.time, gen.rate)]
                if gen.type in ['DELG', 'DMAK', 'DMAT', 'DELT', 'DELW']:
                    ltab = 0 if gen.ltab is None else gen.ltab
                    if ltab > 0:
                        g['deliverability']['productivity'] = {'time': data_table}
                    else:
                        g['deliverability']['pressure'] = {'enthalpy': data_table}
                elif tracer and gen.type in ['COM2', 'TRAC']:
                    g['tracer'] = data_table
                else:
                    if gen.rate: g['rate'] = data_table
                    if gen.enthalpy:
                        g['enthalpy'] = [list(r) for r in zip(gen.time, gen.enthalpy)]
                return g

            if gen.block in geo.block_name_index:
                cell_index = geo.block_name_index[gen.block] - geo.num_atmosphere_blocks
                if cell_index < 0: cell_index = None
            else:
                cell_index = None
            g = {'name': unique_name(gen), 'cell': cell_index}

            if gen.type in mass_component:
                g = specified_injection_generator_json(g, gen)
            elif gen.type == 'DELV':
                g = delv_generator_json(g, gen)
            elif gen.type in ['DELG', 'DELS', 'DELT', 'DELW', 'DMAK', 'DMAT']:
                g = geothermal_deliverability_generator_json(g, gen)
            elif gen.type == 'RECH':
                g = recharge_generator_json(g, gen)
            elif gen.type in ['IMAK', 'XINJ']:
                g = injectivity_generator_json(g, gen)

            if gen.time:
                g = table_generator_json(g, gen)
            return g

        def tmak_json(g, gen, itmak, makeup_inputs):
            """TMAK (total makeup) group with limiter."""
            if gen.name.strip() == '':
                g['name'] = 'makeup %d' % itmak
            del g['cell']
            if gen.hg is None or gen.hg >= 0:
                raise Exception('Unscaled TMAK not supported.')
            elif gen.hg == -1: g['scaling'] = 'uniform'
            else: g['scaling'] = 'progressive'
            limiter = {}
            if gen.gx: limiter['total'] = abs(gen.gx)
            if gen.ex: limiter['steam'] = abs(gen.ex)
            if limiter: g['limiter'] = limiter
            g['in'] = makeup_inputs
            return g

        def reinjector_output_type(gen):
            """For FINJ, PINJ, RINJ and IMAK generators, returns water or steam
            output type."""
            if gen.type in ['FINJ', 'PINJ', 'RINJ']:
                output_type = 'water' if gen.hg > 0 else 'steam'
            elif gen.type == 'IMAK':
                output_type = 'water' if gen.fg > 0 else 'steam'
            else:
                raise Exception('Unrecognised reinjection generator type: %s' % gen.type)
            return output_type

        def reinjector_output_json(g, gen):
            """Returns JSON for reinjector output."""
            output = {'out': g['name']}
            output['enthalpy'] = gen.ex
            if gen.type == 'FINJ':
                output['rate'] = gen.gx
            elif gen.type in ['PINJ', 'RINJ']:
                output['proportion'] = abs(gen.hg)
            return output

        def has_outputs(reinjector):
            """Returns true if JSON has non-empty 'water' or 'steam' properties."""
            return reinjector['water'] or reinjector['steam']

        def prune_reinjector(reinjector):
            """Deletes empty keys from reinjector."""
            for key in ['water', 'steam']:
                if reinjector[key] == []: del reinjector[key]
            return reinjector

        sources, groups, reinjectors = [], [], []
        makeup_inputs, group_inputs = [], []
        itmak, ireinjector = 1, 1
        reinjection = False

        if self.generatorlist:
            for gen in self.generatorlist:

                if gen.type in unsupported_types:
                    raise Exception('Generator type ' + gen.type + ' not supported.')
                else:

                    g = generator_json(gen)
                    if gen.type != 'TMAK': sources.append(g)

                    if gen.type in ['DMAK', 'DMAT']:
                        makeup_inputs.append(g['name'])
                    elif gen.type in reinjection_contributors:
                        group_inputs.append(g['name'])
                    elif gen.type == 'TMAK':
                        tmak_subgroup = tmak_json(g, gen, itmak, makeup_inputs)
                        itmak += 1
                        makeup_inputs = []
                        groups.append(tmak_subgroup)
                        group_inputs.append(tmak_subgroup['name'])
                    elif gen.type in ['FINJ', 'PINJ', 'RINJ', 'IMAK']:
                        if not reinjection:
                            reinjection = True
                            if len(makeup_inputs) == 0 and len(group_inputs) == 1:
                                group_name = group_inputs[0]
                                reinjector_input_group = None
                            else:
                                group_name = 'reinjector group %d' % ireinjector
                                reinjector_input_group = {'name': group_name,
                                                          'in': group_inputs + makeup_inputs}
                            name = 'reinjector %d' % ireinjector
                            reinjector = {'name': name, 'in': group_name,
                                          'water': [], 'steam': []}
                            overflow_outputs = {'water': [], 'steam': []}
                        if reinjection:
                            output_json = reinjector_output_json(g, gen)
                            output_type = reinjector_output_type(gen)
                            if gen.type == 'RINJ':
                                overflow_outputs[output_type].append(output_json)
                            else:
                                reinjector[output_type].append(output_json)
                            if gen.type in ['FINJ', 'PINJ', 'RINJ'] and gen.fg != 0:
                                outputs = has_outputs(reinjector)
                                overflow = has_outputs(overflow_outputs)
                                if outputs or overflow:
                                    if reinjector_input_group: groups.append(reinjector_input_group)
                                    reinjector = prune_reinjector(reinjector)
                                    reinjectors.append(reinjector)
                                    ireinjector += 1
                                if overflow:
                                    name = 'reinjector %d' % ireinjector
                                    overflow_reinjector = {'name': name,
                                                           'water': overflow_outputs['water'],
                                                           'steam': overflow_outputs['steam']}
                                    reinjectors.append(overflow_reinjector)
                                    ireinjector += 1
                                    reinjector['overflow'] = overflow_reinjector['name']
                                reinjection = False
                                makeup_inputs, group_inputs = [], []

            if reinjection:
                # end of generator list without a reinjection reset:
                outputs = has_outputs(reinjector)
                overflow = has_outputs(overflow_outputs)
                if outputs or overflow:
                    if reinjector_input_group: groups.append(reinjector_input_group)
                    reinjectors.append(reinjector)
                if overflow:
                    name = 'reinjector %d' % ireinjector
                    overflow_reinjector = {'name': name,
                                           'water': overflow_outputs['water'],
                                           'steam': overflow_outputs['steam']}
                    reinjector = prune_reinjector(reinjector)
                    reinjectors.append(overflow_reinjector)
                    reinjector['overflow'] = overflow_reinjector['name']

        if sources: jsondata['source'] = sources
        network = {}
        if groups: network['group'] = groups
        if reinjectors: network['reinject'] = reinjectors
        if 'group' in network or 'reinject' in network:
            jsondata['network'] = network
        return jsondata

    def boundaries_json(self, geo, bdy_incons, atmos_volume, eos, mesh_coords, tracer = None):
        """Converts Dirichlet boundary conditions to Waiwera JSON dictionary.
        Currently connections to boundary blocks that are not either horizontal or
        vertical will not be converted correctly.
        """
        jsondata = {}
        vertical_tolerance = 1.e-6
        if eos in primary_to_region_funcs:
            primary_to_region = primary_to_region_funcs[eos]
            num_primary = waiwera_eos_num_primary[eos]
            jsondata['boundaries'] = []
            for blk in self.grid.blocklist:
                if not (0. < blk.volume < atmos_volume):
                    if isinstance(bdy_incons, t2incon):
                        pv = bdy_incons[blk.name].variable
                    else:
                        pv = bdy_incons
                    reg = primary_to_region(pv)
                    bc = {'primary': pv[:num_primary], 'region': reg, 'faces': []}
                    if tracer: bc['tracer'] = pv[num_primary]
                    for conname in blk.connection_name:
                        nz = -self.grid.connection[conname].dircos
                        vertical_connection = abs(nz) > vertical_tolerance
                        names = list(conname)
                        names.remove(blk.name)
                        interior_blkname = names[0]
                        interior_blk = self.grid.block[interior_blkname]
                        if 0. < interior_blk.volume < atmos_volume:
                            cell_index = geo.block_name_index[interior_blkname] - geo.num_atmosphere_blocks
                            if blk.centre is None:
                                if vertical_connection:
                                    normal = np.array([0., 0., nz])
                                else:
                                    raise Exception("Can't find normal vector for connection: " +
                                                    str(conname))
                            else:
                                normal = blk.centre - interior_blk.centre
                            normal /= np.linalg.norm(normal)
                            if mesh_coords != 'xyz':
                                if vertical_connection:
                                    if mesh_coords in ['xz', 'yz', 'rz']:
                                        normal = normal[[0,2]]
                                    elif mesh_coords == 'xy': normal = None
                                else: normal = normal[[0,1]]
                            if normal is not None:
                                bc['faces'].append({"cells": [cell_index],
                                                    "normal": list(normal)})
                    normals = np.array([spec['normal'] for spec in bc['faces']])
                    if np.isclose(normals, normals[0], rtol = 1.e-8).all():
                        allcells = []
                        for spec in bc['faces']:
                            allcells += spec['cells']
                        bc['faces'] = {"cells": allcells,
                                       "normal": bc['faces'][0]["normal"]}
                    if bc['faces']:
                        if isinstance(bc['faces'], list) and \
                           len(bc['faces']) == 1: bc['faces'] = bc['faces'][0]
                        jsondata['boundaries'].append(bc)

            if jsondata['boundaries']:
                # collapse down to one boundary if possible:
                primaries = np.array([bc['primary'] for bc in jsondata['boundaries']])
                if np.isclose(primaries, primaries[0], rtol = 1.e-8).all():
                    regions = np.array([bc['region'] for bc in jsondata['boundaries']])
                    if np.isclose(regions, regions[0]).all():
                        if tracer:
                            tracers = np.array([bc['tracer'] for bc in jsondata['boundaries']])
                            homog_tracer = np.isclose(tracers, tracers[0], rtol = 1.e-8).all()
                        else: homog_tracer = True
                        if homog_tracer:
                            normals = []
                            for bc in jsondata['boundaries']:
                                if isinstance(bc['faces'], dict):
                                    normals.append(bc['faces']['normal'])
                                else:
                                    for face in bc['faces']:
                                        normals.append(face['normal'])
                            normals = np.array(normals)
                            if np.isclose(normals, normals[0], rtol = 1.e-8).all():
                                allcells = []
                                for bc in jsondata['boundaries']:
                                    if isinstance(bc['faces'], dict):
                                        allcells += bc['faces']['cells']
                                    else:
                                        for face in bc['faces']:
                                            allcells += face['cells']
                                normal = list(normals[0, :])
                                primary = list(primaries[0,:])
                                region = int(regions[0])
                                jsondata['boundaries'] = [{"primary": primary, "region": region,
                                                           "faces": {"normal": normal,
                                                                     "cells": allcells}}]
                                if tracer:
                                    jsondata['boundaries'][0]['tracer'] = tracers[0]
        else:
            raise Exception("Finding thermodynamic region from primary variables not yet supported for EOS:" + eos)
        return jsondata

    def output_json(self):
        """Converts output specifications to Waiwera JSON dictionary."""
        datbase, ext = splitext(self.filename)
        jsondata = {}
        if self.parameter['print_interval'] is not None and \
           self.parameter['max_timesteps'] is not None and \
           self.parameter['print_interval'] >= self.parameter['max_timesteps']:
            print_interval = 0
        else:
            print_interval = self.parameter['print_interval']
        jsondata['output'] = {
            'filename': datbase + '.h5',
            'frequency': print_interval,
            'final': True,
            'initial': int(self.parameter['option'][24]) > 0}
        if self.output_times:
            time_tol = 1.e-8
            checkpoint = {'repeat': False}
            if 'num_times_specified' in self.output_times:
                num_times_specified = self.output_times['num_times_specified']
            else:
                num_times_specified = len(self.output_times['time'])
            if 'num_times' in self.output_times:
                num_times = self.output_times['num_times']
            else:
                num_times = num_times_specified
            if num_times_specified >= 0:
                times = self.output_times['time'][:num_times_specified]
                if 'time_increment' in self.output_times:
                    dt = self.output_times['time_increment']
                    if num_times_specified == 1 and abs(times[0] - dt) <= time_tol:
                        checkpoint['repeat'] = num_times
                    else:
                        for i in range(num_times - num_times_specified):
                            times.append(times[-1] + dt)
                checkpoint['time'] = times
            else: # time steps
                steps = self.output_times['time']
                if 'time_increment' in self.output_times:
                    dt = self.output_times['time_increment']
                    for i in range(num_times - abs(num_times_specified)):
                        steps.append(dt)
                checkpoint['step'] = steps
            if self.type == 'AUTOUGH2': checkpoint['tolerance'] = 0.1
            else: checkpoint['tolerance'] = 0.
            jsondata['output']['checkpoint'] = checkpoint
        return jsondata

    def json(self, geo, mesh_filename, atmos_volume = 1.e25, incons = None,
                    eos = None, bdy_incons = None, mesh_coords = 'xyz'):
        """Returns a Waiwera JSON dictionary representing the t2data object
        (with associated mulgrid geometry)."""

        jsondata = {}
        jsondata['title'] = self.title.strip()
        jsondata['gravity'] = self.parameter['gravity']
        jsondata['thermodynamics'] = 'ifc67'
        jsondata.update(self.mesh_json(geo, mesh_filename))
        eos_data, tracer_data = self.eos_json(eos)
        jsondata.update(eos_data)
        if tracer_data: jsondata.update(tracer_data)
        jsondata.update(self.timestepping_json())
        jsondata.update(self.output_json())
        jsondata.update(self.rocks_json(geo, atmos_volume, mesh_coords))
        jsondata['rock'].update(self.relative_permeability_json())
        jsondata['rock'].update(self.capillary_pressure_json())
        if isinstance(incons, str):
            effective_incs = incons
        else:
            effective_incs = self.effective_incons(incons)
        jsondata.update(self.initial_json(geo, effective_incs,
                                          jsondata['eos']['name'],
                                          tracer_data))
        if bdy_incons is None:
            bdy_incons = effective_incs
        jsondata.update(self.boundaries_json(geo, bdy_incons, atmos_volume,
                                             jsondata['eos']['name'],
                                             mesh_coords, tracer_data))
        jsondata.update(self.generators_json(geo, jsondata['eos']['name'],
                                             tracer_data))
        return jsondata
