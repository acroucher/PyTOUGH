"""For reading, parsing and writing fixed format text files."""

"""
Copyright 2013 University of Auckland.

This file is part of PyTOUGH.

PyTOUGH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PyTOUGH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PyTOUGH.  If not, see <http://www.gnu.org/licenses/>."""

class fixed_format_file(file):

    """Class for fixed format text file.  Values from the file may be parsed into variables, 
    according to a specification dictionary.  The keys of the specification dictionary are
    arbitrary and may be assigned for convenience, e.g. referring to specific sections or lines
    in the file.  Each value in the specification dictionary is a list of two lists: first a list
    of the names of variables in the specification, then a list of the corresponding format
    specifications.  The individual format specifications are like those in Python formats,
    consisting of an integer width value followed by a type ('d' for integer, 'f' for float etc.).
    The default conversion functions also allow an 'x' specifier for blanks (like fortran), which
    returns None."""

    default_conversion_function = {'d':int, 'f':float, 'e':float, 'g':float,
                                   's':lambda x:x.rstrip('\n'), 'x':lambda x:None}

    def __init__(self, filename, mode, specification, conversion_function = default_conversion_function):
        self.specification = specification
        self.setup_conversion_functions(conversion_function)
        self.preprocess_specification()
        super(fixed_format_file, self).__init__(filename, mode)

    def setup_conversion_functions(self, conversion_function):
        """Sets up conversion functions for parsing text.  The functions are wrapped with an exception handler
        to return None on ValueError."""
        def value_error_none(f):
            def fn(x):
                try: return f(x)
                except ValueError: return None
            return fn
        self.conversion_function = dict([(typ, value_error_none(f)) for typ,f in conversion_function.iteritems()])

    def preprocess_specification(self):
        """Pre-process specifications to speed up parsing."""
        self.line_spec, self.spec_width={}, {}
        for section, [names,specs] in self.specification.iteritems():
            self.line_spec[section] = []
            pos = 0
            for spec in specs:
                fmt, typ=spec[:-1], spec[-1]
                w = int(fmt.partition('.')[0])
                nextpos = pos+w
                self.line_spec[section].append(((pos,nextpos),typ))
                pos = nextpos
                self.spec_width[fmt] = w
        
    def parse_string(self, line, linetype):
        """Parses a string into values according to specified input format (d,f,s, or x for integer, float, string or skip).
        Blanks are converted to None."""
        return [self.conversion_function[typ](line[i1:i2]) for (i1,i2),typ in self.line_spec[linetype]]

    def write_values_to_string(self, vals, linetype):
        """Inverse of parse_string()."""
        fmt = self.specification[linetype][1]
        strs = []
        for val,f in zip(vals,fmt):
            if (val is not None) and (f[-1]<>'x'): valstr = ('%%%s'%f) % val
            else: valstr = ' '*self.spec_width[f[0:-1]] # blank
            strs.append(valstr)
        return ''.join(strs)

    def read_values(self, linetype):
        """Reads a line from the file, parses it and returns the values."""
        line = self.readline()
        return self.parse_string(line,linetype)

    def write_values(self, vals, linetype):
        """Inverse of read_values()."""
        line = self.write_values_to_string(vals,linetype)
        self.write('%s\n'%line)

    def read_value_line(self, variable, linetype):
        """Reads a line of parameter values from the file into a dictionary variable.
        Null values are ignored."""
        spec = self.specification[linetype]
        vals = self.read_values(linetype)
        for var,val in zip(spec[0],vals):
            if val is not None: variable[var] = val

    def write_value_line(self, variable, linetype):
        """Inverse of read_value_line()."""
        spec = self.specification[linetype]
        vals = []
        for name in spec[0]:
            if name in variable: val = variable[name]
            else: val = None
            vals.append(val)
        self.write_values(vals,linetype)
