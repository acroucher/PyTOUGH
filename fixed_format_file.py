"""For reading, parsing and writing fixed format text files.

Copyright 2013 University of Auckland.

This file is part of PyTOUGH.

PyTOUGH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PyTOUGH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PyTOUGH.  If not, see <http://www.gnu.org/licenses/>."""

from numpy import nan

def fortran_float(s, blank_value = 0.0):
    """Returns float of a string written by Fortran.
    Its behaviour is different from the float() function in the following ways:
    - a blank string will return the specified blank value (default zero)
    - embedded spaces are ignored
    - 'd' specifier will be treated the same as 'e'
    - underflow or overflow in exponent, with the 'e' omitted,
      are treated as if the 'e' was present
    If any other errors are encountered, np.nan is returned."""
    try: return float(s)
    except ValueError:
        s = s.strip()
        if not s: return blank_value
        else:
            try:
                s = s.lower().replace('d', 'e').replace(' ', '')
                return float(s)
            except:
                try:
                    return float(''.join([s[0], s[1:].replace('-', 'e-')]))
                except ValueError:
                    try:
                        return float(''.join([s[0], s[1:].replace('+', 'e')]))
                    except ValueError: return nan
    except: return nan

def fortran_int(s, blank_value = 0):
    """Returns float of a string written by Fortran.
    Its behaviour is different from the float() function in the following ways:
    - a blank string will return the specified blank value (default zero)
    - embedded spaces are ignored
    If any other errors are encountered, None is returned."""
    try: return int(s)
    except ValueError:
        s = s.strip()
        if not s: return blank_value
        else:
            try:
              s = s.replace(' ', '')
              return int(s)
            except: return None

def value_error_none(f):
    """Wraps a function with a handler to return None on a ValueError
    exception."""
    def fn(x):
        try: return f(x)
        except ValueError: return None
    return fn

default_read_float = value_error_none(float)
default_read_int = value_error_none(int)
default_read_str = value_error_none(lambda x: x.rstrip('\n'))
default_read_space = lambda x: None

from functools import partial
fortran_read_float = partial(fortran_float, blank_value = None)
fortran_read_int = partial(fortran_int, blank_value = None)

def read_function_dict(floatfn = default_read_float, intfn = default_read_int,
                             strfn = default_read_str, spacefn = default_read_space):
    """Returns a conversion function dictionary using the specified functions for float,
    int, string and space."""
    result = {'s': strfn, 'x': spacefn, 'd': intfn}
    for typ in ['f','e','g']: result[typ] = floatfn
    return result

default_read_function = read_function_dict()
fortran_read_function = read_function_dict(fortran_read_float, fortran_read_int)

class fixed_format_file(object):
    """Class for fixed format text file.  Values from the file may be
    parsed into variables, according to a specification dictionary.
    The keys of the specification dictionary are arbitrary and may be
    assigned for convenience, e.g. referring to specific sections or
    lines in the file.  Each value in the specification dictionary is
    a list of two lists: first a list of the names of variables in the
    specification, then a list of the corresponding format
    specifications.  The individual format specifications are like
    those in Python formats, consisting of an integer width value
    followed by a type ('d' for integer, 'f' for float etc.).  The
    default conversion functions also allow an 'x' specifier for
    blanks (like fortran), which returns None.
    """

    def __init__(self, filename, mode, specification,
                 read_function = default_read_function):
        self.specification = specification
        self.read_function = read_function
        self.preprocess_specification()
        self.file = open(filename, mode)

    def readline(self):
        """Returns next line from file."""
        return self.file.readline()

    def write(self, s):
        """Writes string s to file."""
        self.file.write(s)

    def close(self):
        """Closes file."""
        self.file.close()

    def preprocess_specification(self):
        """Pre-process specifications to speed up parsing."""
        self.line_spec, self.spec_width={}, {}
        for section, [names,specs] in self.specification.items():
            self.line_spec[section] = []
            pos = 0
            for spec in specs:
                fmt, typ=spec[:-1], spec[-1]
                w = int(fmt.partition('.')[0])
                nextpos = pos + w
                self.line_spec[section].append(((pos, nextpos), typ))
                pos = nextpos
                self.spec_width[fmt] = w
        
    def parse_string(self, line, linetype):
        """Parses a string into values according to specified input format
        (d,f,s, or x for integer, float, string or skip).  Blanks are
        converted to None.
        """
        return [self.read_function[typ](line[i1:i2]) for
                (i1, i2) , typ in self.line_spec[linetype]]

    def write_values_to_string(self, vals, linetype):
        """Inverse of parse_string()."""
        fmt = self.specification[linetype][1]
        strs = []
        for val , f in zip(vals , fmt):
            if (val is not None) and (f[-1] != 'x'): valstr = ('%%%s'%f) % val
            else: valstr = ' ' * self.spec_width[f[0:-1]] # blank
            strs.append(valstr)
        return ''.join(strs)

    def read_values(self, linetype):
        """Reads a line from the file, parses it and returns the values."""
        line = self.file.readline()
        return self.parse_string(line, linetype)

    def write_values(self, vals, linetype):
        """Inverse of read_values()."""
        line = self.write_values_to_string(vals , linetype)
        self.write('%s\n' % line)

    def read_value_line(self, variable, linetype):
        """Reads a line of parameter values from the file into a dictionary variable.
        Null values are ignored."""
        spec = self.specification[linetype]
        vals = self.read_values(linetype)
        for var , val in zip(spec[0] , vals):
            if val is not None: variable[var] = val

    def write_value_line(self, variable, linetype):
        """Inverse of read_value_line()."""
        spec = self.specification[linetype]
        vals = []
        for name in spec[0]:
            if name in variable: val = variable[name]
            else: val = None
            vals.append(val)
        self.write_values(vals, linetype)
