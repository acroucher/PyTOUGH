"""Class for TOUGH2 data"""

"""
Copyright 2011 University of Auckland.

This file is part of PyTOUGH.

PyTOUGH is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PyTOUGH is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PyTOUGH.  If not, see <http://www.gnu.org/licenses/>."""

from t2grids import *
from t2incons import *
from math import ceil

class t2datafile(file):
    """Class for TOUGH2 data file"""
    specification={
        'title':[['title'],['80s']],
        'simulator':[['simulator'],['80s']],
        'rocks1':[['name','nad','density','porosity','k1','k2','k3','conductivity','specific_heat'],
                  ['5s','5d']+['10.4e']*7],
        'rocks1.1':[['compressibility','expansivity','dry_conductivity','tortuosity','klinkenberg','xkd3','xkd4'], ['10.4e']*7],
        'rocks1.2':[['type','']+['parameter']*7,['5d','5x']+['10.3e']*7],
        'rocks1.3':[['type','']+['parameter']*7,['5d','5x']+['10.3e']*7],
        'param1_autough2':[['max_iterations','print_level','max_timesteps','max_duration','print_interval',
                   '_option_str','diff0','texp','be'],
                  ['2d']*2+['4d']*3+['24s']+['10.3e']*3],
        'param1':[['max_iterations','print_level','max_timesteps','max_duration','print_interval','_option_str','texp','be'],
                  ['2d']*2+['4d']*3+['24s']+['10.3e']*2],
        'param2':[['tstart','tstop','const_timestep','max_timestep','print_block','','gravity','timestep_reduction','scale'],
                  ['10.3e']*4+['5s','5x']+['10.4e']*3],
        'param3':[['relative_error','absolute_error','pivot','upstream_weight','newton_weight','derivative_increment'],
                  ['10.4e']*6],
        'timestep': [['timestep']*8,['10.4e']*8],
        'multi': [['num_components','num_equations','num_phases','num_secondary_parameters','num_inc'], ['5d']*5],
        'multi_autough2': [['num_components','num_equations','num_phases','num_secondary_parameters','eos'], ['5d']*4+['4s']],
        'lineq':[['type','epsilon','max_iterations','gauss','num_orthog'],['2d','10.4e','4d','1d','4d']],
        'default_incons':[['incon']*4,['20.14e']*4],
        'output_times1':[['num_times_specified','num_times','max_timestep','time_increment'],['5d']*2+['10.4e']*2],
        'output_times2':[['time']*8,['10.4e']*8],
        'relative_permeability':[['type','']+['parameter']*7,['5d','5x']+['10.3e']*7],
        'capillarity':[['type','']+['parameter']*7,['5d','5x']+['10.3e']*7],
        'blocks':[['name','nseq','nadd','rocktype','volume','ahtx','pmx','x','y','z'],['5s','5d','5d','5s']+['10.4e']*3+
                  ['10.3e']*3],
        'connections':[['block1','block2','nseq','nad1','nad2','direction','distance1','distance2','area','dircos','sigma'],
                       ['5s']*2+['5d']*4+['10.4e']*3+['10.7f','10.4e']],
        'generator':[['block','name','nseq','nadd','nads','ltab','','type','itab','gx','ex','hg','fg'],
                     ['5s']*2+['5d']*3+['5d','5x','4s','1s']+['10.3e']*4],
        'generation_times':[['time']*4,['14.7e']*4],
        'generation_rates':[['rate']*4,['14.7e']*4],
        'generation_enthalpy':[['enthalpy']*4,['14.7e']*4],
        'short':[['','frequency'],['5x','2d']],
        'incon1':[['block','','','porosity'],['5s']+['5x']*2+['15.9e']],
        'incon2':[['incon']*3,['20.14e']*3],
        'solver':[['type','','z_precond','','o_precond','relative_max_iterations','closure'],
                 ['1d','2x','2s','3x','2s']+['10.4e']*2],
        'indom2':[['indom']*4,['20.13e']*4],
        'diffusion': [['diff']*8,['10.4e']*8],
        'selec1': [['int_selec']*16,['5d']*16],
        'selec2': [['float_selec']*8,['10.4e']*8],
        'radii1': [['nrad'],['5d']],
        'radii2': [['radius']*8,['10.4e']*8],
        'equid' : [['nequ','','dr'],['5d','5x','10.4e']],
        'logar' : [['nlog','','rlog','dr'],['5d','5x']+['10.4e']*2],
        'layer1': [['nlay'],['5d']],
        'layer2': [['layer']*8,['10.4e']*8],
        'xyz1'  : [['deg'],['10.4e']],
        'xyz2'  : [['ntype','','no','del'],['2s','3x','5d','10.4e']],
        'xyz3'  : [['deli']*8,['10.4e']*8],
        'minc'  : [['part','type','','dual'],['5s']*2+['5x','5s']],
        'part1' : [['num_continua','nvol','where']+['spacing']*7,['3d']*2+['4s']+['10.4e']*7],
        'part2' : [['vol']*8,['10.4e']*8]
        }
    conversion_function={'d':int,'f':float,'e':float,'g':float,'s':str}
    def parse_string(self,line,linetype):
        """Parses a string into values according to specified input format (d,f,s, or x for integer, float, string or skip).
        Blanks are converted to None."""
        fmt=self.specification[linetype][1]
        result,pos=[],0
        for f in fmt:
            spec,typ=f[0:-1],f[-1]
            width=int(spec.split('.')[0])
            if typ=='x': val=None
            else:
                try: val=self.conversion_function[typ](line[pos:pos+width])
                except ValueError: val=None
            if typ=='s': val=val.replace('\n',' ')
            result.append(val)
            pos+=width
        return result
    def write_values_to_string(self,vals,linetype):
        """Inverse of parse_string()"""
        fmt=self.specification[linetype][1]
        s=""
        for val,f in zip(vals,fmt):
            if (val<>None) and (f[-1]<>'x'): valstr=('%'+f) % val
            else: # blank
                width=int(f[0:-1].split('.')[0])
                valstr=' '*width
            s=s+valstr
        return s
    def read_values(self,linetype):
        line=self.readline()
        return self.parse_string(line,linetype)
    def write_values(self,vals,linetype):
        line=self.write_values_to_string(vals,linetype)
        self.write(line+'\n')
    def read_value_line(self,variable,linetype):
        """Reads a line of parameter values into dictionary variable. Null values are ignored."""
        spec=self.specification[linetype]
        vals=self.read_values(linetype)
        for var,val in zip(spec[0],vals):
            if val<>None: variable[var]=val
    def write_value_line(self,variable,linetype):
        spec=self.specification[linetype]
        vals=[]
        for name in spec[0]:
            if name in variable: val=variable[name]
            else: val=None
            vals.append(val)
        self.write_values(vals,linetype)

class t2generator(object):
    """TOUGH2 generator (source or sink)"""
    def __init__(self,name='     ',block='     ',nseq=None,nadd=None,nads=None,type='MASS',
                 ltab=0,itab='',gx=0.0,ex=0.0,hg=0.0,fg=0.0,time=[],rate=[],enthalpy=[]):
        self.name=name
        self.block=block
        self.nseq,self.nadd,self.nads=nseq,nadd,nads
        self.type=type
        self.ltab=ltab
        self.itab=itab
        self.gx=gx
        self.ex=ex
        self.hg=hg
        self.fg=fg
        self.time=time
        self.rate=rate
        self.enthalpy=enthalpy
    def __repr__(self): return self.block+':'+self.name

default_parameters={'max_iterations':None, 'print_level':None, 'max_timesteps':None, 'max_duration':None, 'print_interval':None, 
                    '_option_str':'0'*24,'option':np.zeros(25,int8), 'diff0':None, 'texp':None, 'tstart':0.0, 'tstop':None,
                    'const_timestep':0.0,'timestep':[], 'max_timestep':None, 'print_block':None, 'gravity':0.0,
                    'timestep_reduction':None, 'scale':None, 'relative_error':None, 'absolute_error':None, 'pivot':None,
                    'upstream_weight':None, 'newton_weight':None, 'derivative_increment':None, 'default_incons':[]}

t2data_sections=['SIMUL','ROCKS','MESHM','PARAM','START','NOVER','RPCAP','LINEQ','SOLVR','MULTI','TIMES',
                 'SELEC','DIFFU','ELEME','CONNE','GENER','SHORT','FOFT','COFT','GOFT','INCON','INDOM']

class t2data(object):
    """Class for TOUGH2 data"""
    def __init__(self,filename=''):
        from copy import copy
        self.filename=filename
        self.title=''
        self.simulator=''
        self.parameter=copy(default_parameters)
        self.multi={}
        self.start=False
        self.relative_permeability={}
        self.capillarity={}
        self.lineq={}
        self.output_times={}
        self.grid=t2grid()
        self.generatorlist=[]
        self.generator={}
        self.short_output={}
        self.incon={}
        self.solver={}
        self.history_block=[]
        self.history_connection=[]
        self.history_generator=[]
        self.indom={}
        self.noversion=False
        self.diffusion=[]
        self.selection={}
        self.meshmaker=[]
        self._sections=[]
        self.end_keyword='ENDCY'
        if self.filename: self.read(filename)

    def __repr__(self): return self.title

    def run(self,save_filename='',incon_filename='',simulator='AUTOUGH2_2',silent=False):
        """Runs simulation using TOUGH2 or AUTOUGH2.  It's assumed that the data object has been written to file
        using write().  For AUTOUGH2, if the filenames for the save file or initial conditions file are not specified,
        they are constructed by changing the extensions of the data filename.  Set silent to True to suppress screen 
        output."""
        if self.filename:
            from os.path import splitext
            from os import devnull,system,remove
            datbase,ext=splitext(self.filename)
            if (self.type=='AUTOUGH2'):
                if save_filename=='': save_filename=datbase+'.save'
                if incon_filename=='': incon_filename=datbase+'.incon'
                savebase,ext=splitext(save_filename)
                inconbase,ext=splitext(incon_filename)
                # write a file containing the filenames, to pipe into AUTOUGH2:
                runfilename=simulator+'_input.dat'
                f=open(runfilename,'w')
                f.write(savebase+'\n')
                f.write(inconbase+'\n')
                f.write(datbase+'\n')
                f.close()
                if silent: out=' > '+devnull
                else: out=''
                # run AUTOUGH2:
                system(simulator+' < '+runfilename+out)
                remove(runfilename)
            else: # run TOUGH2 (need to specify simulator executable name)
                if silent: out=devnull
                else: out=datbase+'.listing'
                system(simulator+' < '+self.filename+' > '+out)

    def get_type(self):
        """Returns type (TOUGH2 or AUTOUGH2) based on whether the simulator has been set."""
        if self.simulator: return 'AUTOUGH2'
        else: return 'TOUGH2'
    type=property(get_type)

    def get_num_generators(self):
        return len(self.generatorlist)
    num_generators=property(get_num_generators)

    def generator_index(self,blocksourcenames):
        """Returns index of generator with specified tuple of block and source names."""
        if blocksourcenames in self.generator:
            return self.generatorlist.index(self.generator[blocksourcenames])
        else: return None

    def total_generation(self,type='MASS',name=''):
        """Returns array containing total generation in each block of the specified generator type and name.  The 
        name parameter specifies a regular expression to be matched."""
        import re
        tg=np.zeros(self.grid.num_blocks,float64)
        gens=[g for g in self.generatorlist if ((type==g.type) and re.search(name,g.name))]
        for g in gens: tg[self.grid.block_index(g.block)]+=g.gx
        return tg

    def specific_generation(self,type='MASS',name=''):
        """Returns array containing total specific generation (i.e. generation per unit volume) in each block of the
        specified generator type and name.  The name parameter specifies a regular expression to be matched."""
        import re
        tg=np.zeros(self.grid.num_blocks,float64)
        gens=[g for g in self.generatorlist if ((g.type==type) and (re.search(name,g.name)))]
        for g in gens:
            blkindex=self.grid.block_index(g.block)
            tg[blkindex]+=g.gx/self.grid.blocklist[blkindex].volume
        return tg

    def read_title(self,infile):
        """Reads simulation title"""
        infile.read_value_line(self.__dict__,'title')

    def write_title(self,outfile):
        outfile.write(self.title.strip()+'\n')

    def read_simulator(self,infile):
        """Reads simulator and EOS type"""
        infile.read_value_line(self.__dict__,'simulator')

    def write_simulator(self,outfile):
        if self.simulator:
            outfile.write('SIMUL\n')
            outfile.write(self.simulator.strip()+'\n')

    def read_rocktypes(self,infile):
        """Reads grid rock types"""
        self.grid.rocktypelist=[]
        self.grid.rocktype={}
        line=padstring(infile.readline())
        while line.strip():
            [name,nad,density,porosity,k1,k2,k3,conductivity,specific_heat]=infile.parse_string(line,'rocks1')
            self.grid.add_rocktype(rocktype(name,nad,density,porosity,[k1,k2,k3],conductivity,specific_heat))
            if nad>=1: # additional lines:
                infile.read_value_line(self.grid.rocktype[name].__dict__,'rocks1.1')
                if nad>=2:
                    vals=infile.read_values('rocks1.2')
                    self.grid.rocktype[name].relative_permeability['type']=vals[0]
                    self.grid.rocktype[name].relative_permeability['parameters']=vals[2:-1]
                    vals=infile.read_values('rocks1.3')
                    self.grid.rocktype[name].capillarity['type']=vals[0]
                    self.grid.rocktype[name].capillarity['parameters']=vals[2:-1]
            line=padstring(infile.readline())

    def write_rocktypes(self,outfile):
        outfile.write('ROCKS\n')
        for rt in self.grid.rocktypelist:
            vals=[rt.name,rt.nad,rt.density,rt.porosity]+list(rt.permeability)+[rt.conductivity,rt.specific_heat]
            outfile.write_values(vals,'rocks1')
            if rt.nad>=1:
                outfile.write_value_line(rt.__dict__,'rocks1.1')
                if rt.nad>=2:
                    vals=[rt.relative_permeability['type'],None]+rt.relative_permeability['parameters']
                    outfile.write_values(vals,'rocks1.2')
                    vals=[rt.capillarity['type'],None]+rt.capillarity['parameters']
                    outfile.write_values(vals,'rocks1.2')
        outfile.write('\n')

    def read_parameters(self,infile):
        """Reads simulation parameters"""
        spec=['param1','param1_autough2'][self.type=='AUTOUGH2']
        infile.read_value_line(self.parameter,spec)
        mops=ljust(self.parameter['_option_str'].rstrip(),24).replace(' ','0')
        self.parameter['option']=np.array([0]+[int(mop) for mop in mops],int8)
        infile.read_value_line(self.parameter,'param2')
        if (self.parameter['print_block']<>None) and (self.parameter['print_block'].strip()==''):
            self.parameter['print_block']=None
        self.read_timesteps(infile)
        infile.read_value_line(self.parameter,'param3')
        for val in infile.read_values('default_incons'): self.parameter['default_incons'].append(val)

    def write_parameters(self,outfile):
        outfile.write('PARAM\n')
        from copy import copy
        paramw=copy(self.parameter)
        if paramw['print_block']<>None: paramw['print_block']=unfix_blockname(paramw['print_block'])
        self.parameter['_option_str']=''.join([str(m) for m in self.parameter['option'][1:]])
        spec=['param1','param1_autough2'][self.type=='AUTOUGH2']
        outfile.write_value_line(self.parameter,spec)
        outfile.write_value_line(paramw,'param2')
        self.write_timesteps(outfile)
        outfile.write_value_line(self.parameter,'param3')
        outfile.write_values(self.parameter['default_incons'],'default_incons')

    def read_timesteps(self,infile):
        """Reads time step sizes from file"""
        if self.parameter['const_timestep']>=0.0:
            self.parameter['timestep']=[self.parameter['const_timestep']]
        else:
            nlines=-int(self.parameter['const_timestep'])
            self.parameter['timestep']=[]
            for i in xrange(nlines):
                for val in infile.read_values('timestep'):
                    if val<>None: self.parameter['timestep'].append(val)

    def write_timesteps(self,outfile):
        if self.parameter['const_timestep']<0.0:
            nlines=-int(self.parameter['const_timestep'])
            for i in xrange(nlines):
                i1,i2=i*8,min((i+1)*8,len(self.parameter['timestep']))
                vals=self.parameter['timestep'][i1:i2]
                if len(vals)<8: vals+=[None]*(8-len(vals))
                outfile.write_values(vals,'timestep')

    def read_multi(self,infile):
        """Reads EOS parameters"""
        spec=['multi','multi_autough2'][self.type=='AUTOUGH2']
        infile.read_value_line(self.multi,spec)
        if 'eos' in self.multi: self.multi['eos']=self.multi['eos'].strip()

    def write_multi(self,outfile):
        if self.multi<>{}:
            outfile.write('MULTI\n')
            spec=['multi','multi_autough2'][self.type=='AUTOUGH2']
            outfile.write_value_line(self.multi,spec)

    def read_start(self,infile):
        """Sets start parameter"""
        self.start=True

    def write_start(self,outfile):
        if self.start: outfile.write('START\n')

    def read_rpcap(self,infile):
        """Reads relative permeability and capillarity parameters"""
        vals=infile.read_values('relative_permeability')
        self.relative_permeability['type'],self.relative_permeability['parameters']=vals[0],vals[2:]
        vals=infile.read_values('capillarity')
        self.capillarity['type'],self.capillarity['parameters']=vals[0],vals[2:]

    def write_rpcap(self,outfile):
        if self.relative_permeability:
            outfile.write('RPCAP\n')
            vals=[self.relative_permeability['type'],None]+self.relative_permeability['parameters']
            outfile.write_values(vals,'relative_permeability')
            vals=[self.capillarity['type'],None]+self.capillarity['parameters']
            outfile.write_values(vals,'capillarity')

    def read_lineq(self,infile):
        """Reads linear equation parameters (AUTOUGH2)"""
        infile.read_value_line(self.lineq,'lineq')

    def write_lineq(self,outfile):
        if self.lineq:
            outfile.write('LINEQ\n')
            outfile.write_value_line(self.lineq,'lineq')

    def read_solver(self,infile):
        """Reads linear equation parameters (TOUGH2)"""
        infile.read_value_line(self.solver,'solver')

    def write_solver(self,outfile):
        if self.solver:
            outfile.write('SOLVR\n')
            outfile.write_value_line(self.solver,'solver')

    def read_blocks(self,infile):
        """Reads grid blocks"""
        self.grid.block,self.grid.blocklist={},[]
        line=padstring(infile.readline())
        while line.strip():
            [name,nseq,nadd,rockname,volume,ahtx,pmx,x,y,z]=infile.parse_string(line,'blocks')
            name=fix_blockname(name)
            rocktype=self.grid.rocktype[rockname]
            if (x<>None) and (y<>None) and (z<>None): centre=np.array([x,y,z])
            else: centre=None
            self.grid.add_block(t2block(name,volume,rocktype,centre=centre,ahtx=ahtx,pmx=pmx))
            if nseq>0:
                self.grid.block[name].nseq=nseq
                self.grid.block[name].nadd=nadd
            line=padstring(infile.readline())

    def write_blocks(self,outfile):
        if self.grid.num_blocks>0:
            outfile.write('ELEME\n')
            from copy import copy
            for blk in self.grid.blocklist:
                blkw=copy(blk.__dict__)
                blkw['name']=unfix_blockname(blkw['name'])
                if blk.centre==None: outfile.write_value_line(blkw,'blocks')
                else:
                    vals=[blkw['name'],blk.nseq,blk.nadd,blk.rocktype.name,blk.volume,
                          blk.ahtx,blk.pmx]+list(blk.centre)
                    outfile.write_values(vals,'blocks')
            outfile.write('\n')

    def read_connections(self,infile):
        """Reads grid connections"""
        self.grid.connectionlist,self.grid.connection=[],{}
        line=padstring(infile.readline())
        while line.strip():
            [name1,name2,nseq,nad1,nad2,isot,d1,d2,areax,betax,sigx]=infile.parse_string(line,'connections')
            name1,name2=fix_blockname(name1),fix_blockname(name2)
            self.grid.add_connection(t2connection([self.grid.block[name1],self.grid.block[name2]],isot,[d1,d2],areax,betax,sigx))
            if nseq>=1:
                self.grid.connection[(name1,name2)].nseq=nseq
                self.grid.connection[(name1,name2)].nad1=nad1
                self.grid.connection[(name1,name2)].nad2=nad2
            line=padstring(infile.readline())

    def write_connections(self,outfile):
        if self.grid.num_connections>0:
            outfile.write('CONNE\n')
            for con in self.grid.connectionlist:
                vals=[unfix_blockname(con.block[0].name),unfix_blockname(con.block[1].name),
                      con.nseq,con.nad1,con.nad2,con.direction]+con.distance+[con.area,con.dircos]
                outfile.write_values(vals,'connections')
            outfile.write('\n')

    def add_generator(self,generator=t2generator()):
        """Adds a generator."""
        self.generatorlist.append(generator)
        self.generator[(generator.block,generator.name)]=self.generatorlist[-1]

    def delete_generator(self,blocksourcenames):
        i=self.generator_index(blocksourcenames)
        del self.generator[blocksourcenames]
        del self.generatorlist[i]

    def clear_generators(self):
        self.generator.clear()
        del self.generatorlist[:]

    def delete_orphan_generators(self):
        """Deletes any generators specified in blocks which do not exist in the grid."""
        delgens=[]
        for gen in self.generatorlist:
            if not (gen.block in self.grid.block): delgens.append((gen.block,gen.name))
        for bg in delgens: self.delete_generator(bg)

    def read_generator(self,line,infile):
        """Returns generator read from line in file"""
        from string import rjust
        [block,name,nseq,nadd,nads,ltab,empty,gentype,itab,gx,ex,hg,fg]=infile.parse_string(line,'generator')
        block,name=fix_blockname(block),fix_blockname(name)
        time,rate,enthalpy=[],[],[]
        if ltab and gentype<>'DELV':
            ntimes=abs(ltab)
            if ntimes>1:
                nlines=int(ceil(ntimes/4.))
                for i in xrange(nlines):
                    for val in infile.read_values('generation_times'): 
                        if val<>None: time.append(val)
                for i in xrange(nlines):
                    for val in infile.read_values('generation_rates'):
                        if val<>None: rate.append(val)
                if itab.strip():
                    for i in xrange(nlines):
                        for val in infile.read_values('generation_enthalpy'):
                            if val<>None: enthalpy.append(val)
        return t2generator(name=name,block=block,nseq=nseq,nadd=nadd,nads=nads,type=gentype,ltab=ltab,itab=itab,
                           gx=gx,ex=ex,hg=hg,fg=fg,time=time,rate=rate,enthalpy=enthalpy)

    def write_generator(self,gen,outfile):
        from copy import copy
        genw=copy(gen.__dict__)
        genw['name'],genw['block']=unfix_blockname(genw['name']),unfix_blockname(genw['block'])
        outfile.write_value_line(genw,'generator')
        if gen.ltab and gen.type<>'DELV': ntimes=abs(gen.ltab)
        else: ntimes=1
        if ntimes>1:
            nlines=int(ceil(ntimes/4.))
            for i in xrange(nlines):
                i1,i2=i*4,min((i+1)*4,ntimes)
                vals=gen.time[i1:i2]
                if len(vals)<4: vals+=[None]*(4-len(vals))
                outfile.write_values(vals,'generation_times')
            for i in xrange(nlines):
                i1,i2=i*4,min((i+1)*4,ntimes)
                vals=gen.rate[i1:i2]
                if len(vals)<4: vals+=[None]*(4-len(vals))
                outfile.write_values(vals,'generation_rates')
            if gen.enthalpy:
                for i in xrange(nlines):
                    i1,i2=i*4,min((i+1)*4,ntimes)
                    vals=gen.enthalpy[i1:i2]
                    if len(vals)<4: vals+=[None]*(4-len(vals))
                    outfile.write_values(vals,'generation_enthalpy')

    def read_generators(self,infile):
        """Reads generators from file"""
        self.generatorlist=[]
        line=infile.readline()
        while line.strip():
            self.add_generator(self.read_generator(line,infile))
            line=infile.readline()

    def write_generators(self,outfile):
        if self.generatorlist:
            outfile.write('GENER\n')
            for generator in self.generatorlist:
                self.write_generator(generator,outfile)
            outfile.write('\n')
        
    def read_times(self,infile):
        """Reads output times from file"""
        infile.read_value_line(self.output_times,'output_times1')
        self.output_times['time']=[]
        nlines=int(ceil(self.output_times['num_times_specified']/8.))
        for i in xrange(nlines):
            for val in infile.read_values('output_times2'): 
                if val<>None: self.output_times['time'].append(val)

    def write_times(self,outfile):
        if self.output_times:
            outfile.write('TIMES\n')
            outfile.write_value_line(self.output_times,'output_times1')
            nlines=int(ceil(self.output_times['num_times_specified']/8.))
            for i in xrange(nlines):
                i1,i2=i*8,min((i+1)*8,len(self.output_times['time']))
                vals=self.output_times['time'][i1:i2]
                if len(vals)<8: vals+=[None]*(8-len(vals))
                outfile.write_values(vals,'output_times2')
        
    def read_incons(self,infile):
        """Reads initial conditions from file"""
        line=infile.readline()
        while line.strip():
            [blockname,empty,empty,porosity]=infile.parse_string(line,'incon1')
            blockname=fix_blockname(blockname)
            variables=infile.read_values('incon2')
            self.incon[blockname]=[porosity,variables]
            line=infile.readline()

    def write_incons(self,outfile):
        if self.incon:
            outfile.write('INCON\n')
            for blkname,inc in self.incon.iteritems():
                vals=[unfix_blockname(blkname),inc[0]]
                outfile.write_values(vals,'incon1')
                outfile.write_values(inc[1],'incon2')
            outfile.write('\n')
        
    def read_short_blocks(self,infile):
        """Reads short output blocks"""
        self.short_output['block']=[]
        badblocks=[]
        more=True
        while more:
            line=infile.readline()
            if line.strip():
                if line[0:5] in ['ELEME','CONNE','GENER']: more=False
                else:
                    blockname=fix_blockname(line[0:5])
                    if blockname in self.grid.block: self.short_output['block'].append(self.grid.block[blockname])
                    else: badblocks.append(blockname)
            else: more=False
        if len(badblocks)>0: print 'Short output blocks',badblocks,'do not exist and will be ignored.'
        return line

    def read_short_connections(self,infile):
        """Reads short output connections"""
        self.short_output['connection']=[]
        badcons=[]
        more=True
        while more:
            line=infile.readline()
            if line.strip():
                if line[0:5] in ['ELEME','CONNE','GENER']: more=False
                else:
                    blknames=(fix_blockname(line[0:5]),fix_blockname(line[5:10]))
                    if blknames in self.grid.connection: self.short_output['connection'].append(self.grid.connection[blknames])
                    else: badcons.append(blknames)
            else: more=False
        if len(badcons)>0: print 'Short output connections',badcons,'do not exist and will be ignored.'
        return line

    def read_short_generators(self,infile):
        """Reads short output generators"""
        self.short_output['generator']=[]
        badgens=[]
        more=True
        while more:
            line=infile.readline()
            if line.strip():
                if line[0:5] in ['ELEME','CONNE','GENER']: more=False
                else:
                    blksourcenames=(fix_blockname(line[0:5]),fix_blockname(line[5:10]))
                    if blksourcenames in self.generator: self.short_output['generator'].append(self.generator[blksourcenames])
                    else: badgens.append(blksourcenames)
            else: more=False
        if len(badgens)>0: print 'Short output generators',badgens,'do not exist and will be ignored.'
        return line

    def read_short_output(self,infile,headerline):
        """Reads short output specifications from file.  'headerline' is passed in to read the frequency parameter."""
        vals=infile.parse_string(headerline,'short')
        if len(vals)>1: self.short_output['frequency']=vals[1]
        read_fn={'ELEME':self.read_short_blocks,'CONNE':self.read_short_connections,'GENER':self.read_short_generators}
        more=True
        line=infile.readline()
        while more:
            if line.strip():
                keyword=line[0:5]
                line=read_fn[keyword](infile)
            else: more=False

    def write_short_output(self,outfile):
        if self.short_output:
            outfile.write('SHORT')
            if 'frequency' in self.short_output:
                if self.short_output['frequency']: outfile.write('%2d' % self.short_output['frequency'])
            outfile.write('\n')
            if 'block' in self.short_output:
                outfile.write('ELEME\n')
                for block in self.short_output['block']:outfile.write(unfix_blockname(block.name)+'\n')
            if 'connection' in self.short_output:
                outfile.write('CONNE\n')
                for con in self.short_output['connection']:outfile.write(unfix_blockname(con.block[0].name)+
                                                                         unfix_blockname(con.block[1].name)+'\n')
            if 'generator' in self.short_output:
                outfile.write('GENER\n')
                for gen in self.short_output['generator']:outfile.write(unfix_blockname(gen.block)+
                                                                        unfix_blockname(gen.name)+'\n')
            outfile.write('\n')

    def read_history_blocks(self,infile):
        """Reads history blocks (TOUGH2)"""
        self.history_block=[]
        badblocks=[]
        line=infile.readline()
        while line.strip():
            blockname=fix_blockname(line[0:5])
            if blockname in self.grid.block: self.history_block.append(self.grid.block[blockname])
            else: badblocks.append(blockname)
            line=infile.readline()
        if len(badblocks)>0: print 'History blocks',badblocks,'do not exist and will be ignored.'
        
    def read_history_connections(self,infile):
        """Reads history connections (TOUGH2)"""
        self.history_connection=[]
        badcons=[]
        line=infile.readline()
        while line.strip():
            blknames=(fix_blockname(line[0:5]),fix_blockname(line[5:10]))
            if blknames in self.grid.connection: self.history_connection.append(self.grid.connection[blknames])
            else: badcons.append(blknames)
            line=infile.readline()
        if len(badcons)>0: print 'History connections',badcons,'do not exist and will be ignored.'

    def read_history_generators(self,infile):
        """Reads history generators (TOUGH2)"""
        self.history_generator=[]
        badgens=[]
        line=infile.readline()
        while line.strip():
            blockname=fix_blockname(line[0:5])
            if blockname in self.grid.block: self.history_generator.append(self.grid.block[blockname])
            else: badgens.append(blockname)
            line=infile.readline()
        if len(badgens)>0: print 'History generator blocks',badgens,'do not exist and will be ignored.'

    def write_history_blocks(self,outfile):
        if self.history_block:
            outfile.write('FOFT\n')
            for block in self.history_block: outfile.write(unfix_blockname(block.name)+'\n')
            outfile.write('\n')
        
    def write_history_connections(self,outfile):
        if self.history_connection:
            outfile.write('COFT\n')
            for con in self.history_connection:
                outfile.write(unfix_blockname(con.block[0].name)+unfix_blockname(con.block[1].name)+'\n')
            outfile.write('\n')
        
    def write_history_generators(self,outfile):
        if self.history_generator:
            outfile.write('GOFT\n')
            for blk in self.history_generator: outfile.write(unfix_blockname(blk.name)+'\n')
            outfile.write('\n')

    def read_indom(self,infile):
        """Reads rock-specific initial conditions from file"""
        line=infile.readline()
        while line.strip():
            rockname=line[0:5]
            variables=infile.read_values('incon2')
            self.indom[rockname]=variables
            line=infile.readline()

    def write_indom(self,outfile):
        if self.indom:
            outfile.write('INDOM\n')
            for rockname,inc in self.indom.iteritems():
                outfile.write(rockname+'\n')
                outfile.write_values(inc,'incon2')
            outfile.write('\n')

    def read_noversion(self,infile):
        """Sets noversion parameter"""
        self.noversion=True

    def write_noversion(self,outfile):
        if self.noversion: outfile.write('NOVER\n')

    def read_diffusion(self,infile):
        """Reads diffusion coefficients from file"""
        if ('num_components' in self.multi) and ('num_phases' in self.multi):
            for comp in xrange(self.multi['num_components']):
                diffs=infile.read_values('diffusion')[0:self.multi['num_phases']]
                self.diffusion.append(diffs)
        else: print 'Unable to read DIFFU block: no MULTI block specified.'

    def write_diffusion(self,outfile):
        if self.diffusion:
            outfile.write('DIFFU\n')
            for comp in self.diffusion: outfile.write_values(comp,'diffusion')

    def read_selection(self,infile):
        """Reads selection parameters from file"""
        int_selec=infile.read_values('selec1')
        self.selection['integer']=int_selec
        nlines=int_selec[0]
        float_selec=[]
        for i in xrange(nlines): float_selec+=infile.read_values('selec2')
        self.selection['float']=float_selec
        
    def write_selection(self,outfile):
        if self.selection:
            outfile.write('SELEC\n')
            outfile.write_values(self.selection['integer'],'selec1')
            nlines=self.selection['integer'][0]
            for i in xrange(nlines):
                i1,i2=i*8,min((i+1)*8,len(self.selection['float']))
                vals=self.selection['float'][i1:i2]
                if len(vals)<8: vals+=[None]*(8-len(vals))
                outfile.write_values(vals,'selec2')

    def read_meshmaker(self,infile):
        """Reads meshmaker data"""
        read_fn={'RZ2D' : self.read_meshmaker_rz2d, 'XYZ'  : self.read_meshmaker_xyz,
                 'MINC' : self.read_meshmaker_minc}
        more=True
        while more:
            line=infile.readline()
            if line.strip():
                keyword=line[0:5].strip()
                if keyword in read_fn: read_fn[keyword](infile)
            else: more=False

    def write_meshmaker(self,outfile):
        if self.meshmaker:
            outfile.write('MESHMAKER\n')
            write_fn={'rz2d': self.write_meshmaker_rz2d, 'xyz': self.write_meshmaker_xyz,
                      'minc': self.write_meshmaker_minc}
            for (stype,section) in self.meshmaker: write_fn[stype.lower()](section,outfile)
            outfile.write('\n')

    def read_meshmaker_rz2d(self,infile):
        """Reads RZ2D meshmaker data"""
        section=('rz2d',[])
        more=True
        while more:
            line=infile.readline()
            keyword=line[0:5].strip()
            subsection={}
            if keyword=='RADII':
                nrad=infile.read_values('radii1')[0]
                subsection['radii']=[]
                nlines=int(ceil(nrad/8.))
                for i in xrange(nlines):
                    for val in infile.read_values('radii2'): 
                        if val<>None: subsection['radii'].append(val)
            elif keyword=='EQUID': infile.read_value_line(subsection,'equid')
            elif keyword=='LOGAR': infile.read_value_line(subsection,'logar')
            elif keyword=='LAYER': 
                nlayers=infile.read_values('layer1')[0]
                nlines=int(ceil(nlayers/8.))
                layer=[]
                for i in xrange(nlines): layer+=infile.read_values('layer2')
                subsection['layer']=layer[0:nlayers]
                more=False # LAYER indicates end of RZ2D
            if subsection: section[1].append((keyword.lower(),subsection))
        self.meshmaker.append(section)

    def write_meshmaker_rz2d(self,section,outfile):
        outfile.write('RZ2D\n')
        for stype,subsection in section:
            outfile.write(stype.upper()+'\n')
            if stype=='radii':
                nrad=len(subsection['radii'])
                outfile.write_values([nrad],'radii1')
                nlines=int(ceil(nrad/8.))
                for i in xrange(nlines):
                    i1,i2=i*8,min((i+1)*8,nrad)
                    vals=subsection['radii'][i1:i2]
                    if len(vals)<8: vals+=[None]*(8-len(vals))
                    outfile.write_values(vals,'radii2')
            elif stype=='equid': outfile.write_value_line(subsection,'equid')
            elif stype=='logar': outfile.write_value_line(subsection,'logar')
            elif stype=='layer':
                nlayers=len(subsection['layer'])
                outfile.write_values([nlayers],'layer1')
                nlines=int(ceil(nlayers/8.))
                for i in xrange(nlines):
                    i1,i2=i*8,min((i+1)*8,nrad)
                    vals=subsection['layer'][i1:i2]
                    if len(vals)<8: vals+=[None]*(8-len(vals))
                    outfile.write_values(vals,'layer2')

    def read_meshmaker_xyz(self,infile):
        """Reads XYZ meshmaker data"""
        deg=infile.read_values('xyz1')[0]
        section=('xyz',[deg])
        more=True
        while more:
            subsection={}
            line=infile.readline()
            if line.strip():
                [subsection['ntype'],blank,subsection['no'],subsection['del']]=infile.parse_string(line,'xyz2')
                if subsection['del']==0:
                    nlines=int(ceil(subsection['no']/8.))
                    deli=[]
                    for i in xrange(nlines): deli+=infile.read_values('xyz3')
                    subsection['deli']=deli[0:subsection['no']]
                section[1].append(subsection)
            else: more=False
        self.meshmaker.append(section)

    def write_meshmaker_xyz(self,section,outfile):
        outfile.write('XYZ\n')
        deg=section[0]
        outfile.write_values([deg],'xyz1')
        for subsection in section[1:]:
            outfile.write_value_line(subsection,'xyz2')
            if subsection['del']==0:
                nlines=int(ceil(subsection['no']/8.))
                for i in xrange(nlines):
                    i1,i2=i*8,min((i+1)*8,subsection['no'])
                    vals=subsection['deli'][i1:i2]
                    if len(vals)<8: vals+=[None]*(8-len(vals))
                    outfile.write_values(vals,'xyz3')
        outfile.write('\n')

    def read_meshmaker_minc(self,infile):
        """Reads MINC meshmaker data"""
        line=infile.readline().strip()
        keyword=line[0:5].strip()
        if keyword=='PART':
            subsection={}
            [part,subsection['type'],dummy,subsection['dual']]=infile.parse_string(line,'minc')
            vals=infile.read_values('part1')
            subsection['num_continua'],nvol,subsection['where'],subsection['spacing']=vals[0],vals[1],vals[2],vals[3:]
            nlines=int(ceil(nvol/8.))
            vol=[]
            for i in xrange(nlines): vol+=infile.read_values('part2')
            subsection['vol']=vol[0:nvol]
            self.meshmaker.append(('minc',subsection))

    def write_meshmaker_minc(self,section,outfile):
        outfile.write('MINC\n')
        outfile.write_values(['PART ',section['type'],'',section['dual']],'minc')
        nvol=len(section['vol'])
        outfile.write_values([section['num_continua'],nvol,section['where']]+section['spacing'],'part1')
        nlines=int(ceil(nvol/8.))
        for i in xrange(nlines):
            i1,i2=i*8,min((i+1)*8,section['vol'])
            vals=section['vol'][i1:i2]
            if len(vals)<8: vals+=[None]*(8-len(vals))
            outfile.write_values(vals,'part2')

    def read(self,filename=''):
        """Reads data from file"""
        if filename: self.filename=filename
        infile=t2datafile(self.filename)
        read_fn=dict(zip(t2data_sections,
                         [self.read_simulator, self.read_rocktypes, self.read_meshmaker, self.read_parameters, self.read_start, 
                         self.read_noversion, self.read_rpcap, self.read_lineq, self.read_solver, self.read_multi, self.read_times,
                         self.read_selection, self.read_diffusion, self.read_blocks, self.read_connections,
                         self.read_generators, self.read_short_output, self.read_history_blocks, 
                         self.read_history_connections, self.read_history_generators, self.read_incons, self.read_indom]))
        self.read_title(infile)
        self._sections=[]
        more=True
        while more:
            line=infile.readline()
            if line:
                keyword=line[0:5].strip()
                if keyword in ['ENDCY','ENDFI']:
                    more=False
                    self.end_keyword=keyword
                elif keyword in t2data_sections:
                    fn=read_fn[keyword]
                    if keyword=='SHORT': fn(infile,line)
                    else: fn(infile)
                    self._sections.append(keyword)
            else: more=False
        infile.close()
        return self

    def write(self,filename=''):
        """Writes data to file"""
        if filename: self.filename=filename
        if self.filename=='': self.filename='t2data.dat'
        outfile=t2datafile(self.filename,'w')
        write_fn=dict(zip(t2data_sections,
                          [self.write_simulator, self.write_rocktypes, self.write_meshmaker, self.write_parameters,
                           self.write_start, self.write_noversion, self.write_rpcap, self.write_lineq, self.write_solver,
                           self.write_multi, self.write_times, self.write_selection, self.write_diffusion, self.write_blocks,
                           self.write_connections, self.write_generators, self.write_short_output, self.write_history_blocks,
                           self.write_history_connections, self.write_history_generators, self.write_incons, self.write_indom]))
        self.write_title(outfile)
        for keyword in self._sections: write_fn[keyword](outfile)
        extra_sections=[keyword for keyword in t2data_sections if not (keyword in self._sections)]
        for keyword in extra_sections: write_fn[keyword](outfile)
        outfile.write(self.end_keyword+'\n')
        outfile.close()
    
    def transfer_rocktypes_from(self,source,mapping):
        """Transfers rock types (definitions and assignments) from another t2data object, using the specified block mapping."""
        from copy import deepcopy
        self.grid.rocktypelist=deepcopy(source.grid.rocktypelist)
        self.grid.rocktype=deepcopy(source.grid.rocktype)
        for blk in self.grid.blocklist: blk.rocktype=self.grid.rocktype[source.grid.block[mapping[blk.name]].rocktype.name]

    def transfer_generators_from(self,source,sourcegeo,geo,top_generator=[],bottom_generator=[],mapping={},colmapping={},rename=False):
        """Transfers generators from another t2data object, using the specified top and bottom generator lists and 
        optional block and column mappings.  If the rename parameter is False, generators other than those at the top or
        bottom of the model will keep their original names; if True, these generators will be renamed according to their
        new column names."""
        from copy import deepcopy
        tablegens=[' AIR','COM1','COM2','COM3','COM4','COM5','HEAT','MASS','NACL','TRAC',' VOL']
        if (colmapping=={}) or (mapping=={}): mapping,colmapping=sourcegeo.block_mapping(geo,True)
        bbox=sourcegeo.bounds
        qt=quadtree(bbox,sourcegeo.columnlist)
        incols=[col for col in geo.columnlist if sourcegeo.column_containing_point(col.centre,qtree=qt)<>None]
        self.clear_generators()
        col_generator=top_generator+bottom_generator
        for sourcegen in source.generatorlist:
            sourcecategory=sourcegeo.layer_name(sourcegen.name)
            sourcecolname=sourcegeo.column_name(sourcegen.block)
            if sourcecategory in col_generator:
                mappedcols=[col for col in incols if colmapping[col.name]==sourcecolname]
                mappedcolarea=sum([col.area for col in mappedcols])
                for col in mappedcols:
                    gen=deepcopy(sourcegen)
                    area_ratio=col.area/mappedcolarea
                    if gen.ltab: ntimes=abs(gen.ltab)
                    else: ntimes=1
                    if gen.type in tablegens:
                        if gen.gx: gen.gx*=area_ratio
                        if ntimes>1: gen.rate=[rate*area_ratio for rate in gen.rate]
                    if geo.convention==sourcegeo.convention: category=sourcecategory
                    else: category=['%2d'%col_generator.index(sourcecategory),sourcecategory,sourcecategory][geo.convention]
                    gen.name=geo.block_name(category,col.name)
                    if sourcecategory in top_generator: layername=geo.layerlist[geo.num_layers-col.num_layers].name
                    elif sourcecategory in bottom_generator: layername=geo.layerlist[-1].name
                    gen.block=geo.block_name(layername,col.name)
                    self.add_generator(gen)
            else: # other generators, do block by block:
                sourceblock=source.grid.block[sourcegen.block]
                mappedblocks=[blk for blk in self.grid.blocklist if mapping[blk.name]==sourceblock.name]
                mappedblockvol=sum([blk.volume for blk in mappedblocks])
                for blk in mappedblocks:
                    gen=deepcopy(sourcegen)
                    if gen.ltab: ntimes=abs(gen.ltab)
                    else: ntimes=1
                    vol_ratio=blk.volume/mappedblockvol
                    if gen.type in tablegens:
                        if gen.gx: gen.gx*=vol_ratio
                        if ntimes>1: gen.rate=[rate*vol_ratio for rate in gen.rate]
                    if rename:
                        if geo.convention==sourcegeo.convention: category=sourcecategory
                        else: category=[' 0',sourcecategory,sourcecategory][geo.convention]
                        colname=geo.column_name(blk.name)
                        gen.name=geo.block_name(category,colname)
                    gen.block=blk.name
                    self.add_generator(gen)
                    
    def transfer_from(self,source,sourcegeo,geo,top_generator=[],bottom_generator=[],sourceinconfilename='',inconfilename='',rename_generators=False):
        """Copies parameters, rock types and assignments, generators and initial conditions
        from another t2data object (without altering the grid structure).
        The top_generator and bottom_generator lists specify the generators which are to be kept at the top or
        bottom of the model, respectively.  They contain the 'layer' part of the generator name for top and bottom
        generators.
        If both the inconfilename parameters are specified, a new initial conditions file with filename 'inconfilename'
        is written to disk, with initial conditions transferred from the file 'sourceinconfilename'.
        If the rename_generators parameter is False, generators (other than those at the top and bottom of the model) keep
        their original names- otherwise, they are renamed according to their new column names. """
        mapping,colmapping=sourcegeo.block_mapping(geo,True)
        from copy import copy
        self.grid=t2grid().fromgeo(geo)
        self.simulator=source.simulator
        self.parameter=copy(source.parameter)
        if self.parameter['print_block']<>None:
            mappedblocks=[blk for blk in self.grid.blocklist if mapping[blk.name]==self.parameter['print_block']]
            if len(mappedblocks)>0: self.parameter['print_block']=mappedblocks[0].name
            else: self.parameter['print_block']=None
        self.multi=copy(source.multi)
        self.start=source.start
        self.noversion=source.noversion
        self.relative_permeability=copy(source.relative_permeability)
        self.capillarity=copy(source.capillarity)
        self.lineq=copy(source.lineq)
        self.solver=copy(source.solver)
        self.diffusion=copy(source.diffusion)
        self.selection=copy(source.selection)
        self.output_times=copy(source.output_times)
        self.transfer_rocktypes_from(source,mapping)
        self.transfer_generators_from(source,sourcegeo,geo,top_generator,bottom_generator,mapping,colmapping,rename_generators)
        # short output (these can't really be transferred sensibly):
        self.short_output={}
        self.history_block={}
        self.history_connection={}
        self.history_generator={}
        # incons (within data file):
        self.incon={}
        for blkname,inc in source.incon.items():
            mappedblocks=[blk for blk in self.grid.blocklist if mapping[blk.name]==blkname]
            for blk in mappedblocks: self.incon[blk.name]=inc
        self.indom=copy(source.indom)
        # incon file:
        if (sourceinconfilename<>'') and (inconfilename<>''):
            sourceinc=t2incon(sourceinconfilename)
            inc=t2incon()
            inc.transfer_from(sourceinc,sourcegeo,geo,mapping,colmapping)
            inc.write(inconfilename)

    def convert_mulkom_heat_conductivity(self):
        """Converts MULKOM-style rock heat conductivities to TOUGH2-style.  MULKOM used a formulation based on wet conductivity
        and liquid conductivity, weighted by porosity (though the liquid value appears to have always been zero).  TOUGH2 uses
        a formulation based on wet and dry conductivities, weighted by liquid saturation (though by default the wet and dry
        values are the same).  Here we scale the rock heat conductivities to give the same effective values as the
        MULKOM formulation."""
        for rt in self.grid.rocktypelist: rt.conductivity*=(1.-rt.porosity)

    def convert_to_TOUGH2(self,warn=True,MP=False):
        """Converts an AUTOUGH2 data file to a TOUGH2 data file.  Various MOP parameters are changed to try to make them
        the TOUGH2 simulation give similar results to AUTOUGH2 where possible.  AUTOUGH2-specific input blocks are removed
        and some generator types changed.  If MP is True, the conversion is done to a TOUGH2_MP data file, which treats a 
        few of the parameters differently."""
        if MP: self.filename='INFILE'
        # remove AUTOUGH2-specific input blocks:
        self.simulator=''
        if self.multi:
            if 'eos' in self.multi: del self.multi['eos']
            self.multi['num_inc']=None
        self.lineq={}
        self.short_output={}
        # convert parameters:
        warnings=[]
        if self.parameter['option'][10]==2:
            self.parameter['option'][10]=0
            self.convert_mulkom_heat_conductivity()
            warnings.append('MOP(10)=2: MULKOM rock heat conductivities (values have been converted to TOUGH2 equivalents)')
        if self.parameter['option'][12]==2:
            self.parameter['option'][12]=0
            warnings.append('MOP(12)=2: piecewise linear well table interpolation')
        if self.parameter['option'][21]>0: self.parameter['option'][21]=0 # not used in AUTOUGH2, but used in TOUGH2
        if self.parameter['option'][22]>0:
            self.parameter['option'][22]=0
            warnings.append('MOP(22)>0: USERBC')
        if self.parameter['option'][23]>0:
            isat2=self.simulator.startswith('AUTOUGH2') and (not self.simulator.startswith('AUTOUGH2.2'))
            ismulkom=self.simulator.startswith('MULKOM')
            mulkom_compatibility=self.parameter['option'][23] in [0,1]
            if (isat2 or ismulkom) and mulkom_compatibility:
                self.convert_mulkom_heat_conductivity()
                warnings.append('MOP(23)>0: MULKOM/TOUGH2 backward compatibility')
            self.parameter['option'][23]=0
        if self.parameter['option'][24]>0:
            self.parameter['option'][24]=0
            warnings.append('MOP(24)>0: initial printout of tables')
        if MP: # these MOPs mean different things in TOUGH2_MP:
            if self.parameter['option'][14]>0:
                self.parameter['option'][14]=0
                warnings.append('MOP(14)>0: Pivot failure handling')
            if self.parameter['option'][17]>0:
                self.parameter['option'][17]=0
                warnings.append('MOP(17)>0: Jacobian scaling')
            if self.parameter['option'][20]>0:
                self.parameter['option'][20]=0
                warnings.append('MOP(20)>0: Disabling vapour pressure lowering')
        if warn and len(warnings)>0:
            print 'The following options are not supported in TOUGH2:'
            for warning in warnings: print warning
        # convert generator types or delete any that can't be converted:
        allowed=['HEAT','WATE','AIR ','MASS','DELV']
        convert={'CO2 ':'COM2'}
        delgens=[]
        for gen in self.generatorlist:
            if gen.type in convert.keys(): gen.type=convert[gen.type]
            elif not ((gen.type in allowed) or gen.type.startswith('COM')): delgens.append((gen.block,gen.name))
        if warn and len(delgens)>0:
            print 'The following generators have types not supported by TOUGH2 and have been deleted:'
            print delgens
