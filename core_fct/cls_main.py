"""
Copyright: IIASA (International Institute for Applied Systems Analysis), 2016-2020; CEA (Commissariat a L'Energie Atomique) & UVSQ (Universite de Versailles et Saint-Quentin), 2014-2016
Contributor(s): Thomas Gasser (gasser@iiasa.ac.at), Yann Quilcaille

This software is a computer program whose purpose is to simulate the behavior of the Earth system, with a specific but not exclusive focus on anthropogenic climate change.

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 

As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. 

In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and,  more generally, to use and operate it in the same conditions as regards security. 

The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
"""

import warnings
import numpy as np
import xarray as xr

from time import perf_counter

from core_fct.fct_misc import Int_ExpInt as Int_dflt


##################################################
##   1. MODELS
##################################################

class Model():
    '''
    Class defining a Model, i.e. a collection of linked Process objects.

    Init:
    ------
    (nothing)

    Options:
    --------
    name (str)      name of the model;
                    default = ''
    '''

    ## ------
    ## Basics
    ## ------

    ## initialization
    def __init__(self, name=''):
        assert type(name) is str
        self.name = name
        self._processes = {}

    ## check if process in model
    def __contains__(self, item):
        assert type(item) == str
        return item in self.proc_all
    
    ## nice display
    def __repr__(self):
        out = '<' + str(type(self)).replace("<class '", "").replace("'>", "")  + " '" + self.name + "' " + '(Processes: ' + str(len(self)) + ')>' + '\n'
        for proc in self._processes.values():
            out += proc.Out + ' ' + str(proc.In).replace("'", "") + ': ' + proc.Eq.__repr__() + '\n'
        return out

    ## number of processes
    def __len__(self): return len(self._processes)

    ## may wanna define:
    ## __del__, __add__, __sub__, __iadd__, __isub__, __radd__, __rsub__

    ## ----------
    ## Properties
    ## ----------

    ## lists of variables
    @property
    def var_all(self): return set([proc.Out for proc in self._processes.values()]) | set([var for proc in self._processes.values() for var in proc.In])
    @property
    def var_mid(self): return set([proc.Out for proc in self._processes.values()]) & set([var for proc in self._processes.values() for var in proc.In])
    @property
    def var_out(self): return set([proc.Out for proc in self._processes.values()]) - self.var_mid
    @property
    def var_in(self): return set([var for proc in self._processes.values() for var in proc.In]) - self.var_mid
    @property
    def var_prog(self): return set([proc.Out for proc in self._processes.values() if proc.prog])
    @property
    def var_diag(self): return set([proc.Out for proc in self._processes.values() if not proc.prog])
    @property
    def proc_all(self): return list(self._processes.keys())

    ## get causality tree levels
    def proc_levels(self, var_node=[]):
        ## initialize level 0 with state variables
        levels = {0:list(self.var_prog) + var_node}
        proc_remain = set(self.proc_all) - self.var_prog - set(var_node)
        ## loop through levels
        while proc_remain != set():
            next_level = []
            for proc in proc_remain:
                ## check whether solvable with lower-level variables
                if all([var in sum([val for val in levels.values()], list(self.var_in)) for var in self._processes[proc].In]):
                    next_level.append(proc)
            levels[max(levels.keys())+1] = next_level
            proc_remain = proc_remain - set(next_level)
            ## break if no variables in this level
            if next_level == []:
                levels[np.inf] = list(proc_remain)
                break
        ## return dic of variables by level
        return levels

    ## ----------
    ## Operations
    ## ----------

    ## copy the model to another one
    def copy(self, add_name='_copy', new_name=None, only=None):
        ## name new model
        if new_name is None: new_model = Model(self.name + add_name)
        else: new_model = Model(new_name)
        ## copy processes
        for proc in self._processes.values():
            if only is None or proc.Out in only:
                new_model.process(proc.Out, proc.In, proc.Eq, units=proc.units, core_dims=proc.core_dims)
        ## return new model
        return new_model
    
    ## check whether same model
    ## /!\ TODO
    def issame(self, other):
        assert type(other) == Model
        if not self.proc_all == other.proc_all: return False
        for key in self.proc_all:
            if self._processes[key].In != other._processes[key].In: return False
            if self._processes[key].Out != other._processes[key].Out: return False
            if self._processes[key].Eq != other._processes[key].Eq: return False
        raise Warning('TODO!')

    ## check whether submodel (excl. same model)
    ## /!\ TODO
    def issubmodel(self, other):
        assert type(other) == Model
        if self.issame(other): return False
        raise Warning('TODO!')

    ## shortcut for comparison (note: non-equal is not logical negation of equal)
    def __eq__(self, other): return self.issame(other)
    def __ne__(self, other): return not self.issame(other) and not self.issubmodel(other) and not other.issubmodel(self)
    def __le__(self, other): return self.issubmodel(other) or self.issame(other)
    def __ge__(self, other): return other.issubmodel(self) or self.issame(other)
    def __lt__(self, other): return self.issubmodel(other)
    def __gt__(self, other): return other.issubmodel(self)

    ## --------
    ## Defining
    ## --------

    ## define process
    def process(self, Out, In, Eq, **proc_args):
        self._processes[Out] = Process(Out, In, Eq, self, **proc_args)
        return self._processes[Out]

    ## get process
    def __getitem__(self, key):
        return self._processes[key]

    ## set process
    def __setitem__(self, key, val):
        assert type(key) == str and type(val) == Process
        if key == val.Out: self._processes[key] = val
        else: raise KeyError('key ({0}) and process name ({1}) must be the same'.format(key, val.Out))

    ## delete process
    def __delitem__(self, key):
        del self._processes[key]

    ## --------
    ## Plotting
    ## --------

    ## display a graph of the model
    def display(self, random=True):
        ## check dependencies
        try:
            import networkx as nx
            import matplotlib.pyplot as plt
        except:
            raise ImportError("displaying the model requires 'matplotlib' and 'networkx' libraries installed")
        ## create graph
        DG = nx.DiGraph()
        for proc in self._processes.values():
            DG.add_edges_from([(var, proc.Out) for var in proc.In])
        ## plot
        plt.figure()
        if random: layout = nx.kamada_kawai_layout(DG, pos=nx.random_layout(DG))
        else: layout = nx.kamada_kawai_layout(DG) # nx.spring_layout(DG, k=0.8, iterations=500)
        nx.draw(DG, pos=layout, with_labels=True, alpha=0.8, font_size=8, node_size=2, edge_color='0.8')
        plt.title(self.name, fontsize='small')
        plt.show()

    ## -------
    ## Running
    ## -------

    ## check solvability by iteration (i.e. no loops in diagnostic variables)
    def _check_solvable(self, var_node=[]):
        assert all(var in self.var_all for var in var_node)
        levels = self.proc_levels(var_node)
        if np.inf in levels.keys(): 
            raise RuntimeError("infinite loop to solve diagnostic variables! set at least one of the following as (fake) prognostic variable: {0}".format([var.replace("'", "") for var in levels[np.inf]]))

    ## check initialized variables
    def _check_Ini(self, Ini):
        var_miss = (self.var_prog) - set(Ini.keys())
        if var_miss != set(): raise RuntimeError('missing initialized variables: {0}'.format(str(var_miss).replace("'", "")))

    ## check forcing variables and time axis
    def _check_For(self, For, time_axis):
        var_miss = self.var_in - set(For.keys())
        if var_miss != set(): raise RuntimeError('missing forcing variables: {0}'.format(str(var_miss).replace("'", "")))
        if time_axis not in For.coords: raise RuntimeError("forcing dataset has no time axis: '{0}'".format(time_axis))

    ## running model
    def __call__(self, Ini, Par, For, dtype=np.float32, var_keep=[], keep_prog=True, time_axis='year', Int=Int_dflt, nt=2, no_warnings=True):
        '''
        Input:
        ------
        Ini (xr.Dataset)        initial conditions
        Par (xr.Dataset)        parameters
        For (xr.Dataset)        forcing data

        Output:
        ------
        Var_out (xr.Dataset)    model outputs

        Options:
        --------
        dtype (type)            data type for computation, either np.float32 or float;
                                default = np.float32
        var_keep (list)         variables to be kept as output;
                                default = []
        keep_prog (bool)        whether prognostic and node variables should be kept as output;
                                default = True
        time_axis (str)         name of the time dimension;
                                default = 'year'
        Int (function)          integration function to solve differential equations (solving scheme);
                                default = exponential integrator ('Int_ExpInt' in 'fct_ancillary')
        nt (int)                number of substeps during the solving of the differential system;
                                default = 2
        '''

        ## various checks
        self._check_solvable()
        self._check_Ini(Ini)
        self._check_For(For, time_axis)
        
        ## force data type
        if dtype is not None:
            Ini = Ini.astype(dtype)
            Par = Par.astype(dtype)
            For = For.astype(dtype)

        ## printing time counter
        print(self.name + ' running')
        t0 = perf_counter()

        ## catch warnings (if requested)
        with warnings.catch_warnings():
            if no_warnings: warnings.filterwarnings('ignore')

            ## INITIALIZATION
            ## get time axis
            time = For.coords[time_axis]
            
            ## level in causality tree of variables (excl. non-kept metrics)
            levels = self.proc_levels()
            levels = {lvl:list(set(levels[lvl]) - (self.var_out - set(var_keep))) for lvl in levels}
            levels = {lvl:levels[lvl] for lvl in levels if len(levels[lvl]) > 0}

            ## initialization of all variables
            Var_old = Ini.copy(deep=True)
            for lvl in np.sort(list(levels.keys()))[1:]:
                for var in levels[lvl]:
                    Var_old[var] = self._processes[var](Var_old, Par, For.sel({time_axis:time[0]}, drop=True))

            ## initialization of kept variables
            Var_out = Var_old.drop([var for var in Var_old if var not in (list(self.var_prog)) * keep_prog + var_keep])
            Var_out = [Var_out.assign_coords(**{time_axis:time[0]}).expand_dims(time_axis, 0)]

            ## LOOP ON TIME-STEP
            for t in range(1, len(time)):
                print(time_axis + ' = ' + str(int(time[t])), end='\n' if t+1==len(time) else '\r')
                
                ## get time step and drivers
                dt = float(time[t] - time[t-1])
                For_t = For.sel({time_axis:time[t]}, drop=True)

                ## LOOP ON SUBSTEPS
                for _ in range(nt):
                    Var_new = xr.Dataset()

                    ## 1. prognostic variables
                    for var in self.var_prog:
                        Var_new[var] = self._processes[var](Var_old, Par, For_t, dt=dt/nt, Int=Int)
                    
                    ## 2. diagnostic variables
                    for lvl in np.sort(list(levels.keys()))[1:]:
                        for var in levels[lvl]:
                            Var_new[var] = self._processes[var](Var_new, Par, For_t)

                    ## finalize (iterate variables)
                    Var_old = Var_new.copy(deep=True)

                ## drop unwanted variables and append to final output
                Var_new = Var_new.drop([var for var in Var_new if var not in (list(self.var_prog)) * keep_prog + var_keep])
                Var_out.append(Var_new.assign_coords(**{time_axis:time[t]}).expand_dims(time_axis, 0))

            ## FINALIZATION
            ## concatenate final output on time axis
            Var_out = xr.concat(Var_out, dim=time_axis)

            ## add model info
            Var_out.attrs['model'] = self.name

            ## add units to output
            for var in Var_out: 
                Var_out[var].attrs['units'] = self._processes[var].units

        ## printing time counter
        print('total running time: {:.0f} seconds'.format(perf_counter() - t0))

        ## return
        return Var_out


##################################################
##   2. PROCESSES
##################################################

class Process():
    '''
    Class defining one single process of a Model object.

    Init:
    ------
    Out (str)           name of output variable (only one!)
    In (tuple)          names of input variables (as str); can be an empty tuple
    Eq (callable)       function linking In to Out

    Options:
    --------
    model (Model)       model to which the process belongs;
                        default = Model()
    units (str)         units of Out;
                        default = '?'
    core_dims (list)    dims over which Out must be defined (relevant only for prognostic variables);
                        default = []
    '''

    ## initialization
    def __init__(self, Out, In, Eq, model=Model(), units='?', core_dims=[]):
        assert type(Out) is str and type(In) is tuple and callable(Eq) and type(model) is Model and type(units) is str
        self.Out, self.In, self.Eq, self.model, self.units, self.core_dims = Out, In, Eq, model, units, core_dims
        if Out in In: self.prog = True
        else: self.prog = False

    ## nice display
    def __repr__(self):
        return self.Out + ' ' + str(self.In).replace("'", "") + ': ' + self.Eq.__repr__() + (" [in model: '"  + self.model.name + "']") * (self.model.name != '')

    ## -------

    ## get process drivers (for recursive call)
    def _get_var(self, var, Var, Par, For):
        ## if prescribed
        if var in For.keys(): return For[var]
        ## if available
        elif var in Var.keys(): return Var[var]
        ## otherwise
        elif var != self.Out: return self.model[var](Var, Par, For, recursive=True)
        else: raise RuntimeError('recursive call of process: {0}'.format(var))

    ## solve process
    def __call__(self, Var, Par, For=xr.Dataset(), recursive=False, time_axis='year', **Int_args):
        ## if prescribed
        if self.Out in For.keys(): return For[self.Out]
        ## otherwise get drivers
        if recursive: Var_in = xr.Dataset({var:self._get_var(var, Var, Par, For) for var in self.In})
        else: Var_in = xr.Dataset(dict([(var, For[var]) if var in For else (var, Var[var]) for var in self.In]))
        ## output (time axis first if requested and available)
        if time_axis not in Var.coords: return self.Eq(Var_in, Par, **Int_args)
        else: return 0.*Var[time_axis] + self.Eq(Var_in, Par, **Int_args)

