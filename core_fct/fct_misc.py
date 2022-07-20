"""
Copyright: IIASA (International Institute for Applied Systems Analysis), 2016-2021; CEA (Commissariat a L'Energie Atomique) & UVSQ (Universite de Versailles et Saint-Quentin), 2016
Contributor(s): Thomas Gasser (gasser@iiasa.ac.at), Yann Quilcaille

This software is a computer program whose purpose is to simulate the behavior of the Earth system, with a specific but not exclusive focus on anthropogenic climate change.

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 

As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. 

In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and,  more generally, to use and operate it in the same conditions as regards security. 

The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
"""

##################################################
##################################################

import os
import csv
import warnings
import numpy as np
import xarray as xr
import scipy.odr as odr
import matplotlib.pyplot as plt


##################################################
##   A. ODR FIT ROUTINES
##################################################

## rolling mean function
def rolling_mean(da, time_axis, time_len):
    return da.rolling({time_axis: time_len}, center=True).mean().dropna(time_axis).values


## rolling std function
def rolling_std(da, time_axis, time_len, min_cv):
    min_std = min_cv * rolling_mean(da, time_axis=time_axis, time_len=time_len)
    std = da.rolling({time_axis: time_len}, center=True).std().dropna(time_axis).values    
    return np.maximum(std, min_std)


## wrapping function for ODR fit
def fit_odr(data_x, data_y, par_fg, par_lb, par_ub, y_func, y_jacb=None, y_jacx=None, 
    std_x=True, std_y=True, time_axis='year', time_len=5, min_cv=1E-6, scen_axis='scen', scen_skip=[]):
    '''
    Wrapper to execute an ODR fit on provided data. See scipy.odr for full documentation.
    
    Input:
    ------
    data_x (list)           input predictor variables as list of xr.DataArray
    data_y (xr.DataArray)   input predicted variable as single data array (can be None for implicit fit)
    par_fg (list)           list of first guess parameters values (len = p)
    par_lb (list)           list of lower bound parameters values (len = p)
    par_ub (list)           list of upper bound parameters values (len = p)
    y_func (callable)       function of the type y = f(par, x)
        
    Output:
    -------
    par_out (list)          output parameters values

    Options:
    --------
    y_jacb (callable)       function for the Jacobian of y_func on the par axis;
                            default = None
    y_jacx (callable)       function for the Jacobian of y_func on the x axis;
                            default = None
    std_x (bool)            whether uncertainty in data_x should be used as weight in the fit;
                            default = True
    std_y (bool)            whether uncertainty in data_y should be used as weight in the fit;
                            default = True
    time_axis (str)         name of time axis over which the rolling values are taken;
                            default = 'year'
    time_len (int)          length of the period over which rolling values are calculated;
                            default = 5
    min_cv (float)          value of minimum relative uncertainty when calculating rolling std;
                            default = 1E-6
    scen_axis (str)         name of scenario/experiment axis over which the data is flattened;
                            default = 'scen'
    scen_skip (list)        list of scenarios/experiments to be ignored by the fit;
                            default = []
    '''
    assert len(par_fg) == len(par_lb) == len(par_ub)

    ## conversions between beta and parameters
    p_ = lambda b, pmin, pmax: pmin + (pmax - pmin) / (1 + np.exp(-b))
    b_ = lambda p, pmin, pmax: np.log((p - pmin) / (pmax - p))
    dp_db = lambda b, pmin, pmax: (pmax - pmin) * np.exp(-b) / (1 + np.exp(-b))**2

    ## shortcut functions
    F = lambda f_, alpha: [pmin if pmin==pmax else f_(a, pmin, pmax) for a, pmin, pmax in zip(alpha, par_lb, par_ub)]
    mean = lambda da: rolling_mean(da, time_axis=time_axis, time_len=time_len)
    std = lambda da: rolling_std(da, time_axis=time_axis, time_len=time_len, min_cv=min_cv)

    ## format data
    x = [np.hstack([mean(da.sel({scen_axis: scen}, drop=True)) for scen in da.scen if scen not in scen_skip]) for da in data_x]
    if len(x) == 1: x = [x[0], 0 * x[0]] # weird fix because of odr
    y = np.hstack([mean(data_y.sel({scen_axis: scen}, drop=True)) for scen in data_y.scen if scen not in scen_skip]) if data_y is not None else True
    sx = [np.hstack([std(da.sel({scen_axis: scen}, drop=True)) for scen in da.scen if scen not in scen_skip]) for da in data_x] if std_x else None
    if len(sx) == 1: sx = [sx[0], min_cv * x[0]] # weird fix because of odr
    sy = np.hstack([std(data_y.sel({scen_axis: scen}, drop=True)) for scen in data_y.scen if scen not in scen_skip]) if std_x else None
    data = odr.RealData(x=x, y=y, sx=sx, sy=sy)

    ## define model
    func = lambda beta, x: y_func(F(p_, beta), x)
    jacb = (lambda beta, x: y_jacb(F(p_, beta), x) * np.array(F( dp_db, beta))[:, np.newaxis]) if y_jacb is not None else None
    jacx = (lambda beta, x: y_jacx(F(p_, beta), x)) if y_jacx is not None else None
    model = odr.Model(func, jacb, jacx, implicit=True if data_y is None else False)

    ## ODR fit and return
    odr_out = odr.ODR(data, model, beta0=F(b_, par_fg)).run()
    return F(p_, odr_out.beta)
    

## checking and plotting fit
def check_odr(data_x, data_y, par_bg, y_func,
    time_axis='year', time_len=5, min_cv=1E-6, scen_axis='scen', scen_skip=[], plot_ax=None):

    ## get data
    x_ref = [da.rolling({time_axis: time_len}, center=True).mean() for da in data_x]
    y_ref = data_y.rolling({time_axis: time_len}, center=True).mean() if data_y is not None else 0.
    y_mod = y_func(par_bg, x_ref)

    ## function for select in or out data
    get_in = lambda da: da.sel({scen_axis: [scen for scen in data_x[0].scen.values if scen not in scen_skip]})
    get_out = lambda da: da.sel({scen_axis: [scen for scen in data_x[0].scen.values if scen in scen_skip]})

    ## calculate metrics
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        MSE_in = ((get_in(y_ref) - get_in(y_mod))**2).mean().values
        MSE_out = ((get_out(y_ref) - get_out(y_mod))**2).mean().values
        R2_in = 1 - MSE_in / ((get_in(y_ref) - get_in(y_ref).mean().values)**2).mean().values
        R2_out = 1 - MSE_out / ((get_out(y_ref) - get_out(y_ref).mean().values)**2).mean().values
        NRMSE_in = np.sqrt(MSE_in) / get_in(y_ref).mean().values
        NRMSE_out = np.sqrt(MSE_out) / get_out(y_ref).mean().values

    ## optional plot (to be refined!)
    if plot_ax is not None:
        plt.style.use('seaborn-colorblind')
        ## plot data
        plt.plot(x_ref[0].year, y_ref, marker='+', ls='none', ms=4, alpha=0.8)
        ## plot fit
        for scen in x_ref[0][scen_axis]:
            not_null = y_ref.sel({scen_axis: scen}).notnull().sum() if data_y is not None else x_ref[-1].sel({scen_axis: scen}).notnull().sum()
            if not_null > 0:
                plt.plot(x_ref[0].year, y_mod.sel({scen_axis: scen}), color='k', lw= 1.5 if scen in scen_skip else 1, ls=':' if scen in scen_skip else '-')
        ## axis
        plt.xticks(fontsize='x-small')
        plt.yticks(fontsize='x-small')
        ## legend
        plt.xlabel(', '.join(['p{0} = {1:.3f}'.format(p, par) for p, par in enumerate(par_bg)]), fontsize='x-small')
        plot_ax.xaxis.set_label_position('top')
        ## legend (metrics)
        for name, val in zip(['NRMSE(in)', 'NRMSE(out)', 'R2(in)', 'R2(out)'], [NRMSE_in, NRMSE_out, R2_in, R2_out]): 
            plt.plot([np.nan], [np.nan], label='{0} = {1:.3f}'.format(name, val))
        plt.legend(loc=0, ncol=1, handlelength=0, handletextpad=0, borderpad=0.2, labelspacing=0.2, frameon=False, fontsize='xx-small')

    ## return metrics
    return {name: val for name, val in zip(['NRMSE(in)', 'NRMSE(out)', 'R2(in)', 'R2(out)'], [NRMSE_in, NRMSE_out, R2_in, R2_out])}

    
##################################################
##   B. PROBABILITY DISTRIBUTIONS
##################################################

## function to get lognorm distrib parameters
def lognorm_distrib_param(mean, std):
    mu = np.log(mean / np.sqrt(1. + std**2./mean**2.))
    sigma = np.sqrt(np.log(1. + std**2./mean**2.))
    return mu, sigma


##################################################
##   C. REGIONAL AGGREGATION
##################################################

## aggregate regional data to OSCAR regions
## /!\ WARNING: translated regions that are not in input data will not appear in output data (i.e. might have to use combine_first instead of merge, with the output data)
def aggreg_region(ds_in, mod_region, weight_var={}, old_axis='reg_iso', new_axis='reg_land', time_axis='year', dataset=None):
    '''
    Function to aggregate data onto OSCAR regions. It uses dictionnaries mapping ISO regions to OSCAR regions defined in 'input_data/regions' by user.
    
    Input:
    ------
    ds_in (xr.Dataset)  input dataset to be aggregated
    mod_region (str)    name of regional aggregation (must a valid option)
        
    Output:
    -------
    ds_out (xr.Dataset) output dataset

    Options:
    --------
    weight_var (dict)   keys variables are weighted using values variables when aggregating; 
                        this is necessary for intensive variables (e.g. temperature that needs to be weighted by area); 
                        keys and values are names (str) of ds_in variables;
                        default = {}
    old_axis (str)      name of regional axis that will be aggregated (must a dim of ds_in);
                        default = 'reg_iso'
    new_axis (str)      name of new aggregated regional axis (must NOT be in ds_in, and will be in ds_out);
                        default = 'reg_land'
    time_axis (str)     name of time axis (to ensure it is first dim in ds_out);
                        default = 'year'
    dataset (str)       name of dataset of original aggregation if not ISO;
                        this likely requires changing old_axis as well;
                        see 'input_data/regions/dico_other_datasets.csv' for examples;
                        default = None
    '''
    
    ## check old axis in ds_in and new_axis not in ds_in
    assert old_axis in ds_in.coords and new_axis not in ds_in.coords
    
    ## make deep copy to be safe
    ds_out = ds_in.copy(deep=True)

    ## region mapping files to be loaded
    list_load = [zou for zou in os.listdir('input_data/regions/') if all([_  in zou for _ in ['dico_', '.csv']])]

    ## load and create combined dictionary
    dico = {}
    for zou in list_load:
        with open('input_data/regions/' + zou) as f: TMP = np.array([line for line in csv.reader(f)])
        ## ISO-formated regions
        if 'ISO' in zou: 
            dico = {**dico, **{int(key):int(val) for key, val in zip(TMP[1:,0], TMP[1:,TMP[0,:].tolist().index(mod_region)])}}
        ## dataset-specific regions
        else: 
            dico = {**dico, **{(key0, key1):int(val) for key0, key1, val in zip(TMP[1:,0], TMP[1:,1], TMP[1:,TMP[0,:].tolist().index(mod_region)])}}

    ## load long region names
    with open('input_data/regions/regions_long_name.csv') as f: 
        TMP = np.array([line for line in csv.reader(f)])
    long_name = {n:name for n, name in enumerate(TMP[1:,TMP[0,:].tolist().index(mod_region)])}

    ## apply weights to weighted variables
    for var in weight_var:
        ds_out[var] = ds_out[var] * ds_out[weight_var[var]]

    ## extract variables without regional axis
    ds_non = ds_out.drop([var for var in ds_out if old_axis in ds_out[var].dims] + [old_axis])
    ds_out = ds_out.drop([var for var in ds_out if old_axis not in ds_out[var].dims])

    ## new regional aggregation
    if dataset is None: 
        ds_out.coords[new_axis] = xr.DataArray([dico[reg] for reg in ds_out[old_axis].values], dims=old_axis)
    else: 
        ds_out.coords[new_axis] = xr.DataArray([dico[(dataset, reg)] for reg in ds_out[old_axis].values], dims=old_axis)
    ds_out = ds_out.groupby(new_axis).sum(old_axis, min_count=1, keep_attrs=True)
    ds_out.coords[new_axis + '_long_name'] = xr.DataArray([long_name[reg] for reg in ds_out[new_axis].values], dims=new_axis)

    ## remove weights
    for var in weight_var:
        ds_out[var] = ds_out[var] / ds_out[weight_var[var]]

    ## merge with extracted variables
    ds_out = xr.merge([ds_out, ds_non])

    ## make sure time axis is first
    if time_axis in ds_out.coords: 
        for var in ds_out:
            if time_axis in ds_out[var].coords:
                ds_out[var] = 0.*ds_out[time_axis] + ds_out[var]
    
    ## return
    return ds_out


##################################################
##   D. GENERAL LOADING
##################################################

## general loading function
def load_data(name):
    '''
    Convenience function to load any dataset in the 'input_data' folder.
    
    Input:
    ------
    name (str)      part of file name -- must be unique to work
        
    Output:
    -------
    (xr.Dataset)    loaded dataset
    '''

    ## lists all available files
    list_files = [os.path.join(dp, f) for dp, dn, fn in os.walk('input_data/') for f in fn if f[-3:] == '.nc']
    
    ## check compatible files
    okay_files = [f for f in list_files if name in f]

    ## error if no file or not unique
    if len(okay_files) == 0:
        raise RuntimeError("no files were found (input name = '{0}')".format(name))
    elif len(okay_files) > 1:
        raise RuntimeError("more than one file was found (input name = '{0}')".format(name))

    ## return loaded file otherwise
    with xr.open_dataset(okay_files[0]) as TMP: 
        return TMP.load()


##################################################
##   E. EXTEND TIME-SERIES
##################################################

## extend one timeseries following another timeseries
def extend_timeseries(ref, ext, direction, time_axis='year', juxtaposition=False, ref_length=1, scale_to_global=False, dump_loss_in=None, stable_axis=[]):
    '''
    Function to extend one time-series with another time-series. It takes the former time-series as-is, and extends it following the latter's relative variations.
    
    Input:
    ------
    ref (xr.DataArray)      initial data array to be extended
    ext (xr.DataArray)      second data array to be used for extension
    direction (str)         must be 'forward' or 'backward' (in time)
        
    Output:
    -------
    new (xr.DataArray)      extented data array

    Options:
    --------
    time_axis (str)         name of time axis;
                            default = 'year'
    juxtaposition (bool)    option to force simple juxtaposition of ref and ext (instead of extending following relative changes);
                            default = False
    ref_length (int)        number of time-steps taken for reference to scale ext to ref;
                            ref and ext must have at least ref_length time-steps in common right before time of junction;
                            default = 1
    scale_to_global (bool)  if False, the extension is done on disaggregated (e.g. regional) basis;
                            if True, the extension will be done over the aggregated data arrays;
                            default = False
    dump_loss_in (dict)     dictionnary of ONE coordinate in which the difference between scale_to_global True and False will be dumped;
                            this ensure that the sum will meet scale_to_global and the disaggregated values be continuous (except for the one used as dump);
                            scale_to_global must also be True for this option to work;
                            default = None
    stable_axis (list)      list of dims of ref/ext to be ignored when scale_to_global is True;
                            default = []
    '''

    ## define summed axis
    summed_axis = [var for var in ref.dims if var not in stable_axis and var != time_axis]

    ## pruning time axis and getting overlapping years
    ref = ref.dropna(time_axis, how='all')
    ext = ext.dropna(time_axis, how='all')
    time_overlap = ref[time_axis] & ext[time_axis]
    
    ## new years to cover
    if direction == 'forward': 
        time_extend = np.arange(int(time_overlap.max())+1, int(ext[time_axis].max())+1)
        time_ref = int(time_overlap.max())+1 + np.arange(-ref_length, 0)
    elif direction == 'backward': 
        time_extend = np.arange(int(ext[time_axis].min()), int(time_overlap.min()))
        time_ref = int(time_overlap.min()) + np.arange(0, +ref_length)
    else: 
        raise ValueError("'direction' can only be 'forward' or 'backward'")
    time_extend = xr.DataArray(time_extend, coords={time_axis:time_extend}, dims=[time_axis], name=time_axis)
    time_ref = xr.DataArray(time_ref, coords={time_axis:time_ref}, dims=[time_axis], name=time_axis)

    ## new dataset with extended time axis
    new = ref.dropna(time_axis, how='all')
    new = ref.combine_first(np.nan * time_extend)

    ## add other stable axis
    for var in stable_axis:
        if var in ext.coords:
            new = new + xr.zeros_like(ext[var], dtype=float)
    
    ## check feasibility of dumping scaling losses
    dump_available = False
    if not juxtaposition and scale_to_global and dump_loss_in is not None:
        if set(dump_loss_in.keys()) == set(summed_axis):
            if all([val in new[key] for key, val in dump_loss_in.items()]):
                dump_available = True
        ## error if requested but unavailable
        if not dump_available:
            raise ValueError("'dump_loss_in' requested but inadequate set of coordinates provided: {0}".format(dump_loss_in))

    ## extend following relative trends
    if not juxtaposition:
        ## catch warnings (normally caused by averaging on ref_length=1)
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')

            ## scale global data and relative contributions, and dump losses in specified address
            if scale_to_global and dump_available:
                new.loc[time_extend] = new.loc[time_extend].fillna(0) + ext.loc[time_extend] * ref.loc[time_ref].mean(time_axis) / ext.loc[time_ref].mean(time_axis)
                glob = ext.loc[time_extend] * ref.loc[time_ref].mean(time_axis).sum(summed_axis) / ext.loc[time_ref].mean(time_axis).sum(summed_axis)
                new.loc[dump_loss_in].loc[time_extend] = new.loc[time_extend].loc[dump_loss_in].fillna(0) + glob.sum(summed_axis) - new.loc[time_extend].sum(summed_axis)
                del glob

            ## scale global data but relative contributions are not affected (and will likely be discountinuous)
            elif scale_to_global:
                new.loc[time_extend] = new.loc[time_extend].fillna(0) + ext.loc[time_extend] * ref.loc[time_ref].mean(time_axis).sum(summed_axis) / ext.loc[time_ref].mean(time_axis).sum(summed_axis)
            
            ## scale only local data (i.e. absolute contributions will be continuous, but global value will vary with aggregation)
            else:
                new.loc[time_extend] = new.loc[time_extend].fillna(0) + ext.loc[time_extend] * ref.loc[time_ref].mean(time_axis) / ext.loc[time_ref].mean(time_axis)

    ## or simply juxtapose both timeseries
    else: 
        new.loc[time_extend] = new.loc[time_extend].fillna(0) + ext.loc[time_extend]

    ## return
    return new


##################################################
##   F. GROUP CONSISTENT SCENARIOS
##################################################

## group scenarios on one common axis
def group_scenarios(ds_in, scen_axis='scen', time_axis='year', group_scen=[], default_scen=[]):
    '''
    Function to group several scenario dims of a dataset onto one unique dim. By default, only scenarios common to all variables are taken.
    
    Input:
    ------
    ds_in (xr.Dataset)      dataset in which each variable may have its own scenario dim
        
    Output:
    -------
    ds_out (xr.Dataset)     dataset with united scenario dim

    Options:
    --------
    scen_axis (str)         name of scenario axis that will be used as output dim name;
                            any dim of ds_in containing this str will be considered a scenario dim;
                            default = 'scen'
    time_axis (str)         name of time axis (to clean-up ds_out)
                            default = 'year'
    group_scen (list)       list of scenarios that will be in ds_out, even if not all variables have them;
                            if empty, output scenarios are only those in common in all variables of ds_in;
                            default = []
    default_scen (list)     if a variable does not have a scenario of group_scen, default_scen are used instead (if possible);
                            taken in priority following the list order;
                            default = []
    '''

    ## list of scenario axes
    list_axis = [var for var in ds_in.coords if scen_axis+'_' in var]

    ## get common scenarios + required scenarios
    final_scen = set(ds_in[list_axis[0]].values)
    for axis in list_axis:
        final_scen = final_scen & set(ds_in[axis].values)
    final_scen = final_scen | set(group_scen)
   
    ## group scenarios on one axis
    ds_out = xr.Dataset()
    if final_scen != set():
        for VAR in ds_in:

            ## get scen axis
            axis = [var for var in ds_in[VAR].coords if scen_axis+'_' in var][0]
            
            ## if all scenarios are available
            if final_scen.issubset(set(ds_in[axis].values)):
                ds_out[VAR] = ds_in[VAR].sel(**{axis:list(final_scen)}).rename({axis:scen_axis})

            ## otherwise apply default scenario
            elif any([var in ds_in[axis] for var in default_scen]):
                default = [var for var in default_scen if var in ds_in[axis]][0]
                TMP = []
                for scen in final_scen:
                    if scen in ds_in[axis]:
                        TMP.append(ds_in[VAR].sel(**{axis:scen}, drop=True).assign_coords(**{scen_axis:scen}).expand_dims(scen_axis, -1))
                    else:
                        TMP.append(ds_in[VAR].sel(**{axis:default}, drop=True).assign_coords(**{scen_axis:scen}).expand_dims(scen_axis, -1))
                ds_out[VAR] = xr.concat(TMP, dim=scen_axis)
    
            ## return NaN if impossible
            else:
                ds_out[VAR] = xr.concat([np.nan * ds_in[VAR].sum(axis).assign_coords(**{scen_axis:scen}).expand_dims(scen_axis, -1) for scen in final_scen], dim=scen_axis)

    ## return (cleaned dataset)
    return ds_out.dropna(time_axis, how='all')

