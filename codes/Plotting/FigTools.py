#!/usr/bin/env python

import os
import pickle

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from Sphere import gc_minor

#_TECTONICS_PATH = '$CODEPATH/Plotting/tectonics'
_TECTONICS_PATH = '/home/gcl/BR/barbara/data08/Plotting/tectonics'

def get_tectonics_path():
    return os.path.expandvars(_TECTONICS_PATH)

def plot_hotspots(m, lon360 = False, **kwargs):
    '''Plot hotspots from Steinberger (2000) on the supplied `Basemap` object.
    The `lon360` keyword argument, when set to True (default: False), will
    ensure all hotspots lons are in [0,360) before plotting. All other kwargs
    will be passed to the plot() method of the `Basemap` object.
    '''
    cache = os.path.expandvars('%s/hotspots.pkl' % (_TECTONICS_PATH))
    with open(cache, 'rb') as f:
        names, hotspots = pickle.load(f)
    if lon360:
        hotspots[:,0] = (hotspots[:,0] + 360) % 360.0
    x, y = m(hotspots[:,0], hotspots[:,1])
    if kwargs:
        m.plot(x, y, **kwargs)
    else:
        m.plot(x, y)
    return names, hotspots

def plot_plates(m, lon360 = False, jump_threshold = 1000000, **kwargs):
    '''Plot plate boundaries from the UTIG PLATES collection on the supplied
    `Basemap` object.  The `lon360` keyword argument, when set to True
    (default: False), will ensure all path lons are in [0,360) before plotting.
    All other kwargs will be passed to the plot() method of the `Basemap`
    object.
    '''
    for bound in ['ridge', 'transform', 'trench']:
        cache = os.path.expandvars('%s/%s.pkl' % (_TECTONICS_PATH, bound))
        with open(cache, 'rb') as f:
            name, segs = pickle.load(f)
        ind_nan, = np.nonzero(np.isnan(segs[:,0]))
        segs[ind_nan,0] = 0
        segs[ind_nan,1] = 0
        if lon360:
            segs[:,0] = (segs[:,0] + 360) % 360.0
        x, y = m(segs[:,0], segs[:,1])
        x[ind_nan] = np.nan
        y[ind_nan] = np.nan
        dx = np.abs(x[1:] - x[:-1])
        ind_jump, = np.nonzero(dx > jump_threshold)
        x[ind_jump] = np.nan
        if kwargs:
            m.plot(x, y, '--', **kwargs)
        else:
            m.plot(x, y, '--')

def plot_gc_clean(m, loc0, loc1, dx, **kwargs):
    '''Plot a great-circle path between `loc0` and `loc1` on the `Basemap`
    object `m` with the step size `dx`.  All other kwargs will be passed to the
    plot() method of the `Basemap` object.
    '''
    plot_path_clean(m, *gc_minor(loc0, loc1, dx), **kwargs)

def plot_path_clean(m, lons, lats, **kwargs):
    '''Plot a path, given by `lons` and `lats`, on the `Basemap` object `m`.
    All kwargs will be passed to the plot() method of the `Basemap` object.
    '''
    xgc, ygc = m(lons,lats)
    xgc_nd = np.array(xgc)
    ygc_nd = np.array(ygc)
    dx = np.abs(xgc_nd[1:] - xgc_nd[:-1])
    ind_flipx, = np.nonzero(dx > 1e6)
    xgc_nd[ind_flipx] = np.nan
    xgc = xgc_nd.tolist()
    if kwargs:
        m.plot(xgc, ygc, **kwargs)
    else:
        m.plot(xgc, ygc)
