#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from os.path import basename
import numpy as np
import Sphere as sp

###########
# constants

# quality tolerance for distance from Euler pole
tol = 1e-6

# earth radius
rn = 6371.0

######
# code

def circle(plon, plat, lon, lat, step=100.0):
    """Returns evenly spaced sampling points along a small circle, 
    specified by an Euler pole and an starting "anchor" point.
    Arguments:
        plon : pole lon (degrees)
        plat : pole lat (degrees)
        lon  : anchor lon (degrees)
        lat  : anchor lat (degrees)
        step : approximate step size along small circle (km)
               optional, default = 100 km
    Returns:
        lons : small-circle lons (degrees)
        lats : small-circle lats (degrees)
    """
    # convert to geocentric
    plat = 90 - sp.geocentric(90 - plat)
    lat  = 90 - sp.geocentric(90 - lat)
    # pole-anchor distance, azimuth
    delta, az0 = sp.delaz((plon, plat), (lon, lat))
    # determine (approximate) step size
    dx = np.degrees(step / (np.sin(np.radians(delta)) * rn))
    # check is pole is, um, polar
    if (90 - np.abs(plat)) / 90 < tol:
        # if so, just rotate around at the same lat, with the direction
        # depending on the hemisphere of the pole (s.t. the result is still
        # right-handed)
        if plat > 0:
            wplons = lon + np.arange(0, 360, dx)
        else:
            wplons = lon - np.arange(0, 360, dx)
        wplons[wplons > 360] -= 360
        wplons[wplons < 0]   += 360
        waypoints = np.array([(wplon,lat) for wplon in wplons])
    else:
        # determine azimuths to waypoints (right handed)
        azs = az0 - np.arange(0, 360, dx)
        azs[azs > 360] -= 360
        azs[azs < 0]   += 360
        # find waypoints
        waypoints = np.array([sp.shoot((plon,plat), delta, a) for a in azs])
    # check pole-waypoint distance
    deltas = np.array([
        sp.delaz((plon,plat), (waypoints[i,0],waypoints[i,1]),
            delta_only=True)
        for i in range(waypoints.shape[0])])
    assert np.abs((deltas - delta) / delta).max() < tol
    # to geodetic, return
    lons = waypoints[:,0]
    lats = 90 - sp.geodetic(90 - waypoints[:,1])
    return lons, lats

def circle_dist(plon, plat, dist, step=100.0):
    """Returns evenly spaced sampling points along a small circle, 
    specified by an Euler pole and an angular distance from it
    Arguments:
        plon : pole lon (degrees)
        plat : pole lat (degrees)
        dist : distance from the euler pole (degrees)
        step : approximate step size along small circle (km)
               optional, default = 100 km
    Returns:
        lons : small-circle lons (degrees)
        lats : small-circle lats (degrees)
    """
    # convert to geocentric
    plat = 90 - sp.geocentric(90 - plat)
    # determine (approximate) step size
    dx = np.degrees(step / (np.sin(np.radians(dist)) * rn))
    # check if pole is, um, polar
    if (90 - np.abs(plat)) / 90 < tol:
        # if so, just rotate around at the same lat, with the direction
        # depending on the hemisphere of the pole (s.t. the result is still
        # right-handed)
        if plat > 0:
            wplons = np.arange(0, 360, dx)
            wplat  = 90 - dist
        else:
            wplons = - np.arange(0, 360, dx)
            wplat  = - 90 + dist
        wplons[wplons > 360] -= 360
        wplons[wplons < 0]   += 360
        waypoints = np.array([(wplon,wplat) for wplon in wplons])
    else:
        # determine azimuths to waypoints (right handed)
        azs = np.arange(0, 360, dx)
        # find waypoints
        waypoints = np.array([sp.shoot((plon,plat), dist, a) for a in azs])
    # check pole-waypoint distance
    deltas = np.array([
        sp.delaz((plon,plat), (waypoints[i,0],waypoints[i,1]),
            delta_only=True)
        for i in range(waypoints.shape[0])])
    assert np.abs((deltas - dist) / dist).max() < tol
    # to geodetic, return
    lons = waypoints[:,0]
    lats = 90 - sp.geodetic(90 - waypoints[:,1])
    return lons, lats


######################
# command-line wrapper

def main():
    ap = ArgumentParser(
        description="""
        Calculate sampling along a small circle, given an Euler pole and an anchor point.
        Result is printed to stdout in GMT format.""",
        usage="""%(prog)s -p lon lat -a lon lat [-s step]

        Example: %(prog)s -p 96.2 -64.1 -a -110.0 -10.0 -s 1000.0""")
    ap.add_argument('-p', '--pole', required=True, type=float, nargs=2, help='Euler pole: lon lat')
    ap.add_argument('-a', '--anchor', required=True, type=float, nargs=2, help='Anchor point: lon lat')
    ap.add_argument('-s', '--step', default=100.0, type=float, help='Step size (approximate); default: 100 km')
    args = ap.parse_args()
    plon, plat = args.pole
    alon, alat = args.anchor
    print('> pole: %.1f %.1f / anchor %.1f %.1f' % (
        plon, plat, alon, alat))
    lons, lats = circle(plon, plat, alon, alat, step=args.step)
    np.savetxt(sys.stdout, np.vstack((lons, lats)).T, fmt='%10.3f %10.3f')

if __name__ == '__main__':
    main()
