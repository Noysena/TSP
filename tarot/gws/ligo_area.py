#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 16:34:55 2019

@author: noysena

Skymap area and distance
"""

"""
Obtain area in sq. deg form fits.gz
"""


from . import ArgumentParser, FileType, SQLiteType, figure_parser


def parser():
    parser = ArgumentParser(parents=[figure_parser])
    parser.add_argument(
        '--annotate', default=False, action='store_true',
        help='annotate plot with information about the event')
    parser.add_argument(
        '--contour', metavar='PERCENT', type=float, nargs='+',
        help='plot contour enclosing this percentage of'
        ' probability mass [may be specified multiple times, default: none]')
    parser.add_argument(
        '--colorbar', default=False, action='store_true',
        help='Show colorbar')
    parser.add_argument(
        '--radec', nargs=2, metavar='deg', type=float, action='append',
        default=[], help='right ascension (deg) and declination (deg) to mark')
    parser.add_argument(
        '--inj-database', metavar='FILE.sqlite', type=SQLiteType('r'),
        help='read injection positions from database')
    parser.add_argument(
        '--geo', action='store_true',
        help='Plot in geographic coordinates, (lat, lon) instead of (RA, Dec)')
    parser.add_argument(
        'input', metavar='INPUT.fits[.gz]', type=FileType('rb'),
        default='-', nargs='?', help='Input FITS file')
    return parser


def main(args=None):
    opts = parser().parse_args(args)
    
    # Late imports

    import numpy as np
#    import matplotlib.pyplot as plt
    import healpy as hp
    from ..io import fits
    from .. import postprocess
    from astropy.time import Time

    skymap, metadata = fits.read_sky_map(opts.input.name, nest=None)
    nside = hp.npix2nside(len(skymap))

    # Convert sky map from probability to probability per square degree.
    deg2perpix = hp.nside2pixarea(nside, degrees=True)
    
    if opts.geo:
        obstime = Time(metadata['gps_time'], format='gps').utc.isot
   #     ax = plt.axes(projection='geo degrees mollweide', obstime=obstime)
    #else:
    #    ax = plt.axes(projection='astro hours mollweide')
    #ax.grid()

    # Add contours.
    if opts.contour:
        cls = 100 * postprocess.find_greedy_credible_levels(skymap)
    radecs = opts.radec
    
    if opts.inj_database:
        query = '''SELECT DISTINCT longitude, latitude FROM sim_inspiral AS si
                   INNER JOIN coinc_event_map AS cm1
                   ON (si.simulation_id = cm1.event_id)
                   INNER JOIN coinc_event_map AS cm2
                   ON (cm1.coinc_event_id = cm2.coinc_event_id)
                   WHERE cm2.event_id = ?'''
        (ra, dec), = opts.inj_database.execute(
            query, (metadata['objid'],)).fetchall()
        radecs.append(np.rad2deg([ra, dec]).tolist())
    
    if opts.annotate:
        text = []
        try:
            objid = metadata['objid']
        except KeyError:
            pass
        else:
            text.append('event ID: {}'.format(objid))
        if opts.contour:
            pp = np.round(opts.contour).astype(int)
            ii = np.round(np.searchsorted(np.sort(cls), opts.contour) *
                          deg2perpix).astype(int)
#    print(pp, ii)
    return [pp, ii]
