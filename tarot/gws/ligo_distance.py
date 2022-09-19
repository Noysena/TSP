#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 16:34:55 2019

@author: noysena

Skymap area and distance
"""

"""
Obtain distance of event with erorr form fits.gz
"""


from . import ArgumentParser, FileType, figure_parser


def parser():
    parser = ArgumentParser(parents=[figure_parser])
    parser.add_argument(
        '--annotate', default=False, action='store_true',
        help='annotate plot with information about the event')
    parser.add_argument(
        '--max-distance', metavar='Mpc', type=float,
        help='maximum distance of plot in Mpc')
    parser.add_argument(
        '--contour', metavar='PERCENT', type=float, nargs='+',
        help='plot contour enclosing this percentage of'
        ' probability mass')
    parser.add_argument(
        '--radecdist', nargs=3, type=float, action='append', default=[],
        help='right ascension (deg), declination (deg), and distance to mark')
    parser.add_argument(
        '--chain', metavar='CHAIN.hdf5', type=FileType('rb'),
        help='optionally plot a posterior sample chain')
    parser.add_argument(
        '--projection', type=int, choices=list(range(4)), default=0,
        help='Plot one specific projection, or 0 for all projections')
    parser.add_argument(
        'input', metavar='INPUT.fits[.gz]', type=FileType('rb'),
        default='-', nargs='?', help='Input FITS file')
    parser.add_argument(
        '--align-to', metavar='SKYMAP.fits[.gz]', type=FileType('rb'),
        help='Align to the principal axes of this sky map')
    parser.set_defaults(figure_width='3.5', figure_height='3.5')
    return parser


def main(args=None):
    opts = parser().parse_args(args)

    # Late imports
    from .. import io
    import healpy as hp

    # Read input, determine input resolution.
#    progress.update(-1, 'Loading FITS file')
    (prob, mu, sigma, norm), metadata = io.read_sky_map(
        opts.input.name, distances=True)
    npix = len(prob)
    nside = hp.npix2nside(npix)

#    progress.update(-1, 'Preparing projection')
    distmean = metadata['distmean']
    diststd = metadata['diststd']
#    print("Mean Dist & Err", round(distmean), round(diststd))
    return [round(distmean), round(diststd)]