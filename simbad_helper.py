#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
#    Helper routines for SIMBAD object queries.
#
# Rob Siverd
# Created:       2018-01-02
# Last modified: 2019-01-30
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.1.7"

## Modules:
import os
import sys
import time
import numpy as np
from collections import OrderedDict

## Need astropy:
try:
    import astropy.coordinates as coord
except ImportError:
    sys.stderr.write("Module 'astropy' not found!  Install and retry ...\n")
    sys.exit(1)

## Need astroquery:
try:
    from astroquery.simbad import Simbad
except ImportError:
    sys.stderr.write("Module 'astroquery' not found!  Install and retry ...\n")
    sys.exit(1)

##--------------------------------------------------------------------------##
## NRES lookup class:
class NRES_Lookup(object):

    def __init__(self):
        self.mysearch = Simbad()
        self.mysearch.add_votable_fields('ra(d)', 'dec(d)')
        self.mysearch.add_votable_fields('propermotions', 'parallax')
        self.mysearch.add_votable_fields('rv_value', 'flux(V)', 'sptype')
        self._colspec = {
                'RA'      :    'RA_d',
                'DE'      :    'DEC_d',
              'pmRA'      :    'PMRA',
              'pmDE'      :    'PMDEC',
              'prlx'      :    'PLX_VALUE',
            'radvel'      :    'RV_VALUE',
              'Vmag'      :    'FLUX_V',
            'sptype'      :    'SP_TYPE',
        }
        self._defvals = [('period', 48.0), ('jitter', 12.0), ('active', 1),
                ('nexp', 1), ('shift', 0.0), ('acq', 'bri'), 
                ('note', 'unspecified')]

        self._txtspec =   [(   "name", 25, "%25s"),
                           (     "RA", 13, "%13.8f"),
                           (     "DE", 13, "%+13.8f"),
                           (   "pmRA", 10, "%10.3f"),
                           (   "pmDE", 10, "%10.3f"),
                           (   "prlx",  9, "%9.3f"),
                           ( "radvel",  8, "%8.3f"),
                           (   "Vmag",  8, "%8.3f"),
                           ("exptime",  8, "%8.1f"),
                           (   "nexp",  5, "%5d"),
                           ( "period",  7, "%7.1f"),
                           ( "jitter",  7, "%7.1f"),
                           (  "shift",  6, "%6.1f"),
                           (    "acq",  4, "%4s"),
                           ( "sptype", 16, "%16s"),
                           (   "note", 16, "%16s"),
                           ( "active",  7, "%7d")]
        return

    def _newitem(self, use_name):
        content = [('name', use_name)]
        content.extend(self._defvals)
        return OrderedDict(content)

    @staticmethod
    def _expcalc(Vmag, fallback=300.0):
        try:
            if (Vmag < 3):
                return  200.0
            elif (Vmag < 6):
                return  600.0
            else:
                return 1200.0
        except:
            return fallback

    # Extract desired quantities from a SIMBAD search result:
    def _qparse(self, result):
        content = []
        for kk,vv in self._colspec.items():
            content.append((kk, result[vv].data.data[0]))
        return dict(content)

    # Make a new line for target list entry:
    def fetch_simbad(self, lookup_name, use_name):
        target = self._newitem(use_name)

        # Include results from SIMBAD:
        sys.stderr.write("Querying SIMBAD ... ")
        self.results = self.mysearch.query_object(lookup_name)
        sys.stderr.write("done.\n")
        target.update(self._qparse(self.results))

        # Replacements for troublesome values:
        for datum in ['prlx', 'radvel']:
            if np.isnan(target[datum]):
                target[datum] = 0.0
        if len(target['sptype']) == 0:
            target['sptype'] = 'unknown'
        if hasattr(target['sptype'], 'decode'):
            target['sptype'] = target['sptype'].decode()

        # Add exposure time:
        target['exptime'] = self._expcalc(target['Vmag'])
        return target

    # ---------------------------------------------------------
    # print a CSV header line:
    def targ2header(self, entry, delim=','):
        col, nchars, fmt = zip(*self._txtspec)
        hfmt = ["%%%ds"%x for x in nchars]
        return delim.join([fmt%x for x,fmt in zip(col, hfmt)])

    # format a result for printing:
    def targ2csv(self, entry, delim=','):
        return delim.join([fmt % entry[x] for x,nc,fmt in self._txtspec])


######################################################################
# CHANGELOG (simbad_helper.py):
#---------------------------------------------------------------------
#
#  2018-07-14:
#     -- Increased __version__ to 0.1.5.
#     -- Added placeholders for troublesome values:
#           * 'nan' in prlx or radvel is replaced with 0.0
#           * empty string in 'sptype' field is replaced with 'unknown'
#     -- NOTE: EB target list format differs from master target list. FIX!
#     -- Added/removed columns and adjusted widths to match current target 
#           list contents and formatting. Columns pertaining to observing
#           frequency (period, jitter, shift) are now gone. New columns
#           pertaining to ephemerides (tc_jdutc, period, dur_t14) added.
#
#  2018-05-01:
#     -- Increased __version__ to 0.1.1.
#     -- Updated column widths to reflect current state of NRES target lists.
#
#  2018-01-02:
#     -- Increased __version__ to 0.1.0.
#     -- First created simbad_helper.py.
#
