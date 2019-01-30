#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Look up a target by name and add to object list.
#
# Rob Siverd
# Created:       2018-05-01
# Last modified: 2018-07-14
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.1.1"

## Python version-agnostic module reloading:
try:
    reload                              # Python 2.7
except NameError:
    try:
        from importlib import reload    # Python 3.4+
    except ImportError:
        from imp import reload          # Python 3.0 - 3.3

## Modules:
import argparse
import os
import sys
import time
import numpy as np
#from numpy.lib.recfunctions import append_fields
#import datetime as dt
#from dateutil import parser as dtp

## Look-up assistant:
try:
    import simbad_helper
    reload(simbad_helper)
except ImportError:
    sys.stderr.write("Module 'simbad_helper' not found!\n")
    sys.exit(1)
nrl = simbad_helper.NRES_Lookup()

##--------------------------------------------------------------------------##

## ASCII I/O:
#try:
#    import astropy.io.ascii as aia
#except ImportError:
#    sys.stderr.write("\nError: astropy module not found!\n")
#    sys.exit(1)

##--------------------------------------------------------------------------##
## Colors for fancy terminal output:
NRED    = '\033[0;31m'   ;  BRED    = '\033[1;31m'
NGREEN  = '\033[0;32m'   ;  BGREEN  = '\033[1;32m'
NYELLOW = '\033[0;33m'   ;  BYELLOW = '\033[1;33m'
NBLUE   = '\033[0;34m'   ;  BBLUE   = '\033[1;34m'
NMAG    = '\033[0;35m'   ;  BMAG    = '\033[1;35m'
NCYAN   = '\033[0;36m'   ;  BCYAN   = '\033[1;36m'
NWHITE  = '\033[0;37m'   ;  BWHITE  = '\033[1;37m'
ENDC    = '\033[0m'

## Suppress colors in cron jobs:
if (os.getenv('FUNCDEF') == '--nocolors'):
    NRED    = ''   ;  BRED    = ''
    NGREEN  = ''   ;  BGREEN  = ''
    NYELLOW = ''   ;  BYELLOW = ''
    NBLUE   = ''   ;  BBLUE   = ''
    NMAG    = ''   ;  BMAG    = ''
    NCYAN   = ''   ;  BCYAN   = ''
    NWHITE  = ''   ;  BWHITE  = ''
    ENDC    = ''

## Fancy text:
degree_sign = u'\N{DEGREE SIGN}'

## Dividers:
halfdiv = "----------------------------------------"
fulldiv = halfdiv + halfdiv

##--------------------------------------------------------------------------##
def ldmap(things):
    return dict(zip(things, range(len(things))))

def argnear(vec, val):
    return (np.abs(vec - val)).argmin()

## Robust location/scale estimate using median/MAD:
def calc_ls_med_MAD(a, axis=None):
    """Return median and median absolute deviation of *a* (scaled to normal)."""
    med_val = np.median(a, axis=axis)
    sig_hat = (1.482602218 * np.median(np.abs(a - med_val), axis=axis))
    return (med_val, sig_hat)

## Robust location/scale estimate using median/IQR:
def calc_ls_med_IQR(a, axis=None):
    """Return median and inter-quartile range of *a* (scaled to normal)."""
    pctiles = np.percentile(a, [25, 50, 75], axis=axis)
    med_val = pctiles[1]
    sig_hat = (0.741301109 * (pctiles[2] - pctiles[0]))
    return (med_val, sig_hat)

## Select inliners given specified sigma threshold:
def pick_inliers(data, sig_thresh):
    med, sig = calc_ls_med_IQR(data)
    return ((np.abs(data - med) / sig) <= sig_thresh)

##--------------------------------------------------------------------------##
## Parse arguments and run script:
if __name__ == '__main__':

    # ------------------------------------------------------------------
    descr_txt = """
    Look up object by name using SIMBAD. Append object to target list.
    
    Version: %s
    """ % __version__
    parser = argparse.ArgumentParser(
            prog=os.path.basename(__file__),
            description=descr_txt)
    parser.add_argument('lookup_name', help='name for SIMBAD lookup')
    parser.add_argument('target_name', help='target name in database')
    parser.add_argument('-o', '--output_file', default=None,
            help='Target list object should append to')
    parser.add_argument('--debug', dest='debug', default=False,
             help='Enable extra debugging messages', action='store_true')
    parser.add_argument('-q', '--quiet', action='count', default=0)
    parser.add_argument('-v', '--verbose', action='count', default=0)
    #parser.add_argument('remainder', help='other stuff', nargs='*')
    # ------------------------------------------------------------------
    #ofgroup = parser.add_argument_group('Output format')
    #fmtparse = ofgroup.add_mutually_exclusive_group()
    #fmtparse.add_argument('--python', required=False, dest='output_mode',
    #        help='Return Python dictionary with results [default]',
    #        default='pydict', action='store_const', const='pydict')
    #bash_var = 'ARRAY_NAME'
    #bash_msg = 'output Bash code snippet (use with eval) to declare '
    #bash_msg += 'an associative array %s containing results' % bash_var
    #fmtparse.add_argument('--bash', required=False, default=None,
    #        help=bash_msg, dest='bash_array', metavar=bash_var)
    # ------------------------------------------------------------------

    context = parser.parse_args()
    #context.vlevel = context.verbose - context.quiet
    context.vlevel = 99 if context.debug else (context.verbose-context.quiet)

    result = nrl.fetch_simbad(context.lookup_name, context.target_name)
    fmtobj = nrl.targ2csv(result)
    sys.stderr.write("result:\n")
    sys.stdout.write("%s\n" % fmtobj)

    # If output given, append to target list:


######################################################################
# CHANGELOG (add_target.py):
#---------------------------------------------------------------------
#
#  2018-07-14:
#     -- Increased __version__ to 0.1.1.
#     -- Target list line is now sent to stdout for use/redirection.
#
#  2018-05-01:
#     -- Increased __version__ to 0.1.0.
#     -- First created add_target.py.
#
