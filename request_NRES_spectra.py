#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
#
#    Submit NRES targets for the coming night.
#
# Rob Siverd
# Created:       2016-04-04
# Last modified: 2019-01-30
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "1.1.5"

## Modules:
import re
import ast
import signal
import getopt
import os
import sys
import time
import numpy as np
import httplib
import urllib
import json
#import copy
#import getpass
import requests

import datetime as dt
import ephem
import ast

#import ascii as asc
#import rjs_request as rjsreq

## LCO info module:
try:
    import lco_site_info
    reload(lco_site_info)
except ImportError:
    sys.stderr.write("\nRequired module 'lco_site_info' not found!\n\n")
    sys.exit(1)
lsi = lco_site_info.info

## Time conversion:
try:
    import astropy.time as astt
except ImportError:
    sys.stderr.write("\nError: astropy module not installed!\n"
           "Please install and try again.\n\n")
    sys.exit(1)

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
## Catch interruption cleanly:
def signal_handler(signal, frame):
    sys.stderr.write("\nInterrupted!\n\n")
    sys.exit(1)

signal.signal(signal.SIGINT, signal_handler)

##--------------------------------------------------------------------------##
## Date/time formatter:
def tformat(when):
    #return when.isoformat().replace('T', ' ')
    iso_date = when.isoformat().replace('T', ' ')
    return iso_date.split('.')[0]

## Remove decimal and beyond:
def strip_dots(stime):
    return stime.split('.')[0]

## Microsecond removal:
def strip_msec(when):
    return (when - dt.timedelta(microseconds=when.microsecond))

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Mapping of options to molecule types:
molmap = {       '--ENG':    'ENGINEERING',
            '--SPECTRUM':  'NRES_SPECTRUM',
                '--TEST':      'NRES_TEST',
                }

## Config:
timer           = False
vlevel          = 0
#obs_note        = 'blablabla'
cred_file       = "credentials.txt"
contirmed       = True
yesterday       = False
obs_period      = None
night_hours     = 10.0
#objs_per_night  = 12
#objs_per_night  =  2
objs_to_submit  =  0
submit_cadence  = False
#lsc_lat         =  -30.1673666667
#lsc_lon         =  -70.8049
#min_Vmag        = -5.0          # ignore objects brighter than this
#max_Vmag        = 15.0          # ignore objects fainter than this
Vmag_limits     = [-5.0, 15.0]  # restrict Vmag to sensible range by default
window_nights   =  1            # number of observing nights in request window
acquire_bright  = True          # by default, acquire on brightest
key_project     = False         # if true, use key project propid
nres_type       = None
site_choice     = None
commissioning   = None
target_filter   = None          # to select target(s) by name
python_strategy = False         # enables 'python' strategies (testing)
nres_site_list  = [x for x in lsi.keys() if lsi[x]['nres_spec'] != None]
target_list     = None
constraints = {
            'max_airmass'   :   4.0, 
        'min_lunar_distance':   0.0,
        }
#cadence_horizon_days = 2.0  # schedule this far in advance

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Valhalla, URLs, and related:
req_urlbase = "http://lco.global/observe/request/"

## Load credentials from external file (credentials.txt):
if not os.path.isfile(cred_file):
    sys.stderr.write("Error: credential file not found: %s\n" % cred_file)
    sys.exit(1)
try:
    with open(cred_file, 'r') as f:
       credentials = ast.literal_eval(f.read())
    auth_token = credentials['auth_token']
except:
    sys.stderr.write("Failed to load credentials!!\n")
    sys.exit(1)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Guider mode:
guider_args = {
    #'ag_name'               : 'kb35',  # spectrograph AGU camera
    'ag_mode'               : 'ON',
    #'ag_filter'             : 'LL',    # enforced in site-software
    #'ag_filter'             : 'BL',    # installed but not selectable
    #'ag_filter'             : 'RL',    # installed but not selectable
}

## Global molecule config:
molecule_defaults = {
    #'instrument_name'       : '1m0-NRES-SciCam',
    'instrument_name'       : '1M0-NRES-SCICAM',
    #'instrument_name'       : '1M0-NRES-COMMISSIONING',
    #'exposure_time'         : exp_time,             # exposure time (seconds)
    #'exposure_count'        : num_exps,             # how many exposure to take
    #'spectra_slit'          : 'slit_1.6as',  # The generic filter name
    'fill_window'           : False,
    #'type'                  : 'NRES_SPECTRUM',  # The type of the molecule
    #'type'                  : 'NRES_TEST',  # The type of the molecule
    #'type'                  : 'ENGINEERING',
    'bin_x'                 : 1,
    'bin_y'                 : 1,
    #'ag_strategy'           : 'pattern',
    'ag_strategy'           : 'window_guiding',

    #'acquire_strategy'      : 'catalogue',
    'acquire_mode'          : 'BRIGHTEST',
    'acquire_strategy'      : 'astrometry',
    'acquire_radius_arcsec' : 5,            # for acquire-on-brightest
}
molecule_defaults.update(guider_args)


## Default cadence is daily:
window  = {}
#cadence =  {
#        "period": 24.0,                   # by default, observe daily
#        "jitter": 10.0,                   # --> doesn't matter what time
#}

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Current semester:
#semester_start = astt.Time("2018-06-01 00:00:00")
#semester_end   = astt.Time("2018-12-01 00:00:00")
#semester_start = astt.Time("2018-12-01 00:00:00")
#semester_end   = astt.Time("2019-06-01 00:00:00")

## Determine current semester:
sem_url = 'https://observe.lco.global/api/semesters/'
now_utc = dt.datetime.utcnow()
params = {'semester_contains':now_utc.isoformat()}
try:
    rr = requests.get(sem_url, params=params)
    rr.raise_for_status()
    this_semester = rr.json()['results'][0]
except:
    sys.stderr.write("Failed to retrieve current semester. API offline?\n")
    sys.exit(1)

sys.stderr.write("Current semester is: %s\n" % this_semester['id'])
semester_start = astt.Time(this_semester['start'])
semester_end   = astt.Time(this_semester['end'])
sys.stderr.write("Start: %s\n" % str(semester_start))
sys.stderr.write("End:   %s\n" % str(semester_end))

##--------------------------------------------------------------------------##
## Argument type-checking:
def is_integer(asdf):
    try:
        int(asdf)
        return True
    except ValueError:
        return False

def is_float(asdf):
    try:
        float(asdf)
        return True
    except ValueError:
        return False

##--------------------------------------------------------------------------##
##*********************     Help and options menu:     *********************##
##--------------------------------------------------------------------------##

## Syntax / how to run:
def usage(stream):
    stream.write("\n"
        + "Usage: %s [options] \n" % __file__
        + "Submit a bunch of NRES requests (generalized version).\n"
        + "Version: %s\n" % __version__
        + "\n"
        #+ "Confirm submission:\n"
        #+ "       --CONFIRM        actually submit observations to POND\n"
        #+ "\n"
        + "Target list:\n"
        + "   -f, --filter=REGEX   keep only objects matching REGEX\n"
        + "   -n, --nobjs=N        set maximum number of submissions\n"
        + "   -t, --targets=FILE   path to ASCII target list [REQUIRED]\n"
        + "       --Vmin=MAG       only consider targets with V >= MAG\n"
        + "       --Vmax=MAG       only consider targets with V <= MAG\n"
        + "\n"
        + "Instrument/proposal type (REQ):\n"
        + "   -C, --commissioning  use commissioning proposal and instrument\n"
        + "   -S, --science        use SCICAM instrument\n"
        + "   -K, --keyproj        use NRES key project proposal [optional]\n"
        + "\n"
        + "Site selection (REQ):\n"
        + "       --any            no site restriction for submitted requests\n"
        + "   -s, --site=SITE      restrict request to specific LCO site.\n"
        + "                           options: %s\n" % str(nres_site_list)
        + "\n"
        + "Observation type (REQ):\n"
        + "       --ENG            submit ENGINEERING molecules (mwa ha ha)\n"
        + "       --SPECTRUM       submit NRES_SPECTRUM molecules\n"
        + "       --TEST           submit NRES_TEST molecules\n"
        + "\n"
        + "Testing options:\n"
        + "       --pystrat        invokes 'python' strategies (testing)\n"
        + "\n"
        + "Scheduling time window:\n"
        + "   -N, --nights=N       number of nights in req window "
        +                                 "[def: %d]\n" % window_nights
        + "   -T, --today          fill around *next* midnight [default]\n"
        + "   -Y, --yesterday      schedule around *previous* midnight\n"
        + "\n"
        + "Available options:\n"
        + "       --debug          extra debugging info\n"
        + "   -h, --help           print this page\n"
        + "   -q, --quiet          suppress unnecessary output\n"
        + "   -t, --timer          report program run-time\n"
        + "   -v, --verbose        more status updates\n"
        + "\n")

##--------------------------------------------------------------------------##
##*********************       Parse command line:      *********************##
##--------------------------------------------------------------------------##

## Options:
short_opts = 'f:n:t:CSKs:N:TYhqtv'
long_opts = ['filter=', 'nobjs=', 'targets=', 'Vmin=', 'Vmax=',
            'commissioning', 'science', 'keyproj',
            'any', 'site=', 'ENG', 'SPECTRUM', 'TEST', 'pystrat',
            'nights=', 'today', 'yesterday',
            'debug', 'help', 'quiet', 'timer', 'verbose']

## GNU-style parsing (with exception handling):
try:
    options, remainder = getopt.gnu_getopt(sys.argv[1:], short_opts, long_opts)
except getopt.GetoptError, err:
    sys.stderr.write("%s\n" % str(err))
    usage(sys.stderr)
    sys.exit(2)

## Handle selected options:
for opt, arg in options:
    # ------------------------------------------------
    if (opt == '--debug'):
        debug = True
        sys.stderr.write(BRED + "\nDebugging output enabled!" + ENDC + "\n")
    # ------------------------------------------------
    elif (opt == '--CONFIRM'):
        confirmed = True
        msg = BYELLOW + "\nSubmission confirmed!" + ENDC + "\n"
        sys.stderr.write(msg)
    # ------------------------------------------------
    elif ((opt == '-f') or (opt == '--filter')):
        target_filter = re.compile(arg)
        msg = "Using target name filter: %s" % target_filter.pattern
        sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    elif ((opt == '-n') or (opt == '--nobjs')):
        if not is_integer(arg):
            sys.stderr.write("Error!  Non-integer nobjs: '%s'\n\n" % arg)
            sys.exit(1)
        objs_to_submit = int(arg)
        if (vlevel >= 0):
            msg = "Submission limit: %d targets" % objs_to_submit
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    elif ((opt == '-t') or (opt == '--targets')):
        target_list = arg
        if not os.path.isfile(target_list):
            sys.stderr.write(NRED + "File not found: '%s'\n\n" % target_list)
            sys.exit(1)
        msg = NYELLOW + "Targets from: %s" % target_list + ENDC + "\n"
        sys.stderr.write(msg)
    elif (opt in ['--Vmin', '--Vmax']):
        _Vmap = {'--Vmin':0, '--Vmax':1}
        if not is_float(arg):
            sys.stderr.write("Error!  Non-numeric Vmag: '%s'\n\n" % arg)
            sys.exit(1)
        Vmag_limits[_Vmap[opt]] = float(arg)
    # ------------------------------------------------
    elif ((opt == '-C') or (opt == '--commission')):
        commissioning = True
        msg = NYELLOW + "Selected COMMISSIONING inst/proposal!" + ENDC + "\n"
        sys.stderr.write(msg)
    elif ((opt == '-S') or (opt == '--science')):
        commissioning = False
        msg = NYELLOW + "Selected SCIENCE inst/proposal!" + ENDC + "\n"
        sys.stderr.write(msg)
    elif ((opt == '-K') or (opt == '--keyproj')):
        key_project = True
        msg = NYELLOW + "Key project proposal selected!" + ENDC + "\n"
        sys.stderr.write(msg)
    # ------------------------------------------------
    elif ((opt == '--any')):
        site_choice = 'any'
        if (vlevel >= 0):
            msg = "No site restriction!"
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    elif ((opt == '-s') or (opt == '--site')):
        site_choice = arg.lower()
        if not site_choice in nres_site_list:
            sys.stderr.write("Invalid site: '%s'\n\n" % arg)
            usage(sys.stderr)
            sys.exit(1)
        if (vlevel >= 0):
            msg = "Requests restricted to: " + arg
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    elif (opt in ['--ENG', '--SPECTRUM', '--TEST']):
        nres_type = molmap.get(opt)
        molecule_defaults['type'] = nres_type
        msg = NYELLOW + "Molecule type: %s" % nres_type + ENDC + "\n"
        sys.stderr.write(msg)
    elif ((opt == '--pystrat')):
        python_strategy = True
        msg = NYELLOW + "WARNING: Python strategies invoked!" + ENDC + "\n"
        sys.stderr.write(msg)

    # ------------------------------------------------
    elif ((opt == '-N') or (opt == '--nights')):
        if not is_integer(arg):
            sys.stderr.write("Error! Non-integer nights: '%s'\n\n" % arg)
            sys.exit(1)
        window_nights = int(arg)
        if (window_nights < 1):
            sys.stderr.write("Number of nights >= 1 required!\n")
            sys.exit(1)
        if (vlevel >= 0):
            msg = "Nights in request window: %d" % window_nights
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    elif ((opt == '-T') or (opt == '--today')):
        yesterday = False
        msg = NYELLOW + "Scheduling next midnight!" + ENDC + "\n"
        sys.stderr.write(msg)
    elif ((opt == '-Y') or (opt == '--yesterday')):
        yesterday = True
        msg = NYELLOW + "Scheduling prev midnight!" + ENDC + "\n"
        sys.stderr.write(msg)
    # ------------------------------------------------
    elif ((opt == '-h') or (opt == '--help')):
        usage(sys.stdout)
        sys.exit(0)
    elif ((opt == '-q') or (opt == '--quiet')):
        vlevel -= 1
    elif ((opt == '-t') or (opt == '--timer')):
        timer = True
    elif ((opt == '-v') | (opt == '--verbose')):
        vlevel += 1
        sys.stderr.write(NYELLOW + "Increasing verbosity." + ENDC + "\n")
    # ------------------------------------------------
    else:
        msg = "Unhandled option: %s" % opt
        sys.stderr.write(BRED + "\n" + msg + ENDC + "\n\n")
        sys.exit(1)
    pass

## Verbosity:
if (vlevel >= 1):
    sys.stderr.write("%sVerbosity level: %d%s\n" % (NYELLOW, vlevel, ENDC))

## Full command line if highly verbose:
if (vlevel >= 2):
    sys.stderr.write("%s\nFull command line:%s\n" % (NCYAN, ENDC))
    sys.stderr.write("   %s\n" % sys.argv)

##--------------------------------------------------------------------------##
##                          Input sanity checks:                            ##
##--------------------------------------------------------------------------##

## Target list must be given:
if target_list == None:
    msg = "No target list provided!"
    sys.stderr.write(BRED + msg + ENDC + "\n")
    usage(sys.stderr)
    sys.exit(1)

## Commissioning/science mode must be selected:
if commissioning == None:
    msg = "Need to specify science/commissioning!"
    sys.stderr.write(BRED + msg + ENDC + "\n")
    usage(sys.stderr)
    sys.exit(1)

## Site choice is required:
if site_choice == None:
    msg = "Must specify SITE or use --any option ..."
    sys.stderr.write(BRED + msg + ENDC + "\n")
    usage(sys.stderr)
    sys.exit(1)

## Molecule type must be specified:
if (nres_type == None) or (molecule_defaults['type'] == None):
    msg = BRED + "Error: molecule type not specified!" + ENDC + "\n"
    sys.stderr.write(msg)
    usage(sys.stderr)
    sys.exit(1)

## Announce option choices:
#sys.stderr.write(NYELLOW + "Site choice:   %s" % site_choice + ENDC + "\n")
sys.stderr.write(NYELLOW + "Molecule type: %s" % nres_type + ENDC + "\n")

## Commissioning time overrides key project proposal choice:
if commissioning and key_project:
    sys.stderr.write("COMMISSIONING mode not compatible with key project!\n")
    sys.stderr.write("Proceeding with engineering proposal.\n")
    key_project = False

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Midnight selector (for single-site obs):
def pick_midnight(site_choice, yest):
    # Define site:
    site = ephem.Observer()
    site.lat = np.radians(lsi[site_choice]['latitude'])
    site.lon = np.radians(lsi[site_choice]['longitude'])

    # Determine previous midnight (so window includes tonight):
    sun = ephem.Sun()
    prev_midnight = strip_msec(site.previous_antitransit(sun).datetime())
    next_midnight = strip_msec(site.next_antitransit(sun).datetime())
    return prev_midnight if yest else next_midnight

## Set up request location and time window:
#location = {'telescope_class':'1m0', 'site':'lsc'} #, 'observatory':'domb'}
if site_choice in nres_site_list:
    location = {'telescope_class':'1m0', 'site':site_choice}
    use_midnight = pick_midnight(site_choice, yesterday)
elif (site_choice == 'any'):
    location = {'telescope_class':'1m0'}
    use_midnight = dt.datetime.utcnow() + dt.timedelta(hours=night_hours)
else:
    sys.stderr.write(BRED + "UNHANDLED site_choice??!!" + ENDC + "\n")
    sys.exit(1)


## Select commissioning vs. science instrument:
if commissioning:
    #proposal_id    = 'ENG2017AB-001'     # NRES commissioning
    molecule_defaults['instrument_name'] = '1M0-NRES-COMMISSIONING'
else:
    #proposal_id    = 'KEY2017AB-002a'    # NRES key project
    molecule_defaults['instrument_name'] = '1M0-NRES-SCICAM'

## Choose proposal (engineering or key project):
if key_project:
    proposal_id    = 'KEY2017AB-002a'    # NRES key project
else:
    proposal_id    = 'ENG2017AB-001'     # NRES commissioning

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
#target_list = "./target_list_nres.txt"
#hdr_txt = asc.get_headtext(target_list)
all_object_data = np.genfromtxt(target_list, dtype=None, delimiter=',',
        names=True, autostrip=True)

## Impose magnitude limit:
n_load = len(all_object_data)
sys.stderr.write("Using mag limits: %.2f <= Vmag <= %.2f\n"
        % tuple(Vmag_limits))
too_bright = (all_object_data['Vmag'] < Vmag_limits[0])
too_faint  = (all_object_data['Vmag'] > Vmag_limits[1])
use_object_data = all_object_data[~too_bright & ~too_faint]
n_used = len(use_object_data)
sys.stderr.write("Kept %d of %d objects.\n" % (n_used, n_load))


## Filter target list by name:
if target_filter != None:
    matching = np.bool_([target_filter.search(x) != None \
            for x in all_object_data['name']])
    if (np.sum(matching) == 0):
        sys.stderr.write("No targets matched: %s\n" % target_filter.pattern)
        sys.exit(1)
    use_object_data = use_object_data[matching] # keep subset

##--------------------------------------------------------------------------##
## Per-target molecule updates:
def configure_NRES_molecule(entry, mconfig):
    obj_name          = entry['name']
   #window_offset_hrs = entry['shift']
    acquire_mode      = entry['acq']

    # ------------------------------------------------
    # molecule updates:
    molecule = {}
    molecule.update(mconfig)
    molecule['exposure_time']  = entry['exptime']
    molecule['exposure_count'] = entry['nexp']
    molecule['ag_exp_time'] = 15.0  # test custom a/g exposure time
    mol_type = molecule['type'].upper()
    if mol_type != nres_type:
        sys.stderr.write("WEIRD!\n")
        sys.stderr.write("mol_type:  %s\n" % mol_type)
        sys.stderr.write("nres_type: %s\n" % nres_type)

    # Apply acquisition mode requested in target list:
    if acquire_mode == 'wcs':
        molecule['acquire_mode'] = "WCS"
        molecule['acquire_radius_arcsec'] = 0
    elif acquire_mode == 'bri':
        molecule['acquire_mode'] = "BRIGHTEST"
        molecule['acquire_radius_arcsec'] = 5
    else:
        sys.stderr.write("Unhandled acquisition mode: %s\n" % acquire_mode)
        sys.exit(1)

    # In python-testing mode, use custom strategies:
    if python_strategy:
        sys.stderr.write("-----------------------------\n")
        sys.stderr.write("PYTHON STRATEGIES INVOKED!!!!\n")
        sys.stderr.write("-----------------------------\n")
        molecule['acquire_radius_arcsec'] = 5
        molecule['acquire_strategy'] = 'python'
        molecule[     'ag_strategy'] = 'python'
        #del molecule['ag_strategy']
        #del molecule['acquire_strategy']
        #molecule['acquire_strategy'] = ""
        #molecule[     'ag_strategy'] = ""

    ## Override target list parameters for specific molecule types:
    #if mol_type in ['ENGINEERING', 'NRES_EXPOSE']:
    #    molecule['acquire_mode'] = 'OFF'
    #    molecule['acquire_radius_arcsec'] = 0
    #    molecule['ag_mode']     = 'ON'
    #    molecule['ag_strategy'] = 'pattern'
    #    sys.stderr.write("mol_type: %s\n" % mol_type)

    # ------------------------------------------------
    # Append molecule type tag to object name:
    if molecule['type'].upper() == 'NRES_TEST':
        obj_name += '_nrt'
    if molecule['type'].upper() == 'ENGINEERING':
        obj_name += '_engr'
    if python_strategy:
        obj_name += '_pystrat'

    # ------------------------------------------------
    # Override exposure time and number for NRES_TEST:
    if molecule['type'].upper() == 'NRES_TEST':
        molecule['exposure_time'] = 300.0
        molecule['exposure_count'] = 1
        sys.stderr.write("Enforcing 1x 300-sec exposure for NRES_TEST!\n")

    # ------------------------------------------------
    # Announce final molecule parameters:
    sys.stderr.write("%s\n" % str(molecule))

    # ------------------------------------------------
    # name enhancements:
    if molecule.get('ag_filter'):
        obj_name += '_%s' % molecule['ag_filter']

    # ------------------------------------------------
    # define the target
    target = {
          'name'              : obj_name,
          'ra'                : entry['RA'],
          'dec'               : entry['DE'],
          'proper_motion_ra'  : entry['pmRA'],
          'proper_motion_dec' : entry['pmDE'],
          'parallax'          : entry['prlx'],
          'vmag'              : entry['Vmag'],
          'radvel'            : entry['radvel'],
          'epoch'             : 2000,
          'type'              : 'SIDEREAL',
    }

    return (molecule, target)

##--------------------------------------------------------------------------##
## Outout config:
#stream = sys.stdout
stream = sys.stderr

## Submit stuff:
errors = 0
naccepted = 0
nrejected = 0
for entry in use_object_data:
    stream.write("\n------------------------------------------------\n")
    stream.write("entry: %s\n" % str(entry))
    obj_name          = entry['name']
   #window_offset_hrs = entry['shift']
    acquire_mode      = entry['acq']
    active            = bool(entry['active'])
    obs_note = entry['note'] if 'note' in entry.dtype.names else ''
    if not active:
        continue
 
    nres_molecule, target = configure_NRES_molecule(entry, molecule_defaults)

    # ------------------------------------------------
    # Start/stop of window:
    obj_window_start = use_midnight - dt.timedelta(hours=night_hours)
    obj_window_end   = use_midnight + dt.timedelta(hours=night_hours)

    # Increase window by 24h for each additional night requested:
    extra_days = window_nights - 1
    obj_window_end += dt.timedelta(days=extra_days)
    #sys.stderr.write("obj_window_start: %s\n" % str(obj_window_start))
    #sys.stderr.write("obj_window_end:   %s\n" % str(obj_window_end))

    # Adjust for semester boundaries and current time:
    if obj_window_start < semester_start.datetime:
        obj_window_start = semester_start.datetime
    if obj_window_end > semester_end.datetime:
        obj_window_end = semester_end.datetime - dt.timedelta(seconds=2)
    if obj_window_end < dt.datetime.utcnow():
        stream.write("Window has passed!\n")
        continue
    if obj_window_start < dt.datetime.utcnow():
        obj_window_start = dt.datetime.utcnow() + dt.timedelta(seconds=120)
        stream.write("Fudged start time!\n")

    window["start"] = strip_dots(tformat(obj_window_start))
    window["end"  ] = strip_dots(tformat(obj_window_end))
    if vlevel >= 1:
        stream.write("Final object window:   %s  <-->  %s\n"
                % (window['start'], window['end']))

    # ------------------------------------------------
    # build request:
    request = {
        "constraints" : constraints,
        "location" : location,
        "molecules" : [nres_molecule],
        "observation_note" : obs_note,
        "observation_type" : "NORMAL",
        "target" : target,
        "type" : "request",
        "windows" : [window],
    }

    ipp_value = 1.05
    new_user_request = {
        "group_id"          : "NRES: %s" % target['name'],
        "ipp_value"         :                   ipp_value,
        "operator"          :                    "SINGLE",
        "observation_type"  :                    "NORMAL",
        "requests"          :                   [request],
        "proposal"          :                 proposal_id,
        #"type"              : "compound_request",
    }

    sys.exit(1)

    # ------------------------------------------------
    # Submit:
    stream.write("\n")
    stream.write("Current time: %s\n" % dt.datetime.utcnow().isoformat())
    stream.write("Submitting %s ... " % target['name'])
    rr = requests.post('https://observe.lco.global/api/userrequests/',
            headers={'Authorization': 'Token {}'.format(auth_token)},
            json=new_user_request)

    # Inspect response:
    try:
        rr.raise_for_status()
        stream.write("success!\n")
        tracknum = rr.json()['id']
        stream.write("Tracking number: %s\n" % str(tracknum))
        naccepted += 1
    except requests.exceptions.HTTPError as exc:
        stream.write("FAILED!!\n\n")
        stream.write("Server response: {}\n".format(rr.content))
        nrejected += 1
    stream.write("\n")

    # Early stop based on accepted requests:
    if (objs_to_submit > 0) and (naccepted >= objs_to_submit):
        break

##--------------------------------------------------------------------------##
## Brief summary:
stream.write("\n\n")
stream.write("Requests accepted: %d\n" % naccepted)
stream.write("Requests rejected: %d\n" % nrejected)
stream.write("\n\n")

#sys.exit(0)
sys.exit(errors)

######################################################################
# CHANGELOG (request_NRES_spectra.py):
#---------------------------------------------------------------------
#
#  2019-01-30:
#     -- Increased __version__ to 1.1.5.
#     -- Now programmatically retrieve current semester info from LCO API.
#     -- Updated to current semester.
#     -- Added test for observation note.
#
#  2018-11-15:
#     -- Increased __version__ to 1.1.0.
#     -- Python strategy is now set in the request again.
#
#  2018-07-16:
#     -- Increased __version__ to 1.0.9.
#     -- Disabled setting of python strategies temporarily for debugging.
#
#  2018-05-22:
#     -- Increased __version__ to 1.0.8.
#     -- Now append _pystrat to object name when those strategies are invoked.
#
#  2018-05-01:
#     -- Increased __version__ to 1.0.7.
#     -- Merged min_Vmag and max_Vmag into Vmag_limits. Added --Vmin and
#           --Vmax command-line options to specify these values.
#     -- Added max_Vmag with default value of 15.
#     -- Increased __version__ to 1.0.6.
#     -- Added --pystrat option to invoke 'python' strategies (testing).
#
#  2018-02-12:
#     -- Increased __version__ to 1.0.5.
#     -- Added -K, --keyproj option to separately control which proposal is
#           used for a request. This permits me to use 1M0-NRES-SCICAM time on
#           the ENG2017AB-001 proposal.
#
#  2017-12-22:
#     -- Increased __version__ to 1.0.4.
#     -- Added -n, --nights=N option to specify request window length.
#
#  2017-12-21:
#     -- Increased __version__ to 1.0.3.
#     -- Added -f, --filter=REGEX option for target name matching.
#     -- Now import re module for regex-based target name filtering.
#     -- Now include target V mag and radial velocity in request.
#     -- Changed min_Vmag from 0 to -5 (no immediate effect expected).
#
#  2017-12-20:
#     -- Increased __version__ to 1.0.2.
#     -- Eliminated ascii and rjsreq module imports (no longer needed).
#     -- Now use requests raise_for_status() to inspect Valhalla response.
#     -- Increased __version__ to 1.0.1.
#     -- Acquisition is no longer disabled. Sequencer changes obviate
#           the need. Further, overheads in the scheduler seem to be grossly
#           underestimated with acquisition disabled.
#
#  2017-12-18:
#     -- Increased __version__ to 1.0.0.
#     -- Added -n, --nobjs=N argument to limit number of requests made.
#     -- Added -t, --targets=FILE (required) to specify target list.
#     -- Added --any, -s, --site=SITE option for request site selection.
#     -- Choice of science/commissioning on CL is now required.
#     -- Changed script name to request_NRES_spectra.py (general version).
#
#  2017-12-13:
#     -- Increased __version__ to 0.7.6.
#     -- Now set ag_mode='OFF' for ENGINEERING molecule.
#
#  2017-12-11:
#     -- Increased __version__ to 0.7.5.
#     -- Added -S, --science and -C, --commission command line options.
#           These designate matching proposal IDs and instrument types.
#
#  2017-11-07:
#     -- Increased __version__ to 0.7.0.
#     -- Molecule type is now a required CL input.
#     -- Fixed not-plumbed-through molecule type selection (yikes!).
#     -- Added --ENG option to submit ENGINEERING molecules.
#
#  2017-10-18:
#     -- Increased __version__ to 0.5.9.
#     -- Switched to 1M0-NRES-COMMISSIONING instrument class.
#
#  2017-07-26:
#     -- Increased __version__ to 0.5.8.
#     -- Added min_Vmag parameter to help control object selection.
#     -- Moved back to domb.
#
#  2017-07-19:
#     -- Increased __version__ to 0.5.7.
#     -- Pinned my submissions to domc.
#
#  2017-07-05:
#     -- Increased __version__ to 0.5.6.
#     -- For testing, now request ag_exp_time of 15 seconds.
#
#  2017-06-06:
#     -- Increased __version__ to 0.5.5.
#     -- Added observation_type and ipp_value to user_request and type to
#           target dict to resolve Valhalla deprecation warnings.
#
#  2017-06-03:
#     -- Increased __version__ to 0.5.0.
#     -- Added acquire_mode = brightest to default molecule args.
#     -- Added command-line option parsing with getopt.
#     -- Now use objs_per_night config item (3 submissions currently).
#
#  2017-04-30:
#     -- Increased __version__ to 0.4.1.
#     -- Now explicitly set acquire_radius_arcsec to 0 in WCS mode.
#
#  2017-04-28:
#     -- Increased __version__ to 0.4.0.
#     -- Added objs_per_night config item, start with 3.
#
#  2017-04-10:
#     -- Increased __version__ to 0.3.6.
#     -- Added type-specific object name suffixes for NRES_TEST and
#           ENGINEERING molecules for easier identification in the POND.
#
#  2017-04-04:
#     -- Increased __version__ to 0.3.5.
#     -- Now try different objects until one gets accepted.
#
#  2016-11-13:
#     -- Increased __version__ to 0.3.0.
#     -- acquire_mode is no longer appended to object name prior to submission
#           (now already appended in the target list).
#     -- Target name field is now 'name' (previously 'object').
#
#  2016-10-18:
#     -- Increased __version__ to 0.2.7.
#     -- Changed LCOGT --> LCO.
#
#  2016-10-04:
#     -- Increased __version__ to 0.2.6.
#     -- Updated semester bounds for 2016B.
#
#  2016-04-20:
#     -- Increased __version__ to 0.2.5.
#     -- Now use 'exptime' and 'nexp' keys to access object data instead of
#           array indices. These were overlooked previously and caused trouble
#           when the CSV format changed slightly.
#
#  2016-04-09:
#     -- Increased __version__ to 0.2.0.
#     -- Failed requests are tracked and reported in script exit status.
#     -- Now report UTC time just before each submission for reference.
#     -- Embedded ODIN password for automatic submission.
#
#  2016-04-05:
#     -- Increased __version__ to 0.1.5.
#     -- Now use rjs_request molecule to clean up code.
#
#  2016-04-04:
#     -- Increased __version__ to 0.1.0.
#     -- First created weekly_nres_cadences.py.
#
