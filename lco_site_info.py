#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
#
#    LCO observatory info for scripts.
#
# Rob Siverd
# Created:       2017-05-22
# Last modified: 2018-07-03
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.5.3"

## Modules:
import os
import sys
import time

##--------------------------------------------------------------------------##
## Site information:
_site_list = ['bpl', 'coj', 'cpt', 'elp', 'lsc', 'ogg', 'sqa', 'tfn', 'tlv']
info = dict([(x,{}) for x in _site_list])
#del x

##--------------------------------------------------------------------------##
## Coordinates:
_site_lonlat = {}
_site_lonlat['bpl'] = (-119.863103  ,  34.433161     )
_site_lonlat['coj'] = ( 149.0708466 , -31.2728196    )
_site_lonlat['cpt'] = (  20.8124    , -32.3826       )
_site_lonlat['elp'] = (-104.015173  ,  30.679833     )
_site_lonlat['lsc'] = ( -70.8049    , -30.1673666667 )
_site_lonlat['ogg'] = (-156.2589    ,  34.433161     )
_site_lonlat['sqa'] = (-120.04222167,  34.691453333  )
_site_lonlat['tfn'] = ( -16.511544  ,  28.300433     )
_site_lonlat['tlv'] = (  34.763333  ,  30.595833     )

for lsite,lonlat in _site_lonlat.items():
    if lsite in info.keys():
        info[lsite]['longitude'], info[lsite]['latitude'] = lonlat
del lsite, lonlat

##--------------------------------------------------------------------------##
## Per-site enclosure list:
_encl_list = {
        'bpl': ['aqwa', 'doma', 'igla'],
        'coj': ['clma', 'doma', 'domb'],
        'cpt': ['aqwa', 'doma', 'domb', 'domc', 'igla'],
        'elp': ['aqwa', 'doma', 'igla'],
        'lsc': ['aqwa', 'aqwb', 'doma', 'domb', 'domc', 'igla'],
        'ogg': ['clma'],
        'sqa': ['doma'],
        'tfn': ['aqwa'],
        'tlv': ['doma', 'igla'],
        }
for lsite,encls in _encl_list.items():
    if lsite in info.keys():
        info[lsite]['enclosures'] = encls
del lsite, encls

##--------------------------------------------------------------------------##
## Daily restart timing:
_restart_hr_utc = {
        'coj':  2,
        'tlv': 11,
        'cpt': 11,
        'tfn': 12,
        'lsc': 16,
        'elp': 18,
        'bpl': 19,
        'sqa': 19,
        'ogg': 22,
        }

for lsite,rhour in _restart_hr_utc.items():
    if lsite in info.keys():
        info[lsite]['restart_hour_utc'] = rhour
del lsite, rhour

##--------------------------------------------------------------------------##
## 1-meter calibration timing, expressed in sun altitude (deg):
_calib_sunalt = {}
_calib_sunalt['evening_autofocus'] = (-7.0, -12.0)
_calib_sunalt['morning_autofocus'] = (-10.0, -7.0)

##--------------------------------------------------------------------------##
## NRES config includes valid permutations of (camid, telescope, enclosure):
info['bpl']['nres_spec'] = None
info['coj']['nres_spec'] = None
info['cpt']['nres_spec'] = [
        {
                'cam' :   'nres03',
              'agcam' :   'ak06',
                'tel' :   '1m0a',
                'enc' :   'domb',
          'bias_spec' :  ( 5,   0.0),   # (count, exp_time_sec)
          'dark_spec' :  ( 3, 720.0),   # (count, exp_time_sec)
          'tung_spec' :  (10, 240.0),   # (count, exp_time_sec)
          'thar_spec' :  ( 5, 600.0),   # (count, exp_time_sec)
        },
        {
                'cam' :   'nres03',
              'agcam' :   'ak05',
                'tel' :   '1m0a',
                'enc' :   'domc',
          'bias_spec' :  ( 5,   0.0),   # (count, exp_time_sec)
          'dark_spec' :  ( 3, 720.0),   # (count, exp_time_sec)
          'tung_spec' :  (10, 240.0),   # (count, exp_time_sec)
          'thar_spec' :  ( 5, 600.0),   # (count, exp_time_sec)
        },
        ]
info['elp']['nres_spec'] = [
        {
                'cam' :   'nres02',
              'agcam' :   'ak04',
                'tel' :   '1m0a',
                'enc' :   'doma',
          'bias_spec' :  ( 5,   0.0),   # (count, exp_time_sec)
          'dark_spec' :  ( 3, 720.0),   # (count, exp_time_sec)
          'tung_spec' :  (10, 210.0),   # (count, exp_time_sec)
          'thar_spec' :  ( 5, 600.0),   # (count, exp_time_sec)
        },
        ]
info['lsc']['nres_spec'] = [
        {
                'cam' :   'nres01',
              'agcam' :   'ak01',
                'tel' :   '1m0a',
                'enc' :   'domb',
          'bias_spec' :  ( 5,   0.0),   # (count, exp_time_sec)
          'dark_spec' :  ( 3, 720.0),   # (count, exp_time_sec)
          'tung_spec' :  (10, 120.0),   # (count, exp_time_sec)
          'thar_spec' :  ( 5, 300.0),   # (count, exp_time_sec)
        },
        {
                'cam' :   'nres01',
              'agcam' :   'ak02',
                'tel' :   '1m0a',
                'enc' :   'domc',
          'bias_spec' :  ( 5,   0.0),   # (count, exp_time_sec)
          'dark_spec' :  ( 3, 720.0),   # (count, exp_time_sec)
          'tung_spec' :  (10,  35.0),   # (count, exp_time_sec)
          'thar_spec' :  ( 5, 180.0),   # (count, exp_time_sec)
        },
        ]
info['ogg']['nres_spec'] = None
info['sqa']['nres_spec'] = None
info['tfn']['nres_spec'] = None
#info['tlv']['nres_spec'] = None
info['tlv']['nres_spec'] = [
        {
                'cam' : 'nres04',
              'agcam' :   'ak03',
                'tel' :   '1m0a',
                'enc' :   'doma',
          'bias_spec' :  (5,   0.0),    # (count, exp_time_sec)
          'dark_spec' :  (3, 720.0),    # (count, exp_time_sec)
          'tung_spec' :  (5, 105.0),    # (count, exp_time_sec)
          'thar_spec' :  (5, 300.0),    # (count, exp_time_sec)
        },
        ]


##--------------------------------------------------------------------------##
## NRES validation and look-up class:
class NRES_Check(object):

    # Basic knowledge:
    def __init__(self):
        self.nres_sites = [x for x in info.keys() \
                if info[x]['nres_spec'] != None]
        self.agcam_info, self.nrcam_info = self._enumerate_nspecs()
        pass

    # ------------------------------
    # Setup helpers:
    @staticmethod
    def _tes_token(tel, enc, site):
        return '.'.join((tel, enc, site))

    def _enumerate_nspecs(self):
        agu_list, cam_list, seq_list = [], [], []
        for site in self.nres_sites:
            for nn in info[site]['nres_spec']:
                token = '.'.join((nn['tel'], nn['enc'], site))
                seq_list.append(token)
                agu_list.append(nn['agcam'])
                cam_list.append(nn['cam'])
        agu_map = dict(zip(seq_list, agu_list))
        cam_map = dict(zip(seq_list, cam_list))
        return agu_map, cam_map

    # ------------------------------
    def _site_has_nres(self, site):
        return (site in self.nres_sites)

    def _fetch_nspecs(self, site, encl, tel='1m0a'):
        hits = []
        if self._site_has_nres(site):
            hits.extend([x for x in info[site]['nres_spec'] \
                    if (x['enc']==encl and x['tel']==tel)])
        return hits

    # ------------------------------
    def is_valid_sequencer(self, token):
        return True if (token in self.agcam_info.keys()) else False

    def is_valid_nres(self, site, encl, tel='1m0a'):
        token = '.'.join((tel, encl, site))
        return self.is_valid_sequencer(token)

    # ------------------------------
    def agcam_lookup(self, site, encl, tel='1m0a'):
        token = '.'.join((tel, encl, site))
        return self.agcam_info.get(token, None)

    def seq_to_agcam(self, sequencer):
        return self.agcam_info.get(sequencer, None)

    # ------------------------------
    def nrcam_lookup(self, site, encl, tel='1m0a'):
        token = '.'.join((tel, encl, site))
        return self.nrcam_info.get(token, None)

    def seq_to_nrcam(self, sequencer):
        return self.nrcam_info.get(sequencer, None)

    pass

##--------------------------------------------------------------------------##




######################################################################
# CHANGELOG (lco_site_info.py):
#---------------------------------------------------------------------
#
#  2018-07-03:
#     -- Added TLV to site list.
#
#  2018-02-15:
#     -- ALMOST: changed out-of-date LSC restart time from 16h to 17h UTC.
#           Confusion about the precise timings of events remains.
#
#  2018-01-18:
#     -- Increased __version__ to 0.5.1.
#     -- Added missing 'aqwb' to LSC enclosure list.
#
#  2017-12-18:
#     -- Increased __version__ to 0.5.0.
#     -- Added NRES_Check() class with look-up and validation utilities.
#     -- Increased __version__ to 0.4.4.
#     -- Set BPL nres_spec to None. Moved old nres_spec to TLV.
#
#  2017-12-04:
#     -- Increased __version__ to 0.4.3.
#     -- Added igla to ELP, CPT enclosure lists. Added aqwa to CPT.
#     -- Increased __version__ to 0.4.2.
#     -- CPT now restarts at 11h UTC so that restart matches day-obs roll.
#
#  2017-10-16:
#     -- Increased __version__ to 0.4.1.
#     -- Added aqwa to elp enclosure list.
#
#  2017-10-12:
#     -- Increased __version__ to 0.4.0.
#     -- Added per-site enclosure list (consolidated from POND-URLs.py).
#
#  2017-10-11:
#     -- Increased __version__ to 0.3.4.
#     -- Switched ELP site restart time back to 18h. The site restart was
#           switched to 17h when it was rigorously bound to a UTC hour.
#           This created a mis-match between site restart and day-obs rollover
#           which led to unexpected problems. The site restart has been moved
#           back to 18h which hopefully resolves those issues.
#
#  2017-09-27:
#     -- Increased __version__ to 0.3.3.
#     -- Added nres_spec for ELP.
#
#  2017-09-14:
#     -- Increased __version__ to 0.3.2.
#     -- Added OGG site restart hour to info collection.
#
#  2017-09-09:
#     -- Increased __version__ to 0.3.1.
#     -- Changed NRES AGU camera at BPL from ak08 to ak07.
#     -- BPL now hosts NRES4.
#
#  2017-08-04:
#     -- Increased __version__ to 0.3.0.
#     -- Added sun altitudes marking the start/end of autofocus windows.
#
#  2017-07-26:
#     -- Increased __version__ to 0.2.8.
#     -- Increased bpl thar01 exposure time from 180 to 300 seconds.
#     -- Increased bpl tung01 exposure time from 35 to 105 seconds.
#
#  2017-07-15:
#     -- Increased __version__ to 0.2.7.
#     -- Increased domc.lsc tung exposures from 5 to 10.
#
#  2017-07-11:
#     -- Increased __version__ to 0.2.6.
#     -- Increased domb.lsc tung exposures from 5 to 10.
#
#  2017-07-10:
#     -- Increased __version__ to 0.2.5.
#     -- Added bias_spec and dark_spec to NRES configs.
#
#  2017-07-06:
#     -- Increased __version__ to 0.2.0.
#     -- Added tung_spec and thar_spec for each NRES instance.
#
#  2017-05-22:
#     -- Increased __version__ to 0.1.0.
#     -- First created lco_site_info.py.
#
