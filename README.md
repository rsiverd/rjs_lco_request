# rjs_lco_request
A few Python scripts to submit LCO observing requests (just NRES right now)
through the Valhalla API. Objects are taken from a target list specified on
the command line. A few basic functions to pick targets, specify the request
window, pick a site, etc. are included.

Written by Rob Siverd.

Installation
============

0) Clone the repository. Consider copying files to a separate location.
```
$ git clone https://github.com/rsiverd/rjs_lco_request
```

1) Make sure you have the required dependencies installed. You will need
requests, pyephem, astropy. The astroquery package will be needed if you want
to use the `add_target.py` script. I have included a requirements.txt file with
the versions of modules that I have installed (with which this script has been
tested). It will likely work with other versions, I have just provided this
file for reference. To add these specific packages to e.g., a new virtual env
for running this script, try 
```
$ pip install -r requirements.txt
```

2) Add your own auth token to credentials.txt. You can find it on your profile
page in the observing portal at `https://observe.lco.global/accounts/profile/`.

3) Edit request_NRES_spectra.py and update the key project `proposal_id` 
as appropriate for your submissions. Currently this script is set to use Tim
Brown's NRES key project proposal.

4) Make a list of targets formatted like the sample provided in the repository.
The `add_target.py` script is meant to retrieve most of the required fields
from SIMBAD. The functional state of this script is not currently known, please
check results in case it needs an update (and let me know if it does).

5) Submit things! Example submission command:
```
./request_NRES_spectra.py -n1 -SK --any --SPECTRUM --pystrat -t ./target_list_nres.master.txt --Vmax=7 -v
```

