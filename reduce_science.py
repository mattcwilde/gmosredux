#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
An `in-folder` reduction script CGM^2 Multi-Object Spectra

Note: You should really be running the prepare file to get this into the correct format...
"""

__author__ = "Matthew Wilde <mwilde@uw.edu>"
__version__ = "March, 2018"

try:
    from pyraf import iraf

except ImportError:
    raise ImportError("No pyraf module -- did you forget to do `source activate geminiconda`?")

import os
import copy
from pyraf.iraf import gemini, gmos, gemtools, onedspec
import numpy as np
from astropy.io import fits
from astropy.table import Table
import glob
import matplotlib.pyplot as plt
from gmosutils import gmosutils


path = os.path.abspath('.')

# Set colormap for plotting
cm = plt.get_cmap('Greys')

flatName = 'MC' + 'gcalFlat'
combName = 'MC' + 'gcalFlat' + 'Comb'
biasName = 'MCbiasFull.fits'


def observation_summary(filenames, additional_headers=None):
    """
    Create a table to serve as an observing summary by reading in the fits image headers.

    :param filenames: list of filenames (usually from glob)
    :param additional_headers: fits headers in case you need them
    :return: returns a astropy table of the observations
    """

    # List the headers we want to extract
    headers = ["INSTRUME", "OBJECT", "OBSTYPE", "MASKNAME", "OBSCLASS",
               "CENTWAVE", "GEMPRGID", "OBSID", "CCDBIN", "RAWPIREQ",
               "DATALAB", "UT", "DATE", "TIME-OBS", "GRATING", "EXPOSURE", "DETECTOR", "DATE-OBS"]

    if additional_headers is not None:
        headers += additional_headers

    # read in each header paramater from each file
    rows = []
    for filename in filenames:
        with fits.open(filename) as image:
            rows.append(dict(zip(headers,
                                 [image[0].header.get(h, None) for h in headers])))

        # Add the filename.
        rows[-1]["FILENAME"] = filename

    headers.insert(0, "FILENAME")

    return Table(rows=rows, names=headers)


def show_img(img, z1, z2):
    """
    Displays a fits image
    :param img: fits image
    :param z1:  vmin
    :param z2:  vmax
    :return:
    """
    plt.clf()
    ax = plt.gca()
    ax.imshow(img, vmin=z1, vmax=z2, cmap=cm, origin='lower')
    #
    plt.show()



# find the files
raw_files = glob.glob('*fits')
summary = gmosutils.observation_summary(raw_files)
cent_waves = list(set(summary[(summary['OBSTYPE'] == 'OBJECT')]['CENTWAVE']))

print ("### Begin Processing GMOS/MOS Spectra ###")
print (' ')
print('IN FOLDER: ',path)
print (' ')
print ("=== Creating FLAT MasterCals ===")
if not os.path.exists(flatName+str(cent_waves[1])+'.fits'):
    print (" -- Creating GCAL Spectral Flat-Field MasterCals --")
    print ("  -Full Flat (GCAL) normalization, non-interactive-")

    # Set the task parameters.
    gmos.gireduce.unlearn()
    gmos.gsflat.unlearn()
    flatFlags = {'fl_over':'yes','fl_trim':'yes','fl_bias':'yes','fl_dark':'no',
                'fl_fixpix':'no','fl_oversize':'no','fl_vardq':'yes','fl_fulldq':'yes',
                'fl_inter':'no','fl_detec':'yes',
                'function':'spline3','order':'8',
                'logfile':'gsflatLog.txt','verbose':'no'
    }


    flatFlags.update({'fl_keep':'yes','fl_usegrad':'yes','fl_detec':'no',
                     'fl_seprows':'no','order':53})

    # flat type
    ft = 'gcalFlat'
    # loop over central wavelengths
    flatName=''
    combName=''
    for cw in cent_waves:
        flatName = 'MC' + ft + str(cw)
        combName = 'MC' + ft + 'Comb' + str(cw)
        flats = summary[(summary['OBSTYPE'] == 'FLAT') & (summary['CENTWAVE'] == cw)]['FILENAME']
        flatFull = [os.path.basename(f) for f in flats]
        gmos.gsflat (','.join(str(x) for x in flatFull), flatName,
                     bias='MCbiasFull.fits', combflat=combName, **flatFlags)
else:
    print('SKIPPING FLAT REDUX: '+flatName+' already exists')