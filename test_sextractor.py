#!/usr/bin/env python
# coding: utf-8


import sys, os
import sewpy
import numpy as np
from astropy.table import Table
from astropy.time import Time, TimeDelta
from astropy.io import fits
import glob

# glob.glob("*")
def find_closest_zpt(table, date):
    """ This function finds the zeropoint mag from 
        Gemini-N zpt tables to our image based on the 
        date the image was observed. 
        
        input
        ------------
        table : astropy.Table 
        date : date (string)
        
        returns
        ------------
        zp : np.float64 : the zeropoint 
    """
    
    zpt = table['ZP'][np.argmin(abs(Time(table['DATE'])-Time(date)))]
    return zpt


def load_ztp_table(band='i', tbl_dir='./'):
    """
        Read in the (changing) zeropoint information directly from the gemini website. 
        
        NOTE: (TODO:) Only for GMOS-N
        
        Args:
            band (str): one of the ugriz filter bands
        
        Returns:
            zpt_tbl_North (astropy.table.Table) : data table containing dates and zeropoints
    """
    

    if band in ['g', 'r', 'i', 'z']:
        raw_tbl_name = 'GMOSphotometryMean_{}.dat'.format(band)
        file = tbl_dir+raw_tbl_name
        if os.path.exists(file):
            zpt_tbl_North = Table.read(file, format='ascii')
        else:            
            url = 'http://www.gemini.edu/sciops/instruments/ipm/data-products/gmos-n-and-s/'+raw_tbl_name
            zpt_tbl_North = Table.read(url, format='ascii')
            zpt_tbl_North.write(raw_tbl_name, format='ascii')
    else:
        raise ValueError('Band must be one of: u, g, r, i, or z')
    
    # rename colm to useful names
    zpt_tbl_North.rename_column('col1','DATE')
    zpt_tbl_North.rename_column('col2','ZP')
    
    return zpt_tbl_North



def calc_Mag_zpt(mzero, exptime, kMK, airmass):
    """ Calculates the zeropoint mag to use for input to sextractor.
    
        Gemini gives:
            mstd = mzero - 2.5 log10 (N[-e]/exptime) - kMK (airmass-1.0)
            
        Jess says:
            mzero + 2.5*np.log10(exptime) - kMK(airmass)
        
    """
    
    mag_zpt = mzero + 2.5*np.log10(exptime) - kMK(airmass)
    
    return mag_zpt


def get_filter(image_file):
    """ Read in the filter used on the instrument and get its simple 
        name.
    """
    with fits.open(image_file) as hdu:
        prime_header = hdu['PRIMARY'].header
        sci_header = hdu['SCI'].header
        
        filter2 = prime_header['FILTER2']
        if 'i' in filter2:
            filt = 'i'
        elif 'g' in filter2:
            filt = 'g'
        elif 'r' in filter2:
            filt = 'r'
        elif 'z' in filter2:
            filt = 'z'
        else:
            raise ValueError('You must not have input an image file for GMOS...?')
            
    return filt


def get_filter(image_file):
    """ Read in the filter used on the instrument and get its simple 
        name.
        
    Args:
        image_file (str): the abspath for the image file name
        
    Returns:
        filt (str): stripped human readable filter 
    """
    with fits.open(image_file) as hdu:
        prime_header = hdu['PRIMARY'].header
        sci_header = hdu['SCI'].header
        
        filter2 = prime_header['FILTER2']
        if 'i' in filter2:
            filt = 'i'
        elif 'g' in filter2:
            filt = 'g'
        elif 'r' in filter2:
            filt = 'r'
        elif 'z' in filter2:
            filt = 'z'
        else:
            raise ValueError('You must not have input an image file for GMOS...?')
            
    return filt


def calc_Mag_zpt(mzero, exptime, kMK, airmass):
    """ Calculates the zeropoint mag to use for input to sextractor.
    
    Gemini gives:
        mstd = mzero - 2.5 log10 (N[-e]/exptime) - kMK (airmass-1.0)

    Jess says:
        mzero - 2.5*np.log10(exptime) - kMK(airmass)
        
    Args:
        mzero (str): the measured zpt mag from gemini
        exptime (str): exporsure time (from image header)
        kMK (float): color correction term
        airmass (str): airmass (from image header)
    
    Returns:
        mag_zpt (str): MAG_ZEROPOINT term for sextractor
    """
    mag_zpt = float(mzero) + 2.5*np.log10(float(exptime)) - float(kMK)*float(airmass)
    
    return mag_zpt


def get_color_correction(filt):
    """ Load the filter specific color corrections from:
        
    https://www.gemini.edu/sciops/instruments/gmos/calibration/photometric-stds
        
    Args:
        filt (str): the human readable filter
    
    Returns:
        kMK (float): the color correction term
        
        
    """

    if filt == 'g': 
        kMK = 0.14
    if filt == 'r':
        kMK = 0.11
    if filt == 'i':
        kMK = 0.10
    if filt == 'z':
        kMK = 0.05
    return kMK






# working_path = '/Users/mwilde/Dropbox/COS-Gemini/gemini_data.GN-2014A-Q-1/reduced_data/J1555+3628/mask_11/Imaging/'
working_path = './'
param_file = working_path+'geminiimaging.param'
image_file = working_path+'mrgN20140327S0155_add.fits'

with fits.open(image_file) as hdu:
    prime_header = hdu['PRIMARY'].header
    sci_header = hdu['SCI'].header


# # Load in the zero-point table from gemini
filt = get_filter(image_file)
zpt_tbl_North = load_ztp_table(filt)


# # Load in param file
# 
# These are the paramters that are output by sextractor that jess saved to output previously. I'll just reuse them.
params = list(np.loadtxt(param_file, dtype=str))


# # load the headers from the image file
# 
# these give us specific values needed to calculate zeropoint for example but maybe dont need this yet
def get_detector_params(prime_header, table, date):
    """ work in progress...
    """
    if 'N' in prime_header['INSTRUME']:
        North = True
        mzero = find_closest_zpt(table, date)
        if 'GMOS + e2v DD' in prime_header['DETECTOR']:
            print('Make sure to use the GMOS-N e2v DD params!')
                      
# get_detector_params()


# #  pull the needed params out of the header

# These are the parameters that sextractor needs specific to the instrument.
airmass = prime_header['AIRMASS']
exptime = prime_header['EXPTIME']
date = prime_header['DATE-OBS']
pixel_scale = prime_header['PIXSCALE']
seeing = sci_header['FWHMPSF']
cat_name = prime_header['OBJECT']+'sex_gem_i_pythontest.fits'
gain = sci_header['GAIN']


# # Find the best zeropoint from the online table
mzero = find_closest_zpt(zpt_tbl_North, date)


# # Load in color corrects and maybe avg zero-points
kMK = get_color_correction(filt)


# # Zpt function:
mag_zpt = calc_Mag_zpt(mzero, exptime, kMK, airmass)


# # Try `sewpy`
# working_path = '/Users/mwilde/Dropbox/COS-Gemini/gemini_data.GN-2014A-Q-1/reduced_data/J1555+3628/mask_11/Imaging/'
working_path = './'

import sewpy

config = {"CATALOG_NAME":cat_name, "SATUR_LEVEL": 50000, "MAG_ZEROPOINT":mag_zpt, "GAIN":gain,
         "PIXEL_SCALE":pixel_scale, "SEEING_FWHM":seeing}

# instatiate SE wrapper object
sew = sewpy.SEW(params=params, config=config, sexpath="sex", workdir=working_path)

sew.tmp = False
# set teh working dir ... gets weird otherwise
sew.workdir = working_path

# run sextractor on the image
out = sew('mrgN20140327S0155_add.fits')

print(out["table"])



