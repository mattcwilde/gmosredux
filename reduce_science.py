#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
An `in-folder` reduction script CGM^2 Multi-Object Spectra

Note: run the prepare file to get this into the correct format...
"""
from __future__ import print_function, division

__author__ = "Matthew Wilde <mwilde@uw.edu>"
__version__ = "April 10, 2018"

try:
    from pyraf import iraf

except ImportError:
    raise ImportError("No pyraf module -- did you forget to do `source activate geminiconda`?")

import copy
import glob
import os
import shutil
import argparse

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from pyraf.iraf import gemini, gmos, gemtools, task

# this task allows us to skip bad slits in gswavelenght!
iraf.task(mgswavelength='/Users/mwilde/iraf/mgswavelength.cl')

path = os.path.abspath('.')

# Set colormap for plotting
cm = plt.get_cmap('Greys')

flatPrefix = 'MC' + 'gcalFlat'
combName = 'MC' + 'gcalFlat' + 'Comb'
biasName = 'MCbiasFull.fits'

# Read in user defined slit(s) to skip with mgswavelength
parser = argparse.ArgumentParser(description="run the iraf reduction on cgmsquared data.")
parser.add_argument("--skip_slit", type=int, nargs='+', metavar='N',
                    help="slit extension to skip")
parser.add_argument("--use_jess", action="store_true", help=" use jess's gswavelength params (lower order fitting")
args = parser.parse_args()



# set some flag variables
if args.skip_slit:
    skip_slit = True
    slits_to_skip = []
    for slit in args.skip_slit:
        slits_to_skip.append(slit)
        print("********* NOTE:Skipping slit {} *********".format(slit))
    slits_to_skip = list(set(slits_to_skip))
    slits_to_skip.sort()
    num_slits = len(slits_to_skip)
else:
    skip_slit = False

if args.use_jess:
    use_jess = True
    print("********* NOTE: using Jess' lower order wave fit")
else: 
    use_jess = False

def observation_summary(filenames, additional_headers=None):
    """
    Create a table to serve as an observing summary.
    """

    # List the headers we want to extract
    headers = ["INSTRUME", "OBJECT", "OBSTYPE", "MASKNAME", "OBSCLASS",
               "CENTWAVE", "GEMPRGID", "OBSID", "CCDBIN", "RAWPIREQ",
               "DATALAB", "UT", "DATE", "TIME-OBS", "GRATING", "EXPOSURE", "DETECTOR", "DATE-OBS"]

    if additional_headers is not None:
        headers += additional_headers

    rows = []
    for filename in filenames:
        with fits.open(filename) as image:
            rows.append(dict(zip(headers,
                                 [image[0].header.get(h, None) for h in headers])))

        # Add the filename.
        rows[-1]["FILENAME"] = filename

    headers.insert(0, "FILENAME")
    # headers.append("TIME")

    return Table(rows=rows, names=headers)


def show_img(img, z1, z2):
    plt.clf()
    ax = plt.gca()
    ax.imshow(img, vmin=z1, vmax=z2, cmap=cm, origin='lower')
    #
    plt.show()


def get_max_extension(file):
    """
    Get the last fits extension of the reduced images from the MDF

    -------
    returns: max_ext (int)
    """


    with fits.open(file) as hdu:
        mdf = hdu['MDF'].data
        return mdf['EXTVER'].max()


def delete_tmp_files():
    """Get rid of the intermediate gsreduced files iraf produces

    """
    intermediate_files = glob.glob('*tmp*')
    for file in intermediate_files:
        try:
            os.remove(file)
        except OSError:
            pass


def delete_sci_files():
    """Get rid of the intermediate gsreduced files iraf produces

    """
    intermediate_files = glob.glob('J*fits')
    for file in intermediate_files:
        try:
            os.remove(file)
        except OSError:
            pass


def clean_up():
    """
    Get rid of the all the debris created by previous reduction attempts.
    This should get the folder back to only having the necessary raw data


    """
    intermediate_files = glob.glob('*g*fits')
    intermediate_files += (glob.glob('*.lst'))
    intermediate_files += (glob.glob('*.txt'))
    intermediate_files += (glob.glob('*.log'))
    intermediate_files += (glob.glob('f*.fits'))
    intermediate_files += (glob.glob('p*.fits'))
    intermediate_files += (glob.glob('[t|u|e]*.fits'))
    intermediate_files += (glob.glob('*sci*.fits'))
    intermediate_files += (glob.glob('*.cl'))
    intermediate_files += (glob.glob('*g[S|N]*fits'))
    intermediate_files += (glob.glob('*gs[S|N]*fits'))
    for file in intermediate_files:
        try:
            os.remove(file)
        except OSError:
            pass

    try:
        shutil.rmtree('database/')
    except OSError:
        pass


def get_bias():
    """
    This function tries to look for the bias based on my usual file convention.
    If not it will then look for the local copy before telling the user there is not one to be had.

    It will try and copy a bias file to the folder.

    :return: None
    """

    # try to copy a fresh bias over first
    try:
        shutil.copy('../../../Bias/MCbiasFull.fits', './MCbiasFull.fits')
        print("copying MCbiasFull from from ../../../Bias/")
    except OSError:
        if not os.path.exists('MCbiasFull.fits'):
            print('No bias? maybe you havent made one yet?')
        else:
            print("bias already exists, ready to roll!")


def make_flat(summary, cent_waves):
    """
    This function creates the Master Calibration flats used to reduce Gemini MOS spectra.

    It will create 2 difference flat images:
        1) a normalized flat field image used to be removed from the science images

        2) a non-normalized (or combined if there are more than 1 flat image per central wavelength)
            used to find the slit edges.

    :param summary: astroypy Table object with header data in it for each file
    :param cent_waves: list of floats, central wavelengths from summary
    :return: None
    """
    print("### Begin Processing GMOS/MOS Spectra ###")
    print(' ')
    print('IN FOLDER: ', path)
    print(' ')
    print("=== Creating FLAT MasterCals ===")
    

    pwd = os.path.abspath(".")
    if not os.path.exists(pwd+'/'+flatPrefix + str(cent_waves[1])[:-2] + '.fits'):
        print(" -- Creating GCAL Spectral Flat-Field MasterCals --")
        print("  -Full Flat (GCAL) normalization, non-interactive-")

        # Set the task parameters.
        gmos.gireduce.unlearn()
        gmos.gsflat.unlearn()
        flatFlags = {'fl_over': 'yes', 'fl_trim': 'yes', 'fl_bias': 'yes', 'fl_dark': 'no',
                     'fl_fixpix': 'no', 'fl_oversize': 'no', 'fl_vardq': 'yes', 'fl_fulldq': 'yes',
                     'fl_inter': 'no', 'fl_detec': 'yes',
                     'function': 'spline3', 'order': '8',
                     'logfile': 'gsflatLog.txt', 'verbose': 'no'
                     }


        flatFlags.update({'fl_keep': 'yes', 'fl_usegrad': 'yes', 'fl_detec': 'no',
                          'fl_seprows': 'no', 'order': 53})

        # flat type
        ft = 'gcalFlat'
        # loop over central wavelengths
        # flatName = ''
        # combName = ''
        for cw in cent_waves:
            flatName = 'MC' + ft + str(cw)[:-2]
            combName = 'MC' + ft + 'Comb' + str(cw)[:-2]
            flats = summary[(summary['OBSTYPE'] == 'FLAT') & (summary['CENTWAVE'] == cw)]['FILENAME']
            flatFull = [os.path.basename(f) for f in flats]
            gmos.gsflat(','.join(str(x) for x in flatFull), flatName,
                        bias='MCbiasFull.fits', combflat=combName, **flatFlags)
    else:
        print('SKIPPING FLAT REDUX: ' + flatPrefix + ' already exists')
    return None


def reduce_science(summary, cent_waves):
    """
    Reduce the science and arc images and apply wavelength transformation
    this produces non sky subtracted, non-combined, transformed 2D slits to then
    be loaded into PYPIT.

    :param summary: astroypy Table object with header data in it for each file
    :param cent_waves: list of floats, central wavelengths from summary
    :return: None
    """

    print("=== Processing Science Files ===")
    print(" -- Performing Basic Processing --")

    # Use primarily the default task parameters.
    gmos.gsreduce.unlearn()
    sciFlags = {
        'fl_over': 'yes', 'fl_trim': 'yes', 'fl_bias': 'yes', 'fl_gscrrej': 'no',
        'fl_dark': 'no', 'fl_flat': 'yes', 'fl_gmosaic': 'yes', 'fl_fixpix': 'no',
        'fl_gsappwave': 'yes', 'fl_oversize': 'no',
        'fl_vardq': 'yes', 'fl_fulldq': 'yes',
        'fl_inter': 'no', 'logfile': 'gsreduceLog.txt', 'verbose': 'no'
    }
    arcFlags = copy.deepcopy(sciFlags)
    arcFlags.update({'fl_flat': 'no', 'fl_vardq': 'no', 'fl_fulldq': 'no'})

    gmos.gswavelength.unlearn()
    # USE JESS' PARAMS
    if use_jess:
        
        waveFlags = {
            'function':'chebyshev','order':4,
            'fl_inter':'no','logfile':'gswaveLog.txt','verbose':'yes', 'step':2, 'nlost':10
        }
    else:
        waveFlags = {
            'coordlist': 'gmos$data/CuAr_GMOS.dat', 'fwidth': 6, 'nsum': 50,
            'function': 'chebyshev', 'order': 7,
            'fl_inter': 'no', 'logfile': 'gswaveLog.txt', 'verbose': 'no'
        }
    waveFlags.update({'order': 7, 'nsum': 20, 'step': 2})



    gmos.gstransform.unlearn()
    transFlags = {
        'fl_vardq': 'yes', 'interptype': 'linear', 'fl_flux': 'yes',
        'logfile': 'gstransformLog.txt', 'verbose': 'no'
    }

    print('  - MOS Science and Arc exposures -')
    prefix = 'gs'  #
    ft = 'gcalFlat'  # flat type
    for i, cw in enumerate(cent_waves):
        flatName = 'MC' + ft + str(cw)[:-2]
        gradName = 'MC' + ft + 'Comb' + str(cw)[:-2]

        ################### ARCS: Reduce and find wavelengths ###############################
        # Arcs
        # arcFull = fs.fileListQuery(dbFile, fs.createQuery('arc', qdf), qdf)
        arcs = summary[(summary['OBSTYPE'] == 'ARC') & (summary['CENTWAVE'] == cw)]['FILENAME']
        arcFull = [os.path.basename(f) for f in arcs]
        reduced_arcs = [prefix + str(x) for x in arcFull]
        if not os.path.exists(reduced_arcs[0]):
            gmos.gsreduce(','.join(str(x) for x in arcFull), bias=biasName, gradimage=gradName, **arcFlags)

            # find max extenstion number:
            max_ext = get_max_extension(reduced_arcs[0])

            # import pdb; pdb.set_trace()

            # IF YOU WANT TO SKIP A SLIT USE THESE LINES!
            if skip_slit:
                if num_slits == 1:
                    if slits_to_skip[0] == 1:
                        iraf.mgswavelength(','.join(prefix + str(x) for x in arcFull), firstsciext=2, lastsciext=max_ext, **waveFlags)
                    else:
                        iraf.mgswavelength(','.join(prefix + str(x) for x in arcFull), firstsciext=1, lastsciext=slits_to_skip[0]-1, **waveFlags)
                        iraf.mgswavelength(','.join(prefix + str(x) for x in arcFull), firstsciext=slits_to_skip[0]+1, lastsciext=max_ext, **waveFlags)
                elif num_slits == 2:
                    # check they arent sequential   
                    if slits_to_skip[1]-slits_to_skip[0] != 1: 
                        if slits_to_skip[0] == 1:
                            iraf.mgswavelength(','.join(prefix + str(x) for x in arcFull), firstsciext=2, lastsciext=slits_to_skip[1] - 1, **waveFlags)
                            iraf.mgswavelength(','.join(prefix + str(x) for x in arcFull), firstsciext=slits_to_skip[1]+1, lastsciext=max_ext, **waveFlags)
                        else:
                            iraf.mgswavelength(','.join(prefix + str(x) for x in arcFull), firstsciext=1, lastsciext=slits_to_skip[0]-1, **waveFlags)
                            iraf.mgswavelength(','.join(prefix + str(x) for x in arcFull), firstsciext=slits_to_skip[0]+1, lastsciext=slits_to_skip[1]-1, **waveFlags)
                            iraf.mgswavelength(','.join(prefix + str(x) for x in arcFull), firstsciext=slits_to_skip[1]+1, lastsciext=max_ext, **waveFlags)
                    # if the slits ARE sequential
                    else:
                        if slits_to_skip[0] == 1:
                            iraf.mgswavelength(','.join(prefix + str(x) for x in arcFull), firstsciext=slits_to_skip[1]+1, lastsciext=ax_ext, **waveFlags)
                        else:
                            iraf.mgswavelength(','.join(prefix + str(x) for x in arcFull), firstsciext=1, lastsciext=slits_to_skip[0]-1, **waveFlags)
                            iraf.mgswavelength(','.join(prefix + str(x) for x in arcFull), firstsciext=slits_to_skip[1]+1, lastsciext=max_ext, **waveFlags)
                else:
                    print("you input too many slits")

            else:
                gmos.gswavelength(','.join(prefix + str(x) for x in arcFull), **waveFlags)
            

        ################### SCI: Reduce and find wavelengths ###############################
        # Science Images
        # sciFull = fs.fileListQuery(dbFile, fs.createQuery('sciSpec', qdf), qdf)
        sci = summary[(summary['OBSTYPE'] == 'OBJECT') & (summary['CENTWAVE'] == cw)]['FILENAME']
        sciFull = [os.path.basename(f) for f in sci]

        reduced_imgs = [prefix + str(x) for x in sciFull]

        # get the quasar name convention right: quasar+mask+centwave
        quasar_name = summary[(summary['OBSTYPE'] == 'OBJECT') & (summary['CENTWAVE'] == cw)]['OBJECT'][0]
        maskname = summary[(summary['OBSTYPE'] == 'OBJECT') & (summary['CENTWAVE'] == cw)]['MASKNAME'][0][-2:]
        outFile = quasar_name + '_m' + maskname + '_' + str(cw)[:-2]
        outFile = outFile.replace('+', '_')

        gmos.gsreduce.unlearn()

        # Rollup: 2 central waves and 2 exposures at each.
        # This means we can combine the exposure at the same central wavelength to perform
        # cosmic ray reduction
        if len(cent_waves) < 3:

            # reduce the SCIENCE images
            if not os.path.exists(reduced_imgs[0]):
                gmos.gsreduce(','.join(str(x) for x in sciFull), bias=biasName, flatim=flatName, gradimage=gradName,
                              **sciFlags)

            # Set CR task parameters.
            # NOTE: Turned on crreject. Does a good job of getting rid of cosmic rays
            gemtools.gemcombine.unlearn()
            sciCombFlags = {
                'combine': 'average', 'reject': 'crreject',
                'fl_vardq': 'yes', 'fl_dqprop': 'yes',
                'logfile': 'gemcombineLog.txt', 'verbose': 'no'
            }
            # Rollup data with 2 science exposures and 2 centwaves
            i = i + 1
            combined_outFile = "cr_sci" + str(i) + "w_auto_slit.fits"

            # combine the 2 exposures at the same central wavelength
            gemtools.gemcombine(','.join(prefix + str(x) for x in sciFull), combined_outFile,
                                **sciCombFlags)
            # transform the sci images
            gmos.gstransform(combined_outFile, wavtraname='gs' + arcFull[0], outimages=outFile,
                             **transFlags)
        else:
            # This is the Long Program data with 3 central wavelength with 1 exposure at each.
            # Cosmic ray reduction needs to be performed
            # can be done with fl_crspec=yes in gs reduce...

            sciFlags = {
                'fl_over': 'yes', 'fl_trim': 'yes', 'fl_bias': 'yes', 'fl_gscrrej': 'yes',
                'fl_dark': 'no', 'fl_flat': 'yes', 'fl_gmosaic': 'yes', 'fl_fixpix': 'no',
                'fl_gsappwave': 'yes', 'fl_oversize': 'no',
                'fl_vardq': 'yes', 'fl_fulldq': 'yes',
                'fl_inter': 'no', 'logfile': 'gsreduceLog.txt', 'verbose': 'no'
            }

            # reduce the SCIENCE images
            if not os.path.exists(reduced_imgs[0]):
                gmos.gsreduce(','.join(str(x) for x in sciFull), bias=biasName, flatim=flatName, gradimage=gradName,
                              **sciFlags)
                # LP with 3 central wavelengths and 1 science exposure
            else:
                print("skipping gsreduce: file already exists " + reduced_imgs[0] + ".fits")
            # get the quasar name convention right: quasar+mask+centwave
            quasar_name = summary[(summary['OBSTYPE'] == 'OBJECT') & (summary['CENTWAVE'] == cw)]['OBJECT'][0]
            maskname = summary[(summary['OBSTYPE'] == 'OBJECT') & (summary['CENTWAVE'] == cw)]['MASKNAME'][0][-2:]
            outFile = quasar_name + '_m' + maskname + '_' + str(cw)[:-2]
            outFile = outFile.replace('+', '_')

            if not os.path.exists(outFile + ".fits"):
                # transform the sci images
                gmos.gstransform(reduced_imgs[0], wavtraname='gs' + arcFull[0], outimages=outFile,
                                 **transFlags)
            else:
                print("skipping gstransform: file already exists" + outFile + ".fits")

    return None


def check_if_reduced():
    """ checks to see if the output science file already in the target folder

    """
    pwd = os.path.abspath(".")
    outpath = '/Users/mwilde/Dropbox/COS-Gemini/2Dspecred/'

    project_IDs = ["GN-2014A-Q-1", "GN-2014B-LP-3", "GN-2015A-LP-3",
                   "GS-2014A-Q-2", "GS-2014B-LP-4", "GS-2015A-LP-4"]

    # figure out which project we are working on
    for pid in project_IDs:
        if pid in pwd:
            outfolder = pid + '/'

            # if the folder doesnt exist, make it
            if not os.path.exists(outpath + outfolder):
                os.mkdir(outpath + outfolder)
        else:
            print('project id not in list')
            break

    # copy final 2Dspecred images to that folder
    final_outfiles = glob.glob('J*.fits')
    if len(final_outfiles) > 1:
        pass

    for outfile in final_outfiles:
        if not os.path.exists(outpath + outfolder + outfile):
            shutil.copy(outfile, outpath + outfolder + outfile)
        else:
            print("{} already in 2Dspecred, you've wasted your time!".format(outpath + outfolder + outfile))


def move_2Dspecred(summary):
    """ Move the finished output science files to their final destination

    """

    outpath = '/Users/mwilde/Dropbox/COS-Gemini/2Dspecred/'

    # get the project ID name
    pid = list(set(summary['GEMPRGID']))[0]

    # if the project ID folder doesnt exist, make it
    outfolder = pid + '/'
    if not os.path.exists(outpath + outfolder):
        os.mkdir(outpath + outfolder)

    # copy final 2D reduced spectra to that folder
    final_outfiles = glob.glob('J*.fits')
    for outfile in final_outfiles:
        if not os.path.exists(outpath + outfolder + outfile):
            shutil.copy(outfile, outpath + outfolder + outfile)
        else:
            print("{} already in 2Dspecred, you've wasted your time!".format(outpath + outfolder + outfile))


def make_2Dspec(start_from_scratch=False, clean=False):
    """

    :return:
    """
    # load in the files
    raw_files = glob.glob('[N|S]*fits')

    # make a nice table out of them
    summary = observation_summary(raw_files)

    # find the central wavelengths we are gonna use 
    cent_waves = list(set(summary[(summary['OBSTYPE'] == 'OBJECT')]['CENTWAVE']))


    # copy the bias over or use the one here
    get_bias()
    # shutil.copy('../../../Bias/MCbiasFull.fits', './MCbiasFull.fits')

    # Use to redo the redux
    if start_from_scratch:
        delete_sci_files()
        # clean up all the old debris created by my past failures
        # clean_up()

    # check to make sure the final 2D science file doesnt already exist. If it does then
    # dont bother
    if len(glob.glob('J*.fits')) < 1:

        flats = glob.glob('MCgcal*.fits')
        if len(flats) == 0:
            # make the normalized flat as well as the stacked one used to find slit edges
            make_flat(summary, cent_waves)

        # reduce the science and arc images and apply wavelength transformation
        # this produces non sky subtracted, non-combined, transformed 2D slits to then
        # be loaded into PYPIT
        reduce_science(summary, cent_waves)
    if clean:
        # clean up
        delete_tmp_files()
        clean_up()

    else:
        print('This folder is already reduced. skipping')

    # move files to their resting place even if they already exist just to make sure
    move_2Dspecred(summary)


if __name__ == "__main__":
    # execute only if run as a script
    make_2Dspec()
