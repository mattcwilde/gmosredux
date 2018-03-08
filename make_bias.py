"""
This script selects out the bias if its mixed into other files and generates a MasterCal Bias image
"""
from __future__ import print_function, division
import os
from pyraf import iraf
from pyraf.iraf import gemini, gmos, gemtools
import glob
from gmosutils import gmosutils


rawPath = './raw/'
dataPath = './'

if not os.path.exists("{0}/MCbiasFull.fits".format(dataPath)):
    print(" --Creating Bias MasterCal--")

    # Use primarily the default task parameters.
    gemtools.gemextn.unlearn()  # Disarm a bug in gbias
    gmos.gbias.unlearn()
    biasFlags = {'logfile': 'biasLog.txt', 'rawpath': rawPath,
                 'fl_vardq': 'yes', 'verbose': 'no'}

    # find the raw fits files
    raw_files = glob.glob(rawPath+'/*fits')

    # make an astropy table summarising them
    summary = gmosutils.observation_summary(raw_files)

    # identify only the bias files and get their filename
    biasFileNames = summary[(summary['OBSTYPE'] == 'BIAS')]['FILENAME']
    # convert the filename to a list
    biasFull = [os.path.basename(f) for f in biasFileNames]

    if len(biasFull) > 20:
        # write the filneame list to a file file
        # IRAF cant take more that ~10 comma seperated files as input
        with open('biases.lis', 'w') as f:
            [f.write(x + '\n') for x in biasFull]
        gmos.gbias('@biases.lis', 'MCbiasFull.fits', **biasFlags)

        # Clean up intermediate files
        iraf.imdel('g*.fits')
    else:
        print('Not enough BIAS images. Make sure you are in the correct folder')

    print(" -- DONE: made MCbiasFull.fits")
else:
    print('Bias already exists, dummy.')