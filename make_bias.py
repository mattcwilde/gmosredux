"""
This script selects out the bias if its mixed into other files and generates a MasterCal Bias image
"""
from __future__ import print, division
import os
from pyraf import iraf
from pyraf.iraf import gemini, gmos
from gmosutils import gmosutils as gu


rawPath = 'data/raw'
dataPath = './'

if not os.path.exists("{0}/MCbiasFull.fits".format(dataPath)):
    print(" --Creating Bias MasterCal--")

    # Use primarily the default task parameters.
    gemtools.gemextn.unlearn()  # Disarm a bug in gbias
    gmos.gbias.unlearn()
    biasFlags = {'logfile': 'biasLog.txt', 'rawpath': rawPath,
                 'fl_vardq': 'yes', 'verbose': 'no'}


    raw_files = glob.glob(rawPath+'/*fits')
    summary = gu.observation_summary(raw_files)
    summary[summary['OBJECT'] != 'Bias']

    biasFileNameist = summary[(summary['OBSTYPE'] == 'BIAS')]['FILENAME']
    biasFull = [os.path.basename(f) for f in bias]

    if len(biasFull) > 20:
        # write a list file
        # IRAF cant take more that ~10 files as input
        with open('biases.lis', 'w') as f:
            [f.write(x + '\n') for x in biasFull]
        gmos.gbias('@biases.lis', 'MCbiasFull.fits', **biasFlags)

        # Clean up intermediate files
        iraf.imdel('g*.fits')
    else:
        print('Not enough BIAS images. Make sure you are in the correct folder')

else:
    print('')