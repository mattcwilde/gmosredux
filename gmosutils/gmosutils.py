import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from shutil import copyfile

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

def delete_inter_files():
    """Get rid of the intermediate gsreduced files iraf produces

    """
    intermediate_files = glob.glob('*gs*fits')
    for file in intermediate_files:
        try:
            os.remove(file)
        except OSError:
            pass


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


def get_bias():
    if not os.path.exists('MCbiasFull.fits'):
        print("No bias in this folder. copying from ../../../Bias/")
        try:
            copyfile('../../../Bias/MCbiasFull.fits', './')
        except OSError:
            print('No bias? maybe need to make one')
            pass
    else:
        print("bias already exists, ready to roll!")
