import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table

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
