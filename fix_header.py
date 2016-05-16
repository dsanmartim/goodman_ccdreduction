__author__ = 'davidsanm'

import glob
import os

from astropy.io import fits
from astropy import log

def fix_header(indir, outdir, prefix='', overwrite=False):
    """
    Remove PARAM0, PARAM61, PARAM62, PARAM63 and other inconvenient
    parameters in the header of the Goodman output files. Some of
    these parameters contain non-printable ASCII characters.

    Parameters
    ----------
    indir: str
    outdir: str
    prefix: str, optional
        path to directory containing FITS files
    overwrite: boolean
    """

    for _file in sorted(glob.glob(os.path.join(indir, '*.fits'))):

        hdulist = fits.open(_file, ignore_missing_end=True)
        hdr = hdulist[0].header

        # cleaning header
        _param0 = 'PARAM0' in hdr
        _param61 = 'PARAM61' in hdr
        _param62 = 'PARAM62' in hdr
        _param63 = 'PARAM63' in hdr
        _pcomm = 'COMMENT' in hdr
        _phist = 'HISTORY' in hdr
        _pinst = 'INSTRUME' in hdr


        # noinspection PyBroadException
        try:

            if _pcomm is True:
                hdr.remove(keyword='COMMENT')
            if _phist is True:
                hdr.remove(keyword='HISTORY')
            if _pinst is True:
                hdr.remove(keyword='INSTRUME')
            if _param0 is True:
                hdr.remove(keyword='PARAM0')
            if _param61 is True:
                hdr.remove(keyword='PARAM61')
            if _param62 is True:
                hdr.remove(keyword='PARAM62')
            if _param63 is True:
                hdr.remove(keyword='PARAM63')

            hdulist.writeto(os.path.join(outdir, '') + prefix + os.path.basename(_file), clobber=overwrite)
            hdulist.close()
            log.info('Keywords header of ' + os.path.basename(_file) + ' have been updated')

        except:
            pass

    log.info('\n Done --> All keywords header are updated.')
    return
