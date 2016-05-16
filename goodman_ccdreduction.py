# __author__ = 'davidsanm'

import os
import glob
import numpy as np
from astropy.io import fits
from astropy import log
from astropy import units as u

from ccdproc import ImageFileCollection
from ccdproc import CCDData
import ccdproc


class Main():
    def __init__(self):

        # if len(sys.argv)!=3:
        #   print('Usage:\npython goodman_ccdreduction.py [full_path_to_raw_data] [full_path_to_reduced_data]\n')
        #   exit()

        # indir = sys.argv[1]
        # outdir = sys.argv[2]

        self.indir = '/home/davidsanm/Dados/Soar/SO2014A-028/2014-07-28'
        self.outdir = '/home/davidsanm/Dados/Soar/SO2014A-028/2014-07-28/RED'

        # checking the output dir
        if not os.path.isdir(self.outdir): os.mkdir(self.outdir)
        os.chdir(self.outdir)

    def clean_dir(self, dir):
        """
        Clean up the directory.

        Parameters
        ----------
        dir: str, input directory
        """
        import glob

        if os.path.exists(dir):
            for file in glob.glob(os.path.join(dir, '*')):
                os.remove(file)

    def fix_header(self, indir, outdir, prefix='', overwrite=False):
        '''
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
        '''

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

    def run(self):

        # cleaning up the output dir
        self.clean_dir(self.outdir)

        # Cleaning Header from Goodman files
        self.fix_header(self.indir, self.outdir, prefix='c', overwrite=True)

        # change this to point to your header updates raw data directory
        ic1 = ImageFileCollection(self.outdir)

        return


if __name__ == '__main__':

    main = Main()
    main.run()

else:
    print('goodman_ccdreduction.py is not being executed as main.')

'''
# create the bias frames
bias_list = []
for filename in ic1.files_filtered(obstype='BIAS'):
    print os.path.join(ic1.location, '') + filename
    ccd = CCDData.read(os.path.join(ic1.location, '') + filename, unit=u.adu)
    # this has to be fixed as the bias section does not include the whole section that will be trimmed
    # ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])
    bias_list.append(ccd)
master_bias = ccdproc.combine(bias_list, output = 'master_bias.fits', method='median')
#master_bias.write('master_bias_blue.fits', clobber=True)

red_bias_list = []
for filename in ic1.files_filtered(obstype='Bias', isiarm='Red arm'):
    print ic1.location + filename
    ccd = CCDData.read(ic1.location + filename, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    red_bias_list.append(ccd)
master_bias_red = ccdproc.combine(red_bias_list, method='median')
master_bias_red.write('master_bias_red.fits', clobber=True)

#create the flat fields
red_flat_list = []
for filename in ic1.files_filtered(obstype='Flat', isiarm='Red arm'):
    ccd = CCDData.read(ic1.location + filename, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_red)
    red_flat_list.append(ccd)
master_flat_red = ccdproc.combine(red_flat_list, method='median')
master_flat_red.write('master_flat_red.fits', clobber=True)

blue_flat_list = []
for filename in ic1.files_filtered(obstype='Flat', isiarm='Blue arm'):
    ccd = CCDData.read(ic1.location + filename, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_blue)
    blue_flat_list.append(ccd)
master_flat_blue = ccdproc.combine(blue_flat_list, method='median')
master_flat_blue.write('master_flat_blue.fits', clobber=True)

#reduce the arc frames
for filename in ic1.files_filtered(obstype='Arc', isiarm='Blue arm'):
    hdu = fits.open(ic1.location + filename)
    ccd = CCDData(hdu[1].data, header=hdu[0].header+hdu[1].header, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_blue)
    ccd = ccdproc.flat_correct(ccd, master_flat_blue)
    ccd.data = ccd.data.T
    ccd.write('arc_'+filename, clobber=True)

red_flat_list = []
for filename in ic1.files_filtered(obstype='Arc', isiarm='Red arm'):
    hdu = fits.open(ic1.location + filename)
    ccd = CCDData(hdu[1].data, header=hdu[0].header+hdu[1].header, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_red)
    ccd = ccdproc.flat_correct(ccd, master_flat_red)
    ccd.data = ccd.data.T
    ccd.write('arc_'+filename, clobber=True)

#reduce the sky frames
for filename in ic1.files_filtered(obstype='Sky', isiarm='Blue arm'):
    hdu = fits.open(ic1.location + filename)
    ccd = CCDData(hdu[1].data, header=hdu[0].header+hdu[1].header, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_blue)
    ccd = ccdproc.flat_correct(ccd, master_flat_blue)
    ccd.data = ccd.data.T
    ccd.write('sky_'+filename, clobber=True)

for filename in ic1.files_filtered(obstype='Sky', isiarm='Red arm'):
    hdu = fits.open(ic1.location + filename)
    ccd = CCDData(hdu[1].data, header=hdu[0].header+hdu[1].header, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_red)
    ccd = ccdproc.flat_correct(ccd, master_flat_red)
    ccd.data = ccd.data.T
    ccd.write('sky_'+filename, clobber=True)

#reduce the object frames
for filename in ic1.files_filtered(obstype='TARGET', isiarm='Blue arm'):
    hdu = fits.open(ic1.location + filename)
    ccd = CCDData(hdu[1].data, header=hdu[0].header+hdu[1].header, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_blue)
    ccd = ccdproc.flat_correct(ccd, master_flat_blue)
    ccd.data = ccd.data.T
    ccd.write('obj_'+filename, clobber=True)

for filename in ic1.files_filtered(obstype='Target', isiarm='Red arm'):
    hdu = fits.open(ic1.location + filename)
    ccd = CCDData(hdu[1].data, header=hdu[0].header+hdu[1].header, unit = u.adu)
    #this has to be fixed as the bias section does not include the whole section that will be trimmed
    ccd = ccdproc.subtract_overscan(ccd, median=True,  overscan_axis=0, fits_section='[1:966,4105:4190]')
    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'] )
    ccd = ccdproc.subtract_bias(ccd, master_bias_red)
    ccd = ccdproc.flat_correct(ccd, master_flat_red)
    ccd.data = ccd.data.T
    ccd.write('obj_'+filename, clobber=True)

'''
