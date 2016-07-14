# -*- coding: utf8 -*-

"""
    Goodman CCDRED - Performs ccd reductions of Goodman spectroscopic data..

    This script performs CCD reduction for spectra taken with the Goodman
    High Throughput Spectrograph at SOAR Telescope. The scrip will make
    (in order):

    - BIAS subtraction
    - TRIMMING including slit edges identification
    - FLAT correction
    - COSMIC rays rejection (optional)

    Users can add a flag in order to clean up science data from cosmic rays, which
    are removed by using the LACosmic code (P. G. van Dokkum, 2001, PASP, 113, 1420)

    This script was designed to make CCD reduction for any spectral configuration, but
    the input dir must contains only an unique spectral configuration (binning, grating,
    slit, gain, rdnoise, CCD ROI, etc). The input dir should contain only the following
    frames:

    1 - BIAS frames
    2 - FLAT frames (Flats taken between science exposures during the night
                     will be combined as normal afternoon calibrations)
    3. ARC frames   (data from focus sequence will not be reduced)

    4. SCIENCE and/or STANDARD frames

    # ToDo

    - Consider internal illumination correction

    Documentations for specific functions of the code can be found directly
    in the corresponding function.

    David Sanmartim (dsanmartim at ctio.noao.edu, dsanmartim at gemini.edu)
    July 2016

    Thanks to Bruno Quint for all comments and helping.

"""

import os
import glob
import argparse
import numpy as np
from astropy.io import fits
from astropy import log
from astropy import units as u
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import ccdproc
from ccdproc import ImageFileCollection
from ccdproc import CCDData


class Main:
    def __init__(self):

        self.raw_path = str(os.path.join(args.raw_path[0], ''))
        self.red_path = str(os.path.join(args.red_path[0], ''))
        
        self.clean = args.clean

        # checking the reduction directory
        if not os.path.isdir(self.red_path):
            os.mkdir(self.red_path)
        os.chdir(self.red_path)

        log.propagate = False
        log.disable_warnings_logging()

    @staticmethod
    def fit_spline3(y, order=3, nsum=5):
        y_resampled = [np.median(y[i:i + nsum]) for i in range(0, len(y) - len(y) % nsum, nsum)]
        x_resampled = np.linspace(0, len(y), len(y_resampled))

        # Fitting
        f = interp1d(x_resampled, y_resampled, kind=order, bounds_error=True)

        # Return function to be constructed with any other x array
        return f

    # Local Minima and Maxima
    @staticmethod
    def local_minmax(data, nmin=1, nmax=1):

        # Identifying indices of local minima-maxima points
        id_min = (np.gradient(np.sign(np.gradient(data))) > 0).nonzero()[0]  # index of local min
        id_max = (np.gradient(np.sign(np.gradient(data))) < 0).nonzero()[0]  # index of local max

        # Taking values at min/max points
        list_min, list_max = data[id_min], data[id_max]

        # Sorting minima-maxima values (bigger --> lower)
        list_min, id_min = (list(p) for p in zip(*sorted(zip(list_min, id_min), reverse=False)))
        list_max, id_max = (list(p) for p in zip(*sorted(zip(list_max, id_max), reverse=True)))

        # Taking the desired number of local minima-maxima points
        list_min, list_max, id_min, id_max = list_min[0:nmin], list_max[0:nmax], id_min[0:nmin], id_max[0:nmax]

        return list_min, list_max, id_min, id_max

    @staticmethod
    def clean_path(path):
        """
        Clean up directoy.
        """
        if os.path.exists(path):
            for _file in glob.glob(os.path.join(path, '*')):
                os.remove(_file)

    @staticmethod
    def fix_header_and_shape(input_path, output_path, prefix='', overwrite=False):
        """
        Remove/Update some  inconvenient parameters in the header of the Goodman FITS
        files. Some of these parameters contain non-printable ASCII characters. The ouptut
        files are created in the output_path. Also convert fits from 3D [1,X,Y] to 2D [X,Y].

        Parameters
        ----------
        indir: str
        outdir: str
        prefix: str, optional
            path to directory containing FITS files
        overwrite: boolean
        """

        for _file in sorted(glob.glob(os.path.join(input_path, '*.fits'))):

            ccddata, hdr = fits.getdata(_file, header=True, ignore_missing_end=True)

            # 3D to 2D
            ccddata = ccddata[0]

            # Inspectin True or False
            _param0 = 'PARAM0' in hdr
            _param61 = 'PARAM61' in hdr
            _param62 = 'PARAM62' in hdr
            _param63 = 'PARAM63' in hdr
            # _pcomm = 'COMMENT' in hdr
            _phist = 'HISTORY' in hdr
            _pinst = 'INSTRUME' in hdr
            _pnaxis = 'NAXIS' in hdr
            _pnaxis3 = 'NAXIS3' in hdr

            # cleaning header
            try:

                if _pnaxis is True:
                    hdr['NAXIS'] = 2
                if _pnaxis3 is True:
                    hdr.remove(keyword='NAXIS3')
                # if _pcomm is True:
                #    hdr.remove(keyword='COMMENT')
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

                fits.writeto(os.path.join(output_path, '') + prefix + os.path.basename(_file), ccddata, hdr,
                             clobber=overwrite)
                log.info('Keywords header of ' + os.path.basename(_file) + ' have been updated')

            except:
                pass
        print '\n'
        log.info('Done --> All keywords header are updated.')
        return

    def find_slitedge(self, ccddata):

        # Reading and Collapsing flat in the dispersion direction
        flat_collapsed = np.sum(ccddata, axis=1) / ccddata.shape[1]
        lines = np.arange(0, np.size(flat_collapsed), 1)

        # Excluding first pixels in the spatial direction
        cut = 3
        c_flat = flat_collapsed[cut:-cut]
        c_lines = np.arange(0, c_flat.size, 1)

        # Fitting cubic spline. It's working very well with order=5, nsum=2
        func_splin3 = self.fit_spline3(c_flat, order=5, nsum=2)
        smooth_flat = func_splin3(c_lines)

        # Compute 1st and flat smoothed
        dy = np.gradient(smooth_flat)
        dy2 = np.gradient(dy)

        # Regions to compute local minina-maxima
        # Region one: it represent first 40 percent of all data
        # Region two: ... last 40%
        pixa, pixb = int(len(c_flat) * 0.4), int(len(c_flat) * 0.6)
        dy2_one, dy2_two = dy2[0:pixa], dy2[pixb:]

        # Reg. 1: Compute local min/max of the 2nd derivative
        list_min_1, list_max_1, id_min_1, id_max_1 = self.local_minmax(dy2_one, nmin=1, nmax=1)
        list_min_2, list_max_2, id_min_2, id_max_2 = self.local_minmax(dy2_two, nmin=1, nmax=1)

        # Indice have to be reshifted to the original indices of the function dy2
        id_min_2 = np.array(id_min_2) + pixb

        # Slit edges are the local maxima/minima 1/2 [accounting the cutted pixels]
        slit_1, slit_2 = int(np.array(id_min_1) + cut), int(np.array(id_min_2) + cut)

        plot = 'no'
        if plot is 'yes':
            c_lines += cut
            plt.plot(lines, flat_collapsed, 'k-', label='Flat Collapsed')
            plt.plot(lines[slit_1:slit_2], flat_collapsed[slit_1:slit_2], 'r-', label='Cutted Flat')
            plt.plot(c_lines, dy, 'g-', label="Dy/dx")
            plt.plot(c_lines, dy2, 'y-', label="Dy2/dx")
            plt.plot(slit_1, list_min_1, 'bo', label='Slit Edge 1 ')
            plt.plot(slit_2, list_min_2, 'ro', label='Slit Edge 2')
            plt.xlim(lines.min() - 50, lines.max() + 50)
            plt.legend(loc='best')
            plt.show()

        return slit_1, slit_2

    def create_master_flat(self, image_collection):

        global slit1, slit2, master_flat, master_flat_nogrt

        # Creating dict. of flats. The key values are expected to be: GRATIN_ID and '<NO GRATING>'
        # if there is flat taken w/o grating
        grt_list = image_collection.values('grating', unique=True)
        dic_all_flats = {}
        for grt in sorted(grt_list):
            dic_all_flats[str(grt)] = image_collection.files_filtered(obstype='FLAT', grating=str(grt))

        # Dict. for flats with grating and without gratintg
        dic_flat = {grt: dic_all_flats[grt] for grt in dic_all_flats if grt != "<NO GRATING>"}
        dic_flatnogrt = {grt: dic_all_flats[grt] for grt in dic_all_flats if grt == "<NO GRATING>"}

        if len(dic_flat.values()) != 0:

            for grt in dic_flat.keys():

                flat_list = []
                log.info('Combining and trimming flat frames:')
                for filename in dic_flat[grt]:
                    log.info(filename)
                    ccd = CCDData.read(os.path.join(image_collection.location, '') + filename, unit=u.adu)
                    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])
                flat_list.append(ccd)
                # combinning and trimming slit edges
                master_flat = ccdproc.combine(flat_list, method='median', sigma_clip=True,
                                              sigma_clip_low_thresh=1.0, sigma_clip_high_thresh=1.0)
                print '\n Wait a minute... finding slit edges.'
                slit1, slit2 = self.find_slitedge(master_flat)
                master_flat = ccdproc.trim_image(master_flat[slit1:slit2, :])
                master_flat_name = 'master_flat_' + grt[5:] + '.fits'
                master_flat.write(master_flat_name, clobber=True)

                print '\n'
                log.info('A master flat has been created.')

        else:
            print '\n'
            log.info('Flat files have not been found.')

        if len(dic_flatnogrt.values()) != 0:

            for grt in dic_flatnogrt.keys():
                flatnogrt_list = []
                log.info('Combining and trimming flat frame taken without grating:')
                for filename in dic_flatnogrt[grt]:
                    log.info(filename)
                    ccd = CCDData.read(os.path.join(image_collection.location, '') + filename, unit=u.adu)
                    ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])
                flatnogrt_list.append(ccd)
                # combining and trimming slit edges
                master_flat_nogrt = ccdproc.combine(flatnogrt_list, method='median', sigma_clip=True,
                                                    sigma_clip_low_thresh=1.0, sigma_clip_high_thresh=1.0)
                master_flat_nogrt = ccdproc.trim_image(master_flat_nogrt[slit1:slit2, :])
                master_flat_nogrt_name = 'master_flat_nogrt.fits'
                master_flat_nogrt.write(master_flat_nogrt_name, clobber=True)
                print '\n'
                log.info('A master flat (taken without grating) have been created.')
        else:
            print '\n'
            log.info('Flat files (taken without grating) have not been found.')

        return

    @staticmethod
    def create_master_bias(image_collection):
        global master_bias
        bias_list = []
        log.info('Combining and trimming bias frames:')
        for filename in image_collection.files_filtered(obstype='BIAS'):
            log.info(filename)
            ccd = CCDData.read(os.path.join(image_collection.location, '') + filename, unit=u.adu)
            # Finding overscan regions... getting from header and assuming it is at the right edge...
            # over_start = int((ccd.header['TRIMSEC'].split(':'))[1].split(',')[0]) - 1
            # over_start += 10 / int(ccd.header['CCDSUM'][0])
            # ccd = ccdproc.subtract_overscan(ccd, median=True, overscan_axis=1, overscan=ccd[:, over_start:])
            ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])
            bias_list.append(ccd)
        master_bias = ccdproc.combine(bias_list, method='median', sigma_clip=True, sigma_clip_low_thresh=1.0,
                                      sigma_clip_high_thresh=1.0)
        master_bias = ccdproc.trim_image(master_bias[slit1:slit2, :])
        master_bias.write('master_bias.fits', clobber=True)

        print '\n'
        log.info('A master bias have been created.')
        return

    # TODO create a function to remove bad pixels from master bias

    @staticmethod
    def reduce_arc(image_collection):
        # reduce the arc frames
        log.info('Reduding Arc frames:')
        for filename in image_collection.files_filtered(obstype='COMP'):
            log.info(filename)
            ccd = CCDData.read(os.path.join(image_collection.location, '') + filename, unit=u.adu)
            ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])
            ccd = ccdproc.trim_image(ccd[slit1:slit2, :])
            ccd = ccdproc.subtract_bias(ccd, master_bias)
            ccd = ccdproc.flat_correct(ccd, master_flat)
            ccd.write('arc_' + filename, clobber=True)
        print '\n'
        log.info('Arc frames have been reduced. ')

    @staticmethod
    def reduce_sci(image_collection, clean):
        # reduce the sci frames
        log.info('Reduding Sci (Std) frames:')
        for filename in image_collection.files_filtered(obstype='OBJECT'):
            log.info(filename)
            ccd = CCDData.read(os.path.join(image_collection.location, '') + filename, unit=u.adu)
            ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])
            ccd = ccdproc.trim_image(ccd[slit1:slit2, :])
            ccd = ccdproc.subtract_bias(ccd, master_bias)
            ccd = ccdproc.flat_correct(ccd, master_flat)
            if clean is True:
                log.info('Cleaning cosmic rays... ')
                print '\n'
                nccd = ccdproc.cosmicray_lacosmic(ccd, sigclip=3.0, sigfrac=2.0, objlim=2.0,
                                                  gain=float(ccd.header['GAIN']),
                                                  readnoise=float(ccd.header['RDNOISE']),
                                                  satlevel=np.inf, sepmed=True, fsmode='convolve',
                                                  psfmodel='gaussy', verbose=True)
                nccd.write('c_obj_' + filename, clobber=True)
            ccd.write('obj_' + filename, clobber=True)

        print '\n'
        log.info('Sci (std) frames have been reduced ')

    def run(self):

        # cleaning up the reduction dir
        self.clean_path(self.red_path)

        # Fixing header and shape of raw data
        self.fix_header_and_shape(self.raw_path, self.red_path, prefix='f', overwrite=True)

        # Create image file collection for raw data
        ic = ImageFileCollection(self.red_path)

        # Create master_flats
        self.create_master_flat(ic)

        # Create master bias
        self.create_master_bias(ic)

        # Reduce Arc frames
        self.reduce_arc(ic)

        # Reduce Sci frames
        self.reduce_sci(ic, self.clean)

        return


if __name__ == '__main__':

    # Parsing Arguments ---
    parser = argparse.ArgumentParser(description="Goodman - CCD Data Reduction.")

    parser.add_argument('-c', '--clean', action='store_true',
                        help="Clean cosmic rays from science data")

    parser.add_argument('raw_path', metavar='raw_path', type=str, nargs=1,
                        help="Full path to raw data (e.g. /home/jamesbond/soardata/).")

    parser.add_argument('red_path', metavar='red_path', type=str, nargs=1,
                        help="Full path to reduced data (e.g /home/jamesbond/soardata/RED/).")

    args = parser.parse_args()

    main = Main()
    main.run()

else:
    print('goodman_ccdreduction.py is not being executed as main.')
