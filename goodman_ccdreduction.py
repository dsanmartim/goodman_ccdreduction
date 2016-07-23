# -*- coding: utf8 -*-

"""
    ## PyGoodman CCD Reduction

    Goodman CCDRed performs ccd reductions for Goodman spectroscopic data.

    This script performs CCD reduction for spectra taken with the Goodman
    High Throughput Spectrograph at SOAR Telescope. The script will make
    (in order):

    - BIAS subtraction
    - TRIMMING including slit edges identification
    - FLAT correction
    - COSMIC rays rejection (optional)

    Users can add a flag in order to clean up science data from cosmic rays, which
    are removed by using the LACosmic code (P. G. van Dokkum, 2001, PASP, 113, 1420)

    ### Data Structure

    Documentations for specific functions of the code can be found directly
    in the corresponding function. (To be done...)

    This script was designed to make CCD reduction for any spectral configuration, but
    the input dir must contains only an unique spectral configuration (binning, grating,
    slit, gain, rdnoise, CCD ROI, etc). The input dir should contain only the following
    frames:

    - BIAS frames
    - FLAT frames  (Flats taken between science exposures will be trimmed and bias subtracted.)
    - ARC frames   (data from focus sequence will not be reduced)
    - SCIENCE and/or STANDARD frames

    # ToDo

    - Consider internal illumination correction

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

from astropy.coordinates import EarthLocation
from astropy.time import Time, TimeDelta
from astroplan import Observer

import ccdproc
from ccdproc import ImageFileCollection
from ccdproc import CCDData
import warnings

__author__ = 'David Sanmartim'
__date__ = '2016-07-15'
__version__ = "0.1"
__email__ = "dsanmartim@ctio.noao.edu"


class Main:
    def __init__(self):

        global memlim, Geodetic_Location

        # Soar Geodetic Location
        Geodetic_Location = ['-70d44m01.11s', '-30d14m16.41s', 2748]

        # Memory Limit to be used
        memlim = 16E9

        # Paths
        self.raw_path = str(os.path.join(args.raw_path[0], ''))
        self.red_path = str(os.path.join(args.red_path[0], ''))

        self.clean = args.clean
        self.slit = args.slit

        # checking the reduction directory
        if not os.path.isdir(self.red_path):
            os.mkdir(self.red_path)
        os.chdir(self.red_path)

        warnings.filterwarnings('ignore')

        log.propagate = False

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
    def fix_header_and_shape(input_path, output_path, prefix, overwrite=False):
        """
        Remove/Update some  inconvenient parameters in the header of the Goodman FITS
        files. Some of these parameters contain non-printable ASCII characters. The ouptut
        files are created in the output_path. Also convert fits from 3D [1,X,Y] to 2D [X,Y].
        """
        for _file in sorted(glob.glob(os.path.join(input_path, '*.fits'))):

            ccddata, hdr = fits.getdata(_file, header=True, ignore_missing_end=True)
            # 3D to 2D
            ccddata = ccddata[0]
            # keywords to remove
            key_list_to_remove = ['PARAM0', 'PARAM61', 'PARAM62', 'PARAM63', 'NAXIS3', 'INSTRUME']

            # Keyword to be changed (3 --> 2)
            hdr['N_PARAM'] -= len(key_list_to_remove)
            hdr['NAXIS'] = 2

            # Specific keywords to be removed
            for key in key_list_to_remove:
                try:
                    if (key in hdr) is True:
                        hdr.remove(keyword=key)
                except KeyError:
                    pass

            # Removing duplicated keywords
            key_list = []
            for key in hdr.iterkeys():
                if key in key_list:
                    hdr.remove(keyword=key)
                key_list.append(key)

            hdr.add_history('Header and Shape fixed.')
            fits.writeto(os.path.join(output_path, '') + prefix + os.path.basename(_file), ccddata, hdr,
                         clobber=overwrite)
            log.info('Keywords header of ' + os.path.basename(_file) + ' have been updated --> ' + prefix
                     + os.path.basename(_file))
        log.info('Done: all keywords header have been updated.')
        print('\n')
        return

    def find_slitedge(self, ccddata):

        # Reading and Collapsing flat in the dispersion direction
        flat_collapsed = np.sum(ccddata, axis=1) / ccddata.shape[1]

        # Excluding first pixels in the spatial direction
        cut = 3
        c_flat = flat_collapsed[cut:-cut]
        c_lines = np.arange(0, c_flat.size, 1)

        # Fitting cubic spline. It's working very well with order=5, nsum=2
        func_splin3 = self.fit_spline3(c_flat, order=5, nsum=2)
        smooth_flat = func_splin3(c_lines)

        # Compute 1st and flat smoothed
        dy2 = np.gradient(np.gradient(smooth_flat))

        # Region one: it represent first 40 percent of all data. Region two: ... last 40%
        pixa, pixb = int(len(c_flat) * 0.4), int(len(c_flat) * 0.6)
        dy2_one, dy2_two = dy2[0:pixa], dy2[pixb:]

        # Reg. 1: Compute local min/max of the 2nd derivative
        list_min_1, _, id_min_1, _ = self.local_minmax(dy2_one, nmin=1, nmax=1)
        list_min_2, _, id_min_2, _ = self.local_minmax(dy2_two, nmin=1, nmax=1)

        # Slit edges are the local maxima/minima 1/2 [accounting the cutted pixels]
        slit_1, slit_2 = int(np.array(id_min_1) + cut), int(np.array(np.array(id_min_2) + pixb) + cut)

        return slit_1, slit_2

    @staticmethod
    def get_twilight_time(image_collection, observatory, long, lat, elevation, timezone, description):

        """
        image_collection: ccdproc object
        observatory: str, observatory name (e.g. 'Soar Telescope',
        long: str, dms or deg
        lat: str, dms or deg
        elevation: int, meters (define through ellipsoid WGS84)
        timezone: str, eg. 'UTC'
        description: str, short description of the observatory

        return: str, twilight evening and twilinght morning (format 'YYYY-MM-DDT00:00:00.00')
        """
        soar_loc = EarthLocation.from_geodetic(long, lat, elevation * u.m, ellipsoid='WGS84')

        soar = Observer(name=observatory, location=soar_loc, timezone=timezone,
                        description=description)

        dateobs_list = image_collection.values('date-obs')
        time_first_frame, time_last_frame = Time(min(dateobs_list)).isot, Time(max(dateobs_list)).isot

        twilight_evening = soar.twilight_evening_astronomical(Time(time_first_frame), which='nearest').isot
        twilight_morning = soar.twilight_morning_astronomical(Time(time_last_frame), which='nearest').isot

        return twilight_evening, twilight_morning

    def get_daycal_flat(self, image_collection):
        """
        image_collection: ccdproc object
        return: list of flats
        """
        twi_eve, _ = self.get_twilight_time(image_collection, 'Soar Telescope', long=Geodetic_Location[0],
                                            lat=Geodetic_Location[1], elevation=Geodetic_Location[2],
                                            timezone='UTC', description='Soar Telescope on Cerro Pachon, Chile')

        df = image_collection.summary.to_pandas()
        dfobj = df['file'][(df['date-obs'] < twi_eve) & (df['obstype'] == 'FLAT')]
        dayflat_list = dfobj.values.tolist()

        return dayflat_list

    def create_daymaster_flat(self, image_collection, twilight_evening, slit):

        global slit1, slit2, \
            master_flat, master_flat_nogrt, \
            master_flat_name, master_flat_nogrt_name

        master_flat = []
        master_flat_nogrt = []

        # Creating dict. of flats. The key values are expected to be: GRATIN_ID and '<NO GRATING>'
        # if there is flat taken w/o grating
        grt_list = image_collection.values('grating', unique=True)
        dic_all_flats = {}

        df = image_collection.summary.to_pandas()
        for grt in sorted(grt_list):
            dfobj = df['file'][(df['obstype'] == 'FLAT') & (df['date-obs'] < twilight_evening) & (df['grating'] == grt)]
            dic_all_flats[str(grt)] = dfobj.tolist()

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
                master_flat = ccdproc.combine(flat_list, method='median', mem_limit=memlim, sigma_clip=True,
                                              sigma_clip_low_thresh=1.0, sigma_clip_high_thresh=1.0)
                if slit is True:
                    print('\n Finding slit edges... \n')
                    slit1, slit2 = self.find_slitedge(master_flat)
                    master_flat = ccdproc.trim_image(master_flat[slit1:slit2, :])

                master_flat_name = 'master_flat_' + grt[5:] + '.fits'
                master_flat.write(master_flat_name, clobber=True)

                log.info('Done: master flat has been created --> ' + 'master_flat_' + grt[5:] + '.fits')
                print('\n')

        else:
            log.info('Flat files have not been found.')
            print('\n')

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
                master_flat_nogrt = ccdproc.combine(flatnogrt_list, method='median', mem_limit=memlim, sigma_clip=True,
                                                    sigma_clip_low_thresh=1.0, sigma_clip_high_thresh=1.0)
                if slit is True:
                    master_flat_nogrt = ccdproc.trim_image(master_flat_nogrt[slit1:slit2, :])

                master_flat_nogrt_name = 'master_flat_nogrt.fits'
                master_flat_nogrt.write(master_flat_nogrt_name, clobber=True)

                log.info('Done: master flat (taken without grating) have been created --> master_flat_nogrt.fits')
                print('\n')
        else:
            log.info('Flat files (taken without grating) have not been found.')
            print('\n')

        return

    @staticmethod
    def create_master_bias(image_collection, slit):

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
        master_bias = ccdproc.combine(bias_list, method='median', mem_limit=memlim, sigma_clip=True,
                                      sigma_clip_low_thresh=1.0, sigma_clip_high_thresh=1.0)
        if slit is True:
            master_bias = ccdproc.trim_image(master_bias[slit1:slit2, :])
        else:
            master_bias = master_bias
        master_bias.header['HISTORY'] = "Trimmed."
        master_bias.write('master_bias.fits', clobber=True)

        # Now I obtained bias... subtracting bias from master flat
        # Testing if master_flats are not empty arrays
        if (not master_flat) is False:
            fccd = ccdproc.subtract_bias(master_flat, master_bias)
            fccd.header['HISTORY'] = "Trimmed. Bias subtracted. Flat corrected."
            fccd.write(master_flat_name, clobber=True)

        if (not master_flat_nogrt) is False:
            ngccd = ccdproc.subtract_bias(master_flat_nogrt, master_bias)
            ngccd.header['HISTORY'] = "Trimmed. Bias subtracted. Flat corrected."
            ngccd.write(master_flat_nogrt_name, clobber=True)

        log.info('Done: a master bias have been created --> master_bias.fits')
        print('\n')
        return

    @staticmethod
    def reduce_nightflats(image_collection, twilight_evening, slit, prefix):

        # 40 min = time before twilight evening to be considered as begining of the night
        time_before = (Time(twilight_evening)-TimeDelta(2400.0, format='sec')).isot

        log.info('Reducing flat frames taken during the night...')
        df = image_collection.summary.to_pandas()
        dfobj = df['file'][(df['obstype'] == 'FLAT') & (df['date-obs'] > time_before) & (df['grating'] != '<NO GRATING>')]
        nightflat_list = dfobj.tolist()

        if len(nightflat_list) > 0:
            for filename in sorted(nightflat_list):
                log.info('Trimming and bias subtracting frame ' + filename + ' --> ' + prefix + filename)
                ccd = CCDData.read(os.path.join(image_collection.location, '') + filename, unit=u.adu)
                ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])
                if slit is True:
                    ccd = ccdproc.trim_image(ccd[slit1:slit2, :])
                ccd = ccdproc.subtract_bias(ccd, master_bias)
                ccd.header['HISTORY'] = "Trimmed. Bias subtracted."
                ccd.write(prefix + filename, clobber=True)
            log.info('Done --> Night flat frames have been reduced.')
            print('\n')
        return

    @staticmethod
    def reduce_arc(image_collection, slit, prefix):

        log.info('Reduding Arc frames...')
        arc_list = image_collection.files_filtered(obstype='COMP')
        if len(arc_list) > 0:
            for filename in sorted(arc_list):
                log.info('Reduding Arc frame ' + filename + ' --> ' + prefix + filename)
                ccd = CCDData.read(os.path.join(image_collection.location, '') + filename, unit=u.adu)
                ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])
                if slit is True:
                    ccd = ccdproc.trim_image(ccd[slit1:slit2, :])
                ccd = ccdproc.subtract_bias(ccd, master_bias)
                ccd = ccdproc.flat_correct(ccd, master_flat)
                ccd.header['HISTORY'] = "Trimmed. Bias subtracted. Flat corrected."
                ccd.write(prefix + filename, clobber=True)
            log.info('Done --> Arc frames have been reduced.')
            print('\n')
        return

    @staticmethod
    def reduce_sci(image_collection, slit, clean, prefix):

        log.info('Reduding Sci/Std frames...')
        for filename in image_collection.files_filtered(obstype='OBJECT'):
            log.info('Reduding Sci/Std frame ' + filename + ' --> ' + prefix + filename)
            ccd = CCDData.read(os.path.join(image_collection.location, '') + filename, unit=u.adu)
            ccd = ccdproc.trim_image(ccd, fits_section=ccd.header['TRIMSEC'])
            if slit is True:
                ccd = ccdproc.trim_image(ccd[slit1:slit2, :])
            ccd = ccdproc.subtract_bias(ccd, master_bias)
            ccd = ccdproc.flat_correct(ccd, master_flat)
            # OBS: cosmic ray rejection is working pretty well by defining gain = 1. It's not working
            # when we use the real gain of the image. In this case the sky level changes by a factor
            # equal the gain.
            if clean is True:
                log.info('Cleaning cosmic rays... ')
                nccd, _ = ccdproc.cosmicray_lacosmic(ccd.data, sigclip=2.5, sigfrac=2.0, objlim=2.0,
                                                     gain=float(ccd.header['GAIN']),
                                                     readnoise=float(ccd.header['RDNOISE']),
                                                     satlevel=np.inf, sepmed=True, fsmode='median',
                                                     psfmodel='gaussy', verbose=True)
                log.info('Cosmic rays have been cleaned ' + prefix + filename + ' --> ' + 'c' + prefix + filename)
                print('\n')
                nccd = np.array(nccd, dtype=np.double) / float(ccd.header['GAIN'])
                ccd.header['HISTORY'] = "Trimmed. Bias subtracted. Flat corrected."
                ccd.header['HISTORY'] = "Cosmic rays rejected."
                fits.writeto('c' + prefix + filename, nccd, ccd.header, clobber=True)
            elif clean is False:
                ccd.header['HISTORY'] = "Trimmed, Bias subtracted, Flat corrected."
            ccd.write(prefix + filename, clobber=True)
        log.info('Done: Sci/Std frames have been reduced.')
        print('\n')
        return

    def run(self):

        # cleaning up the reduction dir
        self.clean_path(self.red_path)

        # Fixing header and shape of raw data
        self.fix_header_and_shape(self.raw_path, self.red_path, prefix='h.', overwrite=True)

        # Create image file collection for raw data
        ic = ImageFileCollection(self.red_path)

        # Getting twilight time
        twi_evening, twi_morning  = self.get_twilight_time(ic, observatory='Soar Telescope', long=Geodetic_Location[0],
                                                           lat=Geodetic_Location[1], elevation=Geodetic_Location[2],
                                                           timezone='UTC',
                                                           description='Soar Telescope on Cerro Pachon, Chile')

        # Create master_flats
        self.create_daymaster_flat(ic, twi_evening, self.slit)

        # Create master bias
        self.create_master_bias(ic, self.slit)

        # Reduce Night Flat frames (if they exist)
        self.reduce_nightflats(ic, twi_evening, self.slit, prefix='z')

        # Reduce Arc frames
        self.reduce_arc(ic, self.slit, prefix='fz')

        # Reduce Sci frames
        self.reduce_sci(ic, self.slit, self.clean, prefix='fz')

        return


if __name__ == '__main__':

    # Parsing Arguments ---
    parser = argparse.ArgumentParser(description="Goodman - CCD Data Reduction.")

    parser.add_argument('-c', '--clean', action='store_true',
                        help="Clean cosmic rays from science data.")

    parser.add_argument('-s', '--slit', action='store_true',
                        help="Find slit edge to make an additional trimming (recommended).")

    parser.add_argument('raw_path', metavar='raw_path', type=str, nargs=1,
                        help="Full path to raw data (e.g. /home/jamesbond/soardata/).")

    parser.add_argument('red_path', metavar='red_path', type=str, nargs=1,
                        help="Full path to reduced data (e.g /home/jamesbond/soardata/RED/).")

    args = parser.parse_args()

    main = Main()
    main.run()

else:
    print('goodman_ccdreduction.py is not being executed as main.')
