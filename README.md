# goodman_ccdreduction

    # Goodman CCDRED - Performs ccd reductions for Goodman spectroscopic data..

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
