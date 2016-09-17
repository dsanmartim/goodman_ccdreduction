# PyGoodman CCD Reduction

Warning: work is still in progress

Goodman CCDRed performs ccd reductions for Goodman spectroscopic
data.

This script performs CCD reduction for long-slit spectra taken
with the Goodman High Throughput Spectrograph at SOAR Telescope.
The script will make (in order):

 - BIAS subtraction
 - TRIMMING including slit edges identification (it does not work
   for MOS spectroscopy)
 - FLAT correction
 - COSMIC rays rejection (optional)

Users can add a flag in order to clean up science data from cosmic
rays, which are removed by using the modified version of the LACosmic
code (P. G. van Dokkum, 2001, PASP, 113, 1420) available through the
astropy affiliated package ccdproc (and astrocrappy).

### I/O Data Structure

This script was designed to make CCD reduction for any spectrograph
configuration, but the input directory should contain only an unique
spectral configuration (binning, grating, slit, gain, rdnoise, CCD ROI,
etc). The input dir should contain only the following frames:

- BIAS frames
- FLAT frames  (Flats taken between science exposures will be trimmed
  and bias subtracted.)
- ARC frames   (data from focus sequence will not be reduced)
- SCIENCE and/or STANDARD frames

Please, inspect you calibration and throw it out the bad ones. The output
data has the same filename of the input data, but with a prefix "fzh". It
means data has its header updated (h), bias subtracted (z) and flat corrected
(f). The prefix "c_fzh" means that cosmic ray correction was applied.

### How to use it...

It can be be executed in terminal running 

    $ python goodman_ccdreduction.py [options] raw_path red_path 
    
More information abotu the options and how to use it can be otained by 
using...

    $ python goodman_ccdreduction.py --help

or

    $ python goodman_ccdreduction.py -h

### ToDo (Short List)

- Consider internal illumination correction (in fact this will not be done)
- Disable the flat correction if there isn't a master no_grating flat
- Automatically determine the best input parameters for LACOSMIC

### Suggestions and Questions

If you have any doubt or question, please contact David Sanmartim 

<b>davidsanm at gmail.com</b>
<br>
<b>dsanmartim at gemini.edu</b>
   
July 2016

