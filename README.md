# PyGoodman CCD Reduction

Warning: work is still in progress

PyGoodman CCDRed performs ccd reductions for Goodman spectroscopic data.

This script performs CCD reduction for long-slit spectra taken with the 
Goodman High Throughput Spectrograph at SOAR Telescope. The script will 
carry out (subsequently):

 - BIAS subtraction
 - TRIMMING including slit edges identification (it does not work
   for MOS spectroscopy)
 - FLAT correction
 - COSMIC rays rejection (optional)

Users can run the code adding a flag in order to clean up science data 
from cosmic rays, which are removed by using the modified version of 
the LACosmic code (P. G. van Dokkum, 2001, PASP, 113, 1420) available 
through the astropy affiliated package ccdproc (and astrocrappy).

For a complete repository with tools for spectroscopy, please see 
[Simon Torres](https://github.com/simontorres/goodman)'s github 
repository.


### I/O Data Structure

Although this script was designed to do CCD reduction for any 
spectrograph configuration, the input directory should contain only 
a unique spectral configuration (binning, grating, slit, gain, rdnoise, 
CCD ROI, etc). The input dir should contain only the following frames:

- BIAS frames
- FLAT frames (Flats taken between science exposures will be trimmed
  and bias subtracted.)
- ARC frames (data from focus sequence will not be reduced)
- SCIENCE and/or STANDARD frames

Please, inspect your calibration and throw out the bad ones. The output 
data will have the same filename as the input, but with a prefix "fzh". 
It means data have their header updated (h), bias subtracted (z) and 
flat corrected (f). The prefix "c_fzh" means that cosmic ray correction 
was also applied.

### How to use it...

It can be be executed in terminal by running: 

    $ python goodman_ccdreduction.py [options] raw_path red_path 
    
More information about the options and how to use it can be obtained by 
typing:

    $ python goodman_ccdreduction.py --help

or

    $ python goodman_ccdreduction.py -h

### ToDo (Short List)

- Consider internal illumination correction (in fact this will not be 
done anymore)
- Disable the flat correction if there isn't a master no-grating flat
- Automatically determine the best input parameters for LACOSMIC

### Suggestions and Questions

If you have any doubt or question, please contact David Sanmartim or Simon 
Torres, who is the current maintainer:

<b>davidsanm at gmail.com</b>
<b>dsanmartim at gemini.edu</b>
<b>storres at ctio.noao.edu</b>
   
Oct 2016

### Notes: By Simon

Install astroplan currently on **heavy development** might requiere a lot
of maintenance

By installing astroplan, it downgrades numpy
