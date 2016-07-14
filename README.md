## PyGoodman CCD Reduction

Warning: work is still in progress

Goodman CCDRed performs ccd reductions for Goodman spectroscopic data..

This script performs CCD reduction for spectra taken with the Goodman High Throughput Spectrograph at SOAR Telescope. The scrip will make (in order):

 - BIAS subtraction
 - TRIMMING including slit edges identification
 - FLAT correction
 - COSMIC rays rejection (optional)

Users can add a flag in order to clean up science data from cosmic rays, which are removed by using the LACosmic code (P. G. van Dokkum, 2001, PASP, 113, 1420).

### Data Structure

This script was designed to make CCD reduction for any spectral configuration, but the input dir must contains only an unique spectral configuration (binning, grating, slit, gain, rdnoise, CCD ROI, etc). The input dir should contain only the following frames:

 - BIAS frames
 - FLAT frames (Flats taken between science exposures will be combined together with afternoon calibrations)
 - ARC frames (data from focus sequence will not be reduced)
 - SCIENCE and/or STANDARD frames
 
Documentation for specific functions of the code can be found directly in the corresponding function.

### How to use it...

It can be be executed in terminal running 

    $ python goodman_ccdreduction.py [options] raw_path red_path 
    
More information abotu the options and how to use it can be otained by using...

    $ python goodman_ccdreduction.py --help
or
    $ python goodman_ccdreduction.py --h

### ToDo (Short List)

 - Consider internal illumination correction (using the flats taken without grating
 -
 -
 -

### Suggestions and Questions

If you have any doubt or question, please contact

David Sanmartim 
<p><b>dsanmartim at ctio.noao.edud</b> (before Set 1st, 2016 ) </p> 
<p><b>dsanmartim at gemini.edu</b> (after Set 1st, 2016) </p> 
   
July 2016

