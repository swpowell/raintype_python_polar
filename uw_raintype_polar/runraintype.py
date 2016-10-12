"""
 ****Rain-type Classification code of Powell et al. (2016, JTECH): Driver code in python*****#
#Author: Scott Powell and Stacy Brodzik, University of Washington
#Date: March 4, 2016
#Description: This is the driver code for the updated version of Steiner et al. (1995)
#   convective/stratiform classification code for use with polar coordinate datasets. Adds new
#   categories for echoes of mixed rain-type near convective cores and correctly identifies
#   isolated, often shallow, convection as convective instead of stratiform. For details, see
#   Powell, S.W., R.A. Houze, JR., and S.R. Brodzik, 2016: Rainfall-type categorization of radar
#   echoes using polar coordinate reflectivity data, J. Atmos. Oceanic Technol., 33, 523-538.
#   The variables listed in the left column immediately below are those in the user-input
#   parameters farther below. The variables listed in the right column below are the
#   corresponding variable names in Table 1 of Powell et al. (2016).
#
""" 
from __future__ import division   #For python2 users only. Alternatively, run interpreter with -Q flag.
import netCDF4 as nc4
import numpy as np
import os
import algorithm as alg
import cfrad_io as io
import warnings

"""
The variables listed in the left column immediately below are those in the user-input
  parameters farther below. The variables listed in the right column below are the
  corresponding variable names in Table 1 of Powell et al. (2016).

  Variable name in this code          Variable name in Powell et al.
  --------------------------          ------------------------------
  minZdiff                            a
  deepcoszero                         b
  shallowconvmin                      Z_shallow
  truncZconvthres                     Z_th
  dBZformaxconvradius                 Z_conv
  weakechothres                       Z_weak
  backgrndradius                      R_bg
  maxconvRadius                       R_conv
  minsize                             A_low
  startslope                          A_med
  maxsize                             A_high

  Inputs:
  refl = Reflectivity
  refl_missing_val = missing value in reflectivity data
  refl_dx (km) = horizontal spacing of grid
  minZdiff = factor for comparing echo to background reflectivity; see equation (1) in journal 
     article referenced above
  deepcoszero = see equation (1) in journal article referenced above
  shallowconvmin = minimum dBZ for classification as convective for objects with area less than 
     startslope
  truncZconvthres = reflectivity threshold at or above which echos are classified as convective;  
     The value in Powell et al. (2016) was used for an S-band radar with a beam width of 0.91 degrees. 
     For C-band and/or larger beam width, this value will probably need to be decreased.  Rain
     type classification is most sensitive to this input.
  dBZformaxconvradius = minimum dBZ required for max_conv_radius to apply; should be somewhere close 
     to truncZconvthres (i.e. within 5 dBZ, and usually greater)
  weakechothres = minimum dBZ for classification as not weak echo; don't change this without a good 
     reason.  7 is about as low as we can go without getting into Bragg scatter territory.
  backgrndradius (km) = radius within which background reflectivity is computed
  maxConvRadius (km) = maximum radius around convective core for possible mixed classification; 
     Powell et al. (2016) tested 5, and showed that too much convection was included 
     in stratiform region.  Don't lower this number without a good reason.
  minsize (km^2) = minimum areal coverage a contiguous echo can cover and still receive an ISO_CONV
     classification (See dBZcluster for use)
  startslope (km^2) = any contiguous echo object with areal coverage greater than this but less than 
     maxsize gets a new convective threshold that is linearly interpolated between shallowconvmin and 
     truncZconvthres depending on where between startslope and maxsize its area is (See makedBZcluster)
  maxsize (km^2) = any contiguous echo object greater than this size gets a convective threshold of 
     truncZconvthres (See makedBZcluster)
"""

## ***************** ALGORITHM USER-INPUT PARAMETERS *****************

## rain type input parameters
minZdiff = 20; 
deepcoszero = 40;
shallowconvmin = 28;
truncZconvthres = 42;
dBZformaxconvradius = 45;
mindbzuse = -50
weakechothres = 7;
backgrndradius = 5;       #(in km)
maxConvRadius = 10;       #(in km)
minsize = 8;              #(in km^2)
startslope = 50;          #(in km^2)
maxsize = 2000;           #(in km^2)

#Reflectivity sweep to base rain-type map on (zero-based 0 = 0.5 deg for SPOL or 0.8 deg for Revelle). sweep_used = 0 means use the lowest elevation angle available. 
sweep_used = 0

#Names of key input variables. SPOL during DYNAMO was dual-pol, so LDR is available. 
#A clutter filter was also applied. These aren't necessary to run the algorithm, but
#if not available, make sure you use a QCed dataset. In some datasets, the "corrected"
#reflectivity is named some form of "CZ" and the original data some form of "dBZ".
#If not using LDR or clutter filter, just leave them as default strings.
reflName = 'DBZ_S'
ldrName = 'LDRH_S'      
clutterName = 'CMD_FLAG_S'

## Information about where the reflectivity data is located and where outputs should be written. Make sure your directory names end with a /.
fileDir = '../example/';
fileDirOut = './';

title = 'Rain type classification of DYNAMO SPolKa radar data in polar coordinates';
institution = 'University of Washington';
source = 'Code used https://github.com/swpowell/raintype_python_polar';
references = 'http://www.atmos.uw.edu/MG/PDFs/JTECH16_Powell-etal_RainCat.pdf';

## *****************  END USER INPUT PARAMETERS *****************

## *****************  BEGIN OUTPUT CONSTANTS *****************
# Output constants: Do not change these without a good reason!
CS_CORE = 8;          #For convective core scheme.
ISO_CS_CORE = 9;      #For isolated convection scheme.
NO_SFC_ECHO = 0;
STRATIFORM = 1;
CONVECTIVE = 2;
MIXED = 3;
#Many users may want to set ISO_CONV_CORE and ISO_CONV_FRINGE to the same value because the core and
#fringe categories are closely related. Or such users can leave this code as-is and process the output
#as if the two categories were the same category.
ISO_CONV_CORE = 4;
ISO_CONV_FRINGE = 5;
WEAK_ECHO = 6;

## ***************** END OUTPUT CONSTANTS   ******************

# We need to make sure that the convsf file starts at sweep_number zero
# with start_ray_index = 0 and end_ray_index = num rays in sweep minus one
# This may need modification if we output multi-sweep convsf file.
sweep_number_out = 0;
start_index_out = 0;

sdir = os.listdir(fileDir)
numfiles = len(sdir)

#This for-loop assumes the code is being run on a batch of several input files. Each file will be
#processed in the order it appears in the directory.
for m in range(0,numfiles):

    try:

    #Filename for input
        fname = fileDir + sdir[m]

        #Read in CF/Radial file
        (sls_size,sweep_number1,volume_number1,time_coverage_start1,time_coverage_end1,starttime,lat1,lon1,alt1,sweep_mode1,fixed_angle1,time1,range1,azimuth1,elevation1,ssri1,seri1,meterstoFirstGate,metersBetweenGates,kmToFirstGate,kmBetweenGates,sweep_start_ray_index,sweep_end_ray_index,ins_name,ncid) = io.readcfrad(fname,sweep_used,reflName,ldrName,clutterName)

        #Read designated sweep
        (dBZsweep,numRanges,numTimes) = io.readsweep(ncid,sweep_start_ray_index,sweep_end_ray_index,sweep_used,fixed_angle1,reflName,ldrName,clutterName)

        #Set up masks for background and mixed region + compute background reflectivities
        (maskcell,convcell,background,sectorarea,dBZsweep,minR,maxR) = alg.convsf(kmToFirstGate,kmBetweenGates,numRanges,numTimes,backgrndradius,maxConvRadius,sweep_used,dBZsweep)

        #Run the algorithm. 
        rtfill = -99 #Set fill value for rain-type (mainly for outer ring of unclassified data)
        raintype = alg.convectivecore(background,dBZsweep,minZdiff,CS_CORE,ISO_CS_CORE,CONVECTIVE,STRATIFORM,MIXED,WEAK_ECHO,ISO_CONV_CORE,ISO_CONV_FRINGE,NO_SFC_ECHO,dBZformaxconvradius,maxConvRadius,weakechothres,deepcoszero,minsize,maxsize,startslope,shallowconvmin,truncZconvthres,mindbzuse,sectorarea,convcell,maxR,numRanges,numTimes,rtfill)

        #Write the data to a new CF/Radial file that is viewable in CIDD.  
        io.writecfrad(fileDirOut,sdir[m],raintype,sls_size,volume_number1,time_coverage_start1,time_coverage_end1,lat1,lon1,alt1,sweep_number1,sweep_mode1,sweep_used,fixed_angle1,ssri1,seri1,time1,numTimes,range1,meterstoFirstGate,metersBetweenGates,azimuth1,elevation1,rtfill,starttime,ins_name,title,institution,source,references)

    except:
        Warning(str('Something went wrong processing this file! ' + fileDir + sdir[m])) 
