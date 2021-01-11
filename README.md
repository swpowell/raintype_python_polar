# raintype_python_polar V1.0 (for native polar grids)

Authors: Scott W. Powell and Stacy Brodzik, University of Washington
Date: October 2016
Email: scott.powell@nps.edu, brodzik@atmos.uw.edu
This is the rainfall type categorization (formerly known as convective/stratiform classification) of Powell et al. (2016, JTECH). It is an update to the Steiner et al. (1995) method of classifying convective and stratiform echoes of tropical and subtropical precipitation. This is the version that runs on the native gridding (the polar coordinate system) of radar output. If you would like to run this algorithm on interpolated radar data in Cartesian coordinates, see the other version for this code (raintype_python).

----------------------------------------------------------------

Setting up and using the code:

Before you do anything, you'll need to make sure that you have numpy, scipy, and netcdf4-python (plus dependencies) installed on your machine. netcdf4-python can be found, as of 2016, here: https://github.com/Unidata/netcdf4-python

The polar version of the code has no modules to install locally. Running runraintype.py with the correct input should work immediately. There is a test file and in subdirectory "example" that you can use to ensure the installation went smoothly. Just go to the directory uw_raintype_polar and run

>> python -W ignore runraintype.py 

It should create a CfRadial output file of rain-type classifications for the example file.

----------------------------------------------------------------

Description of files:

runraintype.py: This is the driver/wrapper code. This is the code that should be modified by the user. User input parameters are listed and described at the top of this code. Input files must be in CfRadial format! If your data is not in CfRadial format, you can usually use RadX (https://www.ral.ucar.edu/projects/titan/docs/radial_formats/radx.html) to make the conversion quite easily. A good tutorial on basic RadX use can be found at Steve Nesbitt's website (https://www.ral.ucar.edu/projects/titan/docs/radial_formats/radx.html).

rtfunctions.py: Contains a variety of functions for implementing algorithm.

algorithm.py: The rain-type classification algorithm.
 
cfrad_io.py: Deals with input and output of CfRadial data.

For basic users, after the user input is tuned appropriately (see below), the code can be executed by entering

>> python -W ignore runraintype.py

The -W ignore flag suppresses warnings that will otherwise pop up because numpy is trying to compare NaNs to real numbers. Don't worry about these warnings when running the code. You'll want to suppress them so they don't keep displaying to the terminal and slowing down your code.

----------------------------------------------------------------

Tuning the user input:

Furthermore, and I cannot stress this enough, the values that are included as "default" in the code (specifically in runraintype.py) are probably not appropriate for your purposes. They need to be tuned based on the radar platform used, the beam width used, the convective regime sampled, etc. The "default" parameters in the code were appropriate for maritime tropical convection observed with an S-band (10 cm wavelength) radar with 0.91 deg beam width (S-PolKa during DYNAMO 2011 in the Central Eq. Indian Ocean).

The classification is *most* sensitive to the selection of the convective threshold, truncZconvthres. For radars with larger wavelengths (like C-band), or for radars using larger beam widths, the value of truncZconvthres will probably need to be reduced. I recommend keeping dBZformaxconvradius within 5 dBZ of truncZconvthres.

----------------------------------------------------------------

Output:

The output is in NetCDF format and is written on the same grid as the reflectivity data used as input.

Values are as follows:

0 = No Echo or Discarded Clutter Echo
1 = Stratiform
2 = Convective
3 = Mixed
4 = Isolated Convective Core
5 = Isolated Convective Fringe
6 = Weak Echo 

For analysis, ignore everything classified as "0" or "6"  unless you have a good reason for specifically examining weak echo. 

The "Mixed" category represents echoes surrounding convective cores (Classification 2). In the old convective/stratiform classification algorithm, these echoes were considered convective, but Powell et al. (2016) shows that the heating profile near convective echoes cannot be distinctively classified as either convective or stratiform based on the distance from a convective core. 

The two isolated convective categories largely contain shallow convective elements. In the old algorithm, convective cores of such elements were usually classified as convective, but the echo surrounding the cores was erroneously classified as stratiform. What Steiner et al. (1995) classified as stratiform is mostly now considered Isolated Convective Fringe (Classification 5). The echoes have composite heating profiles that are consistent with shallow/weak convective echoes, but the shape of droplets (based on ZDR profiles in such echoes) are more stratiform in nature, consistent with the idea that hydrometeors in such echoes consist primarily of "fallout" from nearby convection.

If trying to estimate rainfall in mixed echoes with a method that depends on convective/stratiform classification, it is probably best to use a Z-R (or Z-ZDR-R, etc.) relationship that has been derived for all (convective + stratiform) echoes. If you wish to express a conservative range of potential estimates, you may treat the mixed region as all convective in one estimate and all stratiform in another. It is recommended that Isolated Convective Fringe echoes be treated with a Z-R, etc. relationship derived from stratiform regions. Convective and Isolated Convective Core may be treated with convective Z-R, etc. relationships, and Stratiform echoes, obviously, with a stratiform relationship.

**Never try to interpolate the rain-type output in polar coordinates onto a Cartesian grid! If you want to do that, just use the Cartesian version of this code (i.e. interpolate the radar data, then run the algorithm).**

----------------------------------------------------------------

Known Issues: 

**It is not recommended that this algorithm be used if the user cares about the high-frequency variability (time-scales of approximately less than or equal to 1 day) of convection in their radar domain.** Between two temporally consecutive radar volumes, sometimes echo objects will "flip/flop" between being classified as Convective+Stratiform to being classified as Isolated Convective Core+Fringe. This happens because echoes are only considered "Isolated" if the echo object they reside in is less than maxsize (set to 2000 km by default). If echo objects vary in size around this threshold, their classifications may change back and forth between radar volumes. This is only a problem if you care about the high-frequency evolution of convection in your radar domain. A similar, but less obvious issue exists for echoes with reflectivity that are very near truncZconvthres; such echoes may "flip-flop" between Convective and Stratiform or Mixed classifications. The latter was an issue even in Steiner et al. (1995). 

**If you are interested in variability of Isolated echo classes, restrict your analysis to the innermost 75-125 km of the radar domain.** Because isolated echoes are often shallow, such echoes are more likely to be missed farther from a radar site. Consider that the center of a 0.5 deg scan's beam will be around 2 km in altitude at 150 km from the radar. Most, but not all precipitating elements should exceed this height, but the reflectivities of 2-3 km deep objects near their tops might not be large enough to yield a proper Isolated Convective Core classification.

