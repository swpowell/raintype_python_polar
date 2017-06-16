from __future__ import division     #For python2 users only.

def pol2cart(phi,rho):
  #Just a simple code to convert (azimuth, radius) coordinates to (X,Y) coordinates.
  import numpy as np
  x = rho * np.cos(phi)
  y = rho * np.sin(phi)
  return(x, y)

#*********End pol2cart*******************


def radialdistancemask(Xpt,Ypt,X,Y,radius,convradius):
  import numpy as np
  #csmask is for the mask for MIXED classifications.
  csmask = np.empty([X.shape[0],X.shape[1]])
  csmask[:] = np.nan

  #bgmask is for the mask for background reflectivity calculations.
  bgmask = np.empty([X.shape[0],X.shape[1]])
  bgmask[:] = np.nan

  #Distance formula to determine how far away surrounding points are.
  dtemp = np.sqrt((Xpt-X)**2 + (Ypt-Y)**2)
  dtemp = np.maximum(dtemp,0)[:]

  #If distances are close enough to be masked, then indicate such in csmask or bgmask.
  csmask[dtemp <= convradius] = 1
  csmask[dtemp <= convradius-1] = 2
  csmask[dtemp <= convradius-2] = 3
  csmask[dtemp <= convradius-3] = 4
  csmask[dtemp <= convradius-4] = 5
  bgmask[dtemp <= radius] = 1

  return(csmask,bgmask)


#********End radialdistancemask***************


def makedBZcluster(refl,isCore,convsfmat,weakechothres,minsize,maxsize,startslope,shallowconvmin,truncZconvthres,ISO_CONV_FRINGE,WEAK_ECHO,ISO_CS_CORE,CS_CORE,sectorarea,numTimes):
  import numpy as np
  from scipy import ndimage as nd

  #Allocate matrix indicating whether rain is occurring.
  rain = np.zeros((refl.shape),dtype=np.int)

  #If echo is strong enough, rain = 1.
  rain[refl>=weakechothres] = 1

  halfnumTimes = np.int16(np.ceil(0.5*numTimes))

  #Wrap matrices around so 0 and 359 degrees are continuous.
  rain = np.concatenate((rain[:,halfnumTimes:numTimes+1],rain,rain[:,0:halfnumTimes]),axis=1)
  isCore_minsize = np.concatenate((isCore[:,halfnumTimes:numTimes+1],isCore,isCore[:,0:halfnumTimes]),axis=1)
  convsfmat_minsize = np.concatenate((convsfmat[:,halfnumTimes:numTimes+1],convsfmat,convsfmat[:,0:halfnumTimes]),axis=1)

  #Create truncvalue, which has same shape as reflectivity data. It indicates the 
  #reflectivity over which an echo is automatically classified as some sort of 
  #ISOLATED CONVECTIVE echo. See details below.
  truncvalue = np.ones((isCore_minsize.shape),dtype=np.float64)*truncZconvthres

  #This is a blob detector. Detects contiguous areas of raining pixels. Diagonally
  #touching pixels that share a corner don't count. Edges must touch.
  #echoes contains the blob objects, numechoes is just a count of them.
  (echoes,numechoes) = nd.label(rain)

  for i in range(0,numechoes):

    #Find 2D indices of echo object.
    (I,J) = np.where(echoes==i+1)
    
    #Compute the total areal coverage of the echo object.
    clusterarea = np.nansum(sectorarea[I,J])    #In km^2

    #Any echo object with a size between minsize and maxsize is considered 
    #ISOLATED CONVECTION. First, make all of it FRINGE.
    if clusterarea >= minsize and clusterarea <= maxsize:
      convsfmat_minsize[I,J] = ISO_CONV_FRINGE 

    #Very small echo objects are dismissed as WEAK ECHO.  
    if clusterarea < minsize:
      isCore_minsize[I,J] = 0
      convsfmat_minsize[I,J] = WEAK_ECHO
    #Echo objects with size between minsize and startslope get a small truncvalue
    #equal to shallowconvmin.
    elif clusterarea >= minsize and clusterarea < startslope:
      truncvalue[I,J] = shallowconvmin
    #Echo objects with size between startslope and maxsize get a truncvalue that 
    #is linearly interpolated between shallowconvmin and truncZconvthres depending
    #on the size relative to startslope and maxsize.
    elif clusterarea >= startslope and clusterarea <= maxsize:
      truncvalue[I,J] = shallowconvmin + ((clusterarea-startslope)/(maxsize-startslope))*(truncZconvthres-shallowconvmin)

  #Unwrap and send variables back to original size.
  truncvalue = truncvalue[:,halfnumTimes:halfnumTimes+numTimes]
  isCore = isCore_minsize[:,halfnumTimes:halfnumTimes+numTimes]
  convsfmat = convsfmat_minsize[:,halfnumTimes:halfnumTimes+numTimes]

  #Evaluate isCore with size of echo object accounted for.
  #First, if reflectivity exceeds truncvalue, classify it as ISOLATED CONVECTIVE CORE.
  isCore[refl >= truncvalue] = ISO_CS_CORE

  #But if reflectivity exceeds original reflectivity threshold, classify as CONVECTIVE core.
  isCore[ (refl >= truncZconvthres)*(convsfmat != ISO_CONV_FRINGE)==1 ] = CS_CORE

  return (convsfmat,isCore)
