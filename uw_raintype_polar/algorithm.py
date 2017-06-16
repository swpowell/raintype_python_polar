from __future__ import division     #For python2 users only.

def convsf(kmToFirstGate,kmBetweenGates,numRanges,numTimes,backgrndradius,maxConvRadius,sweep_used,dBZsweep,filenum,maskcell):

    #Purpose: To make background and MIXED region masks and to compute background reflectivity.

    import numpy as np
    import rtfunctions as rt

    halfnumTimes = np.int16(np.ceil(0.5*numTimes))
    #Make a map of the (azimuth, radius) coordinate system.
    [phi,R] = np.meshgrid( np.linspace(0,numTimes-1,numTimes)*np.pi/halfnumTimes, np.linspace(kmToFirstGate,numRanges*kmBetweenGates+kmToFirstGate,numRanges))
    #Convert this map to Cartesian coordinates (an irregularly spaced grid) for calculating distances.
    [X,Y] = rt.pol2cart(phi,R)
    del phi

    #Any data within minR of the radar site will get NaNed out.
    minR = int(round(0.125/kmBetweenGates))

    #Any data beyond maxR of the radar size will also get NaNed out.
    if sweep_used == 0:
      maxR = int(round(numRanges-0.5*backgrndradius/kmBetweenGates))
    else:
      minR = int(round(5/kmBetweenGates))
      maxR = int(round(30/kmBetweenGates))

    #Create masks. One mask per ring of data (i.e. one mask per radius from center).
    #Only do this if mask doesn't yet exist. Only needs to occur on first file in batch.
    #maskcell are indices of  points within backgrndradius of a point at given radius from the radar site.
    #convcell are indices of points within maxConvRadius-5 to maxConvRadius of a convective core.
    #convcell[R,0] is the largest mask, and convcell[R,5] is the smallest mask for weaker convective echoes.
    #Any echoes that end up getting masked by convcell (in def convectivecore) will be MIXED.
    if filenum != 0:
      #do nothing
      dummy = 0
      del dummy
    else:
      centerphi = halfnumTimes  #180
      maskcell = np.empty([numRanges,1],dtype=object)
      convcell = np.empty([numRanges,5],dtype=object)
      for R in range(minR,maxR+1):
        [mask,bgmask] = rt.radialdistancemask(X[R-1,centerphi-1],Y[R-1,centerphi-1],X,Y,backgrndradius,maxConvRadius)
        [I,J] = np.where(bgmask==1)
        [I2,J2] = np.where(mask>=1)
        [I3,J3] = np.where(mask>=2)
        [I4,J4] = np.where(mask>=3)
        [I5,J5] = np.where(mask>=4)
        [I6,J6] = np.where(mask>=5)
        #To create the masks, convert the 2D matrices of data to 1D arrays. Find the indices
        #of the 1D arrays that should be masked. After masking in another subroutine
        #(convectivecore), these will be converted back to 2D matrices.
        centerval = np.ravel_multi_index((R-1,centerphi-1),dBZsweep.shape,order="F")-1
        maskcell[R] = [np.ravel_multi_index((I,J),bgmask.shape,order="F")-centerval]
        convcell[R,0] = [np.ravel_multi_index((I2,J2),mask.shape,order="F")-centerval]
        convcell[R,1] = [np.ravel_multi_index((I3,J3),mask.shape,order="F")-centerval]
        convcell[R,2] = [np.ravel_multi_index((I4,J4),mask.shape,order="F")-centerval]
        convcell[R,3] = [np.ravel_multi_index((I5,J5),mask.shape,order="F")-centerval]
        convcell[R,4] = [np.ravel_multi_index((I6,J6),mask.shape,order="F")-centerval]

      #Compute the areal coverage of each data point. (This gets larger farther from radar.)
      sectorarea = np.empty([numRanges,numTimes])
      sectorarea[:] = np.nan
      for i in range (0,numRanges-1):
        #sectorarea[i,:] = 1/360*np.pi*((kmBetweenGates*(i+1))**2-(kmBetweenGates*i)**2)
        sectorarea[i,:] = 1/numTimes*np.pi*((kmBetweenGates*(i+1))**2-(kmBetweenGates*i)**2)
      #sectorarea = np.concatenate((sectorarea[:,180:360],sectorarea,sectorarea[:,0:180]),axis=1)
      sectorarea = np.concatenate((sectorarea[:,halfnumTimes:numTimes+1],sectorarea,sectorarea[:,0:halfnumTimes]),axis=1)
  
    dBZsweep[1:minR+1,:] = np.nan       #NaN out reflectivity close to radar.
    Zsweep = 10**(0.1*dBZsweep)         #Convert dBZ to Z.

    #Compute background reflectivity at each point.

    #Allocate memory.
    background = np.empty([numTimes,numRanges])
    background[:] = np.nan

    #Wrap data around so that data at 0 and 359 degrees are continuous.
    phi = np.array(range(halfnumTimes,halfnumTimes+numTimes))
    Zsweepconcat = np.concatenate((Zsweep[:,halfnumTimes:numTimes+1],Zsweep,Zsweep[:,0:halfnumTimes]),axis=1)
    Zsweepconcat[np.isnan(Zsweepconcat)] = 0
    Zsweepconcat = np.reshape(Zsweepconcat,(np.size(Zsweepconcat),1), order="F")

    #Use the background mask created above to compute background Z.
    for R in range(minR,maxR+1):
      maskuse = maskcell[R][0]
      background[0:numTimes,R] = np.mean(Zsweepconcat[maskuse[:,np.newaxis]+numRanges*(phi)+R],0)[:,0]

    #Just clearing memory.
    del Zsweepconcat

    #Convert background Z to dBZ. 
    background[background == 0] = np.nan
    background = np.transpose(10*np.log10(background))

    if filenum == 0:
      return(maskcell,convcell,background,sectorarea,dBZsweep,minR,maxR)
    else:
      return(background,dBZsweep,minR,maxR)


#*********End mask production and background calculation*********


def convectivecore(background,refl,minZdiff,CS_CORE,ISO_CS_CORE,CONVECTIVE,STRATIFORM,MIXED,WEAK_ECHO,ISO_CONV_CORE,ISO_CONV_FRINGE,NO_ECHO,dBZformaxconvradius,maxConvRadius,weakechothres,deepcoszero,minsize,maxsize,startslope,shallowconvmin,truncZconvthres,mindbzuse,sectorarea,convcell,maxR,numRanges,numTimes,rtfill):

  import numpy as np
  import rtfunctions as rt
 
  #Allocate isCore, a matrix that contains whether a grid point contains a convective core
  #and convsfmat, what will ultimately be the final rain-type classification.
  isCore = np.ones(background.shape)
  convsfmat = 10*np.ones((refl.shape),dtype=np.int)

  #Allocate zDiff, the variable representing the excess over the background dBZ
  #an echo must achieve to be considered a convective core.
  zDiff = np.empty(refl.shape)
  zDiff[:] = np.nan

  #Compute zDiff. 
  zDiff = 2.5 + minZdiff * np.cos((np.pi)*background/(2*deepcoszero))
  zDiff[(background < 0)] = minZdiff 

  #If reflectivity exceeds background dBZ by zDiff, then echo is convective core.
  isCore[(refl-background >= zDiff)] = CS_CORE;

  #No chance of weak echoes being convective cores.
  isCore[(refl < weakechothres)] = 0

  #Run the shallow, isolated convective core algorithm to detect small echoes that were
  #often identified as STRATIFORM by Steiner et al. (1995)
  (convsfmat,isCore) = rt.makedBZcluster(refl,isCore,convsfmat,weakechothres,minsize,maxsize,startslope,shallowconvmin,truncZconvthres,ISO_CONV_FRINGE,WEAK_ECHO,ISO_CS_CORE,CS_CORE,sectorarea,numTimes)

  #Make initial guesses of classifications. There may be some redundancy in this code,
  #later, but these operations are fast, I think. Better safe than sorry.
  convsfmat[(isCore == CS_CORE)] = CONVECTIVE
  convsfmat[(isCore == ISO_CS_CORE)] = ISO_CONV_CORE
  convsfmat[(isCore == 0)] = WEAK_ECHO
  convsfmat[(convsfmat == 10)] = STRATIFORM
  convsfmat[(np.isnan(refl) == True)] = NO_ECHO
  convsfmat[(refl < weakechothres)] = WEAK_ECHO
  convsfmat[(refl < mindbzuse)] = NO_ECHO

  #Now assign MIXED radius to each core. Currently assumes all echoes within 
  #maxConvRadius - 4 km are MIXED classification. Stronger echoes have larger 
  #mixed radius. Mixed radii of 6-10 km appear to be supported by algorithm 
  #testing on WRF output as seen in Powell et al. (2016). 

  #Compute what the mixed radius is as a function of echo intensity.
  convRadiuskm = np.empty(refl.shape)
  convRadiuskm[:] = np.nan
  convRadiuskm[(background <= dBZformaxconvradius - 15 )] = maxConvRadius - 4
  convRadiuskm[(background > dBZformaxconvradius - 15 )] = maxConvRadius - 3 
  convRadiuskm[(background > dBZformaxconvradius - 10 )] = maxConvRadius - 2
  convRadiuskm[(background > dBZformaxconvradius - 5 )] = maxConvRadius - 1
  convRadiuskm[(background >= dBZformaxconvradius)] = maxConvRadius

  ##Assign MIXED classification to pixels near convective cores.
  
  #If data is too close to the edge throw it out. 
  isCore[maxR:numRanges,:] = 0      
 
  #Find 2D indices of convective cores.
  (I,J) = np.where(isCore==CS_CORE)
 
  for k in range(0,len(I)):
    #Find 1D index of convective core.
    indexhold = np.ravel_multi_index((I[k],J[k]),convsfmat.shape,order="F")

    #Add indices from mask (convcell) to indexhold.
    maskind = indexhold + convcell[I[k],int(maxConvRadius-convRadiuskm[I[k],J[k]])][0]

    #Wrap data near 0 or 359 degrees around for continuity.
    maskind[maskind < 0] = convsfmat.size + maskind[maskind < 0]
    maskind[maskind > convsfmat.size-1] = maskind[maskind > convsfmat.size-1] - convsfmat.size

    #Convert masked data points back to 2D indices
    (K,L) = np.unravel_index(maskind,convsfmat.shape,order='F')

    #Save indices for points that were previously CONVECTIVE or ISOLATED CONVECTIVE.
    (K1,L1) = np.where(convsfmat == CONVECTIVE)
    (K2,L2) = np.where(convsfmat == ISO_CONV_CORE)
    (K3,L3) = np.where(convsfmat == ISO_CONV_FRINGE)

    #Make masked points MIXED.
    convsfmat[K,L] = MIXED

    #Lots of data is made MIXED. Return the cores that aren't mixed back to 
    #their previous classifications.
    convsfmat[K1,L1] = CONVECTIVE
    convsfmat[K2,L2] = ISO_CONV_CORE
    convsfmat[K3,L3] = ISO_CONV_FRINGE
    del(maskind,indexhold)
 
  #Reset the 2D indices for next time. Probably not necessary.
  del(I,J)

  #Make sure original convective cores are CONVECTIVE.
  convsfmat[isCore == CS_CORE] = CONVECTIVE

  #If there is no data, or reflectivity is very low, classify as NO ECHO.
  convsfmat[np.isnan(refl)==1] = NO_ECHO
  convsfmat[refl < mindbzuse] = NO_ECHO
  
  #Classify WEAK_ECHO
  convsfmat[refl < weakechothres] = WEAK_ECHO

  #Any classifications on outer ring of data beyond maxR will be considered "missing".
  convsfmat[maxR:numRanges,:] = rtfill

  #Format the output as integers.
  convsfmat.astype(int)

  return convsfmat
