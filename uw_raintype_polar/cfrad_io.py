def readcfrad(fname,sweep_used,reflName,ldrName,clutterName):

    import netCDF4 as nc4
    import numpy as np
    import warnings

    #Check to make sure file isn't an RHI file. 
    ncid = nc4.Dataset(fname,'r')
#    if max(ncid.variables['elevation'][:]) > 40:
#      Warning('This file might be a misplaced RHI file!')
#      print(fname)

    #Read in some variables just for writing out to file. A few will be called in the algorithm.
    sls_size = 32
    starttime = 'None'
    sweep_number1 = ncid.variables['sweep_number'][sweep_used]
    volume_number1 = ncid.variables['volume_number']
    time_coverage_start1 = ncid.variables['time_coverage_start']
    #This for-loop is for CIDD to correctly display time and date. Sometimes this field gets read
    #in as a masked array and sometimes it doesn't.
    for i in range(0,len(time_coverage_start1)):
      if i == 0:
        starttime = time_coverage_start1[i]
      elif i < 20: #20 because starttime should always take the form YYYY-MM-DDTHH:MM:SSZ
        starttime = starttime + time_coverage_start1[i]
      else:
        break
    time_coverage_end1 = ncid.variables['time_coverage_end']
    lat1 = ncid.variables['latitude']
    lon1 = ncid.variables['longitude']
    alt1 = ncid.variables['altitude']
    sweep_mode1 = ncid.variables['sweep_mode']
    fixed_angle1 = ncid.variables['fixed_angle']
    time1 = ncid.variables['time']
    range1 = ncid.variables['range']
    azimuth1 = ncid.variables['azimuth']
    elevation1 = ncid.variables['elevation']
    ssri1 = ncid.variables['sweep_start_ray_index']
    seri1 = ncid.variables['sweep_end_ray_index']
    ins_name = ncid.instrument_name

    #These are some values that we need for the algorithm. 
    meterstoFirstGate = range1.meters_to_center_of_first_gate
    metersBetweenGates = range1.meters_between_gates
    kmToFirstGate = meterstoFirstGate/1000
    kmBetweenGates = metersBetweenGates/1000
    sweep_start_ray_index = ssri1[:]
    sweep_end_ray_index = seri1[:]

    return(sls_size,sweep_number1,volume_number1,time_coverage_start1,time_coverage_end1,starttime,lat1,lon1,alt1,sweep_mode1,fixed_angle1,time1,range1,azimuth1,elevation1,ssri1,seri1,meterstoFirstGate,metersBetweenGates,kmToFirstGate,kmBetweenGates,sweep_start_ray_index,sweep_end_ray_index,ins_name,ncid)


#**********************End readcfrad********************


def readsweep(ncid,sweep_start_ray_index,sweep_end_ray_index,sweep_used,fixed_angle1,reflName,ldrName,clutterName):
    
    import netCDF4 as nc4
    import numpy as np

    #Get some basic information and allocate space for reflectivity field along the chosen sweep.
    starttimes = ncid.variables['ray_start_index'][sweep_start_ray_index[sweep_used]:sweep_end_ray_index[sweep_used]+1]
    lengths = ncid.variables['ray_n_gates'][sweep_start_ray_index[sweep_used]:sweep_end_ray_index[sweep_used]+1]
    numTimes = len(lengths)
    numRanges = len(ncid.variables['range'][:])
    dBZsweep = np.empty([numRanges,numTimes])
    dBZsweep[:] = np.nan
    
    #Get fill values for missing data.
    refl_fill_value = ncid.variables[reflName]._FillValue

    #Get fill values for missing LDR and clutter data, if the fields are present.
    try:
      ldrh_fill_value = ncid.variables[ldrName]._FillValue
    except:
      dummy = 1
      del dummy
    try:
      clutter_fill_value = ncid.variables[clutterName]._FillValue
    except:
      dummy = 1
      del dummy
    
    #Allocate space for LDR and clutter, if the fields are present.
    if 'ldrh_fill_value' in locals():
      ldrsweep = np.empty([numRanges,numTimes])
      ldrsweep[:] = np.nan
    if 'clutter_fill_value' in locals():
      csweep = np.empty([numRanges,numTimes])
      csweep[:] = np.nan

    #Get reflectivity data, and LDR and clutter data if the latter two are present.
    for phi in range(0,len(starttimes)):
      dBZsweep[0:lengths[phi],phi] = ncid.variables[reflName][starttimes[phi]:starttimes[phi]+lengths[phi]]
      if 'ldrsweep' in locals():
        ldrsweep[0:lengths[phi],phi] = ncid.variables[ldrName][starttimes[phi]:starttimes[phi]+lengths[phi]]
      if 'csweep' in locals():
        csweep[0:lengths[phi],phi] = ncid.variables[clutterName][starttimes[phi]:starttimes[phi]+lengths[phi]]
    del phi

    #Change missing reflectivity data to NaN.
    dBZsweep[dBZsweep == refl_fill_value] = np.nan
    #If LDR data is available, NaN out any second-trip echo.
    if 'ldrsweep' in locals():
      dBZsweep[csweep == 1] = np.nan
    #If clutter filter is available, NaN out any echoes flagged as clutter.
    if 'csweep' in locals():
      dBZsweep[ldrsweep > 0] = np.nan

    #Get the intended elevation angle of this sweep.
    theta = fixed_angle1[sweep_used]

    return(dBZsweep,numRanges,numTimes)


#**********************End readsweep********************


def writecfrad(fileDirOut,sdir,raintype,sls_size,volume_number1,time_coverage_start1,time_coverage_end1,lat1,lon1,alt1,sweep_number1,sweep_mode1,sweep_used,fixed_angle1,ssri1,seri1,time1,numTimes,range1,meterstoFirstGate,metersBetweenGates,azimuth1,elevation1,rtfill,starttime,ins_name,title,institution,source,references):

    import netCDF4 as nc4
    import numpy as np
    
    #Set the name of the output file
    ncname = fileDirOut+'raintype.'+sdir

    #Create the file
    ncid = nc4.Dataset(ncname,'w',format='NETCDF4')

    #Set global attributes
    ncid.Conventions = "CF/Radial instrument_parameters radar_parameters radar_calibration"
    ncid.title = title
    ncid.source = source
    ncid.institution = institution 
    ncid.references = references
    ncid.comment = "NO ECHO = 0, STRATIFORM = 1, CONVECTIVE = 2, MIXED = 3, ISOLATED CONVECTIVE CORE = 4, ISOLATED CONVECTIVE FRINGE = 5, WEAK ECH0 = 6"
    ncid.instrument_name = str(ins_name)

    #Create dimensions
    t = ncid.createDimension('time',raintype.shape[1])
    r = ncid.createDimension('range',raintype.shape[0])
    sweep = ncid.createDimension('sweep',1)
    string_length_short = ncid.createDimension('string_length_short',sls_size)

    #Write the data and variable attributes.
    vn = ncid.createVariable('volume_number',np.int32,fill_value=-9999)
    vn[:] = volume_number1[:]
    vn.standard_name = "data_volume_index_number"
    tcs = ncid.createVariable('time_coverage_start','S1',('string_length_short'))
    tcs[:] = time_coverage_start1[:]
    tcs.standard_name = "data_volume_start_time_utc"
    tcs.comment = "ray times are relative to start time in secs"
    tce = ncid.createVariable('time_coverage_end','S1',('string_length_short'))
    tce[:] = time_coverage_end1[:]
    tce.standard_name = "data_volume_end_time_utc"
    lat = ncid.createVariable('latitude',np.double)
    lat[:] = lat1[:]
    lat.standard_name = "latitude"
    lat.units = "degrees_north"
    lon = ncid.createVariable('longitude',np.double)
    lon[:] = lon1[:]
    lon.standard_name = "longitude"
    lon.units = "degrees_east"
    alt = ncid.createVariable('altitude',np.double)
    alt[:] = alt1[:]
    alt.standard_name = "altitude"
    alt.units = "meters"
    alt.positive = "up"
    sweep_number = ncid.createVariable('sweep_number',np.int32,('sweep'),fill_value=-9999)
    sweep_number[:] = sweep_number1
    #sweep_number.standard_name = "sweep_index_number_0_based"
    #Above line throws an error. Don't know why. Maybe give it long_name instead?
    sweep_mode = ncid.createVariable('sweep_mode','S1',('sweep','string_length_short'))
    sweep_mode[:] = sweep_mode1[sweep_used,:]
    sweep_mode.standard_name = "scan_mode_for_sweep"
    sweep_mode.options = "sector, coplane, rhi, vertical_pointing, idle, azimuth_surveillance, elevation_surveillance, sunscan, pointing, calibration, manual_ppi, manual_rhi"
    fixed_angle = ncid.createVariable('fixed_angle',np.float,('sweep'),fill_value=-9999)
    fixed_angle[:] = fixed_angle1[sweep_used]
    fixed_angle.standard_name = "beam_target_fixed_angle"
    fixed_angle.units = "degrees"
    ssri = ncid.createVariable('sweep_start_ray_index',np.int32,('sweep'),fill_value=-9999)
    ssri[:] = ssri1[sweep_used]
    ssri.standard_name = "index_of_first_ray_in_sweep"
    seri = ncid.createVariable('sweep_end_ray_index',np.int32,('sweep'),fill_value=-9999)
    seri[:] = seri1[sweep_used]
    seri.standard_name = "index_of_last_ray_in_sweep"
    timevar = ncid.createVariable('time',np.double,('time'))
    timevar[:] = time1[sweep_used*numTimes:(sweep_used+1)*numTimes]
    timevar.standard_name = "time"
    timevar.long_name = "time in seconds since volume start"
    try: #Python 2?
      timevar.units = "seconds since " + starttime
    except: #Python 3?
      timevar.units = "seconds since " + str(starttime)[2:22]
    timevar.comment = "times are relative to volume start time"
    rangevar = ncid.createVariable('range',np.float,('range'))
    rangevar[:] = range1[:]
    rangevar.standard_name = "range_to_center_of_measurement_volume"
    rangevar.long_name = "Range from instrument to center of gate"
    rangevar.units = "meters"
    rangevar.spacing_is_constant = "True"
    rangevar.meters_to_center_of_first_gate = str(meterstoFirstGate)
    rangevar.meters_between_gates = str(metersBetweenGates)
    azi = ncid.createVariable('azimuth',np.float,('time'),fill_value = -9999)
    azi[:] = azimuth1[sweep_used*numTimes:(sweep_used+1)*numTimes]
    azi.standard_name = "beam_azimuth_angle"
    azi.units = "degrees"
    elev = ncid.createVariable('elevation',np.float,('time'),fill_value = -9999)
    elev[:] = elevation1[sweep_used*numTimes:(sweep_used+1)*numTimes]
    elev.standard_name = "beam_elevation_angle"
    elev.units = "degrees"
    elev.positive = "up"
    finalrt = ncid.createVariable('raintype',np.int32,('time','range'),fill_value=rtfill )
    finalrt[:,:] = (np.transpose(raintype))[:]
    finalrt.long_name = "rain type classification"
    finalrt.units = "unitless"

    ncid.close()
