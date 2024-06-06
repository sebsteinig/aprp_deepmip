import numpy as np
import netCDF4 as nc4

#Main function to run, for an individual model (which model is specified via the "dataPaths" arguments).
#This version assumes variable names and dimensions following CMIP convention.
#
#Inputs:
#   dataPaths1:  dictionary of paths to the netCDF output for time period 1
#   firstMonth1: first month (indexed from beginning of output) for time period 1--note Python indices start with 0
#   lastMonth1:   last month (indexed from beginning of output) for time period 1
#   dataPaths2:  dictionary of paths to the netCDF output for time period 2
#                (if the two states being compared are different times from the same run, make this the same as dataPaths1)
#   firstMonth2: first month (indexed from beginning of output) for time period 2
#   lastMonth2:   last month (indexed from beginning of output) for time period 2
#
#Outputs: a netcdf file contains APRP results

#=================================================================
# Original Author: Rick Russotto
# Download from  : https://github.com/rdrussotto/pyAPRP
# Adapted for DeepMIP simulations by Jiang Zhu (jiangzhu@ucar.edu) and Sebastian Steinig
#=================================================================

def aprp_main(dataPaths1, firstMonth1, lastMonth1, dataPaths2, firstMonth2, lastMonth2, out_file):

    #Load files and run calculations for first time period
    dict1A  = loadNetCDF(dataPaths1, firstMonth1, lastMonth1)
    dict2A  = loadNetCDF(dataPaths2, firstMonth2, lastMonth2)

    # calculate aprp
    dict1B  = parameters(dict1A)
    dict2B  = parameters(dict2A)
    dictC   = d_albedo(dict1A, dict1B, dict2A, dict2B)

    # output results
    out_nc      = nc4.Dataset(out_file, 'w')
    out_nc.set_auto_mask(True)

    lat_in      = dict1A['lat'] [:]
    lon_in      = dict1A['lon'] [:]
    time_in     = dict1A['time'][:]

    dummy       = out_nc.createDimension('lat' , len(lat_in))
    dummy       = out_nc.createDimension('lon' , len(lon_in))
    dummy       = out_nc.createDimension('time', len(time_in))

    lat         = out_nc.createVariable('lat'    , 'f4', ('lat',))
    lon         = out_nc.createVariable('lon'    , 'f4', ('lon',))
    time        = out_nc.createVariable('time'   , 'f4', ('time',))
    tas1        = out_nc.createVariable('tas1'   , 'f4', ('time','lat','lon'))
    tas2        = out_nc.createVariable('tas2'   , 'f4', ('time','lat','lon'))
    cld         = out_nc.createVariable('cld'    , 'f4', ('time','lat','lon'))
    cld_c       = out_nc.createVariable('cld_c'  , 'f4', ('time','lat','lon'))
    cld_ga      = out_nc.createVariable('cld_ga' , 'f4', ('time','lat','lon'))
    cld_mu      = out_nc.createVariable('cld_mu' , 'f4', ('time','lat','lon'))
    clr         = out_nc.createVariable('clr'    , 'f4', ('time','lat','lon'))
    clr_ga      = out_nc.createVariable('clr_ga' , 'f4', ('time','lat','lon'))
    clr_mu      = out_nc.createVariable('clr_mu' , 'f4', ('time','lat','lon'))
    alf         = out_nc.createVariable('alf'    , 'f4', ('time','lat','lon'))
    alf_clr     = out_nc.createVariable('alf_clr', 'f4', ('time','lat','lon'))
    alf_oc      = out_nc.createVariable('alf_oc' , 'f4', ('time','lat','lon'))
    sw_net_toa  = out_nc.createVariable('sw_net_toa' , 'f4', ('time','lat','lon'))
    lw_net_toa  = out_nc.createVariable('lw_net_toa' , 'f4', ('time','lat','lon'))
    lw_cre      = out_nc.createVariable('lw_cre' , 'f4', ('time','lat','lon'))
    lat.units       = 'degrees_north'
    lon.units       = 'degrees_east'
    time.units      = 'days since 0001-01-01 00:00:00'
    lat.long_name   = 'latitude'
    lon.long_name   = 'longitude'
    time.long_name  = 'time'
    lat.axis        = 'Y'
    lon.axis        = 'X'

    time    [:] = [15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5]
    lat     [:] = dict1A['lat']       [:]
    lon     [:] = dict1A['lon']       [:]
    tas1    [:] = dict1A['tas' ]      [:]
    tas2    [:] = dict2A['tas' ]      [:]
    cld     [:] = dictC['cloud']      [:]
    cld_c   [:] = dictC['cloud_c']    [:]
    cld_ga  [:] = dictC['cloud_ga']   [:]
    cld_mu  [:] = dictC['cloud_mu']   [:]
    clr     [:] = dictC['noncloud']   [:]
    clr_ga  [:] = dictC['noncloud_ga'][:]
    clr_mu  [:] = dictC['noncloud_mu'][:]
    alf     [:] = dictC['surface']    [:]
    alf_clr [:] = dictC['surface_clr'][:]
    alf_oc  [:] = dictC['surface_oc'] [:]
    sw_net_toa[:]= dictC['sw_net_toa'] [:]
    lw_net_toa[:]= dictC['lw_net_toa'] [:]
    lw_cre[:]= dictC['lw_cre'] [:]
    out_nc.close()

# read data
def loadNetCDF(dataPaths, firstMonth, lastMonth):

    tas     = nc4.Dataset(dataPaths['tas'   ]).variables['tas'   ][firstMonth:lastMonth+1,:,:]
    clt     = nc4.Dataset(dataPaths['clt'   ]).variables['clt'   ][firstMonth:lastMonth+1,:,:]
    rsds    = nc4.Dataset(dataPaths['rsds'  ]).variables['rsds'  ][firstMonth:lastMonth+1,:,:]
    rsus    = nc4.Dataset(dataPaths['rsus'  ]).variables['rsus'  ][firstMonth:lastMonth+1,:,:]
    rsut    = nc4.Dataset(dataPaths['rsut'  ]).variables['rsut'  ][firstMonth:lastMonth+1,:,:]
    rsdt    = nc4.Dataset(dataPaths['rsdt'  ]).variables['rsdt'  ][firstMonth:lastMonth+1,:,:]
    rsutcs  = nc4.Dataset(dataPaths['rsutcs']).variables['rsutcs'][firstMonth:lastMonth+1,:,:]
    rsdscs  = nc4.Dataset(dataPaths['rsdscs']).variables['rsdscs'][firstMonth:lastMonth+1,:,:]
    rsuscs  = nc4.Dataset(dataPaths['rsuscs']).variables['rsuscs'][firstMonth:lastMonth+1,:,:]
    rlut    = nc4.Dataset(dataPaths['rlut'  ]).variables['rlut'  ][firstMonth:lastMonth+1,:,:]
    rlutcs  = nc4.Dataset(dataPaths['rlutcs'  ]).variables['rlutcs'  ][firstMonth:lastMonth+1,:,:]
    rlus    = nc4.Dataset(dataPaths['rlus'  ]).variables['rlus'  ][firstMonth:lastMonth+1,:,:]
    rlds    = nc4.Dataset(dataPaths['rlds'  ]).variables['rlds'  ][firstMonth:lastMonth+1,:,:]

    Dataset = nc4.Dataset(dataPaths['tas'   ])
    for var in Dataset.variables:
        if var == 'lat'         : lat  = Dataset.variables['lat'][:]
        if var == 'lon'         : lon  = Dataset.variables['lon'][:]
        if var == 'time'        : time = Dataset.variables['time' ][firstMonth:lastMonth+1]
        if var == 'latitude'    : lat  = Dataset.variables['latitude'][:]
        if var == 'longitude'   : lon  = Dataset.variables['longitude'][:]
        if var == 'month'       : time = Dataset.variables['month'][firstMonth:lastMonth+1]
        if var == 't'           : time = Dataset.variables['t'    ][firstMonth:lastMonth+1]
        if var == 'time_counter': time = Dataset.variables['time_counter'    ][firstMonth:lastMonth+1]

   #clt     = np.ma.masked_greater(clt   ,1.e10)
   #rsds    = np.ma.masked_greater(rsds  ,1.e10)
   #rsus    = np.ma.masked_greater(rsus  ,1.e10)
   #rsut    = np.ma.masked_greater(rsut  ,1.e10)
   #rsdt    = np.ma.masked_greater(rsdt  ,1.e10)
   #rsutcs  = np.ma.masked_greater(rsutcs,1.e10)
   #rsdscs  = np.ma.masked_greater(rsdscs,1.e10)
   #rsuscs  = np.ma.masked_greater(rsuscs,1.e10)

    #Calculate the overcast versions of rsds, rsus, rsut from the clear-sky and all-sky data
    #First mask zero values of cloud fraction so you don't calculate overcast values in clear-sky pixels
    # clt     = np.ma.masked_values(clt, 0) * 100.
    clt     = np.ma.masked_values(clt, 0)
    c       = clt/100. #c is cloud fraction. clt was in percentages
    rsdsoc  = (rsds-(1.-c)*(rsdscs))/c  #Can derive this algebraically from Taylor et al., 2007, Eq. 3
    rsusoc  = (rsus-(1.-c)*(rsuscs))/c
    rsutoc  = (rsut-(1.-c)*(rsutcs))/c

    #Mask zero values of the downward SW radiation (I assume this means polar night, for monthly mean)
    rsds    = np.ma.masked_values(rsds  , 0)
    rsdscs  = np.ma.masked_values(rsdscs, 0)
    rsdsoc  = np.ma.masked_values(rsdsoc, 0)
    rsdt    = np.ma.masked_values(rsdt  , 0)

    #Return dictionary with all the variables calculated here (called "dictA" because calculated in first function called)
    dictA = dict()
    dictA['rsds']   = rsds
    dictA['rsus']   = rsus
    dictA['rsut']   = rsut
    dictA['rsdt']   = rsdt
    dictA['rsutcs'] = rsutcs
    dictA['rsdscs'] = rsdscs
    dictA['rsuscs'] = rsuscs
    dictA['clt']    = clt
    dictA['tas']    = tas
    dictA['lat']    = lat
    dictA['lon']    = lon
    dictA['time']   = time
    dictA['rsdsoc'] = rsdsoc
    dictA['rsusoc'] = rsusoc
    dictA['rsutoc'] = rsutoc
    dictA['rlut']   = rlut
    dictA['rlutcs'] = rlutcs
    dictA['rlus']   = rlus
    dictA['rlds']   = rlds
    dictA['c'] = c #Cloud fraction as fraction, not %

    return dictA

#Calculate the tuning parameters for the idealized single-layer radiative transfer model
#for the individual time period (i.e. control or warmed)
#See Figure 1 of Taylor et al., 2007, and other parts of that paper. Equations referenced are from there.
#Inputs: the dictionary output by loadNetCDF
#Outputs: a dictionary of additional outputs
def parameters(dictA):
    #Clear-sky parameters
    a_clr = dictA['rsuscs']/dictA['rsdscs'] #Surface albedo
    Q = dictA['rsdscs']/dictA['rsdt'] #Ratio of incident surface flux to insolation
    mu_clr = dictA['rsutcs']/dictA['rsdt']+Q*(1.-a_clr) #Atmospheric transmittance (Eq. 9)  #"Invalid value in divide"
    ga_clr = (mu_clr-Q)/(mu_clr-a_clr*Q) #Atmospheric scattering coefficient (Eq. 10)

    #Overcast parameters
    a_oc = dictA['rsusoc']/dictA['rsdsoc'] #Surface albedo
    Q = dictA['rsdsoc']/dictA['rsdt'] #Ratio of incident surface flux to insolation
    mu_oc = dictA['rsutoc']/dictA['rsdt']+Q*(1.-a_oc) #Atmospheric transmittance (Eq. 9)
    ga_oc = (mu_oc-Q)/(mu_oc-a_oc*Q) #Atmospheric scattering coefficient (Eq. 10)

    #Calculating cloudy parameters based on clear-sky and overcast ones
    #Difference between _cld and _oc: _cld is due to the cloud itself, as opposed to
    #scattering and absorption from all constituents including clouds in overcast skies.
    mu_cld = mu_oc / mu_clr            #Eq. 14
    ga_cld = (ga_oc-1.)/(1.-ga_clr)+1. #Eq. 13

    #Save the relevant variables to a dictionary for later use
    dictB = dict()
    dictB['a_clr'] = a_clr
    dictB['a_oc'] = a_oc
    dictB['mu_clr'] = mu_clr
    dictB['mu_cld'] = mu_cld
    dictB['ga_clr'] = ga_clr
    dictB['ga_cld'] = ga_cld

    return dictB


#Calculations for the differences between time periods
def d_albedo(dict1A, dict1B, dict2A, dict2B):

    #First, Ting set cloud values that were masked in one time period
    #equal to the value in the other time period, assuming no cloud changes.
    #I'll take these variables out of the dictionary before modifying them.
    a_oc1 = dict1B['a_oc']
    a_oc2 = dict2B['a_oc']
    a_oc2[a_oc2.mask == True] = a_oc1[a_oc2.mask == True]
    a_oc1[a_oc1.mask == True] = a_oc2[a_oc1.mask == True]

    mu_cld1 = dict1B['mu_cld']
    mu_cld2 = dict2B['mu_cld']
    mu_cld2[mu_cld2.mask == True] = mu_cld1[mu_cld2.mask == True]
    mu_cld1[mu_cld1.mask == True] = mu_cld2[mu_cld1.mask == True]

    ga_cld1 = dict1B['ga_cld']
    ga_cld2 = dict2B['ga_cld']
    ga_cld2[ga_cld2.mask == True] = ga_cld1[ga_cld2.mask == True]
    ga_cld1[ga_cld1.mask == True] = ga_cld2[ga_cld1.mask == True]

    #Now a bunch of calls to the "albedo" function to see how the albedo changes as a result of
    #...the changes to each of the radiative components.

    #Retrieve other variables from dictionaries to make calls to albedo shorter/more readable
    c1 = dict1A['c']
    c2 = dict2A['c']
    a_clr1 = dict1B['a_clr']
    a_clr2 = dict2B['a_clr']
    mu_clr1 = dict1B['mu_clr']
    mu_clr2 = dict2B['mu_clr']
    ga_clr1 = dict1B['ga_clr']
    ga_clr2 = dict2B['ga_clr']

    #Base state albedo
    A1 = albedo(c1, a_clr1, a_oc1, mu_clr1, mu_cld1, ga_clr1, ga_cld1)
    A2 = albedo(c2, a_clr2, a_oc2, mu_clr2, mu_cld2, ga_clr2, ga_cld2)

    #Change in albedo due to each component (Taylor et al., 2007, Eq. 12b)
    dA_c =      .5*(albedo(c2, a_clr1, a_oc1, mu_clr1, mu_cld1, ga_clr1, ga_cld1)-A1)+.5*(
                 A2-albedo(c1, a_clr2, a_oc2, mu_clr2, mu_cld2, ga_clr2, ga_cld2))

    dA_a_clr =  .5*(albedo(c1, a_clr2, a_oc1, mu_clr1, mu_cld1, ga_clr1, ga_cld1)-A1)+.5*(
                 A2-albedo(c2, a_clr1, a_oc2, mu_clr2, mu_cld2, ga_clr2, ga_cld2))

    dA_a_oc =   .5*(albedo(c1, a_clr1, a_oc2, mu_clr1, mu_cld1, ga_clr1, ga_cld1)-A1)+.5*(
                 A2-albedo(c2, a_clr2, a_oc1, mu_clr2, mu_cld2, ga_clr2, ga_cld2))

    dA_mu_clr = .5*(albedo(c1, a_clr1, a_oc1, mu_clr2, mu_cld1, ga_clr1, ga_cld1)-A1)+.5*(
                 A2-albedo(c2, a_clr2, a_oc2, mu_clr1, mu_cld2, ga_clr2, ga_cld2))

    dA_mu_cld = .5*(albedo(c1, a_clr1, a_oc1, mu_clr1, mu_cld2, ga_clr1, ga_cld1)-A1)+.5*(
                 A2-albedo(c2, a_clr2, a_oc2, mu_clr2, mu_cld1, ga_clr2, ga_cld2))

    dA_ga_clr = .5*(albedo(c1, a_clr1, a_oc1, mu_clr1, mu_cld1, ga_clr2, ga_cld1)-A1)+.5*(
                 A2-albedo(c2, a_clr2, a_oc2, mu_clr2, mu_cld2, ga_clr1, ga_cld2))

    dA_ga_cld = .5*(albedo(c1, a_clr1, a_oc1, mu_clr1, mu_cld1, ga_clr1, ga_cld2)-A1)+.5*(
                 A2-albedo(c2, a_clr2, a_oc2, mu_clr2, mu_cld2, ga_clr2, ga_cld1))

    #Set changes due to overcast or cloudy sky parameters, or changes to clouds themselves, to zero
    #...if cloud fraction is less than 3% in either time period
    dA_a_oc  [dict1A['c'] < .03] = 0
    dA_a_oc  [dict2A['c'] < .03] = 0
    dA_mu_cld[dict1A['c'] < .03] = 0
    dA_mu_cld[dict2A['c'] < .03] = 0
    dA_ga_cld[dict1A['c'] < .03] = 0
    dA_ga_cld[dict2A['c'] < .03] = 0
    dA_c     [dict1A['c'] < .03] = 0
    dA_c     [dict2A['c'] < .03] = 0

    #Combine different components into changes due to surface albedo, atmospheric clear-sky and atmospheric cloudy-sky
    dA_a = dA_a_clr + dA_a_oc                #Eq. 16a
    dA_cld = dA_mu_cld + dA_ga_cld + dA_c    #Eq. 16b
    dA_clr = dA_mu_clr + dA_ga_clr           #Eq. 16c

    #Set all planetary albedo changes = zero when incoming solar radaition is zero
    #(This will replace NaNs with zeros in the polar night--affects annual means)
    dA_a     [dict2A['rsdt']<0.1] = 0
    dA_clr   [dict2A['rsdt']<0.1] = 0
    dA_cld   [dict2A['rsdt']<0.1] = 0
    dA_a_clr [dict2A['rsdt']<0.1] = 0
    dA_a_oc  [dict2A['rsdt']<0.1] = 0
    dA_mu_cld[dict2A['rsdt']<0.1] = 0
    dA_ga_cld[dict2A['rsdt']<0.1] = 0
    dA_c     [dict2A['rsdt']<0.1] = 0
    dA_mu_clr[dict2A['rsdt']<0.1] = 0
    dA_ga_clr[dict2A['rsdt']<0.1] = 0

    #Calculate radiative effects in W/m^2 by multiplying negative of planetary albedo changes by downward SW radation
    #(This means positive changes mean more downward SW absorbed)
    surface     = -dA_a     *dict2A['rsdt']   #Radiative effect of surface albedo changes
    cloud       = -dA_cld   *dict2A['rsdt']   #Radiative effect of cloud changes
    noncloud    = -dA_clr   *dict2A['rsdt'] #Radiative effect of non-cloud SW changes (e.g. SW absorption)

    surface_clr = -dA_a_clr *dict2A['rsdt']    #Effects of surface albedo in clear-sky conditions
    surface_oc  = -dA_a_oc  *dict2A['rsdt']      #Effects of surface albedo in overcast conditions
    cloud_c     = -dA_c     *dict2A['rsdt']            #Effects of changes in cloud fraction
    cloud_ga    = -dA_ga_cld*dict2A['rsdt']      #Effects of atmospheric scattering in cloudy conditions
    cloud_mu    = -dA_mu_cld*dict2A['rsdt']      #Effects of atmospheric absorption in cloudy conditions
    noncloud_ga = -dA_ga_clr*dict2A['rsdt']   #Effects of atmospheric scattering in clear-sky conditions
    noncloud_mu = -dA_mu_clr*dict2A['rsdt']   #Effects of atmospheric absorption in clear-sky conditions

    surface    [dict2A['rsdt']<0.1] = 0
    cloud      [dict2A['rsdt']<0.1] = 0
    noncloud   [dict2A['rsdt']<0.1] = 0
    surface_clr[dict2A['rsdt']<0.1] = 0
    surface_oc [dict2A['rsdt']<0.1] = 0
    cloud_c    [dict2A['rsdt']<0.1] = 0
    cloud_ga   [dict2A['rsdt']<0.1] = 0
    cloud_mu   [dict2A['rsdt']<0.1] = 0
    noncloud_ga[dict2A['rsdt']<0.1] = 0
    noncloud_mu[dict2A['rsdt']<0.1] = 0

#   surface = np.ma.masked_outside(surface, -100, 100) # Ting called this "boundary for strange output"
#   cloud   = np.ma.masked_outside(cloud, -100, 100) # Ting called this "boundary for strange output"

    #Calculate more useful radiation output
    CRF = dict1A['rsut'] - dict1A['rsutcs'] - dict2A['rsut'] + dict2A['rsutcs'] #Change in cloud radiative effect
    cs = dict1A['rsutcs'] - dict2A['rsutcs']  #Change in clear-sky upward SW flux at TOA
    sw_net_toa = ( dict2A['rsdt'] - dict2A['rsut'] ) -  ( dict1A['rsdt'] - dict1A['rsut'] )  #Change in net SW flux at TOA
    lw_net_toa =  dict2A['rlut'] - dict1A['rlut']  #Change in OLR at TOA
    lw_cre =  ( dict2A['rlut'] - dict2A['rlutcs'] ) - ( dict1A['rlut'] - dict1A['rlutcs'] ) #Change in LW cloud radiative effect

    #Define a dictionary to return all the variables calculated here
    dictC = dict()
    dictC['A1'] = A1
    dictC['A2'] = A2
    dictC['dA_c'] = dA_c
    dictC['dA_a_clr'] = dA_a_clr
    dictC['dA_a_oc'] = dA_a_oc
    dictC['dA_mu_clr'] = dA_mu_clr
    dictC['dA_mu_cld'] = dA_mu_cld
    dictC['dA_ga_clr'] = dA_ga_clr
    dictC['dA_ga_cld'] = dA_ga_cld
    dictC['dA_a'] = dA_a
    dictC['dA_cld'] = dA_cld
    dictC['dA_clr'] = dA_clr
    dictC['surface'] = surface
    dictC['cloud'] = cloud
    dictC['noncloud'] = noncloud
    dictC['surface_clr'] = surface_clr
    dictC['surface_oc'] = surface_oc
    dictC['cloud_c'] = cloud_c
    dictC['cloud_ga'] = cloud_ga
    dictC['cloud_mu'] = cloud_mu
    dictC['noncloud_ga'] = noncloud_ga
    dictC['noncloud_mu'] = noncloud_mu
    dictC['CRF'] = CRF
    dictC['cs'] = cs
    dictC['sw_net_toa'] = sw_net_toa
    dictC['lw_net_toa'] = lw_net_toa
    dictC['lw_cre'] = lw_cre

    return dictC

#Function to calculate the planetary albedo, A.
#Inputs: (see Fig. 1 of Taylor et al., 2007)
#   c: fraction of the region occupied by clouds
#   a_clr: clear sky surface albedo (SW flux up / SW flux down)
#   a_oc: overcast surface albedo
#   mu_clr: clear-sky transmittance of SW radiation
#   mu_cld: cloudy-sky transmittance of SW radiation
#   ga_clr: clear-sky atmospheric scattering coefficient
#   ga_cld: cloudy-sky atmospheric scattering coefficient
def albedo(c, a_clr, a_oc, mu_clr, mu_cld, ga_clr, ga_cld): #Labeled with equation numbers from Taylor et al. 2007
    mu_oc = mu_clr*mu_cld                                                            #Eq. 14
    ga_oc =  1. - (1.-ga_clr)*(1.-ga_cld)                                            #Eq. 13
    A_clr = mu_clr*ga_clr + mu_clr*a_clr*(1.-ga_clr)*(1.-ga_clr)/(1.-a_clr*ga_clr)   #Eq. 7 (clear-sky)
    A_oc = mu_oc*ga_oc + mu_oc*a_oc*(1.-ga_oc)*(1.-ga_oc)/(1.-a_oc*ga_oc)            #Eq. 7 (overcast sky)
    A = (1-c)*A_clr + c*A_oc                                                         #Eq. 15
    return A
