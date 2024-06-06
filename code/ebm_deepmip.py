import numpy as np
import netCDF4 as nc4

#Main function to run the EBM decomposition, for an individual model (which model is specified via the "dataPaths" arguments).
#This version assumes variable names and dimensions following CMIP convention.
# The structure is identical to the APRP programme `aprp_deepmip.py`.
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
#Outputs: a netcdf file contains EBM results

#=================================================================
# Original Author: Sebastian Steinig
# Adapted for DeepMIP simulations by Sebastian Steinig (sebastian.steinig@bristol.ac.uk)
#=================================================================

def ebm_main(dataPaths1, firstMonth1, lastMonth1, dataPaths2, firstMonth2, lastMonth2, out_file, aprp_file):

    #Load files and run calculations for first time period
    dict_input1  = getZonalMeans(dataPaths1, firstMonth1, lastMonth1)
    dict_input2  = getZonalMeans(dataPaths2, firstMonth2, lastMonth2)
    dict_aprp    = getZonalMeanAPRP(aprp_file, firstMonth1, lastMonth1)

    # print(dict_aprp)
    dict_output  = ebm_partitioning(dict_input1, dict_input2, dict_aprp)

    # output results
    out_nc      = nc4.Dataset(out_file, 'w')
    out_nc.set_auto_mask(True)

    lat_in      = dict_input1['lat'] [:]
    time_in      = dict_input1['time'] [:]

    dummy       = out_nc.createDimension('lat' , len(lat_in))
    dummy       = out_nc.createDimension('time', len(time_in))

    lat         = out_nc.createVariable('lat'    , 'f4', ('lat',))
    time        = out_nc.createVariable('time'   , 'f4', ('time',))
    ts1         = out_nc.createVariable('ts1'   , 'f4', ('time','lat'))
    ts2         = out_nc.createVariable('ts2'   , 'f4', ('time','lat'))
    dt_gcm      = out_nc.createVariable('dt_gcm'   , 'f4', ('time','lat'))
    dt_ebm      = out_nc.createVariable('dt_ebm'   , 'f4', ('time','lat'))
    # emmissivity partitioning
    dt_emm      = out_nc.createVariable('dt_emm'   , 'f4', ('time','lat'))
    dt_lwcre    = out_nc.createVariable('dt_lwcre'   , 'f4', ('time','lat'))
    dt_topo     = out_nc.createVariable('dt_topo'   , 'f4', ('time','lat'))
    dt_gg       = out_nc.createVariable('dt_gg'   , 'f4', ('time','lat'))
    # planetary albedo partitioning
    dt_palb      = out_nc.createVariable('dt_palb'   , 'f4', ('time','lat'))
    dt_salb      = out_nc.createVariable('dt_salb'   , 'f4', ('time','lat'))
    dt_swcre     = out_nc.createVariable('dt_swcre'   , 'f4', ('time','lat'))
    # heat transport convergence
    dt_htc     = out_nc.createVariable('dt_htc'   , 'f4', ('time','lat'))
    # APRP
    dt_aprp_alf = out_nc.createVariable('dt_aprp_alf'   , 'f4', ('time','lat'))
    dt_aprp_cld = out_nc.createVariable('dt_aprp_cld'   , 'f4', ('time','lat'))
    dt_aprp_clr = out_nc.createVariable('dt_aprp_clr'   , 'f4', ('time','lat'))

    lat.units       = 'degrees_north'
    time.units      = 'days since 0001-01-01 00:00:00'
    lat.long_name   = 'latitude'
    time.long_name  = 'time'
    lat.axis        = 'Y'

    ts1     [:] = dict_input1['ts' ][:]
    ts2     [:] = dict_input2['ts' ][:]
    dt_gcm  [:] = dict_output['dt_gcm'][:]
    dt_ebm  [:] = dict_output['dt_ebm'][:]
    dt_emm  [:] = dict_output['dt_emm'][:]
    dt_lwcre[:] = dict_output['dt_lwcre'][:]
    # dt_topo [:] = dict_output['dt_topo'][:]
    # dt_gg   [:] = dict_output['dt_gg'][:]
    dt_palb [:] = dict_output['dt_palb'][:]
    dt_salb [:] = dict_output['dt_salb'][:]
    dt_swcre[:] = dict_output['dt_swcre'][:]
    dt_htc  [:] = dict_output['dt_htc'][:]
    dt_aprp_alf[:] = dict_output['dt_aprp_alf'][:]
    dt_aprp_cld[:] = dict_output['dt_aprp_cld'][:]
    dt_aprp_clr[:] = dict_output['dt_aprp_clr'][:]
    time    [:] = [15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5]
    lat     [:] = dict_input1['lat'][:]

    out_nc.close()

def find_longitude_axis(variable):
    """
    Determine the axis index for longitude in the given variable.

    Args:
    - variable: The NetCDF variable to examine.

    Returns:
    - The index of the longitude dimension.
    """
    for i, dim_name in enumerate(variable.dimensions):
        if dim_name in ['lon', 'longitude', 'LONGITUDE']:
            return i
    raise ValueError("Longitude dimension not found in variable dimensions")

def find_time_axis(variable):
    """
    Determine the axis index for longitude in the given variable.

    Args:
    - variable: The NetCDF variable to examine.

    Returns:
    - The index of the longitude dimension.
    """
    for i, dim_name in enumerate(variable.dimensions):
        if dim_name in ['t', 'time', 'time_counter']:
            return i
    raise ValueError("Longitude dimension not found in variable dimensions")

# read input data
def getZonalMeans(dataPaths, firstMonth, lastMonth):
    datasets = {var: nc4.Dataset(dataPaths[var]) for var in dataPaths}

    # assuming latitude and longitude dimensions are consistent across variables
    Dataset = datasets[next(iter(datasets))]  # get any dataset to extract lat/lon
    lat = lon = None
    for var in Dataset.variables:
        if var in ['lat', 'latitude']:
            lat = Dataset.variables[var][:]
        elif var in ['lon', 'longitude']:
            lon = Dataset.variables[var][:]
        elif var in ['t', 'time', 'month', 'time_counter']:
            time = Dataset.variables[var][:]

    zonal_means = {}
    for var, dataset in datasets.items():
        variable_data = dataset.variables[var][firstMonth:lastMonth+1,:,:]
        lon_axis = find_longitude_axis(dataset.variables[var])
        zonal_mean = np.mean(variable_data, axis=lon_axis)
        # print(zonal_mean[time_axis])
        zonal_mean_t = np.mean(zonal_mean, axis=0)
        zonal_means[f"{var}"] = zonal_mean_t

    zonal_means['lat'] = lat
    zonal_means['lon'] = lon
    zonal_means['time'] = time

    return zonal_means


# read APRP results
def getZonalMeanAPRP(aprp_file, firstMonth, lastMonth):
    dataset = nc4.Dataset(aprp_file)

    lat = lon = None
    for var in dataset.variables:
        if var in ['lat', 'latitude']:
            lat = dataset.variables[var][:]
        elif var in ['lon', 'longitude']:
            lon = dataset.variables[var][:]
        elif var in ['t', 'time', 'month', 'time_counter']:
            time = dataset.variables[var][:]

    var_list = ["alf", "cld", "clr"]
    zonal_means = {}
    for var in var_list:
        variable_data = dataset.variables[var]
        # Replace masked values with 0
        variable_data_filled = dataset.variables[var][firstMonth:lastMonth+1,:,:].filled(0)
        lon_axis = find_longitude_axis(variable_data)
        zonal_mean = np.mean(variable_data_filled, axis=lon_axis)
        zonal_means[f"{var}"] = zonal_mean

    zonal_means['lat'] = lat
    zonal_means['lon'] = lon
    zonal_means['time'] = time

    return zonal_means

def ebm_partitioning(input_ctrl, input_sens, dict_aprp):

    # EBM for control experiment
    # all-sky fluxes
    palb_ctrl       = np.where(input_ctrl['rsdt'] != 0, input_ctrl['rsut'] / input_ctrl['rsdt'], 0)  # Avoid division by zero
    emm_ctrl        = input_ctrl['rlut'] / input_ctrl['rlus'] # Emissivity
    heatt_ctrl      = - ( input_ctrl['rsdt'] - input_ctrl['rsut'] - input_ctrl['rlut'] ) # Heat transport convergence
    insolation_ctrl = input_ctrl['rsdt']
    tsurf_ctrl      = calc_ebm_surface_temp(palb_ctrl, emm_ctrl, heatt_ctrl, insolation_ctrl)
    # clear-sky fluxes
    palb_ctrl_cs    = np.where(input_ctrl['rsdt'] != 0, input_ctrl['rsutcs'] / input_ctrl['rsdt'], 0) 
    emm_ctrl_cs     = input_ctrl['rlutcs'] / input_ctrl['rlus']
    heatt_ctrl_cs   = - ( input_ctrl['rsdt'] - input_ctrl['rsutcs'] - input_ctrl['rlutcs'] )


    # EBM for sensitivity experiment
    # all-sky fluxes
    palb_sens       = np.where(input_sens['rsdt'] != 0, input_sens['rsut'] / input_sens['rsdt'], 0)  # Avoid division by zero
    emm_sens        = input_sens['rlut'] / input_sens['rlus'] # Emissivity
    heatt_sens      = - ( input_sens['rsdt'] - input_sens['rsut'] - input_sens['rlut'] ) # Heat transport convergence
    insolation_sens = input_sens['rsdt']
    tsurf_sens  = calc_ebm_surface_temp(palb_sens, emm_sens, heatt_sens, insolation_sens)
    # clear-sky fluxes
    palb_sens_cs       = np.where(input_sens['rsdt'] != 0, input_sens['rsutcs'] / input_sens['rsdt'], 0)
    emm_sens_cs        = input_sens['rlutcs'] / input_sens['rlus']
    heatt_sens_cs      = - ( input_sens['rsdt'] - input_sens['rsutcs'] - input_sens['rlutcs'] )

    # perturb individual parameters to get their temperature response
    dt_emm = calc_ebm_surface_temp(palb_ctrl, emm_sens, heatt_ctrl, insolation_ctrl) - calc_ebm_surface_temp(palb_ctrl, emm_ctrl, heatt_ctrl, insolation_ctrl)
    dt_emm_cs = calc_ebm_surface_temp(palb_ctrl_cs, emm_sens_cs, heatt_ctrl_cs, insolation_ctrl) - calc_ebm_surface_temp(palb_ctrl_cs, emm_ctrl_cs, heatt_ctrl_cs, insolation_ctrl)
    dt_lwcre = dt_emm - dt_emm_cs

    # error here: 
    # dt_lwcre = calc_ebm_surface_temp(palb_ctrl_cs, emm_sens_cs, heatt_sens_cs, insolation_ctrl) 
    # dt_lwcre = input_sens['rlutcs'] / input_sens['rlus']
    # dt_lwcre = calc_ebm_surface_temp(palb_ctrl_cs, emm_sens_cs, heatt_sens_cs, insolation_ctrl) 
    # dt_lwcre = dt_emm_cs


    dt_palb = calc_ebm_surface_temp(palb_sens, emm_ctrl, heatt_ctrl, insolation_ctrl) - calc_ebm_surface_temp(palb_ctrl, emm_ctrl, heatt_ctrl, insolation_ctrl)
    dt_palb_cs = calc_ebm_surface_temp(palb_sens_cs, emm_ctrl_cs, heatt_ctrl_cs, insolation_ctrl) - calc_ebm_surface_temp(palb_ctrl_cs, emm_ctrl_cs, heatt_ctrl_cs, insolation_ctrl)
 
    # Lunt et al. (2021): surface albedo term is calculated via palb + (salb_sens - salb_ctrl)
    salb_ctrl = np.where(input_ctrl['rsds'] >= 0.1, input_ctrl['rsus'] / input_ctrl['rsds'], 0)  # Avoid division by zero 
    salb_sens = np.where(input_sens['rsds'] >= 0.1, input_sens['rsus'] / input_sens['rsds'], 0)  # Avoid division by zero 
    dt_salb = calc_ebm_surface_temp(palb_ctrl + ( salb_sens - salb_ctrl ), emm_ctrl, heatt_ctrl, insolation_ctrl) - calc_ebm_surface_temp(palb_ctrl, emm_ctrl, heatt_ctrl, insolation_ctrl)

    # other possibility: surface albedo term is just change in clear-sky albedo
    # dt_salb = dt_palb_cs

    dt_swcre = dt_palb - dt_palb_cs

# 
    dt_htc = calc_ebm_surface_temp(palb_ctrl, emm_ctrl, heatt_sens, insolation_ctrl) - calc_ebm_surface_temp(palb_ctrl, emm_ctrl, heatt_ctrl, insolation_ctrl)

    # calculate temperature response due to planetary albedo changes from APRP
    # get plaenatry APRP albedo changes by dividing radiative effects (in W/m^2) by insolation (in W/m^2)
    # basically reversing lines 282-292 in aprp_deepmip.py
    palb_aprp_alf = np.where(input_sens['rsdt'] != 0, dict_aprp['alf'] / input_sens['rsdt'] * ( -1.0 ), 0)  # Avoid division by zero
    palb_aprp_cld = np.where(input_sens['rsdt'] != 0, dict_aprp['cld'] / input_sens['rsdt'] * ( -1.0 ), 0)  # Avoid division by zero
    palb_aprp_clr = np.where(input_sens['rsdt'] != 0, dict_aprp['clr'] / input_sens['rsdt'] * ( -1.0 ), 0)  # Avoid division by zero
    # calculate temperature response due to planetary albedo changes from APRP
    dt_aprp_alf = calc_ebm_surface_temp(palb_ctrl + palb_aprp_alf, emm_ctrl, heatt_ctrl, insolation_ctrl) - calc_ebm_surface_temp(palb_ctrl, emm_ctrl, heatt_ctrl, insolation_ctrl)
    dt_aprp_cld = calc_ebm_surface_temp(palb_ctrl + palb_aprp_cld, emm_ctrl, heatt_ctrl, insolation_ctrl) - calc_ebm_surface_temp(palb_ctrl, emm_ctrl, heatt_ctrl, insolation_ctrl)
    dt_aprp_clr = calc_ebm_surface_temp(palb_ctrl + palb_aprp_clr, emm_ctrl, heatt_ctrl, insolation_ctrl) - calc_ebm_surface_temp(palb_ctrl, emm_ctrl, heatt_ctrl, insolation_ctrl)

    #Save the relevant variables to a dictionary for writing to netCDF
    dict_out = dict()
    dict_out['dt_gcm'] = input_sens['ts'] - input_ctrl['ts']
    dict_out['dt_ebm'] = tsurf_sens - tsurf_ctrl
    dict_out['dt_emm'] = dt_emm
    dict_out['dt_lwcre'] = dt_lwcre
    # dict_out['dt_topo'] = dt_topo
    # dict_out['dt_gg'] = dt_gg
    dict_out['dt_palb'] = dt_palb
    # dict_out['dt_salb'] = dt_palb_cs
    dict_out['dt_salb'] = dt_salb

    dict_out['dt_swcre'] = dt_swcre
    dict_out['dt_htc'] = dt_htc
    dict_out['dt_aprp_alf'] = dt_aprp_alf
    dict_out['dt_aprp_cld'] = dt_aprp_cld
    dict_out['dt_aprp_clr'] = dt_aprp_clr

    return dict_out

def calc_ebm_surface_temp(palb, emm, heat, insolation):
    sigma       = 5.670367e-8 # Stefan-Boltzmann constant
    tsurf       = ( 1. / emm / sigma * ( insolation * ( 1. - palb ) + heat ) ) ** 0.25
    return tsurf

