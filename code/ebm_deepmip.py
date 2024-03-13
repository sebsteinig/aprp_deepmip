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

def ebm_main(dataPaths1, firstMonth1, lastMonth1, dataPaths2, firstMonth2, lastMonth2, out_file):

    #Load files and run calculations for first time period
    dict_input1  = getZonalMeans(dataPaths1, firstMonth1, lastMonth1)
    dict_input2  = getZonalMeans(dataPaths2, firstMonth2, lastMonth2)

    dict_output  = ebm_partitioning(dict_input1, dict_input2)

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

    lat.units       = 'degrees_north'
    time.units      = 'days since 0001-01-01 00:00:00'
    lat.long_name   = 'latitude'
    time.long_name  = 'time'
    lat.axis        = 'Y'

    ts1     [:] = dict_input1['ts' ][:]
    ts2     [:] = dict_input2['ts' ][:]
    dt_gcm  [:] = dict_output['dt_gcm'][:]
    dt_ebm  [:] = dict_output['dt_ebm'][:]
    # dt_emm  [:] = dict_output['dt_emm'][:]
    # dt_lwcre[:] = dict_output['dt_lwcre'][:]
    # dt_topo [:] = dict_output['dt_topo'][:]
    # dt_gg   [:] = dict_output['dt_gg'][:]
    # dt_palb [:] = dict_output['dt_palb'][:]
    # dt_salb [:] = dict_output['dt_salb'][:]
    # dt_swcre[:] = dict_output['dt_swcre'][:]
    # dt_htc  [:] = dict_output['dt_htc'][:]
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

# read data
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
        zonal_means[f"{var}"] = zonal_mean

    zonal_means['lat'] = lat
    zonal_means['lon'] = lon
    zonal_means['time'] = time

    return zonal_means

def ebm_partitioning(input_ctrl, input_sens):

    # EBM for control experiment
    palb_ctrl = np.where(input_ctrl['rsdt'] != 0, input_ctrl['rsut'] / input_ctrl['rsdt'], 0)  # Avoid division by zero
    # palb_ctrl       = input_ctrl['rsut'] / input_ctrl['rsdt'] # Planetary albedo
    emm_ctrl        = input_ctrl['rlut'] / input_ctrl['rlus'] # Emissivity
    heatt_ctrl      = - ( input_ctrl['rsdt'] - input_ctrl['rsut'] - input_ctrl['rlut'] ) # Heat transport convergence
    insolation_ctrl = input_ctrl['rsdt']
    tsurf_ctrl  = calc_ebm_surface_temp(palb_ctrl, emm_ctrl, heatt_ctrl, insolation_ctrl)

    # EBM for sensitivity experiment
    palb_sens       = input_sens['rsut'] / input_sens['rsdt'] # Planetary albedo
    emm_sens        = input_sens['rlut'] / input_sens['rlus'] # Emissivity
    heatt_sens      = - ( input_sens['rsdt'] - input_sens['rsut'] ) - input_sens['rlut'] # Heat transport convergence
    insolation_sens = input_sens['rsdt']
    tsurf_sens  = calc_ebm_surface_temp(palb_sens, emm_sens, heatt_sens, insolation_sens)

    #Save the relevant variables to a dictionary for writing to netCDF
    dict_out = dict()
    dict_out['dt_gcm'] = palb_ctrl
    dict_out['dt_ebm'] = tsurf_ctrl
    # dict_out['dt_emm'] = dt_emm
    # dict_out['dt_lwcre'] = dt_lwcre
    # dict_out['dt_topo'] = dt_topo
    # dict_out['dt_gg'] = dt_gg
    # dict_out['dt_palb'] = dt_palb
    # dict_out['dt_salb'] = dt_salb
    # dict_out['dt_swcre'] = dt_swcre
    # dict_out['dt_htd'] = dt_htd

    return dict_out

def calc_ebm_surface_temp(palb, emm, heat, insolation):
    sigma       = 5.670367e-8 # Stefan-Boltzmann constant
    tsurf       = ( 1. / emm / sigma * ( insolation * ( 1. - palb ) + heat ) ) ** 0.25
    # tsurf       = ( 1. / emm / sigma * ( insolation * ( 1. - palb ) + heat ) ) 
    return tsurf

