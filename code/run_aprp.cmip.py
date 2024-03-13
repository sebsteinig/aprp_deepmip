#!/usr/bin/env python
import os
import aprp_deepmip as aprp
from deepmip_dict import deepmip_dict

data_dir    = '/Users/wb19586/Documents/coding_github/aprp_deepmip/data/cmip6/'
out_dir = '/Users/wb19586/Documents/coding_github/aprp_deepmip/aprp_output_data/cmip/'

months      = [0, 11]

# check if the output directory exists
if not os.path.exists(out_dir):
    # Create the directory
    os.makedirs(out_dir)

for model in os.listdir(data_dir + "piControl/"):
    if model == 'ICON-ESM-LR' or model == '.DS_Store':
        continue
    print(model)
    prefix1     = data_dir + "piControl/" + model + "/"
    prefix2     = data_dir + "abrupt-4xCO2/" + model + "/"
    posfix      = 'climatology.nc'

    ctrl_files  = {'clt'    : prefix1 + 'clt_'   + posfix,
                   'tas'    : prefix1 + 'tas_'   + posfix,
                   'rsds'   : prefix1 + 'rsds_'  + posfix,
                   'rsdscs' : prefix1 + 'rsdscs_'+ posfix,
                   'rsdt'   : prefix1 + 'rsdt_'  + posfix,
                   'rsus'   : prefix1 + 'rsus_'  + posfix,
                   'rsuscs' : prefix1 + 'rsuscs_'+ posfix,
                   'rsut'   : prefix1 + 'rsut_'  + posfix,
                   'rsutcs' : prefix1 + 'rsutcs_'+ posfix}
    sens_files  = {'clt'    : prefix2 + 'clt_'   + posfix,
                   'tas'    : prefix2 + 'tas_'   + posfix,
                   'rsds'   : prefix2 + 'rsds_'  + posfix,
                   'rsdscs' : prefix2 + 'rsdscs_'+ posfix,
                   'rsdt'   : prefix2 + 'rsdt_'  + posfix,
                   'rsus'   : prefix2 + 'rsus_'  + posfix,
                   'rsuscs' : prefix2 + 'rsuscs_'+ posfix,
                   'rsut'   : prefix2 + 'rsut_'  + posfix,
                   'rsutcs' : prefix2 + 'rsutcs_'+ posfix}
    out_file    = out_dir + model + '.piControl.to.abrupt-4xCO2.aprp.nc'
    aprp.aprp_main(ctrl_files, months[0], months[1], sens_files, months[0], months[1], out_file)
