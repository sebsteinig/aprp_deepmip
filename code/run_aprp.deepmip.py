#!/usr/bin/env python
import os
import aprp_deepmip as aprp
from deepmip_dict import deepmip_dict

data_dir    = '/Users/wb19586/Documents/data/deepmip_database/deepmip-eocene-p1/'
out_dir = '/Users/wb19586/Documents/coding_github/aprp_deepmip/aprp_output_data/deepmip/'

months      = [0, 11]

# check if the output directory exists
if not os.path.exists(out_dir):
    # Create the directory
    os.makedirs(out_dir)

for model in deepmip_dict.keys():
    print(model)
    for i in range( (deepmip_dict[model]['ncase'] + deepmip_dict[model]['nsens'])-1):
        print(deepmip_dict[model]['contr'][i])
        prefix1     = data_dir + deepmip_dict[model]['group'] + '/' + model + '/' + deepmip_dict[model]['contr'][i] + '/' + deepmip_dict[model]['versn'] + '/climatology/'
        prefix2     = data_dir + deepmip_dict[model]['group'] + '/' + model + '/' + deepmip_dict[model]['sensi'][i] + '/' + deepmip_dict[model]['versn'] + '/climatology/'
        posfix1     = model + '_' + deepmip_dict[model]['contr'][i] + '_' + deepmip_dict[model]['versn'] + '.mean.nc'
        posfix2     = model + '_' + deepmip_dict[model]['sensi'][i] + '_' + deepmip_dict[model]['versn'] + '.mean.nc'

        ctrl_files  = {'clt'    : prefix1 + 'clt_'   + posfix1,
                       'tas'    : prefix1 + 'tas_'   + posfix1,
                       'rsds'   : prefix1 + 'rsds_'  + posfix1,
                       'rsdscs' : prefix1 + 'rsdscs_'+ posfix1,
                       'rsdt'   : prefix1 + 'rsdt_'  + posfix1,
                       'rsus'   : prefix1 + 'rsus_'  + posfix1,
                       'rsuscs' : prefix1 + 'rsuscs_'+ posfix1,
                       'rsut'   : prefix1 + 'rsut_'  + posfix1,
                       'rsutcs' : prefix1 + 'rsutcs_'+ posfix1}
        sens_files  = {'clt'    : prefix2 + 'clt_'   + posfix2,
                       'tas'    : prefix2 + 'tas_'   + posfix2,
                       'rsds'   : prefix2 + 'rsds_'  + posfix2,
                       'rsdscs' : prefix2 + 'rsdscs_'+ posfix2,
                       'rsdt'   : prefix2 + 'rsdt_'  + posfix2,
                       'rsus'   : prefix2 + 'rsus_'  + posfix2,
                       'rsuscs' : prefix2 + 'rsuscs_'+ posfix2,
                       'rsut'   : prefix2 + 'rsut_'  + posfix2,
                       'rsutcs' : prefix2 + 'rsutcs_'+ posfix2}
        out_file    = out_dir + model + '.' + deepmip_dict[model]['contr'][i] +'.to.' + deepmip_dict[model]['sensi'][i] + '.aprp.nc'
        aprp.aprp_main(ctrl_files, months[0], months[1], sens_files, months[0], months[1], out_file)
