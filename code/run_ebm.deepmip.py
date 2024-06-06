#!/usr/bin/env python
import os
import ebm_deepmip as ebm
from deepmip_dict import deepmip_dict
import subprocess

# data_dir    = '/Volumes/WD_Elements/DeepMIP_database/ceda/deepmip-eocene-p1/'
data_dir    = '/Users/wb19586/Documents/coding_github/aprp_deepmip/data/deepmip-eocene-p1/'
out_dir = '/Users/wb19586/Documents/coding_github/aprp_deepmip/ebm_output_data/deepmip/'
out_dir_aprp = '/Users/wb19586/Documents/coding_github/aprp_deepmip/aprp_output_data/deepmip/'

months      = [0, 11]

# check if the output directory exists
if not os.path.exists(out_dir):
    # Create the directory
    os.makedirs(out_dir)

for model in deepmip_dict.keys():
    print(model)
    for i in range( (deepmip_dict[model]['ncase'])):
        print(deepmip_dict[model]['contr'][i])
        prefix1     = data_dir + deepmip_dict[model]['group'] + '/' + model + '/' + deepmip_dict[model]['contr'][i] + '/' + deepmip_dict[model]['versn'] + '/climatology/'
        prefix2     = data_dir + deepmip_dict[model]['group'] + '/' + model + '/' + deepmip_dict[model]['sensi'][i] + '/' + deepmip_dict[model]['versn'] + '/climatology/'
        posfix1     = model + '_' + deepmip_dict[model]['contr'][i] + '_' + deepmip_dict[model]['versn'] + '.mean.nc'
        posfix2     = model + '_' + deepmip_dict[model]['sensi'][i] + '_' + deepmip_dict[model]['versn'] + '.mean.nc'

        ctrl_files  = {'ts'    : prefix1 + 'ts_'   + posfix1,
                       'ps'    : prefix1 + 'ps_'   + posfix1,
                       'rsut'   : prefix1 + 'rsut_'  + posfix1,
                       'rsdt'   : prefix1 + 'rsdt_'  + posfix1,
                       'rlut'   : prefix1 + 'rlut_'  + posfix1,
                       'rlus'   : prefix1 + 'rlus_'  + posfix1,
                       'rsutcs'   : prefix1 + 'rsutcs_'  + posfix1,
                       'rlutcs'   : prefix1 + 'rlutcs_'  + posfix1,
                       'rsds'   : prefix1 + 'rsds_'  + posfix1,
                       'rsus'   : prefix1 + 'rsus_'  + posfix1}
        sens_files  = {'ts'    : prefix2 + 'ts_'   + posfix2,
                       'ps'    : prefix2 + 'ps_'   + posfix2,
                       'rsut'   : prefix2 + 'rsut_'  + posfix2,
                       'rsdt'   : prefix2 + 'rsdt_'  + posfix2,
                       'rlut'   : prefix2 + 'rlut_'  + posfix2,
                       'rlus'   : prefix2 + 'rlus_'  + posfix2,
                       'rsutcs'   : prefix2 + 'rsutcs_'  + posfix2,
                       'rlutcs'   : prefix2 + 'rlutcs_'  + posfix2,
                       'rsds'   : prefix2 + 'rsds_'  + posfix2,
                       'rsus'   : prefix2 + 'rsus_'  + posfix2}
        
        out_file    = out_dir + model + '.' + deepmip_dict[model]['contr'][i] +'.to.' + deepmip_dict[model]['sensi'][i] + '.ebm.nc'
        # pass APRP output to EBM function to calculate temperature response due to APRP components
        aprp_file    = out_dir_aprp + model + '.' + deepmip_dict[model]['contr'][i] +'.to.' + deepmip_dict[model]['sensi'][i] + '.aprp.nc'

        ebm.ebm_main(ctrl_files, months[0], months[1], sens_files, months[0], months[1], out_file, aprp_file)