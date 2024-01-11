#!/usr/bin/env python
import numpy as np

deepmip_dict    = dict()
deepmip_dict['CESM1.2-CAM5']    = { \
        'ncase' : 5,
        'nsens' : 0,
        'versn' : 'v1.0',
        'group' : 'CESM',
        # 'contr' : ["deepmip-eocene-p1-PI", "deepmip-eocene-p1-x1", "deepmip-eocene-p1-x3", "deepmip-eocene-p1-x6"],
        # 'sensi' : ["deepmip-eocene-p1-x1", "deepmip-eocene-p1-x3", "deepmip-eocene-p1-x6", "deepmip-eocene-p1-x9"],
        'contr' : ["deepmip-eocene-p1-PI", "deepmip-eocene-p1-x1", "deepmip-eocene-p1-x3", "deepmip-eocene-p1-x6"],
        'sensi' : ["deepmip-eocene-p1-x1", "deepmip-eocene-p1-x3", "deepmip-eocene-p1-x6", "deepmip-eocene-p1-x9"],
        'pcolor': 'tab:red',
        }

deepmip_dict['GFDL-CM2.1']      = { \
        'ncase' : 6,
        'nsens' : 2,
        'versn' : 'v1.0',
        'group' : 'GFDL',
        # 'contr' : ["deepmip-eocene-p1-PI", "deepmip-eocene-p1-x1", "deepmip-eocene-p1-x2", "deepmip-eocene-p1-x3", "deepmip-eocene-p1-x4", "deepmip-eocene-p1-x1", "deepmip-eocene-p1-x3"],
        # 'sensi' : ["deepmip-eocene-p1-x1", "deepmip-eocene-p1-x2", "deepmip-eocene-p1-x3", "deepmip-eocene-p1-x4", "deepmip-eocene-p1-x6", "deepmip-eocene-p1-x3", "deepmip-eocene-p1-x6"],
        'contr' : ["deepmip-eocene-p1-PI", "deepmip-eocene-p1-x1", "deepmip-eocene-p1-x2", "deepmip-eocene-p1-x3", "deepmip-eocene-p1-x4", "deepmip-eocene-p1-x1", "deepmip-eocene-p1-x3"],
        'sensi' : ["deepmip-eocene-p1-x1", "deepmip-eocene-p1-x2", "deepmip-eocene-p1-x4", "deepmip-eocene-p1-x6", "deepmip-eocene-p1-x6", "deepmip-eocene-p1-x3", "deepmip-eocene-p1-x6"],
        'pcolor': 'tab:blue',
        }

deepmip_dict['HadCM3B-M2.1aN']  = { \
        'ncase' : 4,
        'nsens' : 1,
        'versn' : 'v1.0',
        'group' : 'HadCM3',
        'contr' : ["deepmip-eocene-p1-PI", "deepmip-eocene-p1-x1", "deepmip-eocene-p1-x2", "deepmip-eocene-p1-x1"],
        'sensi' : ["deepmip-eocene-p1-x1", "deepmip-eocene-p1-x2", "deepmip-eocene-p1-x3", "deepmip-eocene-p1-x3"],
        'pcolor': 'tab:brown',
        }

# deepmip_dict['HadCM3BL-M2.1aN']  = { \
#         'ncase' : 4,
#         'nsens' : 1,
#         'versn' : 'v1.0',
#         'group' : 'HadCM3',
#         'contr' : ["deepmip-eocene-p1-PI", "deepmip-eocene-p1-x1", "deepmip-eocene-p1-x2", "deepmip-eocene-p1-x1"],
#         'sensi' : ["deepmip-eocene-p1-x1", "deepmip-eocene-p1-x2", "deepmip-eocene-p1-x3", "deepmip-eocene-p1-x3"],
#         'pcolor': 'tab:brown',
#         }

deepmip_dict['IPSLCM5A2']    = { \
        'ncase' : 3,
        'nsens' : 0,
        'versn' : 'v1.0',
        'group' : 'IPSL',
        'contr' : ["deepmip-eocene-p1-PI", "deepmip-eocene-p1-x1.5"],
        'sensi' : ["deepmip-eocene-p1-x1.5", "deepmip-eocene-p1-x3"],
        'pcolor': 'tab:orange',
        }

deepmip_dict['MIROC4m']   = { \
        'ncase' : 4,
        'nsens' : 1,
        'versn' : 'v1.0',
        'group' : 'MIROC',
        'contr' : ["deepmip-eocene-p1-PI", "deepmip-eocene-p1-x1", "deepmip-eocene-p1-x2", "deepmip-eocene-p1-x1"],
        'sensi' : ["deepmip-eocene-p1-x1", "deepmip-eocene-p1-x2", "deepmip-eocene-p1-x3", "deepmip-eocene-p1-x3"],
        'pcolor': 'tab:green',
        }

deepmip_dict['NorESM1-F']  = { \
        'ncase' : 3,
        'nsens' : 0,
        'versn' : 'v1.0',
        'group' : 'NorESM',
        'contr' : ["deepmip-eocene-p1-PI", "deepmip-eocene-p1-x2"],
        'sensi' : ["deepmip-eocene-p1-x2", "deepmip-eocene-p1-x4"],
        'pcolor': 'tab:purple',
        }
