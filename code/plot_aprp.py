#!/usr/bin/env python
import netCDF4 as nc4
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from deepmip_dict import deepmip_dict as deepmip

matplotlib.rcParams['font.family']      = 'sans-serif'
matplotlib.rcParams['font.sans-serif']  = 'Helvetica Neue'
matplotlib.rcParams['font.size']        = 10
matplotlib.rcParams['pdf.fonttype']     = 42
matplotlib.rcParams['ps.fonttype']      = 42

varname     = 'cloud'
cld_i       = 0
#varname     = 'clear'
#varname = 'alpha'
figname     = 'deepmip.aprp_lam_' + varname + str(cld_i) + '.pdf'


for model in deepmip.keys():
    for i in range(deepmip[model]['ncase']-1):
        out_file= model + '.' + deepmip[model]['contr'][i] +'.to.' + deepmip[model]['sensi'][i] + '.aprp.nc'
        print(out_file)

        ds          = xr.open_dataset(out_file)

        coslat      = np.cos(np.deg2rad(ds.lat))
        gw          = coslat / coslat.mean(dim='lat')
        dts         = ((ds.tas2 - ds.tas1) * gw).mean(dim=('time', 'lon', 'lat'))

        deepmip[model]['gmst'][i]      = (ds.tas1   * gw).mean(dim=('time', 'lon', 'lat')) - 273.15
        deepmip[model]['clear'][i]     = (ds.clr    * gw).mean(dim=('time', 'lon', 'lat')) / dts
        deepmip[model]['alpha'][i]     = (ds.alf    * gw).mean(dim=('time', 'lon', 'lat')) / dts
        deepmip[model]['cloud'][i][0]  = (ds.cld    * gw).mean(dim=('time', 'lon', 'lat')) / dts
        deepmip[model]['cloud'][i][1]  = (ds.cld_c  * gw).mean(dim=('time', 'lon', 'lat')) / dts
        deepmip[model]['cloud'][i][2]  = (ds.cld_ga * gw).mean(dim=('time', 'lon', 'lat')) / dts
        deepmip[model]['cloud'][i][3]  = (ds.cld_mu * gw).mean(dim=('time', 'lon', 'lat')) / dts

    print('cloud', deepmip[model]['cloud'])
    print('clear', deepmip[model]['clear'])
    print('alpha', deepmip[model]['alpha'])


# make plot

fig1        = plt.figure(figsize=(3.5, 3))
ax          = plt.axes()
for model in deepmip.keys():
    ncase   = deepmip[model]['ncase'] - 1
    var     = deepmip[model][varname][:,cld_i]
   #var     = deepmip[model][varname][:]
    pcolor  = deepmip[model]['pcolor']
    gmst    = deepmip[model]['gmst']
    leg     = deepmip[model]['group']
   #plt.plot(gmst[0], var[0], marker="o", color=pcolor, markersize=6, label=leg)
    plt.plot(gmst[1:ncase], var[1:ncase], marker="o", color=pcolor, linestyle='-', linewidth=2, markersize=8, label=leg)

ax.tick_params(bottom=True, top=True, left=True, right=True, labelbottom=True, labeltop=False, labelleft=True, labelright=False)
ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.yaxis.set_major_locator(ticker.MultipleLocator(.1))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(.1))
plt.xlabel("GMST")
#plt.ylabel("$λ_{cld\_sw}$")
plt.ylabel("$λ_{α}$")
plt.xlim(10, 35)
plt.ylim(-0.5, 1.5)
#plt.ylim(-0.1, .3)
plt.legend(loc=0, fontsize=7, ncol=2, handlelength=1.5, columnspacing=0.5)
plt.tight_layout()
plt.savefig(figname, format='pdf', bbox_inches = "tight")
plt.show()

