# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import glob
import os

import pandas as pd
import numpy as np
# import scipy.stats as st
import matplotlib.pyplot as pl
# import warnings
# import h5py
# from tqdm import tqdm_notebook
# from scipy.interpolate import interp1d
# warnings.simplefilter('ignore')

# %%
os.makedirs('../plots')

# %%
pl.rcParams['figure.figsize'] = (16/2.54, 16/2.54)
pl.rcParams['font.size'] = 11
pl.rcParams['font.family'] = 'Arial'
pl.rcParams['xtick.direction'] = 'out'
pl.rcParams['xtick.minor.visible'] = True
pl.rcParams['ytick.minor.visible'] = True
pl.rcParams['ytick.right'] = True
pl.rcParams['xtick.top'] = True
pl.rcParams['figure.dpi'] = 96

# %%
models = []

for path in glob.glob(os.path.join(*os.path.normpath('../data/cmip6/*').split(os.sep))):
    models.append(os.path.split(path)[-1])

models

# %%
# models = [
# 'ACCESS-CM2',     'CanESM5-CanOE',   'CNRM-CM6-1',    'EC-Earth3-LR',     'GISS-E2-1-G',     'INM-CM5-0',       'NESM3',
# 'ACCESS-ESM1-5',  'CAS-ESM2-0',      'CNRM-CM6-1-HR', 'EC-Earth3-Veg',    'GISS-E2-1-G-CC',  'IPSL-CM6A-LR',    'NorCPM1',
# 'AWI-CM-1-1-MR',  'CESM2',           'CNRM-ESM2-1',   'EC-Earth3-Veg-LR', 'GISS-E2-1-H',     'MIROC6',          'NorESM1-F',
# 'AWI-ESM-1-1-LR', 'CESM2-FV2',                        'FGOALS-f3-L',      'GISS-E2-2-G',     'MIROC-ES2L',      'NorESM2-LM',
# 'BCC-CSM2-MR',    'CESM2-WACCM',     'E3SM-1-0',      'FGOALS-g3',        'HadGEM3-GC31-LL', 'MPI-ESM-1-2-HAM', 'NorESM2-MM',
# 'BCC-ESM1',       'CESM2-WACCM-FV2', 'E3SM-1-1',      'FIO-ESM-2-0',      'HadGEM3-GC31-MM', 'MPI-ESM1-2-HR',   'SAM0-UNICON',
# 'CAMS-CSM1-0',    'CIESM',           'E3SM-1-1-ECA',  'GFDL-CM4',         'IITM-ESM',        'MPI-ESM1-2-LR',   'TaiESM1',
# 'CanESM5',        'CMCC-CM2-SR5',    'EC-Earth3',     'GFDL-ESM4',        'INM-CM4-8',       'MRI-ESM2-0',      'UKESM1-0-LL'
# ]

# %%
historical = {}
accepted_models = []
nyears = {}

for model in models:
    historical[model] = {}
    path_hist_tas  = glob.glob(
        os.path.join(
            '..', 'data', 'cmip6', model, 'historical', '*', 'tas.txt'
        )
    )

#    if model=='CanESM5' or model=='GISS-E2-1-G':
#        dirhist  = [x for x in dirhist if 'r1i1p1f1' in x]
    # experiment missing? skip model
    if len(path_hist_tas)==0:
        print(model + ' not provided historical tas')
        continue
    historical[model]['tas'] = np.zeros((165))
    nens = 0
    for ens in path_hist_tas:
        print(ens)
        tas = np.loadtxt(ens)
        if tas.size >= 165:
            historical[model]['tas']  = historical[model]['tas'] + tas[:165]
            nens = nens + 1
    if nens == 0:
        continue
    historical[model]['tas'] = historical[model]['tas'] / nens
    nyears[model]  = 165
    historical[model]['1951-1980'] = np.mean(historical[model]['tas'][101:131]) - np.mean(historical[model]['tas'][0:51])
    historical[model]['1961-1990'] = np.mean(historical[model]['tas'][111:141]) - np.mean(historical[model]['tas'][0:51])
    historical[model]['1995-2014'] = np.mean(historical[model]['tas'][145:165]) - np.mean(historical[model]['tas'][0:51])
    # if we get this far, things have worked out well
    accepted_models.append(model)

# %%
len(accepted_models)
#nyears

# %%
igcc_temp = pd.read_csv('../data/igcc2024/IGCC_GMST_1850-2024.csv', index_col='year')
igcc_temp

# %%
sixtyoneninety=np.ones(len(accepted_models))*np.nan
fiftyoneeighty=np.ones(len(accepted_models))*np.nan
ninetyfivefourteen = np.ones(len(accepted_models))*np.nan
full=np.ones((165, len(accepted_models)))
for i, model in enumerate(accepted_models):
    full[:,i] = historical[model]['tas'][:165] - np.mean(historical[model]['tas'][0:51])
    pl.plot(np.arange(1850.5, 1850+nyears[model]), historical[model]['tas'] - np.mean(historical[model]['tas'][0:51]), lw=0.5)
    sixtyoneninety[i] = historical[model]['1961-1990']
    fiftyoneeighty[i] = historical[model]['1951-1980']
    ninetyfivefourteen[i] = historical[model]['1995-2014']
pl.fill_between(np.arange(1960, 2001), -0.7, 1.8, color='#f0f0f0')
pl.plot(np.arange(1850.5, 2025), igcc_temp, color='k', lw=2, label='Observations')
pl.xlim(1850, 2015)
pl.ylim(-0.7, 1.8)
pl.legend(frameon=False)
pl.ylabel('Temperature anomaly above 1850-1900, °C')
pl.savefig('../plots/historical_warming.png')

# %%
