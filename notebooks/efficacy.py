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
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

# %%
from fair.energy_balance_model import EnergyBalanceModel

# %%
os.makedirs('../plots', exist_ok=True)

# %%
pl.rcParams['font.size'] = 7
pl.rcParams['font.family'] = 'Arial'
pl.rcParams['xtick.direction'] = 'out'
pl.rcParams['xtick.minor.visible'] = True
pl.rcParams['ytick.minor.visible'] = True
pl.rcParams['ytick.right'] = True
pl.rcParams['xtick.top'] = True
pl.rcParams['figure.dpi'] = 150

# %%
df_params = pd.read_csv("../data/calibrated_constrained_parameters_1.4.3.csv", index_col=0)

# %%
params_median = df_params.median()

# %%
params_median['clim_c1']

# %%
# ebm10 = EnergyBalanceModel(
#     n_timesteps=5000,
#     stochastic_run=False,
#     deep_ocean_efficacy=1,
#     ocean_heat_transfer=[params_median['ocean_heat_transfer[0]'], params_median['ocean_heat_transfer[1]'], params_median['ocean_heat_transfer[2]']], 
#     ocean_heat_capacity=[params_median['ocean_heat_capacity[0]'], params_median['ocean_heat_capacity[1]'], params_median['ocean_heat_capacity[2]']],
# )
# ebm15 = EnergyBalanceModel(
#     n_timesteps=5000,
#     stochastic_run=False,
#     deep_ocean_efficacy=1.5,
#     ocean_heat_transfer=[params_median['ocean_heat_transfer[0]'], params_median['ocean_heat_transfer[1]'], params_median['ocean_heat_transfer[2]']], 
#     ocean_heat_capacity=[params_median['ocean_heat_capacity[0]'], params_median['ocean_heat_capacity[1]'], params_median['ocean_heat_capacity[2]']],
# )
# ebm20 = EnergyBalanceModel(
#     n_timesteps=5000,
#     stochastic_run=False,
#     deep_ocean_efficacy=2,
#     ocean_heat_transfer=[params_median['ocean_heat_transfer[0]'], params_median['ocean_heat_transfer[1]'], params_median['ocean_heat_transfer[2]']], 
#     ocean_heat_capacity=[params_median['ocean_heat_capacity[0]'], params_median['ocean_heat_capacity[1]'], params_median['ocean_heat_capacity[2]']],
# )
# ebm50 = EnergyBalanceModel(
#     n_timesteps=5000,
#     stochastic_run=False,
#     deep_ocean_efficacy=5,
#     ocean_heat_transfer=[params_median['ocean_heat_transfer[0]'], params_median['ocean_heat_transfer[1]'], params_median['ocean_heat_transfer[2]']], 
#     ocean_heat_capacity=[params_median['ocean_heat_capacity[0]'], params_median['ocean_heat_capacity[1]'], params_median['ocean_heat_capacity[2]']],
# )


ebm10 = EnergyBalanceModel(
    n_timesteps=5000,
    stochastic_run=False,
    deep_ocean_efficacy=1,
    ocean_heat_transfer=[params_median['clim_kappa1'], params_median['clim_kappa2']], 
    ocean_heat_capacity=[params_median['clim_c1'], params_median['clim_c2']],
)
ebm15 = EnergyBalanceModel(
    n_timesteps=5000,
    stochastic_run=False,
    deep_ocean_efficacy=1.5,
    ocean_heat_transfer=[params_median['clim_kappa1'], params_median['clim_kappa2']], 
    ocean_heat_capacity=[params_median['clim_c1'], params_median['clim_c2']],
)
ebm20 = EnergyBalanceModel(
    n_timesteps=5000,
    stochastic_run=False,
    deep_ocean_efficacy=2,
    ocean_heat_transfer=[params_median['clim_kappa1'], params_median['clim_kappa2']], 
    ocean_heat_capacity=[params_median['clim_c1'], params_median['clim_c2']],
)
ebm50 = EnergyBalanceModel(
    n_timesteps=5000,
    stochastic_run=False,
    deep_ocean_efficacy=5,
    ocean_heat_transfer=[params_median['clim_kappa1'], params_median['clim_kappa2']], 
    ocean_heat_capacity=[params_median['clim_c1'], params_median['clim_c2']],
)

# %%
ebm10.add_forcing(np.ones(5000)*3.93, 1)
ebm15.add_forcing(np.ones(5000)*3.93, 1)
ebm20.add_forcing(np.ones(5000)*3.93, 1)
ebm50.add_forcing(np.ones(5000)*3.93, 1)

# %%
ebm10.run()
ebm15.run()
ebm20.run()
ebm50.run()

# %%
ebm10.temperature

# %%
ebm10.toa_imbalance


# %%
def effective_climate_sensitivity(ebm):
    n = ebm.toa_imbalance
    f = ebm.forcing
    t = ebm.temperature[:,0]
    lambda_eff = (n - f) / t
    sens_eff = -1/lambda_eff
    return (sens_eff)


# %%
fig, ax = pl.subplots(1, 2, figsize=(15/2.54, 8/2.54))

ax[0].plot(ebm10.temperature[:500, 0], ebm10.toa_imbalance[:500], label=r'$\varepsilon=1.0$')
ax[0].plot(ebm15.temperature[:500, 0], ebm15.toa_imbalance[:500], label=r'$\varepsilon=1.5$')
ax[0].plot(ebm20.temperature[:500, 0], ebm20.toa_imbalance[:500], label=r'$\varepsilon=2.0$')
ax[0].plot(ebm50.temperature[:500, 0], ebm50.toa_imbalance[:500], label=r'$\varepsilon=5.0$')

# zero axis
ax[0].spines['left'].set_position('zero')
ax[0].spines['right'].set_color('none')
ax[0].yaxis.tick_left()
ax[0].spines['bottom'].set_position('zero')
ax[0].spines['top'].set_color('none')
ax[0].xaxis.tick_bottom()
ax[0].set_ylim(0, 4.2)

# arrows at the end of axes
ax[0].plot(1, 0, ">k", transform=ax[0].get_yaxis_transform(), clip_on=False)
ax[0].plot(0, 1, "^k", transform=ax[0].get_xaxis_transform(), clip_on=False)

ax[0].set_xticklabels([]);
ax[0].set_yticklabels([]);
ax[0].legend(frameon=False)

ax[0].set_xlabel('Temperature, K')
ax[0].set_ylabel('Top of atmosphere radiation imbalance, W m$^{-2}$')
ax[0].set_title('(a) Gregory plot')

ax[1].plot(effective_climate_sensitivity(ebm10)[:500]);
ax[1].plot(effective_climate_sensitivity(ebm15)[:500]);
ax[1].plot(effective_climate_sensitivity(ebm20)[:500]);
ax[1].plot(effective_climate_sensitivity(ebm50)[:500]);

# zero axis
ax[1].spines['left'].set_position('zero')
ax[1].spines['right'].set_color('none')
ax[1].yaxis.tick_left()
ax[1].spines['bottom'].set_position('zero')
ax[1].spines['top'].set_color('none')
ax[1].xaxis.tick_bottom()

ax[1].set_ylim(0, 0.82)

# arrows at the end of axes
ax[1].plot(1, 0, ">k", transform=ax[1].get_yaxis_transform(), clip_on=False)
ax[1].plot(0, 1, "^k", transform=ax[1].get_xaxis_transform(), clip_on=False)

ax[1].set_xticklabels([]);
ax[1].set_yticklabels([]);
#ax[1].legend(frameon=False)

ax[1].set_xlabel('Time')
ax[1].set_ylabel('Effective climate sensitivity, W m$^{-2}$ K$^{-1}$')

ax[1].set_title('(b) Effective climate sensitivity')

fig.tight_layout()
pl.savefig('../plots/efficacy.png')

# %%
