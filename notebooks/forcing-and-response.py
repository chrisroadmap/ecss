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
df_params = pd.read_csv("../data/calibrated_constrained_parameters_1.4.0.csv", index_col=0)

# %%
df_forcing = pd.read_csv("../data/igcc2024/ERF_best_aggregates_1750-2024.csv")
trans_in = df_forcing["anthro"] + df_forcing["solar"] + df_forcing["volcanic"]

# %%
params_median = df_params.median()

# %%
params_median['deep_ocean_efficacy']

# %%
ebm_const = EnergyBalanceModel(
    n_timesteps=300,
    stochastic_run=True,
    deep_ocean_efficacy=params_median['deep_ocean_efficacy'],
    ocean_heat_transfer=[params_median['ocean_heat_transfer[0]'], params_median['ocean_heat_transfer[1]'], params_median['ocean_heat_transfer[2]']], 
    ocean_heat_capacity=[params_median['ocean_heat_capacity[0]'], params_median['ocean_heat_capacity[1]'], params_median['ocean_heat_capacity[2]']],
    gamma_autocorrelation=params_median['gamma_autocorrelation'],
    sigma_eta=params_median['sigma_eta'],
    sigma_xi=params_median['sigma_xi'],
    seed=115
)

ebm_trans = EnergyBalanceModel(
    n_timesteps=300,
    stochastic_run=True,
    deep_ocean_efficacy=params_median['deep_ocean_efficacy'],
    ocean_heat_transfer=[params_median['ocean_heat_transfer[0]'], params_median['ocean_heat_transfer[1]'], params_median['ocean_heat_transfer[2]']], 
    ocean_heat_capacity=[params_median['ocean_heat_capacity[0]'], params_median['ocean_heat_capacity[1]'], params_median['ocean_heat_capacity[2]']],
    gamma_autocorrelation=params_median['gamma_autocorrelation'],
    sigma_eta=params_median['sigma_eta'],
    sigma_xi=params_median['sigma_xi'],
    seed=116
)

# %%
const = np.zeros(300)
const[25:] = 2.5
ebm_const.add_forcing(const, 1)
ebm_const.run()

# %%
trans = np.zeros(300)
trans[25:] = trans_in
ebm_trans.add_forcing(trans, 1)
ebm_trans.run()

# %%
fig, ax = pl.subplots(3, 2, figsize=(15/2.54, 8/2.54))
ax[0,0].plot(np.arange(-24.5, 275), ebm_const.forcing, color='tab:blue')
ax[1,0].plot(np.arange(-24.5, 275), ebm_const.toa_imbalance, color='tab:red')
ax[2,0].plot(np.arange(-24.5, 275), ebm_const.temperature[:, 0], color='tab:purple')
ax[0,1].plot(np.arange(-24.5, 275), ebm_trans.forcing, color='tab:blue')
ax[1,1].plot(np.arange(-24.5, 275), ebm_trans.toa_imbalance, color='tab:red')
ax[2,1].plot(np.arange(-24.5, 275), ebm_trans.temperature[:, 0], color='tab:purple')

ax[0,0].set_ylim(-5, 3.2)
ax[0,1].set_ylim(-5, 3.2)
ax[1,0].set_ylim(-4, 3.2)
ax[1,1].set_ylim(-4, 3.2)
ax[2,0].set_ylim(-1.4, 2.0)
ax[2,1].set_ylim(-1.4, 2.0)

for i in range(3):
    for j in range(2):
        ax[i,j].set_xlim(-25, 275)
        ax[i,j].axhline(0, ls=':', lw=0.75, color='k')
        ax[i,j].axvline(0, ls=':', lw=0.75, color='k')
        ax[i,j].set_xticklabels([])
        ax[i,j].set_yticklabels([])
#        ax[i,j].text(0.09, 0.05, '$t=0$', ha='left', va='bottom', transform=ax[i,j].transAxes)

ax[0,0].set_ylabel('W m$^{-2}$')
ax[1,0].set_ylabel('W m$^{-2}$')
ax[2,0].set_ylabel('K')

ax[0,0].set_title('(a) Constant forcing')
ax[0,1].set_title('(b) Transient forcing')

ax[0,0].text(1/24, 0.85, 'F', ha='center', va='center', fontsize=10, transform=ax[0,0].transAxes)
ax[0,1].text(1/24, 0.85, 'F', ha='center', va='center', fontsize=10, transform=ax[0,1].transAxes)
ax[1,0].text(1/24, 0.85, 'N', ha='center', va='center', fontsize=10, transform=ax[1,0].transAxes)
ax[1,1].text(1/24, 0.85, 'N', ha='center', va='center', fontsize=10, transform=ax[1,1].transAxes)
ax[2,0].text(1/24, 0.85, 'T', ha='center', va='center', fontsize=10, transform=ax[2,0].transAxes)
ax[2,1].text(1/24, 0.85, 'T', ha='center', va='center', fontsize=10, transform=ax[2,1].transAxes)

ax[2,1].text(0.09, 0.04, '$t=0$', ha='left', va='bottom', transform=ax[2,1].transAxes)

fig.tight_layout()
pl.savefig('../plots/forcing-and-response.png')

# %%
