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
df_params = pd.read_csv("../data/fair-calibrate/calibrated_constrained_parameters_1.4.0.csv", index_col=0)

# %%
df_forcing = pd.read_csv("../data/ar6/table_A3.4x_ssp534-over_ERF_1750-2500_best_estimate.csv", index_col=0)
trans_in = df_forcing["total"]

# %%
params_median = df_params.median()

# %%
ebm_trans = EnergyBalanceModel(
    n_timesteps=751,
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
trans = trans_in.values
ebm_trans.add_forcing(trans, 1)
ebm_trans.run()

# %%
dir(ebm_trans)

# %%
fig, ax = pl.subplots(3, 1, figsize=(8/2.54, 8/2.54))
ax[0].plot(np.arange(1750, 2501), ebm_trans.temperature[:, 0], color='tab:purple')
ax[1].plot(np.arange(1750, 2501), ebm_trans.toa_imbalance, color='tab:red')
ax[2].plot(np.arange(1750, 2501), ebm_trans.ocean_heat_content_change * 0.12 / 1e24, color='tab:orange')

ax[0].set_ylim(-1, 3)
ax[1].set_ylim(-1.5, 2.2)
ax[2].set_ylim(-0.05, 0.4)

ax[0].set_xticklabels([])
ax[1].set_xticklabels([])

for i in range(3):
    ax[i].set_xlim(1900, 2300)
    ax[i].axhline(0, ls=':', lw=0.75, color='k')
#         ax[i,j].axvline(0, ls=':', lw=0.75, color='k')
#         ax[i,j].set_xticklabels([])
#         ax[i,j].set_yticklabels([])
# #        ax[i,j].text(0.09, 0.05, '$t=0$', ha='left', va='bottom', transform=ax[i,j].transAxes)

ax[0].set_ylabel('K')
ax[1].set_ylabel('W m$^{-2}$')
ax[2].set_ylabel('m')

ax[0].text(0.98, 0.05, '(a) Temperature', ha='right', va='bottom', fontsize=9, transform=ax[0].transAxes)
ax[1].text(0.98, 0.05, '(b) TOA energy imbalance', ha='right', va='bottom', fontsize=9, transform=ax[1].transAxes)
ax[2].text(0.98, 0.13, '(c) Thermosteric SLR', ha='right', va='bottom', fontsize=9, transform=ax[2].transAxes)

# ax[2,1].text(0.09, 0.04, '$t=0$', ha='left', va='bottom', transform=ax[2,1].transAxes)

fig.tight_layout()
pl.savefig('../plots/overshoot.png')

# %%
