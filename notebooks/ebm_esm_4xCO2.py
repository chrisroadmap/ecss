# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import os

from fair.energy_balance_model import EnergyBalanceModel
import matplotlib.pyplot as pl
import pandas as pd
import numpy as np

# %%
pl.rcParams['font.size'] = 7
pl.rcParams['ytick.major.right'] = True
pl.rcParams['ytick.right'] = True
pl.rcParams['xtick.major.top'] = True
pl.rcParams['xtick.top'] = True
pl.rcParams['xtick.direction'] = 'in'
pl.rcParams['ytick.direction'] = 'in'
pl.rcParams['figure.dpi'] = 150

# %%
cmip6_output_df = pd.read_csv('../data/fair-calibrate/4xCO2_cmip6.csv')
cmip6_output_df.head()

# %%
# where models submit more than one simulation, we choose the one that looks the most stable
# same as in the fair-calibrate code
multi_runs = {
    "GISS-E2-1-G": "r1i1p1f1",
    "GISS-E2-1-H": "r1i1p3f1",
    "MRI-ESM2-0": "r1i1p1f1",
    "EC-Earth3": "r3i1p1f1",
    "FIO-ESM-2-0": "r1i1p1f1",
    "CanESM5": "r1i1p2f1",
    "FGOALS-f3-L": "r1i1p1f1",
    "CNRM-ESM2-1": "r1i1p1f2",
}

# %%
models = cmip6_output_df.climate_model.unique()

# %%
tas_cmip6 = {}
for model in models:
    if model in multi_runs:
        tas_cmip6[model] = cmip6_output_df.loc[
            (cmip6_output_df['climate_model']==model) & 
            (cmip6_output_df['member_id']==multi_runs[model]) & 
            (cmip6_output_df['variable']=='tas'), 
            'X1850':'X1999'
        ].values.squeeze()
    else:
        tas_cmip6[model] = cmip6_output_df.loc[
            (cmip6_output_df['climate_model']==model) & 
            (cmip6_output_df['variable']=='tas'), 
            'X1850':'X1999'
        ].values.squeeze()

# %%
ebm_fits_df = pd.read_csv('../data/fair-calibrate/4xCO2_cummins_ebm3_cmip6_1.4.1.csv')
ebm_fits_df.head()

# %%
# run 150 year 4xCO2 with the calibrations
# but run deterministic
ebm3 = {}
df_params = pd.DataFrame()

for isp, model in enumerate(models):
    if model in multi_runs:
        params = ebm_fits_df.loc[
            (ebm_fits_df['model']==model) &
            (ebm_fits_df['run']==multi_runs[model]),
            'gamma':
        ]
    else:
        params = ebm_fits_df.loc[ebm_fits_df['model']==model, 'gamma':]
    df_params = pd.concat([df_params, params])
    ebm3[model] = EnergyBalanceModel(
        ocean_heat_capacity=[params["C1"].values[0], params["C2"].values[0], params["C3"].values[0]],
        ocean_heat_transfer=[params["kappa1"].values[0], params["kappa2"].values[0], params["kappa3"].values[0]],
        deep_ocean_efficacy=params["epsilon"].values[0],
        gamma_autocorrelation=params["gamma"].values[0],
        sigma_xi=params["sigma_xi"].values[0],
        sigma_eta=params["sigma_eta"].values[0],
        forcing_4co2=params["F_4xCO2"].values[0],
        stochastic_run=False,
        seed=10000 * isp + 700  # reproducibility, but a different stochastic seed per model
    )
    ebm3[model].add_forcing(np.ones(151) * params["F_4xCO2"].values[0], timestep=1)
    ebm3[model].run()

# %%
fig, ax = pl.subplots(7, 7, figsize=(18/2.54, 18/2.54))
for isp, model in enumerate(sorted(models)):
    row = isp//7
    col = isp%7
    label1 = 'ESM' if model=='UKESM1-0-LL' else ''
    label2 = 'EBM' if model=='UKESM1-0-LL' else ''
    ax[row, col].plot(np.arange(0.5, 150), tas_cmip6[model], color='tab:red', label=label1, lw=1.0)
    ax[row, col].plot(np.arange(151), ebm3[model].temperature[:, 0], color='tab:purple', label=label2, lw=1.0)
    ax[row, col].set_xlim(0,150)
    ax[row, col].set_ylim(0, 10)
    ax[row, col].text(75, 8.5, model, fontsize=6, ha='center', va='bottom', fontweight='bold')
    ax[row, col].set_xticks([0, 50, 100, 150])
    ax[row, col].set_yticks([0, 2, 4, 6, 8, 10])
    if row<6 or col%2==1:
        ax[row, col].set_xticklabels([])
    if col>0 or row%2==1:
        ax[row, col].set_yticklabels([])
ax[6,3].set_xlabel('Year')
ax[3,0].set_ylabel('°C')
ax[6,6].legend(loc='lower right', frameon=False)
fig.tight_layout()
pl.subplots_adjust(wspace=0, hspace=0)

os.makedirs('../plots/', exist_ok=True)
pl.savefig('../plots/ebm_esm_4xCO2.png')

# %%
fig, ax = pl.subplots(3, 3, figsize=(18/2.54, 18/2.54))
ax[0,0].hist(df_params['kappa1'], bins=np.arange(0, 2.6, 0.2))
ax[0,0].set_xlim(0, 2.5)
ax[0,0].set_xlabel("W m$^{-2}$ K$^{-1}$")
ax[0,0].set_title('(a) $\kappa_1$')

ax[0,1].hist(df_params['kappa2'], bins=np.arange(0, 8.5, 0.5))
ax[0,1].set_xlim(0, 8)
ax[0,1].set_xlabel("W m$^{-2}$ K$^{-1}$")
ax[0,1].set_title('(b) $\kappa_2$')

ax[0,2].hist(df_params['kappa3'], bins=np.arange(0, 8.5, 0.5))
ax[0,2].set_xlim(0, 8)
ax[0,2].set_xlabel("W m$^{-2}$ K$^{-1}$")
ax[0,2].set_title('(c) $\kappa_3$')

ax[1,0].hist(df_params['C1'], bins=np.arange(0, 8.5, 0.5))
ax[1,0].set_xlim(0, 8)
ax[1,0].set_xlabel("W m$^{-2}$ K$^{-1}$ yr")
ax[1,0].set_title('(d) $C_1$')

ax[1,1].hist(df_params['C2'], bins=np.arange(0, 160, 10))
ax[1,1].set_xlim(0, 150)
ax[1,1].set_xlabel("W m$^{-2}$ K$^{-1}$ yr")
ax[1,1].set_title('(e) $C_2$')

ax[1,2].hist(df_params['C3'], bins=np.arange(0, 420, 20))
ax[1,2].set_xlim(0, 400)
ax[1,2].set_xlabel("W m$^{-2}$ K$^{-1}$ yr")
ax[1,2].set_title('(f) $C_3$')

ax[2,0].hist(df_params['epsilon'], bins=np.arange(0, 2.1, 0.1))
ax[2,0].set_xlim(0, 2)
ax[2,0].set_title('(g) $\epsilon$')

ax[2,1].hist(df_params['F_4xCO2'], bins=np.arange(5, 10.2, 0.2))
ax[2,1].set_xlim(5, 10)
ax[2,1].set_xlabel("W m$^{-2}$")
ax[2,1].set_title(r'(h) $F_{4\times\mathrm{CO}_2}$')

ax[2,2].hist(df_params['F_4xCO2']/df_params['kappa1']/2, bins=np.arange(0, 8.5, 0.5))
ax[2,2].set_xlim(0, 8)
ax[2,2].set_xlabel("K")
ax[2,2].set_title('(i) ECS')

for row in range(3):
    for col in range(3):
        ax[row, col].set_yticklabels([])
        ax[row, col].set_yticks([])

fig.tight_layout()
pl.savefig('../plots/ebm_calibrations.png')

# %%
