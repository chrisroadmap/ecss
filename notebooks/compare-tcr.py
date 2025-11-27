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

# %% [markdown]
# This notebook is based on https://github.com/chrisroadmap/fair-calibrate-v1.4-comparison/blob/main/notebooks/ebm-compare.ipynb and was part of the fair-calibrate paper, [Smith et al. (2024)](https://gmd.copernicus.org/articles/17/8569/2024/).

# %%
import os

import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
import xarray as xr
from fair.energy_balance_model import multi_ebm
from fair import FAIR
from fair.interface import fill, initialise
from fair.io import read_properties
from scipy.stats import linregress

# %%
pl.rcParams['font.size'] = 7
pl.rcParams['font.family'] = 'Arial'
pl.rcParams['xtick.direction'] = 'out'
pl.rcParams['xtick.minor.visible'] = True
pl.rcParams['ytick.minor.visible'] = True
pl.rcParams['ytick.right'] = True
pl.rcParams['xtick.top'] = True
pl.rcParams['figure.dpi'] = 150

# %% [markdown]
# ## run 1pctCO2 and get actual TCR

# %%
output_ensemble_size=841

# %%
scenarios = ['1pctCO2']

# %%
df_configs_141 = pd.read_csv(
    "../data/fair-calibrate/calibrated_constrained_parameters_1.4.1.csv",
    index_col=0,
)
df_configs_141.columns

# %%
configs = list(range(output_ensemble_size))

# %%
species, properties = read_properties()

# %%
da_concentration = xr.load_dataarray(
    "../data/ar6/1pctCO2_concentration_1850-2060.nc"
)

# %%
f = FAIR()
f.define_time(1850, 1990, 1)
f.define_scenarios(scenarios)
species = ["CO2", "CH4", "N2O"]
properties = {
    "CO2": {
        "type": "co2",
        "input_mode": "concentration",
        "greenhouse_gas": True,
        "aerosol_chemistry_from_emissions": False,
        "aerosol_chemistry_from_concentration": False,
    },
    "CH4": {
        "type": "ch4",
        "input_mode": "concentration",
        "greenhouse_gas": True,
        "aerosol_chemistry_from_emissions": False,
        "aerosol_chemistry_from_concentration": False,
    },
    "N2O": {
        "type": "n2o",
        "input_mode": "concentration",
        "greenhouse_gas": True,
        "aerosol_chemistry_from_emissions": False,
        "aerosol_chemistry_from_concentration": False,
    },
}
f.define_configs(configs)
f.define_species(species, properties)
f.allocate()

# %%
da = da_concentration.loc[dict(config="unspecified", scenario="1pctCO2", timebounds=np.arange(1850, 1991))]
fe = da.expand_dims(dim=["scenario", "config"], axis=(1, 2))
f.concentration = fe.drop_vars("config") * np.ones((1, 1, output_ensemble_size, 1))

# %%
fill(
    f.climate_configs["ocean_heat_capacity"],
    df_configs_141.loc[:, "clim_c1":"clim_c3"].values
)
fill(
    f.climate_configs["ocean_heat_transfer"],
    df_configs_141.loc[:, "clim_kappa1":"clim_kappa3"].values
)
fill(f.climate_configs["deep_ocean_efficacy"], df_configs_141["clim_epsilon"].values.squeeze())
fill(f.climate_configs["gamma_autocorrelation"], df_configs_141["clim_gamma"].values.squeeze())
fill(f.climate_configs["stochastic_run"], False)
fill(f.climate_configs["forcing_4co2"], df_configs_141["clim_F_4xCO2"].values.squeeze())

# species level
f.fill_species_configs()

# carbon cycle
fill(f.species_configs["iirf_0"], df_configs_141["cc_r0"].values.squeeze(), specie="CO2")
fill(f.species_configs["iirf_airborne"], df_configs_141["cc_rA"].values.squeeze(), specie="CO2")
fill(f.species_configs["iirf_uptake"], df_configs_141["cc_rU"].values.squeeze(), specie="CO2")
fill(f.species_configs["iirf_temperature"], df_configs_141["cc_rT"].values.squeeze(), specie="CO2")

# forcing scaling
fill(f.species_configs["forcing_scale"], df_configs_141["fscale_CO2"].values.squeeze(), specie="CO2")
fill(f.species_configs["forcing_scale"], df_configs_141["fscale_CH4"].values.squeeze(), specie="CH4")
fill(f.species_configs["forcing_scale"], df_configs_141["fscale_N2O"].values.squeeze(), specie="N2O")

# initial condition of CO2 concentration (but not baseline for forcing calculations)
fill(f.species_configs["baseline_concentration"], 284.3169988, specie="CO2")
fill(f.species_configs["baseline_concentration"], 808.2490285, specie="CH4")
fill(f.species_configs["baseline_concentration"], 273.021047, specie="N2O")

fill(
    f.species_configs["forcing_reference_concentration"], 284.3169988, specie="CO2"
)
fill(
    f.species_configs["forcing_reference_concentration"], 808.2490285, specie="CH4"
)
fill(f.species_configs["forcing_reference_concentration"], 273.021047, specie="N2O")

# initial conditions
initialise(f.concentration, f.species_configs["baseline_concentration"])
initialise(f.forcing, 0)
initialise(f.temperature, 0)
initialise(f.cumulative_emissions, 0)
initialise(f.airborne_emissions, 0)

# %%
f.run()

# %%
tcr_yr70 = f.temperature[70, 0, :, 0]

# %% [markdown]
# ## theoretical TCR

# %%
ebms = multi_ebm(
    configs,
    ocean_heat_capacity=df_configs_141.loc[:, "clim_c1":"clim_c3"].values,
    ocean_heat_transfer=df_configs_141.loc[:, "clim_kappa1":"clim_kappa3"].values,
    deep_ocean_efficacy=df_configs_141["clim_epsilon"].values.squeeze(),
    gamma_autocorrelation=df_configs_141["clim_gamma"].values.squeeze(),
    forcing_4co2=df_configs_141["clim_F_4xCO2"].values.squeeze(),
    stochastic_run=[False]*841,
    sigma_eta=df_configs_141["clim_sigma_eta"].values.squeeze(),
    sigma_xi=df_configs_141["clim_sigma_xi"].values.squeeze(),
    seed=12,
    use_seed=[False]*841,
    timestep=1,
    timebounds=[0, 1]
)

# %% [markdown]
# ## Compare

# %%
r = linregress(tcr_yr70, ebms.tcr).rvalue

# %%
fig, ax = pl.subplots(1, 1, figsize=(8/2.54, 8/2.54))
ax.scatter(tcr_yr70, ebms.tcr, marker='x')
ax.plot(np.linspace(0.9, 3.5), np.linspace(0.9, 3.5), color='k')
ax.set_xlim(0.8, 3.2)
ax.set_ylim(0.8, 3.2)
ax.set_ylabel('Theoretical TCR (°C)')
ax.set_xlabel(r'1pctCO2 warming year 70 (°C)')
ax.text(0.9, 3.0, f"$r^2 = {r**2:.2f}$")
#ax.set_title('(b) Transient climate response')
fig.tight_layout()
pl.savefig('../plots/compare-tcr.png')

# %%
