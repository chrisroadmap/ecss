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
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

# %%
from fair.energy_balance_model import EnergyBalanceModel

# %%
df_params = pd.read_csv("../data/calibrated_constrained_parameters_1.4.0.csv", index_col=0)

# %%
params_median = df_params.median()

# %%
params_median['ocean_heat_capacity[0]']

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
    ocean_heat_transfer=[params_median['ocean_heat_transfer[0]'], params_median['ocean_heat_transfer[2]']], 
    ocean_heat_capacity=[params_median['ocean_heat_capacity[0]'], params_median['ocean_heat_capacity[2]']],
)
ebm15 = EnergyBalanceModel(
    n_timesteps=5000,
    stochastic_run=False,
    deep_ocean_efficacy=1.5,
    ocean_heat_transfer=[params_median['ocean_heat_transfer[0]'], params_median['ocean_heat_transfer[2]']], 
    ocean_heat_capacity=[params_median['ocean_heat_capacity[0]'], params_median['ocean_heat_capacity[2]']],
)
ebm20 = EnergyBalanceModel(
    n_timesteps=5000,
    stochastic_run=False,
    deep_ocean_efficacy=2,
    ocean_heat_transfer=[params_median['ocean_heat_transfer[0]'], params_median['ocean_heat_transfer[2]']], 
    ocean_heat_capacity=[params_median['ocean_heat_capacity[0]'], params_median['ocean_heat_capacity[2]']],
)
ebm50 = EnergyBalanceModel(
    n_timesteps=5000,
    stochastic_run=False,
    deep_ocean_efficacy=5,
    ocean_heat_transfer=[params_median['ocean_heat_transfer[0]'], params_median['ocean_heat_transfer[2]']], 
    ocean_heat_capacity=[params_median['ocean_heat_capacity[0]'], params_median['ocean_heat_capacity[2]']],
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
fig, ax = pl.subplots(figsize=(15/2.54, 15/2.54))

ax.plot(ebm10.temperature[:, 0], ebm10.toa_imbalance, label=r'$\varepsilon=1.0$')
ax.plot(ebm15.temperature[:, 0], ebm15.toa_imbalance, label=r'$\varepsilon=1.5$')
ax.plot(ebm20.temperature[:, 0], ebm20.toa_imbalance, label=r'$\varepsilon=2.0$')
ax.plot(ebm50.temperature[:, 0], ebm50.toa_imbalance, label=r'$\varepsilon=5.0$')

# zero axis
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.yaxis.tick_left()
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.tick_bottom()

# arrows at the end of axes
ax.plot(1, 0, ">k", transform=ax.get_yaxis_transform(), clip_on=False)
ax.plot(0, 1, "^k", transform=ax.get_xaxis_transform(), clip_on=False)

ax.set_xticklabels([]);
ax.set_yticklabels([]);
ax.legend(frameon=False)

ax.set_xlabel('Temperature')
ax.set_ylabel('Top of atmosphere radiation imbalance')

# %%
