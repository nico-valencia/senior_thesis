# import pymc as pm
# import arviz as az
# import pandas as pd
# import numpy as np
# import xarray as xr
# import matplotlib.pyplot as plt


# # Load your data
# df = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/CMSA_weird/eventdata.txt',sep='|')
# # df = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/complete/australia/AU_events_measurements_apr21.txt',sep='|')
# print(len(df))
# eventdata = df
# snr_threshold = 8
# snrbefore_threshold = 2.0
# snrafter_threshold = 2.0
# sta_lta_threshold = 6
# ellipticity_threshold = 20
# mag_threshold = 7.0
# sierrmax = 0.3

# s_depth = 200 

# # ### truncate sta, meansnr before and after, ellip

# # # # actually truncate the dataset lmao 
# # eventdata = eventdata[eventdata['snr'] > snr_threshold]
# eventdata = eventdata[eventdata['snrbefore'] > snrbefore_threshold]
# eventdata = eventdata[eventdata['snrafter'] > snrafter_threshold]
# # eventdata = eventdata[eventdata['sta_lta'] > sta_lta_threshold]
# eventdata = eventdata[eventdata['evmg'] < mag_threshold]
# eventdata.loc[:, 'SI_calc_error'] = abs(eventdata['wiggleSI_calc'] - eventdata['wiggleSIlow_calc'])

# # eventdata = eventdata[abs(eventdata['wiggleSI_calc'])<2.0]

# # eventdata = eventdata[eventdata['SI_act_error'] < sierrmax]
# eventdata = eventdata[eventdata['SI_calc_error'] < sierrmax]
# # df = eventdata
# data = df[['evdp', 'wiggleSI_calc']].values
# print(len(data))
# # Standardize for modeling
# data = (data - data.mean(axis=0)) / data.std(axis=0)

# with pm.Model() as model:
#     mu = pm.Normal("mu", mu=0, sigma=1, shape=2)

#     chol, corr, stds = pm.LKJCholeskyCov(
#         "chol", n=2, eta=2, sd_dist=pm.HalfNormal.dist(1.0), compute_corr=True
#     )

#     cov = pm.Deterministic("cov", chol @ chol.T)

#     obs = pm.MvNormal("obs", mu=mu, chol=chol, observed=data)

#     trace = pm.sample(1000, tune=1000, target_accept=0.95)

# # Get the samples from the covariance matrix
# cov = trace.posterior["chol_corr"]

# # Extract scalar entries from the matrix without leftover coords
# cov_00 = xr.DataArray(cov.values[:, :, 0, 0], dims=("chain", "draw"), name="cov_00")
# cov_01 = xr.DataArray(cov.values[:, :, 0, 1], dims=("chain", "draw"), name="cov_01")
# cov_11 = xr.DataArray(cov.values[:, :, 1, 1], dims=("chain", "draw"), name="cov_11")

# # Build a clean xarray.Dataset
# cov_ds = xr.Dataset({
#     "cov_00": cov_00,
#     "cov_01": cov_01,
#     "cov_11": cov_11
# })

# # Wrap in InferenceData
# idata_cov = az.InferenceData(posterior=cov_ds)


# print(cov_00.values[:5])  # preview first few samples

# # Plot pairwise KDE
# az.plot_pair(idata_cov, var_names=["cov_00", "cov_01", "cov_11"], kind="kde")
# # az.plot_pair(idata_cov, var_names=["cov_01"], kind="kde")

# # plt.show()
# plt.savefig("./bayesian_covariance_pairplot.png")

# import pymc as pm
# import arviz as az
# import pandas as pd
# import numpy as np
# import xarray as xr
# import matplotlib.pyplot as plt

# # Load data
# df = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/complete/australia/AU_events_measurements_apr21.txt', sep='|')
# print("Initial size:", len(df))

# # Filtering thresholds
# snrbefore_threshold = 2.0
# snrafter_threshold = 2.0
# mag_threshold = 7.0
# sierrmax = 0.3

# # Apply filters
# df = df[df['snrbefore'] > snrbefore_threshold]
# df = df[df['snrafter'] > snrafter_threshold]
# df = df[df['evmg'] < mag_threshold]
# df['SI_calc_error'] = abs(df['wiggleSI_calc'] - df['wiggleSIlow_calc'])
# df = df[df['SI_calc_error'] < sierrmax]

# # Drop rows with NaNs in relevant columns
# df = df[['evdp', 'wiggleSI_calc']].dropna()
# print("Filtered size:", len(df))

# # Standardize data
# data = df.values
# data = (data - data.mean(axis=0)) / data.std(axis=0)

# with pm.Model() as model:
#     mu = pm.Normal("mu", mu=0, sigma=1, shape=2)

#     chol, corr, stds = pm.LKJCholeskyCov(
#         "chol", n=2, eta=2,
#         sd_dist=pm.HalfNormal.dist(1.0),
#         compute_corr=True
#     )

#     obs = pm.MvNormal("obs", mu=mu, chol=chol, observed=data)

#     trace = pm.sample(1000, tune=1000, target_accept=0.95)

# # Extract correlation matrix
# corr_matrix = trace.posterior["chol_corr"]

# # Extract off-diagonal (correlation between depth and SI)
# corr_01 = xr.DataArray(corr_matrix.values[:, :, 0, 1], dims=("chain", "draw"), name="corr_01")

# # Filter out NaNs/Infs just in case
# corr_01 = corr_01.where(np.isfinite(corr_01), drop=True)

# # Create InferenceData
# idata_corr = az.InferenceData(posterior=xr.Dataset({"corr_01": corr_01}))

# # Plot posterior of correlation
# az.plot_posterior(idata_corr, var_names=["corr_01"], hdi_prob=0.95)
# plt.title("Posterior Correlation: Depth vs. Splitting Intensity")
# plt.xlabel("Correlation coefficient")
# plt.tight_layout()
# plt.savefig("./posterior_corr_01.png")
# plt.close()

# # Summary stats
# summary = az.summary(idata_corr, var_names=["corr_01"], hdi_prob=0.95)
# print(summary)

import pymc as pm
import arviz as az
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Load and filter data
# df = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/complete/australia/AU_events_measurements_apr21.txt', sep='|')
df = pd.read_csv('/work/gcl3/BR/nicolasv03/senior_thesis/data/sample/TAM/TAM_events_measurements_may30.txt',sep='|')
xvar = 'ellipticity'
yvar = 'SI_calc_error'
snrbefore_threshold = 1.0
snrafter_threshold = 1.0
mag_threshold = 7.0
sierrmax = 0.5

# df = df[df['snrbefore'] > snrbefore_threshold]
# df = df[df['snrafter'] > snrafter_threshold]
df = df[df['evmg'] < mag_threshold]
df['SI_calc_error'] = abs(df['wiggleSI_calc'] - df['wiggleSIlow_calc'])
df = df[df['SI_calc_error'] < sierrmax]
# df = df[df['ellipticity']<200]

station_name = 'TAM'
df = df[df['evstation']==station_name].copy()

# Drop NaNs
df = df[[xvar, yvar]].dropna()
print("Filtered size:", len(df))

# Standardize data for modeling
data = df.values
data_std = (data - data.mean(axis=0)) / data.std(axis=0)

# Fit correlation model
with pm.Model() as model:
    mu = pm.Normal("mu", mu=0, sigma=1, shape=2)
    chol, corr, stds = pm.LKJCholeskyCov("chol", n=2, eta=2,
                                         sd_dist=pm.HalfNormal.dist(1.0),
                                         compute_corr=True)
    obs = pm.MvNormal("obs", mu=mu, chol=chol, observed=data_std)
    trace = pm.sample(1000, tune=1000, target_accept=0.95)

# Extract correlation coefficient samples
corr_matrix = trace.posterior["chol_corr"]
corr_01 = xr.DataArray(corr_matrix.values[:, :, 0, 1], dims=("chain", "draw"), name="corr_01")
corr_01 = corr_01.where(np.isfinite(corr_01), drop=True)
idata_corr = az.InferenceData(posterior=xr.Dataset({"corr_01": corr_01}))

# --- Plotting ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# (1) Posterior plot of correlation
az.plot_posterior(idata_corr, var_names=["corr_01"], hdi_prob=0.95, ax=axes[0])
axes[0].set_title(f"Posterior Correlation: {xvar} vs. {yvar}")
axes[0].set_xlabel("Correlation coefficient")

# (2) Scatter plot with regression line (raw data)
depth = df[xvar].values
si = df[yvar].values
slope, intercept, r_value, p_value, std_err = linregress(depth, si)
reg_line = slope * depth + intercept

axes[1].scatter(depth, si, alpha=0.4, label="Data")
axes[1].plot(depth, reg_line, color='red', label=f"Linear fit\n$R$ = {r_value:.2f}")
axes[1].set_xlabel(f"{xvar}")
axes[1].set_ylabel(f"{yvar}")
axes[1].set_title(f"{xvar} vs. {yvar} with Linear Fit")
axes[1].legend()

plt.tight_layout()
plt.savefig("./posterior_and_regression.png")
plt.close()