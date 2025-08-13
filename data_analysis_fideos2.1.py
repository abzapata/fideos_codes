import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
#from astropy.table import Table
from astropy.time import Time
from datetime import datetime
#import astropy.io; from astropy.io import fits as pyfits

plt.rcParams['font.size'] = 8  # Fonts size for plots

data_path = '../all_data/'
# Define the date range for analysis
date_i = "2018-01-01"
date_f = "2018-12-31"

# Convert initial and final dates to Julian Date (JD)
date_i = date_i + 'T' + "05:05:00.0"
time_i = Time(date_i, format="isot", scale="utc")
JD_i = time_i.jd

date_f = date_f + 'T' + "00:00:00.0"
time_f = Time(date_f, format="isot", scale="utc")
JD_f = time_f.jd

print(JD_i, JD_f)

# Define thresholds
min_seeing = 1.5
max_seeing = 2.2
low_snr = 50
high_snr = 100

# Load temperature data
f_temp = np.loadtxt(data_path + "fideos_temperatures.dat", skiprows=1, delimiter=" ", dtype="str")
T_date = f_temp[:, 0].astype(str)
T_time = f_temp[:, 1].astype(float)
T1 = f_temp[:, 2].astype(float)
T2 = f_temp[:, 3].astype(float)
T3 = f_temp[:, 4].astype(float)
T4 = f_temp[:, 5].astype(float)

# Filter temperature data by date range
date_temp, T1_fideos, T2_fideos, T3_fideos, T4_fideos = [], [], [], [], []
for n in range(0, len(T_time), 3):
    if JD_i <= T_time[n] <= JD_f:
        date_temp.append(datetime.fromisoformat(T_date[n]))
        T1_fideos.append(T1[n])
        T2_fideos.append(T2[n])
        T3_fideos.append(T3[n])
        T4_fideos.append(T4[n])

# Load RV data
f_rvs = np.loadtxt(data_path + "RVs_AZ_ALT_all_dimm_meteo.dat", skiprows=1, delimiter=" ", dtype="str")
date_RV = f_rvs[:, 0].astype(str)
RV_time = f_rvs[:, 1].astype(float)
BJD = f_rvs[:, 2].astype(float)
OBJ = f_rvs[:, 3].astype(str)
RV = f_rvs[:, 4].astype(float) * 1000
eRV = f_rvs[:, 5].astype(float)
RV_ALT = f_rvs[:, 6].astype(float)
RV_AZ = f_rvs[:, 7].astype(float)
TEXP = f_rvs[:, 8].astype(float)
AIRMASS = f_rvs[:, 9].astype(float)
DRIFT_CO = f_rvs[:, 10].astype(float)
DRIFT_CO_E = f_rvs[:, 11].astype(float)
BS = f_rvs[:, 12].astype(float)
BS_E = f_rvs[:, 13].astype(float)
SNR = f_rvs[:, 14].astype(float)
SNR_R = f_rvs[:, 15].astype(float)
BAR_COR = f_rvs[:, 16].astype(float)
SEEING = f_rvs[:, 17].astype(float)
PRESS = f_rvs[:, 18].astype(float)
ENV_TEMP = f_rvs[:, 19].astype(float)
HUM = f_rvs[:, 20].astype(float)
WIND_SPEED = f_rvs[:, 21].astype(float)
WIND_DIR = f_rvs[:, 22].astype(float)

# Define RV standards
RV_star = [["HD10700", "HIP8102", "HR509"], ["HD32147", "HIP23311", "HR1614"],
           ["HD72673", "HIP41926", "HR3384"], ["HD157347", "HIP85042", "HR6465"]]

def name_rv(rv_star):
    for i in range(4):
        for j in range(3):
            if rv_star == RV_star[i][j]:
                return RV_star[i][0]
    return None

# Initialize lists for each RV standard star
rv_data = {
    "HD10700": {"rv": [], "date": [], "snr": [], "iso_date": []},
    "HD32147": {"rv": [], "date": [], "snr": [], "iso_date": []},
    "HD72673": {"rv": [], "date": [], "snr": [], "iso_date": []},
    "HD157347": {"rv": [], "date": [], "snr": [], "iso_date": []}
}

# Filter and organize RV data
for i in range(len(RV)):
    if SNR[i] > 60 and SEEING[i] < 1.8 and JD_i <= BJD[i] <= JD_f and SNR[i] >= low_snr:
        star_name = name_rv(OBJ[i])
        if star_name:
            rv_data[star_name]["rv"].append(RV[i])
            rv_data[star_name]["date"].append(BJD[i])
            rv_data[star_name]["snr"].append(SNR[i])
            iso = datetime.fromisoformat(date_RV[i])
            rv_data[star_name]["iso_date"].append(iso.isoformat(timespec='seconds'))

# Calculate average and standard deviation for each star
for star in rv_data:
    rv_data[star]["rv_aver"] = np.mean(rv_data[star]["rv"])
    rv_data[star]["rv_std"] = np.std(rv_data[star]["rv"])

# Print average RV and sigma for each star
for star in rv_data:
    print(f"Aver. RV {star}: {rv_data[star]['rv_aver']:.3f} sigma: {rv_data[star]['rv_std']:.3f}")

# Plotting
fig, axs = plt.subplots(3, 1, figsize=(12, 6), sharex=True)
xfmt = DateFormatter("%d/%m/%y %H:%M")
axs[0].xaxis.set_major_formatter(xfmt)
fig.autofmt_xdate()
fig.subplots_adjust(hspace=0)

# Plot RV data
colors = {'HD10700': 'r', 'HD32147': 'b', 'HD72673': 'g', 'HD157347': 'k'}
for star in rv_data:
    axs[0].plot(rv_data[star]["iso_date"], np.array(rv_data[star]["rv"]) - rv_data[star]["rv_aver"],
                f'{colors[star]}.', label=star, markersize=3)
axs[0].set_ylabel('RV - <RV> (km/s)')
axs[0].legend()

# Plot temperature data
axs[1].plot(date_temp, T3_fideos, 'g', label="Service room", markersize=3)
axs[1].plot(date_temp, T2_fideos, 'b', label="Spec.", markersize=3)
axs[1].plot(date_temp, T1_fideos, 'r', label="Chiller", markersize=3)
axs[1].plot(date_temp, T4_fideos, 'm', label="CCD", markersize=3)
axs[1].set_ylim(13, 22)
axs[1].set_ylabel('Temperature ${}^{\circ}$C')
axs[1].legend()

# Plot SNR data
for star in rv_data:
    axs[2].plot(rv_data[star]["iso_date"], rv_data[star]["snr"], f'{colors[star]}.', markersize=3)
axs[2].set_ylabel('SNR')

plt.show()
