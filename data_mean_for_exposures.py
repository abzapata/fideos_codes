################################################
# Programa que promedia Seeing, wind, temp, etc.
# during the time the exposure lasts.
################################################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from astropy.table import Table
from astropy.time import Time
from datetime import datetime

plt.rcParams['font.size'] = 12 # Fonts size for plots

data_path = "../all_data/"
#Datos de rango de fechas para analizar
date_i = "2018-01-01" #input('Write the Initial date to plot in the format "yyyy-mm-dd": ')
date_f = "2023-05-20" #input('Write the Final date to plot in the format "yyyy-mm-dd": ')
#Convierte fecha inicial de busqueda en JD
date_i = date_i+'T'+"05:05:00.0"
time_i = Time(date_i,format="isot",scale="utc")
JD_i = time_i.jd

#Convierte fecha final de búsqueda en JD
date_f = date_f+'T'+"00:00:00.0"
time_f = Time(date_f,format="isot",scale="utc")
JD_f = time_f.jd

print(JD_i, JD_f)
'''
#################
min_seeing = 1.5
max_seeing = 2.2 # size of the fiber in the sky (pinhole size)

low_snr = 50
high_snr = 100
'''
###########################################################################
save = Table([[], [], [], [], [], [], [], [], [], [], \
	[], [], [], [], [], [], [], [], [], [], [], [], []],\
	names=("#Date", "JD_OBS", "BJD", "OBJ", "RV", "eRV", "ALT", "AZ", \
		"TEXP", "AIRMASS", "DRIFT_CO", "DRIFT_CO_E", "BS", "BS_E", "SNR", \
      "SNR_R", "BAR_COR", "SEEING", "ENV_PRESS", "ENV_TEMP", "REL_HUM", "WIND_speed", "WIND_dir"), \
	dtype=('S26', float, float,'S8', float, float, float, float, \
		float, float, float, float, float, float, float, float, float, float, float, float, float, float, float))
#############################################################################
#SNR_min, SNR_max = input('Write the ranges of SNR in the format "SNR_min, SNR_max": ') 
#Airmass_min, Airmass_max = input('Write the ranges of Airmass in the format "Airmass_min, \
#Airmass_max": ') 
#OBJECT = input('Write the RV Standard to plot ("HD10700, HD32147, HD72673, HD157347"): ') 
 
#Carga de datos de Seeing del DIMM La Silla.
f_dimm = np.loadtxt(data_path + "dimm_data_all.dat", skiprows=1, delimiter=" ", dtype="str")
d_date = f_dimm[:,0]
d_date = d_date.astype(str)
d_JD = f_dimm[:,1]
d_JD = d_JD.astype(float)#-2459000
d_seeing = f_dimm[:,2]
d_seeing = d_seeing.astype(float)

##########################################################
# Crea nueva lista de datos de Seeing con los datos contenidos 
# en el rango de fechas dado.
new_d_JD = list()
new_d_seeing = list()
for i in range(0,len(d_date)):
    if (d_JD[i] >= JD_i and d_JD[i] <= JD_f):
        new_d_JD.append(d_JD[i])
        new_d_seeing.append(d_seeing[i])
print("Data seeing: ",len(new_d_seeing))

#Carga de datos de weather La Silla.
f_weather = np.loadtxt(data_path + "weather_wdb_all_data.dat", skiprows=1, delimiter=" ", dtype="str")
w_date = f_weather[:,0]
w_date = w_date.astype(str)
w_jd = f_weather[:,1]
w_jd = w_jd.astype(float)#-2459000
w_press = f_weather[:,2]
w_press = w_press.astype(float)
w_amb_temp = f_weather[:,3]
w_amb_temp = w_amb_temp.astype(float)
w_dew_temp = f_weather[:,4]
w_dew_temp = w_dew_temp.astype(float)
w_rel_hum = f_weather[:,5]
w_rel_hum = w_rel_hum.astype(float)

##########################################################
# Crea nueva lista de datos de Weather con los datos contenidos 
# en el rango de fechas dado.
new_w_jd = list()
new_w_press = list()
new_amb_temp = list()
new_rel_hum = list()
for i in range(0,len(d_date)):
    if (w_jd[i] >= JD_i and w_jd[i] <= JD_f):
        new_w_jd.append(w_jd[i])
        new_w_press.append(w_press[i])
        new_amb_temp.append(w_amb_temp[i])
        new_rel_hum.append(w_rel_hum[i])
print("Data press: ",len(new_w_press))

#Carga de datos de wind de La Silla.
f_weather = np.loadtxt(data_path + "wind_data_all.dat", skiprows=1, delimiter=" ", dtype="str")
wind_date = f_weather[:,0]
wind_date = wind_date.astype(str)#-2459000
wind_JD = f_weather[:,1]
wind_JD = wind_JD.astype(float)#-2459000
wind_dir = f_weather[:,2]
wind_dir = wind_dir.astype(float)
wind_speed = f_weather[:,3]
wind_speed = wind_speed.astype(float)

##########################################################
# Crea nueva lista de datos de Wind con los datos contenidos 
# en el rango de fechas dado.
new_wind_JD = list()
new_wind_dir = list()
new_wind_speed = list()
for i in range(0,len(d_date)):
    if (wind_JD[i] >= JD_i and wind_JD[i] <= JD_f):
        new_wind_JD.append(wind_JD[i])
        new_wind_dir.append(wind_dir[i])
        new_wind_speed.append(wind_speed[i])
print("Data wind: ",len(new_wind_speed))

#Carga de datos de las RVs standard
f_rvs = np.loadtxt(data_path + "RVs_AZ_ALT_all.dat", skiprows=1, delimiter=" ", dtype="str")
date_RV = f_rvs[:,0]
date_RV = date_RV.astype(str)
RV_time = f_rvs[:,1] # tiempo JD de la observacion
RV_time = RV_time.astype(float)
BJD = f_rvs[:,2]
BJD = BJD.astype(float)#-2459000
OBJ = f_rvs[:,3] # Nombre de la RV standard
OBJ = OBJ.astype(str)
RV = f_rvs[:,4]
RV = RV.astype(float)#*1000
eRV = f_rvs[:,5]
eRV = eRV.astype(float)
OBJ_ALT = f_rvs[:,6]
OBJ_ALT = OBJ_ALT.astype(float)
OBJ_AZ = f_rvs[:,7]
OBJ_AZ = OBJ_AZ.astype(float)
TEXP = f_rvs[:,8]
TEXP = TEXP.astype(float)
AIRMASS = f_rvs[:,9]
AIRMASS = AIRMASS.astype(float)
DRIFT_CO = f_rvs[:,10]
DRIFT_CO = DRIFT_CO.astype(float)
DRIFT_CO_E = f_rvs[:,11]
DRIFT_CO_E = DRIFT_CO_E.astype(float)
BS = f_rvs[:,12]
BS= BS.astype(float)
BS_E = f_rvs[:,13]
BS_E = BS_E.astype(float)
SNR = f_rvs[:,14]
SNR = SNR.astype(float)
SNR_R = f_rvs[:,15]
SNR_R = SNR_R.astype(float)
BAR_COR = f_rvs[:,16]
BAR_COR =BAR_COR.astype(float)
print("baaahhhh")
##########################################################
new_range_date = list()
new_TEXP = list()
for i in range(0,len(RV_time)):
    if (RV_time[i] >= JD_i and RV_time[i] <= JD_f):
        new_range_date.append(RV_time[i])
        new_TEXP.append(TEXP[i])
print("texp data: ",len(new_TEXP))
###############################################
# Crear lista con las fechas de inicio y término de exposición
# y para cada uno de estos pares encontrar los valores de seeing
# promedio, wind promedio, etc.
d_start = list()
d_finish = list()
for i in range(len(new_range_date)):
    date_start = new_range_date[i]
    date_finish = date_start+new_TEXP[i]/(3600*24.) #date_start+TEXP[i]/(3600*24.)
    d_start.append(date_start)
    d_finish.append(date_finish)
print(d_start[10],d_finish[10],new_TEXP[10])
print("Lenght new_range_date:",len(new_range_date))

avg_seeing = list() # Seeing during the exposure
avg_w_amb_press = list() # Environmental pressure during the exposure
avg_w_amb_temp = list() # Environmental temperature during the exposure
avg_w_rel_hum = list() # Environmental humidity during the exposure
avg_wind_speed = list() # Wind speed during the exposure
avg_wind_dir = list() # Wind direction during the exposure
print("blerrrrr")
##
#j = 1
######################################
######################################
for i in range(0,len(new_range_date)):#(len(d_start)): #Fecha y tiempo inicial de cada imagen.
    n_seeing = list()
    n_amb_press = list()
    n_amb_hum = list()
    n_amb_temp = list()
    n_wind_speed = list()
    n_wind_dir = list()
    for j in range(len(d_JD)):
#        print(i,j)
        if (d_JD[j] >= (d_start[i])-0.00208) and (d_JD[j] <= d_finish[i]+0.00208):# +- 3 minutes
            n_seeing.append(d_seeing[j]) 
    try:
        av_seeing = sum(n_seeing)/len(n_seeing)
    except:
        av_seeing = 999
    avg_seeing.append(av_seeing)
############################################################################
    for j in range(len(w_jd)):
#       print(j)
        if (w_jd[j] >= (d_start[i])-0.00208) and (w_jd[j] <= d_finish[i]+0.00208):# +- 3 minutes
            n_amb_press.append(w_press[j])
            n_amb_temp.append(w_amb_temp[j])
            n_amb_hum.append(w_rel_hum[j])
    try:
        av_amb_press = sum(n_amb_press)/len(n_amb_press)
        av_amb_hum = sum(n_amb_hum)/len(n_amb_hum)
        av_amb_temp = sum(n_amb_temp)/len(n_amb_temp)
    except:
        av_amb_press = 999
        av_amb_hum = 999
        av_amb_temp = 999
    avg_w_amb_press.append(av_amb_press)
    avg_w_rel_hum.append(av_amb_hum)
    avg_w_amb_temp.append(av_amb_temp)
############################################################################
    for j in range(len(wind_JD)):
        if (wind_JD[j] >= (d_start[i])-0.00208) and (wind_JD[j] <= d_finish[i]+0.00208):# +- 3 minutes
            n_wind_speed.append(wind_speed[j])
            n_wind_dir.append(wind_dir[j])
    try:
        av_wind_speed = sum(n_wind_speed)/len(n_wind_speed)
        av_wind_dir = sum(n_wind_dir)/len(n_wind_dir)
    except:
        av_wind_speed = 999
        av_wind_dir = 999
    avg_wind_speed.append(av_wind_speed)
    avg_wind_dir.append(av_wind_dir)
#    print(OBJ[i],date_RV[i],TEXP[i],'{:.4f}'.format(RV[i]),
#          '{:.2f}'.format(avg_seeing[i]),'{:.2f}'.format(avg_w_amb_temp[i]),
#          '{:.0f}'.format(avg_wind_dir[i]))
    if j == i + i*10:
        print("aqui va: ", i)
##############################################################################
    row = [date_RV[i],'{:.5f}'.format(RV_time[i]),'{:.5f}'.format(BJD[i]),OBJ[i],RV[i],eRV[i],\
        '{:.5f}'.format(OBJ_ALT[i]),'{:.5f}'.format(OBJ_AZ[i]),TEXP[i],AIRMASS[i],DRIFT_CO[i],\
        DRIFT_CO_E[i],BS[i],BS_E[i],SNR[i],SNR_R[i],'{:.4f}'.format(BAR_COR[i]),'{:.2f}'.format(avg_seeing[i]),\
        '{:.1f}'.format(avg_w_amb_press[i]),'{:.2f}'.format(avg_w_amb_temp[i]),'{:.1f}'.format(avg_w_rel_hum[i]),\
        '{:.0f}'.format(avg_wind_speed[i]),'{:.0f}'.format(avg_wind_dir[i])]
    save.add_row(row)
save.write(data_path + "RVs_AZ_ALT_all_try4.dat", format='ascii', overwrite=True)