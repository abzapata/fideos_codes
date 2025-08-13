import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from astropy.table import Table
from astropy.time import Time
from datetime import datetime
import astropy.io; from astropy.io import fits as pyfits

plt.rcParams['font.size'] = 8 # Fonts size for plots
data_path = '../all_data/'

#Datos de rango de fechas para analizar
date_i = "2018-01-01" #"2020-10-01"#input('Write the Initial date to plot in the format "yyyy-mm-dd": ')
date_f = "2018-12-31"#input('Write the Final date to plot in the format "yyyy-mm-dd": ')
#Convierte fecha inicial de busqueda en JD
date_i = date_i+'T'+"05:05:00.0"
time_i = Time(date_i,format="isot",scale="utc")
JD_i = time_i.jd

#Convierte fecha final de búsqueda en JD
date_f = date_f+'T'+"00:00:00.0"
time_f = Time(date_f,format="isot",scale="utc")
JD_f = time_f.jd

print(JD_i, JD_f)
#################
min_seeing = 1.5
max_seeing = 2.2 # size of the fiber in the sky (pinhole size)

low_snr = 50
high_snr = 100
#############################################################################

#SNR_min, SNR_max = input('Write the ranges of SNR in the format "SNR_min, SNR_max": ') 
#Airmass_min, Airmass_max = input('Write the ranges of Airmass in the format "Airmass_min, \
#Airmass_max": ') 
#OBJECT = input('Write the RV Standard to plot ("HD10700, HD32147, HD72673, HD157347"): ') 

#Carga de datos de temperatura de fideos.
f_temp = np.loadtxt(data_path + "fideos_temperatures.dat", skiprows=1, delimiter=" ", dtype="str")
T_date = f_temp[:,0]
T_date = T_date.astype(str)#-2459000
T_time = f_temp[:,1]
T_time = T_time.astype(float)#-2459000
T1 = f_temp[:,2]
T1 = T1.astype(float)
T2 = f_temp[:,3]
T2 = T2.astype(float)
T3 = f_temp[:,4]
T3 = T3.astype(float)
T4 = f_temp[:,5]
T4 = T4.astype(float)
#
date_temp=list()
T1_fideos=list()
T2_fideos=list()
T3_fideos=list()
T4_fideos=list()
n=0
while n < len(T_time):
	if (T_time[n]>=JD_i) and (T_time[n]<=JD_f): #identifica rango de fechas
		date_temp.append(datetime.fromisoformat(T_date[n]))
		T1_fideos.append(T1[n])
		T2_fideos.append(T2[n])
		T3_fideos.append(T3[n])
		T4_fideos.append(T4[n])
	n=n+3

#Carga de datos de las RVs standard
f_rvs = np.loadtxt(data_path + "RVs_AZ_ALT_all_dimm_meteo.dat", skiprows=1, delimiter=" ", dtype="str")
date_RV = f_rvs[:,0]
date_RV = date_RV.astype(str)
RV_time = f_rvs[:,1] # tiempo JD de la observacion
RV_time = RV_time.astype(float)
BJD = f_rvs[:,2]
BJD = BJD.astype(float)#-2459000
OBJ = f_rvs[:,3] # Nombre de la RV standard
OBJ = OBJ.astype(str)
RV = f_rvs[:,4]
RV = RV.astype(float)*1000
eRV = f_rvs[:,5]
eRV = eRV.astype(float)
RV_ALT = f_rvs[:,6]
RV_ALT = RV_ALT.astype(float)
RV_AZ = f_rvs[:,7]
RV_AZ = RV_AZ.astype(float)
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
SEEING = f_rvs[:,17]
SEEING = SEEING.astype(float)
PRESS = f_rvs[:,18]
PRESS = PRESS.astype(float)
ENV_TEMP = f_rvs[:,19]
ENV_TEMP = ENV_TEMP.astype(float)
HUM = f_rvs[:,20]
HUM = HUM.astype(float)
WIND_SPEED = f_rvs[:,21]
WIND_SPEED = WIND_SPEED.astype(float)
WIND_DIR = f_rvs[:,22]
WIND_DIR = WIND_DIR.astype(float)
#############################################################################
#Nombre de RV Standards usadas en FIDEOS
RV_star = [["HD10700","HIP8102","HR509"],["HD32147","HIP23311","HR1614"], \
    ["HD72673","HIP41926","HR3384"],["HD157347","HIP85042","HR6465"]]
###############################################
def name_rv(rv_star):
    for i in range(4):
        for j in range(3):
#            n_rv=list()
            if (rv_star == RV_star[i][j]):
                name = RV_star[i][0]
    return name
###############################################

#save1 = Table([[], [], [], [], [], [], [], [], [], [], \
#	[], [], [], [], [], []],\
#	names=("#JD_OBS", "BJD", "OBJ", "RV", "eRV", "ALT", "AZ", \
#		"TEXP", "AIRMASS", "DRIFT_CO", "DRIFT_CO_E", "BS", "BS_E", "SNR", "SNR_R", "BAR_COR"), \
#	dtype=(float, float, 'S8', float, float, float, float, float, float, float,\
#		float, float, float, float, float, float))
######
# Se crean listas para los datos de cada RV standard star
######
obj_HD10700 = []
rv_HD10700=list()
erv_HD10700=list()
sig_HD10700=list()
date_HD10700=list()
obj_HD32147 = []
rv_HD32147=list()
erv_HD32147=list()
sig_HD32147=list()
date_HD32147=list()
obj_HD72673 = []
rv_HD72673=list()
erv_HD72673=list() 
sig_HD72673=list()
date_HD72673=list()
obj_HD157347 = []
rv_HD157347=list()
erv_HD157347=list()
sig_HD157347=list()
date_HD157347=list()
HD10700_SNR=list()
HD32147_SNR=list()
HD72673_SNR=list()
HD157347_SNR=list()
obj_alt=list()
obj_az=list()
obj_airm=list()
drift_fco=list()
bisec=list()
std_date_HD10700=list()
std_date_HD32147=list()
std_date_HD72673=list()
std_date_HD157347=list()
std_texp=list()
iso_date_HD10700 = list()
iso_date_HD32147 = list()
iso_date_HD72673 = list()
iso_date_HD157347 = list()
for i in range(len(RV)): # RV data to plot
	if SNR[i] > 60:
		if SEEING[i] < 1.8:
			if (BJD[i]>=JD_i) and (BJD[i]<=JD_f): #identifica rango de fechas
				if (SNR[i] >= low_snr): # Filtra por SNR mayor a low_snr
					for j in range(len(RV_star[0])):
						if name_rv(OBJ[i]) == RV_star[0][j]: # Separa RVs stars por nombre
							obj_HD10700 = RV_star[0][0]
							rv_HD10700.append(RV[i])
							date_HD10700.append(BJD[i])
							HD10700_SNR.append(SNR[i])
							obj_alt.append(RV_ALT[i])
							obj_az.append(RV_AZ[i])
							obj_airm.append(AIRMASS[i])
							drift_fco.append(DRIFT_CO[i])
							bisec.append(BS[i])
							std_date_HD10700.append(date_RV[i])
							std_texp.append(TEXP[i])
							iso = datetime.fromisoformat(date_RV[i])
							iso = iso.isoformat(timespec='seconds')
							iso_date_HD10700.append(datetime.fromisoformat(iso))
						elif name_rv(OBJ[i]) == RV_star[1][j]:
							obj_HD32147 = RV_star[1][0]
							rv_HD32147.append(RV[i])
							date_HD32147.append(BJD[i])
							HD32147_SNR.append(SNR[i])
							obj_alt.append(RV_ALT[i])
							obj_az.append(RV_AZ[i])
							obj_airm.append(AIRMASS[i])
							drift_fco.append(DRIFT_CO[i])
							bisec.append(BS[i])
							std_date_HD32147.append(date_RV[i])
							std_texp.append(TEXP[i])
							iso = datetime.fromisoformat(date_RV[i])
							iso = iso.isoformat(timespec='seconds')
							iso_date_HD32147.append(datetime.fromisoformat(iso))
						elif name_rv(OBJ[i]) == RV_star[2][j]:
							obj_HD72673 = RV_star[2][0]
							rv_HD72673.append(RV[i])
							date_HD72673.append(BJD[i])
							HD72673_SNR.append(SNR[i])
							obj_alt.append(RV_ALT[i])
							obj_az.append(RV_AZ[i])
							obj_airm.append(AIRMASS[i])
							drift_fco.append(DRIFT_CO[i])
							bisec.append(BS[i])
							std_date_HD72673.append(date_RV[i])
							std_texp.append(TEXP[i])
							iso = datetime.fromisoformat(date_RV[i])
							iso = iso.isoformat(timespec='seconds')
							iso_date_HD72673.append(datetime.fromisoformat(iso))
						elif name_rv(OBJ[i]) == RV_star[3][j]:
							obj_HD157347 = RV_star[3][0]
							rv_HD157347.append(RV[i])
							date_HD157347.append(BJD[i])
							HD157347_SNR.append(SNR[i])
							obj_alt.append(RV_ALT[i])
							obj_az.append(RV_AZ[i])
							obj_airm.append(AIRMASS[i])
							drift_fco.append(DRIFT_CO[i])
							bisec.append(BS[i])
							std_date_HD157347.append(date_RV[i])
							std_texp.append(TEXP[i])
							iso = datetime.fromisoformat(date_RV[i])
							iso = iso.isoformat(timespec='seconds')
							iso_date_HD157347.append(datetime.fromisoformat(iso))
					rv_aver_HD10700=np.mean(rv_HD10700)
					sig_HD10700=np.std(rv_HD10700)
					rv_aver_HD32147=np.mean(rv_HD32147)
					sig_HD32147=np.std(rv_HD32147)
					rv_aver_HD72673=np.mean(rv_HD72673)
					sig_HD72673=np.std(rv_HD72673)
					rv_aver_HD157347=np.mean(rv_HD157347)
					sig_HD157347=np.std(rv_HD157347)

##
## Debo crear nueva lista para cada espectro con los datos promedios de seeing, wind, etc.
##

print("Aver. RV HD10700: ",'{:.3f}'.format(rv_aver_HD10700)," sigma: ",'{:.3f}'.format(sig_HD10700),"\n",
	"Aver. RV HD32147: ",'{:.3f}'.format(rv_aver_HD32147)," sigma: ",'{:.3f}'.format(sig_HD32147),"\n",
	"Aver. RV HD72673: ",'{:.3f}'.format(rv_aver_HD72673)," sigma: ",'{:.3f}'.format(sig_HD72673),"\n",
	"Aver. RV HD157347: ",'{:.3f}'.format(rv_aver_HD157347)," sigma: ",'{:.3f}'.format(sig_HD157347))

for i in range(len(rv_HD10700)):
	if (rv_HD10700[i] > (rv_aver_HD10700+10*sig_HD10700)) or (rv_HD10700[i] < (rv_aver_HD10700-10*sig_HD10700)):
		print("HD10700",date_HD10700[i],rv_HD10700[i])
for i in range(len(rv_HD32147)):
	if (rv_HD32147[i] > (rv_aver_HD32147+10*sig_HD32147)) or (rv_HD32147[i] < (rv_aver_HD32147-10*sig_HD32147)):
		print("HD32147",date_HD32147[i],rv_HD32147[i])
for i in range(len(rv_HD72673)):
	if (rv_HD72673[i] > (rv_aver_HD72673+10*sig_HD72673)) or (rv_HD72673[i] < (rv_aver_HD72673-10*sig_HD72673)):
		print("HD72673",date_HD72673[i],rv_HD72673[i])
for i in range(len(rv_HD157347)):
	if (rv_HD157347[i] > (rv_aver_HD157347+10*sig_HD157347)) or (rv_HD157347[i] < (rv_aver_HD32147-10*sig_HD32147)):
		print("HD157347",date_HD157347[i],rv_HD157347[i])


####################################
################
# FILTRAR POR FECHAS SEEING Y TEMPERATURA AMBIENTAL, VIENTO, ETC.
################
####################################
###########################################################
###
### INICIO DE GRÁFICOS!!!
###
###########################################################
fig, axs = plt.subplots(3, 1, figsize=(12, 6), sharex=True)
xfmt = DateFormatter("%d/%m/%y %H:%M")
axs[0].xaxis.set_major_formatter(xfmt)
fig.autofmt_xdate()
fig.subplots_adjust(hspace=0)
line1, = axs[0].plot(iso_date_HD10700,rv_HD10700-rv_aver_HD10700,'r.',label='HD10700',markersize=3)#-rv_aver_HD10700
line2, = axs[0].plot(iso_date_HD32147,rv_HD32147-rv_aver_HD32147,'b.',label='HD32147', markersize=3)#-rv_aver_HD32147
line3, = axs[0].plot(iso_date_HD72673,rv_HD72673-rv_aver_HD72673,'g.',label='HD72673', markersize=3)#-rv_aver_HD72673
line4, = axs[0].plot(iso_date_HD157347,rv_HD157347-rv_aver_HD157347,'k.',label='HD157347', markersize=3)#-rv_aver_HD157347
axs[0].set_ylabel('RV - <RV> (km/s)')
axs[0].legend(handles=[line1, line2, line3, line4])
#ax1.set_title(date_i+' to '+date_f)
#axs.set_legend(bbox_to_anchor=(1.00, 0.8), loc='upper left', borderaxespad=0)
#ax1.set_axis([JD_i-1, JD_f+1, -0.10, 0.10])
#plt.legend((p1[0],p2[0],p3[0],p4[0]),("HD10700","HD32147","HD72673","HD157347"))
#plt.savefig('RV_graphs.png',dpi=600)
####################################
####################################
figT3, =axs[1].plot(date_temp,T3_fideos,'g',label="Service room",markersize=3)
figT2, =axs[1].plot(date_temp,T2_fideos,'b',label="Spec.",markersize=3)
figT1, =axs[1].plot(date_temp,T1_fideos,'r',label="Chiller",markersize=3)
figT4, =axs[1].plot(date_temp,T4_fideos,'m',label="CCD",markersize=3)
axs[1].set_ylim(13, 22)
axs[1].set_ylabel('Temperature ${}^{\circ}$C')
axs[1].legend(handles=[figT3, figT2, figT1, figT4])
####################################
SNR1, = axs[2].plot(iso_date_HD10700,HD10700_SNR,'r.', markersize=3)
SNR2, = axs[2].plot(iso_date_HD32147,HD32147_SNR,'r.', markersize=3)
SNR3, = axs[2].plot(iso_date_HD72673,HD72673_SNR,'r.', markersize=3)
SNR4, = axs[2].plot(iso_date_HD157347,HD157347_SNR,'r.', markersize=3)
axs[2].set_ylabel('SNR')
plt.show()
print("Fig1")
##f1.savefig('RV_temp_graphs_all_snr1.png',dpi=600)
#################################################################
#################################################################
###########################################################
'''
fig, axs = plt.subplots(3, 1, figsize=(12, 6), sharex=True)
fig.subplots_adjust(hspace=0)
#ax1.xaxis.set_major_formatter()
print("mark_1")
p1, = axs[0].plot(date_HD10700,rv_HD10700-rv_aver_HD10700,'r.',label='HD10700',markersize=3)#-rv_aver_HD10700
p2, = axs[0].plot(date_HD32147,rv_HD32147-rv_aver_HD32147,'b.',label='HD32147', markersize=3)#-rv_aver_HD32147
print("mark_2")
p3, = axs[0].plot(date_HD72673,rv_HD72673-rv_aver_HD72673,'g.',label='HD72673', markersize=3)#-rv_aver_HD72673
p4, = axs[0].plot(date_HD157347,rv_HD157347-rv_aver_HD157347,'m.',label='HD157347', markersize=3)#-rv_aver_HD157347
axs[0].set_ylabel('RV - <RV> (km/s)')
print("mark_3")
#axs[0].legend(handles=[p1, p2, p3, p4])
####################################
#axs[1]=plt.subplot(2, 2)
#axs[1].plot(dimm_date,dimm_seeing,'y.', markersize=3)
#axs[1].set_ylabel('Seeing arcsec')
####################################
axs[1].plot(d_date,d_seeing,'y.', markersize=3)
axs[1].set_ylabel('Seeing arcsec')
####################################
axs[2].plot(wind_date,wind_speed,'r.', markersize=3)
axs[2].set_ylabel('Wind speed m/s')
print("mark_4")
plt.show()
print("Fig2")
'''