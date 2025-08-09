import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from astropy.table import Table
from astropy.time import Time
import astropy.io; from astropy.io import fits as pyfits

plt.rcParams['font.size'] = 8 # Fonts size for plots

data_path = '../all_data/'
#Datos de rango de fechas para analizar
date_i = "2020-03-10"#input('Write the Initial date to plot in the format "yyyy-mm-dd": ')
date_f = "2021-12-31"#input('Write the Final date to plot in the format "yyyy-mm-dd": ')
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
f_temp = np.loadtxt(data_path + "tempetatures_fideos.dat", skiprows=1, delimiter=" ", dtype="str")
T_date = f_temp[:,0]
T_date = T_date.astype(float)#-2459000
T1 = f_temp[:,1]
T1 = T1.astype(float)
T2 = f_temp[:,2]
T2 = T2.astype(float)
T3 = f_temp[:,3]
T3 = T3.astype(float)
T4 = f_temp[:,4]
T4 = T4.astype(float)


#Carga de datos de Seeing del DIMM La Silla.
f_dimm = np.loadtxt(data_path + "dimm_data.dat", skiprows=1, delimiter=" ", dtype="str")
d_date = f_dimm[:,0]
d_date = d_date.astype(float)#-2459000
d_seeing = f_dimm[:,1]
d_seeing = d_seeing.astype(float)


#Carga de datos de weather La Silla.
f_weather = np.loadtxt(data_path + "weather_wdb_all_data.dat", skiprows=1, delimiter=" ", dtype="str")
w_date = f_weather[:,0]
w_date = w_date.astype(float)#-2459000
w_press = f_weather[:,1]
w_press = w_press.astype(float)
w_amb_temp = f_weather[:,2]
w_amb_temp = w_amb_temp.astype(float)
w_dew_temp = f_weather[:,3]
w_dew_temp = w_dew_temp.astype(float)
w_rel_hum = f_weather[:,4]
w_rel_hum = w_rel_hum.astype(float)

#Carga de datos de wind de La Silla.
f_weather = np.loadtxt(data_path + "wind_data_all.dat", skiprows=1, delimiter=" ", dtype="str")
wind_date = f_weather[:,0]
wind_date = wind_date.astype(float)#-2459000
wind_dir = f_weather[:,1]
wind_dir = wind_dir.astype(float)
wind_speed = f_weather[:,2]
wind_speed = wind_speed.astype(float)


#Carga de datos de las RVs standard
f_rvs = np.loadtxt(data_path + "RVs_AZ_ALT.dat", skiprows=1, delimiter=" ", dtype="str")
RV_time = f_rvs[:,0] # tiempo UTC de la observacion
RV_time = RV_time.astype(str)
BJD = f_rvs[:,1]
BJD = BJD.astype(float)#-2459000
OBJ = f_rvs[:,2] # Nombre de la RV standard
OBJ = OBJ.astype(str)
RV = f_rvs[:,3]
RV = RV.astype(float)
eRV = f_rvs[:,4]
eRV = eRV.astype(float)
RV_ALT = f_rvs[:,5]
RV_ALT = RV_ALT.astype(float)
RV_AZ = f_rvs[:,6]
RV_AZ = RV_AZ.astype(float)
TEXP = f_rvs[:,7]
TEXP = TEXP.astype(float)
AIRMASS = f_rvs[:,8]
AIRMASS = AIRMASS.astype(float)
DRIFT_CO = f_rvs[:,9]
DRIFT_CO = DRIFT_CO.astype(float)
DRIFT_CO_E = f_rvs[:,10]
DRIFT_CO_E = DRIFT_CO_E.astype(float)
BS = f_rvs[:,11]
BS= BS.astype(float)
BS_E = f_rvs[:,12]
BS_E = BS_E.astype(float)
SNR = f_rvs[:,13]
SNR = SNR.astype(float)
SNR_R = f_rvs[:,14]
SNR_R = SNR_R.astype(float)
BAR_COR = f_rvs[:,15]
BAR_COR =BAR_COR.astype(float)

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
obj_SNR=list()
obj_alt=list()
obj_az=list()
obj_airm=list()
drift_fco=list()
bisec=list()
std_date=list()
std_texp=list()
for i in range(len(RV)): # RV data for plot
#	if SNR[i] > 80:
	if (BJD[i]>=JD_i) and (BJD[i]<=JD_f): #identifica rango de fechas
		if (SNR[i] >= low_snr): # Filtra por SNR mayor a low_snr
			for j in range(len(RV_star[0])):
				if name_rv(OBJ[i]) == RV_star[0][j]: # Separa RVs stars por nombre
					obj_HD10700 = RV_star[0][0]
					rv_HD10700.append(RV[i])
					date_HD10700.append(BJD[i])
					obj_SNR.append(SNR[i])
					obj_alt.append(RV_ALT[i])
					obj_az.append(RV_AZ[i])
					obj_airm.append(AIRMASS[i])
					drift_fco.append(DRIFT_CO[i])
					bisec.append(BS[i])
					std_date.append(BJD[i])
					std_texp.append(TEXP[i])
				elif name_rv(OBJ[i]) == RV_star[1][j]:
					obj_HD32147 = RV_star[1][0]
					rv_HD32147.append(RV[i])
					date_HD32147.append(BJD[i])
					obj_SNR.append(SNR[i])
					obj_alt.append(RV_ALT[i])
					obj_az.append(RV_AZ[i])
					obj_airm.append(AIRMASS[i])
					drift_fco.append(DRIFT_CO[i])
					bisec.append(BS[i])
					std_date.append(BJD[i])
					std_texp.append(TEXP[i])
				elif name_rv(OBJ[i]) == RV_star[2][j]:
					obj_HD72673 = RV_star[2][0]
					rv_HD72673.append(RV[i])
					date_HD72673.append(BJD[i])
					obj_SNR.append(SNR[i])
					obj_alt.append(RV_ALT[i])
					obj_az.append(RV_AZ[i])
					obj_airm.append(AIRMASS[i])
					drift_fco.append(DRIFT_CO[i])
					bisec.append(BS[i])
					std_date.append(BJD[i])
					std_texp.append(TEXP[i])
				elif name_rv(OBJ[i]) == RV_star[3][j]:
					obj_HD157347 = RV_star[3][0]
					rv_HD157347.append(RV[i])
					date_HD157347.append(BJD[i])
					obj_SNR.append(SNR[i])
					obj_alt.append(RV_ALT[i])
					obj_az.append(RV_AZ[i])
					obj_airm.append(AIRMASS[i])
					drift_fco.append(DRIFT_CO[i])
					bisec.append(BS[i])
					std_date.append(BJD[i])
					std_texp.append(TEXP[i])
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


f1=plt.figure(figsize=(12, 6))
plt.subplot(311)
p1=plt.plot(date_HD10700,rv_HD10700-rv_aver_HD10700,'r.',label='HD10700', markersize=3)#-rv_aver_HD10700
p2=plt.plot(date_HD32147,rv_HD32147-rv_aver_HD32147,'b.',label='HD32147', markersize=3)#-rv_aver_HD32147
p3=plt.plot(date_HD72673,rv_HD72673-rv_aver_HD72673,'g.',label='HD72673', markersize=3)#-rv_aver_HD72673
p4=plt.plot(date_HD157347,rv_HD157347-rv_aver_HD157347,'m.',label='HD157347', markersize=3)#-rv_aver_HD157347
plt.ylabel('RV - <RV> (km/s)')
plt.title(date_i+' to '+date_f)
plt.legend(bbox_to_anchor=(1.00, 0.8), loc='upper left', borderaxespad=0)
plt.axis([JD_i-1, JD_f+1, -0.20, 0.20])
#plt.legend((p1[0],p2[0],p3[0],p4[0]),("HD10700","HD32147","HD72673","HD157347"))
#plt.savefig('RV_graphs.png',dpi=600)
####################################
date_temp=list()
T1_fideos=list()
T2_fideos=list()
T3_fideos=list()
T4_fideos=list()
n=0
while n < len(T_date):
	if (T_date[n]>=JD_i) and (T_date[n]<=JD_f): #identifica rango de fechas
		date_temp.append(T_date[n])
		T1_fideos.append(T1[n])
		T2_fideos.append(T2[n])
		T3_fideos.append(T3[n])
		T4_fideos.append(T4[n])
	n=n+3
#####################################
plt.subplot(312)
p1=plt.plot(date_temp,T3_fideos,'g', markersize=3, label="Chiller")
p2=plt.plot(date_temp,T2_fideos,'b', markersize=3, label="Spec.")
p3=plt.plot(date_temp,T1_fideos,'r', markersize=3, label="Serv. room")
p4=plt.plot(date_temp,T4_fideos,'m', markersize=3, label="CCD")
plt.ylabel('Temperature deg C')
plt.legend(bbox_to_anchor=(1.00, 0.8), loc='upper left', borderaxespad=0)
plt.axis([JD_i-1, JD_f+1, 13.5, 20.5])
plt.legend((p1[0],p2[0],p3[0],p4[0]),("Chiller","Spec.","Serv. room","CCD"))
####################################
####################################
plt.subplot(313)
plot1=plt.plot(std_date,obj_SNR,'b.', markersize=3)
plt.ylabel('SNR')
plt.xlabel('Julian date')
plt.axis([JD_i-1, JD_f+1, 20, 180])
plt.show()
#
#f1.savefig('RV_temp_graphs_all_snr1.png',dpi=600)


###########################################################
fig, axs = plt.subplots(3, 1, figsize=(12, 6), sharex=True)
fig.subplots_adjust(hspace=0)
#ax1.xaxis.set_major_formatter()
f1=axs[0].plot(date_HD10700,rv_HD10700-rv_aver_HD10700,'r.',label='HD10700',markersize=3)#-rv_aver_HD10700
f2=axs[0].plot(date_HD32147,rv_HD32147-rv_aver_HD32147,'b.',label='HD32147', markersize=3)#-rv_aver_HD32147
f3=axs[0].plot(date_HD72673,rv_HD72673-rv_aver_HD72673,'g.',label='HD72673', markersize=3)#-rv_aver_HD72673
f4=axs[0].plot(date_HD157347,rv_HD157347-rv_aver_HD157347,'m.',label='HD157347', markersize=3)#-rv_aver_HD157347
axs[0].set_ylabel('RV - <RV> (km/s)')
#ax1.set_title(date_i+' to '+date_f)
#ax.set_legend(bbox_to_anchor=(1.00, 0.8), loc='upper left', borderaxespad=0)
#ax1.set_axis([JD_i-1, JD_f+1, -0.10, 0.10])
#plt.legend((p1[0],p2[0],p3[0],p4[0]),("HD10700","HD32147","HD72673","HD157347"))
#plt.savefig('RV_graphs.png',dpi=600)
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
plt.show()