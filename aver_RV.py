import numpy as np
from astropy.time import Time
data_path="../all_data/"
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

print(JD_i, JD_f)
rv_aver = list()
for i in range(0,len(RV_time)):
    for j in range(len(RV_star[0])):
        if name_rv(OBJ[i]) == RV_star[0][j]: # Separa RVs stars por nombre HD10700
            print(name_rv(OBJ[i]))
        if name_rv(OBJ[i]) == RV_star[1][j]: # Separa RVs stars por nombre HD32147
            print(name_rv(OBJ[i]))
        if name_rv(OBJ[i]) == RV_star[2][j]: # Separa RVs stars por nombre HD72673
            print(name_rv(OBJ[i]))
        if name_rv(OBJ[i]) == RV_star[3][j]: # Separa RVs stars por nombre HD157347
            print(name_rv(OBJ[i]),RV_time[i])