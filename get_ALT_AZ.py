'''
Crea tabla con datos de ALT, AZ de RVs a partir de RA y DEC.
para observatorio La Silla
'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import astropy.io; from astropy.io import fits as pyfits
from astropy.table import Table
from astropy.io import ascii
from astropy.visualization import astropy_mpl_style
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import subprocess
#import sys
#import os
import time
start = time.clock()

#your code here    
#########
plt.style.use(astropy_mpl_style)


# Cargo tabla con datos de RV (DATE,OBS,RA,DEC) y convierto a ALT,AZ, y
# convierto fecha UTC a JD

# Fechas de comienzo de meses del 2018
months = ["2018-01-01T00:00:00.0","2018-02-01T00:00:00.0","2018-03-01T00:00:00.0",\
"2018-04-01T00:00:00.0","2018-05-01T00:00:00.0","2018-06-01T00:00:00.0",\
"2018-07-01T00:00:00.0","2018-08-01T00:00:00.0","2018-09-01T00:00:00.0",\
"2018-10-01T00:00:00.0","2018-11-01T00:00:00.0","2018-12-01T00:00:00.0"]
 
month_JD = list()

for i in range(0,len(months)):
	time = Time(months[i], format="isot",scale="utc")
	month_JD.append(time.jd)

# Ahora los meses estan en JD y guardados en "month_JD[i]"


### SITE Coord. (La Silla Obs. in this case)
OBS_ELEV = 2335
OBS_LAT  = -29.2543
OBS_LON  = -70.7346
LSO = EarthLocation(lat=-29.2543*u.deg, lon=-70.7346*u.deg, height=2335*u.m)

# Archivo headers_keyW_sorted obtenido con read_header.py
f_RADEC = np.loadtxt("headers_keyW_sorted.dat", skiprows=1, \
	comments='#', delimiter=" ", dtype="|S")

#Lista de RVs (diferentes catalogos):
# Los nombres ya estan en la tabla de entrada, asi que no es necesario
# volver a hacer el match.
#RVs_names = [['HD10700','HIP8102','HR509'],['HD32147','HIP23311','HR1614'],\
#['HD72673','HIP41926','HR3384'],['HD157347','HIP85042','HR6465']]

RV_time = f_RADEC[:,0] # tiempo UTC de la observacion
RV_time = RV_time.astype(str)
RV_OBJ = f_RADEC[:,1] # Nombre de la RV standard
RV_OBJ = RV_OBJ.astype(str)
RV = f_RADEC[:,2]
RV = RV.astype(np.float)
eRV = f_RADEC[:,3]
eRV = eRV.astype(np.float)
RV_RA = f_RADEC[:,4]
RV_RA = RV_RA.astype(np.float)
RV_DEC = f_RADEC[:,5]
RV_DEC = RV_DEC.astype(np.float)
TEXP = f_RADEC[:,6]
TEXP = TEXP.astype(np.float)
BAR_COR = f_RADEC[:,7]
BAR_COR =BAR_COR.astype(np.float)
DRIFT_CO = f_RADEC[:,8]
DRIFT_CO = DRIFT_CO.astype(np.float)
DRIFT_CO_E = f_RADEC[:,9]
DRIFT_CO_E = DRIFT_CO_E.astype(np.float)
AIRMASS = f_RADEC[:,10]
BS = f_RADEC[:,11]
BS= BS.astype(np.float)
BS_E = f_RADEC[:,12]
BS_E = BS_E.astype(np.float)
AIRMASS = AIRMASS.astype(np.float)
SNR = f_RADEC[:,13]
SNR = SNR.astype(np.float)
SNR_R = f_RADEC[:,14]
SNR_R = SNR_R.astype(np.float)
BJD = f_RADEC[:,15]
BJD = BJD.astype(np.float)

#time = Time(RV_time)

save = Table([[], [], [], [], [], [], \
	[], [], [], [], [], [], [], [], [], []],\
	names=("#JD_OBS", "OBJ", "RV", "eRV", "ALT", "AZ", \
		"TEXP", "BAR_COR", "DRIFT_CO", "DRIFT_CO_E", "AIRMASS", "BS", "BS_E", "SNR", "SNR_R", "BJD"), \
	dtype=('S26','S8', float, float, float, float, \
		float, float, float, float, float, float, float, float, float, float))


for i in range(0,len(RV_time)):
	time = Time(RV_time[i], format="isot",scale="utc")
	JD = time.jd
	#Convierte a coordenadas ALT,AZ desde RA y DEC
	rv_object = SkyCoord(RV_RA[i],RV_DEC[i],unit="deg",frame="icrs")
	object_altaz = rv_object.transform_to(AltAz(obstime=time,location=LSO))
	OBJ_ALT = object_altaz.alt
	OBJ_AZ = object_altaz.az
	#print(i)
	row = [JD,RV_OBJ[i],RV[i],eRV[i],OBJ_ALT,OBJ_AZ,TEXP[i],BAR_COR[i],DRIFT_CO[i],\
	DRIFT_CO_E[i],AIRMASS[i],BS[i],BS_E[i],SNR[i],SNR_R[i],BJD[i]]		
	save.add_row(row)
save.write("RVs_AZ_ALT.dat", format='ascii', overwrite='true')

import time
end = time.clock()
print(end - start)
