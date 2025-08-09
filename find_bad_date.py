import numpy as np
from matplotlib.dates import DateFormatter
from astropy.table import Table
from astropy.time import Time
from datetime import datetime

#Datos de rango de fechas para analizar
date_i = "2018-01-15" #"2020-10-01"#input('Write the Initial date to plot in the format "yyyy-mm-dd": ')
date_f = "2023-07-15"#input('Write the Final date to plot in the format "yyyy-mm-dd": ')
#Convierte fecha inicial de busqueda en JD
date_i = date_i+'T'+"05:05:00.0"
time_i = Time(date_i,format="isot",scale="utc")
JD_i = time_i.jd

#Carga de datos de temperatura de fideos.

f_temp = np.loadtxt("fideos_temperatures.dat", skiprows=1, delimiter=" ", dtype="str")
T_date = f_temp[:,0]
T_date = T_date.astype(str)
Month=list()
for i in range(0,len(T_date)):
    a = T_date[i].split("-")
    b = float(a[1])
    if b > 12.:
        print(b,T_date[i])