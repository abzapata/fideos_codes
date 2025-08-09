# Separate files
import numpy as np
from astropy.table import Table
from astropy.time import Time

#Datos de rango de fechas para analizar
date_i = "2018-01-15"#input('Write the Initial date to plot in the format "yyyy-mm-dd": ')
date_f = "2021-12-31"#input('Write the Final date to plot in the format "yyyy-mm-dd": ')

#Convierte fecha inicial de busqueda en JD
date_i = date_i+'T'+"00:00:00.0"
time_i = Time(date_i,format="isot",scale="utc")

#Convierte fecha final de búsqueda en JD
date_f = date_f+'T'+"00:00:00.0"
time_f = Time(date_f,format="isot",scale="utc")
JD_f = time_f.jd

print(JD_i, JD_f)
#################