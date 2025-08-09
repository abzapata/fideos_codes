import numpy as np
from astropy.table import Table
from astropy.time import Time
# Este script separa los datos de temperatura de desde el 2018 hasta 2023, cerca de 3 mill. datos 
# en tres archivos de 1 mill. de datos para agilizar el script. Estos datos estan en JD y se transforman
# en fechas iso.

data_path = "../all_data/"
f_temp = np.loadtxt(data_path + "fideos_temperatures.dat", skiprows=1, delimiter=" ", dtype="str")
date = f_temp[:,0]
date = date.astype(str)# Date in ISOT format
Chiller = f_temp[:,1]
Chiller = Chiller.astype(float)
Op_Bench = f_temp[:,2]
Op_Bench = Op_Bench.astype(float)
Room = f_temp[:,3]
Room = Room.astype(float)
CCD = f_temp[:,4]
CCD = CCD.astype(float)

save = Table([[], [], [], [], [], []], names=("#DATE_OBS", "JD", "Chiller", "Op_Bench", "Room", "CCD"), \
	dtype=('S20', float, float, float, float, float))
j = 1000000
while j < 2000000:
    print(j)
    T_date = date[j]
    time = Time(T_date, format="isot", scale="utc")
    JD = time.jd
    row = [date[j],JD,Chiller[j],Op_Bench[j],Room[j],CCD[j]]
    save.add_row(row)
    j = j+1
save.write(data_path + "fideos_temp_data2.dat", format="ascii", overwrite=True)