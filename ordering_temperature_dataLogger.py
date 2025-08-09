import numpy as np
from astropy.table import Table
from astropy.time import Time

save = Table([[], [], [], [], []], names=("#DATE_OBS", "Chiller", "Op_Bench", "Room", "CCD"), \
	dtype=('S19', float, float, float, float))
#
file = np.loadtxt("dataLogger.csv", skiprows=1, delimiter=",",comments="#", dtype="str")
date = file[:,0]
date = date.astype(str)
#hour = file[:,1]
#hour = hour.astype(str)
Chiller = file[:,1]
Chiller = Chiller.astype(np.float64)
Op_Bench = file[:,3]
Op_Bench = Op_Bench.astype(np.float64)
Room = file[:,5]
Room = Room.astype(np.float64)
CCD = file[:,7]
CCD = CCD.astype(np.float64)

time = list()
for i in range(0,len(date)):
	a = date[i].split(" ")
	b = a[0].split("/")
	time.append(a[1])
	dd = b[0]
	mm = b[1]
	yy = "20"+b[2] # cambio el formato de YY a YYYY
	isot_date = yy+"-"+mm+"-"+dd+"T"+a[1]
#	n_time = Time(isot_date, format="isot", scale="utc")
#	JD = n_time.jd
	print(isot_date,Chiller[i],Op_Bench[i],Room[i],CCD[i])
	row = [isot_date,Chiller[i],Op_Bench[i],Room[i],CCD[i]]
	save.add_row(row) 
save.write("fideos_temp_dataLogger.dat", format="ascii",overwrite='True')