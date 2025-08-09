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
#mport os
import time
start = time.clock()

file_list = np.loadtxt("files.list", skiprows=0, dtype="|S")

#
# PARECE que es mejor trabajar con cada archivo por separado
# y no juntarlos de inmediato.
#


#files = []
raw_time = []
offset = [-3,-3,-3,-3,-4,-4,0]
offset = offset*u.hour

save = Table([[], [], [], [], []], names=("#JD_OBS", "Chiller", "Op_Bench", "Room", "CCD"), \
	dtype=('S20', float, float, float, float))

#
for i in range(0,len(file_list)): # Coincide que file_list y files tienen el mismo largo.
	file = np.loadtxt(file_list[i], skiprows=1, delimiter=" ",comments="#", dtype="|S")
	date = file[:,0]
	date = date.astype(str)
	hour = file[:,1]
	hour = hour.astype(str)
	Chiller = file[:,2]
	Chiller = Chiller.astype(np.float)
	Op_Bench = file[:,3]
	Op_Bench = Op_Bench.astype(np.float)
	Room = file[:,4]
	Room = Room.astype(np.float)
	CCD = file[:,5]
	CCD = CCD.astype(np.float)
	print(i)
	print("Room ", len(Room))
	print("date ",len(date))
	for j in range(0,len(date)):
		T_time = date[j]+'T'+hour[j]
		#raw_time.append(T_time)
		time = Time(T_time, format="isot", scale="utc") - offset[i]
		JD = time.jd
		#print(JD,Chiller[j],Op_Bench[j],Room[j],CCD[j])
		row = [JD,Chiller[j],Op_Bench[j],Room[j],CCD[j]]
		print(j)
		save.add_row(row)
#
save.write("fideos_temp_data.dat", format="ascii")