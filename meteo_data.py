#
# Usar este programa para modificar los datos de DIMM y WEATHER de La Silla
#
import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.table import Table

#Carga de datos de Seeing del DIMM La Silla.
f_dimm = pd.read_csv("wdb_query_47395_eso_dimm.csv",skiprows=[0])

dimm_date = f_dimm.Date_time.to_numpy()
dimm_seeing = f_dimm.DIMM_Seeing.to_numpy()

# Creo tabla para guardar los datos del DIMM con date en JD
save = Table([[], []],names=("#JD_time","DIMM"),dtype=(float,float))

#Cambio fecha a JD y agrego los elementos de la nueva tabla de DIMM
j=0
for i in range(len(dimm_date)):#len(f_DATE)):
	d = dimm_date[i]
#	print(i)
	time = Time(d,format="isot",scale="utc")
	dimm_JD = time.jd
	row = [dimm_JD,dimm_seeing[i]]
	save.add_row(row)
	if j==i:
		print(i)
		j=j+1000
save.write("dimm_data.dat",format='ascii',overwrite=True)
print('Done *dimm_data.dat*')

#Carga de datos de weather La Silla.
f_weather = pd.read_csv("weather_wdb_query_40466_eso.csv",skiprows=[0])

w_date = f_weather.Date_time.to_numpy()
w_press = f_weather.Air_Pressure.to_numpy()
w_amb_temp = f_weather.Ambient_temp.to_numpy()
w_dew_temp = f_weather.Dew_temp.to_numpy()
w_rel_hum = f_weather.Rel_humidity.to_numpy()

# Creo tabla para guardar los datos del weather con date en JD
save = Table([[], [], [], [], []],names=("#JD_time", "Air_Pressure", "Ambient_temp", "Dew_temp", \
  "Rel_humidity"),dtype=(float,float,float,float,float))

#Cambio fecha a JD y agrego los elementos de la nueva tabla de weather
j=0
for i in range(len(w_date)):
	d = w_date[i]
	time = Time(d,format="isot",scale="utc")
	w_JD = time.jd
	if w_JD < 2458484.50:
		row = [w_JD,w_press[i],w_amb_temp[i],w_dew_temp[i],w_rel_hum[i]]
		save.add_row(row)
		if j==i:
			print(i)
			j=j+1000
save.write("weather_data_2018.dat",format='ascii',overwrite=True)
print('Done *weather_data_2018.dat*')
#
j=0
for i in range(len(w_date)):
	d = w_date[i]
	time = Time(d,format="isot",scale="utc")
	w_JD = time.jd
	if (w_JD > 2458484.50) and (w_JD < 2458849.50):
		row = [w_JD,w_press[i],w_amb_temp[i],w_dew_temp[i],w_rel_hum[i]]
		save.add_row(row)
		if j==i:
			print(i)
			j=j+1000
save.write("weather_data_2019.dat",format='ascii',overwrite=True)
print('Done *weather_data_2019.dat*')
#
j=0
for i in range(len(w_date)):
	d = w_date[i]
	time = Time(d,format="isot",scale="utc")
	w_JD = time.jd
	if (w_JD > 2458849.50) and (w_JD < 2459215.50):
		row = [w_JD,w_press[i],w_amb_temp[i],w_dew_temp[i],w_rel_hum[i]]
		save.add_row(row)
		if j==i:
			print(i)
			j=j+1000
save.write("weather_data_2020.dat",format='ascii',overwrite=True)
print('Done *weather_data_2020.dat*')
#
j=0
for i in range(len(w_date)):
	d = w_date[i]
	time = Time(d,format="isot",scale="utc")
	w_JD = time.jd
	if (w_JD > 2459215.50) and (w_JD < 2459580.50):
		row = [w_JD,w_press[i],w_amb_temp[i],w_dew_temp[i],w_rel_hum[i]]
		save.add_row(row)
		if j==i:
			print(i)
			j=j+1000
save.write("weather_data_2021.dat",format='ascii',overwrite=True)
print('Done *weather_data_2021.dat*')	
#
j=0
for i in range(len(w_date)):
	d = w_date[i]
	time = Time(d,format="isot",scale="utc")
	w_JD = time.jd
	if (w_JD > 2459580.50):
		row = [w_JD,w_press[i],w_amb_temp[i],w_dew_temp[i],w_rel_hum[i]]
		save.add_row(row)
		if j==i:
			print(i)
			j=j+1000
save.write("weather_data_2022.dat",format='ascii',overwrite=True)
print('Done *weather_data_2022.dat*')



#Carga de datos de wind de La Silla (separado de weather).
f_wind = pd.read_csv("wind_wdb_query_42005_eso.csv",skiprows=[0])
wind_date = f_wind.Date_time.to_numpy()
wind_dir = f_wind.Wind_direc.to_numpy()
wind_speed = f_wind.Wind_speed.to_numpy()

# Creo tabla para guardar los datos del weather con date en JD
save = Table([[], [], []],names=("#JD_time", "Wind_direc", "Wind_speed"),dtype=(float,float,float))

#Cambio fecha a JD y agrego los elementos de la nueva tabla de wind
j=0
for i in range(len(wind_date)):
	d = wind_date[i]
	time = Time(d,format="isot",scale="utc")
	wind_JD = time.jd
	if (wind_JD < 2458484.50):		
		row = [wind_JD,wind_dir[i],wind_speed[i]]
		save.add_row(row)
		if j==i:
			print(i)
			j=j+1000
save.write("wind_data_2018.dat",format='ascii',overwrite=True)
print('Done *wind_data_2018.dat*')
#
j=0
for i in range(len(wind_date)):
	d = wind_date[i]
	time = Time(d,format="isot",scale="utc")
	wind_JD = time.jd
	if (wind_JD > 2458484.50) and (wind_JD < 2458849.50):
		row = [wind_JD,wind_dir[i],wind_speed[i]]
		save.add_row(row)
		if j==i:
			print(i)
			j=j+1000
save.write("wind_data_2019.dat",format='ascii',overwrite=True)
print('Done *wind_data_2019.dat*')
#
j=0
for i in range(len(wind_date)):
	d = wind_date[i]
	time = Time(d,format="isot",scale="utc")
	wind_JD = time.jd
	if (wind_JD > 2458849.50) and (wind_JD < 2459215.50):
		row = [wind_JD,wind_dir[i],wind_speed[i]]
		save.add_row(row)
		if j==i:
			print(i)
			j=j+1000
save.write("wind_data_2020.dat",format='ascii',overwrite=True)
print('Done *wind_data_2020.dat*')
#
j=0
for i in range(len(wind_date)):
	d = wind_date[i]
	time = Time(d,format="isot",scale="utc")
	wind_JD = time.jd
	if (wind_JD > 2459215.50) and (wind_JD < 2459580.50):
		row = [wind_JD,wind_dir[i],wind_speed[i]]
		save.add_row(row)
		if j==i:
			print(i)
			j=j+1000
save.write("wind_data_2021.dat",format='ascii',overwrite=True)
print('Done *wind_data_2021.dat*')
#
j=0
for i in range(len(wind_date)):
	d = wind_date[i]
	time = Time(d,format="isot",scale="utc")
	wind_JD = time.jd
	if (wind_JD > 2459580.50):
		row = [wind_JD,wind_dir[i],wind_speed[i]]
		save.add_row(row)
		if j==i:
			print(i)
			j=j+1000
save.write("wind_data_2022.dat",format='ascii',overwrite=True)
print('Done *wind_data_2022.dat*')