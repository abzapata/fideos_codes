import numpy as np
from astropy.table import Table
from astropy.time import Time
from multiprocessing import Pool

data_path = "../all_data/"

f_wind = np.loadtxt(data_path + "wind_wdb_query_2240_eso.csv", skiprows=1, delimiter=",", dtype="str")
iso_date = f_wind[:,0]
iso_date = iso_date.astype(str)# Date in ISOT format
w_direc = f_wind[:,1]
w_direc = w_direc.astype(float)
w_speed = f_wind[:,2]
w_speed = w_speed.astype(float)

save = Table([[], [], [], []], names=("#DATE_OBS", "JD", "Wind_Direction", "Wind_Speed"), \
	dtype=('S20', float, float, float))
j = 0
with Pool(5) as p:
    while j < 2732710:
        print(j)
        T_date = iso_date[j]
        time = Time(T_date, format="isot", scale="utc")
        JD = time.jd
        row = [iso_date[j],JD,w_direc[j],w_speed[j]]
        save.add_row(row)
        j = j+1
save.write(data_path + "wind_data_all.dat", format="ascii", overwrite=True)