import numpy as np
from astropy.time import Time
import astropy.io; from astropy.io import fits as pyfits
from astropy.table import Table
import subprocess
#import sys
#import os

#subprocess.call(['./create_RV_spec2021'])

'''
 Now is read the file create before named "RVs_spec2021.list". It contains the path to
 the reduced spectra "*.sp.fits"
'''
file = "RVs_spec2021.list"
image = np.loadtxt(file, skiprows=0, delimiter=" ", dtype="str")


# Se crea tabla que guarda datos de DATE, OBJ, RA y DEC.
save = Table([[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []], \
        names=("#DATE_OBS", "BJD", "OBJ", "RV", "eRV", "RA", "DEC", \
                "TEXP", "AIRMASS", "DRIFT_CO", "DRIFT_CO_E", "BS", "BS_E", "SNR", "SNR_R", "BAR_COR"), \
        dtype=('S26', float,'S8', float, float, float, float, float, float, float, float, float, float, \
                float, float, float))

#hdr Muestra todo el Header igual que el cmd print abajo
#print(repr(hdr))
# A continuacion se leen los headers de los espectros reducidos
for i in range(0,len(image)):
        h = pyfits.open(image[i])
        img = h[0].data
        hdr = h[0].header
        RA = hdr['HIERARCH RA']
        DEC = hdr['HIERARCH DEC']
        OBJECT = hdr['HIERARCH TARGET NAME']
        RV = hdr['RV']
        eRV = hdr['RV_E']
        DATE_OBS = hdr['HIERARCH SHUTTER START DATE']
        TIME_OBS = hdr['HIERARCH SHUTTER START UT']
        date = DATE_OBS +'T'+TIME_OBS
        print(OBJECT,RV,date)
        # Descomentar los siguientes campos para leer otros datos del HEADER
        #MJD = hdr['HIERARCH MJD']
        #MBJD = hdr['HIERARCH MBJD']
        TEXP = hdr['HIERARCH TEXP (S)']
        BAR_COR = hdr['HIERARCH BARYCENTRIC CORRECTION (KM/S)']
        #EQUINOX = hdr['HIERARCH EQUINOX']
        DRIFT_CO = hdr['HIERARCH DRIFT_CO']
        DRIFT_CO_E = hdr['HIERARCH DRIFT_CO_E']
        AIRMASS = hdr['HIERARCH AIRMASS']
        #BDATE = hdr['BDATE']
        #REFWAV = hdr['HIERARCH REFWAV']
        BS = hdr['BS']
        BS_E = hdr['BS_E']
        #DISP = hdr['DISP']
        SNR = hdr['SNR']
        SNR_R = hdr['SNR_R']
        #INST = hdr['INST']
        #RESOL = hdr['RESOL']
        BJD = hdr['BJD_OUT']
        #DATE_RED = hdr['DATE']

        #
        row = [date,BJD,OBJECT,RV,eRV,RA,DEC,TEXP,AIRMASS,DRIFT_CO,DRIFT_CO_E,BS,BS_E,SNR,SNR_R,BAR_COR]
        save.add_row(row)
#table_dat = save

# Ahora se guarda la tabla
#output_name = raw_input('Write the output file name: ')
save.write("headers_keyW_2021.dat", format='ascii')#, overwrite=True)

# Llamo a un script para ordenar y eliminar datos repetidos. Entrega el archivo de salida
# "headers_keyW_sorted.dat"
subprocess.call(['./data_sort'])

#
# !!! Editar tabla con vi (titulo al final de archivo.dat)
#