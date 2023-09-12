"""This script can be used to create an input file for the persistence method and the 3D initial thermodynamics method added to DALES.
To use them in DALES, they need to be called 'advfield.inp.iexpnr' or 'coldinifield.inp.iexpnr', where iexpnr is the experiment number as
named in namoptions."""

import numpy as np
import netCDF4 as nc

from scipy.io import FortranFile

# Extension for saved file
file_ext = 'test'

# Folder and experiment number from which to create the input files (reference run)
folderpt = r'C:\Users\matss\Documents\AP\MEP\DALES_results\Astex_nudge\H'
iexpnr = '082'

# Time index at which fields should be extracted
timep = 0

# Reference run simulation data
Nprocx = 1
itot = 256
Nprocy = 32
jtot = 256
k = 427
ti = 13
ktot = 427

xt = np.zeros(itot)
yt = np.zeros(jtot)
thl = np.zeros((k, jtot, itot))
qt = np.zeros((k, jtot, itot))

imax = int(itot / Nprocx)
jmax = int(jtot / Nprocy)

for j in range(0, Nprocy):
    for i in range(0, Nprocx):
        if j < 10 and i < 10:
            fn = folderpt + iexpnr + r'\fielddump.00' + str(i) + r'.00' + str(j) + r'.' + iexpnr + r'.nc'
        if j >= 10 and i < 10:
            fn = folderpt + iexpnr + r'\fielddump.00' + str(i) + r'.0' + str(j) + r'.' + iexpnr + r'.nc'
        if j < 10 and i >= 10:
            fn = folderpt + iexpnr + r'\fielddump.0' + str(i) + r'.00' + str(j) + r'.' + iexpnr + r'.nc'
        if j >= 10 and i >= 10:
            fn = folderpt + iexpnr + r'\fielddump.0' + str(i) + r'.0' + str(j) + r'.' + iexpnr + r'.nc'
        
        df = nc.Dataset(fn)

        xt[int(i * imax) : int((i + 1) * imax)] = df['xt'][:imax]
        yt[int(j * jmax) : int((j + 1) * jmax)] = df['yt'][:jmax]
        thl[:, int(j * jmax) : int((j + 1) * jmax), int(i * imax) : int((i + 1) * imax)] = df['thl'][timep, :, :jmax, :imax]
        qt[:, int(j * jmax) : int((j + 1) * jmax), int(i * imax) : int((i + 1) * imax)] = df['qt'][timep, :, :jmax, :imax]
        t = df['time'][:]
        zt = df['zt'][:]
        df.close()

# Save file
f = FortranFile(r'C:\Users\matss\Documents\AP\MEP\DALES_results\advOnly\advfield_082.inp', 'w')
f.write_record(thl)
f.write_record(qt)
f.close()
