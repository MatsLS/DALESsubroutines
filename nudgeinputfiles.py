"""This script can be used to create input files for 3D (nudge3D.inp.iexpnr) and 1D nudging (nudge.inp.iexpnr) in DALES. Here, 'iexpnr' is the experiment number as
named in namoptions."""
# Import packages
import numpy as np
import netCDF4 as nc
import pandas as pd
from scipy.io import FortranFile

# Path for nc files (DALES output) with the data to create the input files.
folderpt = r'C:\Users\matss\Documents\AP\MEP\DALES_results\Astex_nudge\H'
iexpnr = '112'
# File extension
ext = '_multfields112'

# Data about the nc files
Nprocx = 1
itot = 256
Nprocy = 32
jtot = 256
k = 427
ti = 0
tnstart = 1
tnstop = 12

# Read the file for the first processor to read the time and height already.
fn = folderpt + iexpnr + r'\fielddump.000.000.' + iexpnr + r'.nc'
df = nc.Dataset(fn)

# Convert time
t = df['time'][tnstart : tnstop + 1] - 18000
zt = df['zt'][:]

df.close()

# Create arrays with zeros for the variables in the nc files.
xt = np.zeros(itot)
yt = np.zeros(jtot)
thl = np.zeros((tnstop + 1 - tnstart, k, jtot, itot))
qt = np.zeros((tnstop + 1 - tnstart, k, jtot, itot))

# Maximum i and j indexes per processor
imax = int(itot / Nprocx)
jmax = int(jtot / Nprocy)

# Loop over the amount of used processors and open the nc file for each, then add variables to their array.
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
        thl[:, :, int(j * jmax) : int((j + 1) * jmax), int(i * imax) : int((i + 1) * imax)] = df['thl'][tnstart : tnstop + 1, :, :jmax, :imax]
        qt[:, :, int(j * jmax) : int((j + 1) * jmax), int(i * imax) : int((i + 1) * imax)] = df['qt'][tnstart : tnstop + 1, :, :jmax, :imax]
        
        df.close()

# Create the header for the 1D nudge input file
columns = ['height', 't_nudge', 'u_nudge', 'v_nudge', 'w_nudge', 'thl_nudge', 'qt_nudge']
header = columns[0]
i=1
while i in range(len(columns)):
    header = header + '\t' + columns[i]
    i += 1
zeros = np.zeros(len(zt))
tnudge = np.ones(len(zt))

# Print the header and data for each nudging time step
l = 0
with open(r'C:\Users\matss\Documents\AP\MEP\DALES_results\Astex_nudge\nudge'+ ext + r'.inp', 'wb') as myfile:
    myfile.write(b'ASTEX large domain \n')
    for i in range(0, tnstop + 1 - tnstart):
        header1 = header + '\n # ' + str(t[i])
        # This 1D input file uses mean values of the reference run, like below. zeros for u,v,w indicate no nudging for these.        
        data = {columns[0] : zt, columns[1] : tnudge, columns[2] : zeros, columns[3] : zeros, columns[4] : zeros, 
        columns[5] : np.mean(thl[i,:,:,:], axis = (1,2)), columns[6] : np.mean(qt[i,:,:,:], axis = (1,2))}

        nudge1 = pd.DataFrame(data = data)
        
        np.savetxt(myfile, nudge1.values, delimiter = '\t', header = header1, fmt = '%4.7f', comments = '')
        myfile.write(b'\n')
        
    data = {columns[0] : zt, columns[1] : tnudge, columns[2] : zeros, columns[3] : zeros, columns[4] : zeros, 
        columns[5] : zeros, columns[6] : zeros}
    nudge1 = pd.DataFrame(data = data)
    # Add a nudging time step at the end of the run, needed for standard DALES nudging (not used but necessary in input file)
    header2 = header + '\n # 10800.0'
    np.savetxt(myfile, nudge1.values, delimiter = '\t', header = header2, fmt = '%4.7f', comments = '')

# Create arrays for the 3D nudging desired fields
tnudge3D = np.zeros(tnstop - tnstart + 1, dtype = np.float64)
thlnudge3D = np.zeros([tnstop - tnstart + 1, k, jtot, itot], dtype=np.float64)
qtnudge3D = np.zeros([tnstop - tnstart + 1, k, jtot, itot], dtype=np.float64)

# The nudging time points
tnudge3D[:] = t
# The fluctuation fields
thlnudge3D = thl[:, :, :, :] - np.mean(thl[:, : ,:, :], axis = (2,3))[..., np.newaxis, np.newaxis]
qtnudge3D = qt[:, :, :, :] - np.mean(qt[:, : ,:, :], axis = (2,3))[..., np.newaxis, np.newaxis]

# Write the created arrays to the nudge3D input file.
f = FortranFile(r'C:\Users\matss\Documents\AP\MEP\DALES_results\Astex_nudge\nudge3D' + ext + r'_f1.inp', 'w')
f.write_record(tnudge3D)
f.write_record(thlnudge3D[:, :, :, :])
f.write_record(qtnudge3D[:, :, :, :])
f.close()