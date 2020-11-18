#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""""
Python Module for Importing and Reformatting Ocean Reanalysis System 5 (ORA-S5) Data from the Tropical Atlantic Ocean
by Ryan Glaubke
Date Created: 07/19/20

"""""

# Import packages
import numpy as np
import xarray as xr
from datetime import datetime
import scipy.io
import os

#Enter File Directory where you would like the .mat file to be saved
save_dir = os.getcwd()

# Import ORA-S5 -- Potential Temperature values from 1958 - 2018
# NOTE: the period between 1958 - 1978 is a backwards extension; parameters and resolution are the same as ORA-S5 (refer to Zuo et al. [2019] for more details)
filenames = []
yr = np.arange(1958,2019,1)
mth = np.arange(1,13,1)
for ynum in yr:
    if ynum in range(1958,1979):
        for mnum in mth:
            if mnum < 10:
                url = 'http://icdc.cen.uni-hamburg.de/thredds/dodsC/ftpthredds/EASYInit/oras5_backward_extension/r1x1/votemper/opa0/votemper_ORAS5_1m_'+np.str(ynum)+'0'+np.str(mnum)+'_r1x1.nc'
                filenames.append(url)
            else: 
                url = 'http://icdc.cen.uni-hamburg.de/thredds/dodsC/ftpthredds/EASYInit/oras5_backward_extension/r1x1/votemper/opa0/votemper_ORAS5_1m_'+np.str(ynum)+np.str(mnum)+'_r1x1.nc'
                filenames.append(url)
    else:
        for mnum in mth:
            if mnum < 10:
                url = 'http://icdc.cen.uni-hamburg.de/thredds/dodsC/ftpthredds/EASYInit/oras5/r1x1/votemper/opa4/votemper_ORAS5_1m_'+np.str(ynum)+'0'+np.str(mnum)+'_r1x1.nc'
                filenames.append(url)
            else: 
                url = 'http://icdc.cen.uni-hamburg.de/thredds/dodsC/ftpthredds/EASYInit/oras5/r1x1/votemper/opa4/votemper_ORAS5_1m_'+np.str(ynum)+np.str(mnum)+'_r1x1.nc'
                filenames.append(url)
#print(*filenames, sep='\n') # <-- Displays list of all netCDF files in 'filenames'
data = xr.open_mfdataset(filenames,concat_dim='time_counter',combine='by_coords')

# Subsample ORA-S5 to extract data spanning the Tropical Atlantic (30°S - 30°N; 100°W - 15°E)
# Must split datasets then concatenate (xarray does not like slicing across the meridian...)
sub1 = data.sel(lat=slice(-30,30), lon=slice(259.5,359.5))
sub2 = data.sel(lat=slice(-30,30), lon=slice(0.5,15.5))
TA_datasub = xr.concat([sub1,sub2],dim='lon')

# Convert subsampled "Tropical Atlantic" dataset to .mat file
# First, rearrange netCDF file into the following format: (lon, lat, deptht, time_counter) <-- It's just easier to think about it in that format
TA_datasub = TA_datasub.transpose()
# Next, convert to .mat file that can be uploaded into MATLAB
os.chdir(save_dir)
scipy.io.savemat('TA_ORAS5.mat', mdict={'ORAS5': TA_datasub, 'depths': TA_datasub.deptht})


# In[ ]:




