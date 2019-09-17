# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 15:42:43 2019

@author: fli034
"""

import glob as glb
import segyio

dirName = r"C:\Users\Fli034\Documents\firstArrivalPicker\eRFP_Develop\first-arrival-picker\dataFiles\survey1"
files = glb.glob(dirName + '\*.sgy')
if files == []:
    files = glb.glob(dirName + '\*.segy')
    
if files == []:
    print('No files with *.sgy or *.segy exist in this directory')
    
#Column 1: File Name (str)
#Column 2: SX (float)
fileInfo = []

for i in range(0,len(files)):
    file = files[i]
    tmp = files[i].split(sep='\\')
    print(tmp[-1])
    with segyio.open(file,strict = False) as f:
        shotLoc = f.header[0][segyio.TraceField.SourceX]
        print(shotLoc)
    fileInfo.append([tmp[-1],shotLoc])