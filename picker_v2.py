# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 20:39:09 2018

@author: bflinch1
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.widgets import Slider
import time
import os
from scipy import signal
import segyio
import glob as glb

#FUNCTIONS Called in __main__ prior to passing to visualization classes********
def normalizeTraces(data):
    """
    This function normalizes each trace (column of 2d array) to the maximum
    value. This is a common way to visualize seismic, espeically first arrival
    travel-time data.
    
    INPUTS:
    data = a numpy array that is nt x ns (nt = time samples, ns = number of recievers)
    
    OUTPUTS:
    nData = a numpy array of the same size of input with traces normalized
    """
    nData = 0*data
    for i in range(0,data.shape[1]):
        nData[:,i] = data[:,i]/np.max(np.abs(data[:,i]))
    return nData  

def getFileInfo(dirName):
    """
    This function will read all of the *.segy or & *.sgy files in a given 
    directory. It returns a list with the file name and the shot location.
    This information will be passed to the GUI to display the file names. At
    a latter time it might be worth extracting other things from the headers
    and storing them in this list.
    
    DEPENDENCIES:
        GLOB - this is used to get the file names in the directory
        segyio - this is used to read the segy files and extract header info
    INPUTS:
        dirName (str) = this is a string to the directory that contains all of 
        the segy files from the survey.
    OUTPUTS:
        fileInfo is a list that is total Files by 2.
        Column 1 (str) = file name
        Column 2 (float) = shot location (units assumed to be m)
        
    NOTES:
        At this stage I use two if statemetns to check for segy files. If there
        are no segy files fileInfo will be an empty list and the user will get 
        an error. Though I am not sure where error goes in a GUI?
         - It depends, but we will be able to use try-except blocks for them
        
        It might be worth adding columns to this list if we need more info from
        the files later on
    """
    files = glb.glob(os.path.join(dirName, '*.sgy'))
    if files == []:
        files = glb.glob(os.path.join(dirName, '*.segy'))
        
    if files == []:
        print('No files with *.sgy or *.segy exist in this directory')
    #Column 1: File Name (str)
    #Column 2: SX (float)
    fileInfo = []
    
    for file in files:
        filename = os.path.basename(file)
        #print(filename)
        with segyio.open(file,strict = False) as f:
            shotLoc = f.header[0][segyio.TraceField.SourceX]
            #print(shotLoc)
        fileInfo.append([filename,shotLoc])
    return fileInfo


def getData(fileType, file):
    """
    Read data from segy or su file written to read a single file right now. 
    Could modify to extract shot location from a compiled file (or give file 
    list??) Options but segyio made it pretty easy.
    
    INPUTS
    File type = Str with either segy or su
    file = str with file name with path
    
    OUTPUTS
    x = 1D array with reciever locations in m
    t = 1D array with the time values in s
    data = trace data in an np array that is nt x ns
    gx = reciever spacing (calcualted from header) in m
    shotLoc = Shot Location in m
    """
    
    if str(fileType).lower() == 'segy':
        with segyio.open(file,strict = False) as f:
            t = f.samples/1000
            x = f.attributes(segyio.TraceField.GroupX)[:]
            shotLoc = f.header[0][segyio.TraceField.SourceX]
            gx = np.diff(x)[0]
            ngx = len(x)
            data = np.zeros((len(t),ngx))
            for i in range(0,ngx):
                data[:,i] = f.trace[i]
                
    elif str(fileType).lower() == 'su':
        with segyio.su.open(file) as f:
            t = f.samples/1000
            x = f.attributes(segyio.TraceField.GroupX)[:]
            shotLoc = f.header[0][segyio.TraceField.SourceX]
            gx = np.diff(x)[0]            
            ngx = len(x)
            data = np.zeros((len(t),ngx))
            for i in range(0,ngx):
                data[:,i] = f.trace[i]        
    return x, t, data, gx, shotLoc

def bpData(data,lf,hf,nq,order):
    """
    Applies a band-pass filter to each trace (column in 2d array)
    Inputs
    data = a numpy array that is nt x ns (nt = time samples, ns = number of recievers)
    lf = lower corner frequency (Hz)
    hf = upper corner frequency (Hz)
    nq = nyquist frequency (1/2*dt)
    order = order of the bp filter (required for sp.signal.butter)
    
    Outputs: 
    fData = a filtered (along columns) numpy array that is nt x ns (nt = time samples, ns = number of recievers)
    """
    wl = lf/nq
    wh = hf/nq
    b,a = signal.butter(order,[wl,wh],btype='bandpass')
    fData = data*0
    for i in range(0,data.shape[1]):
        fData[:,i] = signal.filtfilt(b,a,data[:,i])
        
    return fData

#END OF FUNCTIONS IN __main__**************************************************

#Main Window class must intialize first because it will check if pick file exists
class mainPickingWindow:

    def __init__(self,x,t,data,shotLoc,pickFileName):
        
        #LOCAL FUNCTIONS (NOT METHODS)****************************************
        def setUpFigLayout(initAmpSliderVal,initTimeSliderVal):
            """
            This will set up the main window layout.
            INPUTS:
                initAmpSliderVal = the initial value for the amplitude slider
                initTimeSliderVal = the intial value for the time slider
                
            OUTPUTS:
                fig1 = the main figure window (matplotLib Figure object)
                mainDataAxis = Main data axes (matplotLib axis object)
                ampSliderAxis = Amplitude slider for main data (matplotLib axis object)
                timeSliderAxis = Time slider for main data (matplotLib axis object)
                ampSlider = The amplitude "Slider" Object
                timeSlider = the time "Slider" object
            """
            fig1 = plt.figure(1,dpi=100,figsize=[8,7]) #Sizes hard-coded...
            gs = gridspec.GridSpec(5,1, height_ratios=[5,0.5,0.25,0.25,0.25])
            mainDataAxis = fig1.add_subplot(gs[0]) #Main data axes (matplotLib axis object)
            ampSliderAxis = fig1.add_subplot(gs[2]) #Amplitude slider for main data (matplotLib axis object)
            timeSliderAxis = fig1.add_subplot(gs[3]) #Time slider for main data (matplotLib axis object)
            
            ampSlider = Slider(ampSliderAxis, 'Amplitude', 0, 1, valinit=initAmpSliderVal, valstep=0.01)
            timeSlider = Slider(timeSliderAxis, 'Max Time',0,1,valinit=initTimeSliderVal,valstep=0.05)
            
            mainDataAxis.set_ylim([0,initAmpSliderVal])
            mainDataAxis.set_ylim([0,initTimeSliderVal])
            mainDataAxis.invert_yaxis()
            mainDataAxis.set_xlabel('Channel')
            mainDataAxis.set_ylabel('Time (s)')
            return fig1,mainDataAxis,ampSliderAxis,timeSliderAxis,ampSlider,timeSlider
        
        def updateFigure(updateFloat):
            #According to documentation "The function must accept a single float as its arguments."
            cSliderVal = self.ampSlider.val
            cTimeVal = self.timeSlider.val
            mainDataAxis.clear()
            mainDataAxis.pcolorfast(np.append(x,x[-1]+gx),np.append(t,t[-1]+dt),data,vmin=-cSliderVal,vmax=cSliderVal ,cmap='gray')
            if not (self.shotLocs == []):
               indShots = np.where(self.shotLocs==shotLoc)
               mainDataAxis.scatter(self.xPicks[indShots],self.tPicks[indShots],marker=1,s=50,c='c')
            mainDataAxis.set_ylim([0,cTimeVal])
            mainDataAxis.invert_yaxis()
            mainDataAxis.set_xlabel('Distance (m)')
            mainDataAxis.set_ylabel('Time (s)')
            plt.draw()
        
        #END LCOAL FUNCTIONS***************************************************        
        
        """
        Variables that need to be accesible, I am calling these properties.
        These will need to have the term self in front of the name
        shotLocs = list of shot locations
        xPicks = list of x-picks at each shot location
        tPicks = list of t-picks at each shot location        
        
        Functions that need to be accesible, these are methods:
        """
        # OBJECT PROPERTIES****************************************************
        #Initialze or read pick data. These are accesible otuside the class.
        if os.path.exists(pickFileName):
            tempData = np.loadtxt(pickFileName)
            self.shotLocs = tempData[:,0]
            self.xPicks = tempData[:,1]
            self.tPicks = tempData[:,2]
        else:
            self.shotLocs = []
            self.xPicks = []
            self.tPicks = []
        #END OBJECT PROPERTIES*************************************************
        
        #LOCAL VARIABLES*******************************************************
        #Calculate dt and gx (gx = geophone spacing in m)
        dt = np.round(np.diff(t)[0],decimals=4)
        gx = np.round(np.diff(x)[0],decimals=1)
        
        #Intial values for sliders
        initAmp4Slider = 0.5
        initTime4Slider = 0.75
        
        #END LOCAL VARIABLES***************************************************

        #Set up the figure
        fig1,mainDataAxis,ampSliderAxis,timeSliderAxis,self.ampSlider,self.timeSlider = setUpFigLayout(initAmp4Slider,initTime4Slider)
        
        #Initialization of first plot
        mainDataAxis.pcolorfast(np.append(x,x[-1]+gx),np.append(t,t[-1]+dt),data,vmin=-initAmp4Slider,vmax=initAmp4Slider ,cmap='gray')
        if not (self.shotLocs == []):
           indShots = np.where(self.shotLocs==shotLoc)
           mainDataAxis.scatter(self.xPicks[indShots],self.tPicks[indShots],marker=1,s=50,c='c')
        
        #ACTIVATE SLIDERS
        #I had to make these properties so that they would self update....
        #I wanted to keep them local but it didn't work.
        self.ampSlider.on_changed(updateFigure)
        self.timeSlider.on_changed(updateFigure)
        
        plt.show() #This has to be the last command.

        
        
##Trace Window class
#class tracePickingWindow:
#    def __init__(self,x,t,data,shotLoc,pickFileName)

if __name__ == '__main__':

    applyBPFilt = True
    #applyBPFilt = False
    lf = 10
    hf = 120
    nq = 500 #Nyquist Frequency
    order = 4    
    pickFile = 'test.txt'
    suFile = 'dataFiles/shot_01068m_stacked.su'
    
    #* HARD CODED FOR NOW. We will need to click and open a browser to select
    #This Path
    dirName = r"C:\Users\Fli034\Documents\firstArrivalPicker\eRFP_Develop\first-arrival-picker\dataFiles\survey1"
    shotLoc = 918
    fileInfo = getFileInfo(dirName)
    
    #Hard coded logic to search for shot value (shoult this be exact or appox?)
    tmpShotLocs = np.zeros((len(fileInfo),1))
    for k in range(0,len(fileInfo)):
        tmpShotLocs[k] = fileInfo[k][1]
    ind = np.argmin(((tmpShotLocs-shotLoc)**2)**0.5)
    #******************************************************
    
    [x,t,data,gx,shotLoc] = getData('segy', os.path.join(dirName, fileInfo[ind][0]))
    if applyBPFilt:
        data = bpData(data,lf,hf,nq,order)
    data = normalizeTraces(data)
    
    a = mainPickingWindow(x,t,data,shotLoc,pickFile)

