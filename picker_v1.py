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

#Function that normalizes each trace to maximum amplitude for plotting
def normalizeTraces(data):
    #INPUTS:
    #data = a numpy array that is nt x ns (nt = time samples, ns = number of recievers)
    
    #OUTPUTS:
    #nData = a numpy array of the same size of input with traces normalized
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
    files = glb.glob(dirName + '\*.sgy')
    if files == []:
        files = glb.glob(dirName + '\*.segy')
        
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

#Read data from segy or su file
#written to read a single file right now. Could modify to extract shot location
#From a compiled file (or give file list??) Options but segyio made it pretty 
#easy :D
def getData(fileType, file):
    #INPUTS
    #File type = Str with either segy or su
    #file = str with file name with path
    
    #OUTPUTS
    #x = 1D array with reciever locations in m
    #t = 1D array with the time values in s
    #data = trace data in an np array that is nt x ns
    #gx = reciever spacing (calcualted from header) in m
    #shotLoc = Shot Location in m
    
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

#Function Name: bpData
def bpData(data,lf,hf,nq,order):
    #Inputs
    #data = a numpy array that is nt x ns (nt = time samples, ns = number of recievers)
    #lf = lower corner frequency (Hz)
    #hf = upper corner frequency (Hz)
    #nq = nyquist frequency (1/2*dt)
    #order = order of the bp filter (required for sp.signal.butter)
    
    #Outputs: 
    #fData = a filtered (along columns) numpy array that is nt x ns (nt = time samples, ns = number of recievers)
    #Description normalizes traces by columns
    wl = lf/nq
    wh = hf/nq
    b,a = signal.butter(order,[wl,wh],btype='bandpass')
    fData = data*0
    for i in range(0,data.shape[1]):
        fData[:,i] = signal.filtfilt(b,a,data[:,i])
        
    return fData
#Not sure I could comment this as I just got it working mostly through trial and
#error. I am still not quite sure how or why it works.
class pickingModule:
    def __init__(self, x, t, data, shotLoc, pickFileName, gx):
        self.dt = np.round(np.diff(t)[0],decimals=4)
        
        #DEFINE VARIABLES
        if os.path.exists(pickFileName):
            tempData = np.loadtxt(pickFileName)
            self.shotLocs = tempData[:,0]
            self.xPicks = tempData[:,1]
            self.tPicks = tempData[:,2]
        else:
            self.shotLocs = []
            self.xPicks = []
            self.tPicks = []
            
        #Check if pick file exists
            #IF YeS
                #Read File
                #Populate shotLoc,xPicks,tPicks
            #ELSE
                #Initialize blank variables and create file
        
        #SET UP FIGURE SPACES
        self.fig1 = plt.figure(1,dpi=100,figsize=[8,7])
        self.fig2 = plt.figure(2,dpi=100,figsize=[4,7])

        gs = gridspec.GridSpec(5,1, height_ratios=[5,0.5,0.25,0.25,0.25])

        self.ax1 = self.fig1.add_subplot(gs[0]) #Main data axes
        self.ax2 = self.fig1.add_subplot(gs[2]) #Amplitude slider for main data
        self.ax3 = self.fig1.add_subplot(gs[3]) #Time slider for main data
        
        self.ax4 = self.fig2.add_subplot(gs[0]) #Trace window axes
        self.ax5 = self.fig2.add_subplot(gs[2]) #Amplitude slider for main data
        self.ax6 = self.fig2.add_subplot(gs[3]) #start time slider     
        self.ax7 = self.fig2.add_subplot(gs[4]) #start time slider
        
        #INITIALIZE VALUES FOR FIRST PLOT
        self.mainAmpSliderVal = 0.5
        self.mainTimeSliderVal = 0.75
        self.traceAmpSliderVal = 0.5
        self.traceTimeSliderVal = 0
        self.traceWindowSizeSliderVal = 0.1
        
        #CREATE SLIDER BARS
        self.mainAmpSlider = Slider(self.ax2, 'Amplitude', 0, 1, valinit=self.mainAmpSliderVal, valstep=0.01)
        self.mainTimeSlider = Slider(self.ax3, 'Max Time',0,1,valinit=self.mainTimeSliderVal,valstep=0.05)
        
        self.traceAmpSlider = Slider(self.ax5, 'Amplitude', 0, 0.8, valinit=self.traceAmpSliderVal, valstep=0.001)
        self.traceTimeStartSlider = Slider(self.ax6, 'Start Time',-0.01,0.1,valinit=self.traceTimeSliderVal,valstep=0.001)
        self.traceWindowSizeSlider = Slider(self.ax7, 'Window Size',0,0.5,valinit=self.traceWindowSizeSliderVal,valstep=0.001)

        #INITIALIZE WITH FIRST DATA
        #self.ax1.pcolorfast((np.min(x),np.max(x)+gx),(np.min(t),np.max(t)),data,vmin=-self.mainAmpSliderVal,vmax=self.mainAmpSliderVal ,cmap='gray')
        self.ax1.pcolorfast(np.append(x,x[-1]+gx),np.append(t,t[-1]+self.dt),data,vmin=-self.mainAmpSliderVal,vmax=self.mainAmpSliderVal ,cmap='gray')
                
        if not (self.shotLocs == []):
            indShots = np.where(self.shotLocs==shotLoc)
            self.ax1.scatter(self.xPicks[indShots],self.tPicks[indShots],marker=1,s=50,c='c')
        self.ax1.set_ylim([0,self.mainTimeSliderVal])
        self.ax1.invert_yaxis()
        self.ax1.set_xlabel('Channel')
        self.ax1.set_ylabel('Time (s)')
        
        #Extract first traces to plot
        #add some code to extract trace nearest the shot location
        #self.traceNum = 5
        self.traceNum = np.argmin(np.sqrt((x-shotLoc)**2))
        self.tData = data[:,self.traceNum]
        
        def findIndRepeat(xPicks,shotLocs,cxPick,cxShot):
            #xPicks = array of xPicks
            #shotLocs = array of corresponding shot locs
            #cxPick = current x pick (single value)
            #cxShot = current shot location (single value)
            indRepeat = -999 #Initialize
            for kk in range(0,np.size(xPicks,0)):
                if (shotLocs[kk] == cxShot) and (xPicks[kk] == cxPick):
                    indRepeat = kk                    
            return indRepeat
            
        #INITIALIZE WITH FIRST TRACE
        self.ax4.plot(self.tData,t,'k')
        self.ax4.fill_betweenx(t, 0, self.tData, where=self.tData<0, color='k', interpolate=True)
        indRepeat = findIndRepeat(self.xPicks,self.shotLocs,x[self.traceNum],shotLoc)
        if not (indRepeat == -999):
            self.ax4.scatter(0,self.tPicks[indRepeat],marker='_',s=200,c='c')
        #Check for recipricals    
        indRepeat = findIndRepeat(self.xPicks,self.shotLocs,shotLoc,x[self.traceNum])
        if not (indRepeat == -999):
            self.ax4.scatter(0,self.tPicks[indRepeat],marker='P',s=100,c='m')
            #self.ax1.scatter(x[self.traceNum],self.tPicks[indRepeat],marker=1,s=50,c='r')

        self.ax4.set_xlim([-self.traceAmpSliderVal,self.traceAmpSliderVal])
        self.ax4.set_ylim([self.traceTimeSliderVal,self.traceTimeSliderVal+self.traceWindowSizeSliderVal])
        self.ax4.invert_yaxis()
        self.ax4.set_xlabel('Amplitude')
        self.ax4.set_ylabel('Time (s)')
        self.ax4.set_title(str(x[self.traceNum]) + ' m')
        

        #Funciton that will control what the sliders on the trace window will do
        def updateAmp_trace(val):            
            self.traceAmpSliderVal = self.traceAmpSlider.val
            #print(amp)
            self.traceTimeSliderVal = self.traceTimeStartSlider.val
            self.traceWindowSizeSliderVal = self.traceWindowSizeSlider.val
            
            self.ax4.clear()
            self.ax4.plot(self.tData,t,'k')
            self.ax4.fill_betweenx(t, 0, self.tData, where=self.tData<0, color='k', interpolate=True)
            indRepeat = findIndRepeat(self.xPicks,self.shotLocs,x[self.traceNum],shotLoc)            
            if not (indRepeat == -999):
                self.ax4.scatter(0,self.tPicks[indRepeat],marker='_',s=200,c='c') 
            #Check for recipricals    
            indRepeat = findIndRepeat(self.xPicks,self.shotLocs,shotLoc,x[self.traceNum])
            if not (indRepeat == -999):
                self.ax4.scatter(0,self.tPicks[indRepeat],marker='P',s=100,c='m')                
            self.ax4.set_ylim([self.traceTimeSliderVal,self.traceTimeSliderVal+self.traceWindowSizeSliderVal])
            self.ax4.set_xlim([-self.traceAmpSliderVal,self.traceAmpSliderVal])
            self.ax4.invert_yaxis()
            self.ax4.set_xlabel('Amplitude')
            self.ax4.set_ylabel('Time (s)')
            self.ax4.set_title(str(x[self.traceNum]) + ' m')
        
        #Activate the sliders
        self.traceAmpSlider.on_changed(updateAmp_trace)
        self.traceTimeStartSlider.on_changed(updateAmp_trace)
        self.traceWindowSizeSlider.on_changed(updateAmp_trace)
        
        #Function that will control whtat the sliders on the main window will do
        def updateAmp_main(val):
            self.mainAmpSliderVal = self.mainAmpSlider.val
            #print(amp)
            self.mainTimeSliderVal = self.mainTimeSlider.val
            
            self.ax1.clear()
            self.ax1.pcolorfast(np.append(x,x[-1]+gx),np.append(t,t[-1]+self.dt),data,vmin=-self.mainAmpSliderVal,vmax=self.mainAmpSliderVal,cmap='gray')
            
            if not (self.shotLocs == []):
                indShots = np.where(self.shotLocs==shotLoc)
                self.ax1.scatter(self.xPicks[indShots],self.tPicks[indShots],marker=1,s=50,c='c')            
            self.ax1.set_ylim([0,self.mainTimeSliderVal])
            self.ax1.invert_yaxis()
            self.ax1.set_xlabel('Distance (m)')
            self.ax1.set_ylabel('Time (s)')
        
        #Activate the sliders    
        self.mainAmpSlider.on_changed(updateAmp_main)
        self.mainTimeSlider.on_changed(updateAmp_main)
        #OPEN/CREATE EVENT HANDLERS FOR PICKIGN AND DELETING
        #define funciton to move the traces
        def switchTraces(event):            
            #Add some error checking here so I don't get an error if I keep pressing
            # right or left
            if event.key=='right':
                self.traceNum = self.traceNum  + 1
                self.tData = data[:,self.traceNum]
                self.ax4.clear()
                self.ax4.plot(self.tData,t,'k')
                self.ax4.fill_betweenx(t, 0, self.tData, where=self.tData<0, color='k', interpolate=True)
                indRepeat = findIndRepeat(self.xPicks,self.shotLocs,x[self.traceNum],shotLoc)                
                if not (indRepeat == -999):
                    self.ax4.scatter(0,self.tPicks[indRepeat],marker='_',s=200,c='c') 
                #Check for recipricals    
                indRepeat = findIndRepeat(self.xPicks,self.shotLocs,shotLoc,x[self.traceNum])
                if not (indRepeat == -999):
                    self.ax4.scatter(0,self.tPicks[indRepeat],marker='P',s=100,c='m')
                self.ax4.set_ylim([self.traceTimeSliderVal,self.traceTimeSliderVal+self.traceWindowSizeSliderVal])
                self.ax4.set_xlim([-self.traceAmpSliderVal,self.traceAmpSliderVal])
                self.ax4.invert_yaxis()
                self.ax4.set_xlabel('Amplitude')
                self.ax4.set_ylabel('Time (s)')
                self.ax4.set_title(str(x[self.traceNum]) + ' m')
                self.fig2.canvas.draw()
                    
            elif event.key=='left':
                self.traceNum  = self.traceNum - 1
                self.tData = data[:,self.traceNum]
                self.ax4.clear()
                self.ax4.plot(self.tData,t,'k')
                self.ax4.fill_betweenx(t, 0, self.tData, where=self.tData<0, color='k', interpolate=True)
                indRepeat = findIndRepeat(self.xPicks,self.shotLocs,x[self.traceNum],shotLoc)                
                if not (indRepeat == -999):
                    self.ax4.scatter(0,self.tPicks[indRepeat],marker='_',s=200,c='c')
                #Check for recipricals    
                indRepeat = findIndRepeat(self.xPicks,self.shotLocs,shotLoc,x[self.traceNum])
                if not (indRepeat == -999):
                    self.ax4.scatter(0,self.tPicks[indRepeat],marker='P',s=100,c='m')
                self.ax4.set_ylim([self.traceTimeSliderVal,self.traceTimeSliderVal+self.traceWindowSizeSliderVal])
                self.ax4.set_xlim([-self.traceAmpSliderVal,self.traceAmpSliderVal])
                self.ax4.invert_yaxis()
                self.ax4.set_xlabel('Amplitude')
                self.ax4.set_ylabel('Time (s)')
                self.ax4.set_title(str(x[self.traceNum]) + ' m')
                self.fig2.canvas.draw()

        
            
        def whenReleased(event):
            MAX_CLICK_LENGTH = 0.2
            
            if event.inaxes == self.ax4: #Checks to make sure click is in the window
                #Check to make sure it was a left click and released wtihin 0.1 s (so zoom funciton stays usable)
                if event.button == 1 and (time.time() - self.ax4.time_onclick) < MAX_CLICK_LENGTH:
                    #Pull out the time of the pick
                    tPick = event.ydata
                    #Get amplitude of pick
                    ampValatPick = self.tData[np.argmin(np.sqrt((t-tPick)**2))]

                    #print(str(shotLoc) + ' ' + str(x[self.traceNum]) + ' ' + str(np.round(tPick,decimals=5)))
                    
                    #Iinitialize variables if no picks exist
                    if self.shotLocs == []:
                        self.shotLocs = np.array([shotLoc])
                        self.xPicks = np.array([x[self.traceNum]])
                        self.tPicks = np.array([tPick])
                    #if picks alread exist
                    else:
                        #This logic checks for an existing pick with sx=sx and gx=gx
                        indRepeat = findIndRepeat(self.xPicks,self.shotLocs,x[self.traceNum],shotLoc)
                        #If Pick doesn't exists then append new value        
                        if indRepeat == -999:
                            #print('made it here')
                            self.shotLocs = np.append(self.shotLocs,shotLoc)
                            self.xPicks = np.append(self.xPicks,x[self.traceNum])
                            self.tPicks = np.append(self.tPicks,tPick)                

                            tempArr = np.column_stack((self.shotLocs,self.xPicks,self.tPicks))
                            lexInd = np.lexsort((tempArr[:,1],tempArr[:,0]))
                            tempArr = tempArr[lexInd]
                            np.savetxt(pickFileName,tempArr,fmt='%10.5f')
                        
                        #If pick exists then overwrite the value    
                        else:
                            self.shotLocs[indRepeat] = shotLoc
                            self.xPicks[indRepeat] = x[self.traceNum]
                            self.tPicks[indRepeat] = tPick
                            tempArr = np.column_stack((self.shotLocs,self.xPicks,self.tPicks))
                            lexInd = np.lexsort((tempArr[:,1],tempArr[:,0]))
                            tempArr = tempArr[lexInd]
                            np.savetxt(pickFileName,tempArr,fmt='%10.5f')

                        
                #If it's a rigth click delete the pick
                if event.button == 3 and (time.time() - self.ax4.time_onclick) < MAX_CLICK_LENGTH:
                    indRepeat = findIndRepeat(self.xPicks,self.shotLocs,x[self.traceNum],shotLoc)
                    #If Pick doesn't exists then append new value        
                    if not (indRepeat == -999):
                        print('Pick Deleted...')
                        self.shotLocs = np.delete(self.shotLocs,indRepeat)
                        self.xPicks = np.delete(self.xPicks,indRepeat)
                        self.tPicks = np.delete(self.tPicks,indRepeat)
                        tempArr = np.column_stack((self.shotLocs,self.xPicks,self.tPicks))
                        np.savetxt(pickFileName,tempArr,fmt='%10.5f')                  
                
                indRepeat = findIndRepeat(self.xPicks,self.shotLocs,x[self.traceNum],shotLoc)
                if not (indRepeat == -999):
                    self.ax4.scatter(ampValatPick,self.tPicks[indRepeat],marker='_',s=200,c='r')
                    self.ax1.scatter(x[self.traceNum],self.tPicks[indRepeat],marker=1,s=50,c='r')

                    self.fig1.canvas.draw()
                    self.fig2.canvas.draw()
                

        def whenClicked(event):
            self.ax4.time_onclick = time.time()

                        
        self.fig2.canvas.mpl_connect('key_press_event',switchTraces)
        self.fig2.canvas.mpl_connect('button_press_event',whenClicked)
        self.fig2.canvas.mpl_connect('button_release_event',whenReleased)
        
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
        
    a = pickingModule(x,t,data,shotLoc,pickFile,gx)

    
#***************THESE ARE PLOTS THAT I WANT to Incorporate into the GUI*******
#    #They are built from the pick data read from a file during the GUI
#    
#    #Plot recipricals
#    def findIndRepeat(xPicks,shotLocs,cxPick,cxShot):
#        #xPicks = array of xPicks
#        #shotLocs = array of corresponding shot locs
#        #cxPick = current x pick (single value)
#        #cxShot = current shot location (single value)
#        indRepeat = -999 #Initialize
#        for kk in range(0,np.size(xPicks,0)):
#            if (shotLocs[kk] == cxShot) and (xPicks[kk] == cxPick):
#                indRepeat = kk                    
#        return indRepeat
#            
#    xPicks = a.xPicks
#    shotLocs = a.shotLocs
#    ttPicks = a.tPicks
#    plt.figure(9999)
#    plt.scatter(xPicks,shotLocs,c=ttPicks*1000,marker='s',s=20,vmin=10,vmax=150)
#    plt.xlabel('Geophone Location (m)')
#    plt.ylabel('Shot Location (m)')
#    plt.figure(99990)
#    plt.scatter(np.abs(xPicks-shotLocs),ttPicks*1000,marker='s',c=shotLocs,s=5)
#    plt.xlabel('Offset (m)')
#    plt.ylabel('Travel-time (ms)')
#    plt.grid()
#    plt.colorbar().set_label('Shot Location (m)')
#    index = 0
#    obs = np.ones((10000,1))*-999
#    recip = np.ones((10000,1))*-999
#    while len(xPicks) > 0:
#        indRepeat = findIndRepeat(xPicks,shotLocs,shotLocs[0],xPicks[0])
#        if not (indRepeat==-999):
#            obs[index] = ttPicks[0]
#            recip[index] = ttPicks[indRepeat]
#            index = index + 1
#        xPicks = np.delete(xPicks,0)
#        shotLocs = np.delete(shotLocs,0)
#        ttPicks = np.delete(ttPicks,0)
#    
#    ind2rmv = np.where(obs==-999)[0]    
#    obs = np.delete(obs,ind2rmv)
#    recip = np.delete(recip,ind2rmv)
#    
#    f = open('recips.xyz', 'w')
#    for k in range(0,len(obs)):
#        f.write(str(obs[k]*1000) + ' ' + str(recip[k]*1000) + '\n')
#    f.close()  
#    fig999 = plt.figure(999)
#    ax999 = fig999.add_subplot(1,1,1)
#    ax999.plot([0,200],[0,200],'r')
#    ax999.plot([0,200],[5,205],'b')
#    ax999.plot([0,200],[-5,195],'b')
#    ax999.plot([0,200],[10,210],'c')
#    ax999.plot([0,200],[-10,190],'c')
#    ax999.scatter(obs*1000,recip*1000,marker='+',s=25,c='k')
#    ax999.set_xlim([0, 180])
#    ax999.set_ylim([0, 180])
#    avgDiff = np.mean(np.abs(obs-recip))*1000
#    avgStd = np.std(np.abs(obs-recip))*1000
#    print(avgDiff)
#    print(avgStd)

