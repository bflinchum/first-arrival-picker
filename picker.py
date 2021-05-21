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

# FUNCTIONS Called in __main__ prior to passing to visualization classes********
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
    nData = 0 * data
    for i in range(0, data.shape[1]):
        nData[:, i] = data[:, i] / np.max(np.abs(data[:, i]))
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
    files = glb.glob(os.path.join(dirName, "*.sgy"))
    if files == []:
        files = glb.glob(os.path.join(dirName, "*.segy"))

    if files == []:
        print("No files with *.sgy or *.segy exist in this directory")
    # Column 1: File Name (str)
    # Column 2: SX (float)
    fileInfo = []

    for file in files:
        filename = os.path.basename(file)
        # print(filename)
        with segyio.open(file, strict=False) as f:
            shotLoc = f.header[0][segyio.TraceField.SourceX]
            # print(shotLoc)
        fileInfo.append([filename, shotLoc])
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
    if str(fileType).lower() == "segy":
        with segyio.open(file, strict=False) as f:
            t = f.samples / 1000
            x = f.attributes(segyio.TraceField.GroupX)[:]
            shotLoc = f.header[0][segyio.TraceField.SourceX]
            gx = np.diff(x)[0]
            ngx = len(x)
            data = np.zeros((len(t), ngx))
            for i in range(0, ngx):
                data[:, i] = f.trace[i]

    elif str(fileType).lower() == "su":
        with segyio.su.open(file) as f:
            t = f.samples / 1000
            x = f.attributes(segyio.TraceField.GroupX)[:]
            shotLoc = f.header[0][segyio.TraceField.SourceX]
            gx = np.diff(x)[0]
            ngx = len(x)
            data = np.zeros((len(t), ngx))
            for i in range(0, ngx):
                data[:, i] = f.trace[i]
    return x, t, data, gx, shotLoc


def bpData(data, lf, hf, nq, order):
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
    wl = lf / nq
    wh = hf / nq
    b, a = signal.butter(order, [wl, wh], btype="bandpass")
    fData = data * 0
    for i in range(0, data.shape[1]):
        fData[:, i] = signal.filtfilt(b, a, data[:, i])

    return fData


# END OF FUNCTIONS IN __main__**************************************************


class PickingWindow:
    """ Parent window class. Contains all common functionality """

    def __init__(self, x, t, data, figure, initAmp=0.5, initTime=0):
        """
        Constructor for PickingMethod
        Variables that need to be accesible, I am calling these properties.
        These will need to have the term self in front of the name
        shotLocs = list of shot locations
        xPicks = list of x-picks at each shot location
        tPicks = list of t-picks at each shot location
        
        Functions that need to be accesible, these are methods:
        """
        # Define attributes
        self.x = x
        self.t = t
        self.data = data
        self.figure = figure

        # Intial values for sliders
        self.initAmp = initAmp
        self.initTime = initTime

        self.xPicks = []
        self.tPicks = []

        self.xPicks2 = []
        self.tPicks2 = []
        
        # Set up the figure
        self.setUpFigLayout()

        # Plot Data
        self.plot_data()

        # ACTIVATE SLIDERS
        self.activate_sliders()

    def activate_sliders(self):
        self.ampSlider.on_changed(self.updateFigure)
        self.timeSlider.on_changed(self.updateFigure)

    def plot_data(self):
        pass

    def setAxisLimits(self):
        self.mainDataAxis.invert_yaxis()

    def setUpSliders(self):
        self.mainDataAxis = self.figure.subplots(1)

    def setUpFigLayout(self):
        """
        This will set up the main window layout.
        INPUTS:
            initAmp = the initial value for the amplitude slider
            initTime = the intial value for the time slider
            initWindowSize = the initial value for the window size (time)
            
        OUTPUTS:
            figure = the main figure window (matplotLib Figure object)
            mainDataAxis = Main data axes (matplotLib axis object)
            ampSliderAxis = Amplitude slider for main data (matplotLib axis object)
            timeSliderAxis = Time slider for main data (matplotLib axis object)
            windowSizeSliderAxis = Time slider for main data (matplotLib axis object)
            ampSlider = The amplitude "Slider" Object
            timeSlider = the time "Slider" object
            windowSizeSlider = the window size "Slider" object
        """
        self.setUpSliders()
        self.setAxisLimits()

    def updateFigure(self, updateFloat):
        """
        According to documentation:
        "The function must accept a single float as its arguments."
        """
        self.mainDataAxis.clear()
        self.plot_data()
        self.setAxisLimits()


class mainPickingWindow(PickingWindow):
    """ Class for the main window. Inherits from Picking Window """

    def __init__(self, x, t, data, figure):
        """ Constructor """
        # UNIQUE PROPERTIES*****************************************************
        # Calculate dt and gx (gx = geophone spacing in m)
        self.dt = np.round(np.diff(t)[0], decimals=4)
        self.gx = np.round(np.diff(x)[0], decimals=1)

        # Run parent class initialisation
        super().__init__(x, t, data, figure, initTime=0.75)

    def plot_data(self):
        """ Plot the data """
        self.mainDataAxis.pcolorfast(
            np.append(self.x, self.x[-1] + self.gx),
            np.append(self.t, self.t[-1] + self.dt),
            self.data,
            vmin=-self.ampSlider.val,
            vmax=self.ampSlider.val,
            cmap="gray",
        )
        self.mainDataAxis.scatter(self.xPicks, self.tPicks, marker=1, s=50, c="c")
        self.mainDataAxis.scatter(self.xPicks2, self.tPicks2, marker=1, s=50, c="r",alpha=0.5)

    def setAxisLimits(self):
        self.mainDataAxis.set_ylim([0, self.timeSlider.val])
        super().setAxisLimits()
        self.mainDataAxis.set_xlabel("Channel")
        self.mainDataAxis.set_ylabel("Time (s)")

    def setUpSliders(self):
        gs = gridspec.GridSpec(5, 1, height_ratios=[5, 0.5, 0.25, 0.25, 0.25])
        self.mainDataAxis = self.figure.add_subplot(
            gs[0]
        )  # Main data axes (matplotLib axis object)
        ampSliderAxis = self.figure.add_subplot(
            gs[2]
        )  # Amplitude slider for main data (matplotLib axis object)
        timeSliderAxis = self.figure.add_subplot(
            gs[3]
        )  # Time slider for main data (matplotLib axis object)

        self.ampSlider = Slider(
            ampSliderAxis, "Amplitude", 0, 1, valinit=self.initAmp, valstep=0.01
        )
        self.timeSlider = Slider(
            timeSliderAxis, "Max Time", 0, 1, valinit=self.initTime, valstep=0.05
        )


class tracePickingWindow(PickingWindow):
    """ Class for the trace picking window. Inherits from Picking Window """

    def __init__(self, x, t, data, figure, shotLoc, traceNum):
        """ Constructor """
        # UNIQUE PROPERTIES*****************************************************
        self.traceNum = traceNum
        self.initWindowSize = 0.1  # Intial values for window size slider

        super().__init__(x, t, data, figure)

    def activate_sliders(self):
        super().activate_sliders()
        self.windowSizeSlider.on_changed(self.updateFigure)

    @property
    def traceData(self):
        """ Dynamically updates its value when traceNum is changed """
        return self.data[:, self.traceNum]

    def plot_data(self):
        # Initialization of first plot
        self.mainDataAxis.plot(self.traceData, self.t, "k")
        self.mainDataAxis.fill_betweenx(
            self.t,
            0,
            self.traceData,
            where=self.traceData < 0,
            color="k",
            interpolate=True,
        )   
            
        self.mainDataAxis.scatter(self.xPicks, self.tPicks, marker="_", s=200, c="r")

        

    def setAxisLimits(self):
        self.mainDataAxis.set_ylim(
            [self.timeSlider.val, self.timeSlider.val + self.windowSizeSlider.val]
        )
        self.mainDataAxis.set_xlim([-self.ampSlider.val, self.ampSlider.val])
        super().setAxisLimits()
        self.mainDataAxis.set_xlabel("Distance (m)")
        self.mainDataAxis.set_ylabel("Time (s)")

    def setUpSliders(self):
        gs = gridspec.GridSpec(5, 1, height_ratios=[5, 0.5, 0.25, 0.25, 0.25])
        self.mainDataAxis = self.figure.add_subplot(
            gs[0]
        )  # Main data axes (matplotLib axis object)
        ampSliderAxis = self.figure.add_subplot(
            gs[2]
        )  # Amplitude slider for main data (matplotLib axis object)
        timeSliderAxis = self.figure.add_subplot(
            gs[3]
        )  # Time slider for main data (matplotLib axis object)
        windowSizeSliderAxis = self.figure.add_subplot(
            gs[4]
        )  # Time slider for main data (matplotLib axis object)

        self.ampSlider = Slider(
            ampSliderAxis, "Amplitude", 0, 1, valinit=self.initAmp, valstep=0.001
        )
        self.timeSlider = Slider(
            timeSliderAxis,
            "Initial Time",
            -0.1,
            0.3,
            valinit=self.initTime,
            valstep=0.001,
        )
        self.windowSizeSlider = Slider(
            windowSizeSliderAxis,
            "Window Size",
            0,
            0.5,
            valinit=self.initWindowSize,
            valstep=0.001,
        )


class doPicks:
    def __init__(self, x, t, data, shotLoc, initTraceNumb, pickFileName,convSPickFile):
        # Define object-level attributes
        self.x = x
        self.t = t
        self.data = data
        self.shotLoc = shotLoc
        self.pickFileName = pickFileName
        self.convSPickFile = convSPickFile

        # Create both windows with initalized values
        self.tracePickWindow = tracePickingWindow(
            self.x,
            self.t,
            self.data,
            plt.figure(2, dpi=100, figsize=[4, 7]),  # Sizes hard-coded...
            self.shotLoc,
            initTraceNumb,
        )
        self.mainWindowObject = mainPickingWindow(
            self.x,
            self.t,
            self.data,
            plt.figure(1, dpi=100, figsize=[8, 7]),  # Sizes hard-coded...
        )

        # Create an attribute to keep track of current trace (this will need travel throughout the class)
        self.cTrace = initTraceNumb
        # Initialize picking objects

        # Local variables to help me calcualte trace index
        self.gx = np.diff(self.x)[0]
        self.x0 = self.x[0]

        # Plot picks if they exists
        self.updatePicksMainWindow(self.shotLoc, self.pickFileName,self.convSPickFile)
        self.updatePicksTraceWindow(self.shotLoc, self.pickFileName)
        

    def connect(self):
        # Connect the TracePickWindow
        self.tracePickWindow.figure.canvas.mpl_connect(
            "button_press_event", self.whenClickedTraceWindow
        )
        self.tracePickWindow.figure.canvas.mpl_connect(
            "button_release_event", self.whenReleasedTraceWindow
        )

        # This (arrow presses) does not work in the gui
        self.tracePickWindow.figure.canvas.mpl_connect(
            "key_press_event", self.switchTraces
        )

        # Connect the MainWindow
        self.mainWindowObject.figure.canvas.mpl_connect(
            "button_press_event", self.whenClickedMainWindow
        )
        self.mainWindowObject.figure.canvas.mpl_connect(
            "button_release_event", self.getTraceMainWindow
        )

    def whenClickedTraceWindow(self, event):
        self.tracePickWindow.mainDataAxis.time_onclick = time.time()

    def whenReleasedTraceWindow(self, event):
        MAX_CLICK_LENGTH = 0.2
        if event.inaxes == self.tracePickWindow.mainDataAxis:
            if (
                event.button == 1
                and (time.time() - self.tracePickWindow.mainDataAxis.time_onclick)
                < MAX_CLICK_LENGTH
            ):
                tPick = event.ydata
                self.writePick(self.shotLoc, tPick, self.pickFileName)
                self.updatePicksMainWindow(self.shotLoc, self.pickFileName,self.convSPickFile)
                self.updatePicksTraceWindow(self.shotLoc, self.pickFileName)
                # print(tPick)
            if (
                event.button == 3
                and (time.time() - self.tracePickWindow.mainDataAxis.time_onclick)
                < MAX_CLICK_LENGTH
            ):
                tPick = event.ydata
                self.deletePick(self.shotLoc, tPick, self.pickFileName)
                self.updatePicksMainWindow(self.shotLoc, self.pickFileName,self.convSPickFile)
                self.updatePicksTraceWindow(self.shotLoc, self.pickFileName)

    def switchTraces(self, event):
        if event.key == "right":
            cSliderVal = self.tracePickWindow.ampSlider.val
            cTimeVal = self.tracePickWindow.timeSlider.val
            cWindowSize = self.tracePickWindow.windowSizeSlider.val
            self.cTrace = self.cTrace + 1
            self.tracePickWindow.traceNum = self.cTrace
            self.updatePicksTraceWindow(self.shotLoc, self.pickFileName)
        elif event.key == "left":
            cSliderVal = self.tracePickWindow.ampSlider.val
            cTimeVal = self.tracePickWindow.timeSlider.val
            cWindowSize = self.tracePickWindow.windowSizeSlider.val
            self.cTrace = self.cTrace - 1
            self.tracePickWindow.traceNum = (
                self.cTrace
            )  # Last minute modification...Update this attribute so sliders work better.
            self.updatePicksTraceWindow(self.shotLoc, self.pickFileName)

    def whenClickedMainWindow(self, event):
        self.mainWindowObject.mainDataAxis.time_onclick = time.time()

    def getTraceMainWindow(self, event):
        MAX_CLICK_LENGTH = 0.2
        if event.inaxes == self.mainWindowObject.mainDataAxis:
            if (
                event.button == 1
                and (time.time() - self.mainWindowObject.mainDataAxis.time_onclick)
                < MAX_CLICK_LENGTH
            ):
                xPick = event.xdata
                traceNum = np.round((xPick - self.x0) / self.gx)
                if traceNum < 0:
                    traceNum = 0
                elif traceNum > len(self.x):
                    traceNum = len(self.x)

                cSliderVal = self.tracePickWindow.ampSlider.val
                cTimeVal = self.tracePickWindow.timeSlider.val
                cWindowSize = self.tracePickWindow.windowSizeSlider.val
                self.cTrace = int(traceNum)
                self.tracePickWindow.traceNum = self.cTrace
                traceData = self.data[:, self.cTrace]
                self.tracePickWindow.mainDataAxis.clear()
                self.tracePickWindow.mainDataAxis.plot(traceData, self.t, "k")
                self.tracePickWindow.mainDataAxis.fill_betweenx(
                    self.t,
                    0,
                    traceData,
                    where=traceData < 0,
                    color="k",
                    interpolate=True,
                )
                self.tracePickWindow.mainDataAxis.set_ylim(
                    [cTimeVal, cTimeVal + cWindowSize]
                )
                self.tracePickWindow.mainDataAxis.set_xlim([-cSliderVal, cSliderVal])
                self.tracePickWindow.mainDataAxis.invert_yaxis()
                self.tracePickWindow.mainDataAxis.set_xlabel("Distance (m)")
                self.tracePickWindow.mainDataAxis.set_ylabel("Time (s)")
                self.tracePickWindow.figure.canvas.draw()

    def findIndRepeat(self, xPicks, shotLocs, cxPick, cxShot):
        # xPicks = array of xPicks
        # shotLocs = array of corresponding shot locs
        # cxPick = current x pick (single value)
        # cxShot = current shot location (single value)
        indRepeat = -999  # Initialize
        for kk in range(0, np.size(xPicks, 0)):
            if (shotLocs[kk] == cxShot) and (xPicks[kk] == cxPick):
                indRepeat = kk
        return indRepeat

    def deletePick(self, c_shotLoc, c_tPick, pickFile):
        c_xPick = self.x[self.cTrace]
        if os.path.exists(pickFile):
            tempData = np.loadtxt(pickFile)
            shotLocs = tempData[:, 0]
            xPicks = tempData[:, 1]
            tPicks = tempData[:, 2]
            indRepeat = self.findIndRepeat(xPicks, shotLocs, c_xPick, c_shotLoc)
            if indRepeat != -999:  # In otherwords no repeat
                print("Pick Deleted...")
                shotLocs = np.delete(shotLocs, indRepeat)
                xPicks = np.delete(xPicks, indRepeat)
                tPicks = np.delete(tPicks, indRepeat)
                tempArr = np.column_stack((shotLocs, xPicks, tPicks))
                np.savetxt(pickFile, tempArr, fmt="%10.5f")

    def writePick(self, c_shotLoc, c_tPick, pickFile):
        c_xPick = self.x[self.cTrace]

        if os.path.exists(pickFile):
            tempData = np.loadtxt(pickFile)
            shotLocs = tempData[:, 0]
            xPicks = tempData[:, 1]
            tPicks = tempData[:, 2]
            indRepeat = self.findIndRepeat(xPicks, shotLocs, c_xPick, c_shotLoc)
            if indRepeat == -999:  # In otherwords no repeat
                shotLocs = np.append(shotLocs, c_shotLoc)
                xPicks = np.append(xPicks, c_xPick)
                tPicks = np.append(tPicks, c_tPick)
                tempArr = np.column_stack((shotLocs, xPicks, tPicks))
                lexInd = np.lexsort((tempArr[:, 1], tempArr[:, 0]))
                tempArr = tempArr[lexInd]
                np.savetxt(pickFile, tempArr, fmt="%10.5f")
            else:
                shotLocs[indRepeat] = c_shotLoc
                xPicks[indRepeat] = c_xPick
                tPicks[indRepeat] = c_tPick
                tempArr = np.column_stack((shotLocs, xPicks, tPicks))
                lexInd = np.lexsort((tempArr[:, 1], tempArr[:, 0]))
                tempArr = tempArr[lexInd]
                np.savetxt(pickFile, tempArr, fmt="%10.5f")
        else:
            tempArr = [c_shotLoc, c_xPick, c_tPick]
            np.savetxt(pickFile, tempArr, fmt="%10.5f")

    def updatePicksMainWindow(self, shotLoc, pickFile,convSPickFile):
        if os.path.exists(pickFile):
            tempData = np.loadtxt(pickFile)
            shotLocs = tempData[:, 0]
            xPicks = tempData[:, 1]
            tPicks = tempData[:, 2]
            indShots = np.where(shotLocs == shotLoc)
            xPicks = xPicks[indShots]
            tPicks = tPicks[indShots]
            
        else:
            xPicks2 = []
            tPicks2 = []
            
            
        if os.path.exists(convSPickFile):
            tempData = np.loadtxt(convSPickFile)
            shotLocs = tempData[:, 0]
            xPicks2 = tempData[:, 1]
            tPicks2 = tempData[:, 2]
            indShots = np.where(shotLocs == shotLoc)
            xPicks2 = xPicks2[indShots]
            tPicks2 = tPicks2[indShots]        
        else:
            xPicks2 = []
            tPicks2 = []
        
        self.mainWindowObject.tPicks = tPicks
        self.mainWindowObject.xPicks = xPicks
        self.mainWindowObject.tPicks2 = tPicks2
        self.mainWindowObject.xPicks2 = xPicks2
        cSliderVal = self.mainWindowObject.ampSlider.val
        cTimeVal = self.mainWindowObject.timeSlider.val
        self.mainWindowObject.mainDataAxis.clear()
        self.mainWindowObject.mainDataAxis.pcolorfast(
            np.append(self.x, self.x[-1] + self.mainWindowObject.gx),
            np.append(self.t, self.t[-1] + self.mainWindowObject.dt),
            self.data,
            vmin=-cSliderVal,
            vmax=cSliderVal,
            cmap="gray",
        )
        self.mainWindowObject.mainDataAxis.scatter(
            xPicks2, tPicks2, marker=1, s=50, c='r',alpha=0.5
        )
        self.mainWindowObject.mainDataAxis.scatter(
            xPicks, tPicks, marker=1, s=50, c='c'
        )

        #print("made it here")
        #            if not (self.shotLocs == []):
        #               indShots = np.where(self.shotLocs==shotLoc)
        #               self.mainDataAxis.scatter(self.xPicks[indShots],self.tPicks[indShots],marker=1,s=50,c='c')
        self.mainWindowObject.mainDataAxis.set_ylim([0, cTimeVal])
        self.mainWindowObject.mainDataAxis.invert_yaxis()
        self.mainWindowObject.mainDataAxis.set_xlabel("Distance (m)")
        self.mainWindowObject.mainDataAxis.set_ylabel("Time (s)")
        self.mainWindowObject.figure.canvas.draw()

    def updatePicksTraceWindow(self, shotLoc, pickFile):
        if os.path.exists(pickFile):
            tempData = np.loadtxt(pickFile)
            shotLocs = tempData[:, 0]
            xPicks = tempData[:, 1]
            tPicks = tempData[:, 2]
            
            #Find Recipricals
            indRepeat = self.findIndRepeat(xPicks,shotLocs,shotLoc,self.x[self.cTrace])
            if not (indRepeat == -999):
                tPickRecip = tPicks[indRepeat]
            #print(indRepeat)
            
            indShots = np.where(shotLocs == shotLoc)
            xPicks = xPicks[indShots]
            tPicks = tPicks[indShots]
            indTrace = np.where(xPicks == self.x[self.cTrace])
            xPicks = xPicks[indTrace] * 0
            tPicks = tPicks[indTrace]
            

        else:
            xPicks = []
            tPicks = []
            
        self.tracePickWindow.tPicks = tPicks
        self.tracePickWindow.xPicks = xPicks
        cSliderVal = self.tracePickWindow.ampSlider.val
        cTimeVal = self.tracePickWindow.timeSlider.val
        cWindowSize = self.tracePickWindow.windowSizeSlider.val
        self.tracePickWindow.traceNum = self.cTrace
        traceData = self.data[:, self.cTrace]
        self.tracePickWindow.mainDataAxis.clear()
        self.tracePickWindow.mainDataAxis.plot(traceData, self.t, "k")
        self.tracePickWindow.mainDataAxis.fill_betweenx(
            self.t, 0, traceData, where=traceData < 0, color="k", interpolate=True
        )
        self.tracePickWindow.mainDataAxis.scatter(
            xPicks, tPicks, marker="_", s=200, c='c'
        )
        
        if not (indRepeat == -999):
            self.tracePickWindow.mainDataAxis.scatter(0,tPickRecip,marker='P',s=100,c='m') 
         
                
        self.tracePickWindow.mainDataAxis.set_ylim([cTimeVal, cTimeVal + cWindowSize])
        self.tracePickWindow.mainDataAxis.set_xlim([-cSliderVal, cSliderVal])
        self.tracePickWindow.mainDataAxis.invert_yaxis()
        self.tracePickWindow.mainDataAxis.set_xlabel("Distance (m)")
        self.tracePickWindow.mainDataAxis.set_ylabel("Time (s)")
        self.tracePickWindow.figure.canvas.draw()


class picker:
    def __init__(self):
        applyBPFilt = True
        applyBPFilt = False
        lf = 10
        hf = 150
        nq = 2000  # Nyquist Frequency
        order = 4
        pickFile = "CEF_2021_0315_Swave_picks_after_allShots.txt"
        convSPickFile = '/Users/bflinch/Dropbox/Clemson/Research/ResearchProjects/SwaveRefraction/Time-Lapse_CEF/ExtractMatchingPicks/After_allShots/CEF_2021_0315_Swave_picks_after_allShots_Swave_modeledTravelTimes.txt'
        #pickFile = '/Users/bflinch/Desktop/test.txt'
        # Add function to check for duplicates

        # * HARD CODED FOR NOW. We will need to click and open a browser to select
        # This Path
        #dirName = "/Users/bflinch/Dropbox/Clemson/Research/ResearchProjects/SwaveRefraction/BW_072613/StackFiles/segyFiles/stacked"
        #dirName = "/Users/bflinch/Dropbox/Clemson/Research/ResearchProjects/SwaveRefraction/BW_072613/Vertical/segyFiles"
        #dirName = "/Users/bflinch/Dropbox/Clemson/Teaching/Spring2021/GEOL3700_CriticalZoneGeophysics/Data/Day1/Seismic/segyFiles/P-Wave"
        dirName = "/Users/bflinch/Dropbox/Clemson/Teaching/Spring2021/GEOL3700_CriticalZoneGeophysics/Data/Day2/Seismic/segyFiles/S-Wave/stacked"
        #dirName = '/Users/bflinch/Dropbox/Clemson/Research/ResearchProjects/FlatwoodSchool_terry/Seismic/P-waves/stacked/'    
        #dirName = '/Users/bflinch/Dropbox/Clemson/Research/ResearchProjects/FlatwoodSchool_terry/Seismic/S-waves/stacked/'
        shotLoc = 100
        fileInfo = getFileInfo(dirName)
        print(fileInfo)
        
        # Hard coded logic to search for shot value (shoult this be exact or appox?)
        tmpShotLocs = np.zeros((len(fileInfo), 1))
        for k in range(0, len(fileInfo)):
            tmpShotLocs[k] = fileInfo[k][1]
        ind = np.argmin(((tmpShotLocs - shotLoc) ** 2) ** 0.5)
        #ind = 7
        # ******************************************************

        [x, t, data, gx, shotLoc] = getData(
            "segy", os.path.join(dirName, fileInfo[ind][0])
        )
        #x = np.linspace(0,237.5,96)
        if applyBPFilt:
            data = bpData(data, lf, hf, nq, order)
        data = normalizeTraces(data)

        self.c = doPicks(x, t, data, shotLoc, 10, pickFile,convSPickFile)
        # tracePickingWindow(
        #     x,
        #     t,
        #     data,
        #     plt.figure(2, dpi=100, figsize=[4, 7]),  # Sizes hard-coded...
        #     shotLoc,
        #     10,
        # )


if __name__ == "__main__":
    a = picker()
    a.c.connect()
    plt.show()
