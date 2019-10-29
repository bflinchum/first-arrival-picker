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

# Main Window class must intialize first because it will check if pick file exists
class mainPickingWindow:
    def __init__(self, x, t, data, shotLoc):
        # LOCAL FUNCTIONS (NOT METHODS)****************************************
        # This has to be local because it can only take a single float as input
        # if I make it a method it requires "self" as an input--which means two
        # inputs and I get a bunch of errors.

        def updateFigure(updateFloat):
            # According to documentation "The function must accept a single float as its arguments."
            cSliderVal = self.ampSlider.val
            cTimeVal = self.timeSlider.val
            self.mainDataAxis.clear()
            self.mainDataAxis.pcolorfast(
                np.append(x, x[-1] + self.gx),
                np.append(t, t[-1] + self.dt),
                data,
                vmin=-cSliderVal,
                vmax=cSliderVal,
                cmap="gray",
            )
            self.mainDataAxis.scatter(self.xPicks, self.tPicks, marker=1, s=50, c="c")
            self.mainDataAxis.set_ylim([0, cTimeVal])
            self.mainDataAxis.invert_yaxis()
            self.mainDataAxis.set_xlabel("Distance (m)")
            self.mainDataAxis.set_ylabel("Time (s)")
            plt.draw()

        # END LCOAL FUNCTIONS***************************************************

        """
        Variables that need to be accesible, I am calling these properties.
        These will need to have the term self in front of the name
        shotLocs = list of shot locations
        xPicks = list of x-picks at each shot location
        tPicks = list of t-picks at each shot location        
        
        Functions that need to be accesible, these are methods:
        """

        # LOCAL VARIABLES*******************************************************
        # Calculate dt and gx (gx = geophone spacing in m)
        self.dt = np.round(np.diff(t)[0], decimals=4)
        self.gx = np.round(np.diff(x)[0], decimals=1)
        self.xPicks = []
        self.tPicks = []

        # Intial values for sliders
        initAmp4Slider = 0.5
        initTime4Slider = 0.75

        # END LOCAL VARIABLES***************************************************

        # Set up the figure
        self.figureObject, self.mainDataAxis, ampSliderAxis, timeSliderAxis, self.ampSlider, self.timeSlider = self.setUpFigLayout(
            initAmp4Slider, initTime4Slider
        )

        # Initialization of first plot
        self.mainDataAxis.pcolorfast(
            np.append(x, x[-1] + self.gx),
            np.append(t, t[-1] + self.dt),
            data,
            vmin=-initAmp4Slider,
            vmax=initAmp4Slider,
            cmap="gray",
        )
        #        if not (self.shotLocs == []):
        #           indShots = np.where(self.shotLocs==shotLoc)
        #           self.mainDataAxis.scatter(self.xPicks[indShots],self.tPicks[indShots],marker=1,s=50,c='c')

        # ACTIVATE SLIDERS
        # I had to make these properties so that they would self update....
        # I wanted to keep them local but it didn't work.
        self.ampSlider.on_changed(updateFigure)
        self.timeSlider.on_changed(updateFigure)

        # plt.show() #This has to be the last command.

    def setUpFigLayout(self, initAmpSliderVal, initTimeSliderVal):
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
        fig1 = plt.figure(1, dpi=100, figsize=[8, 7])  # Sizes hard-coded...
        gs = gridspec.GridSpec(5, 1, height_ratios=[5, 0.5, 0.25, 0.25, 0.25])
        mainDataAxis = fig1.add_subplot(
            gs[0]
        )  # Main data axes (matplotLib axis object)
        ampSliderAxis = fig1.add_subplot(
            gs[2]
        )  # Amplitude slider for main data (matplotLib axis object)
        timeSliderAxis = fig1.add_subplot(
            gs[3]
        )  # Time slider for main data (matplotLib axis object)

        ampSlider = Slider(
            ampSliderAxis, "Amplitude", 0, 1, valinit=initAmpSliderVal, valstep=0.01
        )
        timeSlider = Slider(
            timeSliderAxis, "Max Time", 0, 1, valinit=initTimeSliderVal, valstep=0.05
        )

        mainDataAxis.set_ylim([0, initAmpSliderVal])
        mainDataAxis.set_ylim([0, initTimeSliderVal])
        mainDataAxis.invert_yaxis()
        mainDataAxis.set_xlabel("Channel")
        mainDataAxis.set_ylabel("Time (s)")
        return fig1, mainDataAxis, ampSliderAxis, timeSliderAxis, ampSlider, timeSlider


class tracePickingWindow:
    def __init__(self, x, t, data, shotLoc, traceNum):
        # LOCAL FUNCTIONS (NOT METHODS)****************************************
        # I know this is unconvential but the matplotlib sliders require this.
        def updateFigure(updateFloat):
            # According to documentation "The function must accept a single float as its arguments."
            cSliderVal = self.ampSlider.val
            cTimeVal = self.timeSlider.val
            cWindowSize = self.windowSizeSlider.val
            traceData = data[:, self.traceNum]
            self.mainDataAxis.clear()
            self.mainDataAxis.plot(traceData, t, "k")
            self.mainDataAxis.fill_betweenx(
                t, 0, traceData, where=traceData < 0, color="k", interpolate=True
            )
            self.mainDataAxis.scatter(
                self.xPicks, self.tPicks, marker="_", s=200, c="r"
            )
            # mainDataAxis.pcolorfast(np.append(x,x[-1]+gx),np.append(t,t[-1]+dt),data,vmin=-cSliderVal,vmax=cSliderVal ,cmap='gray')
            #            if not (self.shotLocs == []):
            #               indShots = np.where(self.shotLocs==shotLoc)
            #               mainDataAxis.scatter(self.xPicks[indShots],self.tPicks[indShots],marker=1,s=50,c='c')
            self.mainDataAxis.set_ylim([cTimeVal, cTimeVal + cWindowSize])
            self.mainDataAxis.set_xlim([-cSliderVal, cSliderVal])
            self.mainDataAxis.invert_yaxis()
            self.mainDataAxis.set_xlabel("Distance (m)")
            self.mainDataAxis.set_ylabel("Time (s)")
            plt.draw()

        # END LCOAL FUNCTIONS***************************************************

        """
        Variables that need to be accesible, I am calling these properties.
        These will need to have the term self in front of the name
        shotLocs = list of shot locations
        xPicks = list of x-picks at each shot location
        tPicks = list of t-picks at each shot location        
        
        Functions that need to be accesible, these are methods:
        """
        self.traceNum = traceNum
        self.xPicks = []
        self.tPicks = []
        # LOCAL VARIABLES*******************************************************
        # Calculate dt and gx (gx = geophone spacing in m)
        traceData = data[:, self.traceNum]
        # Intial values for sliders
        initAmp4Slider = 0.5
        initTime4Slider = 0
        initWindowSize4Slider = 0.1
        # END LOCAL VARIABLES***************************************************

        self.figureObject, self.mainDataAxis, ampSliderAxis, timeSliderAxis, windowSizeSliderAxis, self.ampSlider, self.timeSlider, self.windowSizeSlider = self.setUpFigLayout(
            initAmp4Slider, initTime4Slider, initWindowSize4Slider
        )

        self.mainDataAxis.plot(traceData, t, "k")
        self.mainDataAxis.fill_betweenx(
            t, 0, traceData, where=traceData < 0, color="k", interpolate=True
        )

        self.ampSlider.on_changed(updateFigure)
        self.timeSlider.on_changed(updateFigure)
        self.windowSizeSlider.on_changed(updateFigure)

    def setUpFigLayout(self, initAmpSliderVal, initTimeSliderVal, initWinSizeVal):
        """
        This will set up the main window layout.
        INPUTS:
            initAmpSliderVal = the initial value for the amplitude slider
            initTimeSliderVal = the intial value for the time slider
            initWinSizeVal = the initial value for the window size (time)
            
        OUTPUTS:
            fig1 = the main figure window (matplotLib Figure object)
            mainDataAxis = Main data axes (matplotLib axis object)
            ampSliderAxis = Amplitude slider for main data (matplotLib axis object)
            timeSliderAxis = Time slider for main data (matplotLib axis object)
            windowSizeSliderAxis = Time slider for main data (matplotLib axis object)
            ampSlider = The amplitude "Slider" Object
            timeSlider = the time "Slider" object
            windowSizeSlider = the window size "Slider" object
        """
        fig1 = plt.figure(2, dpi=100, figsize=[4, 7])  # Sizes hard-coded...
        gs = gridspec.GridSpec(5, 1, height_ratios=[5, 0.5, 0.25, 0.25, 0.25])
        mainDataAxis = fig1.add_subplot(
            gs[0]
        )  # Main data axes (matplotLib axis object)
        ampSliderAxis = fig1.add_subplot(
            gs[2]
        )  # Amplitude slider for main data (matplotLib axis object)
        timeSliderAxis = fig1.add_subplot(
            gs[3]
        )  # Time slider for main data (matplotLib axis object)
        windowSizeSliderAxis = fig1.add_subplot(
            gs[4]
        )  # Time slider for main data (matplotLib axis object)

        ampSlider = Slider(
            ampSliderAxis, "Amplitude", 0, 1, valinit=initAmpSliderVal, valstep=0.001
        )
        timeSlider = Slider(
            timeSliderAxis,
            "Initial Time",
            -0.1,
            0.3,
            valinit=initTimeSliderVal,
            valstep=0.001,
        )
        windowSizeSlider = Slider(
            windowSizeSliderAxis,
            "Window Size",
            0,
            0.5,
            valinit=initWinSizeVal,
            valstep=0.001,
        )

        mainDataAxis.set_ylim([initTimeSliderVal, initTimeSliderVal + initWinSizeVal])
        mainDataAxis.set_xlim([-initAmpSliderVal, initAmpSliderVal])
        mainDataAxis.invert_yaxis()
        mainDataAxis.set_xlabel("Normalized Amplitude")
        mainDataAxis.set_ylabel("Time (s)")
        return (
            fig1,
            mainDataAxis,
            ampSliderAxis,
            timeSliderAxis,
            windowSizeSliderAxis,
            ampSlider,
            timeSlider,
            windowSizeSlider,
        )


class doPicks:
    def __init__(self, x, t, data, shotLoc, initTraceNumb, pickFileName):
        # Define object-level attributes
        self.x = x
        self.t = t
        self.data = data

        # Create both windows with initalized values
        tracePickWindow = tracePickingWindow(
            self.x, self.t, self.data, shotLoc, initTraceNumb
        )
        mainWindowObject = mainPickingWindow(self.x, self.t, self.data, shotLoc)

        self.tracePickWindow = tracePickWindow
        self.mainWindowObject = mainWindowObject

        # Create an attribute to keep track of current trace (this will need travel throughout the class)
        self.cTrace = initTraceNumb
        # Initialize picking objects

        # Local variables to help me calcualte trace index
        gx = np.diff(self.x)[0]
        x0 = self.x[0]

        # Plot picks if they exists
        self.updatePicksMainWindow(mainWindowObject, shotLoc, pickFileName)
        self.updatePicksTraceWindow(tracePickWindow, shotLoc, pickFileName)

        # I couldn't seem to move these out as a method because I call tracePickWindow
        # inside the funciton so it has to be defined locally....
        def whenClickedTraceWindow(event):
            tracePickWindow.mainDataAxis.time_onclick = time.time()

        def whenReleasedTraceWindow(event):
            MAX_CLICK_LENGTH = 0.2
            if event.inaxes == tracePickWindow.mainDataAxis:
                if (
                    event.button == 1
                    and (time.time() - tracePickWindow.mainDataAxis.time_onclick)
                    < MAX_CLICK_LENGTH
                ):
                    tPick = event.ydata
                    self.writePick(shotLoc, tPick, pickFileName)
                    self.updatePicksMainWindow(mainWindowObject, shotLoc, pickFileName)
                    self.updatePicksTraceWindow(tracePickWindow, shotLoc, pickFileName)
                    # print(tPick)
                if (
                    event.button == 3
                    and (time.time() - tracePickWindow.mainDataAxis.time_onclick)
                    < MAX_CLICK_LENGTH
                ):
                    tPick = event.ydata
                    self.deletePick(shotLoc, tPick, pickFileName)
                    self.updatePicksMainWindow(mainWindowObject, shotLoc, pickFileName)
                    self.updatePicksTraceWindow(tracePickWindow, shotLoc, pickFileName)

        def switchTraces(event):
            if event.key == "right":
                cSliderVal = tracePickWindow.ampSlider.val
                cTimeVal = tracePickWindow.timeSlider.val
                cWindowSize = tracePickWindow.windowSizeSlider.val
                self.cTrace = self.cTrace + 1
                tracePickWindow.traceNum = self.cTrace
                self.updatePicksTraceWindow(tracePickWindow, shotLoc, pickFileName)
            elif event.key == "left":
                cSliderVal = tracePickWindow.ampSlider.val
                cTimeVal = tracePickWindow.timeSlider.val
                cWindowSize = tracePickWindow.windowSizeSlider.val
                self.cTrace = self.cTrace - 1
                tracePickWindow.traceNum = (
                    self.cTrace
                )  # Last minute modification...Update this attribute so sliders work better.
                self.updatePicksTraceWindow(tracePickWindow, shotLoc, pickFileName)

        tracePickWindow.figureObject.canvas.mpl_connect(
            "button_press_event", whenClickedTraceWindow
        )
        tracePickWindow.figureObject.canvas.mpl_connect(
            "button_release_event", whenReleasedTraceWindow
        )
        tracePickWindow.figureObject.canvas.mpl_connect("key_press_event", switchTraces)

        def whenClickedMainWindow(event):
            mainWindowObject.mainDataAxis.time_onclick = time.time()

        def getTraceMainWindow(event):
            MAX_CLICK_LENGTH = 0.2
            if event.inaxes == mainWindowObject.mainDataAxis:
                if (
                    event.button == 1
                    and (time.time() - mainWindowObject.mainDataAxis.time_onclick)
                    < MAX_CLICK_LENGTH
                ):
                    xPick = event.xdata
                    traceNum = np.round((xPick - x0) / gx)
                    if traceNum < 0:
                        traceNum = 0
                    elif traceNum > len(self.x):
                        traceNum = len(self.x)

                    cSliderVal = tracePickWindow.ampSlider.val
                    cTimeVal = tracePickWindow.timeSlider.val
                    cWindowSize = tracePickWindow.windowSizeSlider.val
                    self.cTrace = int(traceNum)
                    tracePickWindow.traceNum = self.cTrace
                    traceData = self.data[:, self.cTrace]
                    tracePickWindow.mainDataAxis.clear()
                    tracePickWindow.mainDataAxis.plot(traceData, self.t, "k")
                    tracePickWindow.mainDataAxis.fill_betweenx(
                        self.t,
                        0,
                        traceData,
                        where=traceData < 0,
                        color="k",
                        interpolate=True,
                    )
                    tracePickWindow.mainDataAxis.set_ylim(
                        [cTimeVal, cTimeVal + cWindowSize]
                    )
                    tracePickWindow.mainDataAxis.set_xlim([-cSliderVal, cSliderVal])
                    tracePickWindow.mainDataAxis.invert_yaxis()
                    tracePickWindow.mainDataAxis.set_xlabel("Distance (m)")
                    tracePickWindow.mainDataAxis.set_ylabel("Time (s)")
                    tracePickWindow.figureObject.canvas.draw()

        mainWindowObject.figureObject.canvas.mpl_connect(
            "button_press_event", whenClickedMainWindow
        )
        mainWindowObject.figureObject.canvas.mpl_connect(
            "button_release_event", getTraceMainWindow
        )

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

    def updatePicksMainWindow(self, mainWindowObject, shotLoc, pickFile):
        if os.path.exists(pickFile):
            tempData = np.loadtxt(pickFile)
            shotLocs = tempData[:, 0]
            xPicks = tempData[:, 1]
            tPicks = tempData[:, 2]
            indShots = np.where(shotLocs == shotLoc)
            xPicks = xPicks[indShots]
            tPicks = tPicks[indShots]
        else:
            xPicks = []
            tPicks = []
        mainWindowObject.tPicks = tPicks
        mainWindowObject.xPicks = xPicks
        cSliderVal = mainWindowObject.ampSlider.val
        cTimeVal = mainWindowObject.timeSlider.val
        mainWindowObject.mainDataAxis.clear()
        mainWindowObject.mainDataAxis.pcolorfast(
            np.append(self.x, self.x[-1] + mainWindowObject.gx),
            np.append(self.t, self.t[-1] + mainWindowObject.dt),
            self.data,
            vmin=-cSliderVal,
            vmax=cSliderVal,
            cmap="gray",
        )
        mainWindowObject.mainDataAxis.scatter(xPicks, tPicks, marker=1, s=50, c="c")
        print("made it here")
        #            if not (self.shotLocs == []):
        #               indShots = np.where(self.shotLocs==shotLoc)
        #               self.mainDataAxis.scatter(self.xPicks[indShots],self.tPicks[indShots],marker=1,s=50,c='c')
        mainWindowObject.mainDataAxis.set_ylim([0, cTimeVal])
        mainWindowObject.mainDataAxis.invert_yaxis()
        mainWindowObject.mainDataAxis.set_xlabel("Distance (m)")
        mainWindowObject.mainDataAxis.set_ylabel("Time (s)")
        mainWindowObject.figureObject.canvas.draw()

    def updatePicksTraceWindow(self, traceWindowObject, shotLoc, pickFile):
        if os.path.exists(pickFile):
            tempData = np.loadtxt(pickFile)
            shotLocs = tempData[:, 0]
            xPicks = tempData[:, 1]
            tPicks = tempData[:, 2]
            indShots = np.where(shotLocs == shotLoc)
            xPicks = xPicks[indShots]
            tPicks = tPicks[indShots]
            indTrace = np.where(xPicks == self.x[self.cTrace])
            xPicks = xPicks[indTrace] * 0
            tPicks = tPicks[indTrace]
        else:
            xPicks = []
            tPicks = []
        traceWindowObject.tPicks = tPicks
        traceWindowObject.xPicks = xPicks
        cSliderVal = traceWindowObject.ampSlider.val
        cTimeVal = traceWindowObject.timeSlider.val
        cWindowSize = traceWindowObject.windowSizeSlider.val
        traceWindowObject.traceNum = self.cTrace
        traceData = self.data[:, self.cTrace]
        traceWindowObject.mainDataAxis.clear()
        traceWindowObject.mainDataAxis.plot(traceData, self.t, "k")
        traceWindowObject.mainDataAxis.fill_betweenx(
            self.t, 0, traceData, where=traceData < 0, color="k", interpolate=True
        )
        traceWindowObject.mainDataAxis.scatter(xPicks, tPicks, marker="_", s=200, c="r")
        traceWindowObject.mainDataAxis.set_ylim([cTimeVal, cTimeVal + cWindowSize])
        traceWindowObject.mainDataAxis.set_xlim([-cSliderVal, cSliderVal])
        traceWindowObject.mainDataAxis.invert_yaxis()
        traceWindowObject.mainDataAxis.set_xlabel("Distance (m)")
        traceWindowObject.mainDataAxis.set_ylabel("Time (s)")
        traceWindowObject.figureObject.canvas.draw()


class picker:
    def __init__(self):
        applyBPFilt = True
        # applyBPFilt = False
        lf = 10
        hf = 120
        nq = 500  # Nyquist Frequency
        order = 4
        pickFile = "test.txt"
        picks = np.loadtxt(pickFile)
        # Add function to check for duplicates

        suFile = "dataFiles/shot_01068m_stacked.su"

        # * HARD CODED FOR NOW. We will need to click and open a browser to select
        # This Path
        dirName = "./dataFiles/survey1"
        shotLoc = 918
        fileInfo = getFileInfo(dirName)

        # Hard coded logic to search for shot value (shoult this be exact or appox?)
        tmpShotLocs = np.zeros((len(fileInfo), 1))
        for k in range(0, len(fileInfo)):
            tmpShotLocs[k] = fileInfo[k][1]
        ind = np.argmin(((tmpShotLocs - shotLoc) ** 2) ** 0.5)
        # ******************************************************

        [x, t, data, gx, shotLoc] = getData(
            "segy", os.path.join(dirName, fileInfo[ind][0])
        )
        if applyBPFilt:
            data = bpData(data, lf, hf, nq, order)
        data = normalizeTraces(data)

        self.c = doPicks(x, t, data, shotLoc, 10, pickFile)


if __name__ == "__main__":
    picker()
    plt.show()
