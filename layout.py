from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx as NavigationToolbar
from matplotlib.figure import Figure
# import matplotlib.pyplot as plt

import numpy as np

import wx
import wx.lib.mixins.inspection as WIT

from picker import picker

###### Example for displaying a MatPlotLib Figure! ######
class FigurePanel(wx.Panel):
    def __init__(self, parent, figure_window):
        super().__init__(parent)

        self.figure = Figure()
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.subplot = self.figure.subplots(1)
        figure_window.plot_data(self.subplot)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, wx.SizerFlags(1).Left().Top().Shaped())
        self.SetSizer(self.sizer)
        self.Fit()

        self.add_toolbar()  # comment this out for no toolbar

    def add_toolbar(self):
        self.toolbar = NavigationToolbar(self.canvas)
        self.toolbar.Realize()
        # By adding toolbar in sizer, we are able to put it at the bottom
        # of the frame - so appearance is closer to GTK version.
        self.sizer.Add(self.toolbar, wx.SizerFlags(0).Left())
        # update the axes menu on the toolbar
        self.toolbar.update()


class LeftSection(wx.BoxSizer):
    def __init__(self, panel):
        super().__init__(wx.VERTICAL)

        # Top
        topsizer = wx.BoxSizer(wx.VERTICAL)
        topsizer.Add(
            wx.StaticText(panel, label="Enter a number"),
            wx.SizerFlags(0).Centre().TripleBorder(wx.ALL),
        )
        topsizer.Add(wx.Button(panel, label="Btn1"), wx.SizerFlags(0).Expand())
        topsizer.Add(wx.Button(panel, label="Btn2"), wx.SizerFlags(0).Centre())
        topsizer.Add(wx.TextCtrl(panel), wx.SizerFlags(1).Expand())

        # Bottom
        bottomsizer = wx.BoxSizer(wx.HORIZONTAL)
        bottomsizer.Add(
            wx.StaticText(panel, label="Label2"),
            wx.SizerFlags(1).Expand().Border(wx.ALL),
        )
        bottomsizer.Add(wx.Button(panel, label="Btn3"), wx.SizerFlags(1).Right())

        # Add Top and Bottom to the left sizer
        self.Add(topsizer, wx.SizerFlags(1).Expand())
        self.Add(bottomsizer, wx.SizerFlags(1).Expand())


class MiddleSection(wx.BoxSizer):
    def __init__(self, panel, figure_window):
        super().__init__(wx.VERTICAL)

        # Top
        topsizer = wx.BoxSizer(wx.VERTICAL)
        topsizer.Add(
            FigurePanel(panel, figure_window), wx.SizerFlags(1).Shaped().Border(wx.BOTTOM)
        )

        # Bottom
        bottomsizer = wx.BoxSizer(wx.HORIZONTAL)
        bottomsizer.Add(wx.CheckBox(panel, label="Normalize"), wx.SizerFlags(0))
        bottomsizer.Add(wx.CheckBox(panel, label="Wiggle"), wx.SizerFlags(0))
        bottomsizer.Add(wx.CheckBox(panel, label="BP Filter"), wx.SizerFlags(0))
        bottomsizer.Add(
            wx.StaticText(panel, label="Low = 10 Hz, high = 120 Hz, Order = 4"),
            wx.SizerFlags(0),
        )

        # Add Top and Bottom to the left sizer
        self.Add(topsizer, wx.SizerFlags(1).Expand())
        self.Add(bottomsizer, wx.SizerFlags(1).Centre().Expand())


class RightSection(wx.BoxSizer):
    def __init__(self, panel, figure_window):
        super().__init__(wx.VERTICAL)
        self.Add(FigurePanel(panel, figure_window), wx.SizerFlags(1).Shaped())


class MainFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.menubar()

        self.picker = picker()
        self.gui_layout()

        self.Centre()
        self.Show()

    def gui_layout(self):
        panel = wx.Panel(self)

        ### LEFT SECTION ###
        leftsection = LeftSection(panel)

        ### MIDDLE SECTION ###
        middlesection = MiddleSection(
            panel, self.picker.c.mainWindowObject
        )

        ### RIGHT SECTION ###
        rightsection = RightSection(panel, self.picker.c.tracePickWindow)

        ### Add left and right sections to the main sizer
        mainsizer = wx.BoxSizer(wx.HORIZONTAL)

        mainsizer.Add(leftsection, wx.SizerFlags(1).Centre().Expand().Border(wx.ALL))
        mainsizer.Add(middlesection, wx.SizerFlags(1).Centre().Expand().Border(wx.ALL))
        mainsizer.Add(rightsection, wx.SizerFlags(1).Centre().Expand().Border(wx.ALL))

        panel.SetSizerAndFit(mainsizer)
        self.SetSizeHints(self.Size)

    def menubar(self):
        # Define Menu Bar
        menu_bar = wx.MenuBar()

        #### File
        file_menu = wx.Menu()

        # Open
        open_submenu = wx.Menu()
        datafile_item = open_submenu.Append(wx.ID_ANY, "Data File")
        file_menu.AppendSubMenu(open_submenu, "Open")

        # Quit
        exit_button = file_menu.Append(wx.ID_EXIT, "Quit\tCtrl+q", "Quit the app")
        self.Bind(wx.EVT_MENU, self.quit, exit_button)

        #### View
        view_menu = wx.Menu()

        # Reciprocals
        reciprocals_submenu = wx.Menu()
        format_item = reciprocals_submenu.Append(wx.ID_ANY, "Format")
        view_menu.AppendSubMenu(reciprocals_submenu, "Reciprocals")

        # Append File and View to the Menu Bar
        menu_bar.Append(file_menu, "&File")
        menu_bar.Append(view_menu, "View")

        # Set the Menu Bar
        self.SetMenuBar(menu_bar)

    def quit(self, event):
        self.Close()


class App(wx.App):
    def __init__(self):
        "Create the main window and insert the custom frame"
        super().__init__()
        MainFrame(
            None, title="Seismic First Arrival Picker Layout v1", size=(1400, 825)
        )


if __name__ == "__main__":
    app = App()
    app.MainLoop()
