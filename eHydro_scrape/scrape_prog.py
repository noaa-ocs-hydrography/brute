# -*- coding: utf-8 -*-

###########################################################################
## Python code generated with wxFormBuilder (version Dec 12 2018)
## http://www.wxformbuilder.org/
##
## PLEASE DO *NOT* EDIT THIS FILE!
###########################################################################

import wx
import wx.xrc


###########################################################################
## Class Form
###########################################################################

class Form(wx.Frame):
    """ """

    def __init__(self, parent):
        wx.Frame.__init__(self, parent, id=wx.ID_ANY, title=u"eHydro_scrape Progress", pos=wx.DefaultPosition,
                          size=wx.Size(800, 500), style=wx.DEFAULT_FRAME_STYLE | wx.TAB_TRAVERSAL)

        self.SetSizeHints(wx.Size(500, 300), wx.DefaultSize)
        self.SetForegroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
        self.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))

        output_box = wx.BoxSizer(wx.VERTICAL)

        self.progress_bar = wx.Gauge(self, wx.ID_ANY, 100, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL)
        self.progress_bar.SetValue(0)
        output_box.Add(self.progress_bar, 0, wx.ALL | wx.EXPAND, 5)

        self.output_text = wx.TextCtrl(self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize,
                                       wx.HSCROLL | wx.TE_MULTILINE | wx.TE_READONLY)
        output_box.Add(self.output_text, 1, wx.ALL | wx.EXPAND, 5)

        self.SetSizer(output_box)
        self.Layout()
        self.status_bar = self.CreateStatusBar(1, wx.STB_SIZEGRIP, wx.ID_ANY)

        self.Centre(wx.BOTH)

    def __del__(self):
        pass
