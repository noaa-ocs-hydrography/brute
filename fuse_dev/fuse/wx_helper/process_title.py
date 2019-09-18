# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 11:37:21 2019

@author: Casiano.Koprowski
"""

import wx
from . import wx_frame


class Form(wx_frame.Form):
    """ """

    def __init__(self, parent):
        wx_frame.Form.__init__(self, parent)

    def close(self):
        self.Close()


class Open_Frame:

    def __init__(self, title: str):
        self.__app__ = wx.App()
        self.__frame__ = Form(None)
        self.__frame__.SetTitle(title)
        self.__frame__.Show()

    def close(self):
        self.__frame__.close()
        self.__app__.MainLoop()
