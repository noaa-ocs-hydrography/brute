# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 15:48:05 2019

@author: Casiano.Koprowski
"""

import os
import wx
import ncei_ui
import nceiBAGs

progLoc = os.getcwd()
print(progLoc)

class Form(ncei_ui.Form):
    '''
    '''
    def __init__(self, parent):
        ncei_ui.Form.__init__(self, parent)
        
    def programQuit(self, event):
        '''Closes GUI, ends program.
        Maps to Cancel button, File->Quit, and CTRL+Q
        '''
        self.Close()

    def programProg(self, event):
        '''Collects the GUI field values for use in running the main
        function
        '''
        self.status_bar.SetStatusText('')
        name = self.text_file.GetValue()
        nx = self.text_west.GetValue()
        sy = self.text_south.GetValue()
        sx = self.text_east.GetValue()
        ny = self.text_north.GetValue()
        if nceiBAGs.main(name,nx,sy,sx,ny) == True:
            self.status_bar.SetStatusText('Done!')
        
app = wx.App()
frame = Form(None)
frame.Show()
app.MainLoop()
