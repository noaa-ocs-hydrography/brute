# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 11:55:07 2019

@author: Casiano.Koprowski
"""

import os
from datetime import datetime as _dt

import autointerp_ui
import wx

import autointerp


class Form(autointerp_ui.Form):
    """
    Load ui and: define tif storage columns, overwrite ui defined fucntions
    with desired function behaviour
    """

    def __init__(self, parent):
        autointerp_ui.Form.__init__(self, parent)
        self.insInd = 0
        #Instantiation of 'GeoTIFF File List' box Columns:
        self.list_tif.InsertColumn(0, 'File', width=200)
        self.list_tif.InsertColumn(1, 'Path', width=500)

    def programQuit(self, event):
        """
        Closes GUI, ends program.
        Maps to Cancel button, File->Quit, and CTRL+Q
        """

        self.Close()

    def itemInsert(self, event):
        """
        Adds files selected from the 'Add GeoTIFF File' to the 'GeoTIFF File
        List' box. 'File' holds the name of the file and 'Path' holds the
        complete file path
        """

        print (self.picker_tif.GetPath())
        tif = self.picker_tif.GetPath()
        self.gettifList()
        if tif not in self.tifList:
            name = os.path.split(self.picker_tif.GetPath())
            self.list_tif.InsertItem(self.insInd, name[1])
            self.list_tif.SetItem(self.insInd, 1, tif)
            self.insInd += 1

    def itemRemove(self, event):
        """
        Removes selected files from the 'GeoTIFF File List' box when the
        'Remove' button is clicked
        """

        selected = self.list_tif.SelectedItemCount
        for x in range(0, selected):
            sel = self.list_tif.GetFirstSelected()
            self.list_tif.DeleteItem(sel)
            if self.insInd > 0:
                self.insInd -= 1
        self.gettifList()

    def main(self):
        """
        Main function run as a thread by programProg(). This function
        collects the input field values and passes them to autointerp.py's
        'main' function interp()
        """

        st = 'Started - ' + str(_dt.now())
        self.bar_status.SetStatusText(st)
        self.progressBar.Pulse()
        bagPath = self.picker_bag.GetPath()
        self.gettifList()
        tifPath = self.tifList
        desPath = self.picker_des.GetPath()
        if desPath == '':
            desPath = os.path.split(bagPath)[0]
        catzoc = self.choice_catzoc.GetString(self.choice_catzoc.GetCurrentSelection())
        ioOut = self.radio_data.GetSelection()
        returned = autointerp.main(bagPath, tifPath, desPath, catzoc, ioOut)
        self.bar_status.SetStatusText(returned)
        self.progressBar.SetValue(100)

    def programProg(self, event):
        """
        Start thread to run main()
        """

        import threading
        th = threading.Thread(target=self.main)
        th.start()

    def gettifList(self):
        """
        Used for checking and updating the contents of current collection
        of .tiff files
        """

        tifCount = self.list_tif.GetItemCount()
        self.tifList = []
        for x in range(0, tifCount):
            self.tifList.append(self.list_tif.GetItemText(x, col=1))
        print (self.tifList)

class Done(autointerp_ui.Done):
    def __init__(self, parent):
        autointerp_ui.Done.__init__(self, parent)

if __name__ == '__main__':
    app = wx.App()
    frame = Form(None)
    icon = wx.Icon()
    icon.CopyFromBitmap(wx.Bitmap("autointerp.ico", wx.BITMAP_TYPE_ANY))
    frame.SetIcon(icon)
    frame.Show()
    app.MainLoop()
