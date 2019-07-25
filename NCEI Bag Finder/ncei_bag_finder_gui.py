# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 15:48:05 2019

@author: Casiano.Koprowski
"""

import os
import webbrowser

import nceiBAGs
import ncei_ui
import wx

progLoc = os.getcwd()
print(progLoc)


class Form(ncei_ui.Form):
    """ """

    def __init__(self, parent):
        ncei_ui.Form.__init__(self, parent)

    def main(self):
        """ """
        self.status_bar.SetStatusText('')
        name = self.text_file.GetValue()
        xmin = float(self.text_west.GetValue())
        ymin = float(self.text_south.GetValue())
        xmax = float(self.text_east.GetValue())
        ymax = float(self.text_north.GetValue())
        sel = self.radio_query.GetSelection()
        bcm = self.bound_check(xmin, ymin, xmax, ymax)
        if not sel:
            qId = 3
        if sel:
            qId = 0
        if bcm is not None:
            self.status_bar.SetStatusText(bcm)
        else:
            msg = nceiBAGs.main(name, xmin, ymin, xmax, ymax, qId, self.progress_bar)
            if msg is None:
                self.status_bar.SetStatusText('Done!')
            else:
                self.progress_bar.SetValue(0)
                self.status_bar.SetStatusText(msg)

    def bound_check(self, xmin, ymin, xmax, ymax):
        xaxis = xmin < xmax
        yaxis = ymin < ymax
        if not xaxis and not yaxis:
            return 'Please check your North ({ymax}), South ({ymin}), East({xmax}), and West ({xmin})'
        elif not xaxis:
            return f'Please check your East ({xmax}) and West ({xmin})'
        elif not yaxis:
            return f'Please check your North ({ymax}) and South ({ymin})'

    def programAbout(self, event):
        webbrowser.open(r'https://vlab.ncep.noaa.gov/web/national-bathymetric-source/blogs/-/blogs/getting-survey-info-via-ncei-s-rest-api', new=2, autoraise=True)

    def programQuit(self, event):
        """
        Closes GUI, ends program.
        Maps to Cancel button, File->Quit, and CTRL+Q

        Parameters
        ----------
        event :


        Returns
        -------

        """
        self.Close()

    def programProg(self, event):
        """
        Collects the GUI field values for use in running the main
        function

        Parameters
        ----------
        event :


        Returns
        -------

        """
        import threading
        th = threading.Thread(target=self.main)
        th.start()


app = wx.App()
frame = Form(None)
frame.Show()
app.MainLoop()
