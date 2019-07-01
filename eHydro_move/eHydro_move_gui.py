# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 11:01:24 2019

@author: Casiano.Koprowski
"""

import eHydro_move_ui
import wx

import eHydro_move


class Form(eHydro_move_ui.Form):
    """ """

    def __init__(self, parent):
        eHydro_move_ui.Form.__init__(self, parent)

    def main(self):
        """ """
        eHydro_move._main(self.text_region, self.progressBar, self.text_output)
        self.progressBar.SetRange(100)
        self.progressBar.SetValue(100)
        self.text_region.SetValue('Done!')

    #        self.Close()

    def programProg(self):
        """
        Collects the GUI field values for use in running the main
        function'

        Parameters
        ----------

        Returns
        -------

        """
        import threading
        th = threading.Thread(target=self.main)
        th.start()


if __name__ == '__main__':
    app = wx.App()
    frame = Form(None)
    frame.Show()
    frame.programProg()
    app.MainLoop()
