# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 13:48:02 2019

@author: Casiano.Koprowski
"""

import wx
import scrape_prog
import eHydro_scrape

class Form(scrape_prog.Form):
    def __init__(self, parent):
        scrape_prog.Form.__init__(self, parent)
    
    def main(self):
        eHydro_scrape.main(self.progress_bar, self.output_text)
#        self.Close()
        
    def programProg(self):
        '''Collects the GUI field values for use in running the main
        function'
        '''
        import threading
        th = threading.Thread(target=self.main)
        th.start()

if __name__ == '__main__':
    app = wx.App()
    frame = Form(None)
    frame.Show()
    frame.programProg()
    app.MainLoop()