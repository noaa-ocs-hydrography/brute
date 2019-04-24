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

class Form ( wx.Frame ):

	def __init__( self, parent ):
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = u"NCEI Bag Finder", pos = wx.DefaultPosition, size = wx.Size( 300,350 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )

		self.SetSizeHints( wx.Size( 300,350 ), wx.DefaultSize )
		self.SetBackgroundColour( wx.SystemSettings.GetColour( wx.SYS_COLOUR_WINDOW ) )

		self.menu_bar = wx.MenuBar( 0 )
		self.menu_file = wx.Menu()
		self.menu_quit = wx.MenuItem( self.menu_file, wx.ID_ANY, u"Quit"+ u"\t" + u"CTRL+Q", wx.EmptyString, wx.ITEM_NORMAL )
		self.menu_file.Append( self.menu_quit )

		self.menu_bar.Append( self.menu_file, u"File" )

		self.SetMenuBar( self.menu_bar )

		box_container = wx.BoxSizer( wx.VERTICAL )

		box_filename = wx.FlexGridSizer( 0, 3, 0, 0 )
		box_filename.AddGrowableCol( 1 )
		box_filename.SetFlexibleDirection( wx.BOTH )
		box_filename.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )

		self.label_name = wx.StaticText( self, wx.ID_ANY, u"File Name:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.label_name.Wrap( -1 )

		box_filename.Add( self.label_name, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

		self.text_file = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		box_filename.Add( self.text_file, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

		self.label_txt = wx.StaticText( self, wx.ID_ANY, u".txt", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.label_txt.Wrap( -1 )

		box_filename.Add( self.label_txt, 0, wx.ALIGN_CENTER_VERTICAL|wx.TOP|wx.BOTTOM|wx.RIGHT, 5 )


		box_container.Add( box_filename, 1, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5 )

		box_latlong = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Latitude and Longitude (Decimal Degrees)" ), wx.VERTICAL )

		fgSizer2 = wx.FlexGridSizer( 0, 2, 0, 0 )
		fgSizer2.AddGrowableCol( 1 )
		fgSizer2.SetFlexibleDirection( wx.BOTH )
		fgSizer2.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )

		self.label_north = wx.StaticText( box_latlong.GetStaticBox(), wx.ID_ANY, u"North:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.label_north.Wrap( -1 )

		fgSizer2.Add( self.label_north, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT, 5 )

		self.text_north = wx.TextCtrl( box_latlong.GetStaticBox(), wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.text_north, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

		self.lebel_west = wx.StaticText( box_latlong.GetStaticBox(), wx.ID_ANY, u"West:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.lebel_west.Wrap( -1 )

		fgSizer2.Add( self.lebel_west, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT, 5 )

		self.text_west = wx.TextCtrl( box_latlong.GetStaticBox(), wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.text_west, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

		self.label_south = wx.StaticText( box_latlong.GetStaticBox(), wx.ID_ANY, u"South:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.label_south.Wrap( -1 )

		fgSizer2.Add( self.label_south, 0, wx.ALL|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL, 5 )

		self.text_south = wx.TextCtrl( box_latlong.GetStaticBox(), wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.text_south, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

		self.label_east = wx.StaticText( box_latlong.GetStaticBox(), wx.ID_ANY, u"East:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.label_east.Wrap( -1 )

		fgSizer2.Add( self.label_east, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT, 5 )

		self.text_east = wx.TextCtrl( box_latlong.GetStaticBox(), wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		fgSizer2.Add( self.text_east, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )


		box_latlong.Add( fgSizer2, 1, wx.ALIGN_CENTER_HORIZONTAL, 5 )


		box_container.Add( box_latlong, 1, wx.ALIGN_CENTER_HORIZONTAL, 5 )

		bSizer2 = wx.BoxSizer( wx.HORIZONTAL )

		self.button_prog = wx.Button( self, wx.ID_ANY, u"Run", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer2.Add( self.button_prog, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 10 )

		self.button_quit = wx.Button( self, wx.ID_ANY, u"Quit", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer2.Add( self.button_quit, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 10 )


		box_container.Add( bSizer2, 1, wx.ALIGN_CENTER_HORIZONTAL, 5 )

		self.progress_bar = wx.Gauge( self, wx.ID_ANY, 100, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL )
		self.progress_bar.SetValue( 0 )
		box_container.Add( self.progress_bar, 0, wx.ALL|wx.EXPAND, 5 )


		self.SetSizer( box_container )
		self.Layout()
		self.status_bar = self.CreateStatusBar( 1, wx.STB_SIZEGRIP, wx.ID_ANY )

		self.Centre( wx.BOTH )

		# Connect Events
		self.Bind( wx.EVT_MENU, self.programQuit, id = self.menu_quit.GetId() )
		self.button_prog.Bind( wx.EVT_BUTTON, self.programProg )
		self.button_quit.Bind( wx.EVT_BUTTON, self.programQuit )

	def __del__( self ):
		pass


	# Virtual event handlers, overide them in your derived class
	def programQuit( self, event ):
		event.Skip()

	def programProg( self, event ):
		event.Skip()



