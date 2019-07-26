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
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = u"Interpolation Tool", pos = wx.DefaultPosition, size = wx.Size( 500,625 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )

		self.SetSizeHints( wx.Size( 500,625 ), wx.DefaultSize )
		self.SetBackgroundColour( wx.SystemSettings.GetColour( wx.SYS_COLOUR_WINDOW ) )

		self.bar_menu = wx.MenuBar( 0 )
		self.bar_menu.SetForegroundColour( wx.SystemSettings.GetColour( wx.SYS_COLOUR_MENU ) )
		self.bar_menu.SetBackgroundColour( wx.SystemSettings.GetColour( wx.SYS_COLOUR_WINDOWFRAME ) )

		self.menu_file = wx.Menu()
		self.menu_quit = wx.MenuItem( self.menu_file, wx.ID_ANY, u"Quit"+ u"\t" + u"CTRL+Q", wx.EmptyString, wx.ITEM_NORMAL )
		self.menu_file.Append( self.menu_quit )

		self.bar_menu.Append( self.menu_file, u"File" )

		self.SetMenuBar( self.bar_menu )

		container_form = wx.BoxSizer( wx.VERTICAL )

		self.label_title = wx.StaticText( self, wx.ID_ANY, u"Interpolation Tool", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.label_title.Wrap( -1 )

		container_form.Add( self.label_title, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

		opts_form = wx.FlexGridSizer( 0, 2, 0, 0 )
		opts_form.AddGrowableCol( 1 )
		opts_form.AddGrowableRow( 2 )
		opts_form.AddGrowableRow( 4 )
		opts_form.SetFlexibleDirection( wx.BOTH )
		opts_form.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )

		self.label_bag = wx.StaticText( self, wx.ID_ANY, u"BAG File:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.label_bag.Wrap( -1 )

		opts_form.Add( self.label_bag, 0, wx.ALL|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL, 5 )

		self.picker_bag = wx.FilePickerCtrl( self, wx.ID_ANY, wx.EmptyString, u"Select a file", u"BAG Files (*.bag)|*.bag", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
		opts_form.Add( self.picker_bag, 1, wx.ALL|wx.EXPAND, 5 )

		self.label_tif = wx.StaticText( self, wx.ID_ANY, u"Add Coverage File:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.label_tif.Wrap( -1 )

		opts_form.Add( self.label_tif, 0, wx.ALL|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL, 5 )

		self.picker_tif = wx.FilePickerCtrl( self, wx.ID_ANY, wx.EmptyString, u"Select a file", u"GeoTIFF Files (*.tiff; *tif)|*.tiff; *.tif|Oultine Files (*.shp; *.gpkg)|*.shp; *.gpkg", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
		opts_form.Add( self.picker_tif, 1, wx.ALL|wx.EXPAND, 5 )

		self.label_tifList = wx.StaticText( self, wx.ID_ANY, u"File List:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.label_tifList.Wrap( -1 )

		opts_form.Add( self.label_tifList, 0, wx.ALL|wx.ALIGN_RIGHT, 5 )

		opts_tifs = wx.BoxSizer( wx.VERTICAL )

		self.list_tif = wx.ListCtrl( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LC_REPORT|wx.HSCROLL )
		opts_tifs.Add( self.list_tif, 1, wx.ALL|wx.EXPAND, 5 )

		self.button_tifRemove = wx.Button( self, wx.ID_ANY, u"Remove", wx.DefaultPosition, wx.DefaultSize, 0 )
		opts_tifs.Add( self.button_tifRemove, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )


		opts_form.Add( opts_tifs, 1, wx.EXPAND, 5 )

		self.label_des = wx.StaticText( self, wx.ID_ANY, u"Destination (Optional):", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.label_des.Wrap( -1 )

		opts_form.Add( self.label_des, 0, wx.ALL|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL, 5 )

		self.picker_des = wx.DirPickerCtrl( self, wx.ID_ANY, wx.EmptyString, u"Select a folder", wx.DefaultPosition, wx.DefaultSize, wx.DIRP_DEFAULT_STYLE )
		opts_form.Add( self.picker_des, 1, wx.ALL|wx.EXPAND, 5 )

		self.label_catzoc = wx.StaticText( self, wx.ID_ANY, u"Interpolated Uncertainty CATZOC Score:", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_RIGHT )
		self.label_catzoc.Wrap( 150 )

		opts_form.Add( self.label_catzoc, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT, 5 )

		choice_catzocChoices = [ u"A2/B", u"C", u"A1" ]
		self.choice_catzoc = wx.Choice( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, choice_catzocChoices, 0 )
		self.choice_catzoc.SetSelection( 0 )
		opts_form.Add( self.choice_catzoc, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

		self.label_depth = wx.StaticText( self, wx.ID_ANY, u"Depth Threshold for Interpolated Data:", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_RIGHT )
		self.label_depth.Wrap( 150 )

		opts_form.Add( self.label_depth, 0, wx.ALL|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL, 5 )

		self.picker_depth = wx.SpinCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.SP_ARROW_KEYS, 0, 100, 0 )
		opts_form.Add( self.picker_depth, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )


		container_form.Add( opts_form, 1, wx.EXPAND|wx.FIXED_MINSIZE, 5 )

		radio_dataChoices = [ u"All Data", u"Interpolated Only" ]
		self.radio_data = wx.RadioBox( self, wx.ID_ANY, u"Output Data", wx.DefaultPosition, wx.DefaultSize, radio_dataChoices, 2, wx.RA_SPECIFY_COLS )
		self.radio_data.SetSelection( 0 )
		container_form.Add( self.radio_data, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )

		opts_prog = wx.StdDialogButtonSizer()
		self.opts_progOK = wx.Button( self, wx.ID_OK )
		opts_prog.AddButton( self.opts_progOK )
		self.opts_progCancel = wx.Button( self, wx.ID_CANCEL )
		opts_prog.AddButton( self.opts_progCancel )
        opts_prog.Realize()

		container_form.Add( opts_prog, 0, wx.BOTTOM|wx.EXPAND|wx.TOP, 5 )

		self.progressBar = wx.Gauge( self, wx.ID_ANY, 100, wx.DefaultPosition, wx.DefaultSize, wx.GA_HORIZONTAL )
		self.progressBar.SetValue( 0 )
		container_form.Add( self.progressBar, 0, wx.ALL|wx.EXPAND, 5 )


		self.SetSizer( container_form )
		self.Layout()
		self.bar_status = self.CreateStatusBar( 1, wx.STB_SIZEGRIP, wx.ID_ANY )

		self.Centre( wx.BOTH )

		# Connect Events
		self.Bind( wx.EVT_MENU, self.programQuit, id = self.menu_quit.GetId() )
		self.picker_tif.Bind( wx.EVT_FILEPICKER_CHANGED, self.itemInsert )
		self.button_tifRemove.Bind( wx.EVT_BUTTON, self.itemRemove )
		self.opts_progCancel.Bind( wx.EVT_BUTTON, self.programQuit )
		self.opts_progOK.Bind( wx.EVT_BUTTON, self.programProg )

	def __del__( self ):
		pass


	# Virtual event handlers, overide them in your derived class
	def programQuit( self, event ):
		event.Skip()

	def itemInsert( self, event ):
		event.Skip()

	def itemRemove( self, event ):
		event.Skip()


	def programProg( self, event ):
		event.Skip()


###########################################################################
## Class Done
###########################################################################

class Done ( wx.Dialog ):

	def __init__( self, parent ):
		wx.Dialog.__init__ ( self, parent, id = wx.ID_ANY, title = u"Interpolation Tool", pos = wx.DefaultPosition, size = wx.Size( 350,150 ), style = wx.DEFAULT_DIALOG_STYLE )

		self.SetSizeHints( wx.DefaultSize, wx.DefaultSize )

		container_dialog = wx.BoxSizer( wx.VERTICAL )

		self.label_done = wx.StaticText( self, wx.ID_ANY, u"Done!", wx.Point( -1,-1 ), wx.DefaultSize, wx.ALIGN_CENTER_HORIZONTAL )
		self.label_done.Wrap( -1 )

		container_dialog.Add( self.label_done, 1, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL|wx.EXPAND, 5 )

		opts_done = wx.StdDialogButtonSizer()
		self.opts_doneOK = wx.Button( self, wx.ID_OK )
		opts_done.AddButton( self.opts_doneOK )
		self.opts_doneCancel = wx.Button( self, wx.ID_CANCEL )
		opts_done.AddButton( self.opts_doneCancel )
        opts_done.Realize()

		container_dialog.Add( opts_done, 0, wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL, 5 )


		self.SetSizer( container_dialog )
		self.Layout()

		self.Centre( wx.BOTH )

	def __del__( self ):
		pass


