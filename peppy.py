#!/usr/bin/python3
import wx
import peprun
import settings


class MainWindow(wx.Frame):
    """Main application window"""

    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title=title, size=(-1, -1))
        self.__current_run = None
        self.CreateStatusBar()

        file_menu = wx.Menu()
        menu_open = file_menu.Append(wx.ID_OPEN, "&Open", "Open saved run")

        self.__menu_save = file_menu.Append(wx.ID_SAVE, "&Save", "Save run")
        file_menu.Enable(wx.ID_SAVE, False)

        menu_settings = file_menu.Append(wx.ID_PREFERENCES, "Se&ttings", "Edit settings")

        file_menu.AppendSeparator()
        menu_exit = file_menu.Append(wx.ID_EXIT, "E&xit", "Exit program")

        menu_bar = wx.MenuBar()
        menu_bar.Append(file_menu, "&File")
        self.SetMenuBar(menu_bar)

        self.Bind(wx.EVT_MENU, self.on_open, menu_open)
        self.Bind(wx.EVT_MENU, self.on_save, self.__menu_save)
        self.Bind(wx.EVT_MENU, self.on_exit, menu_exit)
        self.Bind(wx.EVT_MENU, self.on_edit_settings, menu_settings)
        self.Show(True)

        self.__settings = settings.get_settings(self)
        if not self.__settings:
            self.Close(True)

        sizer = wx.BoxSizer(wx.VERTICAL)
        self.__run_view = peprun.RunView(self, self.__settings,
                                         lambda enable: file_menu.Enable(wx.ID_SAVE, enable),
                                         lambda msg: self.SetStatusText(msg))
        sizer.Add(self.__run_view, flag=wx.ALL, border=5)
        self.SetSizer(sizer)
        menu_bar.Layout()
        minw, minh = self.GetEffectiveMinSize()
        minh = minh + menu_bar.GetSize()[1]
        self.SetInitialSize((minw, minh))

    def on_open(self, event):
        run_path = None
        with wx.FileDialog(self, "Open Run", style=wx.FD_OPEN) as dlg:
            if dlg.ShowModal() == wx.ID_CANCEL:
                return
            run_path = dlg.GetPath()
        self.__run_view.load_run(run_path)
        self.__current_run = run_path

    def on_save(self, event):
        if not self.__current_run:
            with wx.FileDialog(self, "Save Run", style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as dlg:
                if dlg.ShowModal() == wx.ID_CANCEL:
                    return
                self.__current_run = dlg.GetPath()

        self.__run_view.save_run(self.__current_run)

    def on_exit(self, event):
        self.Close(True)

    def on_edit_settings(self, event):
        settings.edit_settings(self, self.__settings)


if __name__ == '__main__':
    app = wx.App(False)
    MainWindow(None, "PEPpy")
    app.MainLoop()
