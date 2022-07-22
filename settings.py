import wx
from pathlib import Path
from configparser import ConfigParser

CONFIG_DIR = Path.home().joinpath('.peppy')
CONFIG_PATH = CONFIG_DIR.joinpath('config.ini')


def get_settings(parent):
    if CONFIG_PATH.exists():
        return Settings.load()
    else:
        return edit_settings(parent, Settings())


def edit_settings(parent, current_settings):
    with SettingsDlg(parent, current_settings) as dlg:
        if dlg.ShowModal() == wx.ID_OK:
            dlg.update_settings(current_settings)
            current_settings.save()
            return current_settings
        else:
            return None


class Settings:
    """Represents the application settings"""

    def __init__(self, device='', filter_control=False, filter_slots=(), data_dir=CONFIG_DIR):
        self.device = device
        self.filter_control = filter_control
        self.filter_slots = list(filter_slots)
        self.data_dir = data_dir

    @staticmethod
    def load():
        """Load settings from configuration file or create new ones"""
        if not CONFIG_PATH.exists():
            return None

        config = ConfigParser()
        config.read(str(CONFIG_PATH))
        dft_sect = config['DEFAULT']
        device = dft_sect['Device']
        filter_control = bool(dft_sect['FilterControl'])
        filter_slots = []
        if filter_control:
            for i in range(6):
                slot_key = 'Slot{:d}'.format(i + 1)
                if slot_key in dft_sect:
                    filter_slots.append(dft_sect[slot_key])
                else:
                    break
        data_dir = Path(dft_sect['DataDir'])

        return Settings(device, filter_control, filter_slots, data_dir)

    def save(self):
        CONFIG_DIR.mkdir(parents=True, exist_ok=True)
        dft_sect = {'Device': self.device,
                    'FilterControl': self.filter_control,
                    'DataDir': self.data_dir}
        if self.filter_control:
            for slot_idx, slot in enumerate(self.filter_slots, 1):
                slot_key = 'Slot{:d}'.format(slot_idx)
                dft_sect[slot_key] = slot
        config = ConfigParser()
        config['DEFAULT'] = dft_sect
        with open(str(CONFIG_PATH), 'w') as config_file:
            config.write(config_file)


class SettingsDlg(wx.Dialog):
    """Dialog for entering application settings"""

    def __init__(self, parent, settings):
        wx.Dialog.__init__(self, parent)
        self.SetTitle('Settings')

        outer_sizer = wx.BoxSizer(wx.VERTICAL)
        border_width = 5

        top_sizer = wx.BoxSizer(wx.VERTICAL)
        top_sizer.Add(wx.StaticText(self, label='Device'), flag=wx.ALIGN_LEFT | wx.BOTTOM, border=border_width)
        existing_devices = [str(p) for p in Path('/dev').iterdir() if
                            p.stem.startswith('ttyS') or p.stem.startswith('ttyUSB')]
        self.__device_ctrl = wx.ComboBox(self, style=wx.CB_DROPDOWN, choices=existing_devices,
                                         value=settings.device)
        self.__device_ctrl.SetInitialSize(
            self.__device_ctrl.GetSizeFromTextSize(self.__device_ctrl.GetFullTextExtent('/dev/ttyUSB999')[0]))
        top_sizer.Add(self.__device_ctrl, flag=wx.ALIGN_LEFT | wx.BOTTOM, border=border_width)

        filter_ctrl_sizer = wx.BoxSizer(wx.HORIZONTAL)
        filter_ctrl_sizer.Add(wx.StaticText(self, label='Filter Control'), flag=wx.ALIGN_CENTER | wx.RIGHT,
                              border=border_width)
        self.__filter_control_ctrl = wx.CheckBox(self)
        self.__filter_control_ctrl.SetValue(settings.filter_control)
        filter_ctrl_sizer.Add(self.__filter_control_ctrl, flag=wx.ALIGN_CENTER)
        top_sizer.Add(filter_ctrl_sizer, flag=wx.ALIGN_LEFT | wx.TOP | wx.BOTTOM, border=border_width)

        top_sizer.Add(wx.StaticText(self, label='Filter Slots'),
                      flag=wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP, border=border_width)

        filter_slots_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.__filter_slots = []
        for slot in range(6):
            filter_slot = wx.TextCtrl(self)
            filter_slot.SetMaxLength(1)
            filter_slot.SetInitialSize(filter_slot.GetSizeFromTextSize(filter_slot.GetFullTextExtent('CC')[0]))
            filter_slots_sizer.Add(filter_slot, flag=wx.ALIGN_CENTER | wx.RIGHT, border=border_width)
            self.__filter_slots.append(filter_slot)
            if slot < len(settings.filter_slots):
                filter_slot.SetValue(settings.filter_slots[slot])
        top_sizer.Add(filter_slots_sizer, flag=wx.ALIGN_LEFT | wx.BOTTOM, border=border_width)

        top_sizer.Add(wx.StaticText(self, label='Data Directory'), flag=wx.ALIGN_LEFT | wx.TOP | wx.BOTTOM,
                      border=border_width)
        config_path = str(settings.data_dir)
        self.__data_dir_ctrl = wx.DirPickerCtrl(self, path=str(config_path), style=wx.DIRP_USE_TEXTCTRL | wx.DIRP_SMALL)
        size_ctrl = self.__data_dir_ctrl.GetTextCtrl()

        size_ctrl.SetInitialSize(size_ctrl.GetSizeFromTextSize(
            size_ctrl.GetFullTextExtent(str(config_path) + 'XX')[0]
        ))
        top_sizer.Add(self.__data_dir_ctrl, flag=wx.ALIGN_LEFT | wx.BOTTOM, border=border_width)

        button_sizer = self.CreateStdDialogButtonSizer(wx.OK | wx.CANCEL)
        top_sizer.Add(button_sizer, flag=wx.ALIGN_CENTER | wx.TOP, border=2 * border_width)

        outer_sizer.Add(top_sizer, flag=wx.ALL, border=2 * border_width)
        self.SetSizer(outer_sizer)
        minx, miny = outer_sizer.CalcMin()
        miny = miny + button_sizer.CalcMin()[1]
        self.SetInitialSize((minx, miny))
        print(top_sizer.CalcMin())

        self.Show(True)

    def update_settings(self, settings):
        settings.device = self.__device_ctrl.GetValue()
        settings.filter_control = self.__filter_control_ctrl.GetValue()
        if settings.filter_control:
            settings.filter_slots = []
            for slot_ctrl in self.__filter_slots:
                settings.filter_slots.append(slot_ctrl.GetValue())
        settings.data_dir = Path(self.__data_dir_ctrl.GetPath())
