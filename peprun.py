import wx
import json
import collections
import execution


STAR_TYPES = (
    'PGM',  # Program star
    'CHK',  # Check star
    'CMP',  # Comparison star
    'STD',  # Standard
    'NPS',  # North Polar Sequence
    'EXT'   # Extinction
)


class RunItem:
    """A VO representing a measurement of a star"""
    def __init__(self, star_id, star_type):
        self.star_id = star_id
        self.star_type = star_type


class Run:
    """Represents a PEP run"""
    def __init__(self):
        self.__filters = set({})
        self.__items = []

    def add_item(self, star_id, star_type, index=-1):
        """Adds an item to the run.

        Returns
        -------
        Zero-based index of item in run"""
        if 0 <= index < len(self.__items):
            self.__items[index:index] = [RunItem(star_id, star_type)]
            return index
        else:
            self.__items.append(RunItem(star_id, star_type))
            return len(self.__items) - 1

    def update_item(self, item_id, star_id, star_type):
        if 0 <= item_id < len(self.__items):
            item = self.__items[item_id]
            item.star_id = star_id
            item.star_type = star_type

    def delete_item(self, item_id):
        if 0 <= item_id < len(self.__items):
            self.__items[item_id:item_id+1] = []

    def len(self):
        return len(self.__items)

    def enable_filter(self, index, enabled):
        if enabled:
            self.__filters.add(index)
        elif index in self.__filters:
            self.__filters.remove(index)

    def get_filters(self):
        return tuple(self.__filters)

    def get_items(self):
        return tuple(self.__items)

    def get_item(self, item_index):
        return self.__items[item_index] if 0 <= item_index < len(self.__items) else None


class RunEncoder(json.JSONEncoder):
    def default(self, o):
        result = collections.OrderedDict()
        result['filters'] = list(o.get_filters())
        items = [{'star_id': item.star_id, 'star_type': item.star_type} for item in o.get_items()]
        result['items'] = items
        return result


def decode_run(dct):
    result = Run()
    if 'filters' in dct:
        for i in dct['filters']:
            result.enable_filter(int(i), True)
        if 'items' in dct:
            for item in dct['items']:
                result.add_item(item['star_id'], item['star_type'])
        return result

    return dct


class RunView(wx.Panel):
    def __init__(self, parent, settings, enable_save_cb=None, status_text_cb=None):
        wx.Panel.__init__(self, parent)
        self.__enable_save_cb = enable_save_cb
        self.__run = Run()
        self.__editing_item = -1

        self.__execution = None
        self.__filter_slots = settings.filter_slots
        self.__device = settings.device
        self.__status_text_cb = status_text_cb

        outer_sizer = wx.BoxSizer(wx.VERTICAL)
        border_width = 10

        outer_sizer.Add(wx.StaticText(self, label='Filters'), flag=wx.ALIGN_LEFT | wx.BOTTOM, border=5)
        self.__filters = []
        filter_sizer = wx.BoxSizer(wx.HORIZONTAL)
        for i in range(6):
            label = settings.filter_slots[i] if i < len(settings.filter_slots) else ''
            filter_subsizer = self.__add_filter(label, i + 1)
            if i > 0:
                filter_sizer.Add(filter_subsizer, flag=wx.ALIGN_CENTER | wx.LEFT, border=border_width)
            else:
                filter_sizer.Add(filter_subsizer, flag=wx.ALIGN_CENTER)
            if i >= len(settings.filter_slots):
                filter_sizer.Hide(filter_subsizer)
            self.__filters.append(filter_subsizer)
        outer_sizer.Add(filter_sizer, flag=wx.ALIGN_LEFT | wx.BOTTOM, border=border_width)

        star_sizer = wx.BoxSizer(wx.HORIZONTAL)
        star_id = wx.TextCtrl(self)
        star_id.SetInitialSize(star_id.GetSizeFromTextSize(star_id.GetFullTextExtent('M' * 20)[0]))
        star_id.Bind(wx.EVT_KEY_DOWN, self.__escape_edit)
        self.__star_id = star_id
        star_sizer.Add(self.__pack_labelled_field('Star ID', star_id),
                       flag=wx.ALIGN_BOTTOM | wx.RIGHT, border=border_width)
        star_type = wx.ComboBox(self, style=wx.CB_READONLY, choices=STAR_TYPES)
        star_type.SetInitialSize(star_type.GetSizeFromTextSize(star_type.GetFullTextExtent('M' * 3)[0]))
        star_type.SetValue(STAR_TYPES[0])
        star_sizer.Add(self.__pack_labelled_field('Type', star_type), flag=wx.ALIGN_BOTTOM)
        self.__star_type = star_type
        height = star_sizer.GetChildren()[0].CalcMin()[1]
        star_sizer.Add(border_width, height, proportion=1)
        self.__add_button = wx.Button(self, label='Add')
        star_sizer.Add(self.__add_button, flag=wx.ALIGN_BOTTOM)
        outer_sizer.Add(star_sizer, flag=wx.ALIGN_LEFT | wx.BOTTOM, border=border_width)
        self.__add_button.Bind(wx.EVT_BUTTON, lambda evt: self.__add_update_item())

        line_height = self.GetFullTextExtent('lj')[1]
        self.__items = wx.ListCtrl(self, size=(-1, 25 * line_height), style=wx.LC_REPORT | wx.BORDER_SUNKEN)
        col_width = self.__items.GetFullTextExtent('M' * 20)[0]
        self.__items.InsertColumn(0, 'Star ID', width=col_width)
        col_width = self.__items.GetFullTextExtent('M' * 4)[0]
        self.__items.InsertColumn(1, 'Type', width=col_width)
        outer_sizer.Add(self.__items, flag=wx.ALIGN_CENTER)
        self.__items.Bind(wx.EVT_KEY_DOWN, self.__delete_item)
        self.__items.Bind(wx.EVT_LIST_ITEM_ACTIVATED, lambda evt: self.__edit_item())

        self.__exec_button = wx.Button(self, label='Execute')
        self.__exec_button.Enable(False)
        self.__exec_button.Bind(wx.EVT_BUTTON, lambda evt: self.__start_execution())
        outer_sizer.Add(self.__exec_button, flag=wx.ALIGN_LEFT | wx.TOP, border=border_width)

        self.SetSizer(outer_sizer)

    def save_run(self, path):
        with open(path, 'w') as out:
            json.dump(self.__run, out, cls=RunEncoder)

    def load_run(self, path):
        with open(path, 'r') as in_file:
            self.__run = json.load(in_file, object_hook=decode_run)
        self.__cancel_edit()
        self.__items.DeleteAllItems()
        filters = self.__run.get_filters()
        for filter_idx, sizer in enumerate(self.__filters, 1):
            sizer.GetChildren()[1].GetWindow().SetValue(filter_idx in filters)
        for item in self.__run.get_items():
            self.__add_item(item.star_id, item.star_type)
        self.__exec_button.Enable(True)

    def __add_filter(self, label_text, index):
        filter_sizer = wx.BoxSizer(wx.HORIZONTAL)
        label = wx.StaticText(self, label=label_text)
        filter_sizer.Add(label, flag=wx.RIGHT | wx.ALIGN_CENTER, border=5)
        checkbox = wx.CheckBox(self)
        enabled = label_text in ('B', 'V')
        checkbox.SetValue(enabled)
        checkbox.Bind(wx.EVT_CHECKBOX, lambda evt: self.__filter_changed(index, evt.IsChecked()))
        self.__run.enable_filter(index, enabled)
        filter_sizer.Add(checkbox, flag=wx.ALIGN_CENTER)

        return filter_sizer

    def __pack_labelled_field(self, label_text, field):
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(wx.StaticText(self, label=label_text), flag=wx.ALIGN_LEFT | wx.BOTTOM, border=5)
        sizer.Add(field, flag=wx.ALIGN_LEFT)
        return sizer

    def __filter_changed(self, index, enabled):
        self.__run.enable_filter(index, enabled)
        self.__enable_save_cb(True)

    def __add_update_item(self):
        star_id = self.__star_id.GetValue().strip()
        if star_id != '':
            star_type = self.__star_type.GetValue()
            if self.__editing_item == -1:
                selected = self.__items.GetFirstSelected()
                self.__add_item(star_id, star_type)
                self.__run.add_item(star_id, star_type, selected)
            else:
                self.__items.SetItem(self.__editing_item, 0, star_id)
                self.__items.SetItem(self.__editing_item, 1, star_type)
                self.__run.update_item(self.__editing_item, star_id, star_type)
                self.__items.Select(self.__editing_item, 0)
                self.__cancel_edit()
            self.__enable_save_cb(True)
            self.__exec_button.Enable(True)

    def __escape_edit(self, event):
        key_code = event.GetKeyCode()
        if key_code == wx.WXK_ESCAPE:
            if self.__editing_item != -1:
                self.__cancel_edit()
        else:
            event.Skip()

    def __cancel_edit(self):
        self.__items.Select(self.__editing_item, 0)
        self.__add_button.SetLabel('Add')
        self.__editing_item = -1

    def __add_item(self, star_id, star_type):
        index = self.__items.GetFirstSelected()
        if index == -1:
            self.__items.Append((star_id, star_type))
        else:
            self.__items.InsertItem(index, star_id)
            self.__items.SetItem(index, 1, star_type)
            self.__items.Select(index + 1, 0)

    def __delete_item(self, event):
        key_code = event.GetKeyCode()
        if key_code in (wx.WXK_DELETE, wx.WXK_BACK) and not self.__execution:
            selected = self.__items.GetFirstSelected()
            if selected != -1:
                self.__items.DeleteItem(selected)
                self.__run.delete_item(selected)
                self.__editing_item = -1
                self.__enable_save_cb(self.__run.len() != 0)
                self.__exec_button.Enable(self.__run.len() != 0)

        elif key_code == wx.WXK_ESCAPE and self.__editing_item != -1:
            self.__cancel_edit()
        else:
            event.Skip()

    def __edit_item(self):
        selected = self.__items.GetFirstSelected()
        if selected != -1:
            self.__editing_item = selected
            self.__star_id.SetValue(self.__items.GetItem(selected, 0).GetText())
            self.__star_type.SetValue(self.__items.GetItem(selected, 1).GetText())
            self.__add_button.SetLabel('Update')

    def __start_execution(self):
        with wx.FileDialog(self, "Select output file", style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as dlg:
            if dlg.ShowModal() == wx.ID_CANCEL:
                return
            path = dlg.GetPath()
        self.__execution = execution.RunExecution(self.__run, self.__filter_slots, self.__device,
                                                  self.__get_operator_confirmation)
        self.__execution.start(path)
        self.__exec_button.Enable(False)

    def __get_operator_confirmation(self, message, need_confirmation, item_index):
        if not need_confirmation:
            wx.CallAfter(self.__finish_execution, item_index - 1)
        else:
            wx.CallAfter(self.__show_confirmation_dlg, message, item_index)

    def __show_confirmation_dlg(self, message, item_index):
        if item_index >= 0:
            self.__items.Select(item_index - 1, 0)
        self.__items.EnsureVisible(item_index)
        self.__items.Select(item_index)
        self.__status_text_cb('Waiting for operator...')
        with wx.MessageDialog(self, message, caption='Move Scope') as dlg:
            dlg.ShowModal()
        self.__status_text_cb('Measuring...')
        self.__execution.operator_ack()

    def __finish_execution(self, item_index):
        self.__items.Select(item_index, 0)
        self.__execution.finish()
        self.__execution = None
        self.__exec_button.Enable(True)
        self.__status_text_cb('Run complete.')
