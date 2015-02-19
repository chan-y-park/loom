import os
import logging
import Tkinter as tk
import Tkconstants
import tkFileDialog
import pdb

from api import generate_spectral_network, load_spectral_network

class Application(tk.Frame):
    def __init__(self, config, master=None):
        tk.Frame.__init__(self, master)
        self.config = config
        self.pack()
        self.createWidgets()

    def createWidgets(self):
        tk.Label(self, text='curve').grid(row=0, column=0)

        self.entry_curve_text = tk.StringVar()
        self.entry_curve = tk.Entry(self, textvariable=self.entry_curve_text) 
        self.entry_curve_text.set(self.config['sw_curve'])
        self.entry_curve.grid(row=0, column=1)

        self.button_generate = tk.Button(
            self, text='generate',
            command=lambda: generate_spectral_network(self.config),
        )
        self.button_generate.grid(row=1, column=2, sticky=tk.E)

def open_gui(opts, config):
    root = tk.Tk()
    app = Application(config, master=root)
    app.mainloop()

def gui_load_spectral_network(opts):
    root = tk.Tk()
    dir_opts = {
        'initialdir': os.curdir,
        'mustexist': False,
        'parent': root,
        'title': 'Select a directory that contains data files.',
    }
    data_dir = tkFileDialog.askdirectory(**dir_opts)
    root.destroy()
    if data_dir == '':
        return None
    else:
        logging.info('Opening data directory "{}"...'.format(data_dir))
        opts['load-data'] = data_dir
        return load_spectral_network(opts)
    
