import os
import logging
import Tkinter as tk
import tkFileDialog
import pdb

from config import LoomConfig
from api import (generate_spectral_network, load_spectral_network,)

class Application(tk.Frame):
    def __init__(self, config, master=None):
        tk.Frame.__init__(self, master)
        
        self.config = config
        self.entry = {}
        self.entry_var = {} 
        self.mb = None
        self.check = {}
        self.pack()
        self.create_widgets()

    def create_widgets(self):
        # Layout variables
        grid_row = 0
        grid_col = 0

        # Menu
        self.mb = tk.Menubutton(self, text='File', relief=tk.RAISED)
        self.mb.grid(row=grid_row, column=grid_col)
        self.mb.menu = tk.Menu(self.mb)
        self.mb['menu'] = self.mb.menu
        self.mb.menu.add_command(
            label='Load',
            command=self.menu_load_action
        )

        # Associate each config option to an Entry
        for option, value in self.config.iteritems():
            self.entry_var[option] = tk.StringVar()
            self.entry_var[option].set(value)
            self.entry[option] = tk.Entry(
                self,
                textvariable=self.entry_var[option]
            )
            
        # Entry & Label layout
        grid_row += 1
        grid_col = 0
        tk.Label(self, text='sw_curve').grid(row=grid_row, column=grid_col)
        grid_col += 1
        self.entry['sw_curve'].grid(row=grid_row, column=grid_col)

        grid_col += 1
        tk.Label(self, text='sw_diff_v').grid(row=grid_row, column=grid_col)
        grid_col += 1
        self.entry['sw_diff_v'].grid(row=grid_row, column=grid_col)

        grid_col += 1
        tk.Label(self, text='root_system').grid(row=grid_row, column=grid_col)
        grid_col += 1
        self.entry['root_system'].grid(row=grid_row, column=grid_col)

        grid_row += 1
        grid_col = 0
        tk.Label(self, text='sw_parameters').grid(row=grid_row, column=grid_col)
        grid_col += 1
        self.entry['sw_parameters'].grid(
            row=grid_row, column=grid_col,
            columnspan=3,
            sticky=tk.EW,
        )

        grid_row += 1
        grid_col = 0
        tk.Label(self, text='phase_range').grid(row=grid_row, column=grid_col)
        grid_col += 1
        self.entry['phase_range'].grid(row=grid_row, column=grid_col)

        # Check plot_on_cylinder
        grid_row += 1
        grid_col = 0
        self.check['plot_on_cylinder'] = tk.IntVar()
        tk.Checkbutton(
            self,
            text='Plot on cyliner',
            variable=self.check['plot_on_cylinder']
        ).grid(row=grid_row, column=grid_col)

        # Set default parameters.
        for option in self.entry_var:
            self.entry_var[option].set(self.config[option])

        # 'Generate' button
        grid_row += 1
        grid_col = 0
        self.button_generate = tk.Button(
            self, text='generate',
            command=self.button_generate_action,
        )
        grid_col += 1
        self.button_generate.grid(row=grid_row, column=grid_col, sticky=tk.E)

    def menu_load_action(self):
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
            check_plot_on_cylinder = self.check['plot_on_cylinder'].get()
            if check_plot_on_cylinder == 1:
                plot_on_cylinder = True
            else: 
                plot_on_cylinder = False
            return load_spectral_network(data_dir, plot_on_cylinder)

    def button_generate_action(self):
        # Read config options from Entries.
        for section in self.config.parser.sections():
            if (section == 'directories'):
                continue
            elif (section == 'Seiberg-Witten parameters'):
                params_input = eval(self.entry_var['sw_parameters'].get())
                self.config['sw_parameters'] = params_input
                for var, val in params_input.iteritems():
                    self.config.parser.set(section, var, str(val))
                continue
            for option in self.config.parser.options(section):
                value = self.entry_var[option].get()
                if (section == 'symbolic expressions'):
                    self.config[option] = value
                elif (section == 'numerical parameters'):
                    self.config[option] = eval(value)
                self.config.parser.set(section, option, value)
        generate_spectral_network(self.config, phase=None, show_plot=True,
                                  plot_on_cylinder=False)

def open_gui(config):
    root = tk.Tk()
    # Put the window at the center
    ws = root.winfo_screenwidth()
    hs = root.winfo_screenheight()
    root.geometry('+{}+{}'.format(ws/2, hs/2))

    app = Application(config, master=root)
    app.mainloop()


#def gui_load_spectral_network(plot_on_cylinder=False):
#    root = tk.Tk()
#    dir_opts = {
#        'initialdir': os.curdir,
#        'mustexist': False,
#        'parent': root,
#        'title': 'Select a directory that contains data files.',
#    }
#    data_dir = tkFileDialog.askdirectory(**dir_opts)
#    root.destroy()
#    if data_dir == '':
#        return None
#    else:
#        logging.info('Opening data directory "{}"...'.format(data_dir))
#        return load_spectral_network(data_dir, plot_on_cylinder)
    
