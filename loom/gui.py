import os
import time
import glob
import zipfile, zlib
import logging
import Tkinter as tk
import tkFileDialog
import pdb

from math import pi
from config import LoomConfig
from api import (generate_spectral_network, load_spectral_network,)
from plotting import SpectralNetworkPlot

class GUILoom:
    def __init__(self, config=None, spectral_networks=[],):
        root = tk.Tk()
        root.wm_title('loom')
        
        # Put the window at the center
        ws = root.winfo_screenwidth()
        hs = root.winfo_screenheight()
        root.geometry('+{}+{}'.format(ws/2, hs/2))
        
        self.root = root
        self.config = config
        self.entry = {}
        self.entry_var = {} 
        self.mb = None
        self.check = {}
        self.spectral_networks = spectral_networks

    def create_widgets(self):
        # Layout variables
        grid_row = 0
        grid_col = 0

        # Menu
        self.mb = tk.Menubutton(self.root, text='File', relief=tk.RAISED)
        self.mb.grid(row=grid_row, column=grid_col)
        self.mb.menu = tk.Menu(self.mb, tearoff=0,)
        self.mb['menu'] = self.mb.menu
        self.mb.menu.add_command(
            label='Load',
            command=self.menu_load_action,
        )
        self.mb.menu.add_command(
            label='Save',
            command=self.menu_save_action,
        )

        # Associate each config option to an Entry
        for option, value in self.config.iteritems():
            self.entry_var[option] = tk.StringVar()
            self.entry_var[option].set(value)
            self.entry[option] = tk.Entry(
                self.root,
                textvariable=self.entry_var[option]
            )
        self.entry_phase = tk.StringVar()
        self.entry_phase.set('None')
        self.entry['phase'] = tk.Entry(
            self.root,
            textvariable=self.entry_phase
        )
            
        # Entry & Label layout
        grid_row += 1
        grid_col = 0
        tk.Label(self.root,
                 text='sw_curve').grid(row=grid_row, column=grid_col)
        grid_col += 1
        self.entry['sw_curve'].grid(
            row=grid_row, column=grid_col, columnspan=3, sticky=tk.EW
        )

        grid_row += 1
        grid_col = 0
        tk.Label(self.root,
                 text='sw_diff_v').grid(row=grid_row, column=grid_col)
        grid_col += 1
        self.entry['sw_diff_v'].grid(row=grid_row, column=grid_col)

        grid_col += 1
        tk.Label(self.root,
                 text='root_system').grid(row=grid_row, column=grid_col)
        grid_col += 1
        self.entry['root_system'].grid(row=grid_row, column=grid_col)

        grid_row += 1
        grid_col = 0
        tk.Label(self.root,
                 text='sw_parameters').grid(row=grid_row, column=grid_col)
        grid_col += 1
        self.entry['sw_parameters'].grid(
            row=grid_row, column=grid_col, columnspan=3, sticky=tk.EW,
        )

        for option in ['mt_params', 'punctures', 'z_range_limits',
                       'num_of_steps', 'num_of_iterations', 
                       'size_of_small_step', 'size_of_large_step',
                       'size_of_neighborhood', 'size_of_puncture_cutoff',
                       'size_of_ramification_pt_cutoff',
                       'size_of_bin', 'accuracy', 'n_processes']:
            grid_row += 1
            grid_col = 0
            tk.Label(self.root, 
                     text=option).grid(row=grid_row, column=grid_col)
            grid_col += 1
            self.entry[option].grid(row=grid_row, column=grid_col)

        grid_row += 1
        grid_col = 0
        tk.Label(self.root,
                 text='phase_range').grid(row=grid_row, column=grid_col)
        grid_col += 1
        self.entry['phase_range'].grid(row=grid_row, column=grid_col)

        grid_col += 1
        tk.Label(self.root,
                 text='phase').grid(row=grid_row, column=grid_col)
        grid_col += 1
        self.entry['phase'].grid(row=grid_row, column=grid_col)

        # Check plot_on_cylinder
        grid_row += 1
        grid_col = 0
        self.check['plot_on_cylinder'] = tk.IntVar()
        tk.Checkbutton(
            self.root,
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
            self.root, 
            text='Generate',
            command=self.button_generate_action,
        )
        grid_col += 1
        self.button_generate.grid(row=grid_row, column=grid_col, sticky=tk.E)

        # 'Plot' button
        grid_col += 1
        self.button_plot = tk.Button(
            self.root,
            text='Plot',
            command=self.button_plot_action,
        )
        grid_col += 1
        self.button_plot.grid(row=grid_row, column=grid_col, sticky=tk.E)

    def check_plot_on_cylinder(self):
        check = self.check['plot_on_cylinder'].get()
        if check == 1:
            return True
        else: 
            return False

    def menu_load_action(self):
        #t = tk.Toplevel(self.root)
        dir_opts = {
            'initialdir': os.curdir,
            'mustexist': False,
            'parent': self.root,
            'title': 'Select a directory that contains data files.',
        }
        data_dir = tkFileDialog.askdirectory(**dir_opts)
        #t.destroy()
        if data_dir == '':
            return None
        else:
            logging.info('Opening data directory "{}"...'.format(data_dir))
            self.config, self.spectral_networks = load_spectral_network(
                data_dir,
            )
            return None

    def menu_save_action(self):
        file_list = []

        # Prepare to save spectral network data to files.
        timestamp = str(int(time.time()))
        data_save_dir = os.path.join(
            self.config['root_dir'], 
            self.config['data_dir'], 
            timestamp
        )

        logging.info('Make a directory {} to save data.'.format(data_save_dir))
        os.makedirs(data_save_dir)

        # Save configuration to a file.
        config_file_name = os.path.join(data_save_dir, 'config.ini')
        logging.info('Save configuration to {}.'.format(config_file_name))
        with open(config_file_name, 'wb') as fp:
            self.config.parser.write(fp)
            file_list.append(config_file_name)

        # Save spectral network data.
        for i, spectral_network in enumerate(self.spectral_networks):
            data_file_name = os.path.join(
                data_save_dir,
                'data_{}.json'.format(
                    str(i).zfill(len(str(len(self.spectral_networks)-1)))
                )
            )
            logging.info('Saving data to {}.'.format(data_file_name))
            with open(data_file_name, 'wb') as fp:
                spectral_network.save_json_data(fp,)

        # Make a compressed data file.
        file_list += glob.glob(os.path.join(data_save_dir, 'data_*.json'))
        zipped_file_name = data_save_dir + '.zip'
        logging.info('Save compressed data to {}.'.format(zipped_file_name))
        with zipfile.ZipFile(zipped_file_name, 'w',
                             zipfile.ZIP_DEFLATED) as fp:
            for a_file in file_list:
                fp.write(a_file, os.path.relpath(a_file, data_save_dir))

        return None

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

        self.spectral_networks = generate_spectral_network(
            self.config,
            phase=eval(self.entry_phase.get()),
        )

        return None

    def button_plot_action(self):
        # Plot spectral networks.
        if (len(self.spectral_networks) > 0):
            spectral_network_plot = SpectralNetworkPlot(
                master=self.root,
                config=self.config,
                plot_on_cylinder=self.check_plot_on_cylinder(),
                #plot_data_points=True,
                #plot_joints=True,
                #plot_bins=True,
                #plot_segments=True,
            )

            for spectral_network in self.spectral_networks:
                logging.info('Generating the plot of a spectral network '
                             '@ theta = {}...'.format(spectral_network.phase))
                spectral_network_plot.draw(spectral_network)

            return spectral_network_plot.show()
        else:
            logging.warning('No spectral network to plot.')
            return None


def open_gui(config, spectral_networks,):

    gui_loom = GUILoom(config, spectral_networks=spectral_networks)
    gui_loom.create_widgets()
    gui_loom.root.mainloop()

    return None 