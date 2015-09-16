import os
import logging
import Tkinter as tk
import tkFileDialog
import pdb

from api import (generate_spectral_network, load_config, load_spectral_network,
                 save_config, save_spectral_network, 
                 make_spectral_network_plot, SpectralNetworkData,)
from trivialization import SWDataWithTrivialization

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
        self.sw_data = None
        self.spectral_networks = spectral_networks


    def create_widgets(self):
        # Layout variables
        grid_row = 0
        grid_col = 0

        # Menu
        self.mb = tk.Menubutton(self.root, text='File', relief=tk.RAISED)
        self.mb.grid(row=grid_row, column=grid_col, sticky=tk.W)
        self.mb.menu = tk.Menu(self.mb, tearoff=0,)
        self.mb['menu'] = self.mb.menu
        self.mb.menu.add_command(
            label='Load configuration',
            command=self.menu_load_config_action,
        )
        self.mb.menu.add_command(
            label='Save configuration',
            command=self.menu_save_config_action,
        )
        self.mb.menu.add_command(
            label='Load data',
            command=self.menu_load_data_action,
        )
        self.mb.menu.add_command(
            label='Save data',
            command=self.menu_save_data_action,
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
        self.entry_phase.set('1.0')
        self.entry['phase'] = tk.Entry(
            self.root,
            textvariable=self.entry_phase
        )

        ### Entry & Label layout
        grid_row += 1
        grid_col = 0
        tk.Label(self.root,
                 text='casimir_diffenrentials').grid(row=grid_row, 
                                                     column=grid_col)
        grid_col += 1
        self.entry['casimir_differentials'].grid(
            row=grid_row, column=grid_col, columnspan=3, sticky=tk.EW
        )

        grid_row += 1
        grid_col = 0
        tk.Label(self.root,
                 text='root_system').grid(row=grid_row, column=grid_col)
        grid_col += 1
        self.entry['root_system'].grid(row=grid_row, column=grid_col)
        grid_col += 1
        tk.Label(self.root,
                 text='representation').grid(row=grid_row, column=grid_col)
        grid_col += 1
        self.entry['representation'].grid(row=grid_row, column=grid_col)

        grid_row += 1
        grid_col = 0
        (tk.Label(self.root, text='differential_parameters')
         .grid(row=grid_row, column=grid_col))
        grid_col += 1
        self.entry['differential_parameters'].grid(
            row=grid_row, column=grid_col, columnspan=3, sticky=tk.EW,
        )

        for option in ['mt_params', 'punctures', 'plot_range',
                       'num_of_steps', 'num_of_iterations',
                       'size_of_small_step', 'size_of_large_step',
                       'size_of_neighborhood', 'size_of_puncture_cutoff',
                       'size_of_ramification_pt_cutoff',
                       'size_of_bin', 'accuracy', 'n_processes',
                       'mass_limit']:
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

        ### Check plot_on_cylinder
        grid_row += 1
        grid_col = 0
        self.check['plot_on_cylinder'] = tk.IntVar()
        tk.Checkbutton(
            self.root,
            text='Plot on cylinder',
            variable=self.check['plot_on_cylinder']
        ).grid(row=grid_row, column=grid_col)

        ### Set default parameters.
        for option in self.entry_var:
            self.entry_var[option].set(self.config[option])

        ### 'Generate' button
        grid_row += 1
        grid_col = 0
        self.button_generate = tk.Button(
            self.root,
            text='Generate',
            command=self.button_generate_action,
        )
        grid_col += 1
        self.button_generate.grid(row=grid_row, column=grid_col, sticky=tk.E)

        ### 'Plot' button
        grid_col += 1
        self.button_plot = tk.Button(
            self.root,
            text='Plot',
            command=self.button_plot_action,
        )
        grid_col += 1
        self.button_plot.grid(row=grid_row, column=grid_col, sticky=tk.E)


    def button_print_config_action(self):
        self.update_config_from_entries()
        print("config in parser:")
        for section in self.config.parser.sections():
            print("section: {}".format(section))
            for option in self.config.parser.options(section):
                print(
                    "{} = {}".format(
                        option, self.config.parser.getstr(section, option)
                    )
                )

        print("config :")
        for option, value in self.config.iteritems():
            print("{} = {}".format(option, value))


    def check_plot_on_cylinder(self):
        check = self.check['plot_on_cylinder'].get()
        if check == 1:
            return True
        else:
            return False


    def menu_load_config_action(self):
        config = load_config()
        if config is None:
            return None
        else:
            self.config = config

        for option, value in self.config.iteritems():
            try:
                self.entry_var[option].set(value)
            except KeyError:
                logging.warning('No entry for "{}".'.format(option))
                pass


    def menu_save_config_action(self):
        self.update_config_from_entries()
        save_config(self.config)


    def menu_load_data_action(self):
        config, spectral_network_data = load_spectral_network()
        if config is None:
            return None
        self.config = config 
        self.spectral_networks = spectral_network_data.spectral_networks
        self.sw_data = spectral_network_data.sw_data
        logging.info('Finished loading spectral network data.')


    def menu_save_data_action(self):
        dir_opts = {
            'initialdir': os.curdir,
            'mustexist': False,
            'parent': self.root,
            'title': 'Select a directory to save data files.',
        }
        data_dir = tkFileDialog.askdirectory(**dir_opts)
        if data_dir == '':
            return None
        else:
            self.update_config_from_entries()
            save_spectral_network(
                self.config, 
                SpectralNetworkData(self.sw_data, self.spectral_networks,),
                data_dir,
                make_zipped_file=False,
            )
        return None


    def update_config_from_entries(self):
        ### Read config options from Entries.
        for section in self.config.parser.sections():
            for option in self.config.parser.options(section):
                value = self.entry_var[option].get()
                if (section == 'numerical parameters'):
                    self.config[option] = eval(value)
                else:
                    self.config[option] = value
                self.config.parser.set(section, option, value)


    def button_generate_action(self):
        self.update_config_from_entries()

        spectral_network_data = generate_spectral_network(
            self.config,
            phase=eval(self.entry_phase.get()),
        )
        self.sw_data = spectral_network_data.sw_data
        self.spectral_networks = spectral_network_data.spectral_networks

        return None

    def button_plot_action(self):
        self.update_config_from_entries()

        if (len(self.spectral_networks) > 0):
            snd = SpectralNetworkData(self.sw_data, self.spectral_networks)
            make_spectral_network_plot(
                snd,
                master=self.root,
                plot_on_cylinder=self.check_plot_on_cylinder,
                plot_range=self.config['plot_range'],
            )
            return None
        else:
            logging.warning('No spectral network to plot.')
            return None


def open_gui(config, spectral_networks,):

    gui_loom = GUILoom(config, spectral_networks=spectral_networks)
    gui_loom.create_widgets()
    gui_loom.root.protocol("WM_DELETE_WINDOW", lambda: quit_gui(gui_loom),)
    gui_loom.root.mainloop()

    return None 

def quit_gui(gui_loom):
    gui_loom.root.destroy()
