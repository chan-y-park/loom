import os
import logging
import tkinter as tk
import tkFileDialog
import multiprocessing

from io import StringIO
from queue import Empty as QueueEmpty

from api import (set_logging, get_logging_handler,
                 generate_spectral_network, load_config, load_spectral_network,
                 save_config,
                 make_spectral_network_plot, SpectralNetworkData,)

GUI_LOOP_DELAY = 100    # in millisec
LOGGING_FILE_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    '..',
    'logs/gui_loom.log',
)


class GUILoom:
    def __init__(self, config_file_path, logging_level,):
        root = tk.Tk()
        root.wm_title('loom')

        # Put the window at the center.
        ws = root.winfo_screenwidth()
        hs = root.winfo_screenheight()
        root.geometry('+{}+{}'.format(ws / 2, hs / 2))

        # Tkinter root
        self.root = root

        # Logging handling
        self.logging_queue = multiprocessing.Queue()
        self.logger_name = 'gui_loom'
        set_logging(
            logger_name=self.logger_name,
            logging_level=logging_level,
            logging_queue=self.logging_queue,
            logging_file_name=LOGGING_FILE_PATH,
        )
        self.logging_stream = StringIO()
        self.logging_stream_handler = get_logging_handler(
            logging_level,
            logging.StreamHandler,
            self.logging_stream,
        )

        # Load default config.
        self.config = load_config(config_file_path,
                                  logger_name=self.logger_name)

        self.entry = {}
        self.entry_var = {}
        self.mb = None
        self.check = {}
        self.button = {}
        self.sw_data = None
        self.spectral_networks = []

        # Array of config options.
        # Entries that will be placed in the same row
        # are in the same row of this array.
        self.entry_array = [
            ['description'],
            ['root_system', 'representation'],
            ['casimir_differentials'],
            ['differential_parameters'],
            ['regular_punctures'],
            ['irregular_punctures'],
            ['branch_points'],
            ['ramification_points'],
            ['ramification_point_finding_method'],
            ['mt_params'],
            ['plot_range'],
            ['num_of_steps'],
            ['num_of_iterations'],
            ['size_of_small_step'],
            ['size_of_large_step'],
            ['size_of_bp_neighborhood'],
            ['size_of_puncture_cutoff'],
            ['accuracy'],
            ['mass_limit'],
            ['phase'],
        ]

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

        # Entry & Label layout
        # TODO: display config file name.
        for entry_row in self.entry_array:
            grid_row += 1
            for column, config_option in enumerate(entry_row):
                grid_col = column * 2
                entry_label_text = self.config.get_label(config_option)
                label = tk.Label(self.root, text=entry_label_text)
                label.grid(row=grid_row, column=grid_col)

                self.entry_var[config_option] = tk.StringVar()
                self.entry[config_option] = tk.Entry(
                    self.root,
                    textvariable=self.entry_var[config_option]
                )

                grid_col += 1
                if config_option in [
                    'description',
                    'casimir_differentials',
                    'differential_parameters',
                    'regular_punctures',
                    'irregular_punctures',
                    'branch_points',
                    'ramification_points',
                ]:
                    self.entry[config_option].grid(
                        row=grid_row, column=grid_col, columnspan=3,
                        sticky=tk.EW,
                    )
                else:
                    self.entry[config_option].grid(
                        row=grid_row, column=grid_col,
                    )

        # n_processes is not a config option, treated separately here.
        self.entry_n_processes_var = tk.StringVar()
        self.entry_n_processes_var.set('4')
        self.entry_n_processes = tk.Entry(
            self.root,
            textvariable=self.entry_n_processes_var
        )
        grid_col += 1
        tk.Label(self.root,
                 text='n_processes').grid(row=grid_row, column=grid_col)
        grid_col += 1
        self.entry_n_processes.grid(row=grid_row, column=grid_col)

        self.update_entries_from_config()

        # Check plot_on_cylinder
#        grid_row += 1
#        grid_col = 0
#        self.check['plot_on_cylinder'] = tk.IntVar()
#        tk.Checkbutton(
#            self.root,
#            text='Plot on cylinder',
#            variable=self.check['plot_on_cylinder']
#        ).grid(row=grid_row, column=grid_col)

        # 'Generate' button
        grid_row += 1
        grid_col = 0
        self.button['generate'] = tk.Button(
            self.root,
            text='Generate',
            command=self.button_generate_action,
        )
        grid_col += 1
        self.button['generate'].grid(
            row=grid_row, column=grid_col, sticky=tk.E
        )

        # Checkbutton for plot_on_cylinder
        # TODO: implement this functionality.
        self.check['plot_on_cylinder'] = tk.IntVar()
        grid_col += 1
        tk.Checkbutton(
            self.root,
            text='Plot on cylinder',
            variable=self.check['plot_on_cylinder']
        ).grid(row=grid_row, column=grid_col)

        # 'Plot' button
        self.button['plot'] = tk.Button(
            self.root,
            text='Plot',
            command=self.button_plot_action,
        )
        grid_col += 1
        self.button['plot'].grid(row=grid_row, column=grid_col, sticky=tk.E)

        # Log text
        grid_row += 1
        grid_col = 0
        self.log_text = tk.Text(self.root)
        self.log_text.config(
            # spacing1=2,
            spacing2=2,
            spacing3=2,
            height=12,
        )
        self.log_text.grid(
            row=grid_row, column=grid_col, columnspan=4, sticky=tk.EW
        )
        grid_col = 4
        log_text_scroll = tk.Scrollbar(self.root)
        log_text_scroll.grid(row=grid_row, column=grid_col, sticky=tk.NS)
        log_text_scroll.config(command=self.log_text.yview)
        self.log_text.config(yscrollcommand=log_text_scroll.set)

        self.root.after(GUI_LOOP_DELAY, self.get_log)
        return None

    def get_log(self):
        try:
            record = self.logging_queue.get_nowait()
            if record is not None:
                self.logging_stream_handler.handle(record)
                logs = self.logging_stream.getvalue()
                self.logging_stream.truncate(0)
                self.log_text.insert(tk.END, logs)
                self.log_text.see(tk.END)
        except QueueEmpty:
            pass
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            import sys
            import traceback
            print >> sys.stderr, 'logging_listener_process:'
            traceback.print_exc(file=sys.stderr)
        self.root.after(GUI_LOOP_DELAY, self.get_log)

    def menu_load_config_action(self):
        self.change_gui_state('off')
        self.log_text.delete('1.0', tk.END)

        # toplevel = tk.Toplevel()
        file_opts = {
            'defaultextension': '.ini',
            'initialdir': os.curdir,
            'initialfile': 'config.ini',
            # 'parent': toplevel,
            'parent': self.root,
            'title': 'Select a configuration file to load.',
        }
        file_path = tkFileDialog.askopenfilename(**file_opts)
        # toplevel.destroy()
        if file_path == '' or file_path is None:
            self.change_gui_state('on')
            return None

        config = load_config(file_path, logger_name=self.logger_name)
        if config is None:
            self.change_gui_state('on')
            return None
        self.config = config
        self.update_entries_from_config()
        self.change_gui_state('on')
        return None

    def menu_save_config_action(self):
        self.change_gui_state('off')

        # toplevel = tk.Toplevel()
        file_opts = {
            'defaultextension': '.ini',
            'initialdir': os.curdir,
            'initialfile': 'config.ini',
            # 'parent': toplevel,
            'parent': self.root,
            'title': 'Save the current configuration to a file.',
        }
        file_path = tkFileDialog.asksaveasfilename(**file_opts)
        # toplevel.destroy()
        if file_path == '' or file_path is None:
            self.change_gui_state('on')
            return None

        self.change_gui_state('off')
        self.update_config_from_entries()
        save_config(self.config, file_path=file_path,
                    logger_name=self.logger_name,)
        self.change_gui_state('on')
        return None

    def menu_load_data_action(self):
        self.change_gui_state('off')
        self.log_text.delete('1.0', tk.END)
        result_queue = multiprocessing.Queue()

        # toplevel = tk.Toplevel()
        dir_opts = {
            'initialdir': os.curdir,
            'mustexist': False,
            # 'parent': toplevel,
            'parent': self.root,
            'title': 'Select a directory that contains data files.',
        }
        data_dir = tkFileDialog.askdirectory(**dir_opts)
        # toplevel.destroy()
        if data_dir == '' or data_dir is None:
            self.change_gui_state('on')
            return (None, None)

        load_data_process = multiprocessing.Process(
            target=load_spectral_network,
            kwargs=dict(
                data_dir=data_dir,
                result_queue=result_queue,
                logger_name=self.logger_name,
            )
        )
        load_data_process.start()
        self.root.after(
            GUI_LOOP_DELAY,
            self.finish_load_data,
            load_data_process,
            result_queue
        )
        return None

    def finish_load_data(self, *args):
        logger = logging.getLogger(self.logger_name)
        load_data_process, result_queue = args
        if result_queue.empty() is True:
            if load_data_process.is_alive():
                self.root.after(GUI_LOOP_DELAY, self.finish_load_data, *args)
                return None
            else:
                logger.warning(
                    'Loading data failed: pid = {}, exitcode= {}'
                    .format(load_data_process.pid,
                            load_data_process.exitcode,)
                )
                self.change_gui_state('on')
                return None

        spectral_network_data = result_queue.get()
        load_data_process.join()

        self.config = spectral_network_data.config
        if self.config is None:
            logger.warning('Data contains no configuration.')
            self.change_gui_state('on')
            return None

        self.update_entries_from_config()
        self.spectral_networks = spectral_network_data.spectral_networks
        self.sw_data = spectral_network_data.sw_data
        logger.info('Spectral network data is ready.')
        self.change_gui_state('on')
        return None

    def menu_save_data_action(self):
        self.change_gui_state('off')

        dir_opts = {
            'initialdir': os.curdir,
            'mustexist': False,
            'parent': self.root,
            'title': 'Select a directory to save data files.',
        }
        data_dir = tkFileDialog.askdirectory(**dir_opts)
        if data_dir == '' or data_dir is None:
            self.change_gui_state('on')
            return None

        self.update_config_from_entries()
        data = SpectralNetworkData(
            sw_data=self.sw_data,
            spectral_networks=self.spectral_networks,
            config=self.config,
            logger_name=self.logger_name,
        )
        save_data_process = multiprocessing.Process(
            target=data.save,
            kwargs=dict(
                data_dir=data_dir,
                make_zipped_file=False,
            )
        )
        save_data_process.start()
        self.root.after(
            GUI_LOOP_DELAY,
            self.finish_save_data,
            save_data_process,
        )
        return None

    def finish_save_data(self, save_data_process):
        if save_data_process.is_alive():
            self.root.after(GUI_LOOP_DELAY, self.finish_save_data,
                            save_data_process,)
            return None
        else:
            save_data_process.join()
            self.change_gui_state('on')
            return None

    def button_generate_action(self):
        self.change_gui_state('off')
        self.update_config_from_entries()
        result_queue = multiprocessing.Queue()

        generate_process = multiprocessing.Process(
            target=generate_spectral_network,
            args=(
                self.config,
            ),
            kwargs=dict(
                n_processes=eval(self.entry_n_processes_var.get()),
                result_queue=result_queue,
                logger_name=self.logger_name,
            ),
        )
        generate_process.start()
        self.root.after(
            GUI_LOOP_DELAY,
            self.finish_generate,
            generate_process,
            result_queue
        )

    def finish_generate(self, generate_process, result_queue):
        logger = logging.getLogger(self.logger_name)
        if result_queue.empty() is True:
            if generate_process.is_alive():
                self.root.after(
                    GUI_LOOP_DELAY,
                    self.finish_generate,
                    generate_process,
                    result_queue
                )
                return None
            else:
                logger.warning(
                    'Generating spectral networks failed: '
                    'pid = {}, exitcode = {}.'
                    .format(generate_process.pid,
                            generate_process.exitcode,)
                )
                self.change_gui_state('on')
                return None
        spectral_network_data = result_queue.get()
        generate_process.join()
        self.sw_data = spectral_network_data.sw_data
        self.spectral_networks = spectral_network_data.spectral_networks
        logger.info('Finished generating spectral network data.')
        self.change_gui_state('on')
        return None

    def button_plot_action(self):
        logger = logging.getLogger(self.logger_name)
        self.change_gui_state('off')
        self.update_config_from_entries()

        if (len(self.spectral_networks) == 0):
            logger.warning('No spectral network to plot.')
            self.change_gui_state('on')
            return None

        snd = SpectralNetworkData(self.sw_data, self.spectral_networks)
        make_spectral_network_plot(
            snd,
            master=self.root,
            plot_on_cylinder=self.check_plot_on_cylinder,
            plot_range=self.config['plot_range'],
            logger_name=self.logger_name,
        )
        self.change_gui_state('on')

        return None

    def change_gui_state(self, state):
        if state == 'on':
            tk_state = tk.NORMAL
        elif state == 'off':
            tk_state = tk.DISABLED
        else:
            tk_state = state
        for name, button_widget in self.button.items():
            button_widget.config(state=tk_state)
        for name, entry_widget in self.entry.items():
            entry_widget.config(state=tk_state)
        self.mb.config(state=tk_state)

    # For debugging purpose only.
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
        for option, value in self.config.items():
            print("{} = {}".format(option, value))

    def check_plot_on_cylinder(self):
        check = self.check['plot_on_cylinder'].get()
        if check == 1:
            return True
        else:
            return False

    def update_entries_from_config(self):
        """
        Update the text values of the entries
        when a new configuration is loaded to self.config.
        """
        logger = logging.getLogger(self.logger_name)
        config_options = self.config.keys()
        for entry_row in self.entry_array:
            for config_option in entry_row:
                value = self.config[config_option]
                config_options.remove(config_option)
                self.entry_var[config_option].set(value)
        if len(config_options) > 0:
            logger.warning(
                'The following options are in the configuration '
                'file but has no corresponding entry in the GUI: {}.'
                .format(config_options)
            )

    def update_config_from_entries(self):
        """
        Read config options from entries.
        """
        logger = logging.getLogger(self.logger_name)
        for section in self.config.parser.sections():
            for option in self.config.parser.options(section):
                try:
                    value = self.entry_var[option].get()
                    self.config.parser.set(section, option, value)
                    if (section == 'numerical parameters'):
                        self.config[option] = eval(value)
                    else:
                        if value == 'None':
                            value = None
                        self.config[option] = value
                except KeyError:
                    logger.warning(
                        "No entry for option '{}', skip it."
                        .format(option)
                    )
                    pass


def open_gui(config_file_path=None, logging_level=None,):
    gui_loom = GUILoom(config_file_path, logging_level)
    gui_loom.create_widgets()
    gui_loom.root.protocol("WM_DELETE_WINDOW", lambda: quit_gui(gui_loom),)
    gui_loom.root.mainloop()

    return None


def quit_gui(gui_loom):
    gui_loom.root.destroy()
