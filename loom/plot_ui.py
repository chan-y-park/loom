import os
import tkinter as tk
import matplotlib

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg as FigureCanvas,
    NavigationToolbar2Tk as NavigationToolbar,
)
from matplotlib.widgets import Button
from matplotlib import pyplot

import mpldatacursor

from plot_api import SpectralNetworkPlot

class SpectralNetworkPlotUIBase(SpectralNetworkPlot):
    def __init__(self, matplotlib_figure=None, plot_range=None,):
        super(SpectralNetworkPlotUIBase, self).__init__(
            matplotlib_figure=matplotlib_figure,
            plot_range=plot_range,
        )
        self.plots = []
        self.data_cursor = None
        self.current_plot_idx = None

    def draw(
        self,
        **kwargs
    ):
        rect = [.1, 0.15, .8, .8]

        axes = self.figure.add_axes(
            rect,
            label="Network #{}".format(len(self.plots)),
            aspect='equal',
        )

        super(SpectralNetworkPlotUIBase, self).draw(
            axes=axes,
            **kwargs
        )

        axes.set_visible(False)
        self.plots.append(axes)

    def set_data_cursor(self):
        if self.current_plot_idx is None:
            return None

        # Use a DataCursor to interactively display the label
        # for artists of the current axes.
        self.data_cursor = mpldatacursor.datacursor(
            axes=self.plots[self.current_plot_idx],
            formatter='{label}'.format,
            tolerance=4,
            # hover=True,
            # display='single',
            display='multiple',
            draggable=True,
        )

    def change_current_plot(self, new_plot_idx):
        if new_plot_idx == self.current_plot_idx:
            return None
        elif new_plot_idx < 0:
            new_plot_idx = 0
        elif new_plot_idx > len(self.plots) - 1:
            new_plot_idx = len(self.plots) - 1

        if self.current_plot_idx is not None:
            self.plots[self.current_plot_idx].set_visible(False)
        self.plots[new_plot_idx].set_visible(True)
        # Update the index variable for the currently displayed plot.
        self.current_plot_idx = new_plot_idx

    def show_prev_plot(self, event):
        if self.data_cursor is not None:
            self.data_cursor.hide().disable()
        self.change_current_plot(self.current_plot_idx - 1)
        self.set_data_cursor()
        self.figure.show()

    def show_next_plot(self, event):
        if self.data_cursor is not None:
            self.data_cursor.hide().disable()
        self.change_current_plot(self.current_plot_idx + 1)
        self.set_data_cursor()
        self.figure.show()


class SpectralNetworkPlotUI(SpectralNetworkPlotUIBase):
    """
    This class implements UIs using matplotlib widgets
    so that it can be backend-independent.
    """
    def __init__(
        self,
        title=None,
        plot_range=None,
    ):
        super(SpectralNetworkPlotUI, self).__init__(
            matplotlib_figure=pyplot.figure(title),
            plot_range=plot_range,
        )

        self.axes_button_prev = None
        self.axes_button_next = None
        self.index_text = None

    def save(self, plot_dir, file_prefix='', svg=False, pdf=False):
        # TODO: change the current figure to plot_id.
        digits = len(str(len(self.plots) - 1))
        for i, axes in enumerate(self.plots):
            self.change_current_plot(i)
            file_types = [".png"]
            if svg is True:
                file_types.append(".svg")
            if pdf is True:
                file_types.append(".pdf")
            for ft in file_types:
                plot_file_path = os.path.join(
                    plot_dir, file_prefix + str(i).zfill(digits) + ft
                )
                self.figure.savefig(plot_file_path)

    def change_current_plot(self, new_plot_idx):
        super(SpectralNetworkPlotUI, self).change_current_plot(new_plot_idx)
        if self.index_text is not None:
            self.index_text.set_text(
                "{}/{}".format(self.current_plot_idx, len(self.plots) - 1)
            )

    def show_prev_plot(self, event):
        super(SpectralNetworkPlotUI, self).show_prev_plot(event)

    def show_next_plot(self, event):
        super(SpectralNetworkPlotUI, self).show_next_plot(event)

    def show(self):
        plot_idx = 0
        self.current_plot_idx = plot_idx
        self.plots[plot_idx].set_visible(True)
        self.set_data_cursor()

        if(len(self.plots) > 1):
            button_width = .05
            index_width = .03 * len(str(len(self.plots) - 1))
            button_height = .05
            button_bottom = .025
            center = .5
            margin = .005

            axes_prev_rect = [center - index_width / 2 - margin - button_width,
                              button_bottom, button_width, button_height]
            axes_prev = self.figure.add_axes(axes_prev_rect)
            self.button_prev = Button(axes_prev, '<')
            self.button_prev.on_clicked(self.show_prev_plot)

            self.index_text = self.figure.text(
                center - index_width / 2, (button_bottom + button_height) / 2,
                "{}/{}".format(self.current_plot_idx, len(self.plots) - 1)
            )

            axes_next_rect = [center + index_width / 2 + margin,
                              button_bottom, button_width, button_height]
            axes_next = self.figure.add_axes(axes_next_rect)
            self.button_next = Button(axes_next, '>')
            self.button_next.on_clicked(self.show_next_plot)

        self.figure.show()


class SpectralNetworkPlotTk(SpectralNetworkPlotUIBase):
    """
    This class implements UIs using Tkinter.

    The content of this class is independent of the parent class.
    It only depends on the grandparent class, which should be
    'NetworkPlotBase'. Therefore this class can inherit any class
    whose parent is 'NetworkPlotBase'; just change the name of the
    parent in the definition of this class.
    """
    def __init__(self, master=None, title=None, plot_range=None,):
        super(SpectralNetworkPlotTk, self).__init__(
            matplotlib_figure=matplotlib.figure.Figure(),
            plot_range=plot_range,
        )

        if master is None:
            root = tk.Tk()
            root.withdraw()
            self.root = root
        else:
            self.root = master
        self.master = master

        # Create a Toplevel widget, which is a child of GUILoom
        # and contains plots,
        self.toplevel = tk.Toplevel(self.root)
        self.toplevel.wm_title(title)
        self.toplevel.protocol("WM_DELETE_WINDOW", self.on_closing)

        self.plot_idx_scale = None

        self.plot_idx_entry = None
        self.plot_idx_entry_var = tk.StringVar()
        self.plot_idx_entry_var.trace('w', self.plot_idx_entry_change)

        self.canvas = FigureCanvas(
            self.figure,
            master=self.toplevel,
            resize_callback=self.canvas_resize_callback
        )
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        toolbar = NavigationToolbar(self.canvas, self.toplevel)
        toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    def on_closing(self):
        self.toplevel.destroy()
        if self.master is None:
            self.root.destroy()

    def scale_action(self, scale_value):
        new_plot_idx = int(scale_value)
        self.update_current_plot(new_plot_idx)
        self.plot_idx_entry_var.set(new_plot_idx)

    def plot_idx_entry_change(self, *args):
        try:
            new_plot_idx = int(self.plot_idx_entry_var.get())

            if new_plot_idx == self.current_plot_idx:
                return None
            elif new_plot_idx < 0:
                new_plot_idx = 0
            elif new_plot_idx > len(self.plots) - 1:
                new_plot_idx = len(self.plots) - 1

            self.plot_idx_scale.set(new_plot_idx)
            self.update_current_plot(new_plot_idx)

        except ValueError:
            pass

        return None

    def update_current_plot(self, new_plot_idx):
        if self.data_cursor is not None:
            self.data_cursor.hide().disable()

        self.plots[self.current_plot_idx].set_visible(False)
        self.plots[new_plot_idx].set_visible(True)
        # Update the index variable for the currently displayed plot.
        self.current_plot_idx = new_plot_idx
        self.set_data_cursor()
        self.canvas.draw_idle()
        self.canvas.get_tk_widget().focus_set()

        return None

    def canvas_resize_callback(self, event):
        self.set_data_cursor()

    def show(self):
        plot_idx = 0
        self.current_plot_idx = plot_idx
        self.plots[plot_idx].set_visible(True)
        self.set_data_cursor()

        if(len(self.plots) > 1):
            tk.Label(
                self.toplevel,
                text='Plot #',
            ).pack(side=tk.LEFT)

            self.plot_idx_entry_var.set(plot_idx)
            self.plot_idx_entry = tk.Entry(
                self.toplevel,
                textvariable=self.plot_idx_entry_var,
                width=len(str(len(self.plots) - 1)),
            )
            self.plot_idx_entry.pack(side=tk.LEFT)

            tk.Label(
                self.toplevel,
                text='/{}'.format(len(self.plots) - 1),
            ).pack(side=tk.LEFT)

            self.plot_idx_scale = tk.Scale(
                self.toplevel,
                command=self.scale_action,
                # length=100*len(self.plots),
                orient=tk.HORIZONTAL,
                showvalue=0,
                to=len(self.plots) - 1,
                variable=self.current_plot_idx,
            )
            self.plot_idx_scale.pack(
                expand=True,
                fill=tk.X,
                side=tk.LEFT,
            )
