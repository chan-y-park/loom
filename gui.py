from Tkinter import Tk, Frame, Button

class Application(Frame):
    def __init__(self, opts, master=None):
        Frame.__init__(self, master)
        self.opts = opts
        self.pack()
        self.createWidgets()

    def createWidgets(self):
        self.QUIT = Button(self)
        self.QUIT['text'] = 'QUIT'
        self.QUIT['fg'] = 'red'
        self.QUIT['command'] = self.quit
        self.QUIT.pack({'side': 'left'})

        self.show_opts = Button(self)
        self.show_opts['text'] = 'Show options'
        self.show_opts['command'] = self.print_opts
        self.show_opts.pack({'side': 'left'})

    def print_opts(self):
        print self.opts

def open_gui(opts):
    root = Tk()
    app = Application(opts, master=root)
    app.mainloop()
    root.destroy()
