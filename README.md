# ``loom``
A framework to generate and analyze spectral networks, including a GUI based on ``Tkinter`` and a web frontend based on ``flask``. ``loom`` is written in Python, Javascript, and C++.

[Overview](#overview)

[How to use ``loom`` web UI](#how-to-use-loom-web-ui)
[How to run ``loom``](#how-to-run-loom)

## Overview
``loom`` has the following functionality.
* Generate a single spectral network for a given Seiberg-Witten data at a phase.
* Generate multiple spectral networks at various phases by utilizing parallel computation.
  * To generate multiple spectral networks, need to specifiy ```phase``` to ```None```.
* Save and load analytical and numerical configurations.
* Save and load generated data in JSON.
* Plot spectral networks with interactive labeling.

In addition, ``loom`` contains a web frontend that drives a WSGI application, which can be loaded from any web server that supports WSGI applications, including Apache. To see how the WSGI application looks like, visit
* http://het-math2.physics.rutgers.edu/loom/ (stable)
* http://chan.physics.rutgers.edu/loom/ (developmental, alpha version)
* http://het-math2.physics.rutgers.edu/loom/ (developmental, beta version)

### Note for users and developers
``stable_*`` is the branch to use for the study of spectral networks, ``master`` is a developmental branch that may contain up-to-date but unstable features.

### Screenshots

![loom screenshot](https://github.com/chan-y-park/loom/blob/master/screeenshots/loom_desktop.png "loom desktop")
![loom menu screenshot](https://github.com/chan-y-park/loom/blob/master/screeenshots/loom_menu.png "loom menu")
![loom plot screenshot](https://github.com/chan-y-park/loom/blob/master/screeenshots/loom_plot.png "loom plot")

## How to use ``loom`` web UI

### Configuration page
### Plot page

## How to run ``loom``
``loom`` is expected to run on a Linux system. Although ``loom`` does not require any platform-specific library and therefore should run on any platform that runs Python 2.X in principle, it usually runs into a platform-specific problem when run on a platform other than Linux; for example, on Windows it has an issue with multiple processes, and on Mac it has an issue with TKinter GUI.

To run ``loom``'s web frontend, the preparation is more involved, which will be posted here soon. In the meanwhile, please contact https://github.com/chan-y-park for the detail.

### Installation
* The following instruction is based on Ubuntu 14.04 LTS.
* Install Anaconda 2.3.0
  * https://store.continuum.io/cshop/anaconda/
* Install Sage 6.8
  * http://doc.sagemath.org/html/en/installation/index.html
* Install CGAL 4.6.2
  * http://doc.cgal.org/latest/Manual/installation.html
* Then clone this repository.

### Running ``loom``
* ```gmain.py``` is an executable Python script that launches a GUI version of ```loom```.
* In the Python interpreter import ```loom``` module.

### From Python interpreter
* Generating a spectral network data
  1. In the root directory that contains ```loom``` directory, start the python interpreter.
  1. ```>>> import loom```
  1. Load a configuration file.
    1. Load it interactively.
      * ```>>> config = loom.load_config()```
      * A file dialog window opens, select a configuration file to use.
    1. Or give the file name as an argument.
      * ```>>> config = loom.load_config('default.ini')```.
    1. The returned value is ```LoomConfig``` that contains all the configuration information.
  1. To get a K-wall network at a fixed phase, run
    * ```>>> data = loom.generate(config, phase=1.0)```
  1. To get multiple K-wall networks according to the configuration, run
    * ```>>> data = loom.generate(config)```
* Plotting a spectral network
  1. ```>>> spectral_network_plot = loom.plot(data)```
  1. According to the setup, you can either:
    * Hover the mouse cover over the object to display its label, or
    * Click an object in each figure to display its label and press ```d``` to delete all the displayed labels.
* Saving data
  1. ```>>> loom.save(config, data, data_dir='data/name/', make_zipped_file=True)```
    * This saves ```config.ini```, ```data_*.json``` files at ```data/name/``` directory, and make a zipped file of the directory.
    * If ``data_dir`` is not given, this opens a directory dialog window to select a directory to save ```config.ini``` and ```data.mose```.
* Loading the data
  1. Load it interactively.
    * ```>>> config, data = loom.load()```
    * Then you first get a directory dialog window to select a directory that contains the configuration & the data file. Select a directory and it returns ```(config, data)```.
  1. Or give the directory name as an argument.
    * ```>>> config, data = loom.load('data/name/')```.
