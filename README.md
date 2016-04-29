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
* http://het-math2.physics.rutgers.edu/dev_loom/ (developmental, beta version)

### Note for users and developers
``stable_*`` is the branch to use for the study of spectral networks, ``master`` is a developmental branch that may contain up-to-date but unstable features.

### Screenshots

#### Web UI
![config](https://github.com/chan-y-park/loom/blob/master/screeenshots/web_ui_config.png "configuration page")
![progress](https://github.com/chan-y-park/loom/blob/master/screeenshots/web_ui_progress.png "showing progress")
![plot](https://github.com/chan-y-park/loom/blob/master/screeenshots/web_ui_plot.png "plot page")


#### GUI
![loom screenshot](https://github.com/chan-y-park/loom/blob/master/screeenshots/loom_desktop.png "loom desktop")
![loom menu screenshot](https://github.com/chan-y-park/loom/blob/master/screeenshots/loom_menu.png "loom menu")
![loom plot screenshot](https://github.com/chan-y-park/loom/blob/master/screeenshots/loom_plot.png "loom plot")

## How to use ``loom`` web UI

* The front page of the web UI of ``loom`` can be found at
  * http://het-math2.physics.rutgers.edu/loom/ (stable)
  * http://chan.physics.rutgers.edu/loom/ (developmental, alpha version)
  * http://het-math2.physics.rutgers.edu/dev_loom/ (developmental, beta version)
* In the following we will denote the root url as ``/``, that is, the full url of ``/config`` is http://het-math2.physics.rutgers.edu/dev_loom/config.

### Configuration page
* Go to ``/config`` to set a configuration to run ``loom`` on the web.
* The first row specifies the Lie algebra and its representation associated with the spectral network to generate. The default configuration is associated with the standard representation of ``A_2``. For example, to generate an ``D_4`` spectral network in a spinor representation, choose the options in the following way:
  * **Lie algebra type** : D
  * **rank**: 4
  * **representation**: 4
* **Description**: a text to describe the setup; it has no effect in the run of ``loom`` but just a convenient note to understand the data without loading it onto ``loom``.
* **Casimir differentials**: a ``Python`` dictionary of the differentials that specify a class S theory. The default configuration corresponds to
```phi_2(z) = v_2 (dz)^2```
```phi_3(z) = (v_3 + z^2) (dz)^3```
where ``v_2`` and ``v_3`` are complex parameters whose numerical values are to be specified in **Parameters of differentials**. The differentials are intepreted as ``SymPy`` expressions, therefore they should follow the syntax of ``SymPy``. For example, an explicit representation of multiplication using ``*`` is necessary. 
* **Parameters of differentials**: a ``Python`` dictionary of the complex parameters that appear in the definition of the differentials. But for the convienience of users an unusual format of using ``=`` instead of ``:`` is allowed, i.e. ``{v_2 = 1.0, v_3 = 1.0}`` works the same as ``{v_2: 1.0, v_3: 1.0}``.
* **Regular punctures** and **Irregular punctures**: the locations of regular punctures, i.e. a puncture without a branch cut attached to it and the locations of irregular punctures, i.e. a puncture with a branch cut attached to it.  They are ``Python`` lists, and the entries are interpreted as ``SymPy`` expressions. Note that the locations should be the same as those specified by the differentials, and changing just the entries here does not change the spectral network itself but results in a spectral network generated in an incorrect way. In the future this information will be automatically obtained from the differentials.
* **Plot rage**: the initial range of spectral network plot to draw. You can change the range interactively in the plot page.
* **Number of steps**: each S-wall of a spectral network is a collection of numerical points evaluted by solving a differential equation, and this specifies the number of points. The larger this number, the longer it takes to finish one iteration.
* **Number of iterations**: at the start of each iteration ``loom`` grows S-walls from their seed, and at the end of each iteration ``loom`` finds the joint of S-walls and the seed of new S-walls from the joints. This parameter specifies how many iterations ``loom`` will go through, therefore a larger number of iterations will result in a longer run of ``loom``. **It's a good practice to start with ``Number of iterations = 1``**, so that no S-walls from the joints will be grown to have a rough picture of a spectral network, because there is a possibility that the number of joints can be huge even after just one iteration and the running time can be uncontrollable.
* **Mass limit**: controls the length of each S-wall. An S-wall has an associated "mass", which is the integral of the Seiberg-Witten differenital along the S-wall. When growing an S-wall it is truncated when it reaches the mass limit, even if the number of points is less then the number of steps. Therefore another way of controlling the running of ``loom`` is **starting with a small mass limit**.
* **Phase**: specifies the phases of spectral networks to generate.
  * To obtain a single spectral network, specify its phase as a real number in radian. 
  * To get multiple spectral networks, input ``[theta_i, theta_f, theta_n]``, which results in ``theta_n`` number of spectral networks from ``theta_i`` to but not including ``theta_f``. That is, ``[0, 3.14, 4]`` generates spectral networks at ``0``, ``3.14/4``, ``3.14/2``, and ``3.14*(3/4)``.
  * The default configuration generates 8 spectral networks between ``\theta = 0.01`` to ``\theta = 3.14``.
* By clicking **Save configuration** button you can save the configuration specified in this page to an ``.ini`` file on your local machine.
* To load the saved configuration, first click **Browse...** button to select the file from your local machine, and click **Load configuration** to apply the configuration, which will show you the configuration from the file on this page.
* Click **Generate spectral networks** button to start running ``loom`` with the specified configuration.
* By clicking **Show/hide advanced options** you can change more options.
  * **Mobius transformation** specifies an ``SL(2, \mathbb{C})`` transformation on the UV curve. **However it is not fully implemented yet**.
  * **Ramification point finding methods** specifies how to find ramification points of the IR curve, therefore the branch points on the UV curve. There are currently four methods available, ``system_of_eqs``, ``discriminant``, ``from_branch_points``, and ``manual``.
    * If none is specifies the default option is ``system_of_eqs``. For the detail of the methods please see [geometry.py](loom/geometry.py).
    * When using ``from_branch_points``, you need to provide the locations of the branch points on the z-plane as a ``Python`` list of ``SymPy`` expressions in **Branch points**, whose example is given in the default configuration. When both ``system_of_eqs`` and ``discriminant`` methods fail to give answers, use this method as this is more robust than those two because the locations of branch points are manually provided by you, but more convenient than ``manual`` method because you don't have to find the ``x`` roots by yourself.
    * When using ``manual``, you need to provide the locations of the ramification points in **Ramification points** as a ``Python`` list of ``(z_i, x_i)``, each of which is a ``sympy`` expression. This is the most robust method, simply because you already did all the calculations for finding the locations of the ramification points and all ``loom`` does is finding the ramification point index, or the multiplicity of the ``x`` root over a branch point.
  * When **Ramification points** is not ``None``, the ramification point finding method is automatically set to ``manual``.
  * When **Branch points** is not ``None``, the ramification point finding method is automatically set to ``from_branch_points``.
  * **Size of a small step** and **Size of a large step** are two step sizes that ``loom`` adpatively chooses as the size of one step when solving the differential equation for S-walls.
  * **Size of a branch point neighborhood** and **Size of a puncture cutoff** specify the radius around the loci that ``loom`` will stop evaluating the differential equation as the points are singular for the differential equation.
  * **Accuracy** is a cutoff below which a numerical value is considered to be zero.

### Plot page

* Place the mouse cursor onto an arrow, a cross, or a circle to display the information of the S-wall, the branch point, or the puncture, respectively. There will be a tooltip showing the label of the object and the root associated to it.
* Weights and roots are shown in an orthonomal basis of the associated Lie algebra.
* Use the mouse wheel to zoom in or out the plot, drag the plot to move it.
* On the upper right corner of the plot there are tools to save the plot in ``.png``, reset the zoon and the displacement, etc.
* On the right side of the plot, the phase of the displayed spectral network is shown. Below the phase are buttons to change the display of the plot.
  * When you zoomed in and lost some arrows on an S-wall, click **Redraw arrows** to redraw arrows within the plot area so that the corresponding tooltips can be shown by placing the mouse cursor on them.
  * Use **Show data points** and **Hide data points** to display/hide the acutall numerical points on S-walls instead of lines connecting them.
  * **Rotate back** makes the plot and the data to be rotated to the poisition where the calculation of ``loom`` is actually done, this is mostly for the purpose of debugging.
* When multiple spectral networks are displayed, use the slider under the plot to change between spectral networks.
* Extending spectral networks
  * Input **Additional steps**, **New mass limit**, or **Additional iterations** to extend spectral networks by those additional parameters.
  * Input **Additional phases** to add more spectral networks of different phases. The data can be one of the following form:
    * a number, such as ``1.0``, in radian.
    * a ``Python`` list, such as ``[0.01, 3.14, 10]``, which specifies the start, the end, and the number of additional phases.
    * a ``Python`` dict, such as ``{'single': [1.0, 2.0], 'range': [0.01, 3.14, 10], [1.53, 1.55, 10]}``, where ``single`` is a ``Python`` ``dict`` of single phases and ``range`` is a ``Python`` list of ``Python`` lists, each of which specifies a phase range. ``loom`` will automatically generate a list of phases out of this ``dict`` and remove any duplicate, so you don't have to worry about excluding duplicate phases.
  * Click **Extend** to start running ``loom`` for the extension. 
* Saving data and plot
  * To save the data on the server, first specify the name of the data in the text box on the right of **Save data to server as**, then click the **Save** button. The data can be retrieved by going to ``\plot?data=data_name`` for data named as ``data_name``.
  * Use **Download data** to save the ``loom`` data on your local machine, but you need ``loom`` to load the data.
  * Use **Download plot** to save the plot in HTML on your local machine, you don't need to install ``loom`` on your local machine to see this plot as long as it is connected to the internet.

## How to run ``loom``
``loom`` is expected to run on a Linux system. Although ``loom`` does not require any platform-specific library and therefore should run on any platform that runs Python 2.X in principle, it usually runs into a platform-specific problem when run on a platform other than Linux; for example, on Windows it has an issue with multiple processes, and on Mac it has an issue with TKinter GUI.

To run ``loom``'s web frontend, the preparation is more involved. Please contact https://github.com/chan-y-park for the detail.

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
* The most well-maintained but user-friendly UI is an web UI. You can launch a local web frontend by running ``webmain.wsgi``.
  * You need [Flask](http://flask.pocoo.org/) to run the web UI, which can be easily installed via any usual way of installing a ``Python`` package.
  * When you execute ``webmain.wsgi`` there will be a local web server running and listening to port ``8888``. Open your favourite web browser and go to ``localhost:8888/config`` for the [configuration page](#configuration-page).
  * Sometimes the port is occupied by another program. In such a case, run ``webmain.wsgi -p xxxx`` where ``xxxx`` is a port number, for example ``9999``. Then go to ``localhost:xxxx/config``.
  * When running a web UI and in the middle of a run ``loom`` crashes, there can still be a ``Python`` process running. You need to manually kill the process by first finding its ``pid`` using for example ``ps -aux`` and ``kill -9 xxxxx`` where ``xxxxx`` is the ``pid`` of the ``Python`` process.
  * Using ``pdb`` while running ``loom`` using the web UI is not recommended, as it usually won't work unless you carefully place ``pdb.set_trace()`` in an appropriate location. The reason is, as all UI frontends do, ``loom`` runs as a child process behind a parent process driving the ``UI``, and using ``pdb`` on ``loom`` only stops the ``loom`` process, which the parent process sees as the child not responding to itself.
* The most debugging-friendly way and also the most convenient way of running ``loom`` for whom knows about ``Python`` is running it using ``IPython`` and ``Jupyter``. As an example, see [how_to_run_loom.ipynb](how_to_run_loom.ipynb) for the ``Python`` code and [how_to_run_loom.html](http://het-math2.physics.rutgers.edu/loom_docs/how_to_run_loom.html) how it looks like when it is run successfully.
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
