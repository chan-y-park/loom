# loom
Python program to generate and analyze spectral networks

## To Do List
* show the label when placing a mouse cursor on an object of a plot.
  * https://github.com/joferkington/mpldatacursor
* embedding matplotlib into tkinter
  * see matplotlib user_interfaces examples.
  * put Tkinter.Scrollbar & Tkinter.Entry
* open a file dialog when saving data
* implement an option to change between the intersection finding routines. 
* logging to file
* GUI
* egg?
* print time after each child process finishes its job
* SciPy warnings
  * ```/usr/local/lib/python2.7/dist-packages/scipy/optimize/zeros.py:150: RuntimeWarning: Tolerance of 0.000513046300877562 reached
  warnings.warn(msg, RuntimeWarning)```
  * divide by zero
* API for getting D-type joints using other method than the Z_2 projection.

## How to run 
* Python > 2.7.6
* Not tested on Python 3.

### Python library requirement
* NumPy > 1.8.2
* SciPy > 0.15.1
* SymPy > 0.7.6
* Matplotlib > 1.4.1

### Linux

#### Scientific Linux
  1. Install Tcl/Tk and their -devel packages
  1. BLAS/LAPACK/ATLAS (we did this after installing Python but I think it makes more sense if we do it before the build & the installation of Python)
  1. Build Python 2.7.6 and install it (possibly need to change CPPFLAGS variables to include header files for Tcl/Tk)
  1. Install pip by getting its Python source code and run it using Python (http://pip.readthedocs.org/en/latest/installing.html)
  1. Install nose (`sudo pip install nose`)
  1. NumPy
  1. SciPy
  1. Matplotlib
  1. SymPy

#### Ubuntu 14.04

1. from Ubuntu package
  1. ```sudo apt-get install python python-dev libatlas-base-dev gcc gfortran g++ python-numpy python-matplotlib ipython ipython-notebook python-pandas python-nose```
1. install SciPy 0.15.1
  1. Get the source code from http://sourceforge.net/projects/scipy/files/scipy/
  1. Unpack ```scipy-<version>.tar.gz```, change to the ```scipy-<version>/``` directory.
  1. ```python setup.py install```
  1. see http://www.scipy.org/install.html for additional help.

1. install SymPy 0.7.6
  1. Get the source code from https://github.com/sympy/sympy/releases 
  1. Unpack ```sympy-<version>.tar.gz```, change to the ```sympy-<version>/``` directory.
  1. ```sudo python setup.py install``` 

