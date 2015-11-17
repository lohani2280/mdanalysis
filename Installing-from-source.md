MDAnalysis is *open source* and you can always build and install the package from source. This gives you the flexibility to use either the latest releases or, if you like living on the cutting edge, you can also use the current [development version of the code](DevelopmentBranch). This page contains notes on how to install from source.

For alternative ways to install MDAnalysis, see the page [Install](Install).

## Installing from source ##
**MDAnalysis** is also distributed in source form. In order to build the library some C code needs to be compiled. Thus you will need
  * python (>=2.7 and <3)
  * a C compiler (GNU `gcc` works),
  * the Python header files (typically in a package `python-dev`),
  * [numpy](http://numpy.scipy.org) to compile the DCD reader and numerical extensions,

See
[Additional non-standard packages](Install#Additional_non-standard_packages) below for what else you need at run
time.

The source code can be obtained via the [Downloads](Downloads) page and as described below under [Getting the source](#Getting_the_source). The primary dependency is [numpy](http://www.numpy.org/).

Please read through this whole document. Installing MDAnalysis is unfortunately not always absolutely straightforward because of the various dependencies on other packages. **[InstallRecipes](InstallRecipes)** collects a number of examples of how MDAnalysis has been successfully installed on various systems; possibly one of the recipes applies to your situation, too, and you can simply copy and paste.


### Getting the source ###
You can either get the source tar files or check out the source code from the [git](git) repository.

#### Source tar files ####
In order to download the source packages, you need to go to the PyPI repositories:
  * https://pypi.python.org/pypi/MDAnalysis (main package)
  * https://pypi.python.org/pypi/MDAnalysisTests (for the [tests](UnitTests))
and download the tar files from there.

#### Git repository ####
Alternatively, you can check out the MDAnalysis directory from the [git](git)
repository at https://github.com/MDAnalysis/mdanalysis. In most cases
simply do
```
git clone https://github.com/MDAnalysis/mdanalysis
cd mdanalysis
```

If you have cloned the repository before, all you need to do is
```
git pull
```
to update your files.

### Standard installation from source ###
If you want to install MDAnalysis from source you will need
a [C compiler](http://en.wikipedia.org/wiki/List_of_compilers#C_compilers).

We recommend (for python >= 2.7) the following way of installing in your home directory (a so called "user" installation) ::
```
python setup.py build
python setup.py install --user
```

If you want to install in the system wide python directory you will probably require administrative privileges, and so do
```
sudo python setup.py install
```

It is also possible to use `--prefix`, `--home`, or `--user` options for
`setup.py` to install in a different (probably your private) python
directory hierarchy.

If you have problem at this stage then have a look at the operating
system specific notes at the end of this file or look in the issue
tracker --- maybe the problem is recognized and a workaround can be
found in the comments

Unfortunately, installing python packages is not always completely straightforward, so please read through these docs and perhaps [ask on the mailing list](http://groups.google.com/group/mdnalysis-discussion).


#### Selecting an installation directory ####

In you should be able to use
```
 python setup.py install --prefix LOCAL_DIRECTOY
```
or any of the other options of [distutil's setup.py to install in alternative directories](http://docs.python.org/install/index.html#alternate-installation).

Please ask on the [mailing list](http://groups.google.com/group/mdnalysis-discussion) for help with installation and/or [file a bug through the issue tracker](https://github.com/MDAnalysis/mdanalysis/issues).


### pip ###
One can also use [pip](https://pip.pypa.io/en/latest/index.html). For all [releases](Release-Notes), you should be able to just do
```
pip --upgrade MDAnalysis 
pip --upgrade MDAnalysisTests
```
and `pip` will download the source distribution from the Python package index and install them, together with any required dependencies.

If you want to install from the checked out source, you can also use `pip`:
```
pip [options] ./mdanalysis
```
for a standard installation.


### Developer installation ###

A _developer installation_ (that immediately reflects changes to the
sources) can be done with
```
cd ./mdanalysis
python setup develop [options]
```

For testing one can simply use
```
python setup.py develop --install-dir=$HOME/python-lib/ --script-dir=$HOME/bin
```
This builds and installs a working version in
`~/python-lib/MDAnalysis` which is linked to the unpacked source. Then add `$HOME/python-lib` to the `PYTHONPATH`
```
export PYTHONPATH=${PYTHONPATH}:$HOME/python-lib
```
However, the developer installation above is probably cleaner.

# Additional non-standard packages #

See the operating system specific notes below for hints how to get the necessary packages through the native package management system. Please add your own (e.g. for RPM based systems, which the developers are not using heavily).

  * _python-dev_ includes `Python.h`, which is required for compiling.
  * _numpy_ is used at the compilation stage to find maths libraries.
  * _biopython_ is used for strict PDB file parsing and the KD-tree library
  * _scipy_, _netcdf4-python_ are only needed when one wants to use all MDAnalysis functions but are not required for compiling.


## Required python packages ##
In order to make _full use_ of the library in your own python code you will need at least
  * [numpy](http://numpy.scipy.org/) of version 1.5 or greater
  * [scipy](http://www.scipy.org/)
  * [matplotlib](http://matplotlib.sourceforge.net/)
  * [BioPython](http://biopython.org/)
  * [networkx](http://networkx.lanl.gov/) --- for analysis of lipid leaflets via MDAnalysis.analysis.leaflet
  * [GridDataFormats](http://pypi.python.org/pypi/GridDataFormats)
  * [netcdf4-python](http://code.google.com/p/netcdf4-python/): If you want to operate on AMBER binary trajectories (NetCDF) then you will need to have working [NetCDF](netcdf) support based on the netcdf4-python package.
    * If you get a _ValueError: did not find HDF5 headers_ during      installation then you should install the _hdf5 development package_ through your package manager (see below for some hints on how to do this for [Linux](#Linux) and [Mac OS X](#Mac_OS_X)).


## Optional python packages ##

In order to run the [UnitTests](UnitTests) you will need
  * numpy >= 1.5
  * nose >= 1.3.7

# OS specific notes #

See also InstallRecipes for "copy & paste" installation instructions for various operating systems and versions.

## Linux ##
  * Install prerequisite packages on **[Ubuntu](http://www.ubuntu.com/)** or **[Debian](http://www.debian.org/)** with

    ```
    sudo apt-get install python-dev python-cython python-numpy g++
    ```

  * Install packages needed for full functionality

    ```
    sudo apt-get install python-scipy python-matplotlib python-biopython  libhdf5-dev 
    ```

    (If you have issues with the netcdf libraries see the [wiki page on netcdf](netcdf)); the `netcdf` and `libhdf5-dev` packages are  required if you want to be able to process AMBER netcdf (binary) trajectories.)

  * Tested with GNU compilers

## Mac OS X ##
### Mac Ports ###
Tested 10.6.8+MacPorts
```
port install  py27-numpy  py27-cython
port install  py27-scipy  py27-matplotlib py27-biopython hdf5 netcdf+dap+netcdf4
```
(The `netcdf` and `libhdf5-de`v packages are only required if you want to be able to process AMBER netcdf (binary) trajectories.)

Other packages such as `networkx` and `GridDataFormats` are automagically installed when running `python setup.py install` (or `pip`).

### fink ###
Tested with Mac OS X 10.4.11+fink (**has not been tested since 2011 — feedback is welcome**)

Install prerequisite packages using [fink](http://www.finkproject.org/)
```
     fink install cython-py25 scipy-core-py25 scipy-py25
```
Install packages needed for full functionality
```
     fink install  matplotlib-py25 biopython-py25
```

**NOTE**: Use the Apple-provided `gcc`, not the fink provided one. Apparently, only Apple's gcc has the `-arch` flag that appears to be required for `python setup.py install`. (Thanks to Justin Lemkul for pointing it out; see [Issues after Installation](http://groups.google.com/group/mdnalysis-discussion/browse_thread/thread/1bf033ed1d3bb915)). In order to install Apple's gcc you will need to install the Apple Developer Tools that come on one of the additional installation disks with the Mac or install [Xcode](http://developer.apple.com/xcode/) from the web (note that Xcode 3.x is free and can be obtained by registering online for the Apple Developer Connection).



