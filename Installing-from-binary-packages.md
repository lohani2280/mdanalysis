For some platforms, MDAnalysis can be installed directly in pre-compiled form. (If no pre-compiled packages are available, please see [installing MDAnalysis from source](Install).)

# Anaconda
**Note**: This currently **experimental**, feedback welcome (at PR [#751](/MDAnalysis/mdanalysis/issues/751))
* *Experimental, since 0.14.1-dev0* (see issue [#608](/MDAnalysis/mdanalysis/issues/608) with PR [#751](/MDAnalysis/mdanalysis/issues/751)).
* Currently only Linux 64 bit supported.
* Currently only Python 2.

## Installing the anaconda2 distribution for the first time
First install *anaconda* (Python 2 for right now) by [downloading the Python 2.7 Linux 64-bit anaconda2 installer](https://www.continuum.io/downloads#_unix) and running the installer
```bash
bash Anaconda2-2.5.0-Linux-x86_64.sh
```
(adjust the version number according to the name of the downloaded installer). Unless you know better, have the installer "prepend the Anaconda2 install location to `PATH` in your ~/.bashrc" (answer *yes*). See [conda test-drive](http://conda.pydata.org/docs/test-drive.html) for an introduction to using `conda`.

Check that your conda environment is working:
```bash
conda --version
```
should give something like 
```
conda 3.19.1
```
perhaps with a higher version number.

Update everything to the latest version (answer `y` when asked)
```bash
conda update conda
```

## Installing MDAnalysis with `conda`
Add the MDAnalysis anaconda channel (only has to be done once):
```bash
conda config --add channels MDAnalysis
```

Install MDAnalysis
```bash
conda install mdanalysis
```



# Distributions
Some [third-party repositories](Downloads#binary-distributions) might have binary packages of MDAnalysis but we are not currently making any binary packages available.


