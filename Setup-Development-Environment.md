# Virtual Environments

We highly recommend using [virtual environments](https://pypi.python.org/pypi/virtualenv) with [virtualenvwrapper](http://virtualenvwrapper.readthedocs.org/en/latest/).
This allows you to have multiple experimental development versions of MDAnalysis which do not interfere with each other or your own stable version.
Since MDAnalysis is split into the actual package and a testsuite you need to install both modules in development mode.

# Checking out and installing MDAnalysis

Firstly, a copy of the raw source files needs to be downloaded from the code repository.  From the command line in a suitable directory run:

```
git clone https://github.com/MDAnalysis/mdanalysis.git
```
In the future, if you want to update your copy of the code to the latest version, this can be done by running these commands from within the source directory:

```
git fetch origin

git pull origin develop
```

To install this version of the code, run the following commands: 

```
cd <to where ever you checked out MDAnalysis>
mkvirtualenv mdanalysis
pip install numpy
pip install cython
pip install -e package/
pip install -e testsuite/
```

The `-e` flag will cause pip to call setup with the `develop` option. This means that any changes on the source code will immediately be reflected in your virtual environment. 

## macOS specific instructions
One more step is required on macOS because of the number of files that a process can open simultaneously is quite low (256). So to increase the number of files that can be accessed, run the following command:

`ulimit -n 4096`

This sets the number of files to 4096.



# Testing

* To run the unit-tests

```
cd testsuite/MDAnalysisTests
./mda_nosetests
```

* Run specific tests

```
./mda_nosetests test_analysis.py
./mda_nosetests test_analysis.py:TestContactMatrix
```

* test with coverage

```
./mda_nosetest --with-coverage --cover-erase --cover-package=MDAnalsysis
```

If you also want to get a HTML output of the coverage report

```
./mda_nosetest --with-coverage --cover-erase --cover-package=MDAnalsysis --cover-html --cover-html-dir=coverage
```
