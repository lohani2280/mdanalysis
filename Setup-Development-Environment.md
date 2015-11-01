# Virtual Environments

We recommend using [virtual environments](https://pypi.python.org/pypi/virtualenv) with [virtualenvwrapper](http://virtualenvwrapper.readthedocs.org/en/latest/). Since MDAnalysis is split into the actual package and a testsuite you need to install both modules in development mode.

```
cd <to where ever you checked out MDAnalysis>
mkvirtualenv mdanalysis
pip install numpy
pip install cython
pip install -e package/
pip install -e testsuite/
```

The `-e` flag will cause pip to call setup with the `develop` option. This means that any changes on the source code will immediately be reflected in your virtual environment. 

# Testing

* To run the unit-tests

```
cd testsuite/MDAnalsysisTests
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