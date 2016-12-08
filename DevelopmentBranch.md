# Introduction #

The main development branch (i.e. the tar balls one can download and the corresponding git "master" branch) is supposed to be  stable and usable. However, in order to test and introduce new features and improves old ones (while possibly breaking immediate backwards compatibility) before a new release we are doing ongoing development on the **development branch** (as described under [[DevelopmentWorkflow]]).

Everyone is more than welcome to use the development version but you should be aware that you are then first and foremost a alpha-tester: Things may break in horrible ways (but your feedback through the [mailing list](http://groups.google.com/group/mdnalysis-discussion) and the [Issue Tracker](https://github.com/MDAnalysis/mdanalysis/issues) is much appreciated). (That said, many of the core developers are using the development code for their daily work.)

If you find yourself using the development branch a lot then consider signing up for the [developer mailing list](http://groups.google.com/group/mdnalysis-devel) where new features and testing of the development code are primarily discussed.

A quick overview over the latest additions to the development branch is provided by the list of [commits to the develop branch](https://github.com/MDAnalysis/mdanalysis/commits/develop).


# Getting the development branch #

The following instructions are for a read-only version of the source code (see [[git]] for further details on the `git clone` command and and [[DistributedDevelopment]] for more help):
```
git clone https://github.com/MDAnalysis/mdanalysis.git
cd mdanalysis
git fetch --all
git checkout -b develop origin/develop              # this should switch you to the development branch
```
You only have to carry out the above commands once.

Let's check that you actually have switched to the development branch: The command
```
git branch
```
should show something similar to
```
* develop
  master
```
The asterisk in front of _develop_ indicates that you have checked out the development version of MDAnalysis.

You can update it to the very [latest version](https://github.com/MDAnalysis/mdanalysis/tree/develop) with a simple
```
git pull
```

If you want to use the latest officially released version of the source code then switch to the _master_ branch:
```
git checkout master
```
You can go back to the development version with
```
git checkout develop
```


# Installing #
The MDAnalysis [source code repository](Source) contains the actual MDAnalysis library and test code. When using the development branch it is a good idea to install both so that you can always [run the tests](UnitTests) and verify that the current version of the code you're using has not been broken inadvertently.

The following uses `python setup.py install` to mean "install the package in whichever way you normally do". (In particular, as a developer `python setup.py develop` is often a very handy way to work on the code as it it takes changes to the code immediately into account.)

## Pre-requisites for development ##
In order to build from develop you need
* working C/C++ compiler
* Python 2.7.x
* Python packages already installed
  * `cython`
  * `numpy`

## Install the MDAnalysis library ##
The library lives in the [mdanalysis/package](https://github.com/MDAnalysis/mdanalysis/tree/develop/package) directory.
```
cd package
python setup.py install
```
(Ideally, this should also install additional dependencies but you might want to install numpy, scipy, biopython, netcdf4/hdf5 libraries through your operating systems package management system or use something like the [Enthought Python Distribution](https://www.enthought.com/products/epd/).)

## Install the test cases ##
The UnitTests are stored together with test data in [mdanalysis/testsuite](https://github.com/MDAnalysis/mdanalysis/tree/develop/testsuite).
```
cd testsuite
python setup.py install
```

# Running tests #
More details can be found under [[UnitTests]]. In order to run _all tests_ (which takes 5-10 Minutes on a single core) you would do (after installing library and test cases):
```
cd testsuite/MDAnalysisTests
./mda_nosetests -v
```

To run individual tests (e.g. you're only interested in the latest addition of the AMBER netcdf trajectory format)
```
./mda_nosetests -v test_coordinates:TestNCDFReader  test_coordinates:TestTRJReader test_coordinates:TestNCDFWriter
```

You want to see all tests pass and no ERROR or FAIL (unless the test is marked with "knownfailure").

You can report any problems encountered through the mailing list (preferrably the ones for development) or the [Issue Tracker](https://github.com/MDAnalysis/mdanalysis/issues).
