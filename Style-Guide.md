# Requirements for code to be included into MDAnalysis
The aim of this guide is to have a uniform and well-tested code base for MDAnalysis in order to improve code quality and maintainability. The **Style Guide** *prescribes* various aspects that new code has to conform to; old code should be refactored to conform to these requirements. 

The Style Guide evolves through a community-driven process of [[proposals|https://github.com/MDAnalysis/mdanalysis/issues?q=is%3Aissue+is%3Aopen+label%3Aproposal]] and discussions among all interested developers.

1. [[Coding Style|Style-Guide#coding-style]]
1. [[Importing modules|Style-Guide#importing-modules]]
1. [[Writing Docstrings|Style-Guide#writing-docstrings]]  
1. [[Code organization|Style-Guide#code-organization]]
1. [[Version control|Style-Guide#version-control]]
1. [[Tests|Style-Guide#tests]]

(For background and history see [Issue #404](/MDAnalysis/mdanalysis/issues/404), which contains the original discussion on the Style Guide.)

## Coding Style

MDAnalysis is a project with a long history and many contributors and hasn't used a consistent coding style. Since version 0.11.0 we are trying to update all the code to conform with [PEP8](http://legacy.python.org/dev/peps/pep-0008/). Our strategy is to update the style every time we touch an old function and thus switch
to [PEP8](http://legacy.python.org/dev/peps/pep-0008/) continuously.

### Important requirements

* keep line length to **79 characters** or less; break long lines sensibly
* indent with **spaces** and use **4 spaces per level**
* naming (follows PEP8): 
  * classes: *CapitalClasses* (i.e. capitalized nouns without spaces)
  * methods and functions: *underscore_methods* (lower case, with underscores for spaces)

### Python 2/3 compatibility

MDAnalysis strives to be compatible with python 2 and 3 at the same time. To deal with the differences
we use [six](https://pypi.python.org/pypi/six) which takes care of loading the appropriate functions for each version of python. So instead of `xrange` or `iterzip` use the `zip` and `range` function provided by [six](https://pypi.python.org/pypi/six) like this.

```python
from six.moves import zip, range

for i,j in zip(range(3), range(3)):
    print(i, j)
```

### Code Linter

We recommend that you either use a Python IDE ([PyCharm](https://www.jetbrains.com/pycharm/) and others). You can also use external tools like [flake8](http://flake8.readthedocs.org/en/latest/). For integration of external tools with emacs and vim check out [elpy (emacs)](https://github.com/jorgenschaefer/elpy) and [python-mode (vim)](https://github.com/klen/python-mode).

### Code Formatter

To apply the code formatting in an automated way you can also use code formatters. As external tools there are [autopep8](https://github.com/hhatto/autopep8) and [yapf](https://github.com/google/yapf). Most IDE's either have their own code formatter or will work with one of the above through plugins.

## Importing modules
We are striving to keep module dependencies small and lightweight (i.e., easily installable with `pip`).

### General rules for importing
* Imports must all happen at the start of a module (not inside classes or functions).  
* Modules must be imported in the following order:
  - [future](https://docs.python.org/2/library/__future__.html) (`from __future__ import absolute_import, print_function, division`)
  - global imports (installed packages)
  - local imports (MDAnalysis modules)
* use **absolute imports** in the library (i.e. relative imports must be explicitly indicated), e.g.,
  ```python
  from __future__ import absolute_import

  import numpy as np

  import .core
  import ..units
  ```

#### Module imports in `MDAnalysis.analysis`

1. In `MDAnalysis.analysis`, all imports must be at the top level (as in the General Rule) — see [[#666|https://github.com/MDAnalysis/mdanalysis/issues/666]]
1. [[Optional modules|Style-Guide#modules-in-mdanalysisanalysis-and-mdanalysisvisualization]] can be imported.
1. No analysis module is imported automatically at the `MDAnalysis.analysis` level to avoid breaking the installation when optional dependencies are missing.

#### Module imports in the test suite
* In the test suite, use the order above, but import `MDAnalysis` modules before `MDAnalysisTests` imports
* Do not use *relative imports* (e.g. ``import .datafiles``) in the test suite because it breaks running the tests from inside the test directory (see [#189](/MDAnalysis/mdanalysis/issues/189))
* Never import the `MDAnalysis` module from the `__init__.py` of `MDAnalysisTests` or from any of its plugins (it's ok to import from test files). Importing `MDAnalysis` from the test setup code will cause severe coverage underestimation.


### List of core module dependencies

Any module from the standard library can be used, as well as the following nonstandard libraries:

   * `numpy`
   * `biopython`
   * `gridDataFormats`
   * `networkx`

because these packages are always installed. (Note that `scipy` is optional and not guaranteed to be installed, see below).

If you must depend on a new external package, first discuss its use on the [developer mailing list](http://developers.mdanalysis.org) or as part of the issue/PR. 


#### Modules in the "core"
The core of MDAnalysis contains all packages that are not in `MDAnalysis.analysis` or `MDAnalysis.visualization`. Only packages in the [[List of core module dependencies|Style-Guide#list-of-core-module-dependencies]] can be imported.

#### Modules in `MDAnalysis.analysis` and `MDAnalysis.visualization`
Modules under `MDAnalysis.analysis` are considered independent. Each can have its own set of dependencies. We strive to keep them small as well but module authors are in principle free to import what they need. These dependencies beyond the core dependencies are **optional dependencies** (and should be listed in `setup.py` under *analysis*).

A user who does not have a specific optional package installed must still be able to import everything else in MDAnalysis. An analysis module may simply raise an `ImportError` if a package is missing but it is recommended that it should print and log an *error message* notifying the user that a specific additional package needs to be installed to use this module.

If a large portion of the code in the module does not depend on a specific optional module then you should guard the import at the top level with a `try/except`, print and log a *warning*, and only raise an `ImportError` in the specific function or method that would depend on the missing module. (Example: see `scipy.sparse` in `MDAnalysis.analysis.distances` in PR [[#708|https://github.com/MDAnalysis/mdanalysis/pull/708]])


## Writing Docstrings

Since 0.11.0 we adopted the use of [numpy-style doc strings](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt). They are nice to read as normal text and are converted by sphinx to normal ReST through [napoleon](http://sphinxcontrib-napoleon.readthedocs.org/en/latest/index.html). All previous doc-strings are using pure ReST, and we use the same approach as for the coding-style here. When you touch a function please update its docstring to follow the numpy-style. All new functions should have numpy-style docstrings.

```python
def func(arg1, arg2):
    """Summary line.

    Extended description of function.

    Parameters
    ----------
    arg1 : int
        Description of arg1
    arg2 : str
        Description of arg2

    Returns
    -------
    bool
        Description of return value

    """
    return True
```

## Code organization
### Package and test suite
The [[source code|Source#source-code-organisation]] is distributed over the `package` directory (the library and documentation) and the `testsuite` (test code and data files). Commits can contain code in both directories, e.g. a new feature and a test case can be committed together.

### Compiled code
MDAnalysis contains compiled code (cython, C, C++) in addition to pure Python. With [#444](/MDAnalysis/mdanalysis/issues/444) we agreed on the following layout:

*Place all source files for compiled shared object files into the same folder as the final shared object file.*

`*.pyx` files and cython-generated `*.c` would be in the same directory as the `*.so`. External dependent C/C++/Fortran libraries are in dedicated `src` and `include` folders. See the following tree as an example.

```
MDAnalysis 
      |--lib
      |   |-- _distances.so
      |   |-- distances.pyx
      |   |-- distances.c
      |-- coordinates
          |-- _dcdmodule.so
          |-- src
              |-- dcd.c
          |-- include
              |-- dcd.h
```


This is standard. See [numpy](https://github.com/numpy/numpy/tree/master/numpy/linalg/lapack_lite), [scipy](https://github.com/scipy/scipy/tree/master/scipy/spatial), [scikit-learn](https://github.com/scikit-learn/scikit-learn/tree/master/sklearn/svm/src), [mdtraj](https://github.com/mdtraj/mdtraj/tree/master/mdtraj/formats/xtc), for other examples.

## Version control
We use [[git]] for version control. The [[DevelopmentWorkflow]] uses *pull requests* against the [[DevelopmentBranch]] (named *develop*) followed by automated testing and code review by project members. Successful PRs are merged into *develop*.

### Commit messages
Follow [[git commits|Style-Guide:-git-commits]].

## Tests

For now, see [[Writing Tests|UnitTests#writing-test-cases]] (but that needs to be cleaned up). In short:

* changes and additions in the **core** (everything except `MDAnalysis.analysis` and `MDAnalysis.visualization`): unit tests are **mandatory**
* changes and additions to 
   * `MDAnalysis.analysis`
   * `MDAnalysis.visualization`

  Tests are **highly encouraged** (and anyone reviewing commits can ask for at least minimal tests)
