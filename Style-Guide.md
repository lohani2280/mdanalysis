# Code style guidelines

(For background and history see [Issue #404](https://github.com/MDAnalysis/mdanalysis/issues/404), which contains the original discussion on the Style Guide.)

## Coding Style

MDAnalysis is an old project and hasn't used a consistent coding style. Since version
0.11.0 we are trying to update all the code to conform with [PEP8](http://legacy.python.org/dev/peps/pep-0008/). Our strategy is to update the style every time we touch an old function and thus switch
to [PEP8](http://legacy.python.org/dev/peps/pep-0008/) continuously.

### Differences to PEP8

 * 80 columns line length instead of 79

### Code Linter

We recommend that you either use a Python IDE ([PyCharm](https://www.jetbrains.com/pycharm/) and others). You can also use external tools like [flake8](http://flake8.readthedocs.org/en/latest/). For integration of external tools with emacs and vim check out [elpy (emacs)](https://github.com/jorgenschaefer/elpy) and [python-mode (vim)](https://github.com/klen/python-mode)

### Code Formater

To apply the code formatting in an automated way you can also use code formatters. As external tools there are [autopep8](https://github.com/hhatto/autopep8) and [yapf](https://github.com/google/yapf). Most IDE's either have their own code formatter or will work with one of the above through plugins.

## Writing Docstrings

Since 0.11.0 we adopted to use the numpy-style doc strings. They are nice to read as normal text and are converted by sphinx to normal ReST through [napoleon](http://sphinxcontrib-napoleon.readthedocs.org/en/latest/index.html). All previous doc-strings are using pure ReST, we use the same approach as for the coding-style here. When you touch a function please update it's docstring to follow the numpy-style. All new function
should have numpy-style docstrings.

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
MDAnalysis contains compiled code (cython, C, C++) in addition to pure Python. With [[#444|/MDAnalysis/mdanalysis/issues/444]] we agreed on the following layout:

*Place all source files for compiled shared object files into the same folder as the final shared object file.*

`*.pyx` files and cython-generated `*.c` would be in the same directory as the `*.so`. While external dependend C/C++/Fortran libraries are in dedicate `src` and `include` folders. See the following tree as an example.

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


This is standard see [numpy](https://github.com/numpy/numpy/tree/master/numpy/linalg/lapack_lite), [scipy](https://github.com/scipy/scipy/tree/master/scipy/spatial), [scikit-learn](https://github.com/scikit-learn/scikit-learn/tree/master/sklearn/svm/src). [mdtraj](https://github.com/mdtraj/mdtraj/tree/master/mdtraj/formats/xtc)

## Commit messages

Follow [[git commits|Style-Guide:-git-commits]].

## Writing Tests

For now, see [[Writing Tests|UnitTests]] (but that needs to be cleaned up). In short:

* changes and additions in the **core** (everything except `MDAnalysis.analysis` and `MDAnalysis.visualization`): unit tests are **mandatory**
* changes and additions to 
   * `MDAnalysis.analysis`
   * `MDAnalysis.visualization`
  Tests are **highly encouraged** (and anyone reviewing commits can ask for at least minimal tests)
