See also
* proposed [[Documentation Standards]]
* [[Style Guide]]

# Overview #
MDAnalysis has a lot of documentation in the Python doc strings. In the Python interpreter one can simply say
```
  import MDAnalysis
  help(MDAnalysis)
  help(MDAnalysis.Universe)
```
In `ipython` one can use the question mark operator
```
  MDAnalysis.Universe?
```
Interactive use is not always convenient, though, and hence [Issue 26](http://issues.mdanalysis.org/26) called for proper documentation.

As of Jan 2011, the **[Online Documentation](http://docs.mdanalysis.org/index.html)** is available, too.

This page should help developers write good documentation that can be easily processed into online docs.

## Guidelines ##

When writing Python code, you should always add a doc string to each public (user-visible)
  * module
  * class
  * method
  * function
(We consider something public if the name does not start with an underscore.)
The module doc string can be a short sentence describing what the module does or a long document including examples and references. The other doc strings should generally allow a user to use the code element, so they need to say what the code does, what input it requires, and what it returns. Additionally, one should document the exceptions raised and known bugs or quirks.

## Location of documentation ##
  * Code is documented through Python doc strings (see [PEP 257](http://www.python.org/dev/peps/pep-0257/) for good conventions); thus each Python file should contain some documentation.
  * The [Online Documentation](http://docs.mdanalysis.org/index.html) is generated from [reStructuredText](http://docutils.sourceforge.net/rst.html) in [doc/sphinx/source/documentation\_pages](https://github.com/MDAnalysis/mdanalysis/tree/develop/package/doc/sphinx/source/documentation_pages).
  * As outlined in the [[Style Guide]], we use NumPy format for the reST strings
  * The wiki should not duplicate documentation; if it's important for a user of the library then it has to be in the docs.

# Online Docs format: sphinx reStructuredText #

## Sphinx basics ##
We are using [reStructuredText](http://docutils.sourceforge.net/rst.html) ("rst" or "reST") in the Python code and in the [Online Documentation](http://docs.mdanalysis.org/index.html). The reST is processed with [sphinx](http://sphinx.pocoo.org/) to generate HTML or PDF. Thus, the docs should use reST with the sphinx autodoc extensions:
  * sphinx [reST Primer](http://sphinx.pocoo.org/rest.html)
  * sphinx [autodoc extension](http://sphinx.pocoo.org/ext/autodoc.html)

Note that each page of the  [Online Documentation](http://docs.mdanalysis.org/index.html) has a link to the _Source_ of the page. You can look at it in order to find out how a particular page has been written in reST and copy the approach for your own documentation.

## Generating the docs ##
The documentation n HTML format lives in `package/doc/html`. HTML is the primary format in which we provide the docs. Documentation for *releases* is hosted on PyPi and is accessible at https://pythonhosted.org/MDAnalysis (and also http://docs.mdanalysis.org). The documentation for the *development branch* is hosted as [GitHub pages](https://pages.github.com/) in the *gh-pages* branch and is available at http://devdocs.mdanalysis.org.

* release docs are built by a maintainer and uploaded to PyPi
* development docs are built automatically by Travis CI (see issue [#386](/MDAnalysis/mdanalysis/issues/386) for details)

### Build the docs with sphinx ###

In order to **generate the documentation**, one has to have [sphinx](http://sphinx-doc.org/) installed. The *current* version of MDAnalysis must be installed and MDAnalysis must be importable. A 'develop' installation (`python setup.py develop`) is generally useful for quick turn-around.

#### Running sphinx #####
Then one generates the html docs with
```
cd package/doc/sphinx
make html
```
This generates and updates the files in `doc/html`.

Alternatively
```
python setup.py build_sphinx
```
might just work.

#### Checking ####
You should point your browser to the file `../html/index.html` and look through the docs, in particular any pages that you tinkered with.

As of this writing (July 2015), there are still lots of reST warnings and it probably takes a lifetime to fix them (and even that isn't straightforward because it's not always clear where the problem is in the file). Check if there are some serious issues, e.g. "document referenced but no heading" – that indicates some import issues or interference from Gremlins or Faeries or Pixies or all of the above. Visually check anything that looks suspicious in the error output.

It is typical to go through multiple cycles of fix, rebuild the docs with `make html`, check and fix again.

If *any fixes in the restructured text* are needed, *put them in their own commit* (and do not include any generated files under `docs/html`) — `git add FILE` and `git commit --amend` is your friend when piling more and more small reST fixes onto a single "fixed reST" commit.

### Distribute docs for releases ###
Currently we manually upload a zip file of the docs to PyPi as described under [PreparingReleases: Update docs on PyPI](PreparingReleases#Update_docs_on_PyPI).

Or the following might also work
```
python setup.py build_sphinx upload_docs
```

## Mathematics ##
We are using MathJax with sphinx so you can write LaTeX code in math tags, e.g.
```
#<SPACE if there is text above equation>
.. math::
   e^{i\pi} = -1
```
or inline
```
We make use of the identity :math:`e^{i\pi} = -1` to show...
```

Note that you should _always_ make doc strings with math code **raw** python strings **by prefixing them with the letter "r"**:
```
def rotate(self, R):
    r"""Apply a rotation matrix *R* to the selection's coordinates.

    :math:`\mathsf{R}` is a 3x3 orthogonal matrix that transforms a vector
    :math:`\mathbf{x} \rightarrow \mathbf{x}'`:

    .. math::

       \mathbf{x}' = \mathsf{R}\mathbf{x}
    """
```
or else you will get problems with backslashes in unexpected places (see [Stackoverflow: Mathjax expression in sphinx python not rendering correctly](http://stackoverflow.com/questions/16468397/mathjax-expression-in-sphinx-python-not-rendering-correclty")).

