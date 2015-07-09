# Contents #


# New Guidelines Proposal (March 2014) #
In general, we are proposing (still open to some debate) to adjust the MDAnalysis documentation to match the formatting standards adopted by numpy / scipy as described in the [documentation standard](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt). Much of the scientific Python community adopts this standard (numpy / scipy / scikit) and matplotlib appears to be moving in this direction as well with [MEP10](https://github.com/matplotlib/matplotlib/wiki/Mep10). Arguably, new MDAnalysis users from the scientific Python community would appreciate documentation that is consistent with other major scientific libraries, and this format seems to work quite well.

There's a lot to be said for clear and concise documentation--it can save the users and developers time and attract new users (and eventually developers). It's quite nice to be able to google a function or class and obtain a single clear page with its documentation in a familiar format.

The format generally involves a single class / function per documentation page, with example usage strongly encouraged (see, for example: [scipy.spatial.ConvexHull](http://docs.scipy.org/doc/scipy-dev/reference/generated/scipy.spatial.ConvexHull.html)). The order and formatting of sections is perhaps best understood by looking at the [example](https://github.com/numpy/numpy/blob/master/doc/example.py) provided in the numpy documentation standard. In general, we want to avoid large / bloated documentation pages with several different functions / classes. If appropriate, common classes and functions can be associated in a summary table as shown for [scipy.spatial.distance](http://docs.scipy.org/doc/scipy/reference/spatial.distance.html).

There are certainly some points to dicuss: 1) will we adopt a strict [doctest](http://docs.python.org/library/doctest.html) format for examples in docstrings (not meant to replace unit tests)? 2) Will we adopt the [numpydoc](https://pypi.python.org/pypi/numpydoc) Sphinx extension? This would be needed for a more precise match to the numpy documentation standard. 3) other issues?

While I've tried to be reasonably consistent in my recently-added modules (i.e., [3D streamplot](http://mdanalysis.googlecode.com/git-history/develop/package/doc/html/documentation_pages/visualization/streamlines_3D.html)), I'm also guilty of cutting corners in a few of the formatting areas. The idea is to make changes that provide increased consistency with the numpy documentation standard if you stumble across a module that you think could use enhanced documentation.


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

As of Jan 2011, the **[Online Documentation](http://mdanalysis.googlecode.com/git/doc/html/index.html)** is available, too.

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
  * The [Online Documentation](http://mdanalysis.googlecode.com/git/doc/html/index.html) is generated from [reStructuredText](http://docutils.sourceforge.net/rst.html) in [doc/sphinx/source/documentation\_pages](http://code.google.com/p/mdanalysis/source/browse/#git%2Fpackage%2Fdoc%2Fsphinx%2Fsource%2Fdocumentation_pages).
  * The googlecode pages also hold documentation in the [MDAnalysis Wiki](http://code.google.com/p/mdanalysis/w/list) (but as far as the user docs go, they are being phased out in favour of the [Online Docs](http://mdanalysis.googlecode.com/git/doc/html/index.html)).


# Online Docs format: sphinx reStructuredText #

## Sphinx basics ##
We are using [reStructuredText](http://docutils.sourceforge.net/rst.html) ("rst" or "reST") in the Python code and in the [Online Documentation](http://mdanalysis.googlecode.com/git/doc/html/index.html). The reST is processed with [sphinx](http://sphinx.pocoo.org/) to generate HTML or PDF. Thus, the docs should use reST with the sphinx autodoc extensions:
  * sphinx [reST Primer](http://sphinx.pocoo.org/rest.html)
  * sphinx [autodoc extension](http://sphinx.pocoo.org/ext/autodoc.html)

Note that each page of the  [Online Documentation](http://mdanalysis.googlecode.com/git/doc/html/index.html) has a link to the _Source_ of the page. You can look at it in order to find out how a particular page has been written in reST and copy the approach for your own documentation.

## Generating the docs ##
The documentation n HTML format lives in `package/doc/html`. HTML is the primary format in which we provide the docs. Documentation for *releases* is hosted on PyPi and is accessible at https://pythonhosted.org/MDAnalysis (and also http://docs.mdanalysis.org). The documentation for the *development branch* is hosted as [GitHub pages](https://pages.github.com/) in the *gh-pages* branch and is available at http://devdocs.mdanalysis.org.

Documentation is revision controlled and part of the git repository. However, follow the following rules:
* *Generated html and other files* are checked into the git repository *as a separate commit* from any doc changes that you made either in the reST or py files. Please do not mix autogenerated files with manual edits. 
* Newly generated files _must_ be added (`git add`) or they don't show up in the online docs! Files that were removed (e.g. if files were renamed) also need to be removed/renamed in the docs. If you know that any files were renamed or removed (or if you just want to make sure) then start with removing the whole html directory:
```
cd package/doc && rm -r html
```

Building the docs is identical for release and development but distribution differs:

### Build the docs with sphinx ###

In order to **generate the documentation**, one has to have [sphinx](http://sphinx-doc.org/) installed. The *current* version of MDAnalysis must be installed and MDAnalysis must be importable. A 'develop' installation (`python setup.py develop`) is generally useful for quick turn-around.

#### Running sphinx #####
Then one generates the html docs with
```
cd package/doc/sphinx
make html
```
This generates and updates the files in `doc/html`.

#### Checking ####
You should point your browser to the file `../html/index.html` and look through the docs, in particular any pages that you tinkered with.

As of this writing (July 2015), there are still lots of reST warnings and it probably takes a lifetime to fix them (and even that isn't straightforward because it's not always clear where the problem is in the file). Check if there are some serious issues, e.g. "document referenced but no heading" – that indicates some import issues or interference from Gremlins or Faeries or Pixies or all of the above. Visually check anything that looks suspicious in the error output.

It is typical to go through multiple cycles of fix, rebuild the docs with `make html`, check and fix again.

If *any fixes in the restructured text* are needed, *put them in their own commit* (and do not include any generated files under `docs/html`) — `git add FILE` and `git commit --amend` is your friend when piling more and more small reST fixes onto a single "fixed reST" commit.

(For testing, you can also use `python setup.py build_sphinx` which creates docs under `build/sphinx/html` but these files are never committed to the repository so remember to use the procedure with `make` to generate online docs.)

#### Committing changes ####
If all is good, 

1. Run `git status`.
2. `git add` any newly added files under `html` and 
3. properly `git rm` any files that were not rebuilt and show up as deleted. If in doubt `rm -r ../html && make html` again to make sure that only newly generated files are present and old ones show up as deleted under `git status`.
4. Commit with a message such as "regenerated docs with sphinx (new frobazz module now included)". 
5. Push the changes.

### Distribute development docs ###
Now update the [developer docs](http://devdocs.mdanalysis.org) in the *gh-pages* branch:
```
git fetch --all   # just in case
git checkout -b gh-pages
git merge develop
git push origin gh-pages
```
The [developer docs](http://devdocs.mdanalysis.org) should update.

### Distribute docs for releases ###
Currently we manually upload a zip file of the docs to PyPi as described under [PreparingReleases: Update docs on PyPI](PreparingReleases#Update_docs_on_PyPI).

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

# Serving documentation #
## Released versions ##
When a new version is released, the [[docs on pythonhosted are updated|PreparingReleases#update-docs-on-pypi]] as part of the release procedure. They can be accessed from http://docs.mdanalysis.org

## Developer version ##
The http://devdocs.mdanalysis.org are generated from the develop branch but because they are served through gh-pages, a manual procedure needs to be carried out to make it work (perhaps run a webhook or something... – and no, readthedocs did not work for us, see [#183](https://github.com/MDAnalysis/mdanalysis/issues/183) and [#235](https://github.com/MDAnalysis/mdanalysis/issues/235)):

```
git checkout gh-pages
git merge develop
git push origin gh-pages
```



