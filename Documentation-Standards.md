# Proposal for new Documentation Guidelines (2014) #

* author: @tylerjreddy
* link to issue: ?

As set out in the [[Style Guide]], we are moving the MDAnalysis documentation to match the formatting standards adopted by numpy / scipy as described in the [documentation standard](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt). Much of the scientific Python community adopts this standard (numpy / scipy / scikit) and matplotlib appears to be moving in this direction as well with [MEP10](https://github.com/matplotlib/matplotlib/wiki/Mep10). Arguably, new MDAnalysis users from the scientific Python community would appreciate documentation that is consistent with other major scientific libraries, and this format seems to work quite well.

There's a lot to be said for clear and concise documentation--it can save the users and developers time and attract new users (and eventually developers). It's quite nice to be able to google a function or class and obtain a single clear page with its documentation in a familiar format.

The format generally involves a single class / function per documentation page, with example usage strongly encouraged (see, for example: [scipy.spatial.ConvexHull](http://docs.scipy.org/doc/scipy-dev/reference/generated/scipy.spatial.ConvexHull.html)). The order and formatting of sections is perhaps best understood by looking at the [example](https://github.com/numpy/numpy/blob/master/doc/example.py) provided in the numpy documentation standard. In general, we want to avoid large / bloated documentation pages with several different functions / classes. If appropriate, common classes and functions can be associated in a summary table as shown for [scipy.spatial.distance](http://docs.scipy.org/doc/scipy/reference/spatial.distance.html).

There are certainly some points to dicuss: 1) will we adopt a strict [doctest](http://docs.python.org/library/doctest.html) format for examples in docstrings (not meant to replace unit tests)? 2) Will we adopt the [numpydoc](https://pypi.python.org/pypi/numpydoc) Sphinx extension? This would be needed for a more precise match to the numpy documentation standard. 3) other issues?

While I've tried to be reasonably consistent in my recently-added modules (i.e., [3D streamplot](http://devdocs.mdanalysis.org/documentation_pages/visualization/streamlines_3D.html)), I'm also guilty of cutting corners in a few of the formatting areas. The idea is to make changes that provide increased consistency with the numpy documentation standard if you stumble across a module that you think could use enhanced documentation.
