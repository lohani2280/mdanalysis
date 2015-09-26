MDAnalysis defines the [NoDataError](http://devdocs.mdanalysis.org/documentation_pages/core/AtomGroup.html?highlight=nodataerror#MDAnalysis.core.AtomGroup.NoDataError) exception for cases when there's insufficient data to operate on. This can happen, for instance, when a selection did not match anything.

Until release 0.7.5.1 a NoDataError was raised for  selections that did not match anything or an AtomGroup built from an empty list ([Issue 12](http://issues.mdanalysis.org/12)) but this has been changed in [release 0.7.6](ReleaseNotes076) to return simply an empty AtomGroup.
