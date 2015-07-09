This is a minor new release of MDAnalysis, which contains a couple of [enhancements](#Enhancements).

## CHANGELOG ##

  * 0.9.2 
  * released 04/20/15  

###  Enhancements ###

  * Can now set velocity from Atom object. (Issue #221)
  * Atom object now has force attribute for getting/setting forces (requires a trajectory with forces) (Issue #221)
  * Added wrap method. Wrap allows the centers of groups of atoms to be  moved within the primary unit cell without breaking up the group. (Issue #190)

###  Changes ###

  * The MDAnalysis project moved from Google Code to GitHub: the new    website is http://www.mdanalysis.org and the new source code    repository is at https://github.com/MDAnalysis/mdanalysis
  * "applications" were removed from the mdanalysis source code    repository and now exist as independent repositories under    https://github.com/MDAnalysis/
  * Using a non-existent atom name as an instant selector now raises    AttributeError instead of SelectionError (Issue #220)

###  Fixes ###

  * trajectory objects are now properly closed in unit tests (Issue #256)

## Authors ##
tyler.je.reddy, richardjgowers, orbeckst
