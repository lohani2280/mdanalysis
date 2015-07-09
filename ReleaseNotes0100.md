This is a new minor release of MDAnalysis, which contains a large number of [enhancements](#Enhancements), [changes to the API](#Changes) and a few [bug fixes](#Fixes).

## CHANGELOG ##

  * 0.10.0
  * released 06/01/15 

###  Enhancements ###

  * Improved performance of PDB Reading.  Up to 3x faster. (Issue #212)
  * Added the 'same ... as' selection keyword (Issue #217)
  * Added guess_bonds keyword argument to Universe creation.  This will attempt to    guess all topology information on Universe creation. (Issue #245)
  * Added guess_bonds method to AtomGroup. (Issue #245)
  * All TopologyObjects (Bond, Angle etc) now have is_guessed attribute
  * TopologyGroup now has alternate constructor method, .from_indices()
  * Added TopologyObject.indices property
  * Amber netCDF4 Reader will now read Forces (Issue #257)
  * Amber netCDF4 Writer will now write Velocities and Forces
  * Added Amber coordinate/restart file reader (Issue #262)
  * Structural superpositions (MDAnalysis.analysis.align) can work    with partial matches of atoms.
  * new path similarity analysis module MDAnalysis.analysis.psa
  * AtomGroup and TopologyGroup can now be indexed by numpy boolean arrays    works identically to numpy masks. (Issue #282)

###  Changes ###

  * TopologyGroup can now have zero length, and will evaluate to False    when empty.
  * Renamed TopologyGroup.dump_contents to "to_indices"
  * Deprecated 'bonds' keyword from Universe and replaced with 'guess_bonds'
  * PrimitivePDBReader now requires the numatoms keyword
  * Structural superpositions (MDAnalysis.analysis.align) use partial    matches of atoms by default (use strict=True for old behavior)
  * Function rmsd() was removed from MDAnalysis.analysis.align name    space and should be accessed as MDAnalysis.analysis.rms.rmsd()

###  Fixes ###

  * bynum selections now work from AtomGroup instances (Issue #275)
  * Cylinder selections now work from AtomGroup instances and honor
    PBC (Issue #274)
  * NetCDFWriter previously always wrote velocities/forces if found
    in timestep, rather than following how the Writer was created.

## Authors ##
richardjgowers, caio s. souza, manuel.nuno.melo, orbeckst, sseyler
