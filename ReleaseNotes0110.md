Release 0.11.0 is a major update for MDAnalysis. It contains a number of *backwards-incompatible API changes*, i.e. **existing code is likely to break**. Additionally, some old functionality has been deprecated and will be removed for the 1.0.0 release.

The [MDAnalysis 0.11 unifying release user guide](MDAnalysis-0.11-unifying-release-user-guide) should help you with porting your code. We have also included a tool `ten2eleven` as part of the package that automates some of the API changes in your existing code, as described in [Migrating MDAnalysis code with ten2eleven](Migrating-MDAnalysis-code-with-ten2eleven.py).

## CHANGELOG

  This release brings large changes to many parts of the API and might
  break many existing scripts. For a full guide on the API changes,
  please see:

  https://github.com/MDAnalysis/mdanalysis/wiki/MDAnalysis-0.11-unifying-release-user-guide

  Migrating old scripts has been made easier with the introduction of
  the ten2eleven tool which is part of the package.
  Details on how to use this are available at:

  https://github.com/MDAnalysis/mdanalysis/wiki/Migrating-MDAnalysis-code-with-ten2eleven.py

###  API Changes

  * Changed AtomGroup counting methods to properties with different
    names: numberOfAtoms() to n_atoms, numberOfResidues() to
    n_residues, numberOfSegments() --> n_segments (Issue #376)
  * Changed trajectory reader numframes to n_frames (Issue #376)
  * Changed Timestep.numatoms to n_atoms (Issue #376)
  * Deprecated the use of the 'fullgroup' selection keyword (Issue #268)
  * Changed atom.number attribute to atom.index (Issue #372)
  * Changed many AtomGroup methods to properties.  These are: indices,
    masses, charges, names, types, radii, resids, resnames, resnums,
    segids (Issue #372)
  * Timestep can now only init using an integer argument (which
    represents the number of atoms) (Issue #250)
  * Added from_timestep and from_coordinates construction methods
    to base.Timestep (Issue #250)
  * Removed KDTree and CoordinateNeighbor from MDAnalaysis. If you
    want to search in cartesian coordinates directly for nighboring
    points use the BioPython KDTree or scikit-learn Neighbors module.
    The AtomNeighborSearch class has been ported to use the BioPython
    KDTree and is now located in MDAnalaysis.lib.NeighborSearch.
    MDAnalaysis.KDTree still exists in this version so load the
    NeighborSearch module but is deprecated and will be removed in
    1.0. (Issue #383)
  * Moved MDAnalysis.core.transformations to
    MDAnalysis.lib.transformations (Issue #287)
  * Moved MDAnalysis.core.util to MDAnalysis.lib.util (Issue #287)
  * Moved MDAnalysis.core.log to MDAnalysis.lib.log (Issue #287)
  * Moved MDAnalysis.core.units to MDAnalysis.units (Issue #287)
  * Moved MDAnalysis.core.distances to MDAnalysis.lib.distances
    (Issue #287)
  * Moved MDAnalysis.core.parallel to MDAnalysis.lib.parallel
    (Issue #287)
  * Moved norm, normal, angle, stp and dihedral from lib.util to
    lib.mdamath (Issue #287)
  * AtomGroup.bond .angle .dihedral and .improper now return the
    corresponding TopologyObject rather than performing the calculation
    (Issue #373)
  * All TopologyObjects now have a "value" method to evaluate them
    (Issue #373)
  * TopologyGroup now has a "values" methods to evaluate all contained
    bonds (Issue #373)
  * MDAnalysis.lib.distances.calc_torsions renamed to calc_dihedrals
    (Issue #373)
  * TopologyGroup.selectBonds renamed to select_bonds (Issue #389)
  * deprecated camelCase AtomGroup methods in favour of underscore_style
    (Issue #389)
  * deprecate lib.distances.applyPBC in favour of apply_PBC (Issue #389)
  * AtomGroup.res[names,ids,nums] and AtomGroup.segids now give arrays
    of equal size to AtomGroup (Issue #385)
  * ResidueGroup.segids now gives arrays of equal size to ResidueGroup
    (Issue #385)
  * AtomGroup setters `set_<property>` now plural for consistency with
    property names (Issue #385)
  * DCDReader no longer supports the "skip" keyword.  Use slicing
    on reader iteration to achieve same affect. (Issue #350)
  * All Readers have a default "dt" of 1.0 ps.  This applies also to
    single frame readers, these would previously raise an error on
    accessing dt. (Issue #350)
  * NCDF Reader no longer has a default skip_timestep of 1 (Issue #350)

###  Enhancements

  * Added the 'global' selection keyword (Issue #268)
  * Added Jmol selection writer (Issue #356)
  * Added reading of Hoomd XML files (Issue #333)
    These can only act as topology information for now
  * Tests can now detect memleaks on a per-test basis (Issue #323)
  * AtomGroups can now be pickled/unpickled (Issue #293)
  * Universes can have a `__del__` method (not actually added) without
    leaking (Issue #297)
  * Added reading of DL_Poly format CONFIG and HISTORY files, these can
    both act as both Topology and Coordinate information (Issue #298)
  * Timestep objects now have `__eq__` method (Issue #294)
  * coordinates.base.Timestep now can handle velocities and forces
    (Issue #294)
  * Waterdynamics analysis module added, including five analysis
    classes: Hydrogen Bond Lifetimes, Water Orientational Relaxation,
    Angular Distribution, Mean Square Displacement and Survival
    Probability. Documentation and test are included. (Issue #300)
  * RMSF class added to rms analysis module
  * ProgressMeter now outputs every *interval* number of ``update``
    calls (Issue #313)
  * Created lib.mdamath for common mathematical functions. (Issue #287)
  * All Timesteps have the has_positions has_velocities and has_forces
    boolean flags (Issue #213)
  * Timesteps can now allocate velocities and forces if they weren't
    originally created with these through the use of the has_ flags.
    (Issue #213)
  * Timesteps now store 'dt' and 'time_offset' if passed to them by
    Reader to calculate time attribute (Issues #306 and #307)
  * MDAnalysis.selection: can also be written to a NamedStream
  * Added function lib.mdamath.make_whole to "unbreak" molecules
    over periodic boundaries. (Issue #355)
  * Added triclinic_dimensions to Timestep, returns representation of
    unit cell as triclinic vectors (Issue #276)
  * Added setter to bfactors property (Issue #372)
  * Added AtomGroup altLocs and serials properties with setters.
    (Issue #372)
  * MDAnalysis.core.AtomGroup.Merge now copies across bonding
    information (Issue #249)

###  Changes

  * numpy >= 1.5 required
  * A ProtoReader class intermediate between IObase and Reader was added
    so specific Readers can be subclassed without `__del__` (the
    ChainReader and SingleFrameReader), thus preventing memleaks
    (Issue #312).
  * Atoms (and all container classes thereof) are now bound to Universes
    only via weakrefs. If Universes are not explicitly kept in scope
    Atoms will become orphaned. (Issue #297)
  * Removed FormatError, now replaced by ValueError (Issue #294)
  * base.Reader now defines `__iter__` and `__iter__` removed from many
    Readers (now use base.Reader implementation) (Issue #350)
  * Timestep._x _y and _z are now read only (Issue #213)
  * moved MDAnalysis.selections.base.get_writer() to
    MDAnalysis.selections.get_writer() to break a circular import. This
    should not affect any code because
    MDAnalysis.selections.get_writer() already existed.
  * distances.contact_matrix now treats the particle distance with
    itself as a contact for sparse matrices and numpy arrays. The
    progress reporting for sparse calculations has been removed.
    (Issue #375)
  * TopologyObjects and TopologyGroup moved to core.topologyobjects
    module (Issue #373)
  * Consolidated coordinates.guess_format and topology.guess_format to
    lib.util.guess_format  (Issue #336)

###  Fixes

  * All Writers now refer to time between written Timesteps as "dt",
    was previously "delta" in some Writers. (Issue #206)
  * Topology files can now be compressed (Issue #336)
  * Fixed PDB Parser and Reader requiring occupancy field (Issue #396)
  * Amber TRJ and NCDF Reader & Writers now use 'dt' instead of 'delta'
    to refer to time elapsed between timesteps. (Issue #350 and #206)
  * Fixed TPRParser considering LJ 1..4 exclusions as bonds (Issue #351)
  * ChainReaders no longer cause memory leak (Issue #312)
  * analysis.hbonds.HydrogenBondAnalysis performs a sanity check for
    static selections (Issue #296)
  * Fixed TRZWriter failing when passed a non TRZTimestep (Issue #302)
  * relative imports are now banned in unit testing modules
    (Issue #189)
  * Fixed bug and added DivisionByZero exception in
    analysis/waterdynamics.py in SurvivalProbability. (Issue #327)
  * Fixed parsing of PDB header data for PrimitivePDBReader (Issue #332)
  * Fixed contact_matrix handles periodic boundary conditions correctly
    for sparse matrices. (Issue #375)
  * Fixed analysis.hole not using CPOINT (Issue #384)
  * Fixed XTC/TRR .dt rewinding the trajectory (Issue #407)
  * Fixed TopologyGroup.from_indices not guessing topology object class
    (Issue #409)
  * Fixed TopologyGroup.`__contains__` failing if different instance of
    same bond was queried. (Issue #409)

## Authors

tyler.je.reddy, richardjgowers, alejob, orbeckst, dotsdl,
        manuel.nuno.melo, cyanezstange, khuston, ivirshup, kain88-de,
        gormanstock