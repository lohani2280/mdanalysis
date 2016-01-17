## CHANGELOG

API Changes

  * ChainReader `delta` keyword deprecated in favor of `dt`. (Issue #522)
  * XYZWriter can now be used as a container format for different protein models
    as well as a normal trajectory. If `n_atoms` is None (default) MDAnalysis
    assumes that it is used as a container and won't give a warning if the
    number of atoms differs between frames.
  * GROWriter.fmt strings updated to use format style (Issue #494)
  * removed MDAnalysis.lib.parallel.distances; use the new backend="OpenMP"
    keyword for the functions in MDAnalysis.lib.distances (Issue #530)

Enhancement

  * ChainReader now reports times properly summed over sub-readers (Issue #522)
  * GRO file reading approximately 50% faster for large files (Issue #212)
  * GRO file writing now will write velocities where possible (Issue #494)
  * Added bonded selection (Issue #362)
  * Spherical layer and spherical zone selections now much faster (Issue #362)
  * new keyword "backend" for accelerated functions in MDAnalysis.lib.distances
    to select "serial" or "OpenMP"-enabled versions of the code; the default
    is "serial" so old code will behave as before (see Issue #530)
  * Lammps data file parsing improved greatly.  Should now support all files,
    and triclinic geometry. (Issue #139)
  * Added analysis.polymer, currently with PersistenceLength tool (Issue #460)
  * Added analysis.rdf, with InterRDF tool. (Issue #460)
  * Made Reader.check_slice_indices a public method (Issue #604)
  * analysis.helanal.helanal_main() now returns results as dict
  * Added keyword to update selection every frame in density calculation (Issue #584)
  * New keywords start, stop, step for density.density_from_Universe()
    to slice a trajectory.
  * MOL2Reader now reads molecule and substructure into ts.data
  * All subclasses of ProtoReader, Writer and TopologyReader are automatically
    added to the MDAnalysis directory of I/O (Issue #431)

Changes

  * built html doc files are no longer version controlled (Issue #491)
  * The lib._distances and lib_distances_openmp libraries now have a
    OPENMP_ENABLED boolean flag which indicates if openmp was used in
    compilation.  (Issue #530)
  * analysis.helanal.helanal_trajectory() and helanal_main() now use a
    logger at level INFO to output all their computed values instead
    of printing to stdout
  * default offset for ProgressMeter was changed from 0 to 1 (to match
    the change from 1- to 0-based ts.frame counting)
  * removed superfluous analysis.density.density_from_trajectory();
    use density_from_Universe(TOPOL, TRAJ) instead.
  * MOL2Writer.write now only writes a single frame (Issue #521)

Fixes

  * Fixed select_atoms requiring a trajectory be loaded (Issue #270)
  * AtomGroup timesteps no longer cached (Issue #606)
  * GROWriter now truncates atom numbers over 99999 (Issue #550)
  * AMBER netcdf writer now correctly uses float32 precision (Issue #518)
  * Fixed a numpy incompatibility in `analysis.leaflet.LeafletFinder`.
    (Issue #533)
  * Cleaned up `MDAnalysis.Writer` docs regarding `dt` usage. (Issue #522)
  * Fixed setup-time dependency on numpy that broke pip installs. (Issue #479)
  * Fixed unpickling errors due to lingering dead universes. (Issue #487)
  * Fixed analysis.density modules requiring the defunct `skip` attribute
    on trajectories. (Issue #489)
  * ten2eleven camelcase fixer now deals with centerOfMass (Issue #470)
  * ten2eleven will now convert numatoms to n_atoms argument
    for writer() functions (Issue #470)
  * Fixed non-compliant Amber NCDFWriter (Issue #488)
  * Fixed many Timestep methods failing when positions weren't present
    (Issue #512)
  * Fixed PointSelection using KDTree (Issue #362)
  * Fixed GROParser getting tripped up by some file (Issue #548)
  * Fixed writing dx files from analysis.density.density_from_Universe()
    (Issue #544 and #410)
  * Fixed base.Reader._check_slice_indices not liking numpy ints
    (Issue #604)
  * Fixed broken analysis.helanal.helanal_trajectory() and
    helanal_main()
  * Fixed lib.util.greedy_splitext() (now returns full path)
  * Fixed MOL2Reader not reading molecule and substructure on init
    (Issue #521)
  * Fixed MOL2Writer rereading frames when writing them (Issue #521)
  * Fixed PDBWriter not writing occupancies from atoms (Issue #620)

## AUTHORS

  tyler.je.reddy, kain88-de, richardjgowers, manuel.nuno.melo,
  orbeckst, Balasubra