## Highlights

 * [Sean Seyler's](https://github.com/sseyler) Path Similarity analysis module got many more features.  Check out the [tutorial](https://github.com/Becksteinlab/PSAnalysisTutorial)
 * Parsing of TPR files was vastly improved, including supporting v5.0+ of [Gromacs](https://github.com/gromacs/gromacs) thanks to [Jonathan Barnoud](https://github.com/jbarnoud)
 * The installation process has been improved, any various minor annoyances removed, credit to [Max Linke](https://github.com/kain88-de)

## CHANGELOG


API Changes

  * PrimitivePDBReader now imports occupancies into the `TimeStep` object.
    (Issue #396)
  * Atoms without a Universe now return NoDataError instead of
    AttributeError
  * AtomGroups of zero length or containing Atoms with no Universe raise
    a NoDataError when trying to access Universe
  * Atoms now keep a strong reference to Universe, meaning they
    do not become orphaned when the Universe goes out of scope (Issue #297)

Enhancements

  * `Atom` and `AtomGroup` now expose occupancy value as `occupancy` and
    `occupancies` properties (Issue #396)
  * XYZReader now supports frame indexing (Issue #428)
  * Reader objects can now be sliced using lists and arrays of indices
    (Issue #437)
  * `PSAnalysis` now includes Hausdorff pairs analysis and associated nearest
    neighbor plotting method (Issue #438)
  * New class `PSAPair` added to MDAnalysis.analysis.psa for handling
    Hausdorff pairs analysis (Issue #438)
  * `PSAnalysis` can now generate annotated heat maps (Issue #438)
  * Added three new distance functions to MDAnalysis.analysis.psa (Issue #438)
  * Additional getters added to `Path` and `PSAnalysis` (Issue #438)
  * MSD matrix function now globally available in MDAnalysis.analysis.psa
    (Issue #438)
  * Function for obtaining coordinate axes from numpy trajectories now
    globally available in MDAnalysis.analysis.psa (Issue #438)
  * TPR parser updated for Gromacs 5.0.x and 5.1 (Issue #456)
  * Setup.py now looks for some configuration values in a config file. Each
    config option can also be changed via environment variables if they are
    prefixed with 'MDA_'. Current options are 'use_cython', 'use_openmp', 'debug_cflags'

Changes
  * An AtomGroup with 0 atoms now yields an `IndexError` on call to
    `AtomGroup.write` (Issue #434)
  * `PSA` changed to `PSAnalysis` to reduce namespace clutter (Issue #438)
  * To build with debug-symbols use 'MDA_DEBUG_CFLAGS' instead of 'MDA_DEBUG_CFLAGS'

Fixes
  * Fixed minor issue in lib.mdamath.make_whole where if all bonds
    were correctly sized, it wouldn't notice that multiple fragments
    had been given. (Issue #445)
  * Fixed issue with PDB Topology parsing where if serials went
    over 100k, they wrapped to '***', breaking the parser (Issue #446)
  * Fixed AtomGroup.sequence() (Issue #451)
  * Fixed PrimitivePDBParser not detecting when resids had looped over
    10,000. The original resid is stored as Atom.resnum (Issue #454)
  * Fixed TPR topology parser to treat all bonded interactions available in
    Gromacs 5.1 (Issue #222 and #352, pull request #463).

## AUTHORS

kain88-de, richardjgowers, dotsdl, sseyler, orbeckst, jbarnoud