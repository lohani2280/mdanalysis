MDAnalysis 0.15.0 is a new release with many new features and fixes; see also the [blog post](http://blog.mdanalysis.org).

* release: 0.15.0
* release date: 05/15/16 

# CHANGELOG
## Metadata

  * link download_url to GitHub releases so that Depsy recognizes
    contributors (issue #749)
  * a __version__ variable is now exposed; it is built by setup.py from the
    AUTHORS file (Issue #784)

## API Changes

  * rmsd doesn't superimpose by default anymore. The superposition
    is controlled by the 'superposition' keyword now. (see issue #562, #822)

## Enhancements

  * Add conda build scripts (Issue #608)
  * Added read-only property giving Universe init kwargs (Issue #292)
  * Added 'crdbox' as AMBER Trj format extension (Issue #846)
  * Iteration and seeking in PDB files made faster (Issue #848)

## Fixes

  * ENT file format added to PDB Readers/Writers/Parsers (Issue #834)
  * rmsd now returns proper value when given array of weights (Issue #814)
  * change_release now finds number and dev (Issue #776)
  * units.py now correctly prints errors for unknown units.
  * test_shear_from_matrix doesn't fail for MKL builds anymore (Issue #757)
  * HEADER and TITLE now appear just once in the PDB. (Issue #741) (PR #761)
  * MOL2 files without substructure section can now be read (Issue #816)
  * MOL2 files can be written without substructure section (Issue #816)
  * GRO files with an incomplete set of velocities can now be read (Issue #820)
  * Fixed Atom.position/velocity/force returning a view onto Timestep array
    (Issue #755)
  * PDB files can now read a CRYST entry if it happens before model headers
    (Issue #849)
  * Fixed HistoryReader returning 1 based frame indices (Issue #851)

## Changes

  * Added zero_based indices for HBondsAnalysis. (Issue #807)
  * Generalized contact analysis class `Contacts` added. (Issue #702)
  * Removed Bio.PDBParser and sloppy structure builder and all of
    MDAnalysis.coordinates.pdb (Issue #777)
  * PDB parsers/readers/writers replaced by "permissive"/"primitive"
    counterparts (formerly known as PrimitivePDBReader); the
    'permissive' keyword for Universe is now ignored and only the
    native MDAnalysis PDBReader is being used (Issue #777)
  * PDBReader only opens a single file handle in its lifetime,
    previously opened & closed handle each frame (Issue #850)

## Deprecations (Issue #599)

  * Use of PrimitivePDBReader/Writer/Parser deprecated in favor of PDBReader/
    Writer/Parser (Issue #777)
  * Deprecated all `get_*` and `set_*` methods of Groups.
  * Deprecation warnings for accessing atom attributes from Residue,
    ResidueGroup, Segment, SegmentGroup. Will not be present or will
    give per-level results.
  * Deprecation warnings for accessing plural residue attributes from
    Residue or Segment (will disappear), or from SegmentGroup (will give
    per-Segment results).
  * Deprecation warnings for accessing plural segment attributes from Segment
    (will disappear).
  * Deprecated Atom number, pos, centroid, universe setter
  * Deprecated AtomGroup serials, write_selection
  * Deprecated Residue name, id
  * Deprecated Segment id, name
  * Deprecated as_Universe function; not needed
  * Deprecated ContactAnalysis and ContactAnalysis1 classes

# Authors
jandom, abhinavgupta94, orbeckst, kain88-de, hainm, jbarnoud, dotsdl, richardjgowers, BartBruininks, jdetle
