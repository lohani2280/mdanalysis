MDAnalysis 0.14.0 is a new release with many new features and fixes; see also the [blog post](http://blog.mdanalysis.org).

* release: 0.14.0
* release date: 02/28/16 

## Authors
tyler.je.reddy, kain88-de, jbarnoud, richardjgowers, orbeckst, manuel.nuno.melo, Balasubra, Saxenauts, mattihappy

## API Changes

  * Offsets files for Gromacs trajectory formats have been changed to a numpy style format '.npz'. Offsets files will be regenerated when you load a    xtc/trr trajectory again. (Issue #441)
  * rotation_matrix now accepts array-likes as input

## Enhancement

  * XDR file seeking errors now report system errno. (PR #678)
  * Offsets reading for xtc/trr files has been sped up. (Issue #441)
  * select_atoms now implicitly ORs multiple values after a keyword for    many types of selections (Issue #345)
  * Performance improvements for the PDBReader of about 10%
  * LinearDensity analysis module added, which allows to compute linear mass    and charge density profiles along the three cartesian axes of the cell.    Works for orthorombic, fixed volume cells only. (Issue #670)
  * Trajectories can now be sliced using a boolean array (Issue #725)

## Changes

  * xdrlib was rebranded libmdaxdr. (Issue #679)
  * xdrlib has been ported to cython. (Issue #441)
  * util.NamedStream no longer inherits from basestring (Issue #649)
  * Short TRZ titles are striped from trailing spaces. A friendlier error    message is raised when the TRZ writer is asked to write a title longer    than 80 characters. (Issue #689)
  * PDB doesn't save charge information anymore
  * coordinates.core.get_writer_for uses the user defined format if provided     before trying to deduce the format from file extension. (Issue #712)

## Fixes
  * Syntax error corrected in psa.py (Issue #738)
  * XDR file seeking and telling working again for large files (Issue #677).
  * ContactAnalysis1 run method now starts at frame index 0 by default (Issue #624)
  * Fixed PrimitivePDBWriter alignment of the atom name. (Issue #639 and #647)
  * The 'atom' selection keyword returns an empty selection rather than an    error when no atoms are are found. (Issue #644)
  * nucleic selection will now detect nucleic residue names suffixed with 3 or 5    (Issue #461)
  * Fixed Reader returning all frames with stop in slice was 0 (Issue #672)
  * Fixed NCDFReader not reading dt (Issue #676)
  * Fixed PDB-Topology read bonds for atom ids larger then 10000 (Issue #693)
  * Fixed Type Error in qcprot.pyx when no rotation can be fond (Issue #705)
  * Fixed cyzone selection failing in orthogonal systems (Issue #710)
  * Fixed Error in calculation of average grid density (Issue #716)
  * Fixed indexing an AtomGroup using a list of bools now working (Issue #729)