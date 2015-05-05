# Changes #

  * 0.6.3 release
  * minimum requirement is <strike>python 2.4</strike> python 2.5 (using with\_statement in the analysis module and we have not tested on 2.3 in a while)
  * requires Biopython (www.biopython.org)
  * analysis modules (MDAnalysis.analysis):
    * lipid bilayer leaflet detection
    * native contact analysis ("q1-q2")
    * rms-fitting based on sequence alignment
  * write selections for other codes from AtomGroups (VMD, pymol, CHARMM,  Gromacs ndx)
  * gro reader ([Issue 31](http://issues.mdanalysis.org/31))
  * better API for loading a topology and a coordinate file in Universe()
  * trajectory reader: DCDReader can reverse trajectory with negative step increment; XTC/TRRReader can do simple (forward) slices by doing (slow!) sequential iteration
  * deprecated principleAxes() and introduced principalAxes() with less confusing return values ([Issue 33](http://issues.mdanalysis.org/33)).
  * fixed wrong unitcell dimensions for XTC/TRR ([Issue 34](http://issues.mdanalysis.org/34))
  * added basic XYZ reader with compression support ([Issue 35](http://issues.mdanalysis.org/35))
  * PDB reader guesses masses (unknown elements are set to 0)

# Authors #

orbeckst, denniej0, danny.parton, philipwfowler
