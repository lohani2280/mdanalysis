Release 0.12.1 brings a few minor fixes, including better OpenMP detection on Linux/OSX

## CHANGELOG

Fixes
  * Fixed OpenMP detection on Linux/OSX #459
  * Fixed reading of LAMMPS trajectory times: default unit ought
    to be fs and not ps
  * Fixed setting of `dt` for `DCDReader` (and LAMMPS DCDReader) with
    keyword argument `Universe(..., dt=<dt>)`
  * Fixed a bug in `topology.core.guess_atom_element` where a
    single digit atom name would raise an IndexError (#476)
  * Fixed numpy -> np in `LeafletFinder`

## AUTHORS

kain88-de, orbeckst, richardjgowers