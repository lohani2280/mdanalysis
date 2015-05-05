This is a minor new release of MDAnalysis, which contains a couple of [enhancements](#Enhancements) and an important [bug fix](#Fixes).

## CHANGELOG ##

  * 0.9.1
  * released 27 March 2015

### Enhancements ###

  * XYZ file format can be used without an associated topology file.
  * GAMESS output files can be read as trajectories for calculations of type ```SURFACE'' and ```OPTIMIZE'' (work wit both GAMESS-US and Firefly)

### Changes ###

  * removed undocumented MDAnalysis.builder module

### Fixes ###

  * TRR coordinate access via _y and_z now works properly ([Issue 224](http://issues.mdanalysis.org/224))

## Authors ##
manuel.nuno.melo, richardjgowers, comconadin
