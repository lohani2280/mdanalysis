MDAnalysis 0.11.0 is a substantial release that unifies much of the code base at the cost of some necessary changes to the user API. The key points for adjusting your code accordingly are outlined below. The details are in the [Release Notes for 0.11.0](ReleaseNotes0110).

Note that we also plan to have a conversion script `ten2eleven.py` available (see [#377](/MDAnalysis/mdanalysis/issues/377)) that should do some of the necessary changes.

## How to change your code
### Removed contact matrix progress meter
The `MDAnalysis.analysis.distances.contact_matrix` function does not show progress anymore. Remove the keywords `progress_meter_freq` and `suppress_progmet` from your code:
```python
import MDAnalysis
import MDAnalysis.analysis.distances
import numpy as np
#300 rows of random xyz array data, to simulate coordinates:
random_array_coords = np.float32(np.random.rand(300,3)) 

#before 0.11 release:
#contact_matrix = MDAnalysis.analysis.distances.contact_matrix(random_array_coords, returntype = "sparse", progress_meter_freq=10, suppress_progmet=True) 
#after 0.11 release:
contact_matrix = MDAnalysis.analysis.distances.contact_matrix(random_array_coords, returntype = "sparse") 
```

### New Timestep behaviour
Previously, Timesteps could be initiated with either
 - integer (allocated to this size)
 - another Timestep (copied it)
 - coordinates (created a Timestep of this size and filled it)

Now Timesteps can **only** be initiated with an integer argument.

To create a Timestep from another Timestep, or from an array of coordinates, use the new `from_timestep` and `from_coordinates` class methods:

``` python
from MDAnalysis.coordinates.base import Timestep

ts = Timestep(10)
ts = Timestep.from_timestep(other_ts)
ts = Timestep.from_coordinates(coordinates)
# with optional velocities and/or forces
ts = Timestep.from_coordinates(coordinates, velocities=velocities, forces=forces)
```

### atomgroup methods to properties and select_atoms replacing selectAtoms

```python
import MDAnalysis
from MDAnalysis.tests.datafiles import GRO, XTC
universe = MDAnalysis.Universe(GRO, XTC)

#before 0.11
#all_selection = universe.selectAtoms('all')
#after 0.11
all_selection = universe.select_atoms('all')

#before 0.11:
#all_selection.residues()
#after 0.11:
all_selection.residues

#before 0.11:
#all_selection.charges()
#after 0.11:
all_selection.charges
```

### `Atom.number` renamed to `Atom.index`
```python
import MDAnalysis
from MDAnalysis.tests.datafiles import GRO, XTC
universe = MDAnalysis.Universe(GRO, XTC)

atomgroup = universe.select_atoms('all')
first_atom = atomgroup[0]

#before 0.11:
#first_atom.number
#after 0.11:
first_atom.index
```

### Migration of some code to `MDAnalysis.lib`

```python
import MDAnalysis
import numpy

random_coord_array_1 = numpy.float32(numpy.random.rand(20,3))
random_coord_array_2 = numpy.float32(numpy.random.rand(20,3))

#before 0.11:
#import MDAnalysis.core.distances
#distance_array = MDAnalysis.core.distances.distance_array(random_coord_array_1, random_coord_array_2)
#after 0.11: 
import MDAnalysis.lib.distances
distance_array = MDAnalysis.lib.distances.distance_array(random_coord_array_1, random_coord_array_2)
```

### `Timestep._x` `_y` and `_z` are now read only

`Timestep._x` cannot be assigned to, but can be changed in place.  This is to ensure that `_x` remains strictly a view of the x coordinate of the position data.

``` python
# previously:
ts._x = stuff  # would break view onto `_pos`

ts._x[:] = stuff  # still works to assign in place
ts._pos[:,0] = stuff
```


### The `fullgroup` selection keyword is now deprecated
0.11 introduces the `global` selection modifier keyword, which, among other cases, can be used to replace `fullgroup` when combined with `group`. You need only change `fullgroup` to `global group` in your selections:

```python
#before 0.11:
solvent = universe.selectAtoms("name SOL")
solvating = solvent.selectAtoms("around 5 fullgroup ref", ref=some_complex_selection)
#after 0.11: 
solvent = universe.select_atoms("name SOL")
solvating = solvent.select_atoms("around 5 global group ref", ref=some_complex_selection)
```

### `numberOfAtoms()` to `n_atoms` and others
```python
import MDAnalysis
from MDAnalysis.tests.datafiles import GRO, XTC
universe = MDAnalysis.Universe(GRO, XTC)
all_selection = universe.select_atoms('all')

#before 0.11:
#atom_count = all_selection.numberOfAtoms()
#after 0.11:
atom_count = all_selection.n_atoms

#before 0.11:
#residue_count = all_selection.numberOfResidues()
#after 0.11:
residue_count = all_selection.n_residues

#before 0.11:
#segment_count = all_selection.numberOfSegments()
#after 0.11:
segment_count = all_selection.n_segments

#before 0.11:
#frame_count = universe.trajectory.numframes
#after 0.11:
frame_count = universe.trajectory.n_frames
```

### Renamed topology methods

#### Working with topology from AtomGroups

Manipulating AtomGroups as items of topology (bonds, angles or torsions) has been reworked.  `AtomGroup.bond` is now a property which returns a `Bond` object.

``` python
ag = u.atoms[:2]  # 2 size AtomGroup

## Before!
# ag.bond()  # returned the length of the bond
## Now!
ag.bond.value()

```

This now allows for groups of `Bond`s to be collected together in a `TopologyGroup`

``` python
ag1 = u.atoms[0] + u.atoms[10]
ag2 = u.atoms[11] + u.atoms[20]
ag3 = u.atoms[21] + u.atoms[30]

tg = ag1.bond + ag2.bond + ag3.bond
tg.values()  # returns the length between each of the three pairs of atoms
```

The above is identical for `angle` `dihedral` and `improper_dihedral`

Working with dihedrals

``` python
# Previously
ag = u.atoms[:4]  # make a size 4 AtomGroup
ag.dihedral()  # returned the angle

# Now
ag = u.atoms[:4]  # select 4 atoms as usual
di = ag.dihedral  # convert the AtomGroup to a Dihedral object
di.value()  # returns the size of the angle
```

#### Torsions are now called dihedrals throughout MDAnalysis

To avoid confusion between the names dihedral and torsion, the term dihedral is used throughout the package.

``` python
# Previously
ag = u.atoms[:100]
tg = ag.torsions  # returned a TopologyGroup of all dihedrals/torsions

# Now
ag = u.atoms[:100]
tg = ag.dihedrals  # returns a TopologyGroup of all dihedrals/torsions


```

Also, the function "calc_torsions" was deprecated and renamed to "calc_dihedrals"

``` python
# Previously
from MDAnalysis.lib.distances import calc_torsions
t = calc_torsions(ag1.positions, ag2.positions, ag3.positions, ag4.positions)

# Now
from MDAnalysis.lib.distances import calc_dihedrals
t = calc_dihedrals(ag1.positions, ag2.positions, ag3.positions, ag4.positions)

```

### Frame numbering is now 0-based
```python
import MDAnalysis
from MDAnalysis.tests.datafiles import GRO, XTC
universe = MDAnalysis.Universe(GRO, XTC)

#before 0.11:
if 1 == universe.trajectory.frame:
    print 'currently on first frame'

#after 0.11:
if 0 == universe.trajectory.frame:
    print 'currently on first frame'
```