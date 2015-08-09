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