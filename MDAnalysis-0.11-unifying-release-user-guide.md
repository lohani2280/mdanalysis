MDAnalysis 0.11.0 is a substantial release that unifies much of the code base at the cost of some necessary changes to the user API. The key points for adjusting your code accordingly are outlined below. The details are in the [Release Notes for 0.11.0](ReleaseNotes0110).

## How to change your code
### Suppressing contact matrix progress meter (new `quiet` keyword)
```python
import MDAnalysis
import MDAnalysis.analysis.distances
import numpy as np
#300 rows of random xyz array data, to simulate coordinates:
random_array_coords = np.float32(np.random.rand(300,3)) 

#before 0.11 release:
#contact_matrix = MDAnalysis.analysis.distances.contact_matrix(random_array_coords, returntype = "sparse", progress_meter_freq=10, suppress_progmet=True) 
#after 0.11 release:
contact_matrix = MDAnalysis.analysis.distances.contact_matrix(random_array_coords, returntype = "sparse", progress_meter_freq=10, quiet=True) 
```

### New Timestep behaviour

Previously, Timesteps could be initiated with either
 - integer (allocated to this size)
 - another Timestep (copied it)
 - coordinates (created a Timestep of this size and filled it)

Now Timesteps can **only** be initiated with an integer argument.

To create a Timestep from another Timestep, or from an array of coordinates, use the new `from_timestep` and `from_coordinates` class methods

``` python
from MDAnalysis.coordinates.base import Timestep

ts = Timestep(10)

ts = Timestep.from_timestep(other_ts)

ts = Timestep.from_coordinates(coords [, velocities=velos, forces=forces])

```