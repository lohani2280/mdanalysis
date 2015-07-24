MDAnalysis 0.11 is a substantial release that unifies much of the code base at the cost of some necessary changes to the user API. The key points for adjusting your code accordingly are outlined below.

### Suppressing contact matrix progress meter (new `quiet` keyword)
```
import MDAnalysis
import MDAnalysis.analysis.distances
import numpy as np
#300 rows of random xyz array data, to simulate coordinates:
random_array_coords = np.float32(np.random.rand(300,3)) 

#before 0.11 release:
#contact_matrix = MDAnalysis.analysis.distances.contact_matrix(random_array_coords, suppress_progmet=False) 
#after 0.11 release:
contact_matrix = MDAnalysis.analysis.distances.contact_matrix(random_array_coords, quiet=False) 
```