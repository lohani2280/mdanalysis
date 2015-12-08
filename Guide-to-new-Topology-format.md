The new (as of Dec 2015) topology format requires that topology parsers return a single Topology object, where previously they returned a list of Atom instances.

Below is an annotated version of the parse for GRO files.

The `Topology` object contains many tables of data, called `Attributes`, which `*Groups` later refer to.
For example, `Atomnames` is an `Attribute`.
`Attributes` facilitate access to arrays of data.  For an `AtomGroup` accessing the names, this simply applies `AtomGroup.indices` as a mask to the array, but it also handles the `ResidueGroup` and `SegmentGroup` cases.

``` python
import numpy as np

from ..lib.util import openany
from ..core.topologyattrs import (
    Resids,
    Resnames,
    Atomids,
    Atomnames,
)
from ..core.topology import Topology
from .base import TopologyReader, squash_by
```

All Parsers must still inherit from mda.topology.base.TopologyReader
and define a parse method
``` python
class GROParser(TopologyReader)
    def parse(self):
        """Return the *Topology* object for this file"""
        # Gro has the following columns
        # resid, resname, name, index, (x,y,z)
        with openany(self.filename, 'r') as inf:
            inf.readline()
            n_atoms = int(inf.readline())

            # Allocate arrays
            # This isn't necessary, appending to lists is also
            # possible where the number of atoms isn't known beforehand.
            # But ultimately all data must be stored in np.arrays
            resids = np.zeros(n_atoms, dtype=np.int32)
            # String data must use dtype=object so that any size string
            # can be stored.
            resnames = np.zeros(n_atoms, dtype=object)
            names = np.zeros(n_atoms, dtype=object)
            indices = np.zeros(n_atoms, dtype=np.int32)

            for i in xrange(n_atoms):
                line = inf.readline()
                resids[i] = int(line[:5])
                resnames[i] = line[5:10].strip()
                names[i] = line[10:15].strip()
                indices[i] = int(line[15:20])
```

Whereas previously each `Atom` kept a record of its `resname`, all Residue properties are now stored on a per-Residue basis, meaning that only **one** record of Residue name is ever kept.

Considering a trivial case of 4 atoms:
 - Atom 1, Resid 3, Resname 'A'
 - Atom 2, Resid 3, Resname 'A'
 - Atom 3, Resid 4, Resname 'B'
 - Atom 4, Resid 4, Resname 'B'

The `squash_by` function reduces this information to the following:

 - Atom 1, Resindex 0
 - Atom 2, Resindex 0
 - Atom 3, Resindex 1
 - Atom 4, Resindex 1

Represented in the array `residx`

And then a separate record of information per-Residue

 - Resindex 0, Resid 3, Resname 'A'
 - Resindex 1, Resid 4, Resname 'B'

Represented in the arrays `new_resids` and `new_resnames`

``` python
        residx, new_resids, (new_resnames,) = squash_by(resids, resnames)

        # new_resids is len(residues)
        # so resindex 0 has resid new_resids[0]
        atomnames = Atomnames(names)
        atomids = Atomids(indices)
        residueids = Resids(new_resids)
        residuenames = Resnames(new_resnames)
```

Finally, the `Topology` is created by declaring its size in `n_atoms`, `n_residues` and `n_segments`,
passing a list of the `Attributes`, and the relationship between atom indices and residue indices, and residue indices and segment indices.

``` python
        top = Topology(n_atoms, len(new_resids), 0,
                       attrs=[atomnames, atomids, residueids, residuenames],
                       atom_resindex=residx,
                       residue_segindex=None)

        return top
```