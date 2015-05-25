# Introduction #

The topology data defines the structure of the MDAnalysis Universe, and is the first structure generated on initialising the Universe.
The generation of this structure is done by a TopologyReader object, which generally reads this information from a file.  The first argument in Universe `init` is the topology file.

## Using a new parser

It is possible to "inject" a new topology parser into MDAnalysis by supplying it as a keyword argument to
the Universe.  This is done via the keyword `topology_format`.  
This keyword is typically used to override the format detection based on file suffix, but if a subclass
of `TopologyReader` is passed, this is used as the parser.
This is designed to allow the development and testing of new parsers without having to edit the source
of MDAnalysis.

Of course, MDAnalysis prides itself on handling many formats, and so any new topology parser 
would be a very welcome addition into the core of MDAnalysis!

# Rules #

A topology reader subclasses `MDAnalysis.topology.base.TopologyReader`, and then must define the method ``parse``.
The parse method must return a dictionary of fields, as defined in [Topology Data Structures](https://github.com/MDAnalysis/mdanalysis/wiki/TopologyDataStructures)

It is not required to define an `__init__`, this is done in the base class.

The base class (via `MDAnalysis.coordinates.base.IOBase`) also provides context management too.
It is important that file access is done via context managers, so that the file handle is closed
on both success and failure.

Parsers should only raise either `IOError` or `ValueError`.  `IOError` is typically for 
problems with file access, while `ValueError` caters for nonsensical values once the file
has been read.

# Examples #

``` python
from MDAnalysis.topology.base import TopologyParser
from MDAnalysis.core.AtomGroup import Atom

Class MyParser(TopologyParser):
    def parse(self):
        # the filename to parse is accessible via self.filename
        # any kwargs were saved into self.kwargs

        # the return type needs to be a dict!
        struc = dict()
        # required:
        struc['atoms'] = self._parse_atoms()  # required
        # optionally:
        struc['bonds'] = self._parse_bonds()
        # etc for all fields you want to capture within the Topology

        return struc

    def _parse_atoms(self):
        # atoms need to be a list of Atom instances
        # for convenience, the Universe that called the parser is
        # accessible via self._u within the TopologyReader
        a = Atom(universe=self._u)

        return (A,)

    def _parse_bonds(self):
        # bonds will be a tuple of 2 length tuples
        # this example is a bond between Atoms 0 and 1, 1 and 2, and 2 and 3
        # Note that this is zero based!
        return ((0, 1), (1, 2), (2, 3))


```

For an example of a very minimalistic TopologyReader, see the XYZParser

For an example of a TopologyReader populating many fields, see the PSFParser