** Introduction **

The topology data defines the structure of the MDAnalysis Universe, and is the first structure generated on initialising the Universe.
The generation of this structure is done by a TopologyReader object, which generally reads this information from a file.  The first argument in Universe init is the topology file.

** Rules **

A topology reader subclasses the class TopologyReader, and then must define the method ``parse``
The parse method must return a dictionary of fields, as defined in [Topology Data Structures](https://github.com/MDAnalysis/mdanalysis/wiki/TopologyDataStructures)


** Examples **

For an example of a very minimalistic TopologyReader, see the XYZParser

For an example of a TopologyReader populating many fields, see the PSFParser