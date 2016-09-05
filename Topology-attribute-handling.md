As per the discussion started in [#942](https://github.com/MDAnalysis/mdanalysis/issues/942), missing attributes from a topology will be handled according to different categories (a few more categories were added in this page). These are:

Classes
-------

1. Left out of the Topology, with a general `AttributeError` on accession;
2. Left out of the Topology, but with a specific `MDAnalysis.MissingTopologyAttributeError` on accession, perhaps suggesting how to set/guess values;
3. Added to the topology, set to a default;
4. Added to the topology, guessed — possibly with a warning that it is guessed. If guess fails, assign some default;
5. Added to the topology, guessed — possibly with a warning that it is guessed. If guess fails, warn and leave unassigned (like category 2 above);
6. Absolutely necessary; abort Topology/Universe loading if absent;
7. Internal: attribute generated internally that does not depend on topology availability or guessing;
8. *Attribute that you think should be removed from consideration.*

(Thanks to @jbarnoud for suggesting some of the possibilities above)

Poll
----

We need to decide which attributes go into which category. This will be implemented as part of the topology refactor ([#363](https://github.com/MDAnalysis/mdanalysis/issues/363)). Please edit the table below to vote your opinion, and if needed add your justification to [#942](https://github.com/MDAnalysis/mdanalysis/issues/942). Factors you might want to take into account are
- consistency,
- ease for the user in case of missing attributes,
- little boilerplate for the most common cases,
- and backwards compatibility (though this last factor might be of less importance since compatibility will already be broken by the topology refactor).

|Attribute    |@mnmelo| 
|------------:|:-----:| 
|residues*    |   3   |
|segments*    |   3   |
|indices      |   7   | 
|resindices   |   7   | 
|segindices   |   7   | 
|ids          |   3   | 
|names        |   6   | 
|types        |   1   | 
|elements     |   5   | 
|radii        |   5   | 
|chainIDs     |   1   | 
|icodes       |   1   | 
|tempfactors  |   1   | 
|masses       |   5   | 
|charges      |   2   | 
|bfactors     |   1   | 
|occupancies  |   1   | 
|altLocs      |   1   | 
|resids       |   3   | 
|resnames     |   2   | 
|resnums      |   3   | 
|segids       |   3   | 
|bonds        |   5   | 
|angles       |   5   | 
|dihedrals    |   5   |
|impropers    |   1   | 
|*other attrs*|   1   |

\* While not strictly attributes, Residue/Segment assignment falls pretty much in the same categories, and their handling is relevant for related attributes (resindices, resids, etc.).