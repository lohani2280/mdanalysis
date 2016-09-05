As per the discussion started in [#942](https://github.com/MDAnalysis/mdanalysis/issues/942), missing attributes from a topology will be handled according to different categories (a few more categories were added in this page). These are:

# Categories of topology attributes

1. Left out of the Topology, with a general `AttributeError` on accession;
2. Left out of the Topology, but with a specific `MDAnalysis.MissingTopologyAttributeError` on accession, perhaps suggesting how to set/guess values;
3. Added to the topology, set to a default;
4. Added to the topology, guessed — possibly with a warning that it is guessed. If guess fails, assign some default;
5. Added to the topology, guessed — possibly with a warning that it is guessed. If guess fails, warn and leave unassigned (like category 2 above);
6. Absolutely necessary; abort Topology/Universe loading if absent;
7. Internal: attribute generated internally that does not depend on topology availability or guessing;
8. *Attribute that you think should be removed from consideration.*

(Thanks to @jbarnoud for suggesting some of the possibilities above)

# Poll

We need to decide which attributes go into which category. This will be implemented as part of the topology refactor ([#363](https://github.com/MDAnalysis/mdanalysis/issues/363)). Please edit the table below to vote your opinion, and add any comments below. Factors you might want to take into account are
- consistency,
- ease for the user in case of missing attributes,
- little boilerplate for the most common cases,
- and backwards compatibility (though this last factor might be of less importance since compatibility will already be broken by the topology refactor).

|Attribute    |@mnmelo|@orbeckst|@richardjgowers|
|------------:|:-----:|:-------:|:-------------:|
|residues*    |   3   |   3     | 3             |
|segments*    |   3   |   3     | 3             |
|indices      |   7   |   7     | 7             |
|resindices   |   7   |   7     | 7             |
|segindices   |   7   |   7     | 7             |
|ids          |   3   |   3     | 1             |
|names        |   2   |   2     | 2             |
|types        |   1   |   2     | 2             |
|elements     |   5   |   5     | 2             |
|radii        |   5   |   2     | 2             |
|chainIDs     |   1   |   1     | 1             |
|icodes       |   1   |   1     | 1             |
|tempfactors‡ |   1   |   1     | 1             |
|masses       |   5   |   5     | 5             |
|charges      |   2   |   2     | 2             |
|bfactors‡    |   1   |   1     | 1             |
|occupancies  |   1   |   1     | 1             |
|altLocs      |   1   |   1     | 1             |
|resids       |   3   |   3     | 3             |
|resnames     |   2   |   2     | 2             |
|resnums      |   3   |   3     | 2             |
|segids       |   3   |   3     | 3             |
|bonds        |   5   |   5     | 2             |
|angles       |   5   |   5     | 2             |
|dihedrals    |   5   |   5     | 2             |
|impropers    |   1   |   1     | 2             |
|*other attrs*|   1   |   1     | 1             |

\* While not strictly attributes, Residue/Segment assignment falls pretty much in the same categories, and their handling is relevant for related attributes (*resindices*, *resids*, etc.).

‡ *tempfactors* and *bfactors* are the same thing, we only need one.


# Comments

## @orbeckst
- I would also like a specific error message (category **2**) for anything that has a keyword in the search syntax, such as *type*.
- I am wavering on *radii* and I am not sure if they should really be guessed (cat **5**). We could assign Bondi van der Waals radii [A. Bondi. van der waals volumes and radii. The Journal of Physical Chemistry, 68(3):441–451, 1964. doi: 10.1021/j100785a001. URL http://pubs.acs.org/doi/abs/10.1021/j100785a001.] and some ionic radii (or perhaps [Born radii](https://github.com/Becksteinlab/BornProfiler/blob/develop/bornprofiler/templates/bornions.dat)) for ions but ultimately I am not sure how useful this is going to be. What use will people have of such ad-hoc radii? The main use for *radii* is as when reading/writing PQR files, where the radii have well-defined meaning. So maybe cat **2** would be better for *radii*.
- It would be nice of the "PDB" attributes (*chainID*, *iCodes*, *occupancy*, *bfactor* (which is the same as *tempfactor* so we should only have one and I vote for *bfactor* because of continuity) could be elevated to cat **2** because they are pretty common and we might reduce questions regarding "Why does this not work" if the error message just says right away "Load from a PDB file". I left them at **1**, though, in order to keep things simple.

## @richardjgowers

Bond guessing is currently not optimised for large systems, so guessing this by default is a **bad** idea for larger systems.  Maybe we could implement a size cutoff for auto-guessing (~10k atoms or whatever keeps Universe load times below ~1s?)

I'm not sure how an Atom id is different from an Atom index?  I understand that **some** systems may number irregularly, but this isn't necessary everywhere.

Maybe it's just my background, but I'd like things to work with coarse-grained systems just as well, so things like element, resnum and other atomistic based things aren't mandatory there.

## @mnmelo
My rationale was to have a few attributes with VIP status, like mass, charge, and element, which we can expect in a format-independent manner. These are either guessed or left missing if guessing fails. Bond guessing also follows this.

VIP attributes could also be made to be always present as an object, even if guessing fails. This way `u.atoms.masses` is always a present TopologyAttr, albeit one that can raise a `MissingTopologyAttributeError` on access if it couldn't be obtained from the topology nor guessed. 

Most other attributes can either be safely set to defaults, or, for the more format-specific ones, left missing.

I was under the impression that atom names are central to our machinery. I'm happy to demote them to **2.**
