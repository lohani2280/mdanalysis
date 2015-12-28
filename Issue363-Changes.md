Note: Issue-363 is a work in progress, nothing is final yet, everything is up for debate!

[Issue #363](https://github.com/MDAnalysis/mdanalysis/issues/363) brings many changes, these are summarised here:

## All Groups are no longer subclassed from AtomGroup.  

### Explanation

Previously, every object except Atom subclassed from AtomGroup.  This meant that calling `.positions` of would give you the positions of the Atoms contained within that group.

Previous class structure:
```
Atom

AtomGroup  -> Residue
           -> ResidueGroup -> Segment
                           -> SegmentGroup
```
New class structure:
```
Group    -> AtomGroup
         -> ResidueGroup
         -> SegmentGroup

Atom
Residue
Segment
```

Now each object only contains information pertaining to that particular object.  A `Residue` object only yields information about the Residue, to get to the atoms, use `Residue.atoms`

### Why this was changed

Previously everything inheriting from AtomGroup made it unclear at what level of topology a given method or attribute was working on.  Ie does `ResidueGroup.charges` give the charge of the Residues or the Atoms?
Also, it was unclear what size a given output would be (see [Issue-411](https://github.com/MDAnalysis/mdanalysis/issues/411))

### How to work around this

To access Atom level information from anything that isn't an AtomGroup, use the `.atoms` level accessor.
For example, changing all `.positions` calls on anything that isn't an `AtomGroup` to `.atoms.positions`.

## Each Atom is a member of exactly one Residue. Each Residue is a member of exactly one Segment.
The new `Topology` object keeps an array giving the residue membership of each atom. Getting the resname of the residue of a group of atoms, then, is achieved by taking the indices of these atoms to fancy-index the `Atoms->Residues` array, and then using the result of this to fancy-index the `Resnames` array. For example, if the  `Topology` has 5 atoms and 3 residues, with membership (`Atoms->Residues`) and `Resnames` arrays as below:

```
       Atoms->Residues           Resnames
 index ---------------     index --------
     0 0                       0 GLU
     1 2                       1 LYS
     2 1                       2 ALA
     3 1
     4 2
```

calling `AtomGroup.resnames` for an `AtomGroup` with atoms [2, 0, 1, 2] will yield (pseudocode):

```
"Atoms->Residues"[[2, 0, 1, 2]] --> [1, 0, 2, 1]
"Resnames"[[1, 0, 2, 1]]        --> ['LYS', 'GLU', 'ALA', 'LYS']
```

This scheme only works if each atom is a member of one and only one residue. Likewise, residues are members of one and only one segment. Furthermore, `AtomGroup`s, `ResidueGroup`s, and `SegmentGroup`s are very thin, storing only the indices of their members as a `numpy` array. This gives a number of advantages:

1. **Performance**. We get at least an 8x speedup over the old scheme when accessing attributes. Setting attributes can give up to a 40x speedup.
2. **Memory**. We don't store, for example, a resname for each atom, but instead store attributes at the level they make sense for.
3. **Consistency**. Since attributes are stored in one place, we avoid cases where the topology is in an inconsistent state, e.g. two atoms in the same residue give a different resname.
4. **No staleness**. Because e.g. `ResidueGroups` are only an array of indices, not a list of `Residue` objects generated upon creation of the group, changes of resiude-level properties by another `ResidueGroup` are always reflected consistently by every other one. Data is not duplicated anywhere in this scheme, and is all contained in the `Topology` object.

For further performance comparisons, check out this [notebook](https://gist.github.com/dotsdl/0e0fbd409e3e102d0458).