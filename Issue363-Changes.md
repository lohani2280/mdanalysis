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