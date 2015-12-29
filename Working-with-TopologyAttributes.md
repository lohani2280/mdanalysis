TopologyAttributes (TA) are used to define static information in the Universe.



## Creating a new TopologyAttribute

Should a Topology require a TopologyAttribute that doesn't currently exist, these can be easily made....


## The transplant system

All attributes added to the Topology object are automatically added to the namespace of all Groups (AtomGroup, ResidueGroup and SegmentGroup) and components (Atom, Residue & Segment).
It is also possible to make a TA supply arbitrary methods and properties to these classes, done via the transplant system.

To do this, the desired method is defined inside the TopologyAttribute.
A **defaultdict** called transplants is defined, and a tuple containing the method's name and the method is appended to the list of methods.
Transplants may be defined at either the `atom`, `residue`, `segment`, `atomgroup`, `residuegroup` or `segmentgroup` levels.

``` python
class MyAttribute(TopologyAttribute):
    # <normal TA code>
    transplants = defaultdict(list)

    def mass_multiplier(ag, value):
        """Multiply the mass of each item by *value*"""
        return ag.masses * value

    transplants['atomgroup'].append(('multimass', mass_multiplier))


# Creates the following method:
u.atoms.multimass(5)
```

It is also possible to add in properties as well as methods:

``` python
class MyAttribute(TopologyAttribute):
    # <normal TA code>
    transplants = defaultdict(list)

    def mass_doubler(ag):
        """Returns double the mass of each item"""
        return ag.masses * 2

    transplants['atomgroup'].append(('doublemass', property(mass_doubler, None, None, mass_multiplier.__doc__)))


# Creates the following method:
u.atoms.doublemass
```