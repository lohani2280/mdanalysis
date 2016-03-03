<img src="https://developers.google.com/open-source/gsoc/images/gsoc2016-sun-373x373.png" title="Google Summer of Code 2016" alt="Google Summer of Code 2016" align="right"/>
A list of projects ideas for [[Google Summer of Code 2016|Google-Summer-Of-Code]].

The project ideas can be roughly categorized as

1. [[New analysis functionality|GSoC-2016-Project-Ideas#new-analysis-functionality]]
2. [[Increasing performance|GSoC-2016-Project-Ideas#increasing-performance]]
3. [[New input formats|GSoC-2016-Project-Ideas#new-input-formats]]
4. [[Increase platform availability|GSoC-2016-Project-Ideas#increase-platform-availability]]
5. [[Increase ease-of-use|GSoC-2016-Project-Ideas#increase-ease-of-use]]
6. [[Improve the library core|GSoC-2016-Project-Ideas#improve-the-library-core]]

**Or work on your your own idea!** Get in contact with us to propose an idea and we will work with you to flesh it out into a full project. Raise an issue in the [Issue Tracker](/MDAnalysis/mdanalysis/issues) or contact us via the [developer Google group](http://developers.mdanalysis.org).

------

# New analysis functionality
## Implement a general dimension reduction algorithm

**Difficulty**: Hard

**Mentors**: Max, Richard, Manuel

MDAnalysis already comes with a range of different standard analysis tools but currently lacks an implementation of a general dimension reduction algorithm, that can select an arbitrary number of dimensions of interest. 3 common general techniques are

- [Time Independent Component Analysis](http://arxiv.org/abs/1302.6614)
- [Diffusion Maps] (http://arxiv.org/abs/1506.06259)
- [Principle Component Analysis](https://en.wikipedia.org/wiki/Principal_component_analysis)

There are python implementations for all of these algorithms but none of them currently work with MDAnalysis out of the box. This is because the current python implementations work on normal numpy arrays that stores a complete trajectory in memory, but MDAnalysis never loads the whole trajectory but only one frame at a time. This approach allows MDAnalysis to treat very large system on a normal laptop or workstation. A new dimension reduction should be implemented as a class and inherit from [analysis.base](https://github.com/MDAnalysis/mdanalysis/blob/5b6471d93a36581d06ec73a1a0bddc8a460d4213/package/MDAnalysis/analysis/base.py#L35).

Of course you can also suggest us another dimension reduction algorithm that you would like to implement.

# Increasing performance

## Improve distance search 

**Difficulty**: Hard

**Mentors**: Max, Richard, Manuel

To analyze molecular simulations it is often helpful which atoms are close to each other. For this we calculate distance matrices where the distances between every atom pair is calculated. This is a very expensive operation that grows quadratic with the number of atoms involved.

Since we are only interested in atoms that are close to each other we can use some algorithms run faster
after some initial analysis of the coordinates. One class of these algorithms are domain-decomposition algorithms. The basic idea of this type of algorithms is to decompose the volume occupied by the atoms into different cells and then only calculate distances for atoms in neighboring cells. If atoms are not in neighboring cells we already know that the distance is to big for us to be interesting. A theoretical description of these algorithm can be found in [this book Appendix F](http://www.amazon.de/Understanding-Molecular-Simulation-Applications-Computational/dp/0122673514%3FSubscriptionId%3DAKIAILSHYYTFIVPWUY6Q%26tag%3Dduckduckgo-ffnt-de-21%26linkCode%3Dxm2%26camp%3D2025%26creative%3D165953%26creativeASIN%3D0122673514)

One domain decomposition algorithm is [cell grids](https://github.com/richardjgowers/cellgrid).

In this project you would integrate the cell grid algorithm into MDAnalysis. 

# New input formats

## Add new MD-Formats

**Dificulty**: Medium

**Mentors**: Max, Richard, Manuel

One of the strengths of MDAnalysis is its ability to support a wide range of different MD-formats. But we are still missing some like the new [TNG file format](http://onlinelibrary.wiley.com/doi/10.1002/jcc.23495/abstract) from Gromacs or [H5MD](https://github.com/pdebuyl/pyh5md). Alternatively, you can also add a format that you want to use personally in MDAnalysis.
This project will familiarize you with working with and connecting different APIs,
as well as giving insight into how modern portable data storage file formats work.

## Random Walk Trajectory Backend

**Difficulty**: Hard

**Mentors**: Max, Richard, Manuel

To check if a new analysis-method works as intended it is often a good idea to use it with a random walk in different simple energy landscapes (A flat energy, harmonic well, double well). In this project you would develop a 'Reader' that produces random trajectories. 

For analysis of molecular data a comparison against random data can be very useful for several reasons. The first is that we want to test if our analysis can distinguish between a simulation and random noise.
It can also be interesting to see what general analysis methods like Principle Component Analysis produce with [random data](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.62.8438).

The first random trajectory generator would just be a random walk in 3N dimensions (N is the number of particles in the simulation to compare to). The second would be to implement langevin dynamics in either predefined energy landscapes and/or arbitrary ones. Langevin dynamics in a energy landscape are close to the conformational dynamics of proteins, see [1]. As a first start you could implement a integrator for Langevin dynamics and later have the trajectory 'reader' use the integrator to dynamically generate the trajectory.

Please note that this project does require a background in statistical phyics or mathematics.

[1] Robert Zwanzig. Nonequilibrium statistical mechanics. Oxford University Press,
2001

# Increase platform availability

## Help port MDAnalysis to Python 3

**Difficulty**: Easy

**Mentors**: Max, Richard, Manuel

Python 3 is getting adopted by a wider range of users and unix distributions are starting to switch.
MDAnalysis can't run right now under Python 3 mostly due to it's [C/Cython extensions](https://github.com/MDAnalysis/mdanalysis/wiki/List-of-extensions), we currently try to move our C-extensions to cython which supports Python 2 and 3 with one source. See also [#260](https://github.com/MDAnalysis/mdanalysis/issues/260).

Missing here right now is the DCD trajectory readers. There exists an [incomplete work](https://github.com/MDAnalysis/mdanalysis/pull/682) to enable Python 2/3 of the DCD reader. In this project you would finish this work by either writing finishing this work or by rewriting the DCD interface in cython.

The second part of this project is to remove all other incompatibilities with Python 3 we currently have. For this you should work that our test-suite passes on Python 3. 

# Increase ease-of-use

## Create a command line interface for MDAnalysis tasks

**Difficulty**: Medium

**Mentors**: Max, Richard, Manuel

Currently MDAnalysis exists only as a framework, however making common tasks available via the command line would make certain work flows easier.  As an example, the conversion of trajectories between formats could take the form: 

``` bash
mda convert --topology adk.psf -i adk_dims.dcd -o adk_dims.xtc
```

This project would involve creating a template for these command line utilities to follow and implementing a foolproof user interface for navigating them [using a popular command line parsing library](https://realpython.com/blog/python/comparing-python-command-line-parsing-libraries-argparse-docopt-click/).

# Improve the library core

## Switch from pure ndarray's to unit aware nd-arrays.

**Difficulty**: Hard

**Mentors**: Max, Richard, Manuel

MDAnalysis is using Anström and picoseconds as default units. Our Reader/Writer objects are only aware of units to
the extend that they convert other MD-formats to our default units. But we can also read the coordinates in the native units. This can make it hard to remember what units the coordinates of an AtomGroup have, to fix this you should switch from pure numpy arrays to a unit aware numpy-ndarray wrapper. See [Issue #596](https://github.com/MDAnalysis/mdanalysis/issues/596)

