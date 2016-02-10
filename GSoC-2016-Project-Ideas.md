<img src="https://developers.google.com/open-source/gsoc/images/gsoc2016-sun-373x373.png" title="Google Summer of Code 2016" alt="Google Summer of Code 2016" align="right"/>
A list of projects ideas for [[Google Summer of Code 2016|Google-Summer-Of-Code]]. Each with a short sentence describing it. Please also note that not every project we suggest here will take you all summer. Some will be short (eg TNG format) and with others you might only be able to implement parts (eg [[dask backend|#Experiment-with-different-backends-for-the-trajectory-classes]]). 

The project ideas can be roughly categorized as

1. [[New analysis functionality|#new-analysis-functionality]]
2. [[Increasing performance|#increasing-performance]]
3. [[New input formats|#new-input-formats]]
4. [[Increase platform availability|#increase-platform-availability]]
5. [[Increase ease-of-use|#Increase-ease-of-use]]
6. [[Improve the library core|#Improve the library core]]

------

# New analysis functionality
## Implement a general dimension reduction algorithm

MDAnalysis already comes with a range of different standard analysis tools but currently lacks an implementation of a general dimension reduction algorithm, that can select an arbitrary number of dimensions of interest. 3 common general techniques are

- [Time Independent Component Analysis](http://arxiv.org/abs/1302.6614)
- [Diffusion Maps] (http://arxiv.org/abs/1506.06259)
- [Principle Component Analysis](https://en.wikipedia.org/wiki/Principal_component_analysis)

There are python implementations for all of these algorithms but none of them currently work with MDAnalysis out of the box. This is because the current python implementations work on normal numpy arrays that stores a complete trajectory in memory, but MDAnalysis never loads the whole trajectory but only one frame at a time. This approach allows MDAnalysis to treat very large system on a normal laptop or workstation.

Of course you can also suggest us another dimension reduction algorithm that you would like to implement.


## Implement flexible search of atoms inside volume elements defined by densities. 

Combine spatial densities (e.g. from time averaged quantities or experimental data such as electron densities) with atom-based queries in order to aid multiscaling approaches and comparisons between experiment and simulation.

# Increasing performance

## Experiment with different backends for the trajectory classes

Most MD-simulations produce way more data that we could fit into the RAM, even with a modern computer.
To cope with this, MDAnalysis never loads a full trajectory but only one frame at a time. This comes with a performance penalty. There are new Python packages like [Dask](http://dask.pydata.org/en/latest/) and  [Blaze](http://blaze.pydata.org/) that can potentially help us here. You should look into the different distributed computation numerical array libraries in python and implement a reader using it during the summer.


## Develop a analysis pipeline framework for multi-core CPUs

Help us implement a general pipeline to use multiple CPU-cores for analysis tasks ([Dask](http://dask.pydata.org/en/latest/), [MPI](http://pythonhosted.org/mpi4py/usrman/index.html), or even a hybrid approach)
For an example of the direction this is currently taking, [see here](https://github.com/MDAnalysis/mdanalysis/pull/618)

## Improve distance search 

Work with domain-decomposition algorithms to improve our distance search algorithms ([cell grids](https://github.com/richardjgowers/cellgrid)) and/or implement distance search on **GPUs** using CUDA/OpenCL.

# New input formats

## Add new MD-Formats

One of the strengths of MDAnalysis is its ability to support a wide range of different MD-formats. But we are still missing some like the new [TNG file format](http://onlinelibrary.wiley.com/doi/10.1002/jcc.23495/abstract) from Gromacs or [H5MD](https://github.com/pdebuyl/pyh5md). Alternatively, you can also add a format that you want to use personally in MDAnalysis.
This project will familiarise you with working with and connecting different APIs,
as well as giving insight into how modern portable data storage file formats work.

## Random Walk Trajectory Backend

To check if a new analysis-method works as intended it is often a good idea to use it with a random walk in different simple energy landscapes (A flat energy, harmonic well, double well). In this project you would develop a 'Reader' that produces random trajectories. 

Langevin dynamics in a energy landscape are close to the conformational dynamics of proteins, see [1]. As a first
start you could implement a integrator for langevin dynamics and later have the trajectory 'reader' use the integrator to dynamically generate the trajectory. 

[1] Robert Zwanzig. Nonequilibrium statistical mechanics. Oxford University Press,
2001

# Increase platform availability

## Help port MDAnalysis to Python 3

Python 3 is getting adopted by a wider range of users and unix distributions are starting to switch.
MDAnalysis can't run right now under Python 3 mostly due to it's [C/Cython extensions](https://github.com/MDAnalysis/mdanalysis/wiki/List-of-extensions), we currently try to move our C-extensions to cython which supports Python 2 and 3 with one source. See also [#260](https://github.com/MDAnalysis/mdanalysis/issues/260)

## Port MDAnalysis to Windows. 

None of the current devs has a Windows environment. But some research groups do use Windows and it would be nice they could use MDAnalysis as well. Since neither of us has experience with python extensions on windows we don't know what exactly is needed to make this happen.

# Increase ease-of-use

## Create a command line interface for MDAnalysis tasks

Currently MDAnalysis exists only as a framework, however making common tasks available via the command line would make certain work flows easier.  As an example, the conversion of trajectories between formats could take the form: 

``` bash
mda convert --topology adk.psf -i adk_dims.dcd -o adk_dims.xtc
```

This project would involve creating a template for these command line utilities to follow and implementing a foolproof user interface for navigating them [using a popular command line parsing library](https://realpython.com/blog/python/comparing-python-command-line-parsing-libraries-argparse-docopt-click/).

# Improve the library core

## Formalize our Atom selection parser

Implement a formal flexible parser for atom selections (using [pyparsing](https://pyparsing.wikispaces.com/); see also [our discussion](https://github.com/MDAnalysis/mdanalysis/issues/371) on it).

## Switch from pure ndarray's to unit aware nd-arrays.

MDAnalysis is using Anström and picoseconds as default units. Our Reader/Writer objects are only aware of units to
the extend that they convert other MD-formats to our default units. But we can also read the coordinates in the native units. This can make it hard to remember what units the coordinates of an AtomGroup have, to fix this you should switch from pure numpy arrays to a unit aware numpy-ndarray wrapper. See [Issue #596](https://github.com/MDAnalysis/mdanalysis/issues/596)

## Or your idea here! Get in contact with us to propose an idea.
Raise an issue in the [Issue Tracker](/MDAnalysis/mdanalysis/issues) or contact us via the [developer Google group](http://developers.mdanalysis.org).