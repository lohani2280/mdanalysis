<img src="https://developers.google.com/open-source/gsoc/images/gsoc2015-300x270.jpg" title="Google Summer of Code 2016" alt="Google Summer of Code 2016" align="right"/>
A list of projects ideas. Each with a short sentence describing it. Please also note that not every project we suggest here will take you all summer. Some will be short (eg TNG format) and with others you might only be able to implement parts (eg dask backend). 


# Implement a general dimension reduction algorithm

MDAnalysis already comes with a range of different standard analysis tools but currently lacks an implementation of a general dimension reduction algorithm, that can select an arbitrary number of dimensions of interest. 3 common general techniques are

- [Principle Component Analysis](https://en.wikipedia.org/wiki/Principal_component_analysis)
- [Time Independent Component Analysis](http://arxiv.org/abs/1302.6614)
- [Diffusion Maps] (http://arxiv.org/abs/1506.06259)

There are python implementations for all of these algorithms but none of them currently work with MDAnalysis out of the box. This is because the current python impementatoins work on normal numpy arrays that stores a complete trajectory in memory, but MDAnalysis never loads the whole trajectory but only one frame at a time. This approach allows MDAnalysis to treat very large system on a normal laptop or workstation.

Of course you can also suggest us another dimension reduction algorithm that you would like to implement.

# Add new MD-Formats

MDAnalysis already supports a wide range of different MD-formats. But we are still missing some like the new [TNG file format](http://onlinelibrary.wiley.com/doi/10.1002/jcc.23495/abstract) from Gromacs. You can also a format that you want to use personally in MDAnalysis.

# Random Walk Trajectory Backend

To check if a new analysis-method works as intended it is often a good idea to use it with a random walk in different simple energy landscapes (A flat energy, harmonic well, double well). In this project you would develop a 'Reader' that produces random trajectories. 

# Increase the availability of this code by helping port MDAnalysis to Python 3

Python 3 is getting adopted by a wider range of users and unix distributions are starting to switch.
MDAnalysis can't run right now under Python 3 mostly due to it's [C/Cython extensions](https://github.com/MDAnalysis/mdanalysis/wiki/List-of-extensions), we currently try to move our C-extensions to cython which supports Python 2 and 3 with one source. See also [#260](https://github.com/MDAnalysis/mdanalysis/issues/260)

# Port MDAnalysis to Windows. 

None of the current devs has a Windows environment. But some research groups do use Windows and it would be nice they could use MDAnalysis as well. Since neither of us has experience with python extensions on windows we don't know what exactly is needed to make this happen.

# Experiment with different backends for the trajectory classes

Most MD-simulations produce way more data that we could fit into the RAM, even with a modern computer.
To cope with this, MDAnalysis never loads a full trajectory but only one frame at a time. This comes with a performance penalty. There are new Python packages like [Dask](http://dask.pydata.org/en/latest/) and  [Blaze](http://blaze.pydata.org/) that can potentially help us here. You should look into the different distributed computation numerical array libraries in python and implement a reader using it during the summer.

# Develop a analysis pipeline framework for multi-core CPUs

Help us implement a general pipeline to use multiple CPU-cores for analysis tasks ([Dask](http://dask.pydata.org/en/latest/), [MPI](http://pythonhosted.org/mpi4py/usrman/index.html), or even a hybrid approach)

# Improve distance search 

Work with domain-decomposition algorithms to improve our distance search algorithms ([cell grids](https://github.com/richardjgowers/cellgrid)) and/or implement distance search on GPUs.

# Formalize our Atom selection parser

Implement a formal flexible parser for atom selections (using [pyparsing](https://pyparsing.wikispaces.com/); see also [our discussion](https://github.com/MDAnalysis/mdanalysis/issues/371) on it).

# Implement flexible search of atoms inside volume elements defined by densities. 

Combine spatial densities (e.g. from time averaged quantities or experimental data such as electron densities) with atom-based queries in order to aid multiscaling approaches and comparisons between experiment and simulation.

# Or your idea here! Get in contact with us to propose an idea.
Raise an issue in the [Issue Tracker](/MDAnalysis/mdanalysis/issues) or contact us via the [developer Google group](http://developers.mdanalysis.org).