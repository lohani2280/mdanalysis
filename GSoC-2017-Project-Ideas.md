<img src="https://developers.google.com/open-source/gsoc/images/gsoc2016-sun-373x373.png" title="Google Summer of Code 2017" alt="Google Summer of Code 2017" align="right"/>
A list of projects ideas for [[Google Summer of Code 2017|Google-Summer-Of-Code]].

The current proposed projects are:

1. [[Implement efficient parallel analysis of trajectories|GSoC-2017-Project-Ideas#implement-efficient-parallel-analysis-of-trajectories]]
2. [[Improve distance search|GSoC-2017-Project-Ideas#improve-distance-search]]
3. [[Add new MD-Formats|GSoC-2017-Project-Ideas#add-new-md-formats]]

**Or work on your your own idea!** Get in contact with us to propose an idea and we will work with you to flesh it out into a full project. Raise an issue in the [Issue Tracker](/MDAnalysis/mdanalysis/issues) or contact us via the [developer Google group](http://developers.mdanalysis.org).

------

# Implement efficient parallel analysis of trajectories 

**Difficulty**: Hard

**Mentors**: Manuel

Molecular simulation trajectories are very often analyzed frame-by-frame. This is frequently an embarrassingly parallel procedure, in which work can be efficiently divided simply by splitting the trajectory and letting each worker process one of the chunks. The goal of this project is to implement a parallelization framework that automates all the trajectory splitting, work distribution, and eventual result collection.

A parallelization framework should put the least burden possible on the end-user, so that minimal changes are required to turn serial code into parallel. Likewise, the parallelization framework must blend naturally with the analysis API of MDAnalysis. In this way, analyses written using [analysis.base](https://github.com/MDAnalysis/mdanalysis/blob/5b6471d93a36581d06ec73a1a0bddc8a460d4213/package/MDAnalysis/analysis/base.py#L35) will automatically become parallelizable.

Implementing parallelization in Python code can be done in [many ways](https://wiki.python.org/moin/ParallelProcessing). Aspects to consider when choosing one or several approaches are:
- Most users will primarily have access to SMP parallelization;
- Notwithstanding the above point, many users also typically have access to multi-node HPC clusters, and we should be able to leverage their use;
- In an analysis context, being able to write results to shared memory will improve the memory usage footprint and simplify result collection;
- GPU parallelization is attractive for its wide availability (though possibly more complex to implement in a meaningful way).

# Improve distance search 

**Difficulty**: Hard

**Mentors**: Manuel, Richard

Analysis of molecular dynamics simulations typically involves calculations of based upon atoms which are spatially close to each other.  For example a radial distribution function is often only interesting up to distances of around 1.6 nm.
The naive approach to calculate this is to calculate the distance between each pair of atoms, however as the size of the system grows the number of useful pair distances decreases while the computational cost scales as N^2.

To greatly improve the efficiency of this operation, we can first decompose the total simulation volume into smaller cells.  We can then calculate the distances between atom pairs in neighbouring cells. If atoms are not in neighbouring cells we already know that the distance is to large to be interesting. A theoretical description of this algorithm can be found in [this book Appendix F](http://www.amazon.de/Understanding-Molecular-Simulation-Applications-Computational/dp/0122673514%3FSubscriptionId%3DAKIAILSHYYTFIVPWUY6Q%26tag%3Dduckduckgo-ffnt-de-21%26linkCode%3Dxm2%26camp%3D2025%26creative%3D165953%26creativeASIN%3D0122673514)

One domain decomposition algorithm is [cell grids](https://github.com/richardjgowers/cellgrid).

In this project you would integrate the cell grid algorithm into MDAnalysis. 

# Add new MD-Formats

**Dificulty**: Medium

**Mentors**: Richard

One of the strengths of MDAnalysis is its ability to support a wide range of different MD-formats. But we are still missing some like the new [TNG file format](http://onlinelibrary.wiley.com/doi/10.1002/jcc.23495/abstract) from Gromacs , [H5MD](https://github.com/pdebuyl/pyh5md) or the [HALMD](http://halmd.org/) format. Alternatively, you can also add a format that you want to use personally in MDAnalysis.

This project will familiarize you with working with and connecting different APIs,
as well as giving insight into how modern portable data storage file formats work.
It is vitally important that data is read correctly, otherwise analysis will fail at the very first step.
For this reason, there will be a heavy emphasis on the testing for any code written,
and so the project will also teach good practice in software testing.

