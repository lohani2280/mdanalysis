<img src="https://developers.google.com/open-source/gsoc/images/gsoc2016-sun-373x373.png" title="Google Summer of Code 2017" alt="Google Summer of Code 2017" align="right"/>
A list of projects ideas for [[Google Summer of Code 2017|Google-Summer-Of-Code]].

The project ideas can be roughly categorized as

2. [[Increasing performance|GSoC-2017-Project-Ideas#increasing-performance]]
3. [[New input formats|GSoC-2017-Project-Ideas#new-input-formats]]
6. [[Improve the library core|GSoC-2017-Project-Ideas#improve-the-library-core]]

**Or work on your your own idea!** Get in contact with us to propose an idea and we will work with you to flesh it out into a full project. Raise an issue in the [Issue Tracker](/MDAnalysis/mdanalysis/issues) or contact us via the [developer Google group](http://developers.mdanalysis.org).

------

# Increasing performance

## Implement efficient parallel analysis of trajectories 

**Difficulty**: Hard

**Mentors**: 

Molecular simulation trajectories are very often analyzed frame-by-frame. This is frequently an embarrassingly parallel procedure, in which work can be efficiently divided simply by splitting the trajectory and letting each worker process one of the chunks. The goal of this project is to implement a parallelization framework that automates all the trajectory splitting, work distribution, and eventual result collection.

A parallelization framework should put the least burden possible on the end-user, so that minimal changes are required to turn serial code into parallel. Likewise, the parallelization framework must blend naturally with the analysis API of MDAnalysis. In this way, analyses written using [analysis.base](https://github.com/MDAnalysis/mdanalysis/blob/5b6471d93a36581d06ec73a1a0bddc8a460d4213/package/MDAnalysis/analysis/base.py#L35) will automatically become parallelizable.

Implementing parallelization in Python code can be done in [many ways](https://wiki.python.org/moin/ParallelProcessing). Aspects to consider when choosing one or several approaches are:
- Most users will primarily have access to SMP parallelization;
- Notwithstanding the above point, many users also typically have access to multi-node HPC clusters, and we should be able to leverage their use;
- In an analysis context, being able to write results to shared memory will improve the memory usage footprint and simplify result collection;
- GPU parallelization is attractive for its wide availability (though possibly more complex to implement in a meaningful way).

## Improve distance search 

**Difficulty**: Hard

**Mentors**: 

To analyze molecular simulations it is often helpful which atoms are close to each other. For this we calculate distance matrices where the distances between every atom pair is calculated. This is a very expensive operation that grows quadratic with the number of atoms involved.

Since we are only interested in atoms that are close to each other we can use some algorithms run faster
after some initial analysis of the coordinates. One class of these algorithms are domain-decomposition algorithms. The basic idea of this type of algorithms is to decompose the volume occupied by the atoms into different cells and then only calculate distances for atoms in neighboring cells. If atoms are not in neighboring cells we already know that the distance is to big for us to be interesting. A theoretical description of these algorithm can be found in [this book Appendix F](http://www.amazon.de/Understanding-Molecular-Simulation-Applications-Computational/dp/0122673514%3FSubscriptionId%3DAKIAILSHYYTFIVPWUY6Q%26tag%3Dduckduckgo-ffnt-de-21%26linkCode%3Dxm2%26camp%3D2025%26creative%3D165953%26creativeASIN%3D0122673514)

One domain decomposition algorithm is [cell grids](https://github.com/richardjgowers/cellgrid).

In this project you would integrate the cell grid algorithm into MDAnalysis. 

# New input formats

## Add new MD-Formats

**Dificulty**: Medium

**Mentors**: 

One of the strengths of MDAnalysis is its ability to support a wide range of different MD-formats. But we are still missing some like the new [TNG file format](http://onlinelibrary.wiley.com/doi/10.1002/jcc.23495/abstract) from Gromacs , [H5MD](https://github.com/pdebuyl/pyh5md) or the [HALMD](http://halmd.org/) format. Alternatively, you can also add a format that you want to use personally in MDAnalysis.
This project will familiarize you with working with and connecting different APIs,
as well as giving insight into how modern portable data storage file formats work.


# Increase platform availability

## Help port MDAnalysis to Python 3

**Difficulty**: Easy

**Mentors**: 

Python 3 is getting adopted by a wider range of users and unix distributions are starting to switch.
MDAnalysis can't run right now under Python 3 mostly due to it's [C/Cython extensions](https://github.com/MDAnalysis/mdanalysis/wiki/List-of-extensions), we currently try to move our C-extensions to cython which supports Python 2 and 3 with one source. See also [#260](https://github.com/MDAnalysis/mdanalysis/issues/260).

Missing here right now is the DCD trajectory readers. There exists an [incomplete work](https://github.com/MDAnalysis/mdanalysis/pull/682) to enable Python 2/3 of the DCD reader. In this project you would finish this work by either writing finishing this work or by rewriting the DCD interface in cython.

The second part of this project is to remove all other incompatibilities with Python 3 we currently have. For this you should work that our test-suite passes on Python 3. 


