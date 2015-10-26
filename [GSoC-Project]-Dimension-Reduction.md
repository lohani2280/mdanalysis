MDAnalysis already comes with a range of different standard analysis tools but currently lacks an implementation of a general dimension reduction algorithm, that can select an arbitrary number of dimensions of interest. 3 common general techniques are

- [Principle Component Analysis](https://en.wikipedia.org/wiki/Principal_component_analysis)
- [Time Independent Component Analysis](http://arxiv.org/abs/1302.6614)
- [Diffusion Maps] (http://arxiv.org/abs/1506.06259)

There are python implementations for all of these algorithms but none of them currently work with MDAnalysis out of the box. This is because the current python impementatoins work on normal numpy arrays that stores a complete trajectory in memory, but MDAnalysis never loads the whole trajectory but only one frame at a time. This approach allows MDAnalysis to treat very large system on a normal laptop or workstation.

Of course you can also suggest us another dimension reduction algorithm that you would like to implement.