For some platforms, MDAnalysis can be installed directly in pre-compiled form. (If no pre-compiled packages are available, please see [installing MDAnalysis from source](Install).)

# Anaconda
**Note**: This currently **experimental**, feedback welcome (at PR [#751](/MDAnalysis/mdanalysis/issues/751))
* *Experimental, since 0.14.0* (see issue [#608](/MDAnalysis/mdanalysis/issues/608) with PR [#751](/MDAnalysis/mdanalysis/issues/751)).
* Currently only Linux 64 bit supported.
* Currently only Python 2.

Install the current anaconda package according to the [conda docs](https://docs.continuum.io/anaconda/install). Users can either install by one step:

```bash
conda install -c https://conda.anaconda.org/mdanalysis mdanalysis
```

Or two steps, first add the official MDAnalysis anaconda channel, (this has to be done once):

```bash
conda config --add channels MDAnalysis
```

Then install MDAnalysis
```bash
conda install mdanalysis
```
Update MDAnalysis
```bash
conda update mdanalysis
```


# Distributions
Some [third-party repositories](Downloads#binary-distributions) might have binary packages of MDAnalysis but we are not currently making any distribution-specific binary packages available. 

(If you are package maintainer for a distribution and you want your distribution listed here, please [/MDAnalysis/mdanalysis/issues](raise an issue) to notify us and provide information and links to the distribution and the MDAnalysis package itself.)


