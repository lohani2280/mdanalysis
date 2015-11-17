For installing the **latest versions** of MDAnalysis see the [[Install]] notes.

Older versions of MDAnalysis had different requirements and for historical reasons they are listed here:

## up to 0.7.2 ##

A LAPACK library is required (standard [LAPACK/BLAS](http://www.netlib.org/lapack), [ATLAS](http://www.netlib.org/atlas/), or the [Intel Math Kernel Library](http://software.intel.com/en-us/intel-mkl/) (MKL); on Mac OS X we simply use the native fast [VecLib framework](http://developer.apple.com/hardwaredrivers/ve/vector_libraries.html) so you don't have to install anything else).

### Fast math libraries ###

On Mac OS X we use the systems vecLib for fast linear algebra and
nothing needs to be configured. On Linux you will probably want to use
something such as ATLAS or the Intel Math Kernel Libraries for better
performance for the rms fitting.

The fast math library paths are set in `setup.cfg`. By default we ship a template file `setup.cfg.template` that you can
```
cp setup.cfg.template setup.cfg
```
and edit. See the file itself for examples and the ones below.

#### Example MKL ####

Linux section of `setup.cfg` (your paths are probably different from the
ones in this example)::

```
  [linux]
  fast_numeric_include = /opt/intel/cmkl/10.0.5.025/include
  fast_numeric_linkpath = /opt/intel/cmkl/10.0.5.025/lib/em64t
  fast_numeric_libs = mkl_lapack mkl guide
```


#### Example ATLAS LAPACK ####

If you want to use the ATLAS LAPACK libraries use something such as ::

```
[linux]
fast_numeric_include = /usr/include
fast_numeric_linkpath = /usr/lib/atlas
fast_numeric_libs = lapack
```

Your paths could be different. Use a command such as
```
locate lapack
```
or your package manager to get an idea of where the files are to be found.

On some versions of Ubuntu we had to use
```
[linux]
# Ubuntu 9.04 i686
fast_numeric_include = /usr/include
fast_numeric_linkpath = /usr/lib/sse2/atlas
fast_numeric_libs = lapack
```

### Adding lapack packages ###
```
 sudo apt-get install liblapack-dev   # see below for notes on fast_numeric_*
```
