MDAnalysis uses lots of extensions, below is a list of them and their readiness for Python 3

| Extension                 | Description                  | Language   | Python 3 ready?          | Solved in |
|---------------------------|------------------------------|------------|--------------------------|-----------|
| coordinates._dcd          | DCD format reader            | CPython    | No (Cpython api changed) | |
| coordinates.dcdtimeseries | Timeseries for DCD           | Cpython    | No  (as above)           | |
| lib._distances            | Distance calculations        | C / Cython | Yes?                     | |
| lib._distances_openmp     | Parallel of above            | C / Cython | Yes?                     | |
| ~~lib.parallel.distances~~    | More parallel distances      | Cython     | Yes?? | [#530](https://github.com/MDAnalysis/mdanalysis/issues/530) |
| lib.qcprot                | RMSD via QCP                 | Cython     | Yes?                     | |
| lib._transformations      | Coordinate transformations   | C/CPython  | Probably ([#401](https://github.com/MDAnalysis/mdanalysis/issues/401))  | |
| lib.KDTree._KDTree        | KDTree for distance searches | CPP/SWIG   | No (#383)                | [395](https://github.com/MDAnalysis/mdanalysis/pull/395) |
| coordinates.xdrfile       | Gromacs format reader        | C/Swig     | No ([#211](https://github.com/MDAnalysis/mdanalysis/issues/211))                | [441](https://github.com/MDAnalysis/mdanalysis/pull/441#event-517264871) |