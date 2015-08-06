MDAnalysis uses lots of extensions, below is a list of them and their readiness for Python 3

| Extension                 | Description                  | Language   | Python 3 ready?          |
|---------------------------|------------------------------|------------|--------------------------|
| coordinates._dcd          | DCD format reader            | CPython    | No (Cpython api changed) |
| coordinates.dcdtimeseries | Timeseries for DCD           | Cpython    | No  (as above)           |
| lib._distances            | Distance calculations        | C / Cython | Yes?                     |
| lib._distances_openmp     | Parallel of above            | C / Cython | Yes?                     |
| lib.parallel.distances    | More parallel distances      | Cython     | Yes??                    |
| lib.qcprot                | RMSD via QCP                 | Cython     | Yes?                     |
| lib._transformations      | Coordinate transformations   | C/CPython  | No (Rewrite interface)   |
| lib.KDTree._KDTree        | KDTree for distance searches | CPP/SWIG   | No                       |
| coordinates.xdrfile       | Gromacs format reader        | C/Swig     | ???                      |