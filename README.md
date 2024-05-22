# QAP

A TaiC quadratic assignment problem sover.

Sources:

- `src/cb.cc` kb model
- `src/toff.cc` intuitive model
- `src/custom.cc` kb model with ad-hoc solver
- `instance/*.py` instance generation
- `solution/*.py` data processing
- `plot/heatmap.py` display solutions


## Build

Using CMake:

```
mkdir src/build && cd src/build && cmake -DCMAKE_BUILD_TYPE=Release -DCPLEX_ROOT_DIR=/path/to/cplex/root .. && make
```


## Run

```
./src/build/cb [parameters] < instance/tai??c_?x?_??
./src/build/toff [parameters] < instance/tai??c_?x?_??
./src/build/custom [parameters] < instance/tai??c_?x?_??
```

Paramters details at `src/lib/params.cc`.
