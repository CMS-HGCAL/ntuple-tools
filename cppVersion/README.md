# Ntuple tools C++

The C++ version of the ntuple-tools containing ImagingAlgo class running 2D clasterization over ntuples stored in the ROOT files. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

* ROOT (the tool was tested with v6.12.06_2)

### Building

To compile siply run:
```
make
```

## Running

Make will produce an executable:

* algoBenchmark

Configuration if at the moment hardcoded in the *algoBenchmark.cpp*.  Set the `inputPath`, `outDir`, range of ntuples to run over (`minNtuple` and `maxNtuple`). Adjust algorithm settings if needed. Then simply run:
```
./algoBenchmark
```

## Contributing

The project is basically split in two parts:
* core, which is in *src/* and *include/* directories,
* user code, directly in the top directory (e.g. *algoBenchmark*).

If you wish to add another test of the core functionalities:
* add a new file in the main directory,
* add you program to the Makefile, copying parts regarding *algoBenchmark* and changing dependencies accordingly.

If you change the core, you hopefully know very well what to do.
