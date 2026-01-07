# CFD25/26 Project 1

- Author: Jed Nelson
- Email: jn1422@ic.ac.uk
- CID: 02208047
- Github repo: https://github.com/Jedsonofnel/CFD25-project1

## Description

Simple 1D diffusion-convection solver with openFOAM style operators and other
bits.

## Requirements

- Go 1.21+ installed ([see here](go.dev))
- Unix-like system (Linux, macOS or Windows with WSL/Git Bash)

## Building and Running

Compile and generate all data (runs a bunch of predetermined cases):
```bash
make peclet-sweep
make grid-sweep
```

Run a single case:
```bash
make
./bin/cli -n 20 -u 5.0
```

Run a single case and save the data to a .csv file:
```bash
./bin/cli -n 20 -u 5.0 > data.csv
```

See the following for information on CLI parameters:
```bash
./bin/cli --help
```

Clean build artefacts:
```bash
make clean
```

## Output

If you run `make peclet-sweep` or `make grid-sweep` the data is written to the
`data/` directory.
