# POLDIS Source

This directory contains the curated POLDIS source snapshot used for the DIS
comparison studies in this repository.

Requirements:

- a Fortran compiler, defaulting to `gfortran`
- LHAPDF with `lhapdf-config` available in `PATH`

The checked-in helper wrappers build the three supported executables in place
and write `poldis.x` in this directory:

```bash
./compile_dijet
./compile_dijet_rivetplots
./compile_singlejet
```

Useful notes:

- set `FC=ifort` (or another compiler) if you do not want to use `gfortran`
- extra compiler flags can be passed directly, for example
  `./compile_dijet_rivetplots -O2`
- the canonical DIS reference driver for the Rivet comparison is
  `user_dijet_rivetplots.f`
