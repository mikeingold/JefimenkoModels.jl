# JefimenkoModels.jl

`JefimenkoModels.jl` is a time-domain solver for the electromagnetic near-fields produced by some distribution of source charges and currents.

This solver implements a generalized variant of the Jefimenko equations that allows for the consideration of magnetic charges and currents, which are often a useful analytical tool in electromagnetics modeling.

## Status

This package remains in development status. Not all planned solver methods are yet implemented, and results have not yet been validated.

Public Function | Implemented | Tested | Validated
:--|-:-|-:-|-:-
`E(LineSource_General)` | :white_check_mark: | :x: | :x:
`H(LineSource_General)` | :white_check_mark: | :x: | :x:
`E(SurfaceSource_DiskGeneral)` | :white_check_mark: | :x: | :x:
`H(SurfaceSource_DiskGeneral)` | :white_check_mark: | :x: | :x:
