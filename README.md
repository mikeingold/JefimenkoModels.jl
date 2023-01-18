# JefimenkoModels.jl

`JefimenkoModels.jl` is a time-domain solver for the electromagnetic near-fields produced by
some distribution of source charges and currents.

This solver implements a generalized variant of the Jefimenko equations that allows for the
consideration of magnetic charges and currents, which are often a useful analytical tool in
electromagnetics modeling.

## Status

This package remains in development status. Multiple dispatch is used to select the solver
method appropriate for a particular source type. The implementation status of these methods
is detailed in the following table.

| Solver Method | Implemented | Tested | Validated |
|:---|:---:|:---:|:---:|
| `LineSource_Straight`       | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `SurfaceSource_Disk`        | :white_check_mark: | :white_check_mark: | :x: |
| `SurfaceSource_Rectangle`   | :white_check_mark: | :x: | :x: |
| `VolumeSource_Rectangular`  | :white_check_mark: | :x: | :x: |
| `VolumeSource_Cylinder`     | :x: | :x: | :x: |
| `VolumeSource_Sphere`       | :x: | :x: | :x: |

## Usage (TODO)

- How to define a media
- How to define a source signal
- How to construct a model
- How to calculate the fields