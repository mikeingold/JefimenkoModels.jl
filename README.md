# JefimenkoModels.jl

`JefimenkoModels.jl` is a time-domain solver for the electromagnetic near-fields produced by some distribution of source charges and currents.

This solver implements a generalized variant of the Jefimenko equations that allows for the consideration of magnetic charges and currents, which are often a useful analytical tool in electromagnetics modeling.

## Status

This package remains in development status. Not all planned solver methods are yet implemented, and results have not yet been validated.

| Public Function | Implemented | Tested | Validated |
|:---|:---:|:---:|:---:|
| `E(LineSource_Straight)`       | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `E(SurfaceSource_Disk)`        | :white_check_mark: | :white_check_mark: | :x: |
| `E(SurfaceSource_Rectangle)`   | :x: | :x: | :x: |
| `E(VolumeSource_Rectangular)`  | :x: | :x: | :x: |
| `E(VolumeSource_Cylinder)`     | :x: | :x: | :x: |
| `E(VolumeSource_Sphere)`       | :x: | :x: | :x: |
| `H(LineSource_Straight)`       | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `H(SurfaceSource_Disk)`        | :white_check_mark: | :white_check_mark: | :x: |
| `H(SurfaceSource_Rectangle)`   | :x: | :x: | :x: |
| `H(VolumeSource_Rectangular)`  | :x: | :x: | :x: |
| `H(VolumeSource_Cylinder)`     | :x: | :x: | :x: |
| `H(VolumeSource_Sphere)`       | :x: | :x: | :x: |

## Usage (TODO)

- How to construct a model
    - `LineSource_Straight`
    - `SurfaceSource_Disk`
- How to calculate the fields
    - `E(r,t,model)`
    - `H(r,t,model)`
    - `P(r,t,model)`
