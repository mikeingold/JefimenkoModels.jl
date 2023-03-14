# JefimenkoModels.jl

`JefimenkoModels.jl` is a time-domain solver for the electromagnetic near-fields produced by
an arbitrary distribution of charges and currents, both electric and magnetic. This solver
implements a generalized version of the Jefimenko equations that enables the consideration of
magnetic charges and currents, which are often a useful analytical tool in electromagnetics
modeling. The solution process operates purely in the time-domain, enabling the study of
wideband sources without artifacts caused by frequency-domain analysis and with reduced
memory usage compared to FDTD methods.

This package leverages the
[UnitfulCoordinateSystems.jl](https://gitlab.com/mike.ingold/unitfulcoordinatesystems.jl)
package to provide a handy and performant way to deal with `Unitful` coordinate data.

## Status

This package remains in development status. Multiple dispatch is used to select the solver
method appropriate for a particular source type. The implementation status of these methods
is detailed in the following table.

| Solver Method | Implemented | Tested |
|:---|:---:|:---:|
| `LineSource_Straight`       | :white_check_mark: | :white_check_mark: |
| `SurfaceSource_Disk`        | :white_check_mark: | :white_check_mark: |
| `SurfaceSource_Rectangle`   | :white_check_mark: | :white_check_mark: |
| `VolumeSource_Rectangular`  | :white_check_mark: | :white_check_mark: |
| `VolumeSource_Cylinder`     |         :x:        |         :x:        |
| `VolumeSource_Sphere`       |         :x:        |         :x:        |

The `LineSource_Straight` solver methods have been validated against a major commercial
software package's Method of Moments (MoM) solver for electric current line sources. For a
single-frequency (CW) source signal, `JefimenkoModels` produced identical results as the
competitor MoM solver. However, when the source signal was defined as a wideband transient
pulse, the `JefimenkoModels` solver was substantially faster and more accurate: the MoM
solver uses a discrete-frequency-domain transfer function that introduces artifacts/error
when solving for wideband signals.
