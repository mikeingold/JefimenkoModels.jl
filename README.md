# JefimenkoModels.jl

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://mikeingold.github.io/JefimenkoModels.jl/dev/)

`JefimenkoModels.jl` is a time-domain solver for the electromagnetic near-fields produced by an arbitrary distribution of charges and currents, both electric and magnetic. This solver implements a generalized version of the Jefimenko equations that enables the consideration of magnetic charges and currents, which are often a useful analytical tool in electromagnetics modeling. The solution process operates purely in the time-domain, enabling the study of wideband sources without artifacts caused by frequency-domain analysis and with reduced memory usage compared to FDTD methods.

This package leverages the [UnitfulCoordinateSystems.jl](https://github.com/mikeingold/UnitfulCoordinateSystems.jl) package to provide a handy and performant way to deal with `Unitful` coordinate data.

## Status

The solver methods have been validated using a major commercial software package's Method of Moments (MoM) solver as a benchmark. For a single-frequency (CW) source signal, `JefimenkoModels` produced identical results as the competitor MoM solver. However, when the source signal was defined as a wideband transient pulse, the `JefimenkoModels` solver was substantially faster and more accurate: the MoM solver uses a discrete-frequency-domain transfer function that introduces artifacts/error when solving for wideband signals.
