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

The `LineSource_Straight` solver methods have been validated against a major commercial
software package's Method of Moments (MoM) solver for electric current line sources. For a
single-frequency (CW) source signal, `JefimenkoModels` produced identical results as the
competitor MoM solver. However, when the source signal was defined as a wideband transient
pulse, the `JefimenkoModels` solver was substantially faster and more accurate: the MoM
solver uses a discrete-frequency-domain transfer function that introduces artifacts/error
when solving for wideband signals.

## Usage

### Define a propagation media

The propagation media of a model is assumed to be linear, time-invariant, and
spatially-homogeneous.

When designing a model that will propagate in vacuum (free space), a pre-defined media is
provided.
```julia
media::PropagationMedia_Simple = JefimenkoModels.CLASSICAL_VACUUM
```

Alternatively, a propagation media with Real-valued permittivity ($\varepsilon$) and
permeability ($\mu$) can be specified using `Unitful` values. Each term can be defined in
your preferred choice of units, so long as they are dimensionally-equivalent to the reference
units: $\varepsilon$ in [F/m] or [As/Vm], and $\mu$ in [N/A$^2$] or [Vs/Am].
```julia
epsilon = 8.854_188e-15 * (u"kA"*u"s")/(u"V"*u"m")
mu = 1.256_637e-3  * (u"mV"*u"s")/(u"A"*u"m")
c = 2.997_925e5 * u"km"/u"s"
PropagationMedia_Simple(epsilon, mu, c)
```

### Define a source (TODO)

### Construct a model (TODO)

### Calculate the electromagnetic fields (TODO)