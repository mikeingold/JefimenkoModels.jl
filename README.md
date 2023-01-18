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

## Tutorial

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

### Define a source

When any of the source components is neglected, e.g. a source with only currents ($J$) or
charges ($\rho$), a pair of pre-defined null sources are provided for convenience. The
`JefimenkoModels.NULL_CHARGE` function can be used in place of either $\rho(\bar{r},t)$
function, and the `JefimenkoModels.NULL_CURRENT` function can be used in place of either
$J(\bar{r},t)$ function.

In the current version of `JefimenkoModels`, source charge and current functions must be
defined in a fairly-specific format. The functions should take two arguments: a
`UnitfulCoordinateSystem.AbstractCoordinate` indicating the spatial position evaluated, and
the `Real`-typed time in implied units of seconds. The functions should return a Real-valued
number with implied units according to the following tables.

An update is planned that will require a `Unitful` time argument and return types. This will
hopefully simplify the source design process and identify potential dimensional errors.

**Table: Line Source Functions**

| Function | Arg 1 | Arg 2 [Units] | Returns [Units] |
|---|---|---|---|
| Electric charge density $\rho_e(\bar{r},t)$ | `r::AbstractCoordinate` | `t::Real` [s] | `rho_e::Real` [C/m] |
| Magnetic charge density $\rho_h(\bar{r},t)$ | `r::AbstractCoordinate` | `t::Real` [s] | `rho_h::Real` [Wb/m] |
| Electric current density $J_e(\bar{r},t)$   | `r::AbstractCoordinate` | `t::Real` [s] | `J_e::SVector{3,Real}` [A] |
| Magnetic current density $J_h(\bar{r},t)$   | `r::AbstractCoordinate` | `t::Real` [s] | `J_h::SVector{3,Real}` [V] |

**Table: Surface Source Functions**

| Function | Arg 1 | Arg 2 [Units] | Returns [Units] |
|---|---|---|---|
| Electric charge density $\rho_e(\bar{r},t)$ | `r::AbstractCoordinate` | `t::Real` [s] | `rho_e::Real` [C/m$^2$] |
| Magnetic charge density $\rho_h(\bar{r},t)$ | `r::AbstractCoordinate` | `t::Real` [s] | `rho_h::Real` [Wb/m$^2$] |
| Electric current density $J_e(\bar{r},t)$   | `r::AbstractCoordinate` | `t::Real` [s] | `J_e::SVector{3,Real}` [A/m] |
| Magnetic current density $J_h(\bar{r},t)$   | `r::AbstractCoordinate` | `t::Real` [s] | `J_h::SVector{3,Real}` [V/m] |

**Table: Volume Source Functions**

| Function | Arg 1 | Arg 2 [Units] | Returns [Units] |
|---|---|---|---|
| Electric charge density $\rho_e(\bar{r},t)$ | `r::AbstractCoordinate` | `t::Real` [s] | `rho_e::Real` [C/m$^3$] |
| Magnetic charge density $\rho_h(\bar{r},t)$ | `r::AbstractCoordinate` | `t::Real` [s] | `rho_h::Real` [Wb/m$^3$] |
| Electric current density $J_e(\bar{r},t)$   | `r::AbstractCoordinate` | `t::Real` [s] | `J_e::SVector{3,Real}` [A$^2$] |
| Magnetic current density $J_h(\bar{r},t)$   | `r::AbstractCoordinate` | `t::Real` [s] | `J_h::SVector{3,Real}` [V$^2$] |

### Construct a model

`JefimenkoModel`s have a `metadata::Dict` provision. This dictionary is not currently used
by the solver. Rather, it provides the user with a convenient place to store any desired
metadata.

The following example produces a `JefimenkoModel` with a single one-meter line source on
the x-axis. This source is stimulated by a spatially-uniform electric current.
```julia
model_cw = let
    # Single line source on x-axis from -0.5m to +0.5m
    # Electric current only: spatially-uniform, x-directed, driven by 100 MHz CW sinusoid
    a = CoordinateCartesian(-0.5u"m", 0.0u"m", 0.0u"m")
    b = CoordinateCartesian( 0.5u"m", 0.0u"m", 0.0u"m")
    Je(r̄::AbstractCoordinate, t_s::Real) = x̂ .* cos(2π*100e6*t_s)     # t in s -> Je in A
    source = LineSource_Straight{Float64}(a, b, NULL_CHARGE, NULL_CHARGE, Je, NULL_CURRENT)

    metadata = Dict(:name => "Tutorial Example",
                    :charges => "None",
                    :currents => "Electric-Only"
                    :spatial_distribution => "Uniform",
                    :source_length => 1.0u"m",
                    :signal_type => "100 MHz CW")

    JefimenkoModel{Float64}(CLASSICAL_VACUUM, [source], metadata)
end
```

### Calculate the electromagnetic fields (TODO)
