# JefimenkoModels.jl

This package is a time-domain solver for the electromagnetic near-fields produced by
an arbitrary distribution of charges and currents, both electric and magnetic. This solver
implements a generalized version of the Jefimenko equations that enables the consideration of
magnetic charges and currents, which are often a useful analytical tool in electromagnetics
modeling. The solution process operates purely in the time-domain, enabling the study of
wideband sources without artifacts caused by frequency-domain analysis and with reduced
memory usage compared to FDTD methods.

## Dependencies

### UnitfulCoordinateSystems.jl

This package leverages the
[UnitfulCoordinateSystems.jl](https://gitlab.com/mike.ingold/unitfulcoordinatesystems.jl)
package as a handy and performant way of dealing with `Unitful` coordinate data. However,
`UnitfulCoordinateSystems.jl` is not yet registered in Julia's package manager.

This
dependency can be installed in a Julia environment via the `Pkg` REPL-mode, as demonstrated
with the following

```julia
julia> ]
(@v1.X) pkg> add https://github.com/mikeingold/UnitfulCoordinateSystems.jl.git
```

or programmatically via

```julia
import Pkg
Pkg.add(url="https://github.com/mikeingold/UnitfulCoordinateSystems.jl.git")
```
