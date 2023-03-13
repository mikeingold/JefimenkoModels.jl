# JefimenkoModels.jl

This package is a time-domain solver for the electromagnetic near-fields produced by
an arbitrary distribution of charges and currents, both electric and magnetic. This solver
implements a generalized version of the Jefimenko equations that enables the consideration of
magnetic charges and currents, which are often a useful analytical tool in electromagnetics
modeling. The solution process operates purely in the time-domain, enabling the study of
wideband sources without artifacts caused by frequency-domain analysis and with reduced
memory usage compared to FDTD methods.

## Status

This package remains in development status.  Multiple dispatch is used to select the solver
method appropriate for a particular source type, but not all source types have implemented
solver methods yet.

