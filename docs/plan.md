# TODO

## High Priority
- Implement H-field integrand functions
- Implement H-field calculations for LineSource and SurfaceSource types
- Compare results of E-field calculations with prior animation script
- Validate results of calculations against FEKO

## Medium Term
- Determine if performance benefits justify specialized integral functions for non-general sources
    - Consider either adding functions to generalize them, or just removing them totally
- Fully populate the README Status with all struct types and implementation status

## Longer-Term Vision
- Add a CITATION.bib
- Evaluate whether Automatic Differentiation can be made to function through Integral solvers
- Consider consolidating the integrand functions using ComponentArray-parameterized source values
    - This would add complexity to the E/H functions, but would reduce code duplication here
    - The main current difference between R1/R2/R3 is in commented dimensional analysis
    - If Unitful evaluation is a serious performance penalty, then R1/R2/R3 could ustrip source
      values into implied units and then call a consolidated/abstract integrand function

## Performance Improvement
- Re-benchmark all of the solver code
- Compare 3D volume source performance with JefiGPU
- Determine why H function is taking 3x as long to compute as E for same arguments
    - Also requires about 3x the memory/allocations, probably something to do with integral solver
