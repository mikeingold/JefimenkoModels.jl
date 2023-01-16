# TODO

## Short Term
- Rename internal functions with double-underscore
- Remove unicode from function names
- Benchmark the solvers
    - HCubature vs QuadGK for LineSource's
    - Compare 3D volume source performance with JefiGPU

## Medium Term
- Implement E/H(VolumeSource_Cylinder)
- Implement E/H(VolumeSource_Sphere)
- Develop constructor methods for sources and models

## Longer-Term Vision
- Add a CITATION.bib
- Evaluate whether Automatic Differentiation can be made to operate through solutions
- Consider permitting sources to have a variable center/orientation
- Consider consolidating the integrand functions using ComponentArray-parameterized source values
    - This would add complexity to the E/H functions, but would reduce code duplication here
    - The main current difference between R1/R2/R3 is in commented dimensional analysis
    - If Unitful evaluation is a serious performance penalty, then R1/R2/R3 could ustrip source
      values into implied units and then call a consolidated/abstract integrand function
