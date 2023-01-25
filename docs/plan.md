# TODO

## In Progress
-

## Short Term
- Benchmark a 3D volume source against JefiGPU

## Medium Term
- Develop constructor methods for sources and models
- Evaluate type stability of coordinates in integrand functions
- Consider adding default to argument media=CLASSICAL_VACUUM
- Implement solvers for
    - VolumeSource_Cylinder
    - VolumeSource_Sphere

## Longer-Term Vision
- Add a CITATION.bib
- Re-assess the need for solver type parameterization
    - Does it even work as intended?
    - Is there a performance benefit?
    - Should __DEFAULT_RTOL be defined as a function __DEFAULT_RTOL(T)?
- Evaluate whether Automatic Differentiation can be made to operate through solutions
- Consider permitting sources to have a variable center/orientation
- Consider consolidating the integrand functions using ComponentArray-parameterized source values
    - This would add complexity to the E/H functions, but would reduce code duplication here
    - The main current difference between R1/R2/R3 is in commented dimensional analysis
    - If Unitful evaluation is a serious performance penalty, then R1/R2/R3 could ustrip source
      values into implied units and then call a consolidated/abstract integrand function
