# TODO

## In Progress
- Condensing integrand functions into single for each (E vs H), reducing type instability
    - Combine both into single file when done.
    - Document expected units and dimensional analysis somewhere
- Should __DEFAULT_RTOL be defined as a function __DEFAULT_RTOL(T)?

## Short Term
- Figure out how to get E.(r,t,model) broadcasting working
- Either remove Unicode from struct field naming OR provide non-Unicode accessors
- Benchmark a 3D volume source against JefiGPU

## Medium Term
- Develop/document constructor methods for sources and models
- Implement solvers for
    - VolumeSource_Cylinder
    - VolumeSource_Sphere

## Longer-Term Vision
- Add a CITATION.bib
- Re-assess the need for solver type parameterization
    - Does it even work as intended?
    - Is there a performance benefit?
- Evaluate whether Automatic Differentiation can be made to operate through solutions
- Consider permitting sources to have a variable center/orientation
- Consider consolidating the integrand functions using ComponentArray-parameterized source values
    - This would add complexity to the E/H functions, but would reduce code duplication here
    - The main current difference between R1/R2/R3 is in commented dimensional analysis
    - If Unitful evaluation is a serious performance penalty, then R1/R2/R3 could ustrip source
      values into implied units and then call a consolidated/abstract integrand function
