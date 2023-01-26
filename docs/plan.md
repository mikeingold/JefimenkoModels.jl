# TODO

## In Progress
- Should __DEFAULT_RTOL be defined as a function __DEFAULT_RTOL(T)?

## Short Term
- Implement custom `Base.show` pretty-printing for structs
- Either remove Unicode from struct field naming OR provide non-Unicode accessors
- Benchmark a 3D volume source against JefiGPU

## Medium Term
- Make better documentation for source function definitions (units, types, etc)
    - Create a test/inspect/validate function for users to look for issues in their definitions?
- Address type stability of source functions by extend struct parameterization
    - [Reference to conversation on Julia Zulip](https://julialang.zulipchat.com/#narrow/stream/225542-helpdesk/topic/.E2.9C.94.20High.20GC.20Time.20in.20HCubature/near/323730178)
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
