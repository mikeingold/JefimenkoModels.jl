# TODO

## Short Term
- Implement custom `Base.show` pretty-printing for structs
- Add a `runtests.jl` suite
- Register package in General

## Medium Term
- Update `CITATION.bib` when thesis is published
- Make better documentation for source function definitions (units, types, etc)
    - Create a test/inspect/validate function for users to look for issues in their definitions?
- Address type stability of source functions by extend struct parameterization
    - [Reference to conversation on Julia Zulip](https://julialang.zulipchat.com/#narrow/stream/225542-helpdesk/topic/.E2.9C.94.20High.20GC.20Time.20in.20HCubature/near/323730178)
- Document structs and constructor methods
- Implement solvers for
    - VolumeSource_Cylinder
    - VolumeSource_Sphere

## Longer-Term Vision
- Re-assess the need for solver type parameterization
    - Does it even work as intended? Is there a performance benefit?
    - Should __DEFAULT_RTOL be defined as a function __DEFAULT_RTOL(T)?
- Implement sources not centered on the origin
    - Probably use a coordinate system transform to map observation point onto source-centric coordinates
- Consider `Unitful` integral functions
