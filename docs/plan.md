# TODO

## Near Term
- Implement `Base.show` for structs
- Document structs and constructor methods
- Extend the `runtests.jl` suite
- Update `CITATION.bib` when thesis is published
- Implement solvers for 
    - VolumeSource_Cylinder
    - VolumeSource_Sphere

## Longer-Term / Vision
- Improve documentation for source function definitions (units, types, etc)
    - Create a test/inspect/validate function for users to look for issues in their definitions?
- Address possible type stability and parameterization issues
    - Source `::Function`s may be overly broad ([ref conversation on Julia Zulip](https://julialang.zulipchat.com/#narrow/stream/225542-helpdesk/topic/.E2.9C.94.20High.20GC.20Time.20in.20HCubature/near/323730178))
    - Consider applying DynamicQuantities.jl in lieu of Unitful
- Re-assess the need for solver numerical type parameterization
    - Does it even work as intended? Is there a performance benefit?
    - Should __DEFAULT_RTOL be defined as a function __DEFAULT_RTOL(T)?
- Implement a source type for triangular-facet surfaces
- Implement sources not centered on the origin
    - Probably use a coordinate system transform to map observation point onto source-centric coordinates
