# TODO

## Near Term
- Implement `Base.show` for structs
- Document structs and constructor methods

## Longer-Term / Vision
- Improve documentation for source function definitions (units, types, etc)
    - Create a test/inspect/validate function for users to look for issues in their definitions?
- Re-assess the need for solver numerical type parameterization
    - Does it even work as intended? Is there a performance benefit?
    - Should __DEFAULT_RTOL be defined as a function __DEFAULT_RTOL(T)?
- Implement a source type for triangular-facet surfaces
