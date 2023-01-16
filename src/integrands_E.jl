"""
    __integrand_E_R1(r̄′::Coordinate; source::AbstractLineSource{T}, media::PropagationMedia
                    r̄::Coordinate, t::Time)  where {T<:AbstractFloat}

Calculate the integrand function for the electric Jefimenko equation of a one-dimensional
`line source` at the `source point r̄′`. Parameterize the integrand function for a particular
`field source`, `propagation media`, and an observer positioned at space-time point (`r̄`,`t`).

# Arguments
- `r̄′::Coordinate`: spatial coordinate of the source point

# Parameters
- `r̄::UnitfulCoordinateSystems.Coordinate`: spatial coordinate at the observation point
- `t::Unitful.Time`: time at the observation point
- `source::JefimenkoSource`: the source model generating the electric field
- `media::PropagationMedia`: properties of the propagation media

# Returns
- `Vector{Quantity}`: the predicted vector-valued integrand value with appropriate units
"""
function __integrand_E_R1(r̄′::Coordinate; r̄::Coordinate, t::Unitful.Time,
                         source::AbstractLineSource{T}, media::PropagationMedia_Simple
                         ) where {T<:AbstractFloat}
    # Get spatial properties, in implicit units of meters
    r̄′ = CoordinateCartesian(r̄′)                 # ensure r̄′ is in Cartesian format
    r̄ = CoordinateCartesian(r̄)                   # ensure r̄′ is in Cartesian format
    Δr̄_m = ustrip.(T, m, SVector(r̄ - r̄′))        #  vector r̄-r̄′
    r_m = norm(Δr̄_m)                             #  magnitude |r̄-r̄′|

    # Get media properties, in implicit units as specified
    c = ustrip(T, m/s, media.c)               #  speed of light in [m/s]
    ε = ustrip(T, A*s/(V*m), media.ε)         #  permittivity in [(A s)/(V m)]

    # Calculate source-observer retarded time, in implicit units of seconds
    t′_s::T = ustrip(T, s, t′(r̄,t,r̄′,media))

    # Evaluate source function aliases, in implicit units as specified
    ρₑ = source.ρₑ(r̄′, t′_s)                                          # in [C m^-1]
    ∂ρₑ_∂t = ForwardDiff.derivative(t_s -> source.ρₑ(r̄′,t_s), t′_s)   # in [C m^-1 s^-1]
    ∂Jₑ_∂t = ForwardDiff.derivative(t_s -> source.Jₑ(r̄′,t_s), t′_s)   # in [A s^-1]
    Jₕ = source.Jₕ(r̄′, t′_s)                                           # in [V]
    ∂Jₕ_∂t = ForwardDiff.derivative(t_s -> source.Jₕ(r̄′,t_s), t′_s)    # in [V s^-1]

    # Calculate first term, dimensional analysis of implied units commented on right
    term1a = ( (Δr̄_m ./ r_m^3) .* ρₑ )                  # [m/m³ * C/m]        -> [As/m³]
    term1b = ( (Δr̄_m ./ r_m^2) .* (c^-1) .* ∂ρₑ_∂t )    # [m/m² * s/m * C/ms] -> [As/m³]
    term1c = ( (1 / r_m) .* (c^-2) .* ∂Jₑ_∂t )          # [1/m * s²/m² * A/s] -> [As/m³]
    term1  = ( (ε^-1) .* (term1a + term1b - term1c) )   # [Vm/As * As/m³]     -> [V/m²]
    
    # Calculate second term, dimensional analysis of implied units commented on right
    term2a = ( Jₕ ./ r_m^3 )                          # [V/m³] -> [V/m³]
    term2b = ( (1 / r_m^2) .* (c^-1) .* ∂Jₕ_∂t )      # [1/m² * s/m * V/s] -> [V/m³]
    term2  = cross((term2a + term2b), Δr̄_m)           # [V/m³ * m] -> [V/m²]

    # Combine terms, apply appropriate units
    return (term1 - term2) * (V/m^2)
end

"""
    __integrand_E_R2(r̄′::Coordinate; source::AbstractSurfaceSource{T}, media::PropagationMedia
                    r̄::Coordinate, t::Time) where {T<:AbstractFloat}

Calculate the integrand function for the electric Jefimenko equation of a two-dimensional
`surface source` at the `source point r̄′`. Parameterize the integrand function for a particular
`field source`, `propagation media`, and an observer positioned at space-time point (`r̄`,`t`).

# Arguments
- `r̄′::Coordinate`: spatial coordinate of the source point

# Parameters
- `r̄::UnitfulCoordinateSystems.Coordinate`: spatial coordinate at the observation point
- `t::Unitful.Time`: time at the observation point
- `source::JefimenkoSource`: the source model generating the electric field
- `media::PropagationMedia`: properties of the propagation media

# Returns
- `Vector{T}`: the predicted vector-valued integrand value, in implicit units [V/m³]. (The
`HCubature` solver does not yet support `Unitful` types.)
"""
function __integrand_E_R2(r̄′::Coordinate; r̄::Coordinate, t::Unitful.Time,
                         source::AbstractSurfaceSource{T}, media::PropagationMedia_Simple
                         ) where {T<:AbstractFloat}
    # Get spatial properties, in implicit units of meters
    r̄′ = CoordinateCartesian(r̄′)                 # ensure r̄′ is in Cartesian format
    r̄ = CoordinateCartesian(r̄)                   # ensure r̄′ is in Cartesian format
    Δr̄_m = ustrip.(T, m, SVector(r̄ - r̄′))        #  vector r̄-r̄′
    r_m = norm(Δr̄_m)                             #  magnitude |r̄-r̄′|

    # Get media properties, in implicit units as specified
    c = ustrip(T, m/s, media.c)               #  speed of light in [m/s]
    ε = ustrip(T, A*s/(V*m), media.ε)         #  permittivity in [(A s)/(V m)]

    # Calculate source-observer retarded time, in implicit units of seconds
    t′_s::T = ustrip(T, s, t′(r̄,t,r̄′,media))

    # Source functions
    ρₑ = source.ρₑ(r̄′, t′_s)                                          # in [C m^-2]
    ∂ρₑ_∂t = ForwardDiff.derivative(t_s -> source.ρₑ(r̄′,t_s), t′_s)   # in [C m^-2 s^-1]
    ∂Jₑ_∂t = ForwardDiff.derivative(t_s -> source.Jₑ(r̄′,t_s), t′_s)   # in [A m^-1 s^-1]
    Jₕ = source.Jₕ(r̄′, t′_s)                                           # in [V m^-1]
    ∂Jₕ_∂t = ForwardDiff.derivative(t_s -> source.Jₕ(r̄′,t_s), t′_s)    # in [V m^-1 s^-1]

    # Calculate first term
    term1a = ( (Δr̄_m ./ r_m^3) .* ρₑ )                  # [m/m³ * C/m²]     -> [As/m⁴]
    term1b = ( (Δr̄_m ./ r_m^2) .* (c^-1) .* ∂ρₑ_∂t )    # [m/m² * s/m * C/m²s] -> [As/m⁴]
    term1c = ( (1 / r_m) .* (c^-2) .* ∂Jₑ_∂t )          # [1/m * s²/m² * A/sm] -> [As/m⁴]
    term1  = ( (ε^-1) .* (term1a + term1b - term1c) )   # [Vm/As * As/m⁴] -> [V/m³]
    
    # Calculate second term
    term2a = ( Jₕ ./ r_m^3 )                            # [V/m / m³] -> [V/m⁴]
    term2b = ( (1 / r_m^2) .* (c^-1) .* ∂Jₕ_∂t )        # [1/m² * s/m * V/sm] -> [V/m⁴]
    term2  = cross((term2a + term2b), Δr̄_m)             # [V/m⁴ * m] -> [V/m³]

    # Combine terms and return, in implicit units of [V/m³] 
    return (term1 - term2)
end

"""
    __integrand_E_R3(r̄′::Coordinate; source::AbstractVolumeSource{T}, media::PropagationMedia
                    r̄::Coordinate, t::Time) where {T<:AbstractFloat}

Calculate the integrand function for the electric Jefimenko equation of a three-dimensional
`volume source` at the `source point r̄′`. Parameterize the integrand function for a particular
`field source`, `propagation media`, and an observer positioned at space-time point (`r̄`,`t`).

# Arguments
- `r̄′::Coordinate`: spatial coordinate of the source point

# Parameters
- `r̄::UnitfulCoordinateSystems.Coordinate`: spatial coordinate at the observation point
- `t::Unitful.Time`: time at the observation point
- `source::JefimenkoSource`: the source model generating the electric field
- `media::PropagationMedia`: properties of the propagation media

# Returns
- `Vector{T}`: the predicted vector-valued integrand value, in implicit units [V/m⁴]. (The
`HCubature` solver does not yet support `Unitful` types.)
"""
function __integrand_E_R3(r̄′::Coordinate; r̄::Coordinate, t::Unitful.Time,
                         source::AbstractVolumeSource{T}, media::PropagationMedia_Simple
                         ) where {T<:AbstractFloat}
    # Get spatial properties, in implicit units of meters
    r̄′ = CoordinateCartesian(r̄′)                 # ensure r̄′ is in Cartesian format
    r̄ = CoordinateCartesian(r̄)                   # ensure r̄′ is in Cartesian format
    Δr̄_m = ustrip.(T, m, SVector(r̄ - r̄′))        #  vector r̄-r̄′
    r_m = norm(Δr̄_m)                             #  magnitude |r̄-r̄′|

    # Get media properties, in implicit units as specified
    c = ustrip(T, m/s, media.c)               #  speed of light in [m/s]
    ε = ustrip(T, A*s/(V*m), media.ε)         #  permittivity in [(A s)/(V m)]

    # Calculate source-observer retarded time, in implicit units of seconds
    t′_s::T = ustrip(T, s, t′(r̄,t,r̄′,media))

    # Source functions
    ρₑ = source.ρₑ(r̄′, t′_s)                                          # in [C m^-3]
    ∂ρₑ_∂t = ForwardDiff.derivative(t_s -> source.ρₑ(r̄′,t_s), t′_s)   # in [C m^-3 s^-1]
    ∂Jₑ_∂t = ForwardDiff.derivative(t_s -> source.Jₑ(r̄′,t_s), t′_s)   # in [A m^-2 s^-1]
    Jₕ = source.Jₕ(r̄′, t′_s)                                           # in [V m^-2]
    ∂Jₕ_∂t = ForwardDiff.derivative(t_s -> source.Jₕ(r̄′,t_s), t′_s)    # in [V m^-2 s^-1]

    # Calculate first term
    term1a = ( (Δr̄_m ./ r_m^3) .* ρₑ )                  # [m/m³ * C/m³]     -> [As/m⁵]
    term1b = ( (Δr̄_m ./ r_m^2) .* (c^-1) .* ∂ρₑ_∂t )    # [m/m² * s/m * C/m³s] -> [As/m⁵]
    term1c = ( (1 / r_m) .* (c^-2) .* ∂Jₑ_∂t )          # [1/m * s²/m² * A/m²s] -> [As/m⁵]
    term1  = ( (ε^-1) .* (term1a + term1b - term1c) )   # [Vm/As * As/m⁵] -> [V/m⁴]
    
    # Calculate second term
    term2a = ( Jₕ ./ r_m^3 )                            # [V/m² / m³] -> [V/m⁵]
    term2b = ( (1 / r_m^2) .* (c^-1) .* ∂Jₕ_∂t )        # [1/m² * s/m * V/m²s] -> [V/m⁵]
    term2  = cross((term2a + term2b), Δr̄_m)             # [V/m⁵ * m] -> [V/m⁴]

    # Combine terms and return, in implicit units of [V/m⁴]
    return (term1 - term2)
end
