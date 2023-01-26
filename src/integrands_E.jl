"""
    __integrand_E(r̄′::CoordinateCartesian; r̄::CoordinateCartesian, t::Time,
                  source::AbstractJefimenkoSource{T}, media::PropagationMedia)::SVector{3,T}
                        where {T<:AbstractFloat}

Calculate the integrand function for the electric Jefimenko equation at the `source point
r̄′`. Parameterize the integrand function according to a particular `field source`,
`propagation media`, and for an observer positioned at space-time point (`r̄`,`t`).

# Arguments
- `r̄′::AbstractCoordinate`: spatial coordinate of the source point

# Parameters
- `r̄::UnitfulCoordinateSystems.AbstractCoordinate`: spatial coordinate at the observation point
- `t::Unitful.Time`: time at the observation point
- `source::JefimenkoSource`: the source model generating the electric field
- `media::PropagationMedia`: properties of the propagation media

# Returns
- `SVector{3,T}`: the predicted vector-valued integrand value
"""
function __integrand_E(r̄′::CoordinateCartesian; r̄::CoordinateCartesian, t::Unitful.Time,
                         source::AbstractJefimenkoSource{T}, media::PropagationMedia_Simple
                         )::SVector{3,T} where {T<:AbstractFloat}
    # Get spatial properties, in implicit units of meters
    Δr̄_m::SVector{3,T} = ustrip.(T, m, SVector(r̄ - r̄′))      #  vector r̄-r̄′
    r_m::T = norm(Δr̄_m)                                      #  magnitude |r̄-r̄′|

    # Get media properties, in implicit units as specified
    c::T = ustrip(T, m/s, media.c)               #  speed of light
    ε::T = ustrip(T, A*s/(V*m), media.ε)         #  permittivity

    # Calculate source-observer retarded time, in implicit units of seconds
    t′_s::T = ustrip(T, s, t′(r̄,t,r̄′,media))

    # Evaluate source function aliases, in implicit units as specified
    ρₑ::T = source.ρₑ(r̄′, t′_s)
    ∂ρₑ_∂t::T = ForwardDiff.derivative(t_s -> source.ρₑ(r̄′,t_s), t′_s)
    ∂Jₑ_∂t::SVector{3,T} = ForwardDiff.derivative(t_s -> source.Jₑ(r̄′,t_s), t′_s)
    Jₕ::SVector{3,T} = source.Jₕ(r̄′, t′_s)
    ∂Jₕ_∂t::SVector{3,T} = ForwardDiff.derivative(t_s -> source.Jₕ(r̄′,t_s), t′_s)

    # Calculate first term [V/m²]
    term1::SVector{3,T} = (ε^-1) .* (
                                        ((Δr̄_m ./ r_m^3) .* ρₑ)
                                        + ((Δr̄_m ./ r_m^2) .* (c^-1) .* ∂ρₑ_∂t)
                                        - ((1 / r_m) .* (c^-2) .* ∂Jₑ_∂t)
                                    )

    # Calculate second term [V/m²]
    term2::SVector{3,T} = LinearAlgebra.cross(((Jₕ ./ r_m^3) + ((1 / r_m^2) .* (c^-1) .* ∂Jₕ_∂t)), Δr̄_m)

    return (term1 - term2)
end
