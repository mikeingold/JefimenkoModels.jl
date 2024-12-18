"""
    _integrand_E(r̄′, r̄, t, source, media) -> SVector{3, Quantity}

Calculate the integrand function for the electric Jefimenko equation at the `source point
r̄′`. Parameterize the integrand function according to a particular `field source`,
`propagation media`, and for an observer positioned at space-time point (`r̄`,`t`).

# Arguments
- `r̄′::Meshes.Point`: coordinate of the source point

# Parameters
- `r̄::Meshes.Point`: coordinate of the observation point
- `t::Unitful.Time`: time at the observation point
- `source::JefimenkoSource`: the source model generating the field
- `media::PropagationMedia_Simple`: properties of the propagation media

# Returns
- `SVector{3, Quantity}`: the predicted vector-valued integrand value
"""
function _integrand_E(r̄′::Meshes.Point, r̄::Meshes.Point, t::Unitful.Time, source::RadiationSource, media::PropagationMedia_Simple)
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

"""
    _integrand_H(r̄′, r̄, t, source, media) -> SVector{3, Quantity}

Calculate the integrand function for the magnetic Jefimenko equation at the `source point
r̄′`. Parameterize the integrand function according to a particular `field source`,
`propagation media`, and for an observer positioned at space-time point (`r̄`,`t`).

# Arguments
- `r̄′::Meshes.Point`: coordinate of the source point

# Parameters
- `r̄::Meshes.Point`: coordinate of the observation point
- `t::Unitful.Time`: time at the observation point
- `source::JefimenkoSource`: the source model generating the field
- `media::PropagationMedia_Simple`: properties of the propagation media

# Returns
- `SVector{3, Quantity}`: the predicted vector-valued integrand value
"""
function _integrand_H(r̄′::Meshes.Point, r̄::Meshes.Point, t::Unitful.Time, source::RadiationSource, media::PropagationMedia_Simple)
    # Get spatial properties, in implicit units of meters
    Δr̄_m::SVector{3,T} = ustrip.(T, m, SVector(r̄ - r̄′))      #  SVector r̄-r̄′
    r_m::T = norm(Δr̄_m)                                      #  magnitude |r̄-r̄′|

    # Get media properties, in implicit units as specified
    c::T = ustrip(T, m/s, media.c)               #  speed of light in [m/s]
    μ::T = ustrip(T, (V*s)/(A*m), media.μ)       #  permeability in [Vs/Am]

    # Calculate source-observer retarded time, in implicit units of seconds
    t′_s::T = ustrip(T, s, t′(r̄,t,r̄′,media))

    # Source functions
    ρₕ::T = source.ρₕ(r̄′, t′_s)
    ∂ρₕ_∂t::T = ForwardDiff.derivative(t_s -> source.ρₕ(r̄′,t_s), t′_s)
    Jₑ::SVector{3,T} = source.Jₑ(r̄′, t′_s)
    ∂Jₑ_∂t::SVector{3,T} = ForwardDiff.derivative(t_s -> source.Jₑ(r̄′,t_s), t′_s)
    ∂Jₕ_∂t::SVector{3,T} = ForwardDiff.derivative(t_s -> source.Jₕ(r̄′,t_s), t′_s)

    # Calculate first term
    term1::SVector{3,T} = (μ^-1) .* (
                                        ((Δr̄_m ./ r_m^3) .* ρₕ)
                                        + ((Δr̄_m ./ r_m^2) .* (c^-1) .* ∂ρₕ_∂t)
                                        - ((1 / r_m) .* (c^-2) .* ∂Jₕ_∂t)
                                    )
    
    # Calculate second term
    term2::SVector{3,T} = LinearAlgebra.cross((Jₑ ./ r_m^3) + ((1 / r_m^2) .* (c^-1) .* ∂Jₑ_∂t), Δr̄_m)

    return (term1 + term2)
end
