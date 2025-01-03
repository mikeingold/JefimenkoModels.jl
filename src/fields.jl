###########################################################################
#                     RETARDED-TIME CALCULATIONS
###########################################################################

"""
    _t′(r̄, t, r̄′, c)

Calculate the retarded-time at a source point `r̄′` for an observer at the space-time
point (`r̄`,`t`) through a medium with speed of light `c`.

# Arguments
- `r̄::Meshes.Point`: spatial location of the observation point
- `t::Unitful.Time`: time at the observation point
- `r̄′::Meshes.Point`: spatial location of the source point
- `c::Quantity`: Unitful speed of light in the medium between r̄′ and r̄
"""
function _t′(r̄::Meshes.Point, t::Unitful.Time, r̄′::Meshes.Point, c::Quantity)
    return (t - (LinearAlgebra.norm(r̄ - r̄′) / c))
end

"""
    _t′(r̄, t, r̄′, media)

Calculate the retarded-time at a source point `r̄′` for an observer at the space-time
point (`r̄`,`t`) through a `propagation medium`.

# Arguments
- `r̄::Meshes.Point`: spatial location of the observation point
- `t::Unitful.Time`: time at the observation point
- `r̄′::Meshes.Point`: spatial location of the source point
- `media::SimpleMedia`: properties of the medium between r̄′ and r̄
"""
function _t′(r̄::Meshes.Point, t::Unitful.Time, r̄′::Meshes.Point, media::SimpleMedia)
    return t′(r̄, t, r̄′, media.c)
end

#=
function _t′(r̄::Meshes.Point, t::Unitful.Time, r̄′::Meshes.Point, media::PropagationMedia_DiagonallyAnisotropic)
    Δr̄ = r̄ - r̄′
    Δt = norm(media.c^-1 * Δr̄) |> unit(t)
    return (t - Δt)
end
=#

###########################################################################
#                     EM FIELD CALCULATIONS
###########################################################################

"""
    E(r̄, t, model)

Calculate the predicted electric field 𝐇 observed at space-time point (`r̄`,`t`) using
the electric Jefimenko equation for a particular `model`. Calculate the integral using
a specified `relative tolerance`.

# Arguments
- `r̄::UnitfulCoordinateSystems.AbstractCoordinate`: spatial location of the observation point
- `t::Unitful.Time`: time at which the electric field is observed
- `model::JefimenkoModel`: model of the transmitting source and propagation media
"""
function E(
    r̄::Meshes.Point,
    t::Unitful.Time,
    model::JefimenkoModel,
    rule::MeshIntegrals.IntegrationRule = MeshIntegrals.HAdaptiveCubature()
)
    # Superimpose the contributions of the E(r̄,t) produced by each source in model
    return mapreduce(source -> _E(r̄, t, source, model.media, rule), +, model.sources)
end

"""
    H(r̄, t, model)

Calculate the predicted magnetic field 𝐇 observed at space-time point (`r̄`,`t`) using
the magnetic Jefimenko equation for a particular `model`. Calculate the integral using
a specified `relative tolerance`.

# Arguments
- `r̄::UnitfulCoordinateSystems.AbstractCoordinate`: spatial location of the observation point
- `t::Unitful.Time`: time at which the field is observed
- `model::JefimenkoModel`: model of the transmitting source and propagation media
"""
function H(
    r̄::Meshes.Point,
    t::Unitful.Time,
    model::JefimenkoModel,
    rule::MeshIntegrals.IntegrationRule = MeshIntegrals.HAdaptiveCubature()
)
    # Superimpose the contributions of the 𝐇(r̄,t) produced by each source in model
    return mapreduce(source -> _H(r̄, t, source, model.media, rule), +, model.sources) 
end

"""
    P(r̄, t, model)

Calculate the predicted Poynting vector 𝐏 observed at space-time point (`r̄`,`t`) using
the electric and magnetic Jefimenko equations for a particular `model`. Calculate the
integrals using a specified `relative tolerance`.

# Arguments
- `r̄::UnitfulCoordinateSystems.AbstractCoordinate`: spatial location of the observation point
- `t::Unitful.Time`: time at which the field is observed
- `model::JefimenkoModel`: model of the transmitting source and propagation media
"""
function P(
    r̄::Meshes.Point,
    t::Unitful.Time,
    model::JefimenkoModel,
    rule::MeshIntegrals.IntegrationRule = MeshIntegrals.HAdaptiveCubature()
)
    Ert = E(r̄, t, model, rule)
    Hrt = H(r̄, t, model, rule)
    return cross(Ert, Hrt) .|> W/m^2
end

###########################################################################
#                              WORKERS
###########################################################################

"""
    _E(r̄, t, source, media, rule = HAdaptiveCubature())

Calculate the electric field at (`r̄`,`t`) using the electric Jefimenko equation due to a
particular `source`, transmitted through a particular homogeneous `propagation media`.
Calculate the integral using a specified `integration rule`.

# Arguments
- `r̄::Meshes.Point`: spatial location of the observation point
- `t::Unitful.Time`: time at which the electric field is observed
- `source::RadiationSource`: source of the electric field
- `media::PropagationMedia`: properties of the propagation media
- `rule::MeshIntegrals.IntegrationRule`: rule to use for numerical integration
"""
function _E(
    r̄::Meshes.Point,
    t::Unitful.Time,
    source::RadiationSource,
    media::SimpleMedia,
    rule::MeshIntegrals.IntegrationRule = MeshIntegrals.HAdaptiveCubature()
)
    integrand(r̄′) = _integrand_E(r̄′, r̄, t, source, media)
    return (1/4π) .* MeshIntegrals.integral(integrand, source.geometry, rule)
end

"""
    _H(r̄, t, source, media, rule)

Calculate the magnetic field at (`r̄`,`t`) using the magnetic Jefimenko equation due to a
particular `source`, transmitted through a particular homogeneous `propagation media`.
Calculate the integral using a specified `integration rule`.

# Arguments
- `r̄::Meshes.Point`: spatial location of the observation point
- `t::Unitful.Time`: time at which the magnetic field is observed
- `source::RadiationSource`: source of the magnetic field
- `media::PropagationMedia`: properties of the propagation media
- `rule::MeshIntegrals.IntegrationRule`: rule to use for numerical integration
"""
function _H(
    r̄::Meshes.Point,
    t::Unitful.Time,
    source::RadiationSource,
    media::SimpleMedia,
    rule::MeshIntegrals.IntegrationRule = MeshIntegrals.HAdaptiveCubature()
)
    integrand(r̄′) = _integrand_H(r̄′, r̄, t, source, media)
    return (1/4π) .* MeshIntegrals.integral(integrand, source.geometry, rule)
end

###########################################################################
#                             INTEGRANDS
###########################################################################

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
- `media::SimpleMedia`: properties of the propagation media

# Returns
- `SVector{3, Quantity}`: the predicted vector-valued integrand value
"""
function _integrand_E(
    r̄′::Meshes.Point,
    r̄::Meshes.Point,
    t::Unitful.Time,
    source::RadiationSource,
    media::SimpleMedia
)
    Δr̄ = r̄ - r̄′
    r = LinearAlgebra.norm(Δr̄)
    t′ = _t′(r̄, t, r̄′, media)
    ε = media.permittivity
    c = media.c

    # Source functions
    ρₑ = source.rho_e(r̄′, t′)
    ∂ρₑ_∂t = ForwardDiff.derivative(t -> source.rho_e(r̄′, t), t′)
    ∂Jₑ_∂t = ForwardDiff.derivative(t -> source.J_e(r̄′, t), t′)
    Jₕ = source.J_h(r̄′, t′)
    ∂Jₕ_∂t = ForwardDiff.derivative(t -> source.J_h(r̄′, t), t′)

    # Calculate integrand terms
    term1 = (ε^-1) .* ( ((Δr̄ ./ r^3) .* ρₑ) + ((Δr̄ ./ r^2) .* (c^-1) .* ∂ρₑ_∂t) - ((1 / r) .* (c^-2) .* ∂Jₑ_∂t) )
    term2 = LinearAlgebra.cross(((Jₕ ./ r^3) + ((1 / r^2) .* (c^-1) .* ∂Jₕ_∂t)), Δr̄)

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
- `media::SimpleMedia`: properties of the propagation media

# Returns
- `SVector{3, Quantity}`: the predicted vector-valued integrand value
"""
function _integrand_H(
    r̄′::Meshes.Point,
    r̄::Meshes.Point,
    t::Unitful.Time,
    source::RadiationSource,
    media::SimpleMedia
)
    Δr̄ = r̄ - r̄′
    r = LinearAlgebra.norm(Δr̄)
    t′ = _t′(r̄, t, r̄′, media)
    μ = media.permeability
    c = media.c

    # Source functions
    ρₕ = source.rho_h(r̄′, t′)
    ∂ρₕ_∂t = ForwardDiff.derivative(t -> source.rho_h(r̄′, t), t′)
    Jₑ = source.J_e(r̄′, t′)
    ∂Jₑ_∂t = ForwardDiff.derivative(t -> source.J_e(r̄′, t), t′)
    ∂Jₕ_∂t = ForwardDiff.derivative(t -> source.J_h(r̄′, t), t′)

    # Calculate integrand terms
    term1 = (μ^-1) .* ( ((Δr̄_m ./ r_m^3) .* ρₕ) + ((Δr̄_m ./ r_m^2) .* (c^-1) .* ∂ρₕ_∂t) - ((1 / r_m) .* (c^-2) .* ∂Jₕ_∂t) )
    term2 = LinearAlgebra.cross((Jₑ ./ r_m^3) + ((1 / r_m^2) .* (c^-1) .* ∂Jₑ_∂t), Δr̄_m)

    return (term1 + term2)
end
