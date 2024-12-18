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
function E(r̄::Meshes.Point, t::Unitful.Time, model::JefimenkoModel)
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
function H(r̄::Meshes.Point, t::Unitful.Time, model::JefimenkoModel)
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
function P(r̄::Meshes.Point, t::Unitful.Time, model::JefimenkoModel)
    Ert = E(r̄, t, model)
    Hrt = H(r̄, t, model)
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
    media::PropagationMedia_Simple,
    rule::MeshIntegrals.IntegrationRule = MeshIntegrals.HAdaptiveCubature()
)
    integrand = _integrand_E(r̄′; source=source, media=media, r̄=CoordinateCartesian(r̄), t=t)
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
    media::PropagationMedia_Simple,
    rule::MeshIntegrals.IntegrationRule = MeshIntegrals.HAdaptiveCubature()
)
    function integrand_Am4(u, p)
        r̄′ = CoordinateCartesian(u[1]*m, u[2]*m, u[3]*m)
        return _integrand_H(r̄′; source=source, media=media, r̄=CoordinateCartesian(r̄), t=t)
    end

    return (1/4π) .* MeshIntegrals.integral(integrand, source.geometry, rule)
end
