module JefimenkoModels
    import Meshes
    
    using Unitful, UnitfulCoordinateSystems
    using Unitful.DefaultSymbols: W, A, V, C, m, s, rad
    using PhysicalConstants.CODATA2018: c_0, ε_0, μ_0
    using ForwardDiff, Integrals, LinearAlgebra, StaticArrays

    __DEFAULT_RTOL = sqrt(eps())

    ###########################################################################
    #                     DATA STRUCTURES & COMMON DEFAULTS
    ###########################################################################

    # Data structures
    include("structs.jl")
    export RadiationSource

    include("accessors.jl")

    CLASSICAL_VACUUM = let
        ε₀ = uconvert((A*s)/(V*m), float(ε_0))
        μ₀ = uconvert((V*s)/(A*m), float(μ_0))
        c₀ = uconvert(m/s, float(c_0))
        PropagationMedia_Simple(ε₀, μ₀, c₀)
    end

    NULL_CHARGE(r̄::AbstractCoordinate, t_s::Real) = 0
    NULL_CURRENT(r̄::AbstractCoordinate, t_s::Real) = StaticArrays.SVector(0, 0, 0)
    
    export CLASSICAL_VACUUM, NULL_CHARGE, NULL_CURRENT

    ###########################################################################
    #                     RETARDED-TIME CALCULATIONS
    ###########################################################################

    """
        t′(r̄::AbstractCoordinate, t:Time, r̄′::Coordinate, c::Quantity)

    Calculate the retarded-time at a source point `r̄′` for an observer at the space-time
    point (`r̄`,`t`) through a medium with speed of light `c`.

    # Arguments
    - `r̄::UnitfulCoordinateSystems.AbstractCoordinate`: spatial location of the observation point
    - `t::Unitful.Time`: time at the observation point
    - `r̄′::UnitfulCoordinateSystems.AbstractCoordinate`: spatial location of the source point
    - `c::Quantity`: Unitful speed of light in the medium between r̄′ and r̄
    """
    function t′(r̄::AbstractCoordinate, t::Unitful.Time, r̄′::AbstractCoordinate, c::Quantity)::Unitful.Time
        return (t - (norm(r̄-r̄′)/c))
    end

    """
        t′(r̄::Coordinate, t:Time, r̄′::Coordinate, media::PropagationMedia)

    Calculate the retarded-time at a source point `r̄′` for an observer at the space-time
    point (`r̄`,`t`) through a `propagation medium`.

    # Arguments
    - `r̄::UnitfulCoordinateSystems.Coordinate`: spatial location of the observation point
    - `t::Unitful.Time`: time at the observation point
    - `r̄′::UnitfulCoordinateSystems.Coordinate`: spatial location of the source point
    - `media::PropagationMedia`: properties of the medium between r̄′ and r̄
    """
    function t′(r̄::AbstractCoordinate, t::Unitful.Time, r̄′::AbstractCoordinate, media::PropagationMedia_Simple)::Unitful.Time
        return t′(r̄, t, r̄′, media.c)
    end

    function t′(r̄::AbstractCoordinate, t::Unitful.Time, r̄′::AbstractCoordinate, media::PropagationMedia_DiagonallyAnisotropic)::Unitful.Time
        Δr̄ = SVector(r̄ - r̄′)
        Δt = norm(media.c^-1 * Δr̄) |> unit(t)
        return (t - Δt)
    end

    export t′

    ###########################################################################
    #                     EM FIELD CALCULATIONS
    ###########################################################################

    """
        H(r̄::AbstractCoordinate, t::Time, model::JefimenkoModel; rtol=sqrt(eps))

    Calculate the predicted electric field 𝐇 observed at space-time point (`r̄`,`t`) using
    the electric Jefimenko equation for a particular `model`. Calculate the integral using
    a specified `relative tolerance`.

    # Arguments
    - `r̄::UnitfulCoordinateSystems.AbstractCoordinate`: spatial location of the observation point
    - `t::Unitful.Time`: time at which the electric field is observed
    - `model::JefimenkoModel`: model of the transmitting source and propagation media

    # Keywords
    - `rtol::Real`: relative tolerance at which to solve the integral (optional)
    """
    function E(r̄::AbstractCoordinate, t::Unitful.Time, model::JefimenkoModel; rtol=__DEFAULT_RTOL)
        # Superimpose the contributions of the E(r̄,t) produced by each source in model
        E_contrib(source) = __E(r̄, t, source, model.media; rtol=rtol)
        return mapreduce(E_contrib, +, model.sources)
    end

    """
        H(r̄::AbstractCoordinate, t::Time, model::JefimenkoModel; rtol=sqrt(eps))

    Calculate the predicted magnetic field 𝐇 observed at space-time point (`r̄`,`t`) using
    the magnetic Jefimenko equation for a particular `model`. Calculate the integral using
    a specified `relative tolerance`.

    # Arguments
    - `r̄::UnitfulCoordinateSystems.AbstractCoordinate`: spatial location of the observation point
    - `t::Unitful.Time`: time at which the field is observed
    - `model::JefimenkoModel`: model of the transmitting source and propagation media

    # Keywords
    - `rtol::Real`: relative tolerance at which to solve the integral (optional)
    """
    function H(r̄::AbstractCoordinate, t::Unitful.Time, model::JefimenkoModel; rtol=__DEFAULT_RTOL)
        # Superimpose the contributions of the 𝐇(r̄,t) produced by each source in model
        H_contrib(source) = __H(r̄, t, source, model.media; rtol=rtol)
        return mapreduce(H_contrib, +, model.sources) 
    end

    """
        P(r̄::AbstractCoordinate, t::Time, model::JefimenkoModel; rtol=sqrt(eps))

    Calculate the predicted Poynting vector 𝐏 observed at space-time point (`r̄`,`t`) using
    the electric and magnetic Jefimenko equations for a particular `model`. Calculate the
    integrals using a specified `relative tolerance`.

    # Arguments
    - `r̄::UnitfulCoordinateSystems.AbstractCoordinate`: spatial location of the observation point
    - `t::Unitful.Time`: time at which the field is observed
    - `model::JefimenkoModel`: model of the transmitting source and propagation media

    # Keywords
    - `rtol::Real`: relative tolerance at which to solve the integral (optional)
    """
    function P(r̄::AbstractCoordinate, t::Unitful.Time, model::JefimenkoModel; rtol=__DEFAULT_RTOL)
        Ert = E(r̄,t,model; rtol=rtol)
        Hrt = H(r̄,t,model; rtol=rtol)
        return cross(Ert,Hrt) .|> W/m^2
    end

    """
        __P(r̄::AbstractCoordinate, t::Time, source::JefimenkoSource, media::PropagationMedia; rtol)

    Calculate the predicted Poynting vector 𝐏 observed at space-time point (`r̄`,`t`) due to
    a particular `source`, transmitted through a particular `propagation media`. Calculate
    the integral using a specified `relative tolerance`.

    # Arguments
    - `r̄::UnitfulCoordinateSystems.AbstractCoordinate`: spatial location of the observation point
    - `t::Unitful.Time`: time at which the electric field is observed
    - `source::JefimenkoSource`: source of the electric field
    - `media::PropagationMedia`: properties of the propagation media

    # Keywords
    - `rtol::Real`: relative tolerance at which to solve the integral (optional)
    """
    function __P(r̄::AbstractCoordinate, t::Unitful.Time, source::AbstractJefimenkoSource{T},
                media::AbstractPropagationMedia; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
        Ert = __E(r̄,t,source,media; rtol=rtol)
        Hrt = __H(r̄,t,source,media; rtol=rtol)
        return cross(Ert,Hrt) .|> W/m^2
    end

    include("integrands.jl")
    include("fields_E.jl")
    include("fields_H.jl")

    export E, H, P
end
