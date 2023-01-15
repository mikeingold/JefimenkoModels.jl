module JefimenkoModels
    using ForwardDiff, HCubature, Integrals
    using LinearAlgebra
    using StaticArrays
    using Unitful
    using Unitful.DefaultSymbols: W, A, V, C, m, s, rad
    using UnitfulCoordinateSystems

    __DEFAULT_RTOL = sqrt(eps())

    # Data structures
    include("structs.jl")

    ###########################################################################
    #                     RETARDED-TIME CALCULATIONS
    ###########################################################################

    tᵣ(r̄::Coordinate, t::Unitful.Time, r̄′::Coordinate, c::Quantity)::Unitful.Time = t - (norm(r̄-r̄′)/c)

    tᵣ(r̄::Coordinate, t::Unitful.Time, r̄′::Coordinate, media::PropagationMedia_Simple)::Unitful.Time = t - (norm(r̄-r̄′)/media.c)

    function tᵣ(r̄::Coordinate, t::Unitful.Time, r̄′::Coordinate, media::PropagationMedia_DiagonallyAnisotropic)::Unitful.Time
        Δr̄ = SVector(r̄ - r̄′)
        Δt = norm(media.c^-1 * Δr̄) |> unit(t)
        return (t - Δt)
    end

    ###########################################################################
    #                     EM FIELD CALCULATIONS
    ###########################################################################

    """
        𝐄(r̄::Coordinate, t::Time, model::JefimenkoModel; rtol=sqrt(eps))

    Calculate the predicted electric field 𝐇 observed at space-time point (`r̄`,`t`) using
    the electric Jefimenko equation for a particular `model`. Calculate the integral using
    a specified `relative tolerance`.

    # Arguments
    - `r̄::UnitfulCoordinateSystems.Coordinate`: spatial location of the observation point
    - `t::Unitful.Time`: time at which the electric field is observed
    - `model::JefimenkoModel`: model of the transmitting source and propagation media

    # Keywords
    - `rtol::Real`: relative tolerance at which to solve the integral (optional)
    """
    function 𝐄(r̄::Coordinate, t::Unitful.Time, model::JefimenkoModel; rtol=__DEFAULT_RTOL)
        # Superimpose the contributions of the 𝐄(r̄,t) produced by each source in model
        E_contrib(source) = _𝐄(r̄, t, source, model.media; rtol=rtol)
        return mapreduce(E_contrib, +, model.sources)
    end

    """
        𝐇(r̄::Coordinate, t::Time, model::JefimenkoModel; rtol=sqrt(eps))

    Calculate the predicted magnetic field 𝐇 observed at space-time point (`r̄`,`t`) using
    the magnetic Jefimenko equation for a particular `model`. Calculate the integral using
    a specified `relative tolerance`.

    # Arguments
    - `r̄::UnitfulCoordinateSystems.Coordinate`: spatial location of the observation point
    - `t::Unitful.Time`: time at which the field is observed
    - `model::JefimenkoModel`: model of the transmitting source and propagation media

    # Keywords
    - `rtol::Real`: relative tolerance at which to solve the integral (optional)
    """
    function 𝐇(r̄::Coordinate, t::Unitful.Time, model::JefimenkoModel; rtol=__DEFAULT_RTOL)
        # Superimpose the contributions of the 𝐇(r̄,t) produced by each source in model
        H_contrib(source) = _𝐇(r̄, t, source, model.media; rtol=rtol)
        return mapreduce(H_contrib, +, model.sources) 
    end

    """
        𝐏(r̄::Coordinate, t::Time, model::JefimenkoModel; rtol=sqrt(eps))

    Calculate the predicted Poynting vector 𝐏 observed at space-time point (`r̄`,`t`) using
    the electric and magnetic Jefimenko equations for a particular `model`. Calculate the
    integrals using a specified `relative tolerance`.

    # Arguments
    - `r̄::UnitfulCoordinateSystems.Coordinate`: spatial location of the observation point
    - `t::Unitful.Time`: time at which the field is observed
    - `model::JefimenkoModel`: model of the transmitting source and propagation media

    # Keywords
    - `rtol::Real`: relative tolerance at which to solve the integral (optional)
    """
    function 𝐏(r̄::Coordinate, t::Unitful.Time, model::JefimenkoModel; rtol=__DEFAULT_RTOL)
        E = 𝐄(r̄,t,model; rtol=rtol)
        H = 𝐇(r̄,t,model; rtol=rtol)
        return cross(E,H) .|> W/m^2
    end

    """
        _𝐏(r̄::Coordinate, t::Time, source::JefimenkoSource, media::PropagationMedia; rtol)

    Calculate the predicted Poynting vector 𝐏 observed at space-time point (`r̄`,`t`) due to
    a particular `source`, transmitted through a particular `propagation media`. Calculate
    the integral using a specified `relative tolerance`.

    # Arguments
    - `r̄::UnitfulCoordinateSystems.Coordinate`: spatial location of the observation point
    - `t::Unitful.Time`: time at which the electric field is observed
    - `source::JefimenkoSource`: source of the electric field
    - `media::PropagationMedia`: properties of the propagation media

    # Keywords
    - `rtol::Real`: relative tolerance at which to solve the integral (optional)
    """
    function _𝐏(r̄::Coordinate, t::Unitful.Time, source::JefimenkoSource{T},
                media::PropagationMedia; rtol=__DEFAULT_RTOL) where {T<:AbstractFloat}
        E = _𝐄(r̄,t,source,media; rtol=rtol)
        H = _𝐇(r̄,t,source,media; rtol=rtol)
        return cross(E,H) .|> W/m^2
    end

    include("integrands_E.jl")
    include("fields_E.jl")

    include("integrands_H.jl")
    include("fields_H.jl")

    export 𝐄, 𝐇, 𝐏
end
