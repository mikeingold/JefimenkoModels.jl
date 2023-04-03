################################################################################
#                              ACCESSOR FUNCTIONS
################################################################################

function Base.getproperty(s::VolumeSource_Rectangular, sym::Symbol)
    if sym in (:xlims, :ylims, :zlims)  # included
        return getfield(s, sym)
    elseif sym in (:rho_e, :rho_h, :J_e, :J_h)  # included
        return getfield(s, sym)
    elseif sym in (:ρₑ, :ρe)  # aliases
        return getfield(s, :rho_e)
    elseif sym in (:ρₕ, :ρh)  # aliases
        return getfield(s, :rho_h)
    elseif sym in (:Jₑ, :Je)  # aliases
        return getfield(s, :J_e)
    elseif sym in (:Jₕ, :Jh)  # aliases
        return getfield(s, :J_h)
    else  # fallback
        return getfield(s, sym)
    end
end

function Base.getproperty(s::VolumeSource_Cylinder, sym::Symbol)
    if sym in (:r, :philims, :zlims)  # included
        return getfield(s, sym)
    elseif sym in (:rho_e, :rho_h, :J_e, :J_h)  # included
        return getfield(s, sym)
    elseif sym in (:ρₑ, :ρe)  # aliases
        return getfield(s, :rho_e)
    elseif sym in (:ρₕ, :ρh)  # aliases
        return getfield(s, :rho_h)
    elseif sym in (:Jₑ, :Je)  # aliases
        return getfield(s, :J_e)
    elseif sym in (:Jₕ, :Jh)  # aliases
        return getfield(s, :J_h)
    else  # fallback
        return getfield(s, sym)
    end
end

function Base.getproperty(s::VolumeSource_Sphere, sym::Symbol)
    if sym in (:r, :thetalims, :philims)  # included
        return getfield(s, sym)
    elseif sym in (:rho_e, :rho_h, :J_e, :J_h)  # included
        return getfield(s, sym)
    elseif sym in (:ρₑ, :ρe)  # aliases
        return getfield(s, :rho_e)
    elseif sym in (:ρₕ, :ρh)  # aliases
        return getfield(s, :rho_h)
    elseif sym in (:Jₑ, :Je)  # aliases
        return getfield(s, :J_e)
    elseif sym in (:Jₕ, :Jh)  # aliases
        return getfield(s, :J_h)
    else  # fallback
        return getfield(s, sym)
    end
end

function Base.getproperty(s::SurfaceSource_Rectangle, sym::Symbol)
    if sym in (:xlims, :ylims)  # included
        return getfield(s, sym)
    elseif sym in (:rho_e, :rho_h, :J_e, :J_h)  # included
        return getfield(s, sym)
    elseif sym in (:ρₑ, :ρe)  # aliases
        return getfield(s, :rho_e)
    elseif sym in (:ρₕ, :ρh)  # aliases
        return getfield(s, :rho_h)
    elseif sym in (:Jₑ, :Je)  # aliases
        return getfield(s, :J_e)
    elseif sym in (:Jₕ, :Jh)  # aliases
        return getfield(s, :J_h)
    else  # fallback
        return getfield(s, sym)
    end
end

function Base.getproperty(s::SurfaceSource_Disk, sym::Symbol)
    if sym in (:r)  # included
        return getfield(s, sym)
    elseif sym in (:rho_e, :rho_h, :J_e, :J_h)  # included
        return getfield(s, sym)
    elseif sym in (:ρ₀, :ρ_0, :rho_0)  # aliases
        return getfield(s, :r)
    elseif sym in (:ρₑ, :ρe)  # aliases
        return getfield(s, :rho_e)
    elseif sym in (:ρₕ, :ρh)  # aliases
        return getfield(s, :rho_h)
    elseif sym in (:Jₑ, :Je)  # aliases
        return getfield(s, :J_e)
    elseif sym in (:Jₕ, :Jh)  # aliases
        return getfield(s, :J_h)
    else  # fallback
        return getfield(s, sym)
    end
end

function Base.getproperty(s::LineSource_Straight, sym::Symbol)
    if sym in (:a, :b)  # included
        return getfield(s, sym)
    elseif sym in (:rho_e, :rho_h, :J_e, :J_h)  # included
        return getfield(s, sym)
    elseif sym in (:ā)  # aliases
        return getfield(s, :a)
    elseif sym in (:b̄)  # aliases
        return getfield(s, :b)
    elseif sym in (:ρₑ, :ρe)  # aliases
        return getfield(s, :rho_e)
    elseif sym in (:ρₕ, :ρh)  # aliases
        return getfield(s, :rho_h)
    elseif sym in (:Jₑ, :Je)  # aliases
        return getfield(s, :J_e)
    elseif sym in (:Jₕ, :Jh)  # aliases
        return getfield(s, :J_h)
    else  # fallback
        return getfield(s, sym)
    end
end

function Base.getproperty(m::AbstractPropagationMedia, sym::Symbol)
    if sym in (:epsilon, :mu, :c)  # included
        return getfield(s, sym)
    elseif sym in (:ε, :ϵ)  # aliases
        return getfield(s, :epsilon)
    elseif sym in (:μ)  # aliases
        return getfield(s, :mu)
    else  # fallback
        return getfield(s, sym)
    end
end
