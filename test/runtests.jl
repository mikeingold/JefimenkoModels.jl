using Unitful, UnitfulCoordinateSystems
using Unitful.DefaultSymbols: m, ns
using JefimenkoModels, Test

@testset "Accessors for LineSource_Straight" begin
    a = CoordinateCartesian(-0.5m, 0.0m, 0.0m)
    b = CoordinateCartesian( 0.5m, 0.0m, 0.0m)
    rho(r̄::AbstractCoordinate, t_s::Real) = cos(2π*100e6*t_s)^2
    J(r̄::AbstractCoordinate, t_s::Real) = x̂ .* cos(2π*100e6*t_s)
    source = LineSource_Straight{Float64}(a, b, rho, rho, J, J)

	@test source.a.x ≈ source.ā.x ≈ -0.5m
    @test source.b.x ≈ source.b̄.x ≈  0.5m
    @test source.rho_e === source.ρₑ === source.ρe
    @test source.rho_h === source.ρₕ === source.ρh
    @test source.J_e === source.Je === source.Jₑ
    @test source.J_h === source.Jh === source.Jₕ
end
