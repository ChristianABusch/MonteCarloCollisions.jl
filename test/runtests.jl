using MonteCarloCollisions
using Test

@testset "MonteCarloCollision.jl" begin
    @test include("energy_conservation/isotropic.jl")
    @test include("energy_conservation/backscatter.jl")
    @test include("energy_conservation/excitation.jl")
    @test include("energy_conservation/ionization.jl")
end
