using MonteCarloCollision
using Test

@testset "MonteCarloCollision.jl" begin
    @test include("energy_conservation.jl")
end
