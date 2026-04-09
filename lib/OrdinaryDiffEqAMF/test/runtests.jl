using Test

@testset "OrdinaryDiffEqAMF" begin
    include("test_pollu.jl")
    include("test_fd2d.jl")
end
