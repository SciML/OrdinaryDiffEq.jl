using OrdinaryDiffEqExplicitTableaus
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(OrdinaryDiffEqExplicitTableaus)
end

@testset "JET" begin
    JET.test_package(OrdinaryDiffEqExplicitTableaus; target_defined_modules = true)
end
