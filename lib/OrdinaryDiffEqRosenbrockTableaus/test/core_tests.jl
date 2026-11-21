using OrdinaryDiffEqRosenbrockTableaus
using OrdinaryDiffEqRosenbrock: RodasTableau
using Test

const RT = OrdinaryDiffEqRosenbrockTableaus

@testset "OrdinaryDiffEqRosenbrockTableaus" begin
    @testset "Tableau construction" begin
        tab = RT.Rodas4Tableau(Float64, Float64)
        @test tab isa RodasTableau
        @test tab.gamma isa Float64

        tab = RT.Rodas42Tableau(Float64, Float64)
        @test tab isa RodasTableau

        tab = RT.Rodas4PTableau(Float64, Float64)
        @test tab isa RodasTableau

        tab = RT.Rodas4P2Tableau(Float64, Float64)
        @test tab isa RodasTableau

        tab = RT.Rodas5Tableau(Float64, Float64)
        @test tab isa RodasTableau

        tab = RT.ROS3PRodasTableau(Float64, Float64)
        @test tab isa RodasTableau

        tab = RT.Rodas3RodasTableau(Float64, Float64)
        @test tab isa RodasTableau

        tab = RT.Rodas3PRodasTableau(Float64, Float64)
        @test tab isa RodasTableau

        tab = RT.RosShamp4RodasTableau(Float64, Float64)
        @test tab isa RodasTableau

        tab = RT.ROS34PW2RodasTableau(Float64, Float64)
        @test tab isa RodasTableau

        tab = RT.RosenbrockW6S4OSRodasTableau(Float64, Float64)
        @test tab isa RodasTableau
    end

    @testset "BigFloat construction" begin
        tab = RT.Rodas4Tableau(BigFloat, BigFloat)
        @test tab.gamma isa BigFloat
        @test eltype(tab.A) == BigFloat
    end
end

@time @testset "ODE Rosenbrock Convergence Tests" begin
    include("ode_rosenbrock_tests.jl")
end
