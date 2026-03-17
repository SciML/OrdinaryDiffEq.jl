using StochasticDiffEq, DiffEqCallbacks, DiffEqNoiseProcess, Test, Random
import SciMLBase

Random.seed!(100)
# Definitions according to vector/matrix representations of operators from QuantumOptics
u0 = ComplexF64[0.0, 1.0]
A = zeros(ComplexF64, 2, 2)
A[2, 1] = 1.0
A[1, 2] = 1.0
As = zeros(ComplexF64, 2, 2)
As[1, 1] = 0.5
As[2, 2] = -0.5

# SDE functions
function f_qm(du, u, p, t)
    return du .= -1.0im * A * u
end
function g_qm(du, u, p, t)
    return du .= -1.0im * As * u
end

# Set up two identical instances of saving
T = [0:0.1:2;]
fout(u, t, integrator) = abs2(u[1]) - abs2(u[2])

# Define the problems -- note the difference in noise
prob1 = SDEProblem{true}(
    f_qm, g_qm, u0, (T[1], T[end]),
    noise = StochasticDiffEq.RealWienerProcess(
        0.0, 0.0, 0.0,
        rswm = DiffEqNoiseProcess.RSWM(adaptivealg = :RSwM1)
    )
)

prob2 = SDEProblem{true}(
    f_qm, g_qm, u0, (T[1], T[end]),
    noise = StochasticDiffEq.RealWienerProcess(
        0.0, 0.0, 0.0,
        rswm = DiffEqNoiseProcess.RSWM(adaptivealg = :RSwM3)
    )
)

# Simple averaging
Ntraj = 1000
avg1 = zeros(size(T)...)
avg2 = zeros(size(T)...)
avg3 = zeros(size(T)...)

let
    global avg1, avg2, avg3
    for i in 1:Ntraj
        out1 = SavedValues(Float64, ComplexF64)
        scb1 = SavingCallback(
            fout, out1, saveat = T, save_everystep = false, save_start = false
        )

        solve(
            prob1, RKMil{SciMLBase.AlgorithmInterpretation.Stratonovich}(),
            dt = 1.0e-4, callback = scb1, seed = i, adaptive = false
        )
        avg1 .+= out1.saveval ./ Ntraj

        out1 = SavedValues(Float64, ComplexF64)
        scb1 = SavingCallback(
            fout, out1, saveat = T, save_everystep = false, save_start = false
        )

        solve(
            prob1, RKMil{SciMLBase.AlgorithmInterpretation.Stratonovich}(),
            tstops = T, callback = scb1,
            save_everystep = false, save_start = false
        )
        avg2 .+= out1.saveval ./ Ntraj

        out1 = SavedValues(Float64, ComplexF64)
        scb1 = SavingCallback(
            fout, out1, saveat = T, save_everystep = false, save_start = false
        )

        solve(
            prob2, RKMil{SciMLBase.AlgorithmInterpretation.Stratonovich}(),
            tstops = T, callback = scb1,
            save_everystep = false, save_start = false
        )
        avg3 .+= out1.saveval ./ Ntraj
    end
end

#=
using Plots; gr()
plot(T, avg1)
plot!(T, avg2)
plot!(T, avg3)
=#

@test maximum(avg1 - avg2) < 0.02
@test maximum(avg1 - avg3) < 0.03

#Adding Similar Tests for RKMilGeneral
# Simple averaging
Ntraj = 1000
avg1 = zeros(size(T)...)
avg2 = zeros(size(T)...)
avg3 = zeros(size(T)...)

let
    global avg1, avg2, avg3
    for i in 1:Ntraj
        out1 = SavedValues(Float64, ComplexF64)
        scb1 = SavingCallback(
            fout, out1, saveat = T, save_everystep = false, save_start = false
        )

        solve(
            prob1,
            RKMilGeneral(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
            dt = 1.0e-4,
            callback = scb1, seed = i, adaptive = false
        )
        avg1 .+= out1.saveval ./ Ntraj

        out1 = SavedValues(Float64, ComplexF64)
        scb1 = SavingCallback(
            fout, out1, saveat = T, save_everystep = false, save_start = false
        )

        solve(
            prob1,
            RKMilGeneral(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
            tstops = T, callback = scb1,
            save_everystep = false, save_start = false
        )
        avg2 .+= out1.saveval ./ Ntraj

        out1 = SavedValues(Float64, ComplexF64)
        scb1 = SavingCallback(
            fout, out1, saveat = T, save_everystep = false, save_start = false
        )

        solve(
            prob2,
            RKMilGeneral(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich),
            tstops = T, callback = scb1,
            save_everystep = false, save_start = false
        )
        avg3 .+= out1.saveval ./ Ntraj
    end
end

#=
using Plots; gr()
plot(T, avg1)
plot!(T, avg2)
plot!(T, avg3)
=#

@test maximum(avg1 - avg2) < 0.03
@test maximum(avg1 - avg3) < 0.04
