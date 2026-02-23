using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqBDF
using OrdinaryDiffEqFIRK
using OrdinaryDiffEqNonlinearSolve
using CUDA
using LinearAlgebra
using Adapt
using SparseArrays
using Test
using CUDSS
using Printf

#=
du[1] = -u[1]
du[2] = -0.5*u[2]
    0 =  u[1] + u[2] - u[3]
    0 = -u[1] + u[2] - u[4]
=#

function dae!(du, u, p, t)
    return mul!(du, p, u)
end

p = [
    -1 0 0 0
    1 -0.5 0 0
    1 1 -1 0
    -1 1 0 -1
]

mass_matrix = Diagonal([1, 1, 0, 0])
jac_prototype = sparse(map(x -> iszero(x) ? 0.0 : 1.0, p))

u0 = [1.0, 1.0, 0.5, 0.5] # force init
tspan = (0.0, 5.0)

# CPU reference solution
odef = ODEFunction(dae!, mass_matrix = mass_matrix, jac_prototype = jac_prototype)
prob = ODEProblem(odef, u0, tspan, p)
sol = solve(prob, Rodas5P())

# GPU data -- we use F64 for higher accuracy for comparision
u0_d = adapt(CuArray{Float64}, u0)
p_d = adapt(CuArray{Float64}, p)

# dense or spares mass matrix does not work yet!
mass_matrix_d = cu(mass_matrix)

function test_gpu_dae(jac_prototype_d, solver)
    sol_ref = solve(prob, solver)
    odef_d = ODEFunction(dae!, mass_matrix = mass_matrix_d, jac_prototype = jac_prototype_d)
    prob_d = ODEProblem(odef_d, u0_d, tspan, p_d)
    sol_d = solve(prob_d, solver)

    for t in sol_d.t
        u = Vector(sol_d(t))
        @test isapprox(u[1] + u[2], u[3]; atol = 1.0e-6)
        @test isapprox(-u[1] + u[2], u[4]; atol = 1.0e-6)
    end

    max_abserr = 0.0
    max_relerr = 0.0
    for t in tspan[begin]:0.1:tspan[end]
        diff = abs.(Vector(sol_d(t)) - sol_ref(t))
        ref = abs.(sol_ref(t))
        max_abserr = max(max_abserr, maximum(diff))
        max_relerr = max(max_relerr, maximum(diff ./ max.(ref, eps())))
    end
    cond = (max_abserr < 1e-5 || max_relerr < 1e-5)
    cond || println("Solver $sol abstol=$max_abserr retlol=$max_relerr")
    @test (max_abserr < 1e-5 || max_relerr < 1e-5)
end

# Jacobian prototype options
jac_prots = [
    "none" => nothing,
    "CSC" => CUDA.CUSPARSE.CuSparseMatrixCSC(jac_prototype),
    "CSR" => CUDA.CUSPARSE.CuSparseMatrixCSR(jac_prototype),
]

# Solver options
solvers = [
    # Rosenbrock (diagonal mass matrix only)
    "Rosenbrock23" => Rosenbrock23(),           # ✅
    "Rosenbrock32" => Rosenbrock32(),           # 📉 poor fit for this DAE, 💥 CSC catastrophic
    # Rosenbrock 2nd order
    "ROS2" => ROS2(),                           # ✅
    "ROS2PR" => ROS2PR(),                       # ✅
    "ROS2S" => ROS2S(),                         # ✅
    # Rosenbrock 3rd order
    "ROS3" => ROS3(),                           # ✅
    "ROS3PR" => ROS3PR(),                       # 📉 poor fit for this DAE
    "ROS3PRL" => ROS3PRL(),                     # ⚠️ CSC accuracy loss
    "ROS3PRL2" => ROS3PRL2(),                   # ⚠️ CSC accuracy loss
    "ROS3P" => ROS3P(),                         # 📉 poor fit for this DAE
    "Rodas3" => Rodas3(),                       # ⚠️ CSC accuracy loss
    # "Rodas23W" => Rodas23W(),                   # 🚫 scalar indexing, requires large cahnges to `calculate_interpoldiff!`
    # "Rodas3P" => Rodas3P(),                     # 🚫 scalar indexing
    "Scholz4_7" => Scholz4_7(),                 # ✅
    # Rosenbrock 4th order
    "ROS34PW1a" => ROS34PW1a(),                 # 📉 poor fit for this DAE
    "ROS34PW1b" => ROS34PW1b(),                 # 📉 poor fit for this DAE
    "ROS34PW2" => ROS34PW2(),                   # ⚠️ CSC accuracy loss
    "ROS34PW3" => ROS34PW3(),                   # ✅
    "ROS34PRw" => ROS34PRw(),                   # ✅
    "RosShamp4" => RosShamp4(),                 # ✅
    "Veldd4" => Veldd4(),                       # ✅
    "Velds4" => Velds4(),                       # ✅
    "GRK4T" => GRK4T(),                         # ✅
    "GRK4A" => GRK4A(),                         # ✅
    "Ros4LStab" => Ros4LStab(),                 # ✅
    "Rodas4" => Rodas4(),                       # ✅
    "Rodas42" => Rodas42(),                     # ✅
    "Rodas4P" => Rodas4P(),                     # ✅
    "Rodas4P2" => Rodas4P2(),                   # ✅
    "ROK4a" => ROK4a(),                         # ✅
    # Rosenbrock 5th order
    "Rodas5" => Rodas5(),                       # ✅
    "Rodas5P" => Rodas5P(),                     # ✅
    "Rodas5Pe" => Rodas5Pe(),                   # ✅
    "Rodas5Pr" => Rodas5Pr(),                   # ✅
    # Rosenbrock 6th order
    "Rodas6P" => Rodas6P(),                     # ✅
    # SDIRK (don't include fixed time step which need explicit dt)
    "ImplicitEuler" => ImplicitEuler(),         # ✅
    "Trapezoid" => Trapezoid(),                 # ✅
    "SDIRK2" => SDIRK2(),                       # ✅
    "Cash4" => Cash4(),                         # ⚠️ CSC accuracy loss
    "Hairer4" => Hairer4(),                     # ⚠️ CSC accuracy loss
    "Hairer42" => Hairer42(),                   # ⚠️ CSC accuracy loss
    # BDF
    "ABDF2" => ABDF2(),                         # ⚠️ CSC accuracy loss (💥 catastrophic)
    "QNDF1" => QNDF1(),                        # ✅
    "QNDF2" => QNDF2(),                        # ⚠️ CSC accuracy loss
    # "QNDF" => QNDF(),                          # 🔧 DeviceMemory error in LinAlg
    "QBDF1" => QBDF1(),                        # ✅
    "QBDF2" => QBDF2(),                        # ⚠️ accuracy loss (all jac_prots)
    # "QBDF" => QBDF(),                          # 🔧 DeviceMemory error in LinAlg
    # "FBDF" => FBDF(),                           # 🚫 scalar indexing -> needs extensive work on reinitFBDF!
    # FIRK -> all need subtential changes to `perform_step!` for FIRK methods
    # "RadauIIA3" => RadauIIA3(),                 # 🚫 scalar indexing, ComplexF64 sparse unsupported
    # "RadauIIA5" => RadauIIA5(),                 # 🚫 scalar indexing, ComplexF64 sparse unsupported
    # "RadauIIA9" => RadauIIA9(),                 # 🚫 scalar indexing, ComplexF64 sparse unsupported
    # "AdaptiveRadau" => AdaptiveRadau(),         # 🚫 scalar indexing, ComplexF64 sparse unsupported
]

function _maxerrs(sol_a, sol_b)
    max_abs = 0.0
    max_rel = 0.0
    for t in tspan[begin]:0.1:tspan[end]
        a = isa(sol_a(t), CuArray) ? Vector(sol_a(t)) : Vector(sol_a(t))
        b = isa(sol_b(t), CuArray) ? Vector(sol_b(t)) : Vector(sol_b(t))
        diff = abs.(a - b)
        ref = abs.(b)
        max_abs = max(max_abs, maximum(diff))
        max_rel = max(max_rel, maximum(diff ./ max.(ref, eps())))
    end
    return max_abs, max_rel
end

function _printval(val; threshold = 1e-3)
    s = @sprintf("%.2e", val)
    if val > threshold
        printstyled(s; color = :red)
    else
        print(s)
    end
end

function debug_gpu_dae(jac_prototype_d, solver, name)
    # CPU: this solver vs Rodas5P reference
    sol_cpu = solve(prob, solver)
    cpu_abs, cpu_rel = _maxerrs(sol_cpu, sol)

    # GPU solve
    odef_d = ODEFunction(dae!, mass_matrix = mass_matrix_d, jac_prototype = jac_prototype_d)
    prob_d = ODEProblem(odef_d, u0_d, tspan, p_d)
    sol_d = solve(prob_d, solver)

    # GPU vs CPU same solver
    gpu_abs, gpu_rel = _maxerrs(sol_d, sol_cpu)

    # Print row
    print(rpad(name, 30))
    _printval(cpu_abs); print("  ")
    _printval(cpu_rel); print("  ")
    _printval(gpu_abs); print("  ")
    _printval(gpu_rel)
    println()
end

using Printf

println(rpad("", 30), "cpu_abs  cpu_rel  gpu_abs  gpu_rel")
println("-"^72)

for (sn, sv) in solvers, (jn, jp) in jac_prots
    label = "$sn / $jn"
    try
        debug_gpu_dae(jp, sv, label)
    catch e
        printstyled(rpad(label, 30), "ERROR: ", sprint(showerror, e; context=:limit=>true), "\n"; color = :yellow)
    end
end

println("-"^72)

#= Uncomment for actual testing
@testset "End-to-end GPU compat of mass matrix DAE solvers" begin
    for (sn, sv) in solvers
        @testset "GPU DAE: $sn" begin
            for (jn, jp) in jac_prots
                @testset "Jacboain prototype: $jn" begin
                    test_gpu_dae(jp, sv)
                end
            end
        end
    end
end
=#
