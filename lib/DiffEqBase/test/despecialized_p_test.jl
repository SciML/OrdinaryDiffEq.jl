using DiffEqBase
using SciMLBase
using SymbolicIndexingInterface: SymbolCache
using Test

struct DynamicP
    rate::Float64
end

struct OtherDynamicP
    rate::Float64
    unused::Int
end

struct DynamicSDEParameters
    rate::Float64
    sigma::Float64
end

struct DynamicDAEResidual{ID} end
struct DynamicInitializationResidual{ID} end

function (::DynamicDAEResidual)(resid, du, u, p, t)
    seen_dae_parameter[] = typeof(p)
    resid[1] = du[1] + p.rate * u[1]
    return nothing
end

function (::DynamicInitializationResidual)(resid, u, p)
    resid[1] = u[1]
    return nothing
end

function dynamic_initialization_data(::Val{ID}, p) where {ID}
    initprob = NonlinearProblem{true}(DynamicInitializationResidual{ID}(), [0.0], p)
    return SciMLBase.OverrideInitData(initprob, nothing, nothing, nothing)
end

const seen_rhs_parameter = Ref{DataType}()
const seen_jac_parameter = Ref{DataType}()
const seen_tgrad_parameter = Ref{DataType}()
const seen_dae_parameter = Ref{DataType}()
const seen_dde_parameter = Ref{DataType}()
const seen_sde_drift_parameter = Ref{DataType}()
const seen_sde_noise_parameter = Ref{DataType}()
const seen_sdde_drift_parameter = Ref{DataType}()
const seen_sdde_noise_parameter = Ref{DataType}()

function dynamic_rhs!(du, u, p, t)
    seen_rhs_parameter[] = typeof(p)
    du[1] = -p.rate * u[1]
    return nothing
end

function dynamic_jac!(J, u, p, t)
    seen_jac_parameter[] = typeof(p)
    J[1, 1] = -p.rate
    return nothing
end

function dynamic_tgrad!(dT, u, p, t)
    seen_tgrad_parameter[] = typeof(p)
    dT[1] = 0
    return nothing
end

function dynamic_sde_drift!(du, u, p, t)
    seen_sde_drift_parameter[] = typeof(p)
    du[1] = -p.rate * u[1]
    return nothing
end

function dynamic_sde_noise!(du, u, p, t)
    seen_sde_noise_parameter[] = typeof(p)
    du[1] = p.sigma * u[1]
    return nothing
end

function dynamic_sdde_drift!(du, u, h, p, t)
    seen_sdde_drift_parameter[] = typeof(p)
    du[1] = -p.rate * u[1]
    return nothing
end

function dynamic_sdde_noise!(du, u, h, p, t)
    seen_sdde_noise_parameter[] = typeof(p)
    du[1] = p.sigma * u[1]
    return nothing
end

function dynamic_problem(p; sys = nothing)
    f = ODEFunction{true, SciMLBase.AutoDespecialize}(
        dynamic_rhs!; jac = dynamic_jac!, tgrad = dynamic_tgrad!, sys
    )
    return ODEProblem(f, [1.0], (0.0, 1.0), p)
end

concretize(prob, alg = nothing) = DiffEqBase.get_concrete_problem(prob, true; alg)

@testset "AutoDespecialize uses the dynamic parameter barrier" begin
    first_prob = concretize(dynamic_problem(DynamicP(0.5)))
    second_prob = concretize(dynamic_problem(OtherDynamicP(1.5, 1)))

    @test first_prob.p isa SciMLBase.DespecializedParameters
    @test SciMLBase.unwrap_parameters(first_prob.p) isa DynamicP
    @test SciMLBase.unwrap_parameters(second_prob.p) isa OtherDynamicP
    @test typeof(first_prob.p) === typeof(second_prob.p)
    @test typeof(first_prob.f) === typeof(second_prob.f)
    @test typeof(first_prob) === typeof(second_prob)
    @test SciMLBase.unwrapped_f(first_prob.f).f === dynamic_rhs!

    du = zeros(1)
    first_prob.f(du, first_prob.u0, first_prob.p, 0.0)
    @test du == [-0.5]
    @test seen_rhs_parameter[] === DynamicP

    barrier = DiffEqBase.ParameterDespecializationWrapper(dynamic_rhs!)
    barrier(du, first_prob.u0, SciMLBase.unwrap_parameters(first_prob.p), 0.0)
    @test du == [-0.5]
    @test seen_rhs_parameter[] === DynamicP

    J = zeros(1, 1)
    first_prob.f.jac(J, first_prob.u0, first_prob.p, 0.0)
    @test J == [-0.5;;]
    @test seen_jac_parameter[] === DynamicP

    dT = zeros(1)
    first_prob.f.tgrad(dT, first_prob.u0, first_prob.p, 0.0)
    @test dT == [0.0]
    @test seen_tgrad_parameter[] === DynamicP

    remade = concretize(first_prob)
    @test typeof(remade) === typeof(first_prob)
    @test remade.p === first_prob.p
end

@testset "other specialization policies are unchanged" begin
    p = DynamicP(0.5)
    auto = ODEProblem{true, SciMLBase.AutoSpecialize}(
        dynamic_rhs!, [1.0], (0.0, 1.0), p
    )
    respecialized = ODEProblem{true, SciMLBase.AutoRespecialize}(
        dynamic_rhs!, [1.0], (0.0, 1.0), p
    )
    full = ODEProblem{true, SciMLBase.FullSpecialize}(
        dynamic_rhs!, [1.0], (0.0, 1.0), p
    )
    @test concretize(auto).p === p
    @test respecialized.p === p
    @test concretize(full).p === p
end

@testset "symbolic systems retain the dynamic path" begin
    sys = SymbolCache([:x], [:rate], :t)
    prob = concretize(dynamic_problem(DynamicP(0.5); sys))
    @test prob.p isa SciMLBase.DespecializedParameters
    du = zeros(1)
    seen_rhs_parameter[] = Nothing
    prob.f(du, prob.u0, prob.p, 0.0)
    @test du == [-0.5]
    @test seen_rhs_parameter[] === DynamicP
end

@testset "other differential equation functions use the dynamic path" begin
    p = DynamicP(0.5)

    dae_initdata = dynamic_initialization_data(Val(1), p)
    dae_f = DAEFunction{true, SciMLBase.AutoDespecialize}(
        DynamicDAEResidual{1}(); initialization_data = dae_initdata
    )
    dae_prob = DAEProblem(dae_f, [0.0], [1.0], (0.0, 1.0), p)
    dae_concrete = concretize(dae_prob)
    other_p = OtherDynamicP(1.5, 1)
    other_dae_initdata = dynamic_initialization_data(Val(2), other_p)
    other_dae_f = DAEFunction{true, SciMLBase.AutoDespecialize}(
        DynamicDAEResidual{2}(); initialization_data = other_dae_initdata
    )
    other_dae_prob = DAEProblem(
        other_dae_f, [0.0], [1.0], (0.0, 1.0), other_p
    )
    other_dae_concrete = concretize(other_dae_prob)
    @test dae_concrete.p isa SciMLBase.DespecializedParameters
    @test typeof(dae_concrete.f) === typeof(other_dae_concrete.f)
    @test typeof(dae_concrete) === typeof(other_dae_concrete)
    @test dae_concrete.f.initialization_data === dae_initdata
    @test other_dae_concrete.f.initialization_data === other_dae_initdata
    @test SciMLBase.unwrapped_f(dae_concrete.f.f) isa DynamicDAEResidual{1}
    resid = zeros(1)
    dae_concrete.f(resid, [0.0], [1.0], dae_concrete.p, 0.0)
    @test resid == [0.5]
    @test seen_dae_parameter[] === DynamicP

    function dde_rhs!(du, u, h, p, t)
        seen_dde_parameter[] = typeof(p)
        du[1] = -p.rate * u[1]
        return nothing
    end
    history(p, t) = [1.0]
    dde_f = DDEFunction{true, SciMLBase.AutoDespecialize}(dde_rhs!)
    dde_prob = DDEProblem(dde_f, [1.0], history, (0.0, 1.0), p)
    dde_concrete = concretize(dde_prob)
    @test dde_concrete.p isa SciMLBase.DespecializedParameters
    du = zeros(1)
    dde_concrete.f(du, [1.0], history, dde_concrete.p, 0.0)
    @test du == [-0.5]
    @test seen_dde_parameter[] === DynamicP

    sde_p = DynamicSDEParameters(0.5, 0.1)
    sde_f = SDEFunction{true, SciMLBase.AutoDespecialize}(
        dynamic_sde_drift!, dynamic_sde_noise!
    )
    sde_prob = SDEProblem(sde_f, [1.0], (0.0, 1.0), sde_p)
    sde_concrete = concretize(sde_prob)
    @test sde_concrete.p isa SciMLBase.DespecializedParameters
    sde_concrete.f(du, [1.0], sde_concrete.p, 0.0)
    @test du == [-0.5]
    @test seen_sde_drift_parameter[] === DynamicSDEParameters
    sde_concrete.f.g(du, [1.0], sde_concrete.p, 0.0)
    @test du == [0.1]
    @test seen_sde_noise_parameter[] === DynamicSDEParameters
    @test SciMLBase.unwrapped_f(sde_concrete.f).g === dynamic_sde_noise!

    sdde_f = SDDEFunction{true, SciMLBase.AutoDespecialize}(
        dynamic_sdde_drift!, dynamic_sdde_noise!
    )
    sdde_prob = SDDEProblem(sdde_f, dynamic_sde_noise!, [1.0], history, (0.0, 1.0), sde_p)
    sdde_concrete = concretize(sdde_prob)
    @test sdde_concrete.p isa SciMLBase.DespecializedParameters
    sdde_concrete.f(du, [1.0], history, sdde_concrete.p, 0.0)
    @test du == [-0.5]
    @test seen_sdde_drift_parameter[] === DynamicSDEParameters
    sdde_concrete.f.g(du, [1.0], history, sdde_concrete.p, 0.0)
    @test du == [0.1]
    @test seen_sdde_noise_parameter[] === DynamicSDEParameters
    @test SciMLBase.unwrapped_f(sdde_concrete.f.g) === dynamic_sdde_noise!
end
