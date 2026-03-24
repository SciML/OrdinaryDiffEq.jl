using DiffEqBase, ForwardDiff, Test, InteractiveUtils
using ReverseDiff, SciMLStructures

u0 = 2.0
p = 2.0
t0 = 1.0

@test DiffEqBase.promote_u0(u0, p, t0) isa Float64
@test DiffEqBase.promote_u0(u0, p, t0) == 2.0
@test DiffEqBase.promote_u0(cis(u0), p, t0) isa ComplexF64
@test DiffEqBase.promote_u0(cis(u0), p, t0) == cis(2.0)

struct MyStruct{T, T2} <: Number
    x::T
    y::T2
end

struct MyStruct2{T, T2}
    x::T
    y::T2
    MyStruct2(x) = new{typeof(x), Any}(x)
end

struct MyStruct3{T, T2}
    x::T
    y::T2
    MyStruct3(x) = new{typeof(x), Float64}(x)
end

module Mod end

p_possibilities = [
    ForwardDiff.Dual(2.0), (ForwardDiff.Dual(2.0), 2.0),
    [ForwardDiff.Dual(2.0)], ([ForwardDiff.Dual(2.0)], 2.0),
    (2.0, ForwardDiff.Dual(2.0)), (; x = 2.0, y = ForwardDiff.Dual(2.0)),
    (; x = 2.0, y = [ForwardDiff.Dual(2.0)]), (; x = 2.0, y = [[ForwardDiff.Dual(2.0)]]),
    Set([2.0, ForwardDiff.Dual(2.0)]), (SciMLBase.NullParameters(), ForwardDiff.Dual(2.0)),
    ((), ForwardDiff.Dual(2.0)), ForwardDiff.Dual{Nothing}(ForwardDiff.Dual{MyStruct}(2.0)),
    (ForwardDiff.Dual(2.0)), [(1.0, ForwardDiff.Dual(1.0, (1.0,)))],
]

for p in p_possibilities
    @show p
    @test DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual
    local u0 = 2.0
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @test DiffEqBase.promote_u0([cis(u0)], p, t0) isa
        AbstractArray{<:Complex{<:ForwardDiff.Dual}}
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @test DiffEqBase.promote_u0([cis(u0)], p, t0) isa
        AbstractArray{<:Complex{<:ForwardDiff.Dual}}
    @inferred DiffEqBase.anyeltypedual(p)
end

higher_order_p_possibilities = [
    ForwardDiff.Dual{Nothing}(ForwardDiff.Dual{MyStruct}(2.0)),
    (
        ForwardDiff.Dual{Nothing}(ForwardDiff.Dual{MyStruct}(2.0)),
        SciMLBase.NullParameters(),
    ),
    (
        ForwardDiff.Dual{Nothing}(ForwardDiff.Dual{MyStruct}(2.0)),
        ForwardDiff.Dual{Nothing}(2.0),
    ),
    (
        ForwardDiff.Dual{Nothing}(2.0),
        ForwardDiff.Dual{Nothing}(ForwardDiff.Dual{MyStruct}(2.0)),
    ),
]

for p in higher_order_p_possibilities
    @show p
    @test DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual
    @test DiffEqBase.anyeltypedual(p) <:
    ForwardDiff.Dual{Nothing, ForwardDiff.Dual{MyStruct, Float64, 0}, 0}
    local u0 = 2.0
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @test DiffEqBase.promote_u0([cis(u0)], p, t0) isa
        AbstractArray{<:Complex{<:ForwardDiff.Dual}}
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @test DiffEqBase.promote_u0([cis(u0)], p, t0) isa
        AbstractArray{<:Complex{<:ForwardDiff.Dual}}
    @inferred DiffEqBase.anyeltypedual(p)
end

p_possibilities17 = [
    MyStruct(2.0, ForwardDiff.Dual(2.0)), [MyStruct(2.0, ForwardDiff.Dual(2.0))],
    [MyStruct(2.0, [2.0, ForwardDiff.Dual(2.0)])],
    [MyStruct(2.0, (2.0, ForwardDiff.Dual(2.0)))],
    ((;), ForwardDiff.Dual(2.0)), MyStruct3(ForwardDiff.Dual(2.0)),
    (Mod, ForwardDiff.Dual(2.0)), (() -> 2.0, ForwardDiff.Dual(2.0)),
    (Base.pointer([2.0]), ForwardDiff.Dual(2.0)),
]
push!(p_possibilities17, Returns((a = 2, b = 1.3, c = ForwardDiff.Dual(2.0f0))))

for p in p_possibilities17
    @show p
    @test DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual
    local u0 = 2.0
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @test DiffEqBase.promote_u0([cis(u0)], p, t0) isa
        AbstractArray{<:Complex{<:ForwardDiff.Dual}}
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @test DiffEqBase.promote_u0([cis(u0)], p, t0) isa
        AbstractArray{<:Complex{<:ForwardDiff.Dual}}

    @inferred DiffEqBase.anyeltypedual(p)
    ci = InteractiveUtils.@code_typed DiffEqBase.anyeltypedual(p)
    @show filter(!=(Expr(:code_coverage_effect)), ci.first.code)
    #@test count(x -> (x != (Expr(:code_coverage_effect))) &&
    #                (x != GlobalRef(DiffEqBase, :Any)), ci.first.code) == 1
end

p_possibilities_uninferrred = [
    Dict(:x => 2.0, :y => ForwardDiff.Dual(2.0)),
    Dict(:x => 2.0, :y => [ForwardDiff.Dual(2.0)]),
    Dict(:x => 2.0, :y => [(; x = (ForwardDiff.Dual(2.0), 2.0), y = 2.0)]),
    Dict(:x => 2.0, :y => [(; x = [MyStruct(2.0, [2.0, ForwardDiff.Dual(2.0)])], y = 2.0)]),
    [MyStruct("2", [2.0, ForwardDiff.Dual(2.0)])],
    Dict(
        :x => [MyStruct("2", [2.0, MyStruct(ForwardDiff.Dual(2.0), 2.0)])],
        :y => ForwardDiff.Dual{MyStruct}(2.0)
    ),
    ((Dict(:x => nothing)), ForwardDiff.Dual(2.0)),
    MyStruct2(ForwardDiff.Dual(2.0)),
    [MyStruct2(ForwardDiff.Dual(2.0)), 2.0],

    # Vectors of non-number types won't infer
    [MyStruct(2.0, ForwardDiff.Dual(2.0))],
    (; x = 2.0, y = [[MyStruct3(ForwardDiff.Dual(2.0))]]),
    (; x = Vector{Float64}(undef, 2), y = [[MyStruct3(ForwardDiff.Dual(2.0))]]),
    (; x = Matrix{Any}(undef, 2, 2), y = [[MyStruct3(ForwardDiff.Dual(2.0))]]),
]

for p in p_possibilities_uninferrred
    @show p
    @test DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual
    local u0 = 2.0
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @test DiffEqBase.promote_u0([cis(u0)], p, t0) isa
        AbstractArray{<:Complex{<:ForwardDiff.Dual}}
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @test DiffEqBase.promote_u0([cis(u0)], p, t0) isa
        AbstractArray{<:Complex{<:ForwardDiff.Dual}}
end

p_possibilities_missed = [
    Set([2.0, "s", ForwardDiff.Dual(2.0)]),
    Set([2.0, ForwardDiff.Dual(2.0), SciMLBase.NullParameters()]),
    Set([Matrix{Float64}(undef, 2, 2), ForwardDiff.Dual(2.0)]),
]

for p in p_possibilities_missed
    @show p
    @test_broken DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual
    local u0 = 2.0
    @test_broken DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @test_broken DiffEqBase.promote_u0([cis(u0)], p, t0) isa
        AbstractArray{<:Complex{<:ForwardDiff.Dual}}
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @test DiffEqBase.promote_u0([cis(u0)], p, t0) isa
        AbstractArray{<:Complex{<:ForwardDiff.Dual}}
end

p_possibilities_notdual = [
    (), (;), [2.0], [2.0, 2], [2.0, (2.0)], [2.0, MyStruct(2.0, 2.0f0)], pairs((;)),
]

for p in p_possibilities_notdual
    @show p
    @test !(DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual)
    local u0 = 2.0
    @test !(DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual)
    @test !(
        DiffEqBase.promote_u0([cis(u0)], p, t0) isa
            AbstractArray{<:Complex{<:ForwardDiff.Dual}}
    )
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @test DiffEqBase.promote_u0([cis(u0)], p, t0) isa
        AbstractArray{<:Complex{<:ForwardDiff.Dual}}
    @inferred DiffEqBase.anyeltypedual(p)
end

p_possibilities_notdual_uninferred = [
    [],

    # Undefs cause inference loss
    [2.0, MyStruct3(2.0)], [2.0, MyStruct2(2.0)], [2.0, MyStruct2(2.0), []],
    [Dict(:x => 2, "y" => 5), MyStruct2(2.0)],

    # Dictionaries can have inference issues
    Dict(:x => 2, :y => 5), Dict(:x => 2, "y" => 5),
]

# Also check circular references
# https://github.com/SciML/DiffEqBase.jl/issues/784

x = Any[[1.0, 2.0]]
push!(x, x)
push!(p_possibilities_notdual_uninferred, x)

struct X
    x::Any
end
x = Any[[1.0, 2.0]]
push!(x, X(x))
push!(p_possibilities_notdual_uninferred, x)

mutable struct Y
    x::Any
end
x = Y(1)
x.x = x
push!(p_possibilities_notdual_uninferred, x)

for p in p_possibilities_notdual_uninferred
    @test !(DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual)
    local u0 = 2.0
    @test !(DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual)
    @test !(
        DiffEqBase.promote_u0([cis(u0)], p, t0) isa
            AbstractArray{<:Complex{<:ForwardDiff.Dual}}
    )
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @test DiffEqBase.promote_u0([cis(u0)], p, t0) isa
        AbstractArray{<:Complex{<:ForwardDiff.Dual}}
end

f(du, u, p, t) = du .= u
config = ForwardDiff.JacobianConfig(f, ones(5))

p_possibilities_configs = [
    (config, config), (config, 2.0), config, (; x = config, y = 2.0),
]

for p in p_possibilities_configs
    @show p
    @test !(DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual)
    local u0 = 2.0
    @test !(DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual)
    @test !(
        DiffEqBase.promote_u0([cis(u0)], p, t0) isa
            AbstractArray{<:Complex{<:ForwardDiff.Dual}}
    )
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @test DiffEqBase.promote_u0([cis(u0)], p, t0) isa
        AbstractArray{<:Complex{<:ForwardDiff.Dual}}
    @inferred DiffEqBase.anyeltypedual(p)
end

p_possibilities_configs_not_inferred = [
    [2.0, (2.0,), config], [2.0, config, MyStruct(2.0, 2.0f0)],
]

for p in p_possibilities_configs_not_inferred
    @show p
    @test !(DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual)
    local u0 = 2.0
    @test !(DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual)
    @test !(
        DiffEqBase.promote_u0([cis(u0)], p, t0) isa
            AbstractArray{<:Complex{<:ForwardDiff.Dual}}
    )
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @test DiffEqBase.promote_u0([cis(u0)], p, t0) isa
        AbstractArray{<:Complex{<:ForwardDiff.Dual}}
end

# use `getfield` on `Pairs`, see https://github.com/JuliaLang/julia/pull/39448
@test_nowarn DiffEqBase.DualEltypeChecker(pairs((;)), 0)(Val(:data))

# https://discourse.julialang.org/t/type-instability-with-differentialequations-jl-when-using-nested-structs/109764/5
struct Fit
    m₁::Float64
    c₁::Float64
    m₂::Float64
    c₂::Float64

    function Fit()
        m₁ = 1.595
        c₁ = 3.438
        m₂ = 1.075
        c₂ = 3.484

        return new(m₁, c₁, m₂, c₂)
    end
end

struct EOS
    fit::Fit

    function EOS()
        fit = Fit()
        return new(fit)
    end
end

p = EOS()
@test !(DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual)
@inferred DiffEqBase.anyeltypedual(p)

# Check methods used for prevention of Dual-detection when using
# DiffResults.DiffResult in a wrapper.
# https://github.com/SciML/DiffEqBase.jl/issues/1009

struct OutsideWrapper{T}
    a::Float64
    b::T
end

struct InsideWrapper{T, S}
    du::T
    dual_du::S
end

f(x) = 2 * x[1] + 3 * x[2]^2
xdual = ones(
    ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}, 2
)
x = [1.0, 1.0]
diffresult = ForwardDiff.DiffResults.GradientResult(x)
diffresult_dual = ForwardDiff.DiffResults.GradientResult(xdual)
iw = InsideWrapper(diffresult, diffresult_dual)
ow = OutsideWrapper(1.0, iw)

@test !(DiffEqBase.anyeltypedual(iw) <: ForwardDiff.Dual)
@test !(DiffEqBase.anyeltypedual(ow) <: ForwardDiff.Dual)
@inferred DiffEqBase.anyeltypedual(iw)
@inferred DiffEqBase.anyeltypedual(ow)

# Issue https://github.com/SciML/ModelingToolkit.jl/issues/2717
u0 = [1.0, 2.0, 3.0]
p = [1, 2]
t = ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}(1.0)
@test DiffEqBase.promote_u0(u0, p, t) isa AbstractArray{<:ForwardDiff.Dual}
u0 = [1.0 + 1im, 2.0, 3.0]
@test DiffEqBase.promote_u0(u0, p, t) isa AbstractArray{<:Complex{<:ForwardDiff.Dual}}

# Issue https://github.com/SciML/NonlinearSolve.jl/issues/440
f_nlsolve(u, p, t) = [u[2], 1.5u[1]^2]
ode = ODEProblem(f_nlsolve, [0.0, 0.0], (0, 1))
@inferred DiffEqBase.anyeltypedual(ode)
ode = NonlinearProblem(f_nlsolve, [0.0, 0.0], (0, 1))
@inferred DiffEqBase.anyeltypedual(ode)

# Issue https://github.com/SciML/DiffEqBase.jl/issues/1021
f_1021(u, p, t) = 1.01 * u
struct Foo{T}
    sol::T
end
u0 = 1 / 2
tspan = (0.0, 1.0)
prob = ODEProblem{false}(f_1021, u0, tspan)
foo = SciMLBase.build_solution(
    prob, DiffEqBase.InternalEuler.FwdEulerAlg(), [u0, u0], [0.0, 1.0]
)
DiffEqBase.anyeltypedual((; x = foo))
DiffEqBase.anyeltypedual((; x = foo, y = prob.f))

@test DiffEqBase.anyeltypedual(ReverseDiff.track(ones(3))) == Any
@test DiffEqBase.anyeltypedual(typeof(ReverseDiff.track(ones(3)))) == Any
@test DiffEqBase.anyeltypedual(ReverseDiff.track(ones(ForwardDiff.Dual, 3))) ==
    eltype(ones(ForwardDiff.Dual, 3))
@test DiffEqBase.anyeltypedual(typeof(ReverseDiff.track(ones(ForwardDiff.Dual, 3)))) ==
    eltype(ones(ForwardDiff.Dual, 3))

struct FakeParameterObject{T}
    tunables::T
end

SciMLStructures.isscimlstructure(::FakeParameterObject) = true
function SciMLStructures.canonicalize(::SciMLStructures.Tunable, f::FakeParameterObject)
    return f.tunables, x -> FakeParameterObject(x), true
end

@test DiffEqBase.promote_u0(
    ones(3), FakeParameterObject(ReverseDiff.track(ones(3))), 0.0
) isa
    ReverseDiff.TrackedArray
@test DiffEqBase.promote_u0(1.0, FakeParameterObject(ReverseDiff.track(ones(3))), 0.0) isa
    ReverseDiff.TrackedReal
@test DiffEqBase.promote_u0(
    ones(3), FakeParameterObject(ReverseDiff.track(ones(ForwardDiff.Dual, 3))), 0.0
) isa
    ReverseDiff.TrackedArray{<:ForwardDiff.Dual}
@test DiffEqBase.promote_u0(
    1.0, FakeParameterObject(ReverseDiff.track(ones(ForwardDiff.Dual, 3))), 0.0
) isa
    ReverseDiff.TrackedReal{<:ForwardDiff.Dual}
@test DiffEqBase.promote_u0(NaN, [NaN], 0.0) isa Float64
@test DiffEqBase.promote_u0([1.0], [NaN], 0.0) isa Vector{Float64}

# totallength
val = rand(10)
par = rand(10)
u = ForwardDiff.Dual.(val, par)
@test DiffEqBase.totallength(val[1]) == 1
@test DiffEqBase.totallength(val) == length(val)
@test DiffEqBase.totallength(par) == length(par)
@test DiffEqBase.totallength(u[1]) ==
    DiffEqBase.totallength(val[1]) + DiffEqBase.totallength(par[1])
@test DiffEqBase.totallength(u) == sum(DiffEqBase.totallength, u)
