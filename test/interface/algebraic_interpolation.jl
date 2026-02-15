using OrdinaryDiffEq
using Test

function rober_ip(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    return nothing
end
function rober_op(u, p, t)
    du = similar(u)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    return du
end
M = [
    1.0 0 0
    0 1.0 0
    0 0 0
]
fip = ODEFunction(rober_ip, mass_matrix = M)
prob_ip = ODEProblem(fip, [1.0, 0.0, 0.0], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))

fop = ODEFunction(rober_op, mass_matrix = M)
prob_op = ODEProblem(fop, [1.0, 0.0, 0.0], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))

ref_ip = solve(prob_ip, Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-8)
ref_op = solve(prob_op, Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-8)

sol_ip = solve(prob_ip, FBDF(), reltol = 1.0e-8, abstol = 1.0e-8)
sol_op = solve(prob_op, FBDF(), reltol = 1.0e-8, abstol = 1.0e-8)

# make sure interpolation changes don't accidentally break this test suite
# both ref (Rodas5P) and sol (FBDF) use stiffness-aware interpolation
@test occursin("stiffness-aware", SciMLBase.interp_summary(sol_ip))
@test occursin("stiffness-aware", SciMLBase.interp_summary(sol_op))
@test occursin("stiffness-aware", SciMLBase.interp_summary(ref_ip))
@test occursin("stiffness-aware", SciMLBase.interp_summary(ref_op))

reltol = 1.0e-4
abstol = 1.0e-4
t = 1
tv = [1, 10, 100]
idxs = 3
idxsv = [2, 3]

# primal, no index
@test isapprox(ref_ip(t), sol_ip(t), rtol = reltol, atol = abstol) # ip, t
@test isapprox(ref_ip(tv), sol_ip(tv), rtol = reltol, atol = abstol) # ip, tv
@test isapprox(ref_op(t), sol_op(t), rtol = reltol, atol = abstol) # op, t
@test isapprox(ref_op(tv), sol_op(tv), rtol = reltol, atol = abstol) # op, tv

# primal, scalar index
@test isapprox(ref_ip(t, idxs = idxs), sol_ip(t, idxs = idxs), rtol = reltol, atol = abstol) # ip, t
@test isapprox(
    ref_ip(tv, idxs = idxs), sol_ip(tv, idxs = idxs), rtol = reltol, atol = abstol
) # ip, tv
@test isapprox(ref_op(t, idxs = idxs), sol_op(t, idxs = idxs), rtol = reltol, atol = abstol) # op, t
@test isapprox(
    ref_op(tv, idxs = idxs), sol_op(tv, idxs = idxs), rtol = reltol, atol = abstol
) # op, tv

# primal, vector index
@test isapprox(
    ref_ip(t, idxs = idxsv), sol_ip(t, idxs = idxsv), rtol = reltol, atol = abstol
)
@test isapprox(
    ref_ip(tv, idxs = idxsv), sol_ip(tv, idxs = idxsv), rtol = reltol, atol = abstol
)
@test isapprox(
    ref_op(t, idxs = idxsv), sol_op(t, idxs = idxsv), rtol = reltol, atol = abstol
)
@test isapprox(
    ref_op(tv, idxs = idxsv), sol_op(tv, idxs = idxsv), rtol = reltol, atol = abstol
)

abstol = 1.0e-3
# derivative, no index
@test isapprox(ref_ip(t, Val{1}), sol_ip(t, Val{1}), rtol = reltol, atol = abstol)
@test isapprox(ref_ip(tv, Val{1}), sol_ip(tv, Val{1}), rtol = reltol, atol = abstol)
@test isapprox(ref_op(t, Val{1}), sol_op(t, Val{1}), rtol = reltol, atol = abstol)
@test isapprox(ref_op(tv, Val{1}), sol_op(tv, Val{1}), rtol = reltol, atol = abstol)

# derivative, scalar index
@test isapprox(
    ref_ip(t, Val{1}, idxs = idxs),
    sol_ip(t, Val{1}, idxs = idxs), rtol = reltol, atol = abstol
)
@test isapprox(
    ref_ip(tv, Val{1}, idxs = idxs),
    sol_ip(tv, Val{1}, idxs = idxs), rtol = reltol, atol = abstol
)
@test isapprox(
    ref_op(t, Val{1}, idxs = idxs),
    sol_op(t, Val{1}, idxs = idxs), rtol = reltol, atol = abstol
)
@test isapprox(
    ref_op(tv, Val{1}, idxs = idxs),
    sol_op(tv, Val{1}, idxs = idxs), rtol = reltol, atol = abstol
)

# derivative, vector index
@test isapprox(
    ref_ip(tv, Val{1}, idxs = idxsv),
    sol_ip(tv, Val{1}, idxs = idxsv), rtol = reltol, atol = abstol
)
@test isapprox(
    ref_ip(t, Val{1}, idxs = idxsv),
    sol_ip(t, Val{1}, idxs = idxsv), rtol = reltol, atol = abstol
)
@test isapprox(
    ref_op(t, Val{1}, idxs = idxsv),
    sol_op(t, Val{1}, idxs = idxsv), rtol = reltol, atol = abstol
)
@test isapprox(
    ref_op(tv, Val{1}, idxs = idxsv),
    sol_op(tv, Val{1}, idxs = idxsv), rtol = reltol, atol = abstol
)

# higher derivatives should be zero
# second derivative, no index
@test (sol_ip(t, Val{2}) .== 0) == [false, false, true]
@test all(map(==([false, false, true]), sol_ip(tv, Val{2}) .== 0))
@test (sol_op(t, Val{2}) .== 0) == [false, false, true]
@test all(map(==([false, false, true]), sol_op(tv, Val{2}) .== 0))

# second derivative, scalar index
@test sol_ip(t, Val{2}, idxs = idxs) == 0
@test all(sol_ip(tv, Val{2}, idxs = idxs) .== 0)
@test sol_op(t, Val{2}, idxs = idxs) == 0
@test all(sol_op(tv, Val{2}, idxs = idxs) .== 0)

# second derivative, vector index
@test (sol_ip(t, Val{2}, idxs = idxsv) .== 0) == [false, true]
@test all(map(==([false, true]), sol_ip(tv, Val{2}, idxs = idxsv) .== 0))
@test (sol_op(t, Val{2}, idxs = idxsv) .== 0) == [false, true]
@test all(map(==([false, true]), sol_op(tv, Val{2}, idxs = idxsv) .== 0))

# third derivative, no index
@test (sol_ip(t, Val{3}) .== 0) == [false, false, true]
@test all(map(==([false, false, true]), sol_ip(tv, Val{3}) .== 0))
@test (sol_op(t, Val{3}) .== 0) == [false, false, true]
@test all(map(==([false, false, true]), sol_op(tv, Val{3}) .== 0))

# third derivative, scalar index
@test sol_ip(t, Val{3}, idxs = idxs) == 0
@test all(sol_ip(tv, Val{3}, idxs = idxs) .== 0)
@test sol_op(t, Val{3}, idxs = idxs) == 0
@test all(sol_op(tv, Val{3}, idxs = idxs) .== 0)

# third derivative, vector index
@test (sol_ip(t, Val{3}, idxs = idxsv) .== 0) == [false, true]
@test all(map(==([false, true]), sol_ip(tv, Val{3}, idxs = idxsv) .== 0))
@test (sol_op(t, Val{3}, idxs = idxsv) .== 0) == [false, true]
@test all(map(==([false, true]), sol_op(tv, Val{3}, idxs = idxsv) .== 0))
