using OrdinaryDiffEq, LinearAlgebra, LinearSolve
n = 8
dt = 1/16
u0 = ones(n)
tspan = (0.0,1.0)

M1 = 2ones(n) |> Diagonal |> Array
M2 = 2ones(n) |> Diagonal |> Array

f1 = M1 |> DiffEqArrayOperator
f2 = M2 |> DiffEqArrayOperator
prob = SplitODEProblem(f1,f2,u0,tspan)

#=
for algname in (
                :SBDF2,
#               :SBDF3,
#               :KenCarp47,
               )
    println("##########################################################")
    println("Testing out $algname")
    println("##########################################################")

    # prepare_alg() sets alg.linsovle=GenericFactorization()
    alg0 = @eval $algname()
    alg1 = @eval $algname(linsolve=GenericFactorization())
#   alg2 = DiffEqBase.prepare_alg(alg0, prob.u0, prob.p, prob)

    kwargs = (dt=dt,)

#   @show sol  =              solve(prob, alg0; kwargs...).retcode # fails
    @show sol0 = DiffEqBase.__solve(prob, alg0; kwargs...).retcode # passes
    @show sol1 = DiffEqBase.__solve(prob, alg1; kwargs...).retcode # fails
#   @show sol2 = DiffEqBase.__solve(prob, alg2; kwargs...).retcode # fails
end
=#

#####
alg0 = KenCarp47()                                # passing case
alg1 = KenCarp47(linsolve=GenericFactorization()) # failing case

ig0 = SciMLBase.init(prob, alg0; dt=dt)
ig1 = SciMLBase.init(prob, alg1; dt=dt)

lincache0 = ig0.cache.nlsolver.cache.linsolve
lincache1 = ig1.cache.nlsolver.cache.linsolve

# check operators
W0 = lincache0.A
W1 = lincache1.A

for field in W0 |> typeof |> fieldnames
    isequal = getfield(W0, field) == getfield(W1, field)
    println("field $field is $isequal")
end

@show W0._concrete_form == W1._concrete_form # true
@show W0._func_cache    == W1._func_cache    # true

# check operator application
b = ones(n)
@show W0 * b == W1 * b                             # true
@show mul!(rand(n), W0, b) == mul!(rand(n), W1, b) # true
@show W0 \ b == W1 \ b                             # true

# check solve
lincache0.b .= 1.0
lincache1.b .= 1.0

solve(lincache0)
solve(lincache1)

@show lincache0.u == lincache1.u # false

@show lincache0.u'
@show lincache1.u'

#import OrdinaryDiffEq.dolinsolve
#linres0 = dolinsolve(ig0, lincache0; A = W0, b = b, linu = ones(n), reltol = 1e-8)
#linres1 = dolinsolve(ig1, lincache1; A = W1, b = b, linu = ones(n), reltol = 1e-8)
#
#@show linres0 == linres1

nothing
