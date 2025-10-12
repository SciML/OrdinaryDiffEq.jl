using OrdinaryDiffEqVerner, OrdinaryDiffEqCore, DiffEqBase, Test
using LinearAlgebra
using OrdinaryDiffEqSSPRK, DiffEqDevTools, Test, Random
import OrdinaryDiffEqLowStorageRK
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear
using Plots

testTol = 0.3
println(BigFloat(exp(-64.0)))
function print_solution_dts(prob, algorithms; abstol=1e-12, reltol=1e-12, show_solution=true, kwargs...)
    println("\n" * "="^70)
    println("Comparing dt usage across algorithms")
    println("="^70)

    # println("\nProblem Details:")
    # println("  Time span: ", prob.tspan)
    # println("  Initial condition u0: ", prob.u0)
    # println("  Problem type: ", typeof(prob.f).name.name)

    for alg in algorithms
        alg_name = string(typeof(alg).name.name)

        try
            sol = solve(prob, alg; abstol=abstol, reltol=reltol, kwargs...)
            final_dt = sol.t[end] - sol.t[end-1]
            t_end = sol.t[end]

            println("\nAlgorithm: $alg_name")
            if show_solution
                println("  Solution at t=$t_end: ", sol[end])
            end
            println("  Final dt: ", final_dt)
        catch e
            println("\nAlgorithm: $alg_name")
            println("  ERROR: ", e)
        end
    end

    println("\n" * "="^70)
end

# ODE function definitions
f_1 = (u, p, t) -> cos(t)
prob_ode_sin = ODEProblem(ODEFunction(f_1; analytic = (u0, p, t) -> sin(t)), 0.0, (0.0, 1.0))

f_1 = (du, u, p, t) -> du[1] = cos(t)
prob_ode_sin_inplace = ODEProblem(ODEFunction(f_1; analytic = (u0, p, t) -> [sin(t)]), [0.0],
    (0.0, 1.0))

f_2 = (u, p, t) -> sin(u)
prob_ode_nonlinear = ODEProblem(
    ODEFunction(f_2;
        analytic = (u0, p, t) -> 2 * acot(exp(-t) *
                                          cot(0.5))), 1.0,
    (0.0, 0.5))

f_2 = (du, u, p, t) -> du[1] = sin(u[1])
prob_ode_nonlinear_inplace = ODEProblem(
    ODEFunction(f_2;
        analytic = (u0, p, t) -> [
            2 * acot(exp(-t) * cot(0.5))
        ]),
    [1.0], (0.0, 0.5))

f_ssp = (u, p, t) -> begin
    sin(10t) * u * (1 - u)
end
test_problem_ssp = ODEProblem(f_ssp, 0.1, (0.0, 8.0))
test_problem_ssp_long = ODEProblem(f_ssp, 0.1, (0.0, 1.e3))

function f!(du, u, p, t)
    du[1] = -u[1]
end

function f(u, p, t)
    -u
end

test_problems_only_time = [prob_ode_sin, prob_ode_sin_inplace]


t_end=1.0
alg = OrdinaryDiffEqSSPRK.SSPRK22()
t_end = 64.0

setprecision(256)
prob_oop = ODEProblem(f, 1.0, (0.0, t_end))
println("**** exp(-64) ****")
algorithms = [RKV76IIa(), Vern7()]
print_solution_dts(prob_oop, algorithms; abstol=1e-40, reltol=1e-40)
print_solution_dts(prob_oop, [alg]; dt=OrdinaryDiffEqSSPRK.ssp_coefficient(alg), dense=false,abstol=1e-40, reltol=1e-40)
println("Expected value: ", exp(BigFloat(-t_end)))
println("**** exp(-64) ****")


println("***************** sin in and out of place *********************")
algorithms = [RKV76IIa(), Vern7()]
for prob in test_problems_only_time
    print_solution_dts(prob, algorithms; abstol=1e-12, reltol=1e-12)
    print_solution_dts(prob, [alg]; dt=OrdinaryDiffEqSSPRK.ssp_coefficient(alg), dense=false)
end
println("************** sin in and out of place ************************")


test_problems_nonlinear = [prob_ode_nonlinear, prob_ode_nonlinear_inplace]
t_end = 1.e3

sol_oop = solve(test_problem_ssp_long, RKV76IIa(), abstol=1e-12, reltol=1e-12)
println("***************** test_problem_ssp_long *********************")
algorithms = [RKV76IIa(), Vern7(), alg]
print_solution_dts(test_problem_ssp_long, algorithms; abstol=1e-12, reltol=1e-12)
print_solution_dts(test_problem_ssp_long, [alg]; dt=OrdinaryDiffEqSSPRK.ssp_coefficient(alg), dense=false)
println("************** test_problem_ssp_long ************************")

t_end = 64.0

setprecision(256)
prob_oop = ODEProblem(f, 1.0, (0.0, t_end))
println("**** exp(-64) ****")
algorithms = [RKV76IIa(), Vern7()]
print_solution_dts(prob_oop, algorithms; abstol=1e-40, reltol=1e-40)
print_solution_dts(prob_oop, [alg]; dt=OrdinaryDiffEqSSPRK.ssp_coefficient(alg), dense=false)
println("Expected value: ", exp(BigFloat(-t_end)))
println("**** exp(-64) ****")


println("*** Testing Convergence of diff algorithm *** ")

#alg=Vern7()
alg=RKV76IIa()
# dts = BigFloat(1) ./ 2 .^ (1:6)
dts = [8, 6, 4, 2, 1, 0.5, 0.25, 0.125]

errors = zeros(BigFloat, length(dts))
println("Testing order 7 for RKV76IIa()")
for (i, dt) in enumerate(dts)
    sol = solve(prob_oop, alg, dt=dt, adaptive=false)
    # Use BigFloat for error calculation
    errors[i] = abs(BigFloat(sol[end]) - exp(BigFloat(-t_end)))
    println("Computed Solution ", sol[end], " for dt = ", dt, ", error = ", errors[i])
end

for i in 2:length(errors)
    order = log(BigFloat(errors[i-1])/BigFloat(errors[i])) / log(2)
    println("Order between dt=", dts[i-1], " and dt=", dts[i], ": ", order)
end

plot(
    dts, errors;
    xscale = :log10, yscale = :log10,
    marker = :o, linewidth = 2,
    xlabel = "dt", ylabel = "Error",
    title = "Convergence of RKV76IIa",
    label = "Observed Error"
)
# Make reference line pass through the 0.125 dt point
ref_idx = findfirst(x -> x == 0.125, dts)
ref_errors = errors[ref_idx] * (BigFloat.(dts) ./ BigFloat(0.125)).^7
plot!(float.(dts), ref_errors; linestyle = :dash, label = "Order 7 Reference")
display(current())
savefig("convergence_rkv76iia.png")