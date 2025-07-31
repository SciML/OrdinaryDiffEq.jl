using OrdinaryDiffEq, BenchmarkTools, DiffEqBase
using LinearAlgebra, SparseArrays, StaticArrays, StableRNGs

const SUITE = BenchmarkGroup()

# =============================================================================
# Non-Stiff ODE Problems
# =============================================================================

SUITE["nonstiff"] = BenchmarkGroup()

# Lotka-Volterra (Predator-Prey) Problem
function lotka_volterra!(du, u, p, t)
    α, β, γ, δ = p
    x, y = u
    du[1] = α*x - β*x*y  # dx/dt
    du[2] = δ*x*y - γ*y  # dy/dt
    nothing
end

function lotka_volterra_prob()
    p = (1.5, 1.0, 3.0, 1.0)  # α, β, γ, δ
    u0 = [1.0, 1.0]
    tspan = (0.0, 10.0)
    ODEProblem(lotka_volterra!, u0, tspan, p)
end

# Pleiades Problem (7-body celestial mechanics)
function pleiades!(du, u, p, t)
    # u[1:7] = x positions, u[8:14] = y positions
    # u[15:21] = x velocities, u[22:28] = y velocities
    x = @view u[1:7]
    y = @view u[8:14]
    vx = @view u[15:21]
    vy = @view u[22:28]

    # Copy velocities to position derivatives
    du[1:7] .= vx
    du[8:14] .= vy

    # Calculate accelerations
    fill!(du[15:21], 0.0)
    fill!(du[22:28], 0.0)

    for i in 1:7
        for j in 1:7
            if i != j
                dx = x[j] - x[i]
                dy = y[j] - y[i]
                r = sqrt(dx^2 + dy^2)
                r3 = r^3
                du[14 + i] += j * dx / r3  # mass j = j
                du[21 + i] += j * dy / r3
            end
        end
    end
    nothing
end

function pleiades_prob()
    # Initial conditions from literature
    u0 = [3.0, 3.0, -1.0, -3.0, 2.0, -2.0, 2.0, 3.0, -3.0, 2.0, 0, 0, -4.0, 4.0,
        0, 0, 0, 0, 0, 1.75, -1.5, 0, 0, 0, -1.25, 1, 0, 0]
    tspan = (0.0, 3.0)
    ODEProblem(pleiades!, u0, tspan)
end

# FitzHugh-Nagumo Model
function fitzhugh_nagumo!(du, u, p, t)
    a, b, c = p
    v, w = u
    du[1] = c * (v - v^3/3 + w)    # dv/dt
    du[2] = -(v - a + b*w) / c     # dw/dt
    nothing
end

function fitzhugh_nagumo_prob()
    p = (0.7, 0.8, 12.5)
    u0 = [-1.0, 1.0]
    tspan = (0.0, 20.0)
    ODEProblem(fitzhugh_nagumo!, u0, tspan, p)
end

# =============================================================================
# Stiff ODE Problems
# =============================================================================

SUITE["stiff"] = BenchmarkGroup()

# ROBER Problem (Robertson chemical kinetics)
function rober!(du, u, p, t)
    k1, k2, k3 = p
    y1, y2, y3 = u
    du[1] = -k1*y1 + k3*y2*y3         # dy1/dt
    du[2] = k1*y1 - k2*y2^2 - k3*y2*y3   # dy2/dt
    du[3] = k2*y2^2                   # dy3/dt
    nothing
end

function rober_prob()
    p = (0.04, 3e7, 1e4)  # k1, k2, k3
    u0 = [1.0, 0.0, 0.0]
    tspan = (0.0, 1e5)
    ODEProblem(rober!, u0, tspan, p)
end

# Van der Pol Oscillator (stiff)
function van_der_pol!(du, u, p, t)
    μ = p[1]
    x, y = u
    du[1] = y                      # dx/dt
    du[2] = μ * ((1 - x^2)*y - x)  # dy/dt
    nothing
end

function van_der_pol_prob()
    p = [1e6]  # very stiff
    u0 = [1.0, 1.0]
    tspan = (0.0, 6.3)
    ODEProblem(van_der_pol!, u0, tspan, p)
end

# Pollution Problem (atmospheric chemistry, 20D stiff system)
const k1=.35e0
const k2=.266e2
const k3=.123e5
const k4=.86e-3
const k5=.82e-3
const k6=.15e5
const k7=.13e-3
const k8=.24e5
const k9=.165e5
const k10=.9e4
const k11=.22e-1
const k12=.12e5
const k13=.188e1
const k14=.163e5
const k15=.48e7
const k16=.35e-3
const k17=.175e-1
const k18=.1e9
const k19=.444e12
const k20=.124e4
const k21=.21e1
const k22=.578e1
const k23=.474e-1
const k24=.178e4
const k25=.312e1

function pollution!(dy, y, p, t)
    r1 = k1 * y[1]
    r2 = k2 * y[2] * y[4]
    r3 = k3 * y[5] * y[2]
    r4 = k4 * y[7]
    r5 = k5 * y[7]
    r6 = k6 * y[7] * y[6]
    r7 = k7 * y[9]
    r8 = k8 * y[9] * y[6]
    r9 = k9 * y[11] * y[2]
    r10 = k10*y[11]*y[1]
    r11 = k11*y[13]
    r12 = k12*y[10]*y[2]
    r13 = k13*y[14]
    r14 = k14*y[1]*y[6]
    r15 = k15*y[3]
    r16 = k16*y[4]
    r17 = k17*y[4]
    r18 = k18*y[16]
    r19 = k19*y[16]
    r20 = k20*y[17]*y[6]
    r21 = k21*y[19]
    r22 = k22*y[19]
    r23 = k23*y[1]*y[4]
    r24 = k24*y[19]*y[1]
    r25 = k25*y[20]

    dy[1] = -r1-r10-r14-r23-r24 +
            r2 + r3 + r9 + r11 + r12 + r22 + r25
    dy[2] = -r2-r3-r9-r12+r1+r21
    dy[3] = -r15+r1+r17+r19+r22
    dy[4] = -r2-r16-r17-r23+r15
    dy[5] = -r3+r4+r4+r6+r7+r13+r20
    dy[6] = -r6-r8-r14-r20+r3+r18+r18
    dy[7] = -r4-r5-r6+r13
    dy[8] = r4+r5+r6+r7
    dy[9] = -r7-r8
    dy[10] = -r12+r7+r9
    dy[11] = -r9-r10+r8+r11
    dy[12] = r9
    dy[13] = -r11+r10
    dy[14] = -r13+r12
    dy[15] = r14
    dy[16] = -r18-r19+r16
    dy[17] = -r20
    dy[18] = r20
    dy[19] = -r21-r22-r24+r23+r25
    dy[20] = -r25+r24
    nothing
end

function pollution_prob()
    u0 = zeros(20)
    u0[2] = 0.2
    u0[4] = 0.04
    u0[7] = 0.1
    u0[8] = 0.3
    u0[9] = 0.01
    u0[17] = 0.007
    tspan = (0.0, 60.0)
    ODEProblem(pollution!, u0, tspan)
end

# =============================================================================
# Scaling Problems (different dimensions)
# =============================================================================

SUITE["scaling"] = BenchmarkGroup()

# Linear ODE with varying dimensions
function linear_ode!(du, u, p, t)
    A = p
    mul!(du, A, u)
    nothing
end

function create_linear_prob(N::Int)
    rng = StableRNG(123)  # Fixed seed for reproducibility
    A = randn(rng, N, N)
    A = A - A'  # Make skew-symmetric for bounded solutions
    u0 = randn(rng, N)
    tspan = (0.0, 1.0)
    ODEProblem(linear_ode!, u0, tspan, A)
end

# Brusselator PDE (2D reaction-diffusion from advanced ODE example)
# Forcing function - creates a localized disturbance
brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0

# Periodic boundary condition helper
limit(a, N) = a == N + 1 ? 1 : a == 0 ? N : a

function brusselator_2d!(du, u, p, t)
    A, B, α, dx, N = p
    α = α / dx^2

    # Create coordinate arrays for this N
    xyd = range(0, stop = 1, length = N)

    @inbounds for I in CartesianIndices((N, N))
        i, j = Tuple(I)
        x, y = xyd[i], xyd[j]
        ip1, im1, jp1,
        jm1 = limit(i + 1, N), limit(i - 1, N), limit(j + 1, N), limit(j - 1, N)

        # u equation: ∂u/∂t = 1 + u²v - 4.4u + α∇²u + f(x,y,t)
        du[i, j, 1] = α * (u[im1, j, 1] + u[ip1, j, 1] + u[i, jp1, 1] + u[i, jm1, 1] -
                       4*u[i, j, 1]) +
                      B + u[i, j, 1]^2 * u[i, j, 2] - (A + 1) * u[i, j, 1] +
                      brusselator_f(x, y, t)

        # v equation: ∂v/∂t = 3.4u - u²v + α∇²v  
        du[i, j, 2] = α * (u[im1, j, 2] + u[ip1, j, 2] + u[i, jp1, 2] + u[i, jm1, 2] -
                       4*u[i, j, 2]) +
                      A * u[i, j, 1] - u[i, j, 1]^2 * u[i, j, 2]
    end
    nothing
end

function init_brusselator_2d(N::Int)
    xyd = range(0, stop = 1, length = N)
    u = zeros(N, N, 2)
    for I in CartesianIndices((N, N))
        x = xyd[I[1]]
        y = xyd[I[2]]
        u[I, 1] = 22 * (y * (1 - y))^(3/2)  # u initial condition
        u[I, 2] = 27 * (x * (1 - x))^(3/2)  # v initial condition  
    end
    u
end

function create_brusselator_2d_prob(N::Int)
    A, B, α = 3.4, 1.0, 10.0
    xyd = range(0, stop = 1, length = N)
    dx = step(xyd)
    u0 = init_brusselator_2d(N)
    tspan = (0.0, 11.5)
    ODEProblem(brusselator_2d!, u0, tspan, (A, B, α, dx, N))
end

# =============================================================================
# Benchmark Definitions
# =============================================================================

# Non-stiff benchmarks with different solvers
lv_prob = lotka_volterra_prob()
pl_prob = pleiades_prob()
fn_prob = fitzhugh_nagumo_prob()

# Explicit RK methods for non-stiff problems
explicit_solvers = [Tsit5(), Vern6(), Vern7(), DP5(), BS3()]

SUITE["nonstiff"]["lotka_volterra"] = BenchmarkGroup()
SUITE["nonstiff"]["pleiades"] = BenchmarkGroup()
SUITE["nonstiff"]["fitzhugh_nagumo"] = BenchmarkGroup()

for solver in explicit_solvers
    solver_name = string(typeof(solver).name.name)
    SUITE["nonstiff"]["lotka_volterra"][solver_name] = @benchmarkable solve(
        $lv_prob, $solver, reltol = 1e-6, abstol = 1e-8)
    SUITE["nonstiff"]["pleiades"][solver_name] = @benchmarkable solve(
        $pl_prob, $solver, reltol = 1e-6, abstol = 1e-8)
    SUITE["nonstiff"]["fitzhugh_nagumo"][solver_name] = @benchmarkable solve(
        $fn_prob, $solver, reltol = 1e-6, abstol = 1e-8)
end

# Stiff benchmarks with different solvers
rober_prob_instance = rober_prob()
vdp_prob = van_der_pol_prob()
pollution_prob_instance = pollution_prob()

# Stiff solvers
stiff_solvers = [Rosenbrock23(), Rodas4(), TRBDF2(), KenCarp4(), FBDF()]

SUITE["stiff"]["rober"] = BenchmarkGroup()
SUITE["stiff"]["van_der_pol"] = BenchmarkGroup()
SUITE["stiff"]["pollution"] = BenchmarkGroup()

for solver in stiff_solvers
    solver_name = string(typeof(solver).name.name)
    SUITE["stiff"]["rober"][solver_name] = @benchmarkable solve(
        $rober_prob_instance, $solver, reltol = 1e-6, abstol = 1e-8)
    SUITE["stiff"]["van_der_pol"][solver_name] = @benchmarkable solve(
        $vdp_prob, $solver, reltol = 1e-6, abstol = 1e-8)
    SUITE["stiff"]["pollution"][solver_name] = @benchmarkable solve(
        $pollution_prob_instance, $solver, reltol = 1e-6, abstol = 1e-8)
end

# Scaling benchmarks
SUITE["scaling"]["linear"] = BenchmarkGroup()
SUITE["scaling"]["brusselator_2d"] = BenchmarkGroup()

# Linear ODE scaling (different problem sizes)
for N in [10, 50, 100]
    prob = create_linear_prob(N)
    SUITE["scaling"]["linear"]["N$N"] = @benchmarkable solve($prob, Tsit5(), reltol = 1e-6, abstol = 1e-8)
end

# Brusselator 2D scaling (different grid sizes)
for N in [8, 16, 32]
    prob = create_brusselator_2d_prob(N)
    SUITE["scaling"]["brusselator_2d"]["$(N)x$(N)"] = @benchmarkable solve($prob, TRBDF2(),
        reltol = 1e-4, abstol = 1e-6, maxiters = 1000)
end

# =============================================================================
# Problem Construction Benchmarks
# =============================================================================

SUITE["construction"] = BenchmarkGroup()

# Test problem construction overhead
SUITE["construction"]["lotka_volterra"] = @benchmarkable lotka_volterra_prob()
SUITE["construction"]["rober"] = @benchmarkable rober_prob()
SUITE["construction"]["linear_N50"] = @benchmarkable create_linear_prob(50)
