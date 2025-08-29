import LinearAlgebra
using CUDA
using FFTW

using OrdinaryDiffEqRKIP: RKIP
using SciMLBase: SplitODEProblem, SplitFunction, solve
using SciMLOperators: AbstractSciMLOperator

struct DiagonalFourierOperator{T, uType <: AbstractVector{T}, fftType, ifftType} <:
       AbstractSciMLOperator{T}
    diag::uType
    plan_fft!::fftType
    plan_ifft!::ifftType
    tmp::uType
end

function DiagonalFourierOperator(diag, plan_fft!_prototype, plan_ifft!_prototype)
    tmp = similar(diag)
    return DiagonalFourierOperator(
        diag, plan_fft!_prototype(tmp), plan_ifft!_prototype(tmp), tmp)
end

function DiagonalFourierOperator(diag::Vector{T}) where {T <: Complex}
    DiagonalFourierOperator(diag, FFTW.plan_fft!, FFTW.plan_ifft!)
end
function DiagonalFourierOperator(diag::CuArray{T, 1}) where {T <: Complex}
    DiagonalFourierOperator(diag, CUDA.CUFFT.plan_fft!, CUDA.CUFFT.plan_ifft!)
end

function LinearAlgebra.exp(D::DiagonalFourierOperator, t)
    exp_kernel = similar(D.diag)
    D.tmp .= D.diag
    D.tmp .*= t
    exp_kernel .= exp.(D.tmp)
    DiagonalFourierOperator(exp_kernel, D.plan_fft!, D.plan_ifft!, D.tmp)
end

@inline @fastmath function apply_kernel!(du, u, tmp, kernel, op_fft!, op_ifft!)
    copyto!(tmp, u)
    op_fft! * tmp
    tmp .*= kernel
    op_ifft! * tmp
    copyto!(du, tmp)
end

function (D::DiagonalFourierOperator)(du, u, _, _, _)
    apply_kernel!(du, u, D.tmp, D.diag, D.plan_fft!, D.plan_ifft!)
end

gpu_init(u0::CuArray, ::Vector{T}) where {T <: Real} = u0
gpu_init(u0::Vector, ::Vector{T}) where {T <: Real} = CuArray{Complex{T}, 1}(u0)
gpu_init(u0::Function, τ::Vector{T}) where {T <: Real} = CuArray{Complex{T}, 1}(u0.(τ))

cpu_init(u0::Vector, ::Vector{T}) where {T <: Real} = u0
cpu_init(u0::Function, τ::Vector{T}) where {T <: Real} = Vector{Complex{T}}(u0.(τ))

init(u0, τ, ::Val{true}) = gpu_init(u0, τ)
init(u0, τ, ::Val{false}) = cpu_init(u0, τ)

function create_semilinear_problem(
        nb_pts::Int,
        dτ::T,
        integration_limit::T,
        u0,
        semilinear_fourier_part::Function,
        nonlinear_part!::Function,
        use_gpu
) where {T <: Real}
    domain_size = dτ * nb_pts
    τ = Vector{T}(range(-domain_size / 2.0, domain_size / 2.0, nb_pts))

    u0_good_type = init(u0, τ, use_gpu)

    ω = 2π .* Vector{T}(FFTW.fftfreq(nb_pts, nb_pts / domain_size))
    diag = typeof(u0_good_type)(semilinear_fourier_part.(ω))

    D̂ = DiagonalFourierOperator(diag)
    N̂ = (du, u, _, t) -> nonlinear_part!(du, u, t)

    splfc = SplitFunction{true}(D̂, N̂)
    return SplitODEProblem(splfc, u0_good_type, (0.0, integration_limit))
end

@inline @fastmath function nlse_non_linpart(res, u)
    res .= u
    res .= abs2.(res)
    res .*= u
    res .*= 1im
end

@inline @fastmath function lle_non_linpart!(res, u, Δ, S)
    nlse_non_linpart(res, u)
    LinearAlgebra.axpby!(-1im * Δ, u, 1, res)
    res .+= S
end

function create_nlse(
        nb_pts::Int, dτ::T, integration_limit::T, u0, use_gpu) where {T <: Real}
    create_semilinear_problem(nb_pts, dτ, integration_limit, u0, ω -> -0.5im * ω^2,
        (du, u, t) -> nlse_non_linpart(du, u), use_gpu)
end

function create_lle(
        nb_pts::Int, dτ::T, integration_limit::T, u0, use_gpu, Δ::T, S::T) where {T <: Real}
    create_semilinear_problem(nb_pts, dτ, integration_limit, u0, ω -> -1 - 1im * ω^2,
        (du, u, t) -> lle_non_linpart!(du, u, Δ, S), use_gpu)
end

function create_lle_scan(nb_pts::Int, dτ::T, scan_time::T, u0, use_gpu,
        Δ_min::T, Δ_max::T, S::T) where {T <: Real}
    create_semilinear_problem(nb_pts, dτ, scan_time, u0, ω -> -1 - 1im * ω^2,
        (du, u, t) -> lle_non_linpart!(du, u, Δ_min + (Δ_max - Δ_min) * t / scan_time, S), use_gpu)
end

"""
Test function of the NLSE on a conservative soliton
"""
function nlse_test(
        ::Type{T}, use_gpu; B = 5.0, nb_pts = 512, dτ = 5e-3, propag_dist = 10.0,
        save_everystep = false, reltol = 1e-6) where {T <: Real}
    B = T(B)

    problem = create_nlse(
        nb_pts,
        T(dτ),
        T(propag_dist),
        τ -> B * sech(B * τ),
        use_gpu
    )
    return solve(problem, RKIP(T(1e-4), T(1.0)); save_everystep, reltol = T(reltol))
end

function nlse_test(type::Type{T}; kwargs...) where {T <: Real}
    nlse_test(type, Val(false); kwargs...)
end
nlse_test(use_gpu::Val{true}; kwargs...) = nlse_test(Float64, use_gpu; kwargs...)
nlse_test(use_gpu::Val{false}; kwargs...) = nlse_test(Float64, use_gpu; kwargs...)
nlse_test(; kwargs...) = nlse_test(Float64; kwargs...)

"""
Lugiato-Lefever equation test. Useful for benchmark due to its particular dynamic
"""
function lle_scan_test(
        ::Type{T}, use_gpu; nb_pts = 1024, dτ = 5e-2, save_everystep = false,
        reltol = 1e-6, S = 3, Δ_min = -2.0, Δ_max = 12.0, t_max = 10) where {T <: Real}
    problem = create_lle_scan(
        nb_pts,
        T(dτ),
        T(t_max),
        1e-1 * exp.(2im * π * rand(T, nb_pts)),
        use_gpu,
        T(Δ_min),
        T(Δ_max),
        T(S)
    )
    return solve(problem, RKIP(T(1e-4), T(1.0)); save_everystep, reltol = T(reltol))
end

function lle_scan_test(type::Type{T}; kwargs...) where {T <: Real}
    lle_scan_test(type, Val(false); kwargs...)
end
lle_scan_test(use_gpu::Val{true}; kwargs...) = lle_scan_test(Float64, use_gpu; kwargs...)
lle_scan_test(use_gpu::Val{false}; kwargs...) = lle_scan_test(Float64, use_gpu; kwargs...)
lle_scan_test(; kwargs...) = lle_scan_test(Float64; kwargs...)
