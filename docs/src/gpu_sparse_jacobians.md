# Stiff PDEs on the GPU with a Sparse Jacobian

Solving a stiff PDE on a GPU with an implicit method needs three things to line up:
the right-hand side has to execute on the device, the Jacobian prototype has to be a
GPU sparse matrix, and the linear solver has to be one that works on that type. Getting
any of them wrong produces the same message, `Scalar indexing is disallowed`, so this
page works through a complete example and explains what each piece is for.

## A complete example

The Brusselator, solved on a GPU with a sparse Jacobian. Every piece below is explained
in the sections that follow.

```julia
using OrdinaryDiffEq, LinearAlgebra, SparseArrays
using CUDA, CUDA.CUSPARSE
using CUDSS
using SciMLBase: FullSpecialize

CUDA.allowscalar(false)

const N = 32
const xyd = range(0.0, stop = 1.0, length = N)
const dx = step(xyd)
const II = LinearIndices((N, N, 2))

brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0
limit(a, N) = a == N + 1 ? 1 : a == 0 ? N : a

@inline function kernel_u!(du, u, A, B, α, II, I, t)
    i, j = Tuple(I)
    ip1 = limit(i + 1, N); im1 = limit(i - 1, N)
    jp1 = limit(j + 1, N); jm1 = limit(j - 1, N)
    du[II[i, j, 1]] = α * (u[II[im1, j, 1]] + u[II[ip1, j, 1]] + u[II[i, jp1, 1]] +
                           u[II[i, jm1, 1]] - 4u[II[i, j, 1]]) +
        B + u[II[i, j, 1]]^2 * u[II[i, j, 2]] - (A + 1) * u[II[i, j, 1]] +
        brusselator_f(xyd[i], xyd[j], t)
    return nothing
end

@inline function kernel_v!(du, u, A, B, α, II, I, t)
    i, j = Tuple(I)
    ip1 = limit(i + 1, N); im1 = limit(i - 1, N)
    jp1 = limit(j + 1, N); jm1 = limit(j - 1, N)
    du[II[i, j, 2]] = α * (u[II[im1, j, 2]] + u[II[ip1, j, 2]] + u[II[i, jp1, 2]] +
                           u[II[i, jm1, 2]] - 4u[II[i, j, 2]]) +
        A * u[II[i, j, 1]] - u[II[i, j, 1]]^2 * u[II[i, j, 2]]
    return nothing
end

# The index set is a parameter so the same maths can run on either device.
function make_brusselator(idx)
    return function (du, u, p, t)
        A, B, α = p[1], p[2], p[3] / dx^2
        kernel_u!.(Ref(du), Ref(u), A, B, α, Ref(II), idx, t)
        kernel_v!.(Ref(du), Ref(u), A, B, α, Ref(II), idx, t)
        return nothing
    end
end

# The stencil is known, so the sparsity pattern can be written down directly.
function brusselator_pattern()
    rows, cols = Int[], Int[]
    for j in 1:N, i in 1:N
        ip1 = limit(i + 1, N); im1 = limit(i - 1, N)
        jp1 = limit(j + 1, N); jm1 = limit(j - 1, N)
        for c in 1:2
            row = II[i, j, c]
            for (a, b) in ((im1, j), (ip1, j), (i, jp1), (i, jm1), (i, j))
                push!(rows, row); push!(cols, II[a, b, c])
            end
            push!(rows, row); push!(cols, II[i, j, 3 - c])
        end
    end
    return sparse(rows, cols, ones(length(rows)), 2N^2, 2N^2)
end

u0 = vec([c == 1 ? 22 * (xyd[j] * (1 - xyd[j]))^(3 / 2) :
                   27 * (xyd[i] * (1 - xyd[i]))^(3 / 2)
          for i in 1:N, j in 1:N, c in 1:2])
p = (3.4, 1.0, 10.0)
tspan = (0.0, 1.5)

u0_gpu = CuArray(u0)
gpu_rhs! = make_brusselator(CuArray(CartesianIndices((N, N))))
jac_prototype = CuSparseMatrixCSR(CuSparseMatrixCSC(brusselator_pattern()))

prob = ODEProblem(
    ODEFunction{true, FullSpecialize}(gpu_rhs!; jac_prototype),
    u0_gpu, tspan, p
)
sol = solve(prob, Rosenbrock23(); abstol = 1.0e-8, reltol = 1.0e-8)
```

## The right-hand side has to run on the device

This is the part that catches people out, and it has nothing to do with the solver. A
kernel written as a broadcast over an index set only runs on the GPU if the index set
itself lives on the GPU:

```julia
const II = LinearIndices((N, N, 2))

function make_brusselator(idx)
    return function (du, u, p, t)
        A, B, α = p[1], p[2], p[3] / dx^2
        kernel_u!.(Ref(du), Ref(u), A, B, α, Ref(II), idx, t)
        kernel_v!.(Ref(du), Ref(u), A, B, α, Ref(II), idx, t)
        return nothing
    end
end

cpu_rhs! = make_brusselator(CartesianIndices((N, N)))
gpu_rhs! = make_brusselator(CuArray(CartesianIndices((N, N))))
```

Every array argument above is wrapped in `Ref`, which makes it a scalar as far as
broadcasting is concerned. The index set is therefore the only broadcastable argument,
and it alone decides where the loop runs. Pass a host `CartesianIndices` and the whole
right-hand side executes on the CPU and reads the device arrays one element at a time,
which is what raises `Scalar indexing is disallowed` before the solver has done
anything at all.

It is worth checking this directly rather than inferring it from whether the solve
works:

```julia
du_cpu = similar(u0);  cpu_rhs!(du_cpu, u0, p, 0.0)
du_gpu = similar(u0_gpu); gpu_rhs!(du_gpu, u0_gpu, p, 0.0)
@assert Array(du_gpu) ≈ du_cpu
```

`mul!` against a `CuSparseMatrixCSR` is a second way to write a right-hand side that
cannot run here. It has no `ForwardDiff.Dual` kernel, so it falls back to a generic
implementation that indexes scalars as soon as the Jacobian is taken, even though the
same call is fine for plain `Float64`.

## The Jacobian prototype has to be a GPU sparse matrix

The prototype is used as the template for the Jacobian that gets built, so it has to be
the array type you want back:

```julia
using CUDA, CUDA.CUSPARSE

jac_prototype = CuSparseMatrixCSR(CuSparseMatrixCSC(pattern))
f = ODEFunction{true, FullSpecialize}(gpu_rhs!; jac_prototype)
prob = ODEProblem(f, u0_gpu, tspan, p)
```

`FullSpecialize` is used because the default wrapping goes through FunctionWrappers,
whose compiled signature is too narrow for some GPU solver caches.

If you have a sparsity pattern on the host, convert it rather than passing it directly;
a host `SparseMatrixCSC` prototype will produce a host Jacobian and the solve will fall
back to scalar indexing.

## The linear solver

Load [CUDSS.jl](https://github.com/exanauts/CUDSS.jl) for sparse LU on the GPU. Without
it, LinearSolve warns and falls back to a Krylov method:

```julia
using CUDSS

solve(prob, Rosenbrock23())                                        # sparse LU via CUDSS
solve(prob, Rosenbrock23(linsolve = KrylovJL_GMRES(), concrete_jac = true))
solve(prob, FBDF())
```

For a matrix-free preconditioner on the GPU, use
[KrylovPreconditioners.jl](https://github.com/JuliaSmoothOptimizers/KrylovPreconditioners.jl)
(`kp_ilu0`, `kp_ic0`) rather than passing a raw `CUSPARSE.ilu02` result, which stores
the factors but cannot apply them.

## Diagnosing `Scalar indexing is disallowed`

Every mistake above reports the same error, so read the stack trace to tell them apart:

  - frames inside your own kernel or `Base.Broadcast` mean the right-hand side is
    running on the host, and the index set is the thing to check
  - frames inside `jacobian!` or `calc_J!` mean the Jacobian is being built into a host
    matrix, so check the prototype's type
  - frames inside `lu_instance` mean CUDSS is not loaded

Run with `CUDA.allowscalar(false)` so these fail loudly instead of running slowly.
