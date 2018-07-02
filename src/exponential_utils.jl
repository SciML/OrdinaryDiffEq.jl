# exponential_utils.jl
# Contains functions related to the evaluation of scalar/matrix phi functions 
# that are used by the exponential integrators.
#
# TODO: write a version of `expm!` that is non-allocating.

###################################################
# Dense algorithms
const exp! = Base.LinAlg.expm! # v0.7 style

"""
    phi(z,k[;cache]) -> [phi_0(z),phi_1(z),...,phi_k(z)]

Compute the scalar phi functions for all orders up to k.

The phi functions are defined as

```math
\\varphi_0(z) = \\exp(z),\\quad \\varphi_k(z+1) = \\frac{\\varphi_k(z) - 1}{z} 
```

Instead of using the recurrence relation, which is numerically unstable, a 
formula given by Sidje is used (Sidje, R. B. (1998). Expokit: a software 
package for computing matrix exponentials. ACM Transactions on Mathematical 
Software (TOMS), 24(1), 130-156. Theorem 1).
"""
function phi(z::T, k::Integer; cache=nothing) where {T <: Number}
  # Construct the matrix
  if cache == nothing
    cache = zeros(T, k+1, k+1)
  else
    fill!(cache, zero(T))
  end
  cache[1,1] = z
  for i = 1:k
    cache[i,i+1] = one(T)
  end
  P = exp!(cache)
  return P[1,:]
end

"""
    phiv_dense(A,v,k[;cache]) -> [phi_0(A)v phi_1(A)v ... phi_k(A)v]

Compute the matrix-phi-vector products for small, dense `A`. `k`` >= 1.

The phi functions are defined as

```math
\\varphi_0(z) = \\exp(z),\\quad \\varphi_k(z+1) = \\frac{\\varphi_k(z) - 1}{z} 
```

Instead of using the recurrence relation, which is numerically unstable, a 
formula given by Sidje is used (Sidje, R. B. (1998). Expokit: a software 
package for computing matrix exponentials. ACM Transactions on Mathematical 
Software (TOMS), 24(1), 130-156. Theorem 1).
"""
function phiv_dense(A, v, k; cache=nothing)
  w = Matrix{eltype(A)}(length(v), k+1)
  phiv_dense!(w, A, v, k; cache=cache)
end
"""
    phiv_dense!(w,A,v,k[;cache]) -> w

Non-allocating version of `phiv_dense`.
"""
function phiv_dense!(w::AbstractMatrix{T}, A::AbstractMatrix{T}, 
  v::AbstractVector{T}, k::Integer; cache=nothing) where {T <: Number}
  @assert size(w, 1) == size(A, 1) == size(A, 2) == length(v) "Dimension mismatch"
  @assert size(w, 2) == k+1 "Dimension mismatch"
  m = length(v)
  # Construct the extended matrix
  if cache == nothing
    cache = zeros(T, m+k, m+k)
  else
    @assert size(cache) == (m+k, m+k) "Dimension mismatch"
    fill!(cache, zero(T))
  end
  cache[1:m, 1:m] = A
  cache[1:m, m+1] = v
  for i = m+1:m+k-1
    cache[i, i+1] = one(T)
  end
  P = exp!(cache)
  # Extract results
  @views A_mul_B!(w[:, 1], P[1:m, 1:m], v)
  @inbounds for i = 1:k
    @inbounds for j = 1:m
      w[j, i+1] = P[j, m+i]
    end
  end
  return w
end

"""
    phi(A,k[;cache]) -> [phi_0(A),phi_1(A),...,phi_k(A)]

Compute the matrix phi functions for all orders up to k. `k` >= 1.

The phi functions are defined as
  
```math
\\varphi_0(z) = \\exp(z),\\quad \\varphi_k(z+1) = \\frac{\\varphi_k(z) - 1}{z} 
```

Calls `phiv_dense` on each of the basis vectors to obtain the answer.
"""
function phi(A::AbstractMatrix{T}, k; caches=nothing) where {T <: Number}
  m = size(A, 1)
  out = [Matrix{T}(m, m) for i = 1:k+1]
  phi!(out, A, k; caches=caches)
end
"""
    phi!(out,A,k[;caches]) -> out

Non-allocating version of `phi` for matrix inputs.
"""
function phi!(out::Vector{Matrix{T}}, A::AbstractMatrix{T}, k::Integer; caches=nothing) where {T <: Number}
  m = size(A, 1)
  @assert length(out) == k + 1 && all(P -> size(P) == (m,m), out) "Dimension mismatch"
  if caches == nothing
    e = Vector{T}(m)
    W = Matrix{T}(m, k+1)
    C = Matrix{T}(m+k, m+k)
  else
    e, W, C = caches
    @assert size(e) == (m,) && size(W) == (m, k+1) && size(C) == (m+k, m+k) "Dimension mismatch"
  end
  @inbounds for i = 1:m
    fill!(e, zero(T)); e[i] = one(T) # e is the ith basis vector
    phiv_dense!(W, A, e, k; cache=C) # W = [phi_0(A)*e phi_1(A)*e ... phi_k(A)*e]
    @inbounds for j = 1:k+1
      @inbounds for s = 1:m
        out[j][s, i] = W[s, j]
      end
    end
  end
  return out
end

##############################################
# Krylov algorithms
"""
    KrylovSubspace{T}(n,[maxiter=30]) -> Ks

Constructs an uninitialized Krylov subspace, which can be filled by `arnoldi!`.

The dimension of the subspace, `Ks.m`, can be dynamically altered but should 
be smaller than `maxiter`, the maximum allowed arnoldi iterations.

    getV(Ks) -> V
    getH(Ks) -> H

Access methods for the (extended) orthonormal basis `V` and the (extended) 
Gram-Schmidt coefficients `H`. Both methods return a view into the storage 
arrays and has the correct dimensions as indicated by `Ks.m`.

    resize!(Ks, maxiter) -> Ks

Resize `Ks` to a different `maxiter`, destroying its contents.

This is an expensive operation and should be used scarsely.
"""
mutable struct KrylovSubspace{B, T}
  m::Int        # subspace dimension
  maxiter::Int  # maximum allowed subspace size
  beta::B       # norm(b,2)
  V::Matrix{T}  # orthonormal bases
  H::Matrix{T}  # Gram-Schmidt coefficients
  KrylovSubspace{T}(n::Integer, maxiter::Integer=30) where {T} = new{real(T), T}(
    maxiter, maxiter, zero(real(T)), Matrix{T}(n, maxiter + 1), 
    zeros(T, maxiter + 1, maxiter))
end
# TODO: switch to overload `getproperty` in v0.7
getH(Ks::KrylovSubspace) = @view(Ks.H[1:Ks.m + 1, 1:Ks.m])
getV(Ks::KrylovSubspace) = @view(Ks.V[:, 1:Ks.m + 1])
function Base.resize!(Ks::KrylovSubspace{B,T}, maxiter::Integer) where {B,T}
  V = Matrix{T}(size(Ks.V, 1), maxiter + 1)
  H = zeros(T, maxiter + 1, maxiter)
  Ks.V = V; Ks.H = H
  Ks.m = Ks.maxiter = maxiter
  return Ks
end
function Base.show(io::IO, Ks::KrylovSubspace)
  println(io, "$(Ks.m)-dimensional Krylov subspace with fields")
  println(io, "beta: $(Ks.beta)")
  print(io, "V: ")
  println(IOContext(io, limit=true), getV(Ks))
  print(io, "H: ")
  println(IOContext(io, limit=true), getH(Ks))
end

"""
    arnoldi(A,b[;m,tol,norm,cache]) -> Ks

Performs `m` anoldi iterations to obtain the Krylov subspace K_m(A,b).

The n x (m + 1) basis vectors `getV(Ks)` and the (m + 1) x m upper Heisenberg 
matrix `getH(Ks)` are related by the recurrence formula

```
v_1=b,\\quad Av_j = \\sum_{i=1}^{j+1}h_{ij}v_i\\quad(j = 1,2,\\ldots,m)
```

`iop` determines the length of the incomplete orthogonalization procedure [^1]. 
The default value of 0 indicates full Arnoldi. For symmetric/Hermitian `A`, 
`iop` will be ignored and the Lanczos algorithm will be used instead.

Refer to `KrylovSubspace` for more information regarding the output.

Happy-breakdown occurs whenver `norm(v_j) < tol * norm(A, Inf)`, in this case 
the dimension of `Ks` is smaller than `m`.

[^1]: Koskela, A. (2015). Approximating the matrix exponential of an 
advection-diffusion operator using the incomplete orthogonalization method. In 
Numerical Mathematics and Advanced Applications-ENUMATH 2013 (pp. 345-353). 
Springer, Cham.
"""
function arnoldi(A, b; m=min(30, size(A, 1)), tol=1e-7, norm=Base.norm, 
  iop=0, cache=nothing)
  Ks = KrylovSubspace{eltype(b)}(length(b), m)
  arnoldi!(Ks, A, b; m=m, tol=tol, norm=norm, cache=cache, iop=iop)
end
"""
    arnoldi!(Ks,A,b[;tol,m,norm,cache]) -> Ks

Non-allocating version of `arnoldi`.
"""
function arnoldi!(Ks::KrylovSubspace{B, T}, A, b::AbstractVector{T}; tol::Real=1e-7, 
  m::Int=min(Ks.maxiter, size(A, 1)), norm=Base.norm, iop::Int=0, cache=nothing) where {B, T <: Number}
  if ishermitian(A)
    return lanczos!(Ks, A, b; tol=tol, m=m, norm=norm, cache=cache)
  end
  if m > Ks.maxiter
    resize!(Ks, m)
  else
    Ks.m = m # might change if happy-breakdown occurs
  end
  V, H = getV(Ks), getH(Ks)
  vtol = tol * norm(A, Inf)
  if iop == 0
    iop = m
  end
  # Safe checks
  n = size(V, 1)
  @assert length(b) == size(A,1) == size(A,2) == n "Dimension mismatch"
  if cache == nothing
    cache = similar(b)
  else
    @assert size(cache) == (n,) "Dimension mismatch"
  end
  # Arnoldi iterations (with IOP)
  fill!(H, zero(T))
  Ks.beta = norm(b)
  @. V[:, 1] = b / Ks.beta
  @inbounds for j = 1:m
    A_mul_B!(cache, A, @view(V[:, j]))
    @inbounds for i = max(1, j - iop + 1):j
      alpha = dot(@view(V[:, i]), cache)
      H[i, j] = alpha
      Base.axpy!(-alpha, @view(V[:, i]), cache)
    end
    beta = norm(cache)
    H[j+1, j] = beta
    @inbounds for i = 1:n
      @. V[i, j+1] = cache[i] / beta
    end
    if beta < vtol # happy-breakdown
      Ks.m = j
      break
    end
  end
  return Ks
end
"""
    lanczos!(Ks,A,b[;tol,m,norm,cache]) -> Ks

A variation of `arnoldi!` that uses the Lanczos algorithm for Hermitian matrices.
"""
function lanczos!(Ks::KrylovSubspace{B, T}, A, b::AbstractVector{T}; tol=1e-7,
  m=min(Ks.maxiter, size(A, 1)), norm=Base.norm, cache=nothing) where {B, T <: Number}
  if m > Ks.maxiter
    resize!(Ks, m)
  else
    Ks.m = m # might change if happy-breakdown occurs
  end
  V, H = getV(Ks), getH(Ks)
  vtol = tol * norm(A, Inf)
  # Safe checks
  n = size(V, 1)
  @assert length(b) == size(A,1) == size(A,2) == n "Dimension mismatch"
  if cache == nothing
    cache = similar(b)
  else
    @assert size(cache) == (n,) "Dimension mismatch"
  end
  # Lanczos iterations
  fill!(H, zero(T))
  Ks.beta = norm(b)
  @. V[:, 1] = b / Ks.beta
  @inbounds for j = 1:m
    vj = @view(V[:, j])
    A_mul_B!(cache, A, vj)
    alpha = dot(vj, cache)
    H[j, j] = alpha
    Base.axpy!(-alpha, vj, cache)
    if j > 1
      Base.axpy!(-H[j-1, j], @view(V[:, j-1]), cache)
    end
    beta = norm(cache)
    H[j+1, j] = beta
    if j < m
      H[j, j+1] = beta
    end
    @inbounds for i = 1:n
      V[i, j+1] = cache[i] / beta
    end
    if beta < vtol # happy-breakdown
      Ks.m = j
      break
    end
  end
  return Ks
end

# Cache type for expv
mutable struct ExpvCache{T}
  mem::Vector{T}
  ExpvCache{T}(maxiter::Int) where {T} = new{T}(Vector{T}(maxiter^2))
end
function Base.resize!(C::ExpvCache{T}, maxiter::Int) where {T}
  C.mem = Vector{T}(maxiter^2 * 2)
  return C
end
function get_cache(C::ExpvCache, m::Int)
  m^2 > length(C.mem) && resize!(C, m) # resize the cache if needed
  reshape(@view(C.mem[1:m^2]), m, m)
end
"""
    expv(t,A,b; kwargs) -> exp(tA)b

Compute the matrix-exponential-vector product using Krylov.

A Krylov subspace is constructed using `arnoldi` and `expm!` is called 
on the Heisenberg matrix. Consult `arnoldi` for the values of the keyword 
arguments.

    expv(t,Ks; cache) -> exp(tA)b

Compute the expv product using a pre-constructed Krylov subspace.
"""
function expv(t, A, b; m=min(30, size(A, 1)), tol=1e-7, norm=Base.norm, cache=nothing, iop=0)
  Ks = arnoldi(A, b; m=m, tol=tol, norm=norm, iop=iop)
  w = similar(b)
  expv!(w, t, Ks; cache=cache)
end
function expv(t, Ks::KrylovSubspace{B, T}; cache=nothing) where {B, T}
  n = size(getV(Ks), 1)
  w = Vector{T}(n)
  expv!(w, t, Ks; cache=cache)
end
"""
    expv!(w,t,Ks[;cache]) -> w

Non-allocating version of `expv` that uses precomputed Krylov subspace `Ks`.
"""
function expv!(w::Vector{T}, t::Number, Ks::KrylovSubspace{B, T}; 
  cache=nothing) where {B, T <: Number}
  m, beta, V, H = Ks.m, Ks.beta, getV(Ks), getH(Ks)
  @assert length(w) == size(V, 1) "Dimension mismatch"
  if cache == nothing
    cache = Matrix{T}(m, m)
  elseif isa(cache, ExpvCache)
    cache = get_cache(cache, m)
  else
    throw(ArgumentError("Cache must be an ExpvCache"))
  end
  scale!(t, copy!(cache, @view(H[1:m, :])))
  if ishermitian(cache)
    # Optimize the case for symtridiagonal H
    F = eigfact!(SymTridiagonal(cache)) # Note: eigfact! -> eigen! in v0.7
    expHe = F.vectors * (exp.(F.values) .* @view(F.vectors[1, :]))
  else
    expH = exp!(cache)
    expHe = @view(expH[:, 1])
  end
  scale!(beta, A_mul_B!(w, @view(V[:, 1:m]), expHe)) # exp(A) ≈ norm(b) * V * exp(H)e
end

# Cache type for phiv
mutable struct PhivCache{T}
  mem::Vector{T}
  function PhivCache{T}(maxiter::Int, p::Int) where {T}
    numelems = maxiter + maxiter^2 + (maxiter + p)^2 + maxiter*(p + 1)
    new{T}(Vector{T}(numelems))
  end
end
function Base.resize!(C::PhivCache{T}, maxiter::Int, p::Int) where {T}
  numelems = maxiter + maxiter^2 + (maxiter + p)^2 + maxiter*(p + 1)
  C.mem = Vector{T}(numelems * 2)
  return C
end
function get_caches(C::PhivCache, m::Int, p::Int)
  numelems = m + m^2 + (m + p)^2 + m*(p + 1)
  numelems^2 > length(C.mem) && resize!(C, m, p) # resize the cache if needed
  e = @view(C.mem[1:m]); offset = m
  Hcopy = reshape(@view(C.mem[offset + 1:offset + m^2]), m, m); offset += m^2
  C1 = reshape(@view(C.mem[offset + 1:offset + (m+p)^2]), m+p, m+p); offset += (m+p)^2
  C2 = reshape(@view(C.mem[offset + 1:offset + m*(p+1)]), m, p+1)
  return e, Hcopy, C1, C2
end
"""
    phiv(t,A,b,k;correct,kwargs) -> [phi_0(tA)b phi_1(tA)b ... phi_k(tA)b][, errest]

Compute the matrix-phi-vector products using Krylov. `k` >= 1.

The phi functions are defined as

```math
\\varphi_0(z) = \\exp(z),\\quad \\varphi_k(z+1) = \\frac{\\varphi_k(z) - 1}{z} 
```

A Krylov subspace is constructed using `arnoldi` and `phiv_dense` is called 
on the Heisenberg matrix. If `correct=true`, then phi_0 through phi_k-1 are 
updated using the last Arnoldi vector v_m+1 [^1]. If `errest=true` then an 
additional error estimate for the second-to-last phi is also returned. For 
the additional keyword arguments, consult `arnoldi`.

  phiv(t,Ks,k;correct,kwargs) -> [phi_0(tA)b phi_1(tA)b ... phi_k(tA)b][, errest]

Compute the matrix-phi-vector products using a pre-constructed Krylov subspace.

[^1]: Niesen, J., & Wright, W. (2009). A Krylov subspace algorithm for evaluating 
the φ-functions in exponential integrators. arXiv preprint arXiv:0907.4631. 
Formula (10).
"""
function phiv(t, A, b, k; m=min(30, size(A, 1)), tol=1e-7, norm=Base.norm, iop=0, 
  cache=nothing, correct=false, errest=false)
  Ks = arnoldi(A, b; m=m, tol=tol, norm=norm, iop=iop)
  w = Matrix{eltype(b)}(length(b), k+1)
  phiv!(w, t, Ks, k; cache=cache, correct=correct, errest=errest)
end
function phiv(t, Ks::KrylovSubspace{B, T}, k; cache=nothing, correct=false, 
  errest=false) where {B, T}
  n = size(getV(Ks), 1)
  w = Matrix{T}(n, k+1)
  phiv!(w, t, Ks, k; cache=cache, correct=correct, errest=errest)
end
"""
    phiv!(w,t,Ks,k[;cache,correct,errest]) -> w[,errest]

Non-allocating version of 'phiv' that uses precomputed Krylov subspace `Ks`.
"""
function phiv!(w::Matrix{T}, t::Number, Ks::KrylovSubspace{B, T}, k::Integer; 
  cache=nothing, correct=false, errest=false) where {B, T <: Number}
  m, beta, V, H = Ks.m, Ks.beta, getV(Ks), getH(Ks)
  @assert size(w, 1) == size(V, 1) "Dimension mismatch"
  @assert size(w, 2) == k + 1 "Dimension mismatch"
  if cache == nothing
    cache = PhivCache{T}(m, k)
  elseif !isa(cache, PhivCache)
    throw(ArgumentError("Cache must be a PhivCache"))
  end
  e, Hcopy, C1, C2 = get_caches(cache, m, k)
  scale!(t, copy!(Hcopy, @view(H[1:m, :])))
  fill!(e, zero(T)); e[1] = one(T) # e is the [1,0,...,0] basis vector
  phiv_dense!(C2, Hcopy, e, k; cache=C1) # C2 = [ϕ0(H)e ϕ1(H)e ... ϕk(H)e]
  scale!(beta, A_mul_B!(w, @view(V[:, 1:m]), C2)) # f(A) ≈ norm(b) * V * f(H)e
  if correct
    # Use the last Arnoldi vector for correction with little additional cost
    # correct_p = beta * h_{m+1,m} * (em^T phi_p+1(H) e1) * v_m+1
    betah = beta * H[end,end] * t
    vlast = @view(V[:,end])
    @inbounds for i = 1:k
      Base.axpy!(betah * C2[end, i+1], vlast, @view(w[:, i]))
    end
  end
  if errest
    err = abs(beta * H[end, end] * t * C2[end, end])
    return w, err
  else
    return w
  end
end

###########################################
# Krylov phiv with internal time-stepping
"""
    exp_timestep(ts,A,b[;adaptive,tol,kwargs...]) -> U

Evaluates the matrix exponentiation-vector product using time stepping

```math
u = \\exp(tA)b 
```

`ts`` is an array of time snapshots for u, with `U[:,j] ≈ u(ts[j])`. `ts` can 
also be just one value, in which case only the end result is returned and `U` 
is a vector.

The time stepping formula of Niesen & Wright is used [^1]. If the time step 
`tau` is not specified, it is chosen according to (17) of Neisen & Wright. If 
`adaptive==true`, the time step and Krylov subsapce size adaptation scheme of 
Niesen & Wright is used, the relative tolerance of which can be set using the 
keyword parameter `tol`. The delta and gamma parameter of the adaptation 
scheme can also be adjusted.

Set `verbose=true` to print out the internal steps (for debugging). For the 
other keyword arguments, consult `arnoldi` and `phiv`, which are used 
internally.

Note that this function is just a special case of `phiv_timestep` with a more 
intuitive interface (vector `b` instead of a n-by-1 matrix `B`).

[^1]: Niesen, J., & Wright, W. (2009). A Krylov subspace algorithm for 
evaluating the φ-functions in exponential integrators. arXiv preprint 
arXiv:0907.4631.
"""
function expv_timestep(ts::Vector{tType}, A, b; kwargs...) where {tType <: Real}
  U = Matrix{eltype(A)}(size(A, 1), length(ts))
  expv_timestep!(U, ts, A, b; kwargs...)
end
function expv_timestep(t::tType, A, b; kwargs...) where {tType <: Real}
  u = Vector{eltype(A)}(size(A, 1))
  expv_timestep!(u, t, A, b; kwargs...)
end
"""
    expv_timestep!(u,t,A,b[;kwargs]) -> u

Non-allocating version of `expv_timestep`.
"""
function expv_timestep!(u::Vector{T}, t::tType, A, b::Vector{T}; 
  kwargs...) where {T <: Number, tType <: Real}
  expv_timestep!(reshape(u, length(u), 1), [t], A, b; kwargs...)
  return u
end
function expv_timestep!(U::Matrix{T}, ts::Vector{tType}, A, b::Vector{T}; 
  kwargs...) where {T <: Number, tType <: Real}
  B = reshape(b, length(b), 1)
  phiv_timestep!(U, ts, A, B; kwargs...)
end
"""
    phiv_timestep(ts,A,B[;adaptive,tol,kwargs...]) -> U

Evaluates the linear combination of phi-vector products using time stepping

```math
u = \\varphi_0(tA)b_0 + t\\varphi_1(tA)b_1 + \\cdots + t^p\\varphi_p(tA)b_p 
```

`ts`` is an array of time snapshots for u, with `U[:,j] ≈ u(ts[j])`. `ts` can 
also be just one value, in which case only the end result is returned and `U` 
is a vector.

The time stepping formula of Niesen & Wright is used [^1]. If the time step 
`tau` is not specified, it is chosen according to (17) of Neisen & Wright. If 
`adaptive==true`, the time step and Krylov subsapce size adaptation scheme of 
Niesen & Wright is used, the relative tolerance of which can be set using the 
keyword parameter `tol`. The delta and gamma parameter of the adaptation 
scheme can also be adjusted.

Set `verbose=true` to print out the internal steps (for debugging). For the 
other keyword arguments, consult `arnoldi` and `phiv`, which are used 
internally.

[^1]: Niesen, J., & Wright, W. (2009). A Krylov subspace algorithm for 
evaluating the φ-functions in exponential integrators. arXiv preprint 
arXiv:0907.4631.
"""
function phiv_timestep(ts::Vector{tType}, A, B; kwargs...) where {tType <: Real}
  U = Matrix{eltype(A)}(size(A, 1), length(ts))
  phiv_timestep!(U, ts, A, B; kwargs...)
end
function phiv_timestep(t::tType, A, B; kwargs...) where {tType <: Real}
  u = Vector{eltype(A)}(size(A, 1))
  phiv_timestep!(u, t, A, B; kwargs...)
end
"""
    phiv_timestep!(U,ts,A,B[;kwargs]) -> U

Non-allocating version of `phiv_timestep`.
"""
function phiv_timestep!(u::Vector{T}, t::tType, A, B::Matrix{T}; 
  kwargs...) where {T <: Number, tType <: Real}
  phiv_timestep!(reshape(u, length(u), 1), [t], A, B; kwargs...)
  return u
end
function phiv_timestep!(U::Matrix{T}, ts::Vector{tType}, A, B::Matrix{T}; tau::Real=0.0, 
  m::Int=min(10, size(A, 1)), tol::Real=1e-7, norm=Base.norm, iop::Int=0, 
  correct::Bool=false, caches=nothing, adaptive=false, delta::Real=1.2, 
  gamma::Real=0.8, NA::Int=0, verbose=false) where {T <: Number, tType <: Real}
  # Choose initial timestep
  abstol = tol * norm(A, Inf)
  verbose && println("Absolute tolerance: $abstol")
  if iszero(tau)
    Anorm = norm(A, Inf)
    b0norm = norm(@view(B[:, 1]), Inf)
    tau = 10/Anorm * (abstol * ((m+1)/e)^(m+1) * sqrt(2*pi*(m+1)) / 
      (4*Anorm*b0norm))^(1/m)
    verbose && println("Initial time step unspecified, chosen to be $tau")
  end
  # Initialization
  n = size(U, 1)
  sort!(ts); tend = ts[end]
  p = size(B, 2) - 1
  @assert length(ts) == size(U, 2) "Dimension mismatch"
  @assert n == size(A, 1) == size(A, 2) == size(B, 1) "Dimension mismatch"
  if caches == nothing
    u = Vector{T}(n)              # stores the current state
    W = Matrix{T}(n, p+1)         # stores the w vectors
    P = Matrix{T}(n, p+2)         # stores output from phiv!
    Ks = KrylovSubspace{T}(n, m)  # stores output from arnoldi!
    phiv_cache = nothing         # cache used by phiv!
  else
    u, W, P, Ks, phiv_cache = caches
    @assert length(u) == n && size(W) == (n, p+1) && size(P) == (n, p+2) "Dimension mismatch"
  end
  copy!(u, @view(B[:, 1])) # u(0) = b0
  coeffs = ones(tType, p);
  if adaptive # initialization step for the adaptive scheme
    if ishermitian(A)
      iop = 2 # does not have an effect on arnoldi!, just for flops estimation
    end
    if iszero(NA)
      if isa(A, SparseMatrixCSC)
        NA = nnz(A)
      else
        NA = countnz(A) # not constant operation, should be best avoided
      end
    end
  end

  t = 0.0       # current time
  snapshot = 1  # which snapshot to compute next
  while t < tend # time stepping loop
    if t + tau > tend # last step
      tau = tend - t
    end
    # Part 1: compute w0...wp using the recurrence relation (16)
    copy!(@view(W[:, 1]), u) # w0 = u(t)
    @inbounds for l = 1:p-1 # compute cl = t^l/l!
      coeffs[l+1] = coeffs[l] * t / l
    end
    @views @inbounds for j = 1:p
      A_mul_B!(W[:, j+1], A, W[:, j])
      for l = 0:p-j
        Base.axpy!(coeffs[l+1], B[:, j+l+1], W[:, j+1])
      end
    end
    # Part 2: compute ϕp(tau*A)wp using Krylov, possibly with adaptation
    arnoldi!(Ks, A, @view(W[:, end]); tol=tol, m=m, norm=norm, iop=iop, cache=u)
    _, epsilon = phiv!(P, tau, Ks, p + 1; cache=phiv_cache, correct=correct, errest=true)
    verbose && println("t = $t, m = $m, tau = $tau, error estimate = $epsilon")
    if adaptive
      omega = (tend / tau) * (epsilon / abstol)
      epsilon_old = epsilon; m_old = m; tau_old = tau
      q = m/4; kappa = 2.0; maxtau = tend - t
      while omega > delta # inner loop of Algorithm 3
        m_new, tau_new, q, kappa = _phiv_timestep_adapt(
          m, tau, epsilon, m_old, tau_old, epsilon_old, q, kappa, 
          gamma, omega, maxtau, n, p, NA, iop, norm(getH(Ks), 1), verbose)
        m, m_old = m_new, m
        tau, tau_old = tau_new, tau
        # Compute ϕp(tau*A)wp using the new parameters
        arnoldi!(Ks, A, @view(W[:, end]); tol=tol, m=m, norm=norm, iop=iop, cache=u)
        _, epsilon_new = phiv!(P, tau, Ks, p + 1; cache=phiv_cache, correct=correct, errest=true)
        epsilon, epsilon_old = epsilon_new, epsilon
        omega = (tend / tau) * (epsilon / abstol)
        verbose && println("  * m = $m, tau = $tau, error estimate = $epsilon")
      end
    end
    # Part 3: update u using (15)
    scale!(tau^p, copy!(u, @view(P[:, end - 1])))
    @inbounds for l = 1:p-1 # compute cl = tau^l/l!
      coeffs[l+1] = coeffs[l] * tau / l
    end
    @views @inbounds for j = 0:p-1
      Base.axpy!(coeffs[j+1], W[:, j+1], u)
    end
    # Fill out all snapshots in between the current step
    while snapshot <= length(ts) && t + tau >= ts[snapshot]
      tau_snapshot = ts[snapshot] - t
      u_snapshot = @view(U[:, snapshot])
      phiv!(P, tau_snapshot, Ks, p + 1; cache=phiv_cache, correct=correct)
      scale!(tau_snapshot^p, copy!(u_snapshot, @view(P[:, end - 1])))
      @inbounds for l = 1:p-1 # compute cl = tau^l/l!
        coeffs[l+1] = coeffs[l] * tau_snapshot / l
      end
      @views @inbounds for j = 0:p-1
        Base.axpy!(coeffs[j+1], W[:, j+1], u_snapshot)
      end
      snapshot += 1
    end

    t += tau
  end

  return U
end
# Helper functions for phiv_timestep!
function _phiv_timestep_adapt(m, tau, epsilon, m_old, tau_old, epsilon_old, q, kappa, 
  gamma, omega, maxtau, n, p, NA, iop, Hnorm, verbose)
  # Compute new m and tau (Algorithm 4)
  if tau_old > tau
    q = log(tau/tau_old) / log(epsilon/epsilon_old) - 1
  end # else keep q the same
  tau_new = tau * (gamma / omega)^(1/(q + 1))
  tau_new = min(max(tau_new, tau/5), 2*tau, maxtau)
  if m_old < m
    kappa = (epsilon/epsilon_old)^(1/(m_old - m))
  end # else keep kappa the same
  m_new = m + ceil(Int, log(omega / gamma) / log(kappa))
  m_new = min(max(m_new, div(3*m, 4), 1), Int(ceil(4*m / 3)))
  verbose && println("  - Proposed new m: $m_new, new tau: $tau_new")
  # Compare costs of using new m vs new tau (23)
  cost_tau = _phiv_timestep_estimate_flops(m, tau_new, n, p, NA, iop, Hnorm, maxtau)
  cost_m = _phiv_timestep_estimate_flops(m_new, tau, n, p, NA, iop, Hnorm, maxtau)
  verbose && println("  - Cost to use new m: $cost_m flops, new tau: $cost_tau flops")
  if cost_tau < cost_m
    m_new = m
  else
    tau_new = tau
  end
  return m_new, tau_new, q, kappa
end
function _phiv_timestep_estimate_flops(m, tau, n, p, NA, iop, Hnorm, maxtau)
  # Estimate flops for the update of W and u
  flops_W = 2 * (p - 1) * (NA + n)
  flops_u = (2 * p + 1) * n
  # Estimate flops for arnoldi!
  if iop == 0
    iop = m
  end
  flops_matvec = 2 * m * NA
  flops_vecvec = 0
  for i = 1:m
    flops_vecvec += 3 * min(i, iop)
  end
  # Estimate flops for phiv! (7)
  MH = 44/3 + 2 * ceil(max(0.0, log2(Hnorm / 5.37)))
  flops_phiv = round(Int, MH * (m + p)^3)

  flops_onestep = flops_W + flops_u + flops_matvec + flops_vecvec + flops_phiv
  return flops_onestep * Int(ceil(maxtau / tau))
end
function _phiv_timestep_caches(u_prototype, maxiter::Int, p::Int)
  n = length(u_prototype); T = eltype(u_prototype)
  u = similar(u_prototype)                      # stores the current state
  W = Matrix{T}(n, p+1)                         # stores the w vectors
  P = Matrix{T}(n, p+2)                         # stores output from phiv!
  Ks = KrylovSubspace{T}(n, maxiter)            # stores output from arnoldi!
  phiv_cache = PhivCache{T}(maxiter, p+1)       # cache used by phiv! (need +1 for error estimation)
  return u, W, P, Ks, phiv_cache
end
