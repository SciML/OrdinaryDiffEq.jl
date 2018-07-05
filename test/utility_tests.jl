using OrdinaryDiffEq: phi, phi, phiv, phiv_timestep, expv, expv_timestep, arnoldi, getH, getV
using LinearAlgebra

@testset "Exponential Utilities" begin
  # Scalar phi
  K = 4
  z = 0.1
  P = zeros(K+1); P[1] = exp(z)
  for i = 1:K
    P[i+1] = (P[i] - 1/factorial(i-1))/z
  end
  @test phi(z, K) ≈ P

  # Matrix phi
  A = [0.1 0.2; 0.3 0.4]
  P = Vector{Matrix{Float64}}(undef, K+1); P[1] = exp(A)
  for i = 1:K
    P[i+1] = (P[i] - 1/factorial(i-1)*I) / A
  end
  @test phi(A, K) ≈ P

  # Krylov
  n = 20; m = 5
  srand(0)
  A = randn(n, n)
  t = 1e-2
  b = randn(n)
  @test exp(t * A) * b ≈ expv(t, A, b; m=m)
  P = phi(t * A, K)
  W = zeros(n, K+1)
  for i = 1:K+1
    W[:,i] = P[i] * b
  end
  Ks = arnoldi(A, b; m=m)
  W_approx = phiv(t, Ks, K)
  @test W ≈ W_approx

  # Happy-breakdown in Krylov
  v = normalize(randn(n))
  A = v * v' # A is Idempotent
  Ks = arnoldi(A, b)
  @test Ks.m == 2

  # Arnoldi vs Lanczos
  A = Hermitian(randn(n, n))
  Aperm = A + 1e-10 * randn(n, n) # no longer Hermitian
  w = expv(t, A, b; m=m)
  wperm = expv(t, Aperm, b; m=m)
  @test w ≈ wperm

  # Internal time-stepping for Krylov (with adaptation)
  n = 100
  K = 4
  t = 5.0
  tol = 1e-7
  A = spdiagm(-1=>ones(n-1), 0=>-2*ones(n), 1=>ones(n-1))
  B = randn(n, K+1)
  Phi_half = phi(t/2 * A, K)
  Phi = phi(t * A, K)
  uhalf_exact = sum((t/2)^i * Phi_half[i+1] * B[:,i+1] for i = 0:K)
  u_exact = sum(t^i * Phi[i+1] * B[:,i+1] for i = 0:K)
  U = phiv_timestep([t/2, t], A, B; adaptive=true, tol=tol)
  @test norm(U[:,1] - uhalf_exact) / norm(uhalf_exact) < tol
  @test norm(U[:,2] - u_exact) / norm(u_exact) < tol
  # p = 0 special case (expv_timestep)
  u_exact = Phi[1] * B[:, 1]
  u = expv_timestep(t, A, B[:, 1]; adaptive=true, tol=tol)
  @test norm(u - u_exact) / norm(u_exact) < tol
end
