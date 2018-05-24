using OrdinaryDiffEq: phi, phim, _phimv, _expmv, arnoldi

@testset "Phi functions" begin
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
  P = Vector{Matrix{Float64}}(K+1); P[1] = expm(A)
  for i = 1:K
    P[i+1] = (P[i] - 1/factorial(i-1)*I) / A
  end
  @test phim(A, K) ≈ P

  # Krylov
  n = 20; m = 5
  srand(0)
  A = randn(n, n)
  t = 1e-2
  b = randn(n)
  @test expm(t * A) * b ≈ _expmv(t, A, b; m=m)
  P = phim(t * A, K)
  W = zeros(n, K+1)
  for i = 1:K+1
    W[:,i] = P[i] * b
  end
  W_approx = _phimv(t, A, b, K; m=m)
  @test W ≈ W_approx

  # Happy-breakdown in Krylov
  v = normalize(randn(n))
  A = v * v' # A is Idempotent
  Ks = arnoldi(A, b)
  @test Ks.m == 2
end
