# Only diffusion okay.

using SpecialMatrices
A = -full(Strang(11))
g(t,u) = 2-u
u = zeros(11);u[6]=1
nsteps = 10000
tmax = 1.0
h = tmax/nsteps
t = 0

# Lawson-Euler
u = zeros(11);u[6]=1
for k in 1:nsteps
  u = expm(h*A)*(u + h*g(t,u))
  t = k*h
end
@show u

# Norsett-Euler
u = zeros(11);u[6]=1
for k in 1:nsteps
  u = expm(h*A)*u + ((expm(h*A)-I)/A)*g(t,u)
  t = k*h
end
@show u

# Norsett-Euler
u = zeros(11);u[6]=1
for k in 1:nsteps
  u = (I + A*(expm(h*A)-I)/A)*u + ((expm(h*A)-I)/A)*g(t,u)
  t = k*h
end
@show u

# Norsett-Euler
u = zeros(11);u[6]=1
for k in 1:nsteps
  u = u + ((expm(h*A)-I)/A)*(A*u + g(t,u))
  t = k*h
end
@show u


# Norsett-Euler
u = zeros(11);u[6]=1
for k in 1:nsteps
  u = u + ((expm(h*A)-I)/A)*(A*u + g(t,u))
  t = k*h
end
@show u

#  Euler
u = zeros(11);u[6]=1
for k in 1:nsteps
    u += h*(A*u + g(t,u))
    t = k*h
end
@show u
