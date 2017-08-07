# http://www.chimica.unipd.it/antonino.polimeno/pubblica/downloads/JChemPhys_101_4062.pdf

function initialize!(integrator,cache::SymplecticEulerConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  # Do the calculation pre
  # So that way FSAL interpolation
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  integrator.k[2].x[2] = f.f2(integrator.t,uprev,duprev)
  @muladd du = @. duprev + integrator.dt*integrator.k[2].x[2]
  integrator.k[1].x[1] = f.f1(integrator.t,uprev,du)
end

@muladd function perform_step!(integrator,cache::SymplecticEulerConstantCache,f=integrator.f)
  @unpack t,dt = integrator
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  kuprev = integrator.k[1].x[1]
  ku  = integrator.k[2].x[1]
  kdu = integrator.k[2].x[2]
  u = @. uprev + dt*kuprev
  # Now actually compute the step
  # Do it at the end for interpolations!
  kdu = f.f2(t,uprev,duprev)
  du = @. duprev + dt*kdu
  ku = f.f1(t,uprev,du)
  integrator.fsallast = ku
end

function initialize!(integrator,cache::SymplecticEulerCache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  # Do the calculation pre
  # So that way FSAL interpolation
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku = integrator.k[1].x[1]
  kdu = integrator.k[2].x[2]
  f.f2(integrator.t,uprev,duprev,kdu)
  @muladd @. du = duprev + integrator.dt*kdu
  f.f1(integrator.t,uprev,du,ku)
end

@muladd function perform_step!(integrator,cache::SymplecticEulerCache,f=integrator.f)
  @unpack t,dt = integrator
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  kuprev = integrator.k[1].x[1]
  ku  = integrator.k[2].x[1]
  kdu = integrator.k[2].x[2]
  @. u = uprev + dt*kuprev
  # Now actually compute the step
  # Do it at the end for interpolations!
  f.f2(t,uprev,duprev,kdu)
  @. du = duprev + dt*kdu
  f.f1(t,uprev,du,ku)
end

function initialize!(integrator,cache::C,f=integrator.f) where
    {C<:Union{VelocityVerletCache,Symplectic2Cache,Symplectic3Cache,Symplectic4Cache,
              Symplectic45Cache,Symplectic5Cache,Symplectic6Cache,Symplectic62Cache,
              McAte8Cache,KahanLi8Cache,SofSpa10Cache}}
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k

  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast

  uprev,duprev = integrator.uprev.x
  f.f1(integrator.t,uprev,duprev,integrator.k[2].x[1])
  f.f2(integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@muladd function perform_step!(integrator,cache::VelocityVerletCache,f=integrator.f)
  @unpack t,dt = integrator
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # x(t+Δt) = x(t) + v(t)*Δt + 1/2*a(t)*Δt^2
  f.f2(t,uprev,duprev,ku)
  @. u = uprev + dt*duprev + dt^2*(1//2*ku)
  f.f2(t+dt,u,duprev,kdu)
  # v(t+Δt) = v(t) + 1/2*(a(t)+a(t+Δt))*Δt
  @. du = duprev + dt*(1//2*ku + 1//2*kdu)
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@muladd function perform_step!(integrator,cache::Symplectic2Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,b1,b2 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @. u = uprev + dt*b1*duprev
  # update velocity
  f.f2(t,u,duprev,kdu)
  @. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b2*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a2*kdu
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@muladd function perform_step!(integrator,cache::Symplectic3Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,b1,b2,b3 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @. u = uprev + dt*b1*duprev
  # update velocity
  f.f2(integrator.t,u,duprev,kdu)
  @. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b2*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b3*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a3*kdu
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@muladd function perform_step!(integrator,cache::Symplectic4Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,b1,b2,b3,b4 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @. u = uprev + dt*b1*duprev
  # update velocity
  f.f2(t,u,duprev,kdu)
  @. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b2*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b3*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b4*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a4*kdu
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@muladd function perform_step!(integrator,cache::Symplectic45Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,a5,b1,b2,b3,b4,b5 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @. u = uprev + dt*b1*duprev
  # update velocity
  f.f2(t,u,duprev,kdu)
  @. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b2*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b3*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b4*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b5*ku

  f.f2(tnew,u,du,kdu)
  if typeof(integrator.alg) <: McAte42
    @. du = du + dt*a5*kdu
  end
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@muladd function perform_step!(integrator,cache::Symplectic5Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @. u = uprev + dt*b1*duprev
  # update velocity
  f.f2(integrator.t,u,duprev,kdu)
  @. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b2*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b3*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b4*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b5*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a5*kdu

  tnew = tnew + t+a5*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b6*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a6*kdu
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@muladd function perform_step!(integrator,cache::Symplectic6Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @. u = uprev + dt*b1*duprev
  # update velocity
  f.f2(integrator.t,u,duprev,kdu)
  @. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b2*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b3*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b4*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b5*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a5*kdu

  tnew = tnew + a5*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b6*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a6*kdu

  tnew = tnew + a6*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b7*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a7*kdu

  tnew = tnew + a7*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b8*ku

  f.f2(tnew,u,du,kdu)
  # @. du = du + dt*a8*kdu
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@muladd function perform_step!(integrator,cache::Symplectic62Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @. u = uprev + dt*b1*duprev
  # update velocity
  f.f2(integrator.t,u,duprev,kdu)
  @. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b2*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b3*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b4*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b5*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a5*kdu

  tnew = tnew + a5*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b6*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a6*kdu

  tnew = tnew + a6*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b7*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a7*kdu

  tnew = tnew + a7*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b8*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a8*kdu

  tnew = tnew + a8*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b9*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a9*kdu

  tnew = tnew + a9*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b10*ku

  f.f2(tnew,u,du,kdu)
  # @. du = du + dt*a10*kdu
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@muladd function perform_step!(integrator,cache::McAte8Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,
          b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @. u = uprev + dt*b1*duprev
  # update velocity
  f.f2(integrator.t,u,duprev,kdu)
  @. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b2*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b3*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b4*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b5*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a5*kdu

  tnew = tnew + a5*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b6*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a6*kdu

  tnew = tnew + a6*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b7*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a7*kdu

  tnew = tnew + a7*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b8*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a8*kdu

  tnew = tnew + a8*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b9*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a9*kdu

  tnew = tnew + a9*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b10*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a10*kdu

  tnew = tnew + a10*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b11*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a11*kdu

  tnew = tnew + a11*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b12*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a12*kdu

  tnew = tnew + a12*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b13*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a13*kdu

  tnew = tnew + a13*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b14*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a14*kdu

  tnew = tnew + a14*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b15*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a15*kdu

  tnew = tnew + a15*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b16*ku

  f.f2(tnew,u,du,kdu)
  # @. du = du + dt*a16*kdu
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@muladd function perform_step!(integrator,cache::KahanLi8Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,
          b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @. u = uprev + dt*b1*duprev
  # update velocity
  f.f2(integrator.t,u,duprev,kdu)
  @. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b2*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b3*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b4*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b5*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a5*kdu

  tnew = tnew + a5*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b6*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a6*kdu

  tnew = tnew + a6*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b7*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a7*kdu

  tnew = tnew + a7*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b8*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a8*kdu

  tnew = tnew + a8*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b9*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a9*kdu

  tnew = tnew + a9*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b10*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a10*kdu

  tnew = tnew + a10*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b11*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a11*kdu

  tnew = tnew + a11*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b12*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a12*kdu

  tnew = tnew + a12*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b13*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a13*kdu

  tnew = tnew + a13*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b14*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a14*kdu

  tnew = tnew + a14*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b15*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a15*kdu

  tnew = tnew + a15*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b16*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a16*kdu

  tnew = tnew + a16*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b17*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a17*kdu

  tnew = tnew + a17*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b18*ku

  f.f2(tnew,u,du,kdu)
  # @. du = du + dt*a18*kdu
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@muladd function perform_step!(integrator,cache::SofSpa10Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,
          a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,
          a35,a36,
          b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,
          b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,
          b35,b36 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @. u = uprev + dt*b1*duprev
  # update velocity
  f.f2(integrator.t,u,duprev,kdu)
  @. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b2*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b3*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b4*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b5*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a5*kdu

  tnew = tnew + a5*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b6*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a6*kdu

  tnew = tnew + a6*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b7*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a7*kdu

  tnew = tnew + a7*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b8*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a8*kdu

  tnew = tnew + a8*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b9*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a9*kdu

  tnew = tnew + a9*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b10*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a10*kdu

  tnew = tnew + a10*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b11*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a11*kdu

  tnew = tnew + a11*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b12*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a12*kdu

  tnew = tnew + a12*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b13*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a13*kdu

  tnew = tnew + a13*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b14*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a14*kdu

  tnew = tnew + a14*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b15*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a15*kdu

  tnew = tnew + a15*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b16*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a16*kdu

  tnew = tnew + a16*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b17*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a17*kdu

  tnew = tnew + a17*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b18*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a18*kdu

  tnew = tnew + a18*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b19*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a19*kdu

  tnew = tnew + a19*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b20*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a20*kdu

  tnew = tnew + a20*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b21*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a21*kdu

  tnew = tnew + a21*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b22*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a22*kdu

  tnew = tnew + a22*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b23*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a23*kdu

  tnew = tnew + a23*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b24*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a24*kdu

  tnew = tnew + a24*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b25*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a25*kdu

  tnew = tnew + a25*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b26*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a26*kdu

  tnew = tnew + a26*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b27*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a27*kdu

  tnew = tnew + a27*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b28*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a28*kdu

  tnew = tnew + a28*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b29*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a29*kdu

  tnew = tnew + a29*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b30*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a30*kdu

  tnew = tnew + a30*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b31*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a31*kdu

  tnew = tnew + a31*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b32*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a32*kdu

  tnew = tnew + a32*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b33*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a33*kdu

  tnew = tnew + a33*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b34*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a34*kdu

  tnew = tnew + a34*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b35*ku

  f.f2(tnew,u,du,kdu)
  @. du = du + dt*a35*kdu

  tnew = tnew + a35*dt
  f.f1(tnew,u,du,ku)
  @. u = u + dt*b36*ku

  f.f2(tnew,u,du,kdu)
  # @. du = du + dt*a30*kdu
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end
