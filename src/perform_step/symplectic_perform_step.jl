# http://www.chimica.unipd.it/antonino.polimeno/pubblica/downloads/JChemPhys_101_4062.pdf

function initialize!(integrator,cache::SymplecticEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  # Do the calculation pre
  # So that way FSAL interpolation
  duprev,uprev = integrator.uprev.x
  du,u = integrator.u.x
  kdu = integrator.f.f1(duprev,uprev,integrator.p,integrator.t)
  kuprev = integrator.f.f2(duprev,uprev,integrator.p,integrator.t)
  @muladd du = duprev + integrator.dt*kdu
  ku = integrator.f.f2(du,uprev,integrator.p,integrator.t)
  integrator.destats.nf2 += 1
  integrator.destats.nf += 2
  integrator.fsalfirst = ArrayPartition((kdu,kuprev))
  integrator.fsallast = ArrayPartition((zero(kdu),ku))
end

@muladd function perform_step!(integrator,cache::SymplecticEulerConstantCache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  duprev,uprev = integrator.uprev.x
  kuprev = integrator.fsalfirst.x[2]
  u = uprev + dt*kuprev
  # Now actually compute the step
  # Do it at the end for interpolations!
  kdu = f.f1(duprev,u,p,t)
  du = duprev + dt*kdu

  ku = f.f2(du,u,p,t)
  integrator.destats.nf2 += 1
  integrator.destats.nf += 1

  integrator.u = ArrayPartition((du,u))
  integrator.fsallast = ArrayPartition((kdu,ku))
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function initialize!(integrator,cache::SymplecticEulerCache)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  # Do the calculation pre
  # So that way FSAL interpolation
  duprev,uprev = integrator.uprev.x
  du,u = integrator.u.x
  kuprev = integrator.fsalfirst.x[2]
  kdu,ku = integrator.fsallast.x
  integrator.f.f1(kdu,duprev,uprev,integrator.p,integrator.t)
  integrator.f.f2(kuprev,duprev,uprev,integrator.p,integrator.t)
  @muladd @.. du = duprev + integrator.dt*kdu
  integrator.f.f2(ku,du,uprev,integrator.p,integrator.t)
  integrator.destats.nf += 1
  integrator.destats.nf2 += 2
end

@muladd function perform_step!(integrator,cache::SymplecticEulerCache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  duprev,uprev = integrator.uprev.x
  du,u = integrator.u.x
  kuprev = integrator.fsalfirst.x[2]
  kdu,ku = integrator.fsallast.x
  @.. u = uprev + dt*kuprev
  # Now actually compute the step
  # Do it at the end for interpolations!
  integrator.destats.nf2 += 1
  integrator.destats.nf += 1
  f.f1(kdu,duprev,u,p,t)
  @.. du = duprev + dt*kdu
  f.f2(ku,du,u,p,t)
end

const CachesInpHamilton = Union{Symplectic2Cache,Symplectic3Cache,
Symplectic4Cache,Symplectic45Cache,Symplectic5Cache,
Symplectic6Cache,Symplectic62Cache,
McAte8Cache,KahanLi8Cache,SofSpa10Cache,}
const CachesInpNewton = Union{VelocityVerletCache,}

const CachesNipHamilton = Union{Symplectic2ConstantCache,Symplectic3ConstantCache,
Symplectic4ConstantCache,Symplectic45ConstantCache,Symplectic5ConstantCache,
Symplectic6ConstantCache,Symplectic62ConstantCache,
McAte8ConstantCache,KahanLi8ConstantCache,SofSpa10ConstantCache,}
const CachesNipNewton = Union{VelocityVerletConstantCache,}

# some of the algorithms are designed only for the case
# f.f2(p, q, pa, t) = p which is the Newton/Lagrange equations
# If called with different functions (which are possible in the Hamiltonian case)
# an exception is thrown to avoid silently calculate wrong results.
verify_f2(::C, f, p, q, pa, t) where C<:CachesNipHamilton = f(p, q, pa, t)
verify_f2(::C, f, res, p, q, pa, t) where C<:CachesInpHamilton = f(res, p, q, pa, t)

function verify_f2(::C, f, p, q, pa, t) where C<:CachesNipNewton
    res = f(p, q, pa, t)
    res == p ? p : throwex()
end
function verify_f2(::C, f, res, p, q, pa, t) where C<:CachesInpNewton
    f(res, p, q, pa, t)
    res == p ? res : throwex()
end
throwex() = throw(ArgumentError("This algorithm is invalid if f2(p, q, t) != p"))

function initialize!(integrator,cache::C) where C<:Union{CachesInpHamilton,CachesInpNewton}
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k

  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast

  duprev,uprev = integrator.uprev.x
  integrator.f.f1(integrator.k[2].x[1],duprev,uprev,integrator.p,integrator.t)
  verify_f2(cache, integrator.f.f2, integrator.k[2].x[2], duprev, uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
  integrator.destats.nf2 += 1
end

function initialize!(integrator,cache::C) where C<:Union{CachesNipHamilton,CachesNipNewton}
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  duprev,uprev = integrator.uprev.x
  kdu  = integrator.f.f1(duprev,uprev,integrator.p,integrator.t)
  ku = verify_f2(cache, integrator.f.f2, duprev, uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
  integrator.destats.nf2 += 1
  integrator.fsalfirst = ArrayPartition((kdu,ku))
  integrator.k[2] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::VelocityVerletConstantCache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  duprev,uprev = integrator.uprev.x
  # x(t+Δt) = x(t) + v(t)*Δt + 1/2*a(t)*Δt^2
  ku = f.f1(duprev,uprev,p,t)
  dtsq = dt^2
  half = cache.half
  u = uprev + dt*duprev + dtsq*(half*ku)
  kdu = f.f1(duprev,u,p,t+dt)
  integrator.destats.nf += 2
  # v(t+Δt) = v(t) + 1/2*(a(t)+a(t+Δt))*Δt
  du = duprev + dt*(half*ku + half*kdu)

  integrator.u = ArrayPartition((du,u))
  integrator.fsallast = ArrayPartition((kdu,du))
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::VelocityVerletCache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  duprev,uprev = integrator.uprev.x
  du,u = integrator.u.x
  kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # x(t+Δt) = x(t) + v(t)*Δt + 1/2*a(t)*Δt^2
  f.f1(ku,duprev,uprev,p,t)
  dtsq = dt^2
  half = cache.half
  @.. u = uprev + dt*duprev + dtsq*(half*ku)
  f.f1(kdu,duprev,u,p,t+dt)
  integrator.destats.nf += 2
  # v(t+Δt) = v(t) + 1/2*(a(t)+a(t+Δt))*Δt
  @.. du = duprev + dt*(half*ku + half*kdu)
  copyto!(integrator.k[1].x[1],integrator.k[2].x[1])
  copyto!(integrator.k[1].x[2],integrator.k[2].x[2])
  copyto!(integrator.k[2].x[1],kdu)
  copyto!(integrator.k[2].x[2],du)
end

@muladd function perform_step!(integrator,cache::Symplectic2ConstantCache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,b1,b2 = cache
  duprev,uprev = integrator.uprev.x
  kuprev = integrator.k[2].x[2]
  # update position
  u = uprev + dt*b1*kuprev
  # update velocity
  kdu = f.f1(duprev,u,p,t)
  du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b2*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a2*kdu
  kdu = f.f1(du,u,p,tnew)
  ku = f.f2(du,u,p,tnew)
  integrator.destats.nf += 3
  integrator.destats.nf2 += 2

  integrator.u = ArrayPartition((du,u))
  integrator.fsallast = ArrayPartition((kdu,ku))
  integrator.k[1] = integrator.k[2]
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::Symplectic2Cache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,b1,b2 = cache.tab
  duprev,uprev = integrator.uprev.x
  kuprev = integrator.k[2].x[2]
  du,u = integrator.u.x
  kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @.. u = uprev + dt*b1*kuprev
  # update velocity
  f.f1(kdu,duprev,u,p,t)
  @.. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b2*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a2*kdu
  f.f1(kdu,du,u,p,tnew)
  f.f2(ku,du,u,p,tnew)
  integrator.destats.nf += 3
  integrator.destats.nf2 += 2
  copyto!(integrator.k[1].x[1],integrator.k[2].x[1])
  copyto!(integrator.k[1].x[2],integrator.k[2].x[2])
  copyto!(integrator.k[2].x[2],ku)
  copyto!(integrator.k[2].x[1],kdu)
end

@muladd function perform_step!(integrator,cache::Symplectic3ConstantCache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,b1,b2,b3 = cache
  duprev,uprev = integrator.uprev.x
  kuprev = integrator.k[2].x[2]
  # update position
  u = uprev + dt*b1*kuprev
  # update velocity
  kdu = f.f1(duprev,u,p,integrator.t)
  du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b2*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b3*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a3*kdu
  kdu = f.f1(du,u,p,tnew)
  ku = f.f2(du,u,p,tnew)
  integrator.destats.nf += 4
  integrator.destats.nf2 += 3
  integrator.u = ArrayPartition((du,u))
  integrator.fsallast = ArrayPartition((kdu,ku))
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::Symplectic3Cache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,b1,b2,b3 = cache.tab
  duprev,uprev = integrator.uprev.x
  kuprev = integrator.k[2].x[2]
  du,u = integrator.u.x
  kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @.. u = uprev + dt*b1*kuprev
  # update velocity
  f.f1(kdu,duprev,u,p,integrator.t)
  @.. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b2*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b3*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a3*kdu
  f.f1(kdu,du,u,p,tnew)
  f.f2(ku,du,u,p,tnew)
  integrator.destats.nf += 4
  integrator.destats.nf2 += 3
  copyto!(integrator.k[1].x[1],integrator.k[2].x[1])
  copyto!(integrator.k[1].x[2],integrator.k[2].x[2])
  copyto!(integrator.k[2].x[2],ku)
  copyto!(integrator.k[2].x[1],kdu)
end

@muladd function perform_step!(integrator,cache::Symplectic4ConstantCache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,b1,b2,b3,b4 = cache
  duprev,uprev = integrator.uprev.x
  kuprev = integrator.k[2].x[2]
  # update position
  u = uprev + dt*b1*kuprev
  # update velocity
  kdu = f.f1(duprev,u,p,t)
  du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b2*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b3*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b4*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a4*kdu
  kdu = f.f1(du,u,p,tnew)
  ku = f.f2(du,u,p,tnew)
  integrator.destats.nf += 5
  integrator.destats.nf2 += 4
  integrator.u = ArrayPartition((du,u))
  integrator.fsallast = ArrayPartition((kdu,ku))
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::Symplectic4Cache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,b1,b2,b3,b4 = cache.tab
  duprev,uprev = integrator.uprev.x
  du,u = integrator.u.x
  kuprev = integrator.k[2].x[2]
  kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @.. u = uprev + dt*b1*kuprev
  # update velocity
  f.f1(kdu,duprev,u,p,t)
  @.. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b2*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b3*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b4*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a4*kdu
  f.f1(kdu,du,u,p,tnew)
  f.f2(ku,du,u,p,tnew)
  integrator.destats.nf += 5
  integrator.destats.nf2 += 4
  copyto!(integrator.k[1].x[1],integrator.k[2].x[1])
  copyto!(integrator.k[1].x[2],integrator.k[2].x[2])
  copyto!(integrator.k[2].x[2],ku)
  copyto!(integrator.k[2].x[1],kdu)
end

@muladd function perform_step!(integrator,cache::Symplectic45ConstantCache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,a5,b1,b2,b3,b4,b5 = cache
  duprev,uprev = integrator.uprev.x
  kuprev = integrator.k[2].x[2]
  # update position
  u = uprev + dt*b1*kuprev
  # update velocity
  kdu = f.f1(duprev,u,p,t)
  du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b2*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b3*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b4*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b5*ku

  kdu = f.f1(du,u,p,tnew)
  if typeof(integrator.alg) <: McAte42
    du = du + dt*a5*kdu
    kdu = f.f1(du,u,p,tnew)
    integrator.destats.nf += 1
  end
  ku = f.f2(du,u,p,tnew)
  integrator.destats.nf += 5
  integrator.destats.nf2 += 5
  integrator.u = ArrayPartition((du,u))
  integrator.fsallast = ArrayPartition((kdu,ku))
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::Symplectic45Cache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,a5,b1,b2,b3,b4,b5 = cache.tab
  duprev,uprev = integrator.uprev.x
  du,u = integrator.u.x
  kuprev = integrator.k[2].x[2]
  kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @.. u = uprev + dt*b1*kuprev
  # update velocity
  f.f1(kdu,duprev,u,p,t)
  @.. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b2*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b3*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b4*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b5*ku

  f.f1(kdu,du,u,p,tnew)
  if typeof(integrator.alg) <: McAte42
    @.. du = du + dt*a5*kdu
    f.f1(kdu,du,u,p,tnew)
    integrator.destats.nf += 1
  end
  f.f2(ku,du,u,p,tnew)
  integrator.destats.nf += 5
  integrator.destats.nf2 += 5
  copyto!(integrator.k[1].x[1],integrator.k[2].x[1])
  copyto!(integrator.k[1].x[2],integrator.k[2].x[2])
  copyto!(integrator.k[2].x[2],ku)
  copyto!(integrator.k[2].x[1],kdu)
end

@muladd function perform_step!(integrator,cache::Symplectic5ConstantCache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6 = cache
  duprev,uprev = integrator.uprev.x
  kuprev = integrator.k[2].x[2]
  # update position
  u = uprev + dt*b1*kuprev
  # update velocity
  kdu = f.f1(duprev,u,p,integrator.t)
  du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b2*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b3*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b4*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b5*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a5*kdu

  tnew = tnew + t+a5*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b6*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a6*kdu
  kdu = f.f1(du,u,p,tnew)
  ku = f.f2(du,u,p,tnew)
  integrator.destats.nf += 7
  integrator.destats.nf2 += 6
  integrator.u = ArrayPartition((du,u))
  integrator.fsallast = ArrayPartition((kdu,ku))
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::Symplectic5Cache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6 = cache.tab
  duprev,uprev = integrator.uprev.x
  kuprev = integrator.k[2].x[2]
  du,u = integrator.u.x
  kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @.. u = uprev + dt*b1*kuprev
  # update velocity
  f.f1(kdu,duprev,u,p,integrator.t)
  @.. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b2*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b3*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b4*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b5*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a5*kdu

  tnew = tnew + t+a5*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b6*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a6*kdu
  f.f1(kdu,du,u,p,tnew)
  f.f2(ku,du,u,p,tnew)
  integrator.destats.nf += 7
  integrator.destats.nf2 += 6
  copyto!(integrator.k[1].x[1],integrator.k[2].x[1])
  copyto!(integrator.k[1].x[2],integrator.k[2].x[2])
  copyto!(integrator.k[2].x[2],ku)
  copyto!(integrator.k[2].x[1],kdu)
end

@muladd function perform_step!(integrator,cache::Symplectic6ConstantCache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8 = cache
  duprev,uprev = integrator.uprev.x
  kuprev = integrator.k[2].x[2]
  # update position
  u = uprev + dt*b1*kuprev
  # update velocity
  kdu = f.f1(duprev,u,p,integrator.t)
  du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b2*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b3*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b4*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b5*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a5*kdu

  tnew = tnew + a5*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b6*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a6*kdu

  tnew = tnew + a6*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b7*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a7*kdu

  tnew = tnew + a7*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b8*ku

  kdu = f.f1(du,u,p,tnew)
  # @.. du = du + dt*a8*kdu
  ku = f.f2(du,u,p,tnew)
  integrator.destats.nf += 8
  integrator.destats.nf2 += 8
  integrator.u = ArrayPartition((du,u))
  integrator.fsallast = ArrayPartition((kdu,ku))
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::Symplectic6Cache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8 = cache.tab
  duprev,uprev = integrator.uprev.x
  du,u = integrator.u.x
  kuprev = integrator.k[2].x[2]
  kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @.. u = uprev + dt*b1*kuprev
  # update velocity
  f.f1(kdu,duprev,u,p,integrator.t)
  @.. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b2*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b3*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b4*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b5*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a5*kdu

  tnew = tnew + a5*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b6*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a6*kdu

  tnew = tnew + a6*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b7*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a7*kdu

  tnew = tnew + a7*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b8*ku

  f.f1(kdu,du,u,p,tnew)
  # @.. du = du + dt*a8*kdu
  f.f2(ku,du,u,p,tnew)
  integrator.destats.nf += 8
  integrator.destats.nf2 += 8
  copyto!(integrator.k[1].x[1],integrator.k[2].x[1])
  copyto!(integrator.k[1].x[2],integrator.k[2].x[2])
  copyto!(integrator.k[2].x[2],ku)
  copyto!(integrator.k[2].x[1],kdu)
end

@muladd function perform_step!(integrator,cache::Symplectic62ConstantCache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10 = cache
  duprev,uprev = integrator.uprev.x
  kuprev = integrator.k[2].x[2]
  # update position
  u = uprev + dt*b1*kuprev
  # update velocity
  kdu = f.f1(duprev,u,p,integrator.t)
  du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b2*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b3*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b4*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b5*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a5*kdu

  tnew = tnew + a5*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b6*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a6*kdu

  tnew = tnew + a6*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b7*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a7*kdu

  tnew = tnew + a7*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b8*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a8*kdu

  tnew = tnew + a8*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b9*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a9*kdu

  tnew = tnew + a9*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b10*ku

  kdu = f.f1(du,u,p,tnew)
  # @.. du = du + dt*a10*kdu
  ku = f.f2(du,u,p,tnew)
  integrator.destats.nf += 10
  integrator.destats.nf2 += 10
  integrator.u = ArrayPartition((du,u))
  integrator.fsallast = ArrayPartition((kdu,ku))
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::Symplectic62Cache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10 = cache.tab
  duprev,uprev = integrator.uprev.x
  du,u = integrator.u.x
  kuprev = integrator.k[2].x[2]
  kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @.. u = uprev + dt*b1*kuprev
  # update velocity
  f.f1(kdu,duprev,u,p,integrator.t)
  @.. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b2*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b3*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b4*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b5*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a5*kdu

  tnew = tnew + a5*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b6*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a6*kdu

  tnew = tnew + a6*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b7*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a7*kdu

  tnew = tnew + a7*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b8*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a8*kdu

  tnew = tnew + a8*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b9*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a9*kdu

  tnew = tnew + a9*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b10*ku

  f.f1(kdu,du,u,p,tnew)
  # @.. du = du + dt*a10*kdu
  f.f2(ku,du,u,p,tnew)
  integrator.destats.nf += 10
  integrator.destats.nf2 += 10
  copyto!(integrator.k[1].x[1],integrator.k[2].x[1])
  copyto!(integrator.k[1].x[2],integrator.k[2].x[2])
  copyto!(integrator.k[2].x[2],ku)
  copyto!(integrator.k[2].x[1],kdu)
end

@muladd function perform_step!(integrator,cache::McAte8ConstantCache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,
          b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16 = cache
  duprev,uprev = integrator.uprev.x
  kuprev = integrator.k[2].x[2]
  # update position
  u = uprev + dt*b1*kuprev
  # update velocity
  kdu = f.f1(duprev,u,p,integrator.t)
  du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b2*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b3*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b4*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b5*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a5*kdu

  tnew = tnew + a5*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b6*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a6*kdu

  tnew = tnew + a6*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b7*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a7*kdu

  tnew = tnew + a7*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b8*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a8*kdu

  tnew = tnew + a8*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b9*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a9*kdu

  tnew = tnew + a9*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b10*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a10*kdu

  tnew = tnew + a10*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b11*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a11*kdu

  tnew = tnew + a11*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b12*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a12*kdu

  tnew = tnew + a12*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b13*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a13*kdu

  tnew = tnew + a13*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b14*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a14*kdu

  tnew = tnew + a14*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b15*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a15*kdu

  tnew = tnew + a15*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b16*ku

  kdu = f.f1(du,u,p,tnew)
  # @.. du = du + dt*a16*kdu
  ku = f.f2(du,u,p,tnew)
  integrator.destats.nf += 16
  integrator.destats.nf2 += 16
  integrator.u = ArrayPartition((du,u))
  integrator.fsallast = ArrayPartition((kdu,ku))
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::McAte8Cache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,
          b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16 = cache.tab
  duprev,uprev = integrator.uprev.x
  du,u = integrator.u.x
  kuprev = integrator.k[2].x[2]
  kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @.. u = uprev + dt*b1*kuprev
  # update velocity
  f.f1(kdu,duprev,u,p,integrator.t)
  @.. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b2*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b3*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b4*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b5*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a5*kdu

  tnew = tnew + a5*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b6*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a6*kdu

  tnew = tnew + a6*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b7*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a7*kdu

  tnew = tnew + a7*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b8*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a8*kdu

  tnew = tnew + a8*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b9*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a9*kdu

  tnew = tnew + a9*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b10*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a10*kdu

  tnew = tnew + a10*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b11*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a11*kdu

  tnew = tnew + a11*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b12*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a12*kdu

  tnew = tnew + a12*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b13*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a13*kdu

  tnew = tnew + a13*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b14*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a14*kdu

  tnew = tnew + a14*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b15*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a15*kdu

  tnew = tnew + a15*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b16*ku

  f.f1(kdu,du,u,p,tnew)
  # @.. du = du + dt*a16*kdu
  f.f2(ku,du,u,p,tnew)
  integrator.destats.nf += 16
  integrator.destats.nf2 += 16
  copyto!(integrator.k[1].x[1],integrator.k[2].x[1])
  copyto!(integrator.k[1].x[2],integrator.k[2].x[2])
  copyto!(integrator.k[2].x[2],ku)
  copyto!(integrator.k[2].x[1],kdu)
end

@muladd function perform_step!(integrator,cache::KahanLi8ConstantCache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,
          b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18 = cache
  duprev,uprev = integrator.uprev.x
  kuprev = integrator.k[2].x[2]
  # update position
  u = uprev + dt*b1*kuprev
  # update velocity
  kdu = f.f1(duprev,u,p,integrator.t)
  du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b2*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b3*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b4*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b5*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a5*kdu

  tnew = tnew + a5*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b6*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a6*kdu

  tnew = tnew + a6*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b7*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a7*kdu

  tnew = tnew + a7*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b8*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a8*kdu

  tnew = tnew + a8*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b9*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a9*kdu

  tnew = tnew + a9*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b10*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a10*kdu

  tnew = tnew + a10*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b11*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a11*kdu

  tnew = tnew + a11*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b12*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a12*kdu

  tnew = tnew + a12*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b13*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a13*kdu

  tnew = tnew + a13*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b14*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a14*kdu

  tnew = tnew + a14*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b15*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a15*kdu

  tnew = tnew + a15*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b16*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a16*kdu

  tnew = tnew + a16*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b17*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a17*kdu

  tnew = tnew + a17*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b18*ku

  kdu = f.f1(du,u,p,tnew)
  # @.. du = du + dt*a18*kdu
  ku = f.f2(du,u,p,tnew)

  integrator.destats.nf += 18
  integrator.destats.nf2 += 18
  integrator.u = ArrayPartition((du,u))
  integrator.fsallast = ArrayPartition((kdu,ku))
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::KahanLi8Cache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,
          b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18 = cache.tab
  duprev,uprev = integrator.uprev.x
  du,u = integrator.u.x
  kuprev = integrator.k[2].x[2]
  kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @.. u = uprev + dt*b1*kuprev
  # update velocity
  f.f1(kdu,duprev,u,p,integrator.t)
  @.. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b2*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b3*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b4*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b5*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a5*kdu

  tnew = tnew + a5*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b6*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a6*kdu

  tnew = tnew + a6*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b7*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a7*kdu

  tnew = tnew + a7*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b8*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a8*kdu

  tnew = tnew + a8*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b9*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a9*kdu

  tnew = tnew + a9*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b10*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a10*kdu

  tnew = tnew + a10*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b11*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a11*kdu

  tnew = tnew + a11*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b12*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a12*kdu

  tnew = tnew + a12*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b13*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a13*kdu

  tnew = tnew + a13*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b14*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a14*kdu

  tnew = tnew + a14*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b15*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a15*kdu

  tnew = tnew + a15*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b16*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a16*kdu

  tnew = tnew + a16*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b17*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a17*kdu

  tnew = tnew + a17*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b18*ku

  f.f1(kdu,du,u,p,tnew)
  # @.. du = du + dt*a18*kdu
  f.f2(ku,du,u,p,tnew)

  integrator.destats.nf += 18
  integrator.destats.nf2 += 18
  copyto!(integrator.k[1].x[1],integrator.k[2].x[1])
  copyto!(integrator.k[1].x[2],integrator.k[2].x[2])
  copyto!(integrator.k[2].x[2],ku)
  copyto!(integrator.k[2].x[1],kdu)
end

@muladd function perform_step!(integrator,cache::SofSpa10ConstantCache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,
          a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,
          a35,a36,
          b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,
          b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,
          b35,b36 = cache
  duprev,uprev = integrator.uprev.x
  kuprev = integrator.k[2].x[2]
  # update position
  u = uprev + dt*b1*kuprev
  # update velocity
  kdu = f.f1(duprev,u,p,integrator.t)
  du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b2*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b3*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b4*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b5*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a5*kdu

  tnew = tnew + a5*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b6*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a6*kdu

  tnew = tnew + a6*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b7*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a7*kdu

  tnew = tnew + a7*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b8*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a8*kdu

  tnew = tnew + a8*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b9*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a9*kdu

  tnew = tnew + a9*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b10*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a10*kdu

  tnew = tnew + a10*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b11*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a11*kdu

  tnew = tnew + a11*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b12*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a12*kdu

  tnew = tnew + a12*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b13*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a13*kdu

  tnew = tnew + a13*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b14*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a14*kdu

  tnew = tnew + a14*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b15*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a15*kdu

  tnew = tnew + a15*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b16*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a16*kdu

  tnew = tnew + a16*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b17*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a17*kdu

  tnew = tnew + a17*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b18*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a18*kdu

  tnew = tnew + a18*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b19*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a19*kdu

  tnew = tnew + a19*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b20*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a20*kdu

  tnew = tnew + a20*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b21*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a21*kdu

  tnew = tnew + a21*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b22*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a22*kdu

  tnew = tnew + a22*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b23*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a23*kdu

  tnew = tnew + a23*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b24*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a24*kdu

  tnew = tnew + a24*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b25*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a25*kdu

  tnew = tnew + a25*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b26*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a26*kdu

  tnew = tnew + a26*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b27*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a27*kdu

  tnew = tnew + a27*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b28*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a28*kdu

  tnew = tnew + a28*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b29*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a29*kdu

  tnew = tnew + a29*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b30*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a30*kdu

  tnew = tnew + a30*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b31*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a31*kdu

  tnew = tnew + a31*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b32*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a32*kdu

  tnew = tnew + a32*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b33*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a33*kdu

  tnew = tnew + a33*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b34*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a34*kdu

  tnew = tnew + a34*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b35*ku

  kdu = f.f1(du,u,p,tnew)
  du = du + dt*a35*kdu

  tnew = tnew + a35*dt
  ku = f.f2(du,u,p,tnew)
  u = u + dt*b36*ku

  kdu = f.f1(du,u,p,tnew)
  # @.. du = du + dt*a30*kdu
  ku = f.f2(du,u,p,tnew)
  integrator.destats.nf += 36
  integrator.destats.nf2 += 36
  integrator.u = ArrayPartition((du,u))
  integrator.fsallast = ArrayPartition((kdu,ku))
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::SofSpa10Cache,repeat_step=false)
  @unpack t,dt,f,p = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,
          a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,
          a35,a36,
          b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,
          b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,
          b35,b36 = cache.tab
  duprev,uprev = integrator.uprev.x
  du,u = integrator.u.x
  kuprev = integrator.k[2].x[2]
  kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @.. u = uprev + dt*b1*kuprev
  # update velocity
  f.f1(kdu,duprev,u,p,integrator.t)
  @.. du = duprev + dt*a1*kdu
  # update position & velocity
  tnew = t+a1*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b2*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a2*kdu

  # update position & velocity
  tnew = tnew + a2*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b3*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a3*kdu

  # update position & velocity
  tnew = tnew + a3*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b4*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a4*kdu

  # update position & velocity
  tnew = tnew + a4*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b5*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a5*kdu

  tnew = tnew + a5*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b6*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a6*kdu

  tnew = tnew + a6*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b7*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a7*kdu

  tnew = tnew + a7*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b8*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a8*kdu

  tnew = tnew + a8*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b9*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a9*kdu

  tnew = tnew + a9*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b10*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a10*kdu

  tnew = tnew + a10*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b11*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a11*kdu

  tnew = tnew + a11*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b12*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a12*kdu

  tnew = tnew + a12*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b13*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a13*kdu

  tnew = tnew + a13*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b14*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a14*kdu

  tnew = tnew + a14*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b15*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a15*kdu

  tnew = tnew + a15*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b16*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a16*kdu

  tnew = tnew + a16*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b17*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a17*kdu

  tnew = tnew + a17*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b18*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a18*kdu

  tnew = tnew + a18*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b19*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a19*kdu

  tnew = tnew + a19*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b20*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a20*kdu

  tnew = tnew + a20*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b21*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a21*kdu

  tnew = tnew + a21*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b22*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a22*kdu

  tnew = tnew + a22*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b23*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a23*kdu

  tnew = tnew + a23*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b24*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a24*kdu

  tnew = tnew + a24*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b25*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a25*kdu

  tnew = tnew + a25*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b26*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a26*kdu

  tnew = tnew + a26*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b27*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a27*kdu

  tnew = tnew + a27*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b28*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a28*kdu

  tnew = tnew + a28*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b29*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a29*kdu

  tnew = tnew + a29*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b30*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a30*kdu

  tnew = tnew + a30*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b31*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a31*kdu

  tnew = tnew + a31*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b32*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a32*kdu

  tnew = tnew + a32*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b33*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a33*kdu

  tnew = tnew + a33*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b34*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a34*kdu

  tnew = tnew + a34*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b35*ku

  f.f1(kdu,du,u,p,tnew)
  @.. du = du + dt*a35*kdu

  tnew = tnew + a35*dt
  f.f2(ku,du,u,p,tnew)
  @.. u = u + dt*b36*ku

  f.f1(kdu,du,u,p,tnew)
  # @.. du = du + dt*a30*kdu
  f.f2(ku,du,u,p,tnew)
  integrator.destats.nf += 36
  integrator.destats.nf2 += 36
  copyto!(integrator.k[1].x[1],integrator.k[2].x[1])
  copyto!(integrator.k[1].x[2],integrator.k[2].x[2])
  copyto!(integrator.k[2].x[2],ku)
  copyto!(integrator.k[2].x[1],kdu)
end
