@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::DiscreteConstantCache,idxs,T::Type{Val{0}})
  y₀
end

@inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::DiscreteCache,idxs,T::Type{Val{0}})
  recursivecopy!(out,y₀)
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 192
"""
@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::DP5ConstantCache,idxs,T::Type{Val{0}})
  Θ1 = 1-Θ
  y₀ + dt*Θ*(k[1]+Θ1*(k[2]+Θ*(k[3]+Θ1*k[4])))
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 192
"""
@inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::DP5Cache,idxs,T::Type{Val{0}})
  Θ1 = 1-Θ
  if out == nothing
    return y₀[idxs] + dt*(c1*k[1][idxs] + c2*k[2][idxs])
  else
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = y₀[i] + dt*Θ*(k[1][i]+Θ1*(k[2][i]+Θ*(k[3][i]+Θ1*k[4][i])))
    end
  end
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 192
"""
@inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::DP5ThreadedCache,idxs,T::Type{Val{0}})
  Θ1 = 1-Θ
  if out == nothing
    return y₀[idxs] + dt*(c1*k[1][idxs] + c2*k[2][idxs])
  else
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = y₀[i] + dt*Θ*(k[1][i]+Θ1*(k[2][i]+Θ*(k[3][i]+Θ1*k[4][i])))
    end
  end
end

"""
From MATLAB ODE Suite by Shampine
"""
@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23ConstantCache,Rosenbrock32ConstantCache},idxs,T::Type{Val{0}})
  d = cache.d
  c1 = Θ*(1-Θ)/(1-2d)
  c2 = Θ*(Θ-2d)/(1-2d)
  y₀ + dt*(c1*k[1] + c2*k[2])
end

# First Derivative of the dense output
@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23ConstantCache,Rosenbrock32ConstantCache},idxs,T::Type{Val{1}})
  d = cache.d
  c1diff = (1-2*Θ)/(1-2*d)
  c2diff = (2*Θ-2*d)/(1-2*d)
  c1diff*k[1] + c2diff*k[2]
end

"""
From MATLAB ODE Suite by Shampine
"""
@inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23Cache,Rosenbrock32Cache},idxs,T::Type{Val{0}})
  d = cache.tab.d
  c1 = Θ*(1-Θ)/(1-2d)
  c2 = Θ*(Θ-2d)/(1-2d)
  if out == nothing
    return y₀[idxs] + dt*(c1*k[1][idxs] + c2*k[2][idxs])
  else
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = y₀[i] + dt*(c1*k[1][i] + c2*k[2][i])
    end
  end
end

@inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Union{Rosenbrock23Cache,Rosenbrock32Cache},idxs,T::Type{Val{1}})
  d = cache.tab.d
  c1diff = (1-2*Θ)/(1-2*d)
  c2diff = (2*Θ-2*d)/(1-2*d)
  if out == nothing
    return c1diff*k[1][idxs] + c2diff*k[2][idxs]
  else
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = c1diff*k[1][i] + c2diff*k[2][i]
    end
  end
end

"""
Runge–Kutta pairs of order 5(4) satisfying only the first column
simplifying assumption

Ch. Tsitouras
"""
@inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Tsit5Cache,idxs,T::Type{Val{0}})
  b1Θ = -1.0530884977290216Θ * (Θ - 1.3299890189751412)*(Θ^2 - 1.4364028541716351Θ + 0.7139816917074209)
  b2Θ = 0.1017Θ^2 * (Θ^2 - 2.1966568338249754Θ + 1.2949852507374631)
  b3Θ = 2.490627285651252793Θ^2 * (Θ^2 - 2.38535645472061657Θ + 1.57803468208092486)
  b4Θ = -16.54810288924490272*(Θ - 1.21712927295533244)*(Θ - 0.61620406037800089)*Θ^2
  b5Θ = 47.37952196281928122*(Θ - 1.203071208372362603)*(Θ - 0.658047292653547382)*Θ^2
  b6Θ = -34.87065786149660974*(Θ - 1.2)*(Θ - 0.666666666666666667)*Θ^2
  b7Θ = 2.5*(Θ - 1)*(Θ - 0.6)*Θ^2
  if out == nothing
    return y₀[idxs] + dt*(k[1][idxs]*b1Θ + k[2][idxs]*b2Θ + k[3][idxs]*b3Θ + k[4][idxs]*b4Θ + k[5][idxs]*b5Θ + k[6][idxs]*b6Θ + k[7][idxs]*b7Θ)
  else
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = y₀[i] + dt*(k[1][i]*b1Θ + k[2][i]*b2Θ + k[3][i]*b3Θ + k[4][i]*b4Θ + k[5][i]*b5Θ + k[6][i]*b6Θ + k[7][i]*b7Θ)
    end
  end
end

"""
Runge–Kutta pairs of order 5(4) satisfying only the first column
simplifying assumption

Ch. Tsitouras
"""
@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::Tsit5ConstantCache,idxs,T::Type{Val{0}})
  b1Θ = -1.0530884977290216Θ * (Θ - 1.3299890189751412)*(Θ^2 - 1.4364028541716351Θ + 0.7139816917074209)
  b2Θ = 0.1017Θ^2 * (Θ^2 - 2.1966568338249754Θ + 1.2949852507374631)
  b3Θ = 2.490627285651252793Θ^2 * (Θ^2 - 2.38535645472061657Θ + 1.57803468208092486)
  b4Θ = -16.54810288924490272*(Θ - 1.21712927295533244)*(Θ - 0.61620406037800089)*Θ^2
  b5Θ = 47.37952196281928122*(Θ - 1.203071208372362603)*(Θ - 0.658047292653547382)*Θ^2
  b6Θ = -34.87065786149660974*(Θ - 1.2)*(Θ - 0.666666666666666667)*Θ^2
  b7Θ = 2.5*(Θ - 1)*(Θ - 0.6)*Θ^2
  y₀ + dt*(k[1]*b1Θ + k[2]*b2Θ + k[3]*b3Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ)
end

"""
Coefficients taken from RKSuite
"""
@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::BS5ConstantCache,idxs,T::Type{Val{0}})
  @unpack r016,r015,r014,r013,r012,r036,r035,r034,r033,r032,r046,r045,r044,r043,r042,r056,r055,r054,r053,r052,r066,r065,r064,r063,r062,r076,r075,r074,r073,r072,r086,r085,r084,r083,r082,r096,r095,r094,r093,r106,r105,r104,r103,r102,r116,r115,r114,r113,r112 = cache

  b1Θ  = @evalpoly(Θ, 0, 0, r012, r013, r014, r015, r016)
  b3Θ  = @evalpoly(Θ, 0, 0, r032, r033, r034, r035, r036)
  b4Θ  = @evalpoly(Θ, 0, 0, r042, r043, r044, r045, r046)
  b5Θ  = @evalpoly(Θ, 0, 0, r052, r053, r054, r055, r056)
  b6Θ  = @evalpoly(Θ, 0, 0, r062, r063, r064, r065, r066)
  b7Θ  = @evalpoly(Θ, 0, 0, r072, r073, r074, r075, r076)
  b8Θ  = @evalpoly(Θ, 0, 0, r082, r083, r084, r085, r086)
  b9Θ  = @evalpoly(Θ, 0, 0,    0, r093, r094, r095, r096)
  b10Θ = @evalpoly(Θ, 0, 0, r102, r103, r104, r105, r106)
  b11Θ = @evalpoly(Θ, 0, 0, r112, r113, r114, r115, r116)

  y₀ + dt*Θ*k[1] + dt*(k[1]*b1Θ  + k[3]*b3Θ + k[4]*b4Θ  + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ)
end

@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::BS5ConstantCache,idxs,T::Type{Val{1}})
  @unpack r016,r015,r014,r013,r012,r036,r035,r034,r033,r032,r046,r045,r044,r043,r042,r056,r055,r054,r053,r052,r066,r065,r064,r063,r062,r076,r075,r074,r073,r072,r086,r085,r084,r083,r082,r096,r095,r094,r093,r106,r105,r104,r103,r102,r116,r115,r114,r113,r112 = cache
  #TODO change names to <name>diff
  b1Θdiff  = @evalpoly(Θ, 0, 2*r012, 3*r013, 4*r014, 5*r015, 6*r016)
  b3Θdiff  = @evalpoly(Θ, 0, 2*r032, 3*r033, 4*r034, 5*r035, 6*r036)
  b4Θdiff  = @evalpoly(Θ, 0, 2*r042, 3*r043, 4*r044, 5*r045, 6*r046)
  b5Θdiff  = @evalpoly(Θ, 0, 2*r052, 3*r053, 4*r054, 5*r055, 6*r056)
  b6Θdiff  = @evalpoly(Θ, 0, 2*r062, 3*r063, 4*r064, 5*r065, 6*r066)
  b7Θdiff  = @evalpoly(Θ, 0, 2*r072, 3*r073, 4*r074, 5*r075, 6*r076)
  b8Θdiff  = @evalpoly(Θ, 0, 2*r082, 3*r083, 4*r084, 5*r085, 6*r086)
  b9Θdiff  = @evalpoly(Θ, 0,      0, 3*r093, 4*r094, 5*r095, 6*r096)
  b10Θdiff = @evalpoly(Θ, 0, 2*r102, 3*r103, 4*r104, 5*r105, 6*r106)
  b11Θdiff = @evalpoly(Θ, 0, 2*r112, 3*r113, 4*r114, 5*r115, 6*r116)

  k[1] + k[1]*b1Θdiff  + k[3]*b3Θdiff + k[4]*b4Θdiff  + k[5]*b5Θdiff + k[6]*b6Θdiff + k[7]*b7Θdiff + k[8]*b8Θdiff + k[9]*b9Θdiff + k[10]*b10Θdiff + k[11]*b11Θdiff
end

"""
Coefficients taken from RKSuite
"""
@inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::BS5Cache,idxs,T::Type{Val{0}})
  @unpack r016,r015,r014,r013,r012,r036,r035,r034,r033,r032,r046,r045,r044,r043,r042,r056,r055,r054,r053,r052,r066,r065,r064,r063,r062,r076,r075,r074,r073,r072,r086,r085,r084,r083,r082,r096,r095,r094,r093,r106,r105,r104,r103,r102,r116,r115,r114,r113,r112 = cache.tab

  b1Θ  = @evalpoly(Θ, 0, 0, r012, r013, r014, r015, r016)
  b3Θ  = @evalpoly(Θ, 0, 0, r032, r033, r034, r035, r036)
  b4Θ  = @evalpoly(Θ, 0, 0, r042, r043, r044, r045, r046)
  b5Θ  = @evalpoly(Θ, 0, 0, r052, r053, r054, r055, r056)
  b6Θ  = @evalpoly(Θ, 0, 0, r062, r063, r064, r065, r066)
  b7Θ  = @evalpoly(Θ, 0, 0, r072, r073, r074, r075, r076)
  b8Θ  = @evalpoly(Θ, 0, 0, r082, r083, r084, r085, r086)
  b9Θ  = @evalpoly(Θ, 0, 0,    0, r093, r094, r095, r096)
  b10Θ = @evalpoly(Θ, 0, 0, r102, r103, r104, r105, r106)
  b11Θ = @evalpoly(Θ, 0, 0, r112, r113, r114, r115, r116)

  if out == nothing
    return y₀[idxs] + dt*Θ*k[1][idxs] + dt*(k[1][idxs]*b1Θ  + k[3][idxs]*b3Θ + k[4][idxs]*b4Θ  + k[5][idxs]*b5Θ + k[6][idxs]*b6Θ + k[7][idxs]*b7Θ + k[8][idxs]*b8Θ + k[9][idxs]*b9Θ + k[10][idxs]*b10Θ + k[11][idxs]*b11Θ)
  else
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = y₀[i] + dt*Θ*k[1][i] + dt*(k[1][i]*b1Θ  + k[3][i]*b3Θ + k[4][i]*b4Θ  + k[5][i]*b5Θ + k[6][i]*b6Θ + k[7][i]*b7Θ + k[8][i]*b8Θ + k[9][i]*b9Θ + k[10][i]*b10Θ + k[11][i]*b11Θ)
    end
  end
end

"""

"""
@inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Vern6Cache,idxs,T::Type{Val{0}})
  @unpack r011,r012,r013,r014,r015,r016,r042,r043,r044,r045,r046,r052,r053,r054,r055,r056,r062,r063,r064,r065,r066,r072,r073,r074,r075,r076,r082,r083,r084,r085,r086,r092,r093,r094,r095,r096,r102,r103,r104,r105,r106,r112,r113,r114,r115,r116,r122,r123,r124,r125,r126 = cache.tab

  b1Θ  = @evalpoly(Θ, 0, r011, r012, r013, r014, r015, r016)
  b4Θ  = @evalpoly(Θ, 0,    0, r042, r043, r044, r045, r046)
  b5Θ  = @evalpoly(Θ, 0,    0, r052, r053, r054, r055, r056)
  b6Θ  = @evalpoly(Θ, 0,    0, r062, r063, r064, r065, r066)
  b7Θ  = @evalpoly(Θ, 0,    0, r072, r073, r074, r075, r076)
  b8Θ  = @evalpoly(Θ, 0,    0, r082, r083, r084, r085, r086)
  b9Θ  = @evalpoly(Θ, 0,    0, r092, r093, r094, r095, r096)
  b10Θ = @evalpoly(Θ, 0,    0, r102, r103, r104, r105, r106)
  b11Θ = @evalpoly(Θ, 0,    0, r112, r113, r114, r115, r116)
  b12Θ = @evalpoly(Θ, 0,    0, r122, r123, r124, r125, r126)

  if out == nothing
    return y₀[idxs] + dt*(k[1][idxs]*b1Θ + k[4][idxs]*b4Θ + k[5][idxs]*b5Θ + k[6][idxs]*b6Θ + k[7][idxs]*b7Θ + k[8][idxs]*b8Θ + k[9][idxs]*b9Θ + k[10][idxs]*b10Θ + k[11][idxs]*b11Θ + k[12][idxs]*b12Θ)
  else
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = y₀[i] + dt*(k[1][i]*b1Θ + k[4][i]*b4Θ + k[5][i]*b5Θ + k[6][i]*b6Θ + k[7][i]*b7Θ + k[8][i]*b8Θ + k[9][i]*b9Θ + k[10][i]*b10Θ + k[11][i]*b11Θ + k[12][i]*b12Θ)
    end
  end
end

"""

"""
@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::Vern6ConstantCache,idxs,T::Type{Val{0}})
  @unpack r011,r012,r013,r014,r015,r016,r042,r043,r044,r045,r046,r052,r053,r054,r055,r056,r062,r063,r064,r065,r066,r072,r073,r074,r075,r076,r082,r083,r084,r085,r086,r092,r093,r094,r095,r096,r102,r103,r104,r105,r106,r112,r113,r114,r115,r116,r122,r123,r124,r125,r126 = cache

  b1Θ  = @evalpoly(Θ, 0, r011, r012, r013, r014, r015, r016)
  b4Θ  = @evalpoly(Θ, 0,    0, r042, r043, r044, r045, r046)
  b5Θ  = @evalpoly(Θ, 0,    0, r052, r053, r054, r055, r056)
  b6Θ  = @evalpoly(Θ, 0,    0, r062, r063, r064, r065, r066)
  b7Θ  = @evalpoly(Θ, 0,    0, r072, r073, r074, r075, r076)
  b8Θ  = @evalpoly(Θ, 0,    0, r082, r083, r084, r085, r086)
  b9Θ  = @evalpoly(Θ, 0,    0, r092, r093, r094, r095, r096)
  b10Θ = @evalpoly(Θ, 0,    0, r102, r103, r104, r105, r106)
  b11Θ = @evalpoly(Θ, 0,    0, r112, r113, r114, r115, r116)
  b12Θ = @evalpoly(Θ, 0,    0, r122, r123, r124, r125, r126)

  y₀ + dt*(k[1]*b1Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ + k[12]*b12Θ)
end

"""

"""
@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::Vern7ConstantCache,idxs,T::Type{Val{0}})
  @unpack r011,r012,r013,r014,r015,r016,r017,r042,r043,r044,r045,r046,r047,r052,r053,r054,r055,r056,r057,r062,r063,r064,r065,r066,r067,r072,r073,r074,r075,r076,r077,r082,r083,r084,r085,r086,r087,r092,r093,r094,r095,r096,r097,r112,r113,r114,r115,r116,r117,r122,r123,r124,r125,r126,r127,r132,r133,r134,r135,r136,r137,r142,r143,r144,r145,r146,r147,r152,r153,r154,r155,r156,r157,r162,r163,r164,r165,r166,r167 = cache

  b1Θ  = @evalpoly(Θ, 0, r011, r012, r013, r014, r015, r016, r017)
  b4Θ  = @evalpoly(Θ, 0,    0, r042, r043, r044, r045, r046, r047)
  b5Θ  = @evalpoly(Θ, 0,    0, r052, r053, r054, r055, r056, r057)
  b6Θ  = @evalpoly(Θ, 0,    0, r062, r063, r064, r065, r066, r067)
  b7Θ  = @evalpoly(Θ, 0,    0, r072, r073, r074, r075, r076, r077)
  b8Θ  = @evalpoly(Θ, 0,    0, r082, r083, r084, r085, r086, r087)
  b9Θ  = @evalpoly(Θ, 0,    0, r092, r093, r094, r095, r096, r097)
  b11Θ = @evalpoly(Θ, 0,    0, r112, r113, r114, r115, r116, r117)
  b12Θ = @evalpoly(Θ, 0,    0, r122, r123, r124, r125, r126, r127)
  b13Θ = @evalpoly(Θ, 0,    0, r132, r133, r134, r135, r136, r137)
  b14Θ = @evalpoly(Θ, 0,    0, r142, r143, r144, r145, r146, r147)
  b15Θ = @evalpoly(Θ, 0,    0, r152, r153, r154, r155, r156, r157)
  b16Θ = @evalpoly(Θ, 0,    0, r162, r163, r164, r165, r166, r167)

  y₀ + dt*(k[1]*b1Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[11]*b11Θ + k[12]*b12Θ + k[13]*b13Θ + k[14]*b14Θ + k[15]*b15Θ + k[16]*b16Θ)
end

"""

"""
@inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Vern7Cache,idxs,T::Type{Val{0}})
  @unpack r011,r012,r013,r014,r015,r016,r017,r042,r043,r044,r045,r046,r047,r052,r053,r054,r055,r056,r057,r062,r063,r064,r065,r066,r067,r072,r073,r074,r075,r076,r077,r082,r083,r084,r085,r086,r087,r092,r093,r094,r095,r096,r097,r112,r113,r114,r115,r116,r117,r122,r123,r124,r125,r126,r127,r132,r133,r134,r135,r136,r137,r142,r143,r144,r145,r146,r147,r152,r153,r154,r155,r156,r157,r162,r163,r164,r165,r166,r167 = cache.tab

  b1Θ  = @evalpoly(Θ, 0, r011, r012, r013, r014, r015, r016, r017)
  b4Θ  = @evalpoly(Θ, 0,    0, r042, r043, r044, r045, r046, r047)
  b5Θ  = @evalpoly(Θ, 0,    0, r052, r053, r054, r055, r056, r057)
  b6Θ  = @evalpoly(Θ, 0,    0, r062, r063, r064, r065, r066, r067)
  b7Θ  = @evalpoly(Θ, 0,    0, r072, r073, r074, r075, r076, r077)
  b8Θ  = @evalpoly(Θ, 0,    0, r082, r083, r084, r085, r086, r087)
  b9Θ  = @evalpoly(Θ, 0,    0, r092, r093, r094, r095, r096, r097)
  b11Θ = @evalpoly(Θ, 0,    0, r112, r113, r114, r115, r116, r117)
  b12Θ = @evalpoly(Θ, 0,    0, r122, r123, r124, r125, r126, r127)
  b13Θ = @evalpoly(Θ, 0,    0, r132, r133, r134, r135, r136, r137)
  b14Θ = @evalpoly(Θ, 0,    0, r142, r143, r144, r145, r146, r147)
  b15Θ = @evalpoly(Θ, 0,    0, r152, r153, r154, r155, r156, r157)
  b16Θ = @evalpoly(Θ, 0,    0, r162, r163, r164, r165, r166, r167)

  if out == nothing
    return y₀[idxs] + dt*(k[1][idxs]*b1Θ + k[4][idxs]*b4Θ + k[5][idxs]*b5Θ + k[6][idxs]*b6Θ + k[7][idxs]*b7Θ + k[8][idxs]*b8Θ + k[9][idxs]*b9Θ + k[11][idxs]*b11Θ + k[12][idxs]*b12Θ + k[13][idxs]*b13Θ + k[14][idxs]*b14Θ + k[15][idxs]*b15Θ + k[16][idxs]*b16Θ)
  else
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = y₀[i] + dt*(k[1][i]*b1Θ + k[4][i]*b4Θ + k[5][i]*b5Θ + k[6][i]*b6Θ + k[7][i]*b7Θ + k[8][i]*b8Θ + k[9][i]*b9Θ + k[11][i]*b11Θ + k[12][i]*b12Θ + k[13][i]*b13Θ + k[14][i]*b14Θ + k[15][i]*b15Θ + k[16][i]*b16Θ)
    end
  end
end

"""

"""
@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::Vern8ConstantCache,idxs,T::Type{Val{0}})
  @unpack r011,r012,r013,r014,r015,r016,r017,r018,r062,r063,r064,r065,r066,r067,r068,r072,r073,r074,r075,r076,r077,r078,r082,r083,r084,r085,r086,r087,r088,r092,r093,r094,r095,r096,r097,r098,r102,r103,r104,r105,r106,r107,r108,r112,r113,r114,r115,r116,r117,r118,r122,r123,r124,r125,r126,r127,r128,r142,r143,r144,r145,r146,r147,r148,r152,r153,r154,r155,r156,r157,r158,r162,r163,r164,r165,r166,r167,r168,r172,r173,r174,r175,r176,r177,r178,r182,r183,r184,r185,r186,r187,r188,r192,r193,r194,r195,r196,r197,r198,r202,r203,r204,r205,r206,r207,r208,r212,r213,r214,r215,r216,r217,r218 = cache

  b1Θ  = @evalpoly(Θ, 0, r011, r012, r013, r014, r015, r016, r017, r018)
  b6Θ  = @evalpoly(Θ, 0,    0, r062, r063, r064, r065, r066, r067, r068)
  b7Θ  = @evalpoly(Θ, 0,    0, r072, r073, r074, r075, r076, r077, r078)
  b8Θ  = @evalpoly(Θ, 0,    0, r082, r083, r084, r085, r086, r087, r088)
  b9Θ  = @evalpoly(Θ, 0,    0, r092, r093, r094, r095, r096, r097, r098)
  b10Θ = @evalpoly(Θ, 0,    0, r102, r103, r104, r105, r106, r107, r108)
  b11Θ = @evalpoly(Θ, 0,    0, r112, r113, r114, r115, r116, r117, r118)
  b12Θ = @evalpoly(Θ, 0,    0, r122, r123, r124, r125, r126, r127, r128)
  b14Θ = @evalpoly(Θ, 0,    0, r142, r143, r144, r145, r146, r147, r148)
  b15Θ = @evalpoly(Θ, 0,    0, r152, r153, r154, r155, r156, r157, r158)
  b16Θ = @evalpoly(Θ, 0,    0, r162, r163, r164, r165, r166, r167, r168)
  b17Θ = @evalpoly(Θ, 0,    0, r172, r173, r174, r175, r176, r177, r178)
  b18Θ = @evalpoly(Θ, 0,    0, r182, r183, r184, r185, r186, r187, r188)
  b19Θ = @evalpoly(Θ, 0,    0, r192, r193, r194, r195, r196, r197, r198)
  b20Θ = @evalpoly(Θ, 0,    0, r202, r203, r204, r205, r206, r207, r208)
  b21Θ = @evalpoly(Θ, 0,    0, r212, r213, r214, r215, r216, r217, r218)

  y₀ + dt*(k[1]*b1Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ + k[12]*b12Θ + k[14]*b14Θ + k[15]*b15Θ + k[16]*b16Θ + k[17]*b17Θ + k[18]*b18Θ + k[19]*b19Θ + k[20]*b20Θ + k[21]*b21Θ)
end

"""

"""
@inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Vern8Cache,idxs,T::Type{Val{0}})
  @unpack r011,r012,r013,r014,r015,r016,r017,r018,r062,r063,r064,r065,r066,r067,r068,r072,r073,r074,r075,r076,r077,r078,r082,r083,r084,r085,r086,r087,r088,r092,r093,r094,r095,r096,r097,r098,r102,r103,r104,r105,r106,r107,r108,r112,r113,r114,r115,r116,r117,r118,r122,r123,r124,r125,r126,r127,r128,r142,r143,r144,r145,r146,r147,r148,r152,r153,r154,r155,r156,r157,r158,r162,r163,r164,r165,r166,r167,r168,r172,r173,r174,r175,r176,r177,r178,r182,r183,r184,r185,r186,r187,r188,r192,r193,r194,r195,r196,r197,r198,r202,r203,r204,r205,r206,r207,r208,r212,r213,r214,r215,r216,r217,r218 = cache.tab

  b1Θ  = @evalpoly(Θ, 0, r011, r012, r013, r014, r015, r016, r017, r018)
  b6Θ  = @evalpoly(Θ, 0,    0, r062, r063, r064, r065, r066, r067, r068)
  b7Θ  = @evalpoly(Θ, 0,    0, r072, r073, r074, r075, r076, r077, r078)
  b8Θ  = @evalpoly(Θ, 0,    0, r082, r083, r084, r085, r086, r087, r088)
  b9Θ  = @evalpoly(Θ, 0,    0, r092, r093, r094, r095, r096, r097, r098)
  b10Θ = @evalpoly(Θ, 0,    0, r102, r103, r104, r105, r106, r107, r108)
  b11Θ = @evalpoly(Θ, 0,    0, r112, r113, r114, r115, r116, r117, r118)
  b12Θ = @evalpoly(Θ, 0,    0, r122, r123, r124, r125, r126, r127, r128)
  b14Θ = @evalpoly(Θ, 0,    0, r142, r143, r144, r145, r146, r147, r148)
  b15Θ = @evalpoly(Θ, 0,    0, r152, r153, r154, r155, r156, r157, r158)
  b16Θ = @evalpoly(Θ, 0,    0, r162, r163, r164, r165, r166, r167, r168)
  b17Θ = @evalpoly(Θ, 0,    0, r172, r173, r174, r175, r176, r177, r178)
  b18Θ = @evalpoly(Θ, 0,    0, r182, r183, r184, r185, r186, r187, r188)
  b19Θ = @evalpoly(Θ, 0,    0, r192, r193, r194, r195, r196, r197, r198)
  b20Θ = @evalpoly(Θ, 0,    0, r202, r203, r204, r205, r206, r207, r208)
  b21Θ = @evalpoly(Θ, 0,    0, r212, r213, r214, r215, r216, r217, r218)

  if out == nothing
    return y₀[idxs] + dt*(k[1][idxs]*b1Θ + k[6][idxs]*b6Θ + k[7][idxs]*b7Θ + k[8][idxs]*b8Θ + k[9][idxs]*b9Θ + k[10][idxs]*b10Θ + k[11][idxs]*b11Θ + k[12][idxs]*b12Θ + k[14][idxs]*b14Θ + k[15][idxs]*b15Θ + k[16][idxs]*b16Θ + k[17][idxs]*b17Θ + k[18][idxs]*b18Θ + k[19][idxs]*b19Θ + k[20][idxs]*b20Θ + k[21][idxs]*b21Θ)
  else
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = y₀[i] + dt*(k[1][i]*b1Θ + k[6][i]*b6Θ + k[7][i]*b7Θ + k[8][i]*b8Θ + k[9][i]*b9Θ + k[10][i]*b10Θ + k[11][i]*b11Θ + k[12][i]*b12Θ + k[14][i]*b14Θ + k[15][i]*b15Θ + k[16][i]*b16Θ + k[17][i]*b17Θ + k[18][i]*b18Θ + k[19][i]*b19Θ + k[20][i]*b20Θ + k[21][i]*b21Θ)
    end
  end
end

"""

"""
@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::Vern9ConstantCache,idxs,T::Type{Val{0}})
  @unpack r011,r012,r013,r014,r015,r016,r017,r018,r019,r082,r083,r084,r085,r086,r087,r088,r089,r092,r093,r094,r095,r096,r097,r098,r099,r102,r103,r104,r105,r106,r107,r108,r109,r112,r113,r114,r115,r116,r117,r118,r119,r122,r123,r124,r125,r126,r127,r128,r129,r132,r133,r134,r135,r136,r137,r138,r139,r142,r143,r144,r145,r146,r147,r148,r149,r152,r153,r154,r155,r156,r157,r158,r159,r172,r173,r174,r175,r176,r177,r178,r179,r182,r183,r184,r185,r186,r187,r188,r189,r192,r193,r194,r195,r196,r197,r198,r199,r202,r203,r204,r205,r206,r207,r208,r209,r212,r213,r214,r215,r216,r217,r218,r219,r222,r223,r224,r225,r226,r227,r228,r229,r232,r233,r234,r235,r236,r237,r238,r239,r242,r243,r244,r245,r246,r247,r248,r249,r252,r253,r254,r255,r256,r257,r258,r259,r262,r263,r264,r265,r266,r267,r268,r269 = cache

  b1Θ  = @evalpoly(Θ, 0, r011, r012, r013, r014, r015, r016, r017, r018, r019)
  b8Θ  = @evalpoly(Θ, 0,    0, r082, r083, r084, r085, r086, r087, r088, r089)
  b9Θ  = @evalpoly(Θ, 0,    0, r092, r093, r094, r095, r096, r097, r098, r099)
  b10Θ = @evalpoly(Θ, 0,    0, r102, r103, r104, r105, r106, r107, r108, r109)
  b11Θ = @evalpoly(Θ, 0,    0, r112, r113, r114, r115, r116, r117, r118, r119)
  b12Θ = @evalpoly(Θ, 0,    0, r122, r123, r124, r125, r126, r127, r128, r129)
  b13Θ = @evalpoly(Θ, 0,    0, r132, r133, r134, r135, r136, r137, r138, r139)
  b14Θ = @evalpoly(Θ, 0,    0, r142, r143, r144, r145, r146, r147, r148, r149)
  b15Θ = @evalpoly(Θ, 0,    0, r152, r153, r154, r155, r156, r157, r158, r159)
  b17Θ = @evalpoly(Θ, 0,    0, r172, r173, r174, r175, r176, r177, r178, r179)
  b18Θ = @evalpoly(Θ, 0,    0, r182, r183, r184, r185, r186, r187, r188, r189)
  b19Θ = @evalpoly(Θ, 0,    0, r192, r193, r194, r195, r196, r197, r198, r199)
  b20Θ = @evalpoly(Θ, 0,    0, r202, r203, r204, r205, r206, r207, r208, r209)
  b21Θ = @evalpoly(Θ, 0,    0, r212, r213, r214, r215, r216, r217, r218, r219)
  b22Θ = @evalpoly(Θ, 0,    0, r222, r223, r224, r225, r226, r227, r228, r229)
  b23Θ = @evalpoly(Θ, 0,    0, r232, r233, r234, r235, r236, r237, r238, r239)
  b24Θ = @evalpoly(Θ, 0,    0, r242, r243, r244, r245, r246, r247, r248, r249)
  b25Θ = @evalpoly(Θ, 0,    0, r252, r253, r254, r255, r256, r257, r258, r259)
  b26Θ = @evalpoly(Θ, 0,    0, r262, r263, r264, r265, r266, r267, r268, r269)

  y₀ + dt*(k[1]*b1Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ + k[12]*b12Θ + k[13]*b13Θ + k[14]*b14Θ + k[15]*b15Θ + k[17]*b17Θ + k[18]*b18Θ + k[19]*b19Θ + k[20]*b20Θ + k[21]*b21Θ + k[22]*b22Θ + k[23]*b23Θ + k[24]*b24Θ + k[25]*b25Θ + k[26]*b26Θ)
end

"""

"""
@inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::Vern9Cache,idxs,T::Type{Val{0}})
  @unpack r011,r012,r013,r014,r015,r016,r017,r018,r019,r082,r083,r084,r085,r086,r087,r088,r089,r092,r093,r094,r095,r096,r097,r098,r099,r102,r103,r104,r105,r106,r107,r108,r109,r112,r113,r114,r115,r116,r117,r118,r119,r122,r123,r124,r125,r126,r127,r128,r129,r132,r133,r134,r135,r136,r137,r138,r139,r142,r143,r144,r145,r146,r147,r148,r149,r152,r153,r154,r155,r156,r157,r158,r159,r172,r173,r174,r175,r176,r177,r178,r179,r182,r183,r184,r185,r186,r187,r188,r189,r192,r193,r194,r195,r196,r197,r198,r199,r202,r203,r204,r205,r206,r207,r208,r209,r212,r213,r214,r215,r216,r217,r218,r219,r222,r223,r224,r225,r226,r227,r228,r229,r232,r233,r234,r235,r236,r237,r238,r239,r242,r243,r244,r245,r246,r247,r248,r249,r252,r253,r254,r255,r256,r257,r258,r259,r262,r263,r264,r265,r266,r267,r268,r269 = cache.tab

  b1Θ  = @evalpoly(Θ, 0, r011, r012, r013, r014, r015, r016, r017, r018, r019)
  b8Θ  = @evalpoly(Θ, 0,    0, r082, r083, r084, r085, r086, r087, r088, r089)
  b9Θ  = @evalpoly(Θ, 0,    0, r092, r093, r094, r095, r096, r097, r098, r099)
  b10Θ = @evalpoly(Θ, 0,    0, r102, r103, r104, r105, r106, r107, r108, r109)
  b11Θ = @evalpoly(Θ, 0,    0, r112, r113, r114, r115, r116, r117, r118, r119)
  b12Θ = @evalpoly(Θ, 0,    0, r122, r123, r124, r125, r126, r127, r128, r129)
  b13Θ = @evalpoly(Θ, 0,    0, r132, r133, r134, r135, r136, r137, r138, r139)
  b14Θ = @evalpoly(Θ, 0,    0, r142, r143, r144, r145, r146, r147, r148, r149)
  b15Θ = @evalpoly(Θ, 0,    0, r152, r153, r154, r155, r156, r157, r158, r159)
  b17Θ = @evalpoly(Θ, 0,    0, r172, r173, r174, r175, r176, r177, r178, r179)
  b18Θ = @evalpoly(Θ, 0,    0, r182, r183, r184, r185, r186, r187, r188, r189)
  b19Θ = @evalpoly(Θ, 0,    0, r192, r193, r194, r195, r196, r197, r198, r199)
  b20Θ = @evalpoly(Θ, 0,    0, r202, r203, r204, r205, r206, r207, r208, r209)
  b21Θ = @evalpoly(Θ, 0,    0, r212, r213, r214, r215, r216, r217, r218, r219)
  b22Θ = @evalpoly(Θ, 0,    0, r222, r223, r224, r225, r226, r227, r228, r229)
  b23Θ = @evalpoly(Θ, 0,    0, r232, r233, r234, r235, r236, r237, r238, r239)
  b24Θ = @evalpoly(Θ, 0,    0, r242, r243, r244, r245, r246, r247, r248, r249)
  b25Θ = @evalpoly(Θ, 0,    0, r252, r253, r254, r255, r256, r257, r258, r259)
  b26Θ = @evalpoly(Θ, 0,    0, r262, r263, r264, r265, r266, r267, r268, r269)

  if out == nothing
    return y₀[idxs] + dt*(k[1][idxs]*b1Θ + k[8][idxs]*b8Θ + k[9][idxs]*b9Θ + k[10][idxs]*b10Θ + k[11][idxs]*b11Θ + k[12][idxs]*b12Θ + k[13][idxs]*b13Θ + k[14][idxs]*b14Θ + k[15][idxs]*b15Θ + k[17][idxs]*b17Θ + k[18][idxs]*b18Θ + k[19][idxs]*b19Θ + k[20][idxs]*b20Θ + k[21][idxs]*b21Θ + k[22][idxs]*b22Θ + k[23][idxs]*b23Θ + k[24][idxs]*b24Θ + k[25][idxs]*b25Θ + k[26][idxs]*b26Θ)
  else
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = y₀[i] + dt*(k[1][i]*b1Θ + k[8][i]*b8Θ + k[9][i]*b9Θ + k[10][i]*b10Θ + k[11][i]*b11Θ + k[12][i]*b12Θ + k[13][i]*b13Θ + k[14][i]*b14Θ + k[15][i]*b15Θ + k[17][i]*b17Θ + k[18][i]*b18Θ + k[19][i]*b19Θ + k[20][i]*b20Θ + k[21][i]*b21Θ + k[22][i]*b22Θ + k[23][i]*b23Θ + k[24][i]*b24Θ + k[25][i]*b25Θ + k[26][i]*b26Θ)
    end
  end
end

"""

"""
@inline function ode_interpolant(Θ,dt,y₀,y₁,k,cache::DP8ConstantCache,idxs,T::Type{Val{0}})
  Θ1 = 1-Θ
  conpar = k[4] + Θ*(k[5] + Θ1*(k[6]+Θ*k[7]))
  y₀ + dt*Θ*(k[1] + Θ1*(k[2] + Θ*(k[3]+Θ1*conpar)))
end

"""

"""
@inline function ode_interpolant!(out,Θ,dt,y₀,y₁,k,cache::DP8Cache,idxs,T::Type{Val{0}})
  Θ1 = 1-Θ
  if out == nothing
    return y₀[idxs] + dt*Θ*(k[1][idxs] + Θ1*(k[2][idxs] + Θ*(k[3][idxs]+Θ1*(k[4][idxs] + Θ*(k[5][idxs] + Θ1*(k[6][idxs]+Θ*k[7][idxs]))))))
  else
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = y₀[i] + dt*Θ*(k[1][i] + Θ1*(k[2][i] + Θ*(k[3][i]+Θ1*(k[4][i] + Θ*(k[5][i] + Θ1*(k[6][i]+Θ*k[7][i]))))))
    end
  end
end
