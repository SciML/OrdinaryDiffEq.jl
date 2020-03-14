abstract type DAEInitializationAlgorithm end

struct DefaultInit <: DAEInitializationAlgorithm end
struct NoInit <: DAEInitializationAlgorithm end

struct ShampineCollocationInit{T} <: DAEInitializationAlgorithm
	initdt::T
end
ShampineCollocationInit() = ShampineCollocationInit(nothing)

struct BrownFullBasicInit{T} <: DAEInitializationAlgorithm
	abstol::T
end
BrownFullBasicInit() = BrownFullBasicInit(1e-10)

## Notes

#=
differential_vars = [any(!iszero,x) for x in eachcol(M)]

A column should be zero for an algebraic variable, since that means that the
derivative term doesn't show up in any equations (i.e. is an algebraic variable).
The rows are not necessarily non-zero, for example a flux condition between two
differential variables. But if it's a condition that doesn't involve the algebraic
variable, then the system is not Index 1!

=#

## Default algorithms

function initialize_dae!(integrator, prob::ODEProblem, u, du,
						 alg::DefaultInit, x::Val{true})
	initialize_dae!(integrator, prob, u, du,
					ShampineCollocationInit(), x)
end

function initialize_dae!(integrator, prob::ODEProblem, u, du,
						 alg::DefaultInit, x::Val{false})
	initialize_dae!(integrator, prob, u, du,
					ShampineCollocationInit(), x)
end

function initialize_dae!(integrator, prob::DAEProblem, u, du,
						 alg::DefaultInit, x::Val{false})
	initialize_dae!(integrator, prob, u, du,
					BrownFullBasicInit(), x)
end

function initialize_dae!(integrator, prob::DAEProblem, u, du,
						 alg::DefaultInit, x::Val{true})
	initialize_dae!(integrator, prob, u, du,
					BrownFullBasicInit(), x)
end

## NoInit

function initialize_dae!(integrator, prob::ODEProblem, u, du,
						 alg::NoInit, x::Val{true})
end

function initialize_dae!(integrator, prob::ODEProblem, u, du,
						 alg::NoInit, x::Val{false})
end

function initialize_dae!(integrator, prob::DAEProblem, u, du,
						 alg::NoInit, x::Val{false})
end

function initialize_dae!(integrator, prob::DAEProblem, u, du,
						 alg::NoInit, x::Val{true})
end

## ShampineCollocationInit

#=
The method:

du = (u-u0)/h
Solve for `u`

=#

function initialize_dae!(integrator, prob::ODEProblem, u0, du0,
						 alg::ShampineCollocationInit, ::Val{true})

	@unpack p, t, f = integrator
 	M = integrator.f.mass_matrix
	dtmax = integrator.opts.dtmax
 	update_coefficients!(M,u,p,t)
	tmp = first(get_tmp_cache(integrator))

	dt = t != 0 ? min(t/1000,dtmax) : dtmax # Haven't implemented norm reduction

 	nlequation! = function (out,u)
 		#M * (u-u0)/dt - f(u,p,t)
		@. tmp = (u - u0)/dt
		mul!(out,M,tmp)
		f(tmp,u,p,t)
		out .-= tmp
		nothing
 	end

	differential_vars = [any(!iszero,x) for x in eachcol(M)]
	f(tmp,u0,p,t)
	tmp .= (differential_vars .== false) .* tmp

 	integrator.opts.internalnorm(tmp,t) <= integrator.opts.abstol && return

	integrator.u .= nlsolve(nlequation!, u0).zero
	recursivecopy!(integrator.uprev,integrator.u)
	if alg_extrapolates(integrator.alg)
		recursivecopy!(integrator.uprev2,integrator.uprev)
	end

end

function initialize_dae!(integrator, prob::ODEProblem, u0, du0,
						 alg::ShampineCollocationInit, ::Val{false})

	@unpack p, t, f = integrator
	M = integrator.f.mass_matrix
	dtmax = integrator.opts.dtmax
	update_coefficients!(M,u,p,t)

	dt = t != 0 ? min(t/1000,dtmax/10) : dtmax # Haven't implemented norm reduction

	differential_vars = [any(!iszero,x) for x in eachcol(M)]
	du = f(u0,p,t)
	resid = du[!differential_vars]

	integrator.opts.internalnorm(resid,t) <= integrator.opts.abstol && return

	nlequation_oop = function (u)
		M * (u-u0)/dt - f(u,p,t)
	end

	nlequation! = (out,u) -> out .= nlequation_oop(u)

	integrator.u = nlsolve(nlequation!, u0).zero
	integrator.uprev = integrator.u
	if alg_extrapolates(integrator.alg)
		integrator.uprev2 = integrator.uprev
	end
end

function initialize_dae!(integrator, prob::DAEProblem, u0, du0,
						 alg::ShampineCollocationInit, ::Val{true})

	@unpack p, t, f = integrator
	dtmax = integrator.opts.dtmax
 	update_coefficients!(M,u,p,t)
	tmp = get_tmp_cache(integrator)[1]
	resid = get_tmp_cache(integrator)[2]

	dt = t != 0 ? min(t/1000,dtmax) : dtmax # Haven't implemented norm reduction

 	nlequation! = function (out,u)
 		#M * (u-u0)/dt - f(u,p,t)
		@. tmp = (u - u0)/dt
		f(out,tmp,u,p,t)
		nothing
 	end

	nlequation!(tmp,u0)
	f(resid,tmp,u0,p,t)
 	integrator.opts.internalnorm(resid,t) <= integrator.opts.abstol && return

	integrator.u .= nlsolve(nlequation!, u0).zero
	recursivecopy!(integrator.uprev,integrator.u)
	if alg_extrapolates(integrator.alg)
		recursivecopy!(integrator.uprev2,integrator.uprev)
	end

end

function initialize_dae!(integrator, prob::DAEProblem, u0, du0,
						 alg::ShampineCollocationInit, ::Val{false})

	@unpack p, t, f = integrator
	dtmax = integrator.opts.dtmax
	update_coefficients!(M,u,p,t)

	dt = t != 0 ? min(t/1000,dtmax/10) : dtmax # Haven't implemented norm reduction

	nlequation_oop = function (u)
		f((u-u0)/dt,u,p,t)
	end

	nlequation! = (out,u) -> out .= nlequation_oop(u)

	resid = f(du,u0,p,t)
	integrator.opts.internalnorm(resid,t) <= integrator.opts.abstol && return

	integrator.u = nlsolve(nlequation!, u0).zero
	integrator.uprev = integrator.u
	if alg_extrapolates(integrator.alg)
		integrator.uprev2 = integrator.uprev
	end
end

## BrownFullBasic

#=
The method:

Keep differential variables constant
Solve for the algebraic variables

=#

function initialize_dae!(integrator, prob::ODEProblem, u, du,
						 alg::BrownFullBasicInit, ::Val{true})
	@unpack p, t, f = integrator
	differential_vars = [any(!iszero,x) for x in eachcol(M)]

	f(tmp,u0,p,t)
	tmp .= (differential_vars .== false) .* tmp

	integrator.opts.internalnorm(tmp,t) <= alg.abstol && return
	alg_u = @view u[!differential_vars]

	nlequation = (out, x) -> begin
		alg_u .= x
		du = f(u, p, t)
		out .= @view du[!differential_vars]
	end

	r = nlsolve(nlequation, u[!differential_vars])
	alg_u .= r.zero

	recursivecopy!(integrator.uprev,integrator.u)
	if alg_extrapolates(integrator.alg)
		recursivecopy!(integrator.uprev2,integrator.uprev)
	end

	return
end

function initialize_dae!(integrator, prob::ODEProblem, u0, du0,
						 alg::BrownFullBasicInit, ::Val{false})
	@unpack p, t, f = integrator

	du = f(u0,p,t)
	resid = du[!differential_vars]

	integrator.opts.internalnorm(resid,t) <= alg.abstol && return
	alg_u = @view u[!differential_vars]

	if u0 isa Number && du0 isa Number
		# This doesn't fix static arrays!
		u = [u0]
		du = [_du]
	else
		u = u0
		du = _du
	end

	nlequation = (out,x) -> begin
		alg_u .= x
		du = f(u,p,t)
		out .= du[!differential_vars]
	end

	r = nlsolve(nlequation, u0[!differential_vars])
	alg_u .= x

	if _u isa Number && _du isa Number
		# This doesn't fix static arrays!
		integrator.u = first(u)
		integrator.cache.du = first(du)
	else
		integrator.u = u
		integrator.cache.du = du
	end

	integrator.uprev = integrator.u
	if alg_extrapolates(integrator.alg)
		integrator.uprev2 = integrator.uprev
	end

	return
end

function initialize_dae!(integrator, prob::DAEProblem, u, du,
						 alg::BrownFullBasicInit, ::Val{true})
	@unpack p, t, f = integrator
	differential_vars = prob.differential_vars

	tmp = get_tmp_cache(integrator)[1]
	f(tmp, du, u, p, t)

	if integrator.opts.internalnorm(tmp,t) <= alg.abstol
		return
	elseif differential_vars === nothing
		error("differential_vars must be set for DAE initialization to occur. Either set consistent initial conditions or differential_vars")
	end

	nlequation = (out, x) -> begin
		for i in 1:length(x)
			if differential_vars[i]
				du[i] = x[i]
			else
				u[i] = x[i]
			end
		end
		f(out, du, u, p, t)
	end

	r = nlsolve(nlequation, zero(u))

	for i in 1:length(u)
		if differential_vars[i]
			du[i] = r.zero[i]
		else
			u[i] = r.zero[i]
		end
	end

	recursivecopy!(integrator.uprev,integrator.u)
	if alg_extrapolates(integrator.alg)
		recursivecopy!(integrator.uprev2,integrator.uprev)
	end

	return
end

function initialize_dae!(integrator, prob::DAEProblem, _u, _du,
						 differential_vars, alg::BrownFullBasicInit, ::Val{false})
	@unpack p, t, f = integrator
	differential_vars = prob.differential_vars

	if integrator.opts.internalnorm(f(_du, _u, p, t),t) <= alg.abstol
		return
	elseif differential_vars === nothing
		error("differential_vars must be set for DAE initialization to occur. Either set consistent initial conditions or differential_vars")
	end

	if _u isa Number && _du isa Number
		# This doesn't fix static arrays!
		u = [_u]
		du = [_du]
	else
		u = _u
		du = _du
	end

	nlequation = (out,x) -> begin
		for i in 1:length(x)
			if differential_vars[i]
				du[i] = x[i]
			else
				u[i] = x[i]
			end
		end
		out .= f(du, u, p, t)
	end

	r = nlsolve(nlequation, zero(u))

	for i in 1:length(u)
		if differential_vars[i]
			du[i] = r.zero[i]
		else
			u[i] = r.zero[i]
		end
	end

	if _u isa Number && _du isa Number
		# This doesn't fix static arrays!
		integrator.u = first(u)
		integrator.cache.du = first(du)
	else
		integrator.u = u
		integrator.cache.du = du
	end

	integrator.uprev = integrator.u
	if alg_extrapolates(integrator.alg)
		integrator.uprev2 = integrator.uprev
	end

	return
end
