abstract type DAEInitializationAlgorithm end

struct BrownFullBasicInit{T} <: DAEInitializationAlgorithm
	abstol::T
end
BrownFullBasicInit() = BrownFullBasicInit(1e-10)

function initialize_dae!(integrator, u, du, differential_vars, alg::BrownFullBasicInit, ::Val{true})
	@unpack p, t, f = integrator

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

	return
end

function initialize_dae!(integrator, _u, _du, differential_vars, alg::BrownFullBasicInit, ::Val{false})
	@unpack p, t, f = integrator

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

	return
end
