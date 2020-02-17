abstract type DAEInitializationAlgorithm end

struct BrownFullBasicInit <: DAEInitializationAlgorithm end

function initialize_dae!(integrator, u, du, differential_vars, alg::BrownFullBasicInit, ::Val{true})
	@unpack p, t, f = integrator

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
end

function initialize_dae!(integrator, u, du, differential_vars, alg::BrownFullBasicInit, ::Val{false})
	@unpack p, t, f = integrator

	nlequation = (dx,x) -> begin
		for i in 1:length(x)
			if differential_vars[i]
				du[i] = x[i]
			else
				u[i] = x[i]
			end
		end
		du .= f(u, p, t)
	end

	r = nlsolve(nlequation, zero(u))

	for i in 1:length(u)
		if differential_vars[i]
			du[i] = r.zero[i]
		else
			u[i] = r.zero[i]
		end
	end
end
