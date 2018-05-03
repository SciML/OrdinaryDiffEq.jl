function ϕ_and_ϕstar!(cache, dy, next_point, last_idx)
	@unpack grid_points, ϕstar_nm1, k = cache
	ϕ_n = zeros(k)
	ϕstar_n = zeros(k)
	β = zeros(k)
	for i = 0:(k)-1
		if i == 0
			β[(i)+1] = 1
			ϕ_n[(i)+1] = dy
			ϕstar_n[(i)+1] = dy
		else
			β[(i)+1] = β[i] * (next_point - grid_points[last_idx-i+1])/(grid_points[last_idx] - grid_points[last_idx-i])
			ϕ_n[(i)+1] = ϕ_n[i] - ϕstar_nm1[i]
			ϕstar_n[(i)+1] = β[(i)+1] * ϕ_n[(i)+1]
		end
	end
	return ϕ_n, ϕstar_n 
end

function ϕ_np1!(ϕstar_n, dy, k)
	ϕ_np1 = zeros(Float64,k+1)
    for i = 1:(k)+1
    	if i == 1
    		ϕ_np1[i] = dy
    	else
    		ϕ_np1[i] = ϕ_np1[i-1] - ϕstar_n[i-1]
    	end
    end
    return ϕ_np1
end

function g_coefs!(cache, dt, next_point, last_idx)
	@unpack grid_points, k = cache
	c = zeros(Float64, k, k)
	g = zeros(Float64, k)
	for i = 1:k
		for q = 1:(k-(i-1))
			if i == 1
				c[(i),q] = 1/q
			elseif i == 2
				c[(i),q] = 1/q/(q+1)
			else
				c[(i),q] = c[i-1,q] - c[i-1,q+1] * (dt)/(next_point - grid_points[last_idx-i+1+1])
			end
		end
		g[(i)] = c[(i),1]
	end
	return g
end
