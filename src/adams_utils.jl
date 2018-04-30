function ϕ_and_ϕstar!(cache, dy)
	@unpack grid_points, ϕstar_nm1, k = cache
	ϕ_n = zeros(k)
	ϕstar_n = zeros(k)
	β = zeros(k)
	for i = 1:k
		if i == 1
			β[i] = 1
			ϕ_n[i] = dy
			ϕstar_n[i] = dy
		else
			β[i] = β[i-1] * (grid_points[end] - grid_points[end-i+1])/(grid_points[end-1] - grid_points[end-i])
			ϕ_n[i] = ϕ_n[i-1] - ϕstar_nm1[i-1]
			ϕstar_n[i] = β[i] * ϕ_n[i]
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

function g_coefs!(cache)
	@unpack grid_points, k = cache
	c = zeros(Float64, k+1, k+1)
	g = zeros(Float64, k+1)
	for i = 1:(k)+1
		for q = 1:(k-i+1)
			if i == 1
				c[i,q] = 1/q
			elseif i == 2
				c[i,q] = 1/q/(q+1)
			else
				c[i,q] = c[i-1,q] - c[i-1,q+1] * (grid_points[end] - grid_points[end-1])/(grid_points[end] - grid_points[end-i+1])
			end
		end
		g[i] = c[i,1]
	end
	return g
end
