@recipe function f(sim::ConvergenceSimulation)
    if any((x) -> ndims(x) > 1, values(sim.errors)) #Monte Carlo
        vals = [x isa AbstractVector ? x : vec(mean(x, 1))
                for x in values(sim.errors)]
    else #Deterministic
        vals = [x for x in values(sim.errors)]
    end
    seriestype --> :path
    label -->
    reshape([string(key) for key in keys(sim.errors)], 1, length(keys(sim.errors)))
    xguide --> "Convergence Axis"
    yguide --> "Error"
    linewidth --> 3
    xscale --> :log10
    yscale --> :log10
    marker --> :auto
    sim.convergence_axis, vals
end

@recipe function f(shoot::Shootout)
    seriestype --> :bar
    legend := false
    xguide --> "Algorithms"
    yguide --> "Efficiency"
    shoot.names, shoot.effs
end

@recipe function f(wp::WorkPrecision)
    seriestype --> :path
    label --> wp.name
    linewidth --> 3
    yguide --> "Time (s)"
    xguide --> "Error ($(wp.error_estimate))"
    xscale --> :log10
    yscale --> :log10
    markershape --> :auto
    errors = getproperty(wp.errors, wp.error_estimate)
    return errors, wp.times
end

function get_val_from_wp(wp::WorkPrecision, key::Symbol)
    if key == :abstols
        return wp.abstols
    elseif key == :reltols
        return wp.reltols
    elseif key == :times
        return wp.times
    elseif key == :dts
        return wp.dts
    elseif !isnothing(wp.stats) && key in propertynames(wp.stats[1])
        return [getproperty(s, key) for s in wp.stats]
    elseif key in propertynames(wp.errors)
        return getproperty(wp.errors, key)
    else
        throw(ArgumentError("Key $key does not specify a valid property of WorkPrecision"))
    end
end

function key_to_label(key::Symbol)
    if key == :abstols
        return "Absolute tolerance"
    elseif key == :reltols
        return "Relative tolerance"
    elseif key == :times
        return "Time (s)"
    elseif key == :dts
        return "Δt"
    elseif key in ALL_ERRORS
        return "Error ($key)"
        # Some DEStats cases copied from SciMLBase: src/solutions/ode_solutions.jl#L43-L56
    elseif key == :nf
        return "Number of function evaluations"
    elseif key == :nf2
        return "Number of function 2 evaluations"
    elseif key == :nw
        return "Number of W matrix evaluations"
    elseif key == :nsolve
        return "Number of linear solves"
    elseif key == :njacs
        return "Number of Jacobians created"
    elseif key == :nnonliniter
        return "Number of nonlinear solver iterations"
    elseif key == :nnonlinconvfail
        return "Number of nonlinear solver convergence failures"
    elseif key == :ncondition
        return "Number of rootfind condition calls"
    elseif key == :naccept
        return "Number of accepted steps"
    elseif key == :nreject
        return "Number of rejected steps"
    elseif key == :maxeig
        return "Maximum eigenvalue recorded"
    else
        return String(key)
    end
end

@recipe function f(wp_set::WorkPrecisionSet;
        x = wp_set.error_estimate,
        y = :times,
        view = :benchmark,
        color = nothing,
        tags = nothing,
        include_tags = nothing,
        exclude_tags = nothing,
        reference_tags = nothing,
        reference_style = (linestyle = :dash, linewidth = 1, alpha = 0.5))
    # Apply tag-based filtering if specified
    if tags !== nothing || include_tags !== nothing || exclude_tags !== nothing
        filtered = wp_set
        if tags !== nothing
            filtered = filter_by_tags(filtered, tags...)
        end
        if include_tags !== nothing
            extra = filter_by_tags(wp_set, include_tags...)
            seen = Set(wp.name for wp in filtered.wps)
            extra_wps = [wp for wp in extra.wps if wp.name ∉ seen]
            if !isempty(extra_wps)
                all_wps = vcat(filtered.wps, extra_wps)
                all_names = [wp.name for wp in all_wps]
                filtered = WorkPrecisionSet(
                    all_wps, length(all_wps), filtered.abstols, filtered.reltols,
                    filtered.prob, filtered.setups, all_names, filtered.error_estimate,
                    filtered.numruns)
            end
        end
        if exclude_tags !== nothing
            filtered = exclude_by_tags(filtered, exclude_tags...)
        end
        wp_set = filtered
    end

    if view == :benchmark
        if reference_tags !== nothing
            ref_tags = reference_tags isa Symbol ? (reference_tags,) : reference_tags
            ref_indices = findall(wp -> any(t -> t in wp.tags, ref_tags), wp_set.wps)
            main_indices = setdiff(1:length(wp_set.wps), ref_indices)

            if !isempty(ref_indices)
                @series begin
                    seriestype --> :path
                    linewidth --> get(reference_style, :linewidth, 1)
                    linestyle --> get(reference_style, :linestyle, :dash)
                    seriesalpha --> get(reference_style, :alpha, 0.5)
                    xscale --> :log10
                    yscale --> :log10
                    markershape --> :auto
                    xguide --> key_to_label(x)
                    yguide --> key_to_label(y)
                    legend --> :outerright
                    ref_xs = [get_val_from_wp(wp_set.wps[i], x) for i in ref_indices]
                    ref_ys = [get_val_from_wp(wp_set.wps[i], y) for i in ref_indices]
                    label --> reshape([wp_set.wps[i].name for i in ref_indices], 1, length(ref_indices))
                    ref_xs, ref_ys
                end
            end

            if !isempty(main_indices)
                seriestype --> :path
                linewidth --> 3
                xscale --> :log10
                yscale --> :log10
                markershape --> :auto
                xguide --> key_to_label(x)
                yguide --> key_to_label(y)
                legend --> :outerright
                main_xs = [get_val_from_wp(wp_set.wps[i], x) for i in main_indices]
                main_ys = [get_val_from_wp(wp_set.wps[i], y) for i in main_indices]
                label --> reshape([wp_set.wps[i].name for i in main_indices], 1, length(main_indices))
                return main_xs, main_ys
            end
        else
            seriestype --> :path
            linewidth --> 3
            xscale --> :log10
            yscale --> :log10
            markershape --> :auto
            xs = [get_val_from_wp(wp, x) for wp in wp_set.wps]
            ys = [get_val_from_wp(wp, y) for wp in wp_set.wps]
            xguide --> key_to_label(x)
            yguide --> key_to_label(y)
            legend --> :outerright
            label --> reshape(wp_set.names, 1, length(wp_set))
            return xs, ys
        end
    elseif view == :dt_convergence
        idts = filter(i -> haskey(wp_set.setups[i], :dts), 1:length(wp_set))
        length(idts) > 0 ||
            throw(ArgumentError("Convergence with respect to Δt requires runs with fixed time steps"))
        dts = Vector{Any}(undef, 0)
        errors = Vector{Any}(undef, 0)
        ps = Vector{Any}(undef, 0)
        convs = Vector{Any}(undef, 0)
        for i in idts
            push!(dts, wp_set.setups[i][:dts])
            push!(errors, getproperty(wp_set[i].errors, wp_set.error_estimate))
            lc, p = [one.(dts[end]) log.(dts[end])] \ log.(errors[end])
            push!(ps, p)
            push!(convs, exp(lc) * dts[end] .^ p)
        end
        names = wp_set.names[idts] .*
                map(p -> " (Δtᵖ order p=$(round(p, sigdigits=2)))", ps)
        if color === nothing
            color = reshape(1:length(idts), 1, :)
        end
        @series begin
            seriestype --> :path
            linestyle --> :dash
            linewidth --> 1
            color --> color
            xscale --> :log10
            yscale --> :log10
            markersize --> 0
            label --> nothing
            return dts, convs
        end
        @series begin
            seriestype --> :path
            linewidth --> 3
            color --> color
            yguide --> "Error ($(wp_set.error_estimate))"
            xguide --> "Δt"
            xscale --> :log10
            yscale --> :log10
            markershape --> :auto
            label --> reshape(names, 1, length(idts))
            return dts, errors
        end
    else
        throw(ArgumentError("view argument `$view` not implemented"))
    end
end

@recipe function f(tab::ODERKTableau; dx = 1 / 100, dy = 1 / 100, order_star = false,
        embedded = false)
    xlims = get(plotattributes, :xlims, (-6, 1))
    ylims = get(plotattributes, :ylims, (-5, 5))
    x = xlims[1]:dx:xlims[2]
    y = ylims[1]:dy:ylims[2]

    if order_star
        f = (u, v) -> abs(stability_region(u + v * im, tab; embedded = embedded) /
                          exp(u + v * im)) < 1
    else
        f = (u, v) -> abs(stability_region(u + v * im, tab; embedded = embedded)) < 1
    end
    seriestype --> :contour
    fill --> true
    colorbar --> false
    seriescolor --> :grays
    aspect_ratio --> 1
    x, y, f
end
