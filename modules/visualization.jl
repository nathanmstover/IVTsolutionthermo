# Plotting functions for Mg2PPi phase diagrams.

"""Compute a position relative to the current plot limits (0 to 1 scale)."""
function relative(f, r; sp)
    p = plot!()
    lims = f(p[sp])
    return lims[1] + r * (lims[2] - lims[1])
end
relativex(r; sp::Int=1) = relative(Plots.xlims, r; sp=sp)
relativey(r; sp::Int=1) = relative(Plots.ylims, r; sp=sp)

"""
    plotCryscatter(crystaldata; color1, color2)

Scatter plot of experimental crystal/no-crystal observations.
Column 6 of crystaldata is a boolean (1 = precipitation observed).
"""
function plotCryscatter(crystaldata; color1 = :snow4, color2 = :tomato)
    crystal = Bool.(crystaldata[:, 6])
    plt = scatter(1e-3*crystaldata[crystal, 1], 1e-3*crystaldata[crystal, 2],
        xscale=:log, yscale=:log, label="Precipitation Observed", legend=:outertop,
        xlabel="Magnesium (M)", ylabel="PPi (M)",
        color=color1, markershape=:xcross, markerstrokewidth=3)
    scatter!(1e-3*crystaldata[.!crystal, 1], 1e-3*crystaldata[.!crystal, 2],
        label="No Precipitation Observed", color=color2, markerstrokewidth=0)
    return plt
end

"""
    solubilitycurve!(plt, params; ATP, spermadine, color, lw, label, legend)

Draw the Mg2PPi solubility boundary on an existing plot by sweeping PPi (upper arm)
and Mg (lower arm) and solving for the equilibrium phase boundary at each point.
"""
function solubilitycurve!(plt, params; ATP = 0, spermadine = 1e-8, color = :black, lw = 2, label = "", legend = :none)
    n = 100

    # Upper arm: sweep PPi at fixed MgCl2
    MgCl2 = 5e-3
    PPisol = []
    Mgsol = []
    for PPi in 10 .^ LinRange(-3, 0, n)
        (Mgend, PPiend, res) = getsolubilityRO(params, PPi, MgCl2, ATP, spermadinetot=spermadine, tol=1e-6)
        if converged(res)
            push!(PPisol, PPiend)
            push!(Mgsol, Mgend)
        end
    end

    # Lower arm: sweep Mg at fixed PPi
    PPi = 5e-3
    PPisol_lower = []
    Mgsol_lower = []
    for MgCl2 in 10 .^ LinRange(log10(5e-3), 0, n)
        (Mgend, PPiend, res) = getsolubilityRO(params, PPi, MgCl2, ATP, spermadinetot=spermadine)
        if converged(res)
            push!(PPisol_lower, PPiend)
            push!(Mgsol_lower, Mgend)
        end
    end

    plot!(plt, Mgsol, PPisol, yscale=:log, label=label, color=color, linewidth=lw)
    plot!(plt, Mgsol_lower, PPisol_lower, xscale=:log, label="", color=color, linewidth=lw, legend=legend)
end
