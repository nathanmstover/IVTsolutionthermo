
"""
    relative(f, r; sp)

Generate relative coordinates for plotting.""" 
function relative(f, r; sp)
    p = plot!()
    lims = f(p[sp])
    return lims[1] + r * (lims[2]-lims[1])
end
relativex(r; sp::Int=1) = relative(Plots.xlims, r; sp=sp)
relativey(r; sp::Int=1) = relative(Plots.ylims, r; sp=sp)

function plotCryscatter(crystaldata;color1 = :snow4, color2 = :tomato)
    crystal = Bool.(crystaldata[:,6])
    plt = scatter(1e-3*crystaldata[crystal,1],1e-3*crystaldata[crystal,2],xscale = :log,yscale = :log,label = "Precipitation Observed",legend = :outertop,xlabel = "Magnesium (M)", ylabel = "PPi (M)", color = color1, markershape = :xcross,markerstrokewidth = 3)
    scatter!(1e-3*crystaldata[.!crystal,1],1e-3*crystaldata[.!crystal,2],label = "No Precipitation Observed", color = color2,markerstrokewidth = 0)
    return plt
end

function solubilitycurve!(plt,params; ATP = 0, spermadine = 1e-8, color = :black, lw = 2,label = "",reduceorder = true, legend = :none)  
    n = 100
    
    #Drawing "upper arm"
    MgCl2 = 5e-3
    PPisol = []
    Mgsol = []
    PPistart = 10 .^(LinRange(-3,0,n))
    nuguess = 0
    for (i, PPi) in enumerate(PPistart)
        if reduceorder
            (Mgend,PPiend,res) = getsolubilityRO(params,PPi,MgCl2,ATP,spermadinetot = spermadine, tol = 1e-6)
        else
            (Mgend,PPiend,res) = getsolubility(params,PPi,MgCl2,ATP, nuguess = nuguess, tol = 1e-6)
            nuguess = res.zero[6]
        end
        if  converged(res)
            append!(PPisol,[PPiend])
            append!(Mgsol,[Mgend])
        end
    end
    
    #Drawing "lower arm"
    PPi = 5e-3
    PPisol_lower = []
    Mgsol_lower = []
    Mgstart = 10 .^(LinRange(log10(5e-3),-0,n))
    nuguess = 0
    for (i, MgCl2) in enumerate(Mgstart)
        if reduceorder
            (Mgend,PPiend,res) = getsolubilityRO(params,PPi,MgCl2,ATP,spermadinetot = spermadine)
        else
            (Mgend,PPiend,res) = getsolubility(params,PPi,MgCl2,ATP, nuguess = nuguess)
            nuguess = res.zero[6]
        end
        if converged(res)
            append!(PPisol_lower,[PPiend])
            append!(Mgsol_lower,[Mgend])
        end
    end
    
    plot!(plt, Mgsol,PPisol,yscale = :log, label = label, color = color, linewidth = lw)
    plot!(plt, Mgsol_lower,PPisol_lower,xscale = :log, label = "", color = color, linewidth = lw, legend = legend)
end

function solubilitydata!(plt, solubilitymatrix; color = :black)
        scatter!(plt,1e-3 .*solubilitymatrix[:,3],1e-3 .*solubilitymatrix[:,4],label = "",markersize = 3,color = color)
end

function solubilitybool!(plt,Mg_PPi,formedcrystal;color = :black)
    scatter!(plt,1e-3*Mg_PPi[formedcrystal,1],1e-3*Mg_PPi[formedcrystal,2],color = color,markershape = :star5)
    scatter!(plt,1e-3*Mg_PPi[.!formedcrystal,1],1e-3*Mg_PPi[.!formedcrystal,2],color = color,markershape = :circle)
end