function newtempKs(Ks25,ΔHs,T)
    newtempKs = copy(Ks25)
    for paramname in keys(newtempKs)
        if paramname in keys(ΔHs)
            ΔH = ΔHs[paramname]
            K25 = Ks25[paramname]
            newtempKs[paramname] = vantHoff(K25,T,ΔH)
        end
    end
    return newtempKs
end
    
function vantHoff(K25,T,ΔH)
    Tkelvin = T+273.15
    R = 8.314e-3 #Units of kJ/mol K
    newK = K25*exp((ΔH/R)*(1/298.15-1/Tkelvin))
end


"""
    IVTmodel

List of parameter objects."""
struct IVTmodel
    parameters
end

"""
    Parameter

Collects information of name, fitting status, prior (hasprior, prior value, and prior variance), and upper and lower bound of parameter for optimization."""
struct Parameter
    name
    isfitted
    hasprior
    value
    variance
    upperbound
    lowerbound
end

Parameter(name,value) = Parameter(name,false,false,value,0,0,0)
Parameter(name,value,upperbound,lowerbound) = Parameter(name,true,false,value,0,upperbound,lowerbound)
Parameter(name,value,upperbound,lowerbound,standarddev) = Parameter(name,true,true,value,standarddev,upperbound,lowerbound)

function reportparametertypes(model::IVTmodel)
    fittedparametercounter = 0
    bayesianparametercounter = 0
    for i in model.parameters
        if i.isfitted
            fittedparametercounter+=1
            if i.hasprior
                bayesianparametercounter+=1
            end
        end
    end
    return (fittedparametercounter,bayesianparametercounter)
end

function fullparameterset(model::IVTmodel,roundparameters)
    namesofvalues::Vector{String} = []
    parametervalues = []
    fittedparametercounter = 1
    for i in model.parameters
        append!(namesofvalues,[i.name])
        if i.isfitted
            append!(parametervalues,[10^(roundparameters[fittedparametercounter])])
            fittedparametercounter+=1
        else
            append!(parametervalues,[i.value])  
        end
    end
    nt = ComponentArray(namedtuple(namesofvalues, parametervalues))
    return nt
end

addtoinputlist!(list,parameter::Parameter) = append!(list,[parameter])

function setupmodel_thermo()#using adjusted values for temp/IS, used for parameter estimation with new system
    #This is our "interface": each line here represents adding a parameter and some associated information
    #format is Parameter(name, bayesian prior, lower bound (useful for optimization), upper bound, bayesian stdev)
    Ks = []
    ΔHs = []
    #Parameters for equilibria


    addtoinputlist!(Ks,Parameter("K_HPPi",10^(-8.94)))#Kern and Davis 1997
    addtoinputlist!(ΔHs,Parameter("K_HPPi",1.63176))#From Smith and Martell Vol 4

    addtoinputlist!(Ks,Parameter("K_NaPPi",10^(-0.21)))#From Smith and Martell Vol 4
    addtoinputlist!(ΔHs,Parameter("K_NaPPi",-5.8576))#From Smith and Martell Vol 4

    addtoinputlist!(Ks,Parameter("K_Na2PPi",10^(0.8)))#From Smith and Martell Vol 4
    addtoinputlist!(ΔHs,Parameter("K_Na2PPi",0))#From Smith and Martell Vol 4, unknown

    addtoinputlist!(Ks,Parameter("K_MgPPi",10^(-5.42)))#Kern and Davis 1997
    addtoinputlist!(ΔHs,Parameter("K_MgPPi",-12.552))#From Smith and Martell Vol 4

    addtoinputlist!(Ks,Parameter("K_Mg2PPi",10^(-2.33)))#Kern and Davis 1997
    addtoinputlist!(ΔHs,Parameter("K_Mg2PPi",0))#Irani 1961

    addtoinputlist!(Ks,Parameter("K_Mg3PPi",10^(10)))#Made up, setting to insignificant value
    addtoinputlist!(ΔHs,Parameter("K_Mg3PPi",-10.8))#Guess from previous Mg addition

    # addtoinputlist!(Ks,Parameter("K_MgPPi2",10^(-7.80)))#Made up, value guess from previous addition of Mg
    # addtoinputlist!(Ks,Parameter("K_Mg3PPiNTP",10^(-7.80)))#Made up, value guess from previous addition of Mg

    addtoinputlist!(Ks,Parameter("K_HMgPPi",10^(-3.05)))#Kern and Davis 1997      
    addtoinputlist!(ΔHs,Parameter("K_HMgPPi",0))#Irani 1961

    addtoinputlist!(Ks,Parameter("K_H2PPi",10^(-6.13)))#Kern and Davis 1997   
    addtoinputlist!(ΔHs,Parameter("K_H2PPi",0.54))#From Smith and Martell Vol 4  

    addtoinputlist!(Ks,Parameter("K_H2MgPPi",7.7e-3))#Young and Davis 1997   

    addtoinputlist!(Ks,Parameter("K_HNTP",3.426e-7))#Alberty and Goldberg 1992
    addtoinputlist!(ΔHs,Parameter("K_HNTP",-6.3))#Alberty and Goldberg 1992

    addtoinputlist!(Ks,Parameter("K_HMgNTP",1.181e-2))#Alberty and Goldberg 1992
    addtoinputlist!(ΔHs,Parameter("K_HMgNTP",-16.9))#Alberty and Goldberg 1992

    addtoinputlist!(Ks,Parameter("K_MgNTP",1.229e-4))#Alberty and Goldberg 1992
    addtoinputlist!(ΔHs,Parameter("K_MgNTP",-22.9))#Alberty and Goldberg 1992

    addtoinputlist!(Ks,Parameter("K_Mg2NTP",2.785e-2))#Alberty and Goldberg 1992
    addtoinputlist!(ΔHs,Parameter("K_Mg2NTP",-10.8))#Alberty and Goldberg 1992

    # addtoinputlist!(Ks,Parameter("K_Mg3NTP",10^(-2.33)))#Made up, designed to be the same as for PPi
    # addtoinputlist!(ΔHs,Parameter("K_Mg3NTP",-10.8))#Made up, designed to be the same as for PPi

    addtoinputlist!(Ks,Parameter("K_Pi",10^(-12.37)))#From Kern and Davis 1995
    addtoinputlist!(Ks,Parameter("K_MgPi",10^(-2.20)))#From Kern and Davis 1995
    addtoinputlist!(Ks,Parameter("K_HPi",10^(-6.92)))#From Kern and Davis 1995
    addtoinputlist!(Ks,Parameter("K_NaPi",10^(-0.77)))#From Kern and Davis 1995
    addtoinputlist!(Ks,Parameter("K_NaNTP",10^(-1.16)))#From Kern and Davis 1995

    addtoinputlist!(Ks,Parameter("K_spermadineNTP",10^(-2.95424)))#From Gunther et al 1994
    addtoinputlist!(Ks,Parameter("K_spermadine2NTP",10^(-2.44715)))#From Gunther et al 1994

    #MgPi solid formation
    addtoinputlist!(Ks,Parameter("Ksp_Mg2PPi",8.0e-14))#Guess, this parameter will be fitted
    addtoinputlist!(ΔHs,Parameter("Ksp_Mg2PPi",8.0e-14))#Guess, this parameter will be fitted

    # addtoinputlist!(Ks,Parameter("MgPiKsp",6.3e-26))#Rough from Taylor et al 1963

    return IVTmodel(Ks), IVTmodel(ΔHs)
end