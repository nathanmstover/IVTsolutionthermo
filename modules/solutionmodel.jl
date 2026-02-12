function getsolubility(parameters,PPi,MgCl2,ATP; kwargs...)
    TrisBuffer = 40e-3
    NTPtot = ATP
    Mgtot = MgCl2
    Buffertot = TrisBuffer
    PPitot = PPi
    Pitot = 0
    Natot = 3*PPi+4*ATP
    Cltot = 2*Mgtot+0.5*TrisBuffer
    OActot = 0
    spermadinetot = 1e-8
    result = Mg2PPiphaseequilibrium(parameters, NTPtot, Mgtot, Buffertot, PPitot, Pitot, Natot, Cltot, OActot, spermadinetot; kwargs...)
    nu = result.zero[6]
    Mg = Mgtot - 2*nu
    PPi = PPitot-nu
    return Mg, PPi, result
end

"""
    initialconcentrationresidual!(F, x, ptot, param)

Generate residual for nonlinear solving of free ion concentrations. Only used for initialization of DAE.

"""
function initialconcentrationresidual!(F, x, ptot, param;buffer_pka = 8.1, RNAtotalbases = 0)
    H = exp(x[1])
    Mg = exp(x[2])
    Na = exp(x[3])
    spermadine = exp(x[4])
    PPi = exp(x[5])
    (ions,(Chargebalance,Mgbalance,Nabalance,spermadinebalance, PPibalance)) = speciationmodel(param, ptot.NTPtot, ptot.Mgtot, ptot.Buffertot, ptot.PPitot, ptot.Pitot, ptot.Natot, ptot.spermadinetot, H, Mg, Na, ptot.Cltot, ptot.OActot, spermadine, PPi, buffer_pka = buffer_pka)
    #Algebraic Equations
    # H mass balance
    F[1] = Chargebalance
    # Mg mass balance
    F[2] = Mgbalance
    #Na mass balance
    F[3] = Nabalance
    #spermadine mass balance
    F[4] = spermadinebalance
    #PPi mass balance
    F[5] = PPibalance
    nothing
end

"""
    initialconcentrationresidual!(F, x, ptot, param)

Generate residual for nonlinear solving of free ion concentrations. Only used for initialization of DAE.

"""
function Mg2PPiphaseresidual!(F, x, ptot, param;buffer_pka = 8.1, RNAtotalbases = 0)
    H = exp(x[1])
    Mg = exp(x[2])
    Na = exp(x[3])
    spermadine = exp(x[4])
    PPi = exp(x[5])
    nu = x[6]
    Mgtot = max(ptot.Mgtot - 2*nu,1e-9)
    PPitot = max(ptot.PPitot - nu,1e-9)
    (ions,(Chargebalance,Mgbalance,Nabalance,spermadinebalance, PPibalance)) = speciationmodel(param, ptot.NTPtot, Mgtot, ptot.Buffertot, PPitot, ptot.Pitot, ptot.Natot, ptot.spermadinetot, H, Mg, Na, ptot.Cltot, ptot.OActot, spermadine, PPi, buffer_pka = buffer_pka)
    thermodrivingforce = log(ions[4])
    #Algebraic Equations
    # H mass balance
    F[1] = Chargebalance
    # Mg mass balance
    F[2] = Mgbalance
    #Na mass balance
    F[3] = Nabalance
    #spermadine mass balance
    F[4] = spermadinebalance
    #PPi mass balance
    F[5] = PPibalance
    #supersaturation
    F[6] = thermodrivingforce
    nothing
end

"""
    getfreeconcentrations(params::AbstractArray{T}, NTPtot, Mgtot, Buffertot, PPitot, Pitot)  where {T<:Real}

Take parameters and total concentration of ions, return free concentrations of H and Mg. Performs solving of nonlinear system of equations in log space. Only used for initialization of DAE.

"""
function Mg2PPiphaseequilibrium(params::AbstractArray{T}, NTPtot, Mgtot, Buffertot, PPitot, Pitot, Natot, Cltot, OActot, spermadinetot; buffer_pka = 8.1, init = log.([1e-7,max(1e-5,Mgtot-2*PPitot),Natot,spermadinetot/2,max(1e-5,PPitot-0.5*Mgtot),1]), tol = 1e-8, nuguess = 0)  where {T<:Real}
    #Defining Parameters
    initialtotalconcentrations = (NTPtot=NTPtot, Mgtot=Mgtot, Buffertot=Buffertot, PPitot = PPitot, Pitot = Pitot, Natot = Natot, Cltot = Cltot, OActot = OActot, spermadinetot = spermadinetot)
    guessfreeconcentrations = zeros(T,6)
    guessfreeconcentrations += init
    guessfreeconcentrations[end] += nuguess
    logsolvedfreeconcentrations = nlsolve((F,x)->Mg2PPiphaseresidual!(F, x, initialtotalconcentrations, params,buffer_pka = buffer_pka), guessfreeconcentrations,ftol = tol)
    solvedfreeconcentrations = [exp(x) for x in logsolvedfreeconcentrations.zero]
    return logsolvedfreeconcentrations
end

"""
    getfreeconcentrations(params::AbstractArray{T}, NTPtot, Mgtot, Buffertot, PPitot, Pitot)  where {T<:Real}

Take parameters and total concentration of ions, return free concentrations of H and Mg. Performs solving of nonlinear system of equations in log space. Only used for initialization of DAE.

"""
function getfreeconcentrations(params::AbstractArray{T}, NTPtot, Mgtot, Buffertot, PPitot, Pitot, Natot, Cltot, OActot, spermadinetot; buffer_pka = 8.1, init = log.([1e-7,0.0001,Natot,spermadinetot/2,max(1e-5,PPitot-0.5*Mgtot)]))  where {T<:Real}
    #Defining Parameters
    initialtotalconcentrations = (NTPtot=NTPtot, Mgtot=Mgtot, Buffertot=Buffertot, PPitot = PPitot, Pitot = Pitot, Natot = Natot, Cltot = Cltot, OActot = OActot, spermadinetot = spermadinetot)
    guessfreeconcentrations = zeros(T,5)
    guessfreeconcentrations += init
    logsolvedfreeconcentrations = nlsolve((F,x)->initialconcentrationresidual!(F, x, initialtotalconcentrations, params,buffer_pka = buffer_pka), guessfreeconcentrations,ftol = 1e-8)
    solvedfreeconcentrations = [exp(x) for x in logsolvedfreeconcentrations.zero]
    return solvedfreeconcentrations
end


function getfreeconcentrations(params::AbstractArray{Float64}, NTPtot, Mgtot, Buffertot::T, PPitot, Pitot, Natot::Float64, Cltot, OActot, spermadinetot; buffer_pka = 8.1, init = log.([1e-7,0.0001,Natot,spermadinetot/2,max(1e-5,PPitot-0.5*Mgtot)]))  where {T<:Real}
    #Defining Parameters
    initialtotalconcentrations = (NTPtot=NTPtot, Mgtot=Mgtot, Buffertot=Buffertot, PPitot = PPitot, Pitot = Pitot, Natot = Natot, Cltot = Cltot, OActot = OActot, spermadinetot = spermadinetot)
    guessfreeconcentrations = zeros(T,5)
    guessfreeconcentrations += init
    logsolvedfreeconcentrations = nlsolve((F,x)->initialconcentrationresidual!(F, x, initialtotalconcentrations, params,buffer_pka = buffer_pka), guessfreeconcentrations,ftol = 1e-8)
    solvedfreeconcentrations = [exp(x) for x in logsolvedfreeconcentrations.zero]
    return solvedfreeconcentrations
end

"""
    speciationmodel(param, NTPtot, Mgtot, Buffertot, PPitot, Pitot, H, Mg)

Take parameters and total concentration of ions, return concentrations of ionic species and ion equilibrium residuals for use in rate calculations and DAE solving.

"""
function speciationmodel(param, NTPtot, Mgtot, Buffertot, PPitot, Pitot, Natot, spermadinetot, H, Mg, Na, Cl, OAc, spermadine, PPi;buffer_pka = 8.1)

    #Calculating HBuffer concentrations - we do this out of order of the others since it is needed for MgRNA binding
    HBuffer_dimless = H*10^buffer_pka
    Buffer = Buffertot/(1+HBuffer_dimless)
    HBuffer = HBuffer_dimless*Buffer

    #Dimensionless concentrations for NTP
    HNTP_dimless = H*param.K_HNTP^(-1)
    HMgNTP_dimless = HNTP_dimless*Mg*param.K_HMgNTP^(-1)
    MgNTP_dimless = Mg*param.K_MgNTP^(-1)
    Mg2NTP_dimless = MgNTP_dimless*Mg*param.K_Mg2NTP^(-1)
    NaNTP_dimless = Na*param.K_NaNTP^(-1)
    spermadineNTP_dimless = spermadine*param.K_spermadineNTP^(-1)
    spermadine2NTP_dimless = spermadine*spermadineNTP_dimless*param.K_spermadine2NTP^(-1)


    #concentrations for PPi
    HPPi = H*PPi*param.K_HPPi^(-1)
    HMgPPi = HPPi*Mg*param.K_HMgPPi^(-1)
    H2PPi = HPPi*H*param.K_H2PPi^(-1)
    H2MgPPi = H2PPi*Mg*param.K_H2MgPPi^(-1)
    MgPPi = Mg*PPi*param.K_MgPPi^(-1)
    Mg2PPi = MgPPi*Mg*param.K_Mg2PPi^(-1)
    MgPPi2 = Mg*PPi^2*param.K_MgPPi2^(-1)

    #Dimensionless concentrations for Pi
    MgPi_dimless = Mg*param.K_MgPi^(-1)
    HPi_dimless = H*param.K_HPi^(-1)
    NaPi_dimless = Na*param.K_NaPi^(-1)

    #Solve for Anion Concentrations
    NTP = NTPtot/(1 + HNTP_dimless + HMgNTP_dimless + MgNTP_dimless + Mg2NTP_dimless + NaNTP_dimless + spermadineNTP_dimless + spermadine2NTP_dimless)
    Pi = Pitot/(1+MgPi_dimless+HPi_dimless+NaPi_dimless)
    OH = 10^-14/H

    #Redimensionalizing concentrations for NTP
    HNTP = HNTP_dimless*NTP
    HMgNTP = HMgNTP_dimless*NTP
    MgNTP = MgNTP_dimless*NTP
    Mg2NTP = Mg2NTP_dimless*NTP
    NaNTP = NaNTP_dimless*NTP
    spermadineNTP = spermadineNTP_dimless*NTP
    spermadine2NTP = spermadine2NTP_dimless*NTP

    #Redimensionalizing concentrations for Pi
    MgPi = MgPi_dimless*Pi
    HPi = HPi_dimless*Pi
    NaPi = NaPi_dimless*Pi

    #Algebraic Equations
    #Change balance
    #total charge equations (correct but slower)
    negativecharge = OAc+Cl+OH+4*NTP+4*PPi+2*Pi+3*HNTP+HMgNTP+2*MgNTP+3*NaNTP+3*HPPi+HMgPPi+2*H2PPi+2*MgPPi+HPi+NaPi+1*spermadineNTP-6*MgPPi2
    positivecharge = H+HBuffer+Na+2*Mg+1*spermadine

    #totalfreeions = Cl+OH+NTP+PPi+Pi+HNTP+HMgNTP+MgNTP+NaNTP+HPPi+HMgPPi+H2PPi+MgPPi+HPi+NaPi+H+HBuffer+Na+Mg
    Mg2PPisupersaturation = Mg^2*PPi/param.Ksp_Mg2PPi
    MgPisupersaturation = ((Pi/H)/param.K_Pi)^2*Mg^3/param.MgPiKsp
    
    #PPi only (good for pyrophosphate solutions)
    # negativecharge = Cl+OH+4*PPi+3*HPPi+2*H2PPi
    # positivecharge = H+HBuffer+Na

    #NTP only (good for NTP solutions)
    # negativecharge = Cl+OH+4*NTP+3*HNTP+3*NaNTP
    # positivecharge = H+HBuffer+Na

    #Basic (good for strongly buffered solutions where pH doesn't matter much)
    # negativecharge = Cl+OH
    # positivecharge = H+HBuffer

    Chargebalance  = 1-(negativecharge)/(positivecharge)
    # Mg mass balance
    Mgbalance = 1 - (1/(Mgtot)) * (Mg + MgPPi + HMgPPi + MgNTP + 2*Mg2NTP + 2*Mg2PPi + HMgNTP + H2MgPPi + MgPi + MgPPi2)
    # Na mass balance
    Nabalance = 1 - (1/(Natot)) * (Na + NaPi + NaNTP)
    # Spermadine mass balance
    spermadinebalance = 1 - (1/(2*spermadinetot)) * (spermadine + spermadineNTP + spermadine2NTP)
    # PPi mass balance 
    #PPibalance = 1 - (1/(PPitot)) * (HPPi + HMgPPi + H2PPi + H2MgPPi + MgPPi + Mg2PPi + 2*MgPPi2 + MgPPiNTP)
    PPibalance = PPitot - (HPPi + HMgPPi + H2PPi + H2MgPPi + MgPPi + Mg2PPi + 2*MgPPi2)
    ionspecies = (Mg,MgNTP,MgPPi,Mg2PPisupersaturation,H,MgPisupersaturation)
    balances = (Chargebalance,Mgbalance,Nabalance,spermadinebalance,PPibalance)
    return (ionspecies,balances)
end