function getsolubilityRO(parameters,PPi,MgCl2,ATP; spermadinetot = 1e-8, tol = 1e-10)
    TrisBuffer = 40e-3
    NTPtot = ATP
    Mgtot = MgCl2
    Buffertot = TrisBuffer
    PPitot = PPi
    Pitot = 0
    Natot = 3*PPi+4*ATP
    Cltot = 2*Mgtot+0.5*TrisBuffer
    OActot = 0
    result = Mg2PPiphaseequilibriumRO(parameters, NTPtot, Mgtot, Buffertot, PPitot, Pitot, Natot, Cltot, OActot, spermadinetot; tol = tol)
    nu = result.zero[5]
    Mg = Mgtot - 2*nu
    PPi = PPitot-nu
    return Mg, PPi, result
end



"""
    getfreeconcentrations(params::AbstractArray{T}, NTPtot, Mgtot, Buffertot, PPitot, Pitot)  where {T<:Real}

Take parameters and total concentration of ions, return free concentrations of H and Mg. Performs solving of nonlinear system of equations in log space. Only used for initialization of DAE.

"""
function Mg2PPiphaseequilibriumRO(params::AbstractArray{T}, NTPtot, Mgtot, Buffertot, PPitot, Pitot, Natot, Cltot, OActot, spermadinetot; buffer_pka = 8.1, init = log.([1e-7,max(1e-5,Mgtot-2*PPitot),Natot,spermadinetot/2,1]), tol = 1e-8)  where {T<:Real}
    #Defining Parameters
    initialtotalconcentrations = (NTPtot=NTPtot, Mgtot=Mgtot, Buffertot=Buffertot, PPitot = PPitot, Pitot = Pitot, Natot = Natot, Cltot = Cltot, OActot = OActot, spermadinetot = spermadinetot)
    guessfreeconcentrations = zeros(T,5)
    guessfreeconcentrations += init
    logsolvedfreeconcentrations = nlsolve((F,x)->Mg2PPiphaseresidualRO!(F, x, initialtotalconcentrations, params,buffer_pka = buffer_pka), guessfreeconcentrations,ftol = tol)
    solvedfreeconcentrations = [exp(x) for x in logsolvedfreeconcentrations.zero]
    return logsolvedfreeconcentrations
end

"""
    initialconcentrationresidual!(F, x, ptot, param)

Generate residual for nonlinear solving of free ion concentrations. Only used for initialization of DAE.

"""
function Mg2PPiphaseresidualRO!(F, x, ptot, param;buffer_pka = 8.1, RNAtotalbases = 0)
    H = exp(x[1])
    Mg = exp(x[2])
    Na = exp(x[3])
    spermadine = exp(x[4])
    nu = x[5]
    Mgtot = max(ptot.Mgtot - 2*nu,1e-9)
    PPitot = max(ptot.PPitot - nu,1e-9)
    (ions,(Chargebalance,Mgbalance,Nabalance,spermadinebalance)) = speciationmodelRO(param, ptot.NTPtot, Mgtot, ptot.Buffertot, PPitot, ptot.Pitot, ptot.Natot, ptot.spermadinetot, H, Mg, Na, ptot.Cltot, ptot.OActot, spermadine, buffer_pka = buffer_pka)
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
    #spermadine mass balance
    F[5] = thermodrivingforce
    nothing
end

function getsupersaturationRO(parameters,PPi,MgCl2,ATP; spermadinetot = 1e-8, tol = 1e-8)
    TrisBuffer = 40e-3
    NTPtot = ATP
    Mgtot = MgCl2
    Buffertot = TrisBuffer
    PPitot = PPi
    Pitot = 0
    Natot = 3*PPi+4*ATP
    Cltot = 2*Mgtot+0.5*TrisBuffer
    OActot = 0
    (H, Mg, Na, spermadine) = getfreeconcentrationsRO(parameters, NTPtot, Mgtot, Buffertot, PPitot, Pitot, Natot, Cltot, OActot, spermadinetot; tol = tol)
    (ionspecies,balances) = speciationmodelRO(parameters, NTPtot, Mgtot, Buffertot, PPitot, Pitot, Natot, spermadinetot, H, Mg, Na, Cltot, OActot, spermadine)
    return ionspecies[4]
end

"""
    getfreeconcentrations(params::AbstractArray{T}, NTPtot, Mgtot, Buffertot, PPitot, Pitot)  where {T<:Real}

Take parameters and total concentration of ions, return free concentrations of H and Mg. Performs solving of nonlinear system of equations in log space. Only used for initialization of DAE.

"""
function getfreeconcentrationsRO(params::AbstractArray{T}, NTPtot, Mgtot, Buffertot, PPitot, Pitot, Natot, Cltot, OActot, spermadinetot; buffer_pka = 8.1, tol = 1e-8, init = log.([1e-7,0.0001,Natot,spermadinetot/2]))  where {T<:Real}
    #Defining Parameters
    initialtotalconcentrations = (NTPtot=NTPtot, Mgtot=Mgtot, Buffertot=Buffertot, PPitot = PPitot, Pitot = Pitot, Natot = Natot, Cltot = Cltot, OActot = OActot, spermadinetot = spermadinetot)
    guessfreeconcentrations = zeros(T,4)
    guessfreeconcentrations += init
    logsolvedfreeconcentrations = nlsolve((F,x)->initialconcentrationresidualRO!(F, x, initialtotalconcentrations, params,buffer_pka = buffer_pka), guessfreeconcentrations,ftol = tol)
    solvedfreeconcentrations = [exp(x) for x in logsolvedfreeconcentrations.zero]
    return solvedfreeconcentrations
end

function getfreeconcentrationsRO(params::AbstractArray{Float64}, NTPtot, Mgtot, Buffertot::T, PPitot, Pitot, Natot::Float64, Cltot, OActot, spermadinetot; buffer_pka = 8.1, tol = 1e-8, init = log.([1e-7,0.0001,Natot,spermadinetot/2]))  where {T<:Real}
    #Defining Parameters
    initialtotalconcentrations = (NTPtot=NTPtot, Mgtot=Mgtot, Buffertot=Buffertot, PPitot = PPitot, Pitot = Pitot, Natot = Natot, Cltot = Cltot, OActot = OActot, spermadinetot = spermadinetot)
    guessfreeconcentrations = zeros(T,4)
    guessfreeconcentrations += init
    logsolvedfreeconcentrations = nlsolve((F,x)->initialconcentrationresidualRO!(F, x, initialtotalconcentrations, params,buffer_pka = buffer_pka), guessfreeconcentrations,ftol = tol)
    solvedfreeconcentrations = [exp(x) for x in logsolvedfreeconcentrations.zero]
    return solvedfreeconcentrations
end

"""
    initialconcentrationresidual!(F, x, ptot, param)

Generate residual for nonlinear solving of free ion concentrations. Only used for initialization of DAE.

"""
function initialconcentrationresidualRO!(F, x, ptot, param;buffer_pka = 8.1, RNAtotalbases = 0)
    H = exp(x[1])
    Mg = exp(x[2])
    Na = exp(x[3])
    spermadine = exp(x[4])
    (ions,(Chargebalance,Mgbalance,Nabalance,spermadinebalance)) = speciationmodelRO(param, ptot.NTPtot, ptot.Mgtot, ptot.Buffertot, ptot.PPitot, ptot.Pitot, ptot.Natot, ptot.spermadinetot, H, Mg, Na, ptot.Cltot, ptot.OActot, spermadine, buffer_pka = buffer_pka)
    #Algebraic Equations
    # H mass balance
    F[1] = Chargebalance
    # Mg mass balance
    F[2] = Mgbalance
    #Na mass balance
    F[3] = Nabalance
    #spermadine mass balance
    F[4] = spermadinebalance
    nothing
end

"""
    speciationmodel(param, NTPtot, Mgtot, Buffertot, PPitot, Pitot, H, Mg)

Take parameters and total concentration of ions, return concentrations of ionic species and ion equilibrium residuals for use in rate calculations and DAE solving.

"""
function speciationmodelRO(param, NTPtot, Mgtot, Buffertot, PPitot, Pitot, Natot, spermadinetot, H, Mg, Na, Cl, OAc, spermadine;buffer_pka = 8.1)

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

    #Dimensionless concentrations for PPi
    HPPi_dimless = H*param.K_HPPi^(-1)
    HMgPPi_dimless = HPPi_dimless*Mg*param.K_HMgPPi^(-1)
    H2PPi_dimless = HPPi_dimless*H*param.K_H2PPi^(-1)
    H2MgPPi_dimless = H2PPi_dimless*Mg*param.K_H2MgPPi^(-1)
    MgPPi_dimless = Mg*param.K_MgPPi^(-1)
    Mg2PPi_dimless = MgPPi_dimless*Mg*param.K_Mg2PPi^(-1)
    Mg3PPi_dimless = Mg2PPi_dimless*Mg*param.K_Mg3PPi^(-1)
    NaPPi_dimless = Na*param.K_NaPPi^(-1)
    Na2PPi_dimless = NaPPi_dimless*Na*param.K_Na2PPi^(-1)

    #Dimensionless concentrations for Pi
    MgPi_dimless = Mg*param.K_MgPi^(-1)
    HPi_dimless = H*param.K_HPi^(-1)
    NaPi_dimless = Na*param.K_NaPi^(-1)

    #Solve for Anion Concentrations
    NTP = NTPtot/(1 + HNTP_dimless + HMgNTP_dimless + MgNTP_dimless + Mg2NTP_dimless + NaNTP_dimless + spermadineNTP_dimless + spermadine2NTP_dimless)
    PPi = PPitot/(1 + MgPPi_dimless + Mg2PPi_dimless + Mg3PPi_dimless + HPPi_dimless + HMgPPi_dimless + H2PPi_dimless + H2MgPPi_dimless + NaPPi_dimless + Na2PPi_dimless)
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

    #Redimensionalizing concentrations for PPi
    HPPi = HPPi_dimless*PPi
    HMgPPi = HMgPPi_dimless*PPi
    H2PPi = H2PPi_dimless*PPi
    H2MgPPi = H2MgPPi_dimless*PPi
    MgPPi = MgPPi_dimless*PPi
    Mg2PPi = Mg2PPi_dimless*PPi
    Mg3PPi = Mg3PPi_dimless*PPi
    NaPPi = NaPPi_dimless*PPi
    Na2PPi = Na2PPi_dimless*PPi
    
    #Redimensionalizing concentrations for Pi
    MgPi = MgPi_dimless*Pi
    HPi = HPi_dimless*Pi
    NaPi = NaPi_dimless*Pi

    #Algebraic Equations
    #Change balance
    negativecharge = OAc+Cl+OH+4*NTP+4*PPi+2*Pi+3*HNTP+HMgNTP+2*MgNTP+3*NaNTP+3*HPPi+HMgPPi+2*H2PPi+2*MgPPi+HPi+NaPi+1*spermadineNTP+3*NaPPi+2*Na2PPi
    positivecharge = H+HBuffer+Na+2*Mg+3*spermadine+2*Mg3PPi+2*spermadine2NTP

    Mg2PPisupersaturation = Mg^2*PPi/param.Ksp_Mg2PPi
    
    Chargebalance  = 1-(negativecharge)/(positivecharge)
    # Mg mass balance
    Mgbalance = 1 - (1/(Mgtot)) * (Mg + MgPPi + HMgPPi + MgNTP + 2*Mg2NTP + 2*Mg2PPi + 3*Mg3PPi +HMgNTP + H2MgPPi + MgPi)
    # Na mass balance
    Nabalance = 1 - (1/(Natot)) * (Na + NaPi + NaNTP + NaPPi + 2*Na2PPi)
    # Spermadine mass balance
    spermadinebalance = 1 - (1/(2*spermadinetot)) * (spermadine + spermadineNTP + spermadine2NTP)

    ionspecies = (Mg,MgNTP,MgPPi,Mg2PPisupersaturation,H)
    balances = (Chargebalance,Mgbalance,Nabalance,spermadinebalance)
    return (ionspecies,balances)
end