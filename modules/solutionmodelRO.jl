# Solubility model for Mg2PPi in IVT reaction conditions.
#
# Solves for the Mg2PPi phase boundary by finding the set of free ion concentrations
# (H+, Mg2+, Na+, spermidine) and extent of precipitation (nu) that simultaneously
# satisfy mass balances, charge balance, and the solubility product constraint.
#
# "RO" = reduced order: PPi concentration is determined from the other species via
# mass balance rather than being solved as an independent free variable.

"""
    getsolubilityRO(parameters, PPi, MgCl2, ATP; spermadinetot, tol)

Compute the equilibrium Mg and PPi concentrations after Mg2PPi precipitation.
Assumes 40 mM Tris buffer at pH 8.0. Returns (Mg_eq, PPi_eq, result).
"""
function getsolubilityRO(parameters, PPi, MgCl2, ATP; spermadinetot = 1e-8, tol = 1e-10)
    TrisBuffer = 40e-3
    NTPtot = ATP
    Mgtot = MgCl2
    Buffertot = TrisBuffer
    PPitot = PPi
    Pitot = 0
    Natot = 3*PPi + 4*ATP       # Na from Na4PPi and Na4ATP stock solutions
    Cltot = 2*Mgtot + 0.5*TrisBuffer  # Cl from MgCl2 and Tris-HCl
    OActot = 0
    result = Mg2PPiphaseequilibriumRO(parameters, NTPtot, Mgtot, Buffertot, PPitot, Pitot, Natot, Cltot, OActot, spermadinetot; tol = tol)
    nu = result.zero[5]          # moles of Mg2PPi precipitated per liter
    Mg = Mgtot - 2*nu
    PPi = PPitot - nu
    return Mg, PPi, result
end

"""
    Mg2PPiphaseequilibriumRO(params, NTPtot, Mgtot, ...)

Solve the nonlinear system for phase equilibrium. Free variables (solved in log space):
  x[1..4] = log of free [H+], [Mg2+], [Na+], [spermidine]
  x[5]    = nu (extent of precipitation, linear scale)

The system is at equilibrium when all mass/charge balances are satisfied and
the supersaturation sigma = 0 (i.e., Mg^2 * PPi / Ksp = 1).
"""
function Mg2PPiphaseequilibriumRO(params::AbstractArray{T}, NTPtot, Mgtot, Buffertot, PPitot, Pitot, Natot, Cltot, OActot, spermadinetot; buffer_pka = 8.1, init = log.([1e-7, max(1e-5, Mgtot-2*PPitot), Natot, spermadinetot/2, 1]), tol = 1e-8) where {T<:Real}
    initialtotalconcentrations = (NTPtot=NTPtot, Mgtot=Mgtot, Buffertot=Buffertot, PPitot=PPitot, Pitot=Pitot, Natot=Natot, Cltot=Cltot, OActot=OActot, spermadinetot=spermadinetot)
    guessfreeconcentrations = zeros(T, 5)
    guessfreeconcentrations += init
    logsolvedfreeconcentrations = nlsolve((F,x) -> Mg2PPiphaseresidualRO!(F, x, initialtotalconcentrations, params, buffer_pka=buffer_pka), guessfreeconcentrations, ftol=tol)
    return logsolvedfreeconcentrations
end

"""
    Mg2PPiphaseresidualRO!(F, x, ptot, param)

Residual function for the phase equilibrium solver. Equations:
  F[1] = charge balance
  F[2] = Mg mass balance
  F[3] = Na mass balance
  F[4] = spermidine mass balance
  F[5] = log(supersaturation) = 0 at equilibrium
"""
function Mg2PPiphaseresidualRO!(F, x, ptot, param; buffer_pka = 8.1)
    H = exp(x[1])
    Mg = exp(x[2])
    Na = exp(x[3])
    spermadine = exp(x[4])
    nu = x[5]
    # Adjust total concentrations for amount precipitated as Mg2PPi
    Mgtot = max(ptot.Mgtot - 2*nu, 1e-9)
    PPitot = max(ptot.PPitot - nu, 1e-9)
    (ions, (Chargebalance, Mgbalance, Nabalance, spermadinebalance)) = speciationmodelRO(param, ptot.NTPtot, Mgtot, ptot.Buffertot, PPitot, ptot.Pitot, ptot.Natot, ptot.spermadinetot, H, Mg, Na, ptot.Cltot, ptot.OActot, spermadine, buffer_pka=buffer_pka)
    thermodrivingforce = log(ions[4])  # log(supersaturation), = 0 at equilibrium
    F[1] = Chargebalance
    F[2] = Mgbalance
    F[3] = Nabalance
    F[4] = spermadinebalance
    F[5] = thermodrivingforce
    nothing
end

"""
    speciationmodelRO(param, NTPtot, Mgtot, Buffertot, PPitot, Pitot, Natot,
                      spermadinetot, H, Mg, Na, Cl, OAc, spermadine)

Given free ion concentrations and total concentrations, compute:
  1. All complex species concentrations via equilibrium expressions
  2. Supersaturation of Mg2PPi: sigma = [Mg2+]^2 * [PPi4-] / Ksp  (Eq. 1 in paper)
  3. Residuals for charge balance and mass balances (Mg, Na, spermidine)

PPi is determined from its mass balance (reduced-order approach) rather than
being a free variable. Complexes considered:
  - PPi: HPPi, H2PPi, MgPPi, Mg2PPi, Mg3PPi, HMgPPi, H2MgPPi, NaPPi, Na2PPi
  - NTP: HNTP, MgNTP, Mg2NTP, HMgNTP, NaNTP, spermidine-NTP, spermidine2-NTP
  - Pi:  MgPi, HPi, NaPi
  - Buffer: HBuffer (protonated Tris)
"""
function speciationmodelRO(param, NTPtot, Mgtot, Buffertot, PPitot, Pitot, Natot, spermadinetot, H, Mg, Na, Cl, OAc, spermadine; buffer_pka = 8.1)

    # Buffer speciation (Tris)
    HBuffer_dimless = H * 10^buffer_pka
    Buffer = Buffertot / (1 + HBuffer_dimless)
    HBuffer = HBuffer_dimless * Buffer

    # NTP complexes (dimensionless ratios relative to free NTP4-)
    HNTP_dimless = H * param.K_HNTP^(-1)
    HMgNTP_dimless = HNTP_dimless * Mg * param.K_HMgNTP^(-1)
    MgNTP_dimless = Mg * param.K_MgNTP^(-1)
    Mg2NTP_dimless = MgNTP_dimless * Mg * param.K_Mg2NTP^(-1)
    NaNTP_dimless = Na * param.K_NaNTP^(-1)
    spermadineNTP_dimless = spermadine * param.K_spermadineNTP^(-1)
    spermadine2NTP_dimless = spermadine * spermadineNTP_dimless * param.K_spermadine2NTP^(-1)

    # PPi complexes (dimensionless ratios relative to free PPi4-)
    HPPi_dimless = H * param.K_HPPi^(-1)
    HMgPPi_dimless = HPPi_dimless * Mg * param.K_HMgPPi^(-1)
    H2PPi_dimless = HPPi_dimless * H * param.K_H2PPi^(-1)
    H2MgPPi_dimless = H2PPi_dimless * Mg * param.K_H2MgPPi^(-1)
    MgPPi_dimless = Mg * param.K_MgPPi^(-1)
    Mg2PPi_dimless = MgPPi_dimless * Mg * param.K_Mg2PPi^(-1)
    Mg3PPi_dimless = Mg2PPi_dimless * Mg * param.K_Mg3PPi^(-1)
    NaPPi_dimless = Na * param.K_NaPPi^(-1)
    Na2PPi_dimless = NaPPi_dimless * Na * param.K_Na2PPi^(-1)

    # Pi complexes (dimensionless ratios relative to free PO4 3-)
    MgPi_dimless = Mg * param.K_MgPi^(-1)
    HPi_dimless = H * param.K_HPi^(-1)
    NaPi_dimless = Na * param.K_NaPi^(-1)

    # Solve for free anion concentrations from mass balances
    NTP = NTPtot / (1 + HNTP_dimless + HMgNTP_dimless + MgNTP_dimless + Mg2NTP_dimless + NaNTP_dimless + spermadineNTP_dimless + spermadine2NTP_dimless)
    PPi = PPitot / (1 + MgPPi_dimless + Mg2PPi_dimless + Mg3PPi_dimless + HPPi_dimless + HMgPPi_dimless + H2PPi_dimless + H2MgPPi_dimless + NaPPi_dimless + Na2PPi_dimless)
    Pi = Pitot / (1 + MgPi_dimless + HPi_dimless + NaPi_dimless)
    OH = 10^-14 / H

    # Absolute concentrations of NTP complexes
    HNTP = HNTP_dimless * NTP
    HMgNTP = HMgNTP_dimless * NTP
    MgNTP = MgNTP_dimless * NTP
    Mg2NTP = Mg2NTP_dimless * NTP
    NaNTP = NaNTP_dimless * NTP
    spermadineNTP = spermadineNTP_dimless * NTP
    spermadine2NTP = spermadine2NTP_dimless * NTP

    # Absolute concentrations of PPi complexes
    HPPi = HPPi_dimless * PPi
    HMgPPi = HMgPPi_dimless * PPi
    H2PPi = H2PPi_dimless * PPi
    H2MgPPi = H2MgPPi_dimless * PPi
    MgPPi = MgPPi_dimless * PPi
    Mg2PPi = Mg2PPi_dimless * PPi
    Mg3PPi = Mg3PPi_dimless * PPi
    NaPPi = NaPPi_dimless * PPi
    Na2PPi = Na2PPi_dimless * PPi

    # Absolute concentrations of Pi complexes
    MgPi = MgPi_dimless * Pi
    HPi = HPi_dimless * Pi
    NaPi = NaPi_dimless * Pi

    # Charge balance residual
    negativecharge = OAc + Cl + OH + 4*NTP + 4*PPi + 2*Pi + 3*HNTP + HMgNTP + 2*MgNTP + 3*NaNTP + 3*HPPi + HMgPPi + 2*H2PPi + 2*MgPPi + HPi + NaPi + 1*spermadineNTP + 3*NaPPi + 2*Na2PPi
    positivecharge = H + HBuffer + Na + 2*Mg + 3*spermadine + 2*Mg3PPi + 2*spermadine2NTP

    # Supersaturation: sigma = [Mg2+]^2 * [PPi4-] / Ksp  (Eq. 1)
    Mg2PPisupersaturation = Mg^2 * PPi / param.Ksp_Mg2PPi

    # Balance residuals (each should be ~0 at solution)
    Chargebalance = 1 - (negativecharge) / (positivecharge)
    Mgbalance = 1 - (1/Mgtot) * (Mg + MgPPi + HMgPPi + MgNTP + 2*Mg2NTP + 2*Mg2PPi + 3*Mg3PPi + HMgNTP + H2MgPPi + MgPi)
    Nabalance = 1 - (1/Natot) * (Na + NaPi + NaNTP + NaPPi + 2*Na2PPi)
    spermadinebalance = 1 - (1/(2*spermadinetot)) * (spermadine + spermadineNTP + spermadine2NTP)

    ionspecies = (Mg, MgNTP, MgPPi, Mg2PPisupersaturation, H)
    balances = (Chargebalance, Mgbalance, Nabalance, spermadinebalance)
    return (ionspecies, balances)
end
