# Parameter management and thermodynamic constants for the Mg2PPi speciation model.
#
# Equilibrium constants (Ks) are defined at 25°C. Temperature adjustment to 37°C
# is done via the van't Hoff equation using enthalpy values (ΔHs).

"""
    newtempKs(Ks25, ΔHs, T)

Adjust all equilibrium constants from 25°C to temperature T (°C) using the
van't Hoff equation. Only parameters present in both Ks25 and ΔHs are adjusted.
"""
function newtempKs(Ks25, ΔHs, T)
    newtempKs = copy(Ks25)
    for paramname in keys(newtempKs)
        if paramname in keys(ΔHs)
            newtempKs[paramname] = vantHoff(Ks25[paramname], T, ΔHs[paramname])
        end
    end
    return newtempKs
end

"""
    vantHoff(K25, T, ΔH)

Van't Hoff equation: K(T) = K(25°C) * exp[(ΔH/R)(1/298.15 - 1/T)].
ΔH in kJ/mol, T in °C.
"""
function vantHoff(K25, T, ΔH)
    Tkelvin = T + 273.15
    R = 8.314e-3  # kJ/(mol·K)
    return K25 * exp((ΔH / R) * (1/298.15 - 1/Tkelvin))
end

# --- Parameter types for model setup ---

struct IVTmodel
    parameters
end

struct Parameter
    name
    isfitted
    hasprior
    value
    variance
    upperbound
    lowerbound
end

Parameter(name, value) = Parameter(name, false, false, value, 0, 0, 0)
Parameter(name, value, upperbound, lowerbound) = Parameter(name, true, false, value, 0, upperbound, lowerbound)

"""
    fullparameterset(model::IVTmodel, roundparameters)

Extract a named ComponentArray of parameter values from the model.
Fitted parameters are taken from `roundparameters` (log10 scale); others use defaults.
"""
function fullparameterset(model::IVTmodel, roundparameters)
    namesofvalues::Vector{String} = []
    parametervalues = []
    fittedparametercounter = 1
    for i in model.parameters
        append!(namesofvalues, [i.name])
        if i.isfitted
            append!(parametervalues, [10^(roundparameters[fittedparametercounter])])
            fittedparametercounter += 1
        else
            append!(parametervalues, [i.value])
        end
    end
    return ComponentArray(namedtuple(namesofvalues, parametervalues))
end

addtoinputlist!(list, parameter::Parameter) = append!(list, [parameter])

"""
    setupmodel_thermo()

Define all equilibrium constants and enthalpies for the IVT speciation model.
Returns (Ks_model, ΔHs_model) where each is an IVTmodel containing Parameter objects.

Sources for equilibrium constants are noted inline. The Ksp for Mg2PPi·3.5H2O
is set to a placeholder here and overridden in the notebook after fitting.
"""
function setupmodel_thermo()
    Ks = []
    ΔHs = []

    # --- Pyrophosphate (PPi) complexes ---

    # HPPi3- <-> H+ + PPi4-
    addtoinputlist!(Ks, Parameter("K_HPPi", 10^(-8.94)))       # Kern & Davis 1997
    addtoinputlist!(ΔHs, Parameter("K_HPPi", 1.63176))         # Smith & Martell Vol 4

    # NaPPi3- <-> Na+ + PPi4-
    addtoinputlist!(Ks, Parameter("K_NaPPi", 10^(-0.21)))      # Smith & Martell Vol 4
    addtoinputlist!(ΔHs, Parameter("K_NaPPi", -5.8576))        # Smith & Martell Vol 4

    # Na2PPi2- <-> Na+ + NaPPi3-
    addtoinputlist!(Ks, Parameter("K_Na2PPi", 10^(0.8)))       # Smith & Martell Vol 4
    addtoinputlist!(ΔHs, Parameter("K_Na2PPi", 0))             # unknown

    # MgPPi2- <-> Mg2+ + PPi4-
    addtoinputlist!(Ks, Parameter("K_MgPPi", 10^(-5.42)))      # Kern & Davis 1997
    addtoinputlist!(ΔHs, Parameter("K_MgPPi", -12.552))        # Smith & Martell Vol 4

    # Mg2PPi(aq) <-> Mg2+ + MgPPi2-
    addtoinputlist!(Ks, Parameter("K_Mg2PPi", 10^(-2.33)))     # Kern & Davis 1997
    addtoinputlist!(ΔHs, Parameter("K_Mg2PPi", 0))             # Irani 1961

    # Mg3PPi2+ <-> Mg2+ + Mg2PPi(aq)  (set to negligible)
    addtoinputlist!(Ks, Parameter("K_Mg3PPi", 10^(10)))
    addtoinputlist!(ΔHs, Parameter("K_Mg3PPi", -10.8))

    # HMgPPi- <-> Mg2+ + HPPi3-
    addtoinputlist!(Ks, Parameter("K_HMgPPi", 10^(-3.05)))     # Kern & Davis 1997
    addtoinputlist!(ΔHs, Parameter("K_HMgPPi", 0))             # Irani 1961

    # H2PPi2- <-> H+ + HPPi3-
    addtoinputlist!(Ks, Parameter("K_H2PPi", 10^(-6.13)))      # Kern & Davis 1997
    addtoinputlist!(ΔHs, Parameter("K_H2PPi", 0.54))           # Smith & Martell Vol 4

    # H2MgPPi(aq) <-> Mg2+ + H2PPi2-
    addtoinputlist!(Ks, Parameter("K_H2MgPPi", 7.7e-3))        # Young & Davis 1997

    # --- NTP complexes ---

    # HNTP3- <-> H+ + NTP4-
    addtoinputlist!(Ks, Parameter("K_HNTP", 3.426e-7))         # Alberty & Goldberg 1992
    addtoinputlist!(ΔHs, Parameter("K_HNTP", -6.3))

    # HMgNTP- <-> Mg2+ + HNTP3-
    addtoinputlist!(Ks, Parameter("K_HMgNTP", 1.181e-2))       # Alberty & Goldberg 1992
    addtoinputlist!(ΔHs, Parameter("K_HMgNTP", -16.9))

    # MgNTP2- <-> Mg2+ + NTP4-
    addtoinputlist!(Ks, Parameter("K_MgNTP", 1.229e-4))        # Alberty & Goldberg 1992
    addtoinputlist!(ΔHs, Parameter("K_MgNTP", -22.9))

    # Mg2NTP(aq) <-> Mg2+ + MgNTP2-
    addtoinputlist!(Ks, Parameter("K_Mg2NTP", 2.785e-2))       # Alberty & Goldberg 1992
    addtoinputlist!(ΔHs, Parameter("K_Mg2NTP", -10.8))

    # --- Phosphate (Pi) complexes ---

    # PO43- <-> H+ + ...
    addtoinputlist!(Ks, Parameter("K_Pi", 10^(-12.37)))        # Kern & Davis 1995
    addtoinputlist!(Ks, Parameter("K_MgPi", 10^(-2.20)))       # Kern & Davis 1995
    addtoinputlist!(Ks, Parameter("K_HPi", 10^(-6.92)))        # Kern & Davis 1995
    addtoinputlist!(Ks, Parameter("K_NaPi", 10^(-0.77)))       # Kern & Davis 1995

    # --- Na-NTP and spermidine-NTP complexes ---

    addtoinputlist!(Ks, Parameter("K_NaNTP", 10^(-1.16)))      # Kern & Davis 1995
    addtoinputlist!(Ks, Parameter("K_spermadineNTP", 10^(-2.95424)))   # Gunther et al 1994
    addtoinputlist!(Ks, Parameter("K_spermadine2NTP", 10^(-2.44715)))  # Gunther et al 1994

    # --- Solubility product (fitted in notebook) ---

    addtoinputlist!(Ks, Parameter("Ksp_Mg2PPi", 8.0e-14))      # placeholder, overridden after fitting
    addtoinputlist!(ΔHs, Parameter("Ksp_Mg2PPi", 8.0e-14))     # placeholder, overridden after fitting

    return IVTmodel(Ks), IVTmodel(ΔHs)
end
