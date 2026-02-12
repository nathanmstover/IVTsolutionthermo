# IVTsolutionthermo

Ionic speciation model for predicting Mg2PPi solubility in IVT reaction mixtures. Generates the phase diagrams in Figure 4b,c of:

> S. Ahmadi, N. M. Stover, K. Ganko, F. B. Kayser, R. D. Braatz, A. S. Myerson, "DNA nanoflowers assemble via template-induced crystallization during mRNA synthesis."

## What this code does

The model solves for the supersaturation of Mg2PPi·3.5H2O as defined in Eq. 1. The notebook reproduces Figure 4b,c — phase diagrams comparing predicted solubility boundaries against experimental precipitation observations at 4 mM and 20 mM ATP. Free ion concentrations [Mg2+] and [PPi4-] are determined by solving the speciation equilibria listed in Table S1 (SI Section 5), which govern how Mg, PPi, ATP, spermidine, Na+, and H+ partition among complexes in the IVT mixture (see Figure 4a). Equilibrium constants at 37°C are obtained from the 25°C values in Table S1 via the van't Hoff relation. The single fitted parameter is Ksp of Mg2PPi·3.5H2O, fit to the ICP-MS solubility measurements in Figure S3.


## Running

Requires Julia and IJulia (for notebooks). Install dependencies:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Then run `Mg2PPi_phasediagramdemo.ipynb`.

## Contents

- `modules/` — Speciation model, solver, parameter definitions, and plotting functions
- `data/` — Experimental solubility and precipitation observations (CSV)
- `figures/` — Output phase diagrams

## Dependencies

ComponentArrays, NamedTupleTools, NLsolve, Plots, CSV, DataFrames
