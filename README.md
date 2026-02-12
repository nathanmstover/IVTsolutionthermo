# IVTsolutionthermo

Ionic speciation model for predicting Mg2PPi solubility in IVT reaction mixtures. Generates the phase diagrams in Figure 4b,c of:

> S. Ahmadi, N. M. Stover, K. Ganko, F. B. Kayser, R. D. Braatz, A. S. Myerson, "DNA nanoflowers assemble via template-induced crystallization during mRNA synthesis."

## What this code does

The model implements the thermodynamic framework described in Section 2.1 of the paper. It computes the supersaturation of Mg2PPi (Eq. 1) by solving for the free ion concentrations [Mg2+] and [PPi4-] from a speciation network that includes complexes with ATP, spermidine, Na+, H+, and Tris buffer (Figure 4a, Table S1). The single fitted parameter is the solubility product Ksp of Mg2PPi·3.5H2O, fit to ICP-MS solubility measurements at 25 and 37°C. Temperature dependence is handled via the van't Hoff equation.

The notebook reproduces the model predictions in Figure 4b,c — phase diagrams comparing predicted solubility boundaries against experimental precipitation observations at 4 mM and 20 mM ATP.

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
