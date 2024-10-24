
# SOLPS2ctrl.jl 

```@contents
Pages = ["index.md"]
Depth = 3
```

This repository serves as the top most workflow manager with helpful utilities to use other repositories in this project. Following steps are supported/planned in near future:

* Loading SOLPS outputs into IMAS data dictionary format
* Loading equilibrium (that the SOLPS mesh was based on) into IMAS DD format
* Extrapolating results into core and far SOL region if required
* Running synthetic diagnostics on them
* Performing system identification (to be added)
* Designing and tuning linear causal and model predictive controllers (to be added)

## Documentation of other repositories in this project

### [IMASggd.jl](https://projecttorreypines.github.io/IMASggd.jl/stable)

### [SOLPS2imas.jl](https://projecttorreypines.github.io/SOLPS2imas.jl/stable)

### [FusionSyntheticDiagnostics.jl](https://projecttorreypines.github.io/FusionSyntheticDiagnostics.jl/stable)

SOLPS2ctrl is registered with public repository [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl/). For installation:

## Installation

```
using Pkg
Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
Pkg.Registry.add("General")
Pkg.add("SOLPS2ctrl")
```

## Top file handling functions

```@docs
find_files_in_allowed_folders
geqdsk_to_imas!
preparation
```

## Repairing/filling out partial equilibrium files

Tools for repairing/filling out partial equilibrium files.

Some of the added fields may not be totally accurate, so it is recommended to
use this tool mainly for test cases, as a utility. For a real equilibrium,
problems should be fixed properly.

```@docs
add_rho_to_equilibrium!
check_rho_1d
```

## Extrapolations

Utilities for extrapolating profiles

### Core profile extrapolations

```@docs
extrapolate_core
fill_in_extrapolated_core_profile!
```

### Edge profiles extrapolations

These functions have not been fully tested and/or supported yet.

```@docs
mesh_psi_spacing
cached_mesh_extension!
```

## Unit conversion utilities

```@docs
gas_unit_converter
```