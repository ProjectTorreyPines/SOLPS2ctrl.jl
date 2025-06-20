
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

## Installation

```
using Pkg
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

## System identification and modeling

```@docs
offset_scale
unscale_unoffset
system_id
system_id_optimal_input_cond
model_evolve
```

## Controllers

```@docs
state_prediction_matrices
```

### State Prediction Matrix Algebra

At step k:

```math
\begin{split}
y_k &= C x_k + D u_k\\
x_k &= A x_{k-1} + B u_{k-1}\\
x_{k+1} &= A x_{k} + B u_{k}
\end{split}
```

Therefore to h steps in history:

For all output estimates:

```math
\begin{split}
y_k & = C x_k + D u_k \\
    & = C (A x_{k-1} + B u_{k-1}) + D u_k \\
    & = C (A (A x_{k-2} + B u_{k-2}) + B u_{k-1}) + D u_k \\
    & = ... \\
    & = C A^{h-1} x_{k-(h-1)} + C A^{h-2} B u_{k-(h-1)} + C A^{h-3} B u_{k-(h-2)} + ... + CAB u_{k-2} + CB u_{k-1} + D u_k \\
\end{split}
```

```math
\begin{split}
y_k         &= C A^{h-1} x_{k-(h-1)} + C A^{h-2} B u_{k-(h-1)} + C A^{h-3} B u_{k-(h-2)} + ... + CAB u_{k-2} + CB u_{k-1} + D u_k \\
y_{k-1}     &= C A^{h-2} x_{k-(h-1)} + C A^{h-3} B u_{k-(h-1)} + C A^{h-4} B u_{k-(h-2)} + ... + CAB u_{k-3} + CB u_{k-2} + D u_{k-1} \\

...         &= ... \\
y_{k-i}     &= C A^{h-1-i} x_{k-(h-1)} + C A^{h-2-i} B u_{k-(h-1)} + C A^{h-3-i} B u_{k-(h-2)} + ... + CB u_{k-i-1} + D u_{k-i} \\
...         &= ... \\
y_{k-(h-3)} &= C A^2 x_{k-(h-1)} + CAB u_{k-(h-1)}+ CB u_{k-(h-2)} + D u_{k-(h-3)} \\
y_{k-(h-2)} &= C A x_{k-(h-1)} + CB u_{k-(h-1)} + D u_{k-(h-2)} \\
y_{k-(h-1)} &= C x_{k-(h-1)} + D u_{k-(h-1)}
\end{split}
```

For predicted next state:

```math
\begin{split}
x_{k+1} & = A x_{k} + B u_{k} \\
        & = A (A x_{k-1} + B u_{k-1}) + B u_{k} \\
        & = A (A (A x_{k-2} + B u_{k-2})  + B u_{k-1}) + B u_{k} \\
        & = ... \\
        & = A^h x_{k-(h-1)} + A^{h-1} B u_{k-(h-1)} + A^{h-2} B u_{k-(h-2)} + ... + A^2 B u_{k-2} + A B u_{k-1} + B u_{k} \\
\end{split}
```

In terms of mega-matrices, define mega-vectors of inputs and outputs:

```math
\vec{Y}  = \begin{bmatrix}
                y_{k-(h-1)} \\
                y_{k-(h-2)} \\
                ... \\
                y_{k-1} \\
                y_{k}
              \end{bmatrix}
```

Note that for multiple outputs, each output vector will be stacked vertically to create a single column.

```math
\vec{U}  = \begin{bmatrix}
                u_{k-(h-1)} \\
                u_{k-(h-2)} \\
                ... \\
                u_{k-1} \\
                u_{k}
              \end{bmatrix}
```

Note that for multiple inputs, each input vector will be stacked vertically to create a single column.

Then, from $x_{k-(h-1)}$ and $\vec{U}$, we get $\vec{Y}$ and predicted state $x_{k+1}$:
```math
\vec{Y} = \mathcal{L} x_{k-(h-1)} + \mathcal{M} \vec{U}
```
```math
x_{k+1} = \mathcal{N} x_{k-(h-1)} + \mathcal{O} \vec{U}
```

Where the mega-matrices $\mathcal{L}$, $\mathcal{M}$, $\mathcal{N}$, $\mathcal{O}$ are:

``\mathcal{L}`` is a matrix with (h x no. of outputs) rows and state-space order columns:

```math
\mathcal{L} = \begin{bmatrix}
                  C... \\
                  CA... \\
                  CA^2... \\
                  ... \\
                  CA^{h-2}... \\
                  CA^{h-1}...
               \end{bmatrix}
```

``\mathcal{M}`` is a matrix with (h x no. of outputs) rows and (h x no. of inputs) columns:

```math
\mathcal{M} = \begin{bmatrix}
                  D           & 0           & 0           & ... & 0     & 0     & 0   & 0   &     \\
                  CB          & D           & 0           & ... & 0     & 0     & 0   & 0   &     \\
                  CAB         & CB          & D           & ... & 0     & 0     & 0   & 0   &     \\
                  ...         & ...         & ...         & ... & ...   & ...   & ... & ... & ... \\
                  C A^{h-4} B & C A^{h-5} B & C A^{h-6} B & ... & CAB   & CB    & D   & 0   & 0   \\
                  C A^{h-3} B & C A^{h-4} B & C A^{h-5} B & ... & CA^2B & CAB   & CB  & D   & 0   \\
                  C A^{h-2} B & C A^{h-3} B & C A^{h-4} B & ... & CA^3B & CA^2B & CAB & CB  & D   \\
                 \end{bmatrix}
```

``\mathcal{N}`` is a square matrix with state-space order columns and rows:

```math
\mathcal{N} = A^{h}
```

``\mathcal{O}`` is a matrix with state-space order rows and (h x no. of inputs) columns:

```math
\mathcal{O} = \begin{bmatrix} A^{h-1} B & A^{h-2} B & A^{h-3} B & ... & A^2 B & AB & B \end{bmatrix}
```

Then, future state can be predicted by:

```math
x_{k+1} = \mathcal{N} \mathcal{L}^{-1} (\vec{Y} - \mathcal{M} \vec{U}) + \mathcal{O} \vec{U}
```

## Unit conversion utilities

```@docs
gas_unit_converter
```