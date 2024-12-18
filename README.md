# SOLPS2ctrl

![Format Check](https://github.com/ProjectTorreyPines/SOLPS2ctrl.jl/actions/workflows/format_check.yml/badge.svg)
![Docs](https://github.com/ProjectTorreyPines/SOLPS2ctrl.jl/actions/workflows/make_docs.yml/badge.svg)
![Tests](https://github.com/ProjectTorreyPines/SOLPS2ctrl.jl/actions/workflows/test.yml/badge.svg)
[![codecov](https://codecov.io/gh/ProjectTorreyPines/SOLPS2ctrl.jl/graph/badge.svg?token=ZJBRLAXIS1)](https://codecov.io/gh/ProjectTorreyPines/SOLPS2ctrl.jl)

This repository is the top level layer for managing workflow for:

* Loading SOLPS outputs into IMAS data dictionary format
* Loading equilibrium (that the SOLPS mesh was based on) into IMAS DD format
* Extrapolating results into core and far SOL region if required
* Running synthetic diagnostics on them
* Performing system identification (to be added)
* Designing and tuning linear causal and model predictive controllers (to be added)

For installation and usage instructions, see the [online documentation](https://projecttorreypines.github.io/SOLPS2ctrl.jl/stable). For documentation on under development branch, see [dev online documentation](https://projecttorreypines.github.io/SOLPS2ctrl.jl/dev).

## Installation

```
using Pkg
Pkg.add("SOLPS2ctrl")
```

## Examples

Refer to the instructions on this [wiki page](https://github.com/ProjectTorreyPines/SOLPS2ctrl.jl/wiki/Demo) to see how to run `examples/demo.ipynb`.
