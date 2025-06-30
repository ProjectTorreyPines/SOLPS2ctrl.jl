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

## For developers

If you are contributing to this project, you would need to install [dvc](https://dvc.org/) to fetch sample files for testing. Once installed, please configure your ssh so that you can ssh into omega tunneling through cybele without requiring to enter password. This is optional but will make it much easier for you.

Once you have completed above steps, inside the git repo, simply do:
```bash
dvc pull
```

This would download the sample files in the `samples` directory. Then to run tests, you would first need to instantiate the project:
```bash
julia --project
```
Then press `]`:
```julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.11.5 (2025-04-14)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> ]
```
Then type:
```julia
(SOLPS2ctrl) pkg> instantiate
```
Once the package has been instantiated, you can run the tests using:
```julia
(SOLPS2ctrl) pkg> test
```
This would run all the tests though. To run specific tests, you can do following from the command line to see help options (this works after you ahve instantiated the project like mentioned above):
```bash
% julia --project test/test.jl help
Usage (from inside SOLPS2ctrl.jl): 
julia --project test/test.jl [units] [core] [edge] [heavy] [repair] [geqdsk] [prep] [sysid] [state] [controller] [h] [help]

Run tests. Default is all tests.

Optional arguments:
    units      Test unit conversion utilities  
    core       Test core profile extension     
    edge       Test edge profile extension     
    heavy      Test heavy utilities            
    repair     Test repair_eq                  
    geqdsk     Test geqdsk_to_imas             
    prep       Test preparation                
    sysid      Test system id                  
    state      Test state prediction           
    controller Test linear and PVLC controllers
    h          Show this help message and exit 
    help       Show this help message and exit 

```
