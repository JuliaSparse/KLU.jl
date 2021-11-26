# KLU
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://wimmerer.github.io/KLU.jl/dev/)

A linear solver for sparse circuit simulation problems.

## Installation
```julia
using Pkg
Pkg.add(url="https://github.com/Wimmerer/KLU.jl")
```

## Usage

```julia
julia> using KLU
julia> using SparseArrays
julia> A = sprand(15, 15, 0.45);
julia> factor = klu(A);
julia> b = rand(15);
julia> x = factor \ b
15-element Vector{Float64}:
  22.501711574339556
  27.80357152530039
 -37.08488119597194
  61.470604418791346
  -0.49109038415655193
  38.03163985997178
   1.2417021657549396
 -51.154984106483155
   0.9488554483842341
   1.587572050319512
 -31.251127382496147
   1.4508846663863602
 -25.866070230462796
  -5.46661465074706
  10.688731342850245
```
<!---[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://USER_NAME.github.io/PACKAGE_NAME.jl/stable) --->


