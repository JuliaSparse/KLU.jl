# KLU.jl

This package wraps the [SuiteSparse KLU](https://github.com/DrTimothyAldenDavis/SuiteSparse) solver. A description from the publication:

KLU is a set of routines for solving sparse linear systems of equations. It is particularly well-suited to matrices arising in SPICE-like circuit simulation applications. It relies on a permutation to block triangular form (BTF), several methods for finding a fill-reducing ordering (variants of approximate minimum degree, and nested dissection), and a sparse left-looking LU factorization method to factorize each block.


## Installation
```julia
using Pkg
Pkg.add("KLU")
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

## Citing this package

Please cite both [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/dev/CITATION.bib) as well as [KLU](https://github.com/JuliaSparse/KLU.jl/blob/main/CITATION.bib):

## Exported API

```@docs
klu
```
