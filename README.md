# KLU
[![](https://img.shields.io/badge/docs-dev-blue.svg)](http://klu.juliasparse.org/dev/)

A linear solver for sparse circuit simulation problems.
This package wraps the [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) solver KLU.

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
<!---[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://USER_NAME.github.io/PACKAGE_NAME.jl/stable) --->

## Citing this package

Please cite both [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) as well as KLU:

```bibtex
@article{10.1145/1824801.1824814,
author = {Davis, Timothy A. and Palamadai Natarajan, Ekanathan},
title = {Algorithm 907: KLU, A Direct Sparse Solver for Circuit Simulation Problems},
year = {2010},
issue_date = {September 2010},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
volume = {37},
number = {3},
issn = {0098-3500},
url = {https://doi.org/10.1145/1824801.1824814},
doi = {10.1145/1824801.1824814},
journal = {ACM Trans. Math. Softw.},
month = {sep},
articleno = {36},
numpages = {17},
keywords = {LU factorization, circuit simulation, sparse matrices}
}
```