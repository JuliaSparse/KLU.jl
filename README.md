# KLU
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://klu.juliasparse.org/)

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

## Supported types

KLU provides `Float64` and `ComplexF64` values with `Int32` and `Int64` indices, and nothing else is
accepted. Right-hand sides must have the same element type as the factorization. Convert anything
else first, for example with `klu(Float64.(A))` or `factor \ Float64.(b)`.

## Reusing a factorization

If only the values of a matrix change, `klu!` refactors faster than `klu`. It reuses the
pivot order, and throws an `ArgumentError` if the pattern is different.

```julia
julia> klu!(factor, B);          # B has the same pattern as A

julia> klu!(factor, nonzeros(B)); # or just its values
```

If `klu!` throws a `SingularException`, `factor` has no numeric factorization until a later
`klu!` or `klu_factor!` succeeds, and solving with it throws as well.

To avoid using more memory, `factor` shares the vectors it is given: after `klu(A)` or
`klu!(factor, B)`, `nonzeros(factor)` is the vector of `A` or `B`. It is the caller's
responsibility not to change them in a way that breaks `factor`, as described in the docstring of
`klu`.

`klu(A; full_factor = false)` only analyzes the pattern of `A`, and `klu!` then computes the
first numeric factorization. Use `copy(factor)` to get an independent factorization, for
example to refactor it without affecting the original. This costs as much as a new
factorization.

## Checking a factorization

`rcond` is a cheap estimate of the reciprocal condition number, `condest` is a more accurate
estimate of the 1-norm condition number, and `rgrowth` is the reciprocal pivot growth. `det`,
`logabsdet`, and `logdet` give the determinant of the factored matrix.

## Thread safety

KLU's C library keeps scratch workspace *inside* each numeric factorization, so two
simultaneous solves (or a solve and a refactorization) on the **same** `KLUFactorization`
would corrupt each other. KLU.jl therefore guards every libklu call with a lock stored
in the factorization object: sharing one factorization between threads is safe, but those
calls are serialized. Distinct factorizations share no state and run fully in parallel.

To solve many right-hand sides against one matrix, batch them as the columns of a single
`B` in one `solve!`/`ldiv!` call, which runs on one thread. To solve them in parallel, give
each task its own `copy(F)` or `klu(A)`, each of which costs as much time and memory as a
new factorization.

To make a sequence of calls on a factorization `F` atomic, hold its reentrant lock with
`Base.@lock F ...`.

<!---[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://USER_NAME.github.io/PACKAGE_NAME.jl/stable) --->

## Citing this package

Please cite both [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/dev/CITATION.bib) as well as [KLU](https://github.com/JuliaSparse/KLU.jl/blob/main/CITATION.bib).
