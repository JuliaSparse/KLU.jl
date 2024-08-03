var documenterSearchIndex = {"docs":
[{"location":"functions_list/#Functions","page":"Reference","title":"Functions","text":"","category":"section"},{"location":"functions_list/","page":"Reference","title":"Reference","text":"CurrentModule = KLU","category":"page"},{"location":"functions_list/","page":"Reference","title":"Reference","text":"Modules = [KLU]\nOrder   = [:module, :type, :function]","category":"page"},{"location":"functions_list/#KLU.KLUFactorization","page":"Reference","title":"KLU.KLUFactorization","text":"KLUFactorization <: Factorization\n\nMatrix factorization type of the KLU factorization of a sparse matrix A. This is the return type of klu, the corresponding matrix factorization function.\n\nThe factors can be obtained from K::KLUFactorization via K.L, K.U and K.F See the klu docs for more information.\n\nYou typically should not construct this directly, instead use klu.\n\n\n\n\n\n","category":"type"},{"location":"functions_list/#KLU.klu_common","page":"Reference","title":"KLU.klu_common","text":"Data structure for parameters of and statistics generated by KLU functions.\n\nThis is the Int32 version of klu_l_common.\n\n\n\n\n\n","category":"type"},{"location":"functions_list/#KLU.klu_l_common","page":"Reference","title":"KLU.klu_l_common","text":"Data structure for parameters of and statistics generated by KLU functions.\n\nFields\n\ntol::Float64: Partial pivoting tolerance for diagonal preference\nbtf::Int64: If btf != 0 use BTF pre-ordering\nordering::Int64: If ordering == 0 use AMD to permute, if ordering == 1 use COLAMD,\n\nif ordering == 3 use the user provided ordering function.\n\nscale::Int64: If scale == 1 then A[:,i] ./= sum(abs.(A[:,i])), if scale == 2 then\n\nA[:,i] ./= maximum(abs.(A[:,i])). If scale == 0 no scaling is done, and the input is checked for errors if scale >= 0.\n\nSee the KLU User Guide for more information.\n\n\n\n\n\n","category":"type"},{"location":"functions_list/#KLU.condest-Tuple{KLU.KLUFactorization{ComplexF64, Int32}}","page":"Reference","title":"KLU.condest","text":"condest(K::KLUFactorization)\n\nAccurately estimate the 1-norm condition number of the factorization.\n\n\n\n\n\n","category":"method"},{"location":"functions_list/#KLU.condest-Tuple{KLU.KLUFactorization{ComplexF64, Int64}}","page":"Reference","title":"KLU.condest","text":"condest(K::KLUFactorization)\n\nAccurately estimate the 1-norm condition number of the factorization.\n\n\n\n\n\n","category":"method"},{"location":"functions_list/#KLU.condest-Tuple{KLU.KLUFactorization{Float64, Int32}}","page":"Reference","title":"KLU.condest","text":"condest(K::KLUFactorization)\n\nAccurately estimate the 1-norm condition number of the factorization.\n\n\n\n\n\n","category":"method"},{"location":"functions_list/#KLU.condest-Tuple{KLU.KLUFactorization{Float64, Int64}}","page":"Reference","title":"KLU.condest","text":"condest(K::KLUFactorization)\n\nAccurately estimate the 1-norm condition number of the factorization.\n\n\n\n\n\n","category":"method"},{"location":"functions_list/#KLU.klu!","page":"Reference","title":"KLU.klu!","text":"klu!(K::KLUFactorization, A::SparseMatrixCSC; check=true, allowsingular=false) -> K::KLUFactorization\nklu(K::KLUFactorization, nzval::Vector{Tv}; check=true, allowsingular=false) -> K::KLUFactorization\n\nRecompute the KLU factorization K using the values of A or nzval. The pattern of the original matrix used to create K must match the pattern of A. \n\nFor sparse A with real or complex element type, the return type of K is KLUFactorization{Tv, Ti}, with Tv = Float64 or ComplexF64 respectively and Ti is an integer type (Int32 or Int64).\n\nSee also: klu\n\nArguments\n\nK::KLUFactorization: The matrix factorization object, previously created by a call to klu, to be re-factored.\nA::SparseMatrixCSC or n::Integer, colptr::Vector{Ti}, rowval::Vector{Ti}, nzval::Vector{Tv}: The sparse matrix or the zero-based sparse matrix components to be factored.\ncheck::Bool: If true (default) check for errors after the factorization. If false errors must be checked by the user with klu.common.status.\nallowsingular::Bool: If true (default false) allow the factorization to proceed even if the matrix is singular. Note that this will allow for\n\nsilent divide by zero errors in subsequent solve! or ldiv! calls if singularity is not checked by the user with klu.common.status == KLU.KLU_SINGULAR\n\nnote: Note\nklu(A::SparseMatrixCSC) uses the KLU[ACM907] library that is part of SuiteSparse. As this library only supports sparse matrices with Float64 or ComplexF64 elements, lu converts A into a copy that is of type SparseMatrixCSC{Float64} or SparseMatrixCSC{ComplexF64} as appropriate.\n\n[ACM907]: Davis, Timothy A., & Palamadai Natarajan, E. (2010). Algorithm 907: KLU, A Direct Sparse Solver for Circuit Simulation Problems. ACM Trans. Math. Softw., 37(3). doi:10.1145/1824801.1824814\n\n\n\n\n\n","category":"function"},{"location":"functions_list/#KLU.klu-Union{Tuple{Tv}, Tuple{Ti}, Tuple{Any, Vector{Ti}, Vector{Ti}, Vector{Tv}}} where {Ti<:Union{Int32, Int64}, Tv<:AbstractFloat}","page":"Reference","title":"KLU.klu","text":"klu(A::SparseMatrixCSC) -> K::KLUFactorization\nklu(n, colptr::Vector{Ti}, rowval::Vector{Ti}, nzval::Vector{Tv}) -> K::KLUFactorization\n\nCompute the LU factorization of a sparse matrix A using KLU[ACM907].\n\nFor sparse A with real or complex element type, the return type of K is KLUFactorization{Tv, Ti}, with Tv = Float64 or ComplexF64 respectively and Ti is an integer type (Int32 or Int64).\n\nThe individual components of the factorization K can be accessed by indexing:\n\nComponent Description\nL L (lower triangular) part of LU of the diagonal blocks\nU U (upper triangular) part of LU of the diagonal blocks\nF F (upper triangular) part of LU + F, the off-diagonal blocks\np right permutation Vector\nq left permutation Vector\nRs Vector of scaling factors\n\nThe relation between K and A is\n\nK.L * K.U + K.F  == K.Rs \\ A[K.p, K.q]\n\nK further supports the following functions:\n\nLinearAlgebra.\\\n\nArguments\n\nA::SparseMatrixCSC or n::Integer, colptr::Vector{Ti}, rowval::Vector{Ti}, nzval::Vector{Tv}: The sparse matrix or the zero-based sparse matrix components to be factored.\ncheck::Bool: If true (default) check for errors after the factorization. If false errors must be checked by the user with klu.common.status.\nallowsingular::Bool: If true (default false) allow the factorization to proceed even if the matrix is singular. Note that this will allow for\n\nsilent divide by zero errors in subsequent solve! or ldiv! calls if singularity is not checked by the user with klu.common.status == KLU.KLU_SINGULAR\n\nnote: Note\nklu(A::SparseMatrixCSC) uses the KLU[ACM907] library that is part of SuiteSparse. As this library only supports sparse matrices with Float64 or ComplexF64 elements, lu converts A into a copy that is of type SparseMatrixCSC{Float64} or SparseMatrixCSC{ComplexF64} as appropriate.\n\n[ACM907]: Davis, Timothy A., & Palamadai Natarajan, E. (2010). Algorithm 907: KLU, A Direct Sparse Solver for Circuit Simulation Problems. ACM Trans. Math. Softw., 37(3). doi:10.1145/1824801.1824814\n\n\n\n\n\n","category":"method"},{"location":"functions_list/#KLU.klu_factor!","page":"Reference","title":"KLU.klu_factor!","text":"klu_factor!(K::KLUFactorization; check=true, allowsingular=false) -> K::KLUFactorization\n\nFactor K into components K.L, K.U, and K.F. This function will perform both the symbolic and numeric steps of factoriation on an existing KLUFactorization instance.\n\nArguments\n\nK::KLUFactorization: The matrix factorization object to be factored.\ncheck::Bool: If true (default) check for errors after the factorization. If false errors must be checked by the user with klu.common.status.\nallowsingular::Bool: If true (default false) allow the factorization to proceed even if the matrix is singular. Note that this will allow for\n\nsilent divide by zero errors in subsequent solve! or ldiv! calls if singularity is not checked by the user with klu.common.status == KLU.KLU_SINGULAR\n\nThe K.common struct can be used to modify certain options and parameters, see the KLU documentation or klu_common for more information.\n\n\n\n\n\n","category":"function"},{"location":"functions_list/#KLU.rcond-Tuple{KLU.AbstractKLUFactorization{ComplexF64, Int32}}","page":"Reference","title":"KLU.rcond","text":"rcond(K::KLUFactorization)\n\nCheaply estimate the reciprocal condition number.\n\n\n\n\n\n","category":"method"},{"location":"functions_list/#KLU.rcond-Tuple{KLU.AbstractKLUFactorization{ComplexF64, Int64}}","page":"Reference","title":"KLU.rcond","text":"rcond(K::KLUFactorization)\n\nCheaply estimate the reciprocal condition number.\n\n\n\n\n\n","category":"method"},{"location":"functions_list/#KLU.rcond-Tuple{KLU.AbstractKLUFactorization{Float64, Int32}}","page":"Reference","title":"KLU.rcond","text":"rcond(K::KLUFactorization)\n\nCheaply estimate the reciprocal condition number.\n\n\n\n\n\n","category":"method"},{"location":"functions_list/#KLU.rcond-Tuple{KLU.AbstractKLUFactorization{Float64, Int64}}","page":"Reference","title":"KLU.rcond","text":"rcond(K::KLUFactorization)\n\nCheaply estimate the reciprocal condition number.\n\n\n\n\n\n","category":"method"},{"location":"functions_list/#KLU.rgrowth-Tuple{KLU.KLUFactorization{ComplexF64, Int32}}","page":"Reference","title":"KLU.rgrowth","text":"rgrowth(K::KLUFactorization)\n\nCalculate the reciprocal pivot growth.\n\n\n\n\n\n","category":"method"},{"location":"functions_list/#KLU.rgrowth-Tuple{KLU.KLUFactorization{ComplexF64, Int64}}","page":"Reference","title":"KLU.rgrowth","text":"rgrowth(K::KLUFactorization)\n\nCalculate the reciprocal pivot growth.\n\n\n\n\n\n","category":"method"},{"location":"functions_list/#KLU.rgrowth-Tuple{KLU.KLUFactorization{Float64, Int32}}","page":"Reference","title":"KLU.rgrowth","text":"rgrowth(K::KLUFactorization)\n\nCalculate the reciprocal pivot growth.\n\n\n\n\n\n","category":"method"},{"location":"functions_list/#KLU.rgrowth-Tuple{KLU.KLUFactorization{Float64, Int64}}","page":"Reference","title":"KLU.rgrowth","text":"rgrowth(K::KLUFactorization)\n\nCalculate the reciprocal pivot growth.\n\n\n\n\n\n","category":"method"},{"location":"functions_list/#KLU.solve!","page":"Reference","title":"KLU.solve!","text":"solve!(klu::Union{KLUFactorization, AdjointFact{Tv, KLUFactorization}, TransposeFact{Tv, KLUFactorization}}, B::StridedVecOrMat) -> B\n\nSolve the linear system AX = B or A'X = B where using the matrix facotrization klu and right-hand side B.\n\nThis function overwrites B with the solution X, for a new solution vector X use solve(klu, B).\n\nArguments\n\nklu::KLUFactorization: The matrix factorization of A to use in the solution.\nB::StridedVecOrMat: The right-hand side of the linear system.\ncheck::Bool: If true (default) check for errors after the solve. If false errors must be checked using klu.common.status. The return value of this function is always B, so the status of the solve must be checked using the factorization object itself when check=false.\n\nwarning: Solve with a singular factorization object\nIf the factorization object klu has klu.common.status == KLU.KLU_SINGULAR then the solve! or ldiv! will result in a silent divide by zero error.This status should be checked by the user before solve calls if singularity checks were disabled on factorization using check=false or allowsingular=true.\n\n\n\n\n\n","category":"function"},{"location":"#KLU.jl","page":"KLU","title":"KLU.jl","text":"","category":"section"},{"location":"","page":"KLU","title":"KLU","text":"This package wraps the SuiteSparse KLU solver. A description from the publication:","category":"page"},{"location":"","page":"KLU","title":"KLU","text":"KLU is a set of routines for solving sparse linear systems of equations. It is particularly well-suited to matrices arising in SPICE-like circuit simulation applications. It relies on a permutation to block triangular form (BTF), several methods for finding a fill-reducing ordering (variants of approximate minimum degree, and nested dissection), and a sparse left-looking LU factorization method to factorize each block.","category":"page"},{"location":"#Installation","page":"KLU","title":"Installation","text":"","category":"section"},{"location":"","page":"KLU","title":"KLU","text":"using Pkg\nPkg.add(\"KLU\")","category":"page"},{"location":"#Usage","page":"KLU","title":"Usage","text":"","category":"section"},{"location":"","page":"KLU","title":"KLU","text":"julia> using KLU\njulia> using SparseArrays\njulia> A = sprand(15, 15, 0.45);\njulia> factor = klu(A);\njulia> b = rand(15);\njulia> x = factor \\ b\n15-element Vector{Float64}:\n  22.501711574339556\n  27.80357152530039\n -37.08488119597194\n  61.470604418791346\n  -0.49109038415655193\n  38.03163985997178\n   1.2417021657549396\n -51.154984106483155\n   0.9488554483842341\n   1.587572050319512\n -31.251127382496147\n   1.4508846663863602\n -25.866070230462796\n  -5.46661465074706\n  10.688731342850245","category":"page"},{"location":"#Citing-this-package","page":"KLU","title":"Citing this package","text":"","category":"section"},{"location":"","page":"KLU","title":"KLU","text":"Please cite both SuiteSparse as well as KLU:","category":"page"},{"location":"#Exported-API","page":"KLU","title":"Exported API","text":"","category":"section"},{"location":"","page":"KLU","title":"KLU","text":"klu","category":"page"},{"location":"#KLU.klu","page":"KLU","title":"KLU.klu","text":"klu(A::SparseMatrixCSC) -> K::KLUFactorization\nklu(n, colptr::Vector{Ti}, rowval::Vector{Ti}, nzval::Vector{Tv}) -> K::KLUFactorization\n\nCompute the LU factorization of a sparse matrix A using KLU[ACM907].\n\nFor sparse A with real or complex element type, the return type of K is KLUFactorization{Tv, Ti}, with Tv = Float64 or ComplexF64 respectively and Ti is an integer type (Int32 or Int64).\n\nThe individual components of the factorization K can be accessed by indexing:\n\nComponent Description\nL L (lower triangular) part of LU of the diagonal blocks\nU U (upper triangular) part of LU of the diagonal blocks\nF F (upper triangular) part of LU + F, the off-diagonal blocks\np right permutation Vector\nq left permutation Vector\nRs Vector of scaling factors\n\nThe relation between K and A is\n\nK.L * K.U + K.F  == K.Rs \\ A[K.p, K.q]\n\nK further supports the following functions:\n\nLinearAlgebra.\\\n\nArguments\n\nA::SparseMatrixCSC or n::Integer, colptr::Vector{Ti}, rowval::Vector{Ti}, nzval::Vector{Tv}: The sparse matrix or the zero-based sparse matrix components to be factored.\ncheck::Bool: If true (default) check for errors after the factorization. If false errors must be checked by the user with klu.common.status.\nallowsingular::Bool: If true (default false) allow the factorization to proceed even if the matrix is singular. Note that this will allow for\n\nsilent divide by zero errors in subsequent solve! or ldiv! calls if singularity is not checked by the user with klu.common.status == KLU.KLU_SINGULAR\n\nnote: Note\nklu(A::SparseMatrixCSC) uses the KLU[ACM907] library that is part of SuiteSparse. As this library only supports sparse matrices with Float64 or ComplexF64 elements, lu converts A into a copy that is of type SparseMatrixCSC{Float64} or SparseMatrixCSC{ComplexF64} as appropriate.\n\n[ACM907]: Davis, Timothy A., & Palamadai Natarajan, E. (2010). Algorithm 907: KLU, A Direct Sparse Solver for Circuit Simulation Problems. ACM Trans. Math. Softw., 37(3). doi:10.1145/1824801.1824814\n\n\n\n\n\n","category":"function"}]
}