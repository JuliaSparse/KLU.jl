__free_sym(::Type{Int32}, args...) = KLU.klu_free_symbolic(args...)
__free_sym(::Type{Int64}, args...) = KLU.klu_l_free_symbolic(args...)

__free_num(::Type{Float64}, ::Type{Int32}, args...) = KLU.klu_free_numeric(args...)
__free_num(::Type{Float64}, ::Type{Int64}, args...) = KLU.klu_l_free_numeric(args...)
__free_num(::Type{ComplexF64}, ::Type{Int32}, args...) = KLU.klu_z_free_numeric(args...)
__free_num(::Type{ComplexF64}, ::Type{Int64}, args...) = KLU.klu_zl_free_numeric(args...)

__analyze(::Type{Int32}, args...) = KLU.klu_analyze(args...)
__analyze(::Type{Int64}, args...) = KLU.klu_l_analyze(args...)

__analyze!(::Type{Int32}, args...) = KLU.klu_analyze_given(args...)
__analyze!(::Type{Int64}, args...) = KLU.klu_l_analyze_given(args...)

__factor(::Type{Float64}, ::Type{Int32}, args...) = KLU.klu_factor(args...)
__factor(::Type{Float64}, ::Type{Int64}, args...) = KLU.klu_l_factor(args...)
__factor(::Type{ComplexF64}, ::Type{Int32}, args...) = KLU.klu_z_factor(args...)
__factor(::Type{ComplexF64}, ::Type{Int64}, args...) = KLU.klu_zl_factor(args...)

__rcond(::Type{Float64}, ::Type{Int32}, args...) = KLU.klu_rcond(args...)
__rcond(::Type{Float64}, ::Type{Int64}, args...) = KLU.klu_l_rcond(args...)
__rcond(::Type{ComplexF64}, ::Type{Int32}, args...) = KLU.klu_z_rcond(args...)
__rcond(::Type{ComplexF64}, ::Type{Int64}, args...) = KLU.klu_zl_rcond(args...)

__rgrowth(::Type{Float64}, ::Type{Int32}, args...) = KLU.klu_rgrowth(args...)
__rgrowth(::Type{Float64}, ::Type{Int64}, args...) = KLU.klu_l_rgrowth(args...)
__rgrowth(::Type{ComplexF64}, ::Type{Int32}, args...) = KLU.klu_z_rgrowth(args...)
__rgrowth(::Type{ComplexF64}, ::Type{Int64}, args...) = KLU.klu_zl_rgrowth(args...)

__condest(::Type{Float64}, ::Type{Int32}, args...) = KLU.klu_condest(args...)
__condest(::Type{Float64}, ::Type{Int64}, args...) = KLU.klu_l_condest(args...)
__condest(::Type{ComplexF64}, ::Type{Int32}, args...) = KLU.klu_z_condest(args...)
__condest(::Type{ComplexF64}, ::Type{Int64}, args...) = KLU.klu_zl_condest(args...)

__refactor(::Type{Float64}, ::Type{Int32}, args...) = KLU.klu_refactor(args...)
__refactor(::Type{Float64}, ::Type{Int64}, args...) = KLU.klu_l_refactor(args...)
__refactor(::Type{ComplexF64}, ::Type{Int32}, args...) = KLU.klu_z_refactor(args...)
__refactor(::Type{ComplexF64}, ::Type{Int64}, args...) = KLU.klu_zl_refactor(args...)

__solve(::Type{Float64}, ::Type{Int32}, args...) = KLU.klu_solve(args...)
__solve(::Type{Float64}, ::Type{Int64}, args...) = KLU.klu_l_solve(args...)
__solve(::Type{ComplexF64}, ::Type{Int32}, args...) = KLU.klu_z_solve(args...)
__solve(::Type{ComplexF64}, ::Type{Int64}, args...) = KLU.klu_zl_solve(args...)

__symType(::Type{Int32}) = KLU.klu_symbolic
__symType(::Type{Int64}) = KLU.klu_l_symbolic

__numType(::Type{Int32}) = KLU.klu_numeric
__numType(::Type{Int64}) = KLU.klu_l_numeric

# Could rewrite with multiple dispatch, but more awkward because of the 
# varying number of arguments between real vs complex cases: tsolve, extract.
# sort is only used in extract, so also not implemented here.