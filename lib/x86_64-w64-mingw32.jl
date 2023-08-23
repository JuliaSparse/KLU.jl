mutable struct klu_symbolic
    symmetry::Cdouble
    est_flops::Cdouble
    lnz::Cdouble
    unz::Cdouble
    Lnz::Ptr{Cdouble}
    n::Int32
    nz::Int32
    P::Ptr{Int32}
    Q::Ptr{Int32}
    R::Ptr{Int32}
    nzoff::Int32
    nblocks::Int32
    maxblock::Int32
    ordering::Int32
    do_btf::Int32
    structural_rank::Int32
    klu_symbolic() = new()
end

mutable struct klu_l_symbolic
    symmetry::Cdouble
    est_flops::Cdouble
    lnz::Cdouble
    unz::Cdouble
    Lnz::Ptr{Cdouble}
    n::Int64
    nz::Int64
    P::Ptr{Int64}
    Q::Ptr{Int64}
    R::Ptr{Int64}
    nzoff::Int64
    nblocks::Int64
    maxblock::Int64
    ordering::Int64
    do_btf::Int64
    structural_rank::Int64
    klu_l_symbolic() = new()
end

mutable struct klu_numeric
    n::Int32
    nblocks::Int32
    lnz::Int32
    unz::Int32
    max_lnz_block::Int32
    max_unz_block::Int32
    Pnum::Ptr{Int32}
    Pinv::Ptr{Int32}
    Lip::Ptr{Int32}
    Uip::Ptr{Int32}
    Llen::Ptr{Int32}
    Ulen::Ptr{Int32}
    LUbx::Ptr{Ptr{Cvoid}}
    LUsize::Ptr{Csize_t}
    Udiag::Ptr{Cvoid}
    Rs::Ptr{Cdouble}
    worksize::Csize_t
    Work::Ptr{Cvoid}
    Xwork::Ptr{Cvoid}
    Iwork::Ptr{Int32}
    Offp::Ptr{Int32}
    Offi::Ptr{Int32}
    Offx::Ptr{Cvoid}
    nzoff::Int32
    klu_numeric() = new()
end

mutable struct klu_l_numeric
    n::Int64
    nblocks::Int64
    lnz::Int64
    unz::Int64
    max_lnz_block::Int64
    max_unz_block::Int64
    Pnum::Ptr{Int64}
    Pinv::Ptr{Int64}
    Lip::Ptr{Int64}
    Uip::Ptr{Int64}
    Llen::Ptr{Int64}
    Ulen::Ptr{Int64}
    LUbx::Ptr{Ptr{Cvoid}}
    LUsize::Ptr{Csize_t}
    Udiag::Ptr{Cvoid}
    Rs::Ptr{Cdouble}
    worksize::Csize_t
    Work::Ptr{Cvoid}
    Xwork::Ptr{Cvoid}
    Iwork::Ptr{Int64}
    Offp::Ptr{Int64}
    Offi::Ptr{Int64}
    Offx::Ptr{Cvoid}
    nzoff::Int64
    klu_l_numeric() = new()
end

mutable struct klu_common_struct
    tol::Cdouble
    memgrow::Cdouble
    initmem_amd::Cdouble
    initmem::Cdouble
    maxwork::Cdouble
    btf::Cint
    ordering::Cint
    scale::Cint
    user_order::Ptr{Cvoid}
    user_data::Ptr{Cvoid}
    halt_if_singular::Cint
    status::Cint
    nrealloc::Cint
    structural_rank::Int32
    numerical_rank::Int32
    singular_col::Int32
    noffdiag::Int32
    flops::Cdouble
    rcond::Cdouble
    condest::Cdouble
    rgrowth::Cdouble
    work::Cdouble
    memusage::Csize_t
    mempeak::Csize_t
    klu_common_struct() = new()
end

const klu_common = klu_common_struct

mutable struct klu_l_common_struct
    tol::Cdouble
    memgrow::Cdouble
    initmem_amd::Cdouble
    initmem::Cdouble
    maxwork::Cdouble
    btf::Cint
    ordering::Cint
    scale::Cint
    user_order::Ptr{Cvoid}
    user_data::Ptr{Cvoid}
    halt_if_singular::Cint
    status::Cint
    nrealloc::Cint
    structural_rank::Int64
    numerical_rank::Int64
    singular_col::Int64
    noffdiag::Int64
    flops::Cdouble
    rcond::Cdouble
    condest::Cdouble
    rgrowth::Cdouble
    work::Cdouble
    memusage::Csize_t
    mempeak::Csize_t
    klu_l_common_struct() = new()
end

const klu_l_common = klu_l_common_struct

function klu_defaults(Common)
    @ccall libklu.klu_defaults(Common::Ptr{klu_common})::Cint
end

function klu_l_defaults(Common)
    @ccall libklu.klu_l_defaults(Common::Ptr{klu_l_common})::Cint
end

function klu_analyze(n, Ap, Ai, Common)
    @ccall libklu.klu_analyze(n::Int32, Ap::Ptr{Int32}, Ai::Ptr{Int32}, Common::Ptr{klu_common})::Ptr{klu_symbolic}
end

function klu_l_analyze(arg1, arg2, arg3, Common)
    @ccall libklu.klu_l_analyze(arg1::Int64, arg2::Ptr{Int64}, arg3::Ptr{Int64}, Common::Ptr{klu_l_common})::Ptr{klu_l_symbolic}
end

function klu_analyze_given(n, Ap, Ai, P, Q, Common)
    @ccall libklu.klu_analyze_given(n::Int32, Ap::Ptr{Int32}, Ai::Ptr{Int32}, P::Ptr{Int32}, Q::Ptr{Int32}, Common::Ptr{klu_common})::Ptr{klu_symbolic}
end

function klu_l_analyze_given(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_l_analyze_given(arg1::Int64, arg2::Ptr{Int64}, arg3::Ptr{Int64}, arg4::Ptr{Int64}, arg5::Ptr{Int64}, arg6::Ptr{klu_l_common})::Ptr{klu_l_symbolic}
end

function klu_factor(Ap, Ai, Ax, Symbolic, Common)
    @ccall libklu.klu_factor(Ap::Ptr{Int32}, Ai::Ptr{Int32}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Common::Ptr{klu_common})::Ptr{klu_numeric}
end

function klu_z_factor(Ap, Ai, Ax, Symbolic, Common)
    @ccall libklu.klu_z_factor(Ap::Ptr{Int32}, Ai::Ptr{Int32}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Common::Ptr{klu_common})::Ptr{klu_numeric}
end

function klu_l_factor(arg1, arg2, arg3, arg4, arg5)
    @ccall libklu.klu_l_factor(arg1::Ptr{Int64}, arg2::Ptr{Int64}, arg3::Ptr{Cdouble}, arg4::Ptr{klu_l_symbolic}, arg5::Ptr{klu_l_common})::Ptr{klu_l_numeric}
end

function klu_zl_factor(arg1, arg2, arg3, arg4, arg5)
    @ccall libklu.klu_zl_factor(arg1::Ptr{Int64}, arg2::Ptr{Int64}, arg3::Ptr{Cdouble}, arg4::Ptr{klu_l_symbolic}, arg5::Ptr{klu_l_common})::Ptr{klu_l_numeric}
end

function klu_solve(Symbolic, Numeric, ldim, nrhs, B, Common)
    @ccall libklu.klu_solve(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, ldim::Int32, nrhs::Int32, B::Ptr{Cdouble}, Common::Ptr{klu_common})::Cint
end

function klu_z_solve(Symbolic, Numeric, ldim, nrhs, B, Common)
    @ccall libklu.klu_z_solve(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, ldim::Int32, nrhs::Int32, B::Ptr{Cdouble}, Common::Ptr{klu_common})::Cint
end

function klu_l_solve(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_l_solve(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Int64, arg4::Int64, arg5::Ptr{Cdouble}, arg6::Ptr{klu_l_common})::Cint
end

function klu_zl_solve(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_zl_solve(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Int64, arg4::Int64, arg5::Ptr{Cdouble}, arg6::Ptr{klu_l_common})::Cint
end

function klu_tsolve(Symbolic, Numeric, ldim, nrhs, B, Common)
    @ccall libklu.klu_tsolve(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, ldim::Int32, nrhs::Int32, B::Ptr{Cdouble}, Common::Ptr{klu_common})::Cint
end

function klu_z_tsolve(Symbolic, Numeric, ldim, nrhs, B, conj_solve, Common)
    @ccall libklu.klu_z_tsolve(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, ldim::Int32, nrhs::Int32, B::Ptr{Cdouble}, conj_solve::Cint, Common::Ptr{klu_common})::Cint
end

function klu_l_tsolve(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_l_tsolve(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Int64, arg4::Int64, arg5::Ptr{Cdouble}, arg6::Ptr{klu_l_common})::Cint
end

function klu_zl_tsolve(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libklu.klu_zl_tsolve(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Int64, arg4::Int64, arg5::Ptr{Cdouble}, arg6::Cint, arg7::Ptr{klu_l_common})::Cint
end

function klu_refactor(Ap, Ai, Ax, Symbolic, Numeric, Common)
    @ccall libklu.klu_refactor(Ap::Ptr{Int32}, Ai::Ptr{Int32}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_z_refactor(Ap, Ai, Ax, Symbolic, Numeric, Common)
    @ccall libklu.klu_z_refactor(Ap::Ptr{Int32}, Ai::Ptr{Int32}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_l_refactor(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_l_refactor(arg1::Ptr{Int64}, arg2::Ptr{Int64}, arg3::Ptr{Cdouble}, arg4::Ptr{klu_l_symbolic}, arg5::Ptr{klu_l_numeric}, arg6::Ptr{klu_l_common})::Cint
end

function klu_zl_refactor(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_zl_refactor(arg1::Ptr{Int64}, arg2::Ptr{Int64}, arg3::Ptr{Cdouble}, arg4::Ptr{klu_l_symbolic}, arg5::Ptr{klu_l_numeric}, arg6::Ptr{klu_l_common})::Cint
end

function klu_free_symbolic(Symbolic, Common)
    @ccall libklu.klu_free_symbolic(Symbolic::Ptr{Ptr{klu_symbolic}}, Common::Ptr{klu_common})::Cint
end

function klu_l_free_symbolic(arg1, arg2)
    @ccall libklu.klu_l_free_symbolic(arg1::Ptr{Ptr{klu_l_symbolic}}, arg2::Ptr{klu_l_common})::Cint
end

function klu_free_numeric(Numeric, Common)
    @ccall libklu.klu_free_numeric(Numeric::Ptr{Ptr{klu_numeric}}, Common::Ptr{klu_common})::Cint
end

function klu_z_free_numeric(Numeric, Common)
    @ccall libklu.klu_z_free_numeric(Numeric::Ptr{Ptr{klu_numeric}}, Common::Ptr{klu_common})::Cint
end

function klu_l_free_numeric(arg1, arg2)
    @ccall libklu.klu_l_free_numeric(arg1::Ptr{Ptr{klu_l_numeric}}, arg2::Ptr{klu_l_common})::Cint
end

function klu_zl_free_numeric(arg1, arg2)
    @ccall libklu.klu_zl_free_numeric(arg1::Ptr{Ptr{klu_l_numeric}}, arg2::Ptr{klu_l_common})::Cint
end

function klu_sort(Symbolic, Numeric, Common)
    @ccall libklu.klu_sort(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_z_sort(Symbolic, Numeric, Common)
    @ccall libklu.klu_z_sort(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_l_sort(arg1, arg2, arg3)
    @ccall libklu.klu_l_sort(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Ptr{klu_l_common})::Cint
end

function klu_zl_sort(arg1, arg2, arg3)
    @ccall libklu.klu_zl_sort(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Ptr{klu_l_common})::Cint
end

function klu_flops(Symbolic, Numeric, Common)
    @ccall libklu.klu_flops(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_z_flops(Symbolic, Numeric, Common)
    @ccall libklu.klu_z_flops(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_l_flops(arg1, arg2, arg3)
    @ccall libklu.klu_l_flops(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Ptr{klu_l_common})::Cint
end

function klu_zl_flops(arg1, arg2, arg3)
    @ccall libklu.klu_zl_flops(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Ptr{klu_l_common})::Cint
end

function klu_rgrowth(Ap, Ai, Ax, Symbolic, Numeric, Common)
    @ccall libklu.klu_rgrowth(Ap::Ptr{Int32}, Ai::Ptr{Int32}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_z_rgrowth(Ap, Ai, Ax, Symbolic, Numeric, Common)
    @ccall libklu.klu_z_rgrowth(Ap::Ptr{Int32}, Ai::Ptr{Int32}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_l_rgrowth(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_l_rgrowth(arg1::Ptr{Int64}, arg2::Ptr{Int64}, arg3::Ptr{Cdouble}, arg4::Ptr{klu_l_symbolic}, arg5::Ptr{klu_l_numeric}, arg6::Ptr{klu_l_common})::Cint
end

function klu_zl_rgrowth(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_zl_rgrowth(arg1::Ptr{Int64}, arg2::Ptr{Int64}, arg3::Ptr{Cdouble}, arg4::Ptr{klu_l_symbolic}, arg5::Ptr{klu_l_numeric}, arg6::Ptr{klu_l_common})::Cint
end

function klu_condest(Ap, Ax, Symbolic, Numeric, Common)
    @ccall libklu.klu_condest(Ap::Ptr{Int32}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_z_condest(Ap, Ax, Symbolic, Numeric, Common)
    @ccall libklu.klu_z_condest(Ap::Ptr{Int32}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_l_condest(arg1, arg2, arg3, arg4, arg5)
    @ccall libklu.klu_l_condest(arg1::Ptr{Int64}, arg2::Ptr{Cdouble}, arg3::Ptr{klu_l_symbolic}, arg4::Ptr{klu_l_numeric}, arg5::Ptr{klu_l_common})::Cint
end

function klu_zl_condest(arg1, arg2, arg3, arg4, arg5)
    @ccall libklu.klu_zl_condest(arg1::Ptr{Int64}, arg2::Ptr{Cdouble}, arg3::Ptr{klu_l_symbolic}, arg4::Ptr{klu_l_numeric}, arg5::Ptr{klu_l_common})::Cint
end

function klu_rcond(Symbolic, Numeric, Common)
    @ccall libklu.klu_rcond(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_z_rcond(Symbolic, Numeric, Common)
    @ccall libklu.klu_z_rcond(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_l_rcond(arg1, arg2, arg3)
    @ccall libklu.klu_l_rcond(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Ptr{klu_l_common})::Cint
end

function klu_zl_rcond(arg1, arg2, arg3)
    @ccall libklu.klu_zl_rcond(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Ptr{klu_l_common})::Cint
end

function klu_scale(scale, n, Ap, Ai, Ax, Rs, W, Common)
    @ccall libklu.klu_scale(scale::Cint, n::Int32, Ap::Ptr{Int32}, Ai::Ptr{Int32}, Ax::Ptr{Cdouble}, Rs::Ptr{Cdouble}, W::Ptr{Int32}, Common::Ptr{klu_common})::Cint
end

function klu_z_scale(scale, n, Ap, Ai, Ax, Rs, W, Common)
    @ccall libklu.klu_z_scale(scale::Cint, n::Int32, Ap::Ptr{Int32}, Ai::Ptr{Int32}, Ax::Ptr{Cdouble}, Rs::Ptr{Cdouble}, W::Ptr{Int32}, Common::Ptr{klu_common})::Cint
end

function klu_l_scale(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libklu.klu_l_scale(arg1::Cint, arg2::Int64, arg3::Ptr{Int64}, arg4::Ptr{Int64}, arg5::Ptr{Cdouble}, arg6::Ptr{Cdouble}, arg7::Ptr{Int64}, arg8::Ptr{klu_l_common})::Cint
end

function klu_zl_scale(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libklu.klu_zl_scale(arg1::Cint, arg2::Int64, arg3::Ptr{Int64}, arg4::Ptr{Int64}, arg5::Ptr{Cdouble}, arg6::Ptr{Cdouble}, arg7::Ptr{Int64}, arg8::Ptr{klu_l_common})::Cint
end

function klu_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx, P, Q, Rs, R, Common)
    @ccall libklu.klu_extract(Numeric::Ptr{klu_numeric}, Symbolic::Ptr{klu_symbolic}, Lp::Ptr{Int32}, Li::Ptr{Int32}, Lx::Ptr{Cdouble}, Up::Ptr{Int32}, Ui::Ptr{Int32}, Ux::Ptr{Cdouble}, Fp::Ptr{Int32}, Fi::Ptr{Int32}, Fx::Ptr{Cdouble}, P::Ptr{Int32}, Q::Ptr{Int32}, Rs::Ptr{Cdouble}, R::Ptr{Int32}, Common::Ptr{klu_common})::Cint
end

function klu_z_extract(Numeric, Symbolic, Lp, Li, Lx, Lz, Up, Ui, Ux, Uz, Fp, Fi, Fx, Fz, P, Q, Rs, R, Common)
    @ccall libklu.klu_z_extract(Numeric::Ptr{klu_numeric}, Symbolic::Ptr{klu_symbolic}, Lp::Ptr{Int32}, Li::Ptr{Int32}, Lx::Ptr{Cdouble}, Lz::Ptr{Cdouble}, Up::Ptr{Int32}, Ui::Ptr{Int32}, Ux::Ptr{Cdouble}, Uz::Ptr{Cdouble}, Fp::Ptr{Int32}, Fi::Ptr{Int32}, Fx::Ptr{Cdouble}, Fz::Ptr{Cdouble}, P::Ptr{Int32}, Q::Ptr{Int32}, Rs::Ptr{Cdouble}, R::Ptr{Int32}, Common::Ptr{klu_common})::Cint
end

function klu_l_extract(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14, arg15, arg16)
    @ccall libklu.klu_l_extract(arg1::Ptr{klu_l_numeric}, arg2::Ptr{klu_l_symbolic}, arg3::Ptr{Int64}, arg4::Ptr{Int64}, arg5::Ptr{Cdouble}, arg6::Ptr{Int64}, arg7::Ptr{Int64}, arg8::Ptr{Cdouble}, arg9::Ptr{Int64}, arg10::Ptr{Int64}, arg11::Ptr{Cdouble}, arg12::Ptr{Int64}, arg13::Ptr{Int64}, arg14::Ptr{Cdouble}, arg15::Ptr{Int64}, arg16::Ptr{klu_l_common})::Cint
end

function klu_zl_extract(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19)
    @ccall libklu.klu_zl_extract(arg1::Ptr{klu_l_numeric}, arg2::Ptr{klu_l_symbolic}, arg3::Ptr{Int64}, arg4::Ptr{Int64}, arg5::Ptr{Cdouble}, arg6::Ptr{Cdouble}, arg7::Ptr{Int64}, arg8::Ptr{Int64}, arg9::Ptr{Cdouble}, arg10::Ptr{Cdouble}, arg11::Ptr{Int64}, arg12::Ptr{Int64}, arg13::Ptr{Cdouble}, arg14::Ptr{Cdouble}, arg15::Ptr{Int64}, arg16::Ptr{Int64}, arg17::Ptr{Cdouble}, arg18::Ptr{Int64}, arg19::Ptr{klu_l_common})::Cint
end

function klu_malloc(n, size, Common)
    @ccall libklu.klu_malloc(n::Csize_t, size::Csize_t, Common::Ptr{klu_common})::Ptr{Cvoid}
end

function klu_free(p, n, size, Common)
    @ccall libklu.klu_free(p::Ptr{Cvoid}, n::Csize_t, size::Csize_t, Common::Ptr{klu_common})::Ptr{Cvoid}
end

function klu_realloc(nnew, nold, size, p, Common)
    @ccall libklu.klu_realloc(nnew::Csize_t, nold::Csize_t, size::Csize_t, p::Ptr{Cvoid}, Common::Ptr{klu_common})::Ptr{Cvoid}
end

function klu_l_malloc(arg1, arg2, arg3)
    @ccall libklu.klu_l_malloc(arg1::Csize_t, arg2::Csize_t, arg3::Ptr{klu_l_common})::Ptr{Cvoid}
end

function klu_l_free(arg1, arg2, arg3, arg4)
    @ccall libklu.klu_l_free(arg1::Ptr{Cvoid}, arg2::Csize_t, arg3::Csize_t, arg4::Ptr{klu_l_common})::Ptr{Cvoid}
end

function klu_l_realloc(arg1, arg2, arg3, arg4, arg5)
    @ccall libklu.klu_l_realloc(arg1::Csize_t, arg2::Csize_t, arg3::Csize_t, arg4::Ptr{Cvoid}, arg5::Ptr{klu_l_common})::Ptr{Cvoid}
end

const KLU_OK = 0

const KLU_SINGULAR = 1

const KLU_OUT_OF_MEMORY = -2

const KLU_INVALID = -3

const KLU_TOO_LARGE = -4

const KLU_DATE = "June 16, 2023"

const KLU_MAIN_VERSION = 2

const KLU_SUB_VERSION = 0

const KLU_SUBSUB_VERSION = 4

const KLU_VERSION = KLU_VERSION_CODE(KLU_MAIN_VERSION, KLU_SUB_VERSION)

