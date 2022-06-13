mutable struct klu_symbolic
    symmetry::Cdouble
    est_flops::Cdouble
    lnz::Cdouble
    unz::Cdouble
    Lnz::Ptr{Cdouble}
    n::Cint
    nz::Cint
    P::Ptr{Cint}
    Q::Ptr{Cint}
    R::Ptr{Cint}
    nzoff::Cint
    nblocks::Cint
    maxblock::Cint
    ordering::Cint
    do_btf::Cint
    structural_rank::Cint
    klu_symbolic() = new()
end

mutable struct klu_l_symbolic
    symmetry::Cdouble
    est_flops::Cdouble
    lnz::Cdouble
    unz::Cdouble
    Lnz::Ptr{Cdouble}
    n::Clonglong
    nz::Clonglong
    P::Ptr{Clonglong}
    Q::Ptr{Clonglong}
    R::Ptr{Clonglong}
    nzoff::Clonglong
    nblocks::Clonglong
    maxblock::Clonglong
    ordering::Clonglong
    do_btf::Clonglong
    structural_rank::Clonglong
    klu_l_symbolic() = new()
end

mutable struct klu_numeric
    n::Cint
    nblocks::Cint
    lnz::Cint
    unz::Cint
    max_lnz_block::Cint
    max_unz_block::Cint
    Pnum::Ptr{Cint}
    Pinv::Ptr{Cint}
    Lip::Ptr{Cint}
    Uip::Ptr{Cint}
    Llen::Ptr{Cint}
    Ulen::Ptr{Cint}
    LUbx::Ptr{Ptr{Cvoid}}
    LUsize::Ptr{Csize_t}
    Udiag::Ptr{Cvoid}
    Rs::Ptr{Cdouble}
    worksize::Csize_t
    Work::Ptr{Cvoid}
    Xwork::Ptr{Cvoid}
    Iwork::Ptr{Cint}
    Offp::Ptr{Cint}
    Offi::Ptr{Cint}
    Offx::Ptr{Cvoid}
    nzoff::Cint
    klu_numeric() = new()
end

mutable struct klu_l_numeric
    n::Clonglong
    nblocks::Clonglong
    lnz::Clonglong
    unz::Clonglong
    max_lnz_block::Clonglong
    max_unz_block::Clonglong
    Pnum::Ptr{Clonglong}
    Pinv::Ptr{Clonglong}
    Lip::Ptr{Clonglong}
    Uip::Ptr{Clonglong}
    Llen::Ptr{Clonglong}
    Ulen::Ptr{Clonglong}
    LUbx::Ptr{Ptr{Cvoid}}
    LUsize::Ptr{Csize_t}
    Udiag::Ptr{Cvoid}
    Rs::Ptr{Cdouble}
    worksize::Csize_t
    Work::Ptr{Cvoid}
    Xwork::Ptr{Cvoid}
    Iwork::Ptr{Clonglong}
    Offp::Ptr{Clonglong}
    Offi::Ptr{Clonglong}
    Offx::Ptr{Cvoid}
    nzoff::Clonglong
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
    structural_rank::Cint
    numerical_rank::Cint
    singular_col::Cint
    noffdiag::Cint
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
    btf::Clonglong
    ordering::Clonglong
    scale::Clonglong
    user_order::Ptr{Cvoid}
    user_data::Ptr{Cvoid}
    halt_if_singular::Clonglong
    status::Clonglong
    nrealloc::Clonglong
    structural_rank::Clonglong
    numerical_rank::Clonglong
    singular_col::Clonglong
    noffdiag::Clonglong
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
    @ccall libklu.klu_l_defaults(Common::Ptr{klu_l_common})::Clonglong
end

function klu_analyze(n, Ap, Ai, Common)
    @ccall libklu.klu_analyze(n::Cint, Ap::Ptr{Cint}, Ai::Ptr{Cint}, Common::Ptr{klu_common})::Ptr{klu_symbolic}
end

function klu_l_analyze(arg1, arg2, arg3, Common)
    @ccall libklu.klu_l_analyze(arg1::Clonglong, arg2::Ptr{Clonglong}, arg3::Ptr{Clonglong}, Common::Ptr{klu_l_common})::Ptr{klu_l_symbolic}
end

function klu_analyze_given(n, Ap, Ai, P, Q, Common)
    @ccall libklu.klu_analyze_given(n::Cint, Ap::Ptr{Cint}, Ai::Ptr{Cint}, P::Ptr{Cint}, Q::Ptr{Cint}, Common::Ptr{klu_common})::Ptr{klu_symbolic}
end

function klu_l_analyze_given(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_l_analyze_given(arg1::Clonglong, arg2::Ptr{Clonglong}, arg3::Ptr{Clonglong}, arg4::Ptr{Clonglong}, arg5::Ptr{Clonglong}, arg6::Ptr{klu_l_common})::Ptr{klu_l_symbolic}
end

function klu_factor(Ap, Ai, Ax, Symbolic, Common)
    @ccall libklu.klu_factor(Ap::Ptr{Cint}, Ai::Ptr{Cint}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Common::Ptr{klu_common})::Ptr{klu_numeric}
end

function klu_z_factor(Ap, Ai, Ax, Symbolic, Common)
    @ccall libklu.klu_z_factor(Ap::Ptr{Cint}, Ai::Ptr{Cint}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Common::Ptr{klu_common})::Ptr{klu_numeric}
end

function klu_l_factor(arg1, arg2, arg3, arg4, arg5)
    @ccall libklu.klu_l_factor(arg1::Ptr{Clonglong}, arg2::Ptr{Clonglong}, arg3::Ptr{Cdouble}, arg4::Ptr{klu_l_symbolic}, arg5::Ptr{klu_l_common})::Ptr{klu_l_numeric}
end

function klu_zl_factor(arg1, arg2, arg3, arg4, arg5)
    @ccall libklu.klu_zl_factor(arg1::Ptr{Clonglong}, arg2::Ptr{Clonglong}, arg3::Ptr{Cdouble}, arg4::Ptr{klu_l_symbolic}, arg5::Ptr{klu_l_common})::Ptr{klu_l_numeric}
end

function klu_solve(Symbolic, Numeric, ldim, nrhs, B, Common)
    @ccall libklu.klu_solve(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, ldim::Cint, nrhs::Cint, B::Ptr{Cdouble}, Common::Ptr{klu_common})::Cint
end

function klu_z_solve(Symbolic, Numeric, ldim, nrhs, B, Common)
    @ccall libklu.klu_z_solve(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, ldim::Cint, nrhs::Cint, B::Ptr{Cdouble}, Common::Ptr{klu_common})::Cint
end

function klu_l_solve(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_l_solve(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Clonglong, arg4::Clonglong, arg5::Ptr{Cdouble}, arg6::Ptr{klu_l_common})::Clonglong
end

function klu_zl_solve(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_zl_solve(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Clonglong, arg4::Clonglong, arg5::Ptr{Cdouble}, arg6::Ptr{klu_l_common})::Clonglong
end

function klu_tsolve(Symbolic, Numeric, ldim, nrhs, B, Common)
    @ccall libklu.klu_tsolve(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, ldim::Cint, nrhs::Cint, B::Ptr{Cdouble}, Common::Ptr{klu_common})::Cint
end

function klu_z_tsolve(Symbolic, Numeric, ldim, nrhs, B, conj_solve, Common)
    @ccall libklu.klu_z_tsolve(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, ldim::Cint, nrhs::Cint, B::Ptr{Cdouble}, conj_solve::Cint, Common::Ptr{klu_common})::Cint
end

function klu_l_tsolve(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_l_tsolve(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Clonglong, arg4::Clonglong, arg5::Ptr{Cdouble}, arg6::Ptr{klu_l_common})::Clonglong
end

function klu_zl_tsolve(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    @ccall libklu.klu_zl_tsolve(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Clonglong, arg4::Clonglong, arg5::Ptr{Cdouble}, arg6::Clonglong, arg7::Ptr{klu_l_common})::Clonglong
end

function klu_refactor(Ap, Ai, Ax, Symbolic, Numeric, Common)
    @ccall libklu.klu_refactor(Ap::Ptr{Cint}, Ai::Ptr{Cint}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_z_refactor(Ap, Ai, Ax, Symbolic, Numeric, Common)
    @ccall libklu.klu_z_refactor(Ap::Ptr{Cint}, Ai::Ptr{Cint}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_l_refactor(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_l_refactor(arg1::Ptr{Clonglong}, arg2::Ptr{Clonglong}, arg3::Ptr{Cdouble}, arg4::Ptr{klu_l_symbolic}, arg5::Ptr{klu_l_numeric}, arg6::Ptr{klu_l_common})::Clonglong
end

function klu_zl_refactor(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_zl_refactor(arg1::Ptr{Clonglong}, arg2::Ptr{Clonglong}, arg3::Ptr{Cdouble}, arg4::Ptr{klu_l_symbolic}, arg5::Ptr{klu_l_numeric}, arg6::Ptr{klu_l_common})::Clonglong
end

function klu_free_symbolic(Symbolic, Common)
    @ccall libklu.klu_free_symbolic(Symbolic::Ptr{Ptr{klu_symbolic}}, Common::Ptr{klu_common})::Cint
end

function klu_l_free_symbolic(arg1, arg2)
    @ccall libklu.klu_l_free_symbolic(arg1::Ptr{Ptr{klu_l_symbolic}}, arg2::Ptr{klu_l_common})::Clonglong
end

function klu_free_numeric(Numeric, Common)
    @ccall libklu.klu_free_numeric(Numeric::Ptr{Ptr{klu_numeric}}, Common::Ptr{klu_common})::Cint
end

function klu_z_free_numeric(Numeric, Common)
    @ccall libklu.klu_z_free_numeric(Numeric::Ptr{Ptr{klu_numeric}}, Common::Ptr{klu_common})::Cint
end

function klu_l_free_numeric(arg1, arg2)
    @ccall libklu.klu_l_free_numeric(arg1::Ptr{Ptr{klu_l_numeric}}, arg2::Ptr{klu_l_common})::Clonglong
end

function klu_zl_free_numeric(arg1, arg2)
    @ccall libklu.klu_zl_free_numeric(arg1::Ptr{Ptr{klu_l_numeric}}, arg2::Ptr{klu_l_common})::Clonglong
end

function klu_sort(Symbolic, Numeric, Common)
    @ccall libklu.klu_sort(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_z_sort(Symbolic, Numeric, Common)
    @ccall libklu.klu_z_sort(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_l_sort(arg1, arg2, arg3)
    @ccall libklu.klu_l_sort(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Ptr{klu_l_common})::Clonglong
end

function klu_zl_sort(arg1, arg2, arg3)
    @ccall libklu.klu_zl_sort(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Ptr{klu_l_common})::Clonglong
end

function klu_flops(Symbolic, Numeric, Common)
    @ccall libklu.klu_flops(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_z_flops(Symbolic, Numeric, Common)
    @ccall libklu.klu_z_flops(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_l_flops(arg1, arg2, arg3)
    @ccall libklu.klu_l_flops(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Ptr{klu_l_common})::Clonglong
end

function klu_zl_flops(arg1, arg2, arg3)
    @ccall libklu.klu_zl_flops(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Ptr{klu_l_common})::Clonglong
end

function klu_rgrowth(Ap, Ai, Ax, Symbolic, Numeric, Common)
    @ccall libklu.klu_rgrowth(Ap::Ptr{Cint}, Ai::Ptr{Cint}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_z_rgrowth(Ap, Ai, Ax, Symbolic, Numeric, Common)
    @ccall libklu.klu_z_rgrowth(Ap::Ptr{Cint}, Ai::Ptr{Cint}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_l_rgrowth(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_l_rgrowth(arg1::Ptr{Clonglong}, arg2::Ptr{Clonglong}, arg3::Ptr{Cdouble}, arg4::Ptr{klu_l_symbolic}, arg5::Ptr{klu_l_numeric}, arg6::Ptr{klu_l_common})::Clonglong
end

function klu_zl_rgrowth(arg1, arg2, arg3, arg4, arg5, arg6)
    @ccall libklu.klu_zl_rgrowth(arg1::Ptr{Clonglong}, arg2::Ptr{Clonglong}, arg3::Ptr{Cdouble}, arg4::Ptr{klu_l_symbolic}, arg5::Ptr{klu_l_numeric}, arg6::Ptr{klu_l_common})::Clonglong
end

function klu_condest(Ap, Ax, Symbolic, Numeric, Common)
    @ccall libklu.klu_condest(Ap::Ptr{Cint}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_z_condest(Ap, Ax, Symbolic, Numeric, Common)
    @ccall libklu.klu_z_condest(Ap::Ptr{Cint}, Ax::Ptr{Cdouble}, Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_l_condest(arg1, arg2, arg3, arg4, arg5)
    @ccall libklu.klu_l_condest(arg1::Ptr{Clonglong}, arg2::Ptr{Cdouble}, arg3::Ptr{klu_l_symbolic}, arg4::Ptr{klu_l_numeric}, arg5::Ptr{klu_l_common})::Clonglong
end

function klu_zl_condest(arg1, arg2, arg3, arg4, arg5)
    @ccall libklu.klu_zl_condest(arg1::Ptr{Clonglong}, arg2::Ptr{Cdouble}, arg3::Ptr{klu_l_symbolic}, arg4::Ptr{klu_l_numeric}, arg5::Ptr{klu_l_common})::Clonglong
end

function klu_rcond(Symbolic, Numeric, Common)
    @ccall libklu.klu_rcond(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_z_rcond(Symbolic, Numeric, Common)
    @ccall libklu.klu_z_rcond(Symbolic::Ptr{klu_symbolic}, Numeric::Ptr{klu_numeric}, Common::Ptr{klu_common})::Cint
end

function klu_l_rcond(arg1, arg2, arg3)
    @ccall libklu.klu_l_rcond(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Ptr{klu_l_common})::Clonglong
end

function klu_zl_rcond(arg1, arg2, arg3)
    @ccall libklu.klu_zl_rcond(arg1::Ptr{klu_l_symbolic}, arg2::Ptr{klu_l_numeric}, arg3::Ptr{klu_l_common})::Clonglong
end

function klu_scale(scale, n, Ap, Ai, Ax, Rs, W, Common)
    @ccall libklu.klu_scale(scale::Cint, n::Cint, Ap::Ptr{Cint}, Ai::Ptr{Cint}, Ax::Ptr{Cdouble}, Rs::Ptr{Cdouble}, W::Ptr{Cint}, Common::Ptr{klu_common})::Cint
end

function klu_z_scale(scale, n, Ap, Ai, Ax, Rs, W, Common)
    @ccall libklu.klu_z_scale(scale::Cint, n::Cint, Ap::Ptr{Cint}, Ai::Ptr{Cint}, Ax::Ptr{Cdouble}, Rs::Ptr{Cdouble}, W::Ptr{Cint}, Common::Ptr{klu_common})::Cint
end

function klu_l_scale(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libklu.klu_l_scale(arg1::Clonglong, arg2::Clonglong, arg3::Ptr{Clonglong}, arg4::Ptr{Clonglong}, arg5::Ptr{Cdouble}, arg6::Ptr{Cdouble}, arg7::Ptr{Clonglong}, arg8::Ptr{klu_l_common})::Clonglong
end

function klu_zl_scale(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    @ccall libklu.klu_zl_scale(arg1::Clonglong, arg2::Clonglong, arg3::Ptr{Clonglong}, arg4::Ptr{Clonglong}, arg5::Ptr{Cdouble}, arg6::Ptr{Cdouble}, arg7::Ptr{Clonglong}, arg8::Ptr{klu_l_common})::Clonglong
end

function klu_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx, P, Q, Rs, R, Common)
    @ccall libklu.klu_extract(Numeric::Ptr{klu_numeric}, Symbolic::Ptr{klu_symbolic}, Lp::Ptr{Cint}, Li::Ptr{Cint}, Lx::Ptr{Cdouble}, Up::Ptr{Cint}, Ui::Ptr{Cint}, Ux::Ptr{Cdouble}, Fp::Ptr{Cint}, Fi::Ptr{Cint}, Fx::Ptr{Cdouble}, P::Ptr{Cint}, Q::Ptr{Cint}, Rs::Ptr{Cdouble}, R::Ptr{Cint}, Common::Ptr{klu_common})::Cint
end

function klu_z_extract(Numeric, Symbolic, Lp, Li, Lx, Lz, Up, Ui, Ux, Uz, Fp, Fi, Fx, Fz, P, Q, Rs, R, Common)
    @ccall libklu.klu_z_extract(Numeric::Ptr{klu_numeric}, Symbolic::Ptr{klu_symbolic}, Lp::Ptr{Cint}, Li::Ptr{Cint}, Lx::Ptr{Cdouble}, Lz::Ptr{Cdouble}, Up::Ptr{Cint}, Ui::Ptr{Cint}, Ux::Ptr{Cdouble}, Uz::Ptr{Cdouble}, Fp::Ptr{Cint}, Fi::Ptr{Cint}, Fx::Ptr{Cdouble}, Fz::Ptr{Cdouble}, P::Ptr{Cint}, Q::Ptr{Cint}, Rs::Ptr{Cdouble}, R::Ptr{Cint}, Common::Ptr{klu_common})::Cint
end

function klu_l_extract(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14, arg15, arg16)
    @ccall libklu.klu_l_extract(arg1::Ptr{klu_l_numeric}, arg2::Ptr{klu_l_symbolic}, arg3::Ptr{Clonglong}, arg4::Ptr{Clonglong}, arg5::Ptr{Cdouble}, arg6::Ptr{Clonglong}, arg7::Ptr{Clonglong}, arg8::Ptr{Cdouble}, arg9::Ptr{Clonglong}, arg10::Ptr{Clonglong}, arg11::Ptr{Cdouble}, arg12::Ptr{Clonglong}, arg13::Ptr{Clonglong}, arg14::Ptr{Cdouble}, arg15::Ptr{Clonglong}, arg16::Ptr{klu_l_common})::Clonglong
end

function klu_zl_extract(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19)
    @ccall libklu.klu_zl_extract(arg1::Ptr{klu_l_numeric}, arg2::Ptr{klu_l_symbolic}, arg3::Ptr{Clonglong}, arg4::Ptr{Clonglong}, arg5::Ptr{Cdouble}, arg6::Ptr{Cdouble}, arg7::Ptr{Clonglong}, arg8::Ptr{Clonglong}, arg9::Ptr{Cdouble}, arg10::Ptr{Cdouble}, arg11::Ptr{Clonglong}, arg12::Ptr{Clonglong}, arg13::Ptr{Cdouble}, arg14::Ptr{Cdouble}, arg15::Ptr{Clonglong}, arg16::Ptr{Clonglong}, arg17::Ptr{Cdouble}, arg18::Ptr{Clonglong}, arg19::Ptr{klu_l_common})::Clonglong
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

const KLU_DATE = "Mar 12, 2018"

const KLU_MAIN_VERSION = 1

const KLU_SUB_VERSION = 3

const KLU_SUBSUB_VERSION = 9

const KLU_VERSION = KLU_VERSION_CODE(KLU_MAIN_VERSION, KLU_SUB_VERSION)

