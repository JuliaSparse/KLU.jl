using KLU
using KLU: increment!, KLUITypes, decrement, klu!, KLUFactorization,
klu_analyze!, klu_factor!
using Test
using SparseArrays: SparseMatrixCSC, sparse, sprand, nnz, nonzeros
using LinearAlgebra

@testset "KLU Wrappers" begin
    Ap = increment!([0,4,1,1,2,2,0,1,2,3,4,4])
    Ai = increment!([0,4,0,2,1,2,1,4,3,2,1,2])
    Ax = [2.,1.,3.,4.,-1.,-3.,3.,6.,2.,1.,4.,2.]
    A0 = sparse(Ap, Ai, Ax)
    A1 = sparse(increment!([0,4,1,1,2,2,0,1,2,3,4,4]),
                increment!([0,4,0,2,1,2,1,4,3,2,1,2]),
                [2.,1.,3.,4.,-1.,-3.,3.,9.,2.,1.,4.,2.], 5, 5)
    @testset "Core functionality for $Tv elements" for Tv in (Float64, ComplexF64)
        @testset "Core functionality for $Ti indices" for Ti ∈ Base.uniontypes(KLUITypes)
            A = convert(SparseMatrixCSC{Tv, Ti}, A0)
            # test the raw vector construction method.
            klua = klu(size(A, 1), decrement(A.colptr), decrement(A.rowval), A.nzval)
            @test klua.F == klu(A).F
            @test nnz(klua) == 18
            R = Diagonal(Tv == ComplexF64 ? complex.(klua.Rs) : klua.Rs)
            @test R \ A[klua.p, klua.q] ≈ (klua.L * klua.U + klua.F)

            b = [8., 45., -3., 3., 19.]
            x = klua \ b
            @test x ≈ float([1:5;])
            @test A*x ≈ b

            z = complex.(b)
            x = ldiv!(klua, z)
            @test x ≈ float([1:5;])
            @test z === x
            # Can't match UMFPACK's ldiv!(<OUTPUT>, <KLU>, <INPUT>)
            # since klu_solve(<KLU>, <INPUT>) modifies <INPUT>, and has no field for <OUTPUT>.
            @test A*x ≈ b

            b = [8., 20., 13., 6., 17.]
            x = klua'\b
            @test x ≈ float([1:5;])

            @test A'*x ≈ b
            z = complex.(b)
            x = ldiv!(adjoint(klua), z)
            @test x ≈ float([1:5;])
            @test x === z

            @test A'*x ≈ b
            x = transpose(klua) \ b
            @test x ≈ float([1:5;])
            @test transpose(A) * x ≈ b

            x = ldiv!(transpose(klua), complex.(b))
            @test x ≈ float([1:5;])
            @test transpose(A) * x ≈ b

            @inferred klua\fill(1, size(A, 2))

            @testset "Permutation vectors" begin
                #Just to test this works, we'll use the existing permutation vectors.
                @test klu_analyze!(klua, klua.p, klua.q).common.status == 0
                klu_factor!(klua)
                x = klua \ b
                @test A*x ≈ b
            end
            @testset "Utility functions" begin
                K = KLUFactorization(A);
                @test size(K) == (5, 5)
                @test size(K, 3) == 1
                @test_throws ArgumentError K.symbolic
                @test_throws ArgumentError K.numeric
                klu_analyze!(K);
                @test K.symbolic.nz == 12
                @test K.symbolic.nzoff == 4
                klu_factor!(K)
                @test K.numeric.lnz == 7 == K.numeric.unz
                @test K.nblocks == 3
                @test propertynames(K) == (:lnz, :unz, :nzoff, :L, :U, :F, :q, :p, :Rs, :symbolic, :numeric,)
                @test_throws ArgumentError size(K, -1)
                @test size(K, 3) == 1
            end
            @testset "Refactorization" begin
                B = convert(SparseMatrixCSC{Tv, Ti}, A1)
                b = Tv[8., 45., -3., 3., 19.]
                F = klu(A)
                klu!(F, B)
                @test F\b ≈ B\b ≈ Matrix(B)\b
                @test klu!(F, B.nzval)\b ≈ B\b ≈ Matrix(B)\b #test just supplying nzval for recompute
            end

            @testset "Singular matrix" begin
                S = sparse([1 2; 0 0])
                S = convert(SparseMatrixCSC{Tv, Ti}, S)
                @test_throws SingularException klu(S)
            end
        end
    end
    @testset "REPL printing of KLU" begin
        A = sparse([1, 2], [1, 2], Float64[1.0, 1.0])
        F = klu(A)
        facstring = sprint((t, s) -> show(t, "text/plain", s), F)
        lstring = sprint((t, s) -> show(t, "text/plain", s), F.L)
        ustring = sprint((t, s) -> show(t, "text/plain", s), F.U)
        fstring = sprint((t, s) -> show(t, "text/plain", s), F.F)
        @test facstring == "$(summary(F))\nL factor:\n$lstring\nU factor:\n$ustring\nF factor:\n$fstring"
    end
end

@testset "Issue #4" begin
    A = SparseMatrixCSC(15, 15, 
    [1, 8, 12, 16, 21, 22, 27, 34, 37, 39, 41, 44, 47, 50, 55, 58],
    [3, 4, 6, 8, 11, 13, 15, 2, 6, 8, 15, 1, 2, 5, 12, 3, 5, 6, 7, 9, 14, 7, 10, 11, 13, 15, 4, 5, 6, 8, 10, 12, 13, 8, 12, 15, 2, 13, 1, 7, 3, 4, 14, 2, 4, 5, 1, 5, 9, 4, 7, 9, 11, 14, 13, 14, 15],
    [0.2775474841561938, 0.19549953706849155, 0.7221976371086005, 0.4339373082200655, 0.983079431343046, 0.10918778088879799, 0.3728676112188065, 0.9134045061777432, 0.14560891622463457, 0.7715210431553383, 0.2945501295372417, 0.6134722502122134, 0.8777181195348973, 0.6995382425541914, 0.9562490972786235, 0.27001502642215325, 0.8573661029146233, 0.13020432565448115, 0.9221068910751316, 0.17087414970038983, 0.7062975193151109, 0.7668596005709167, 0.46967704631299334, 0.31764226298560083, 0.39054386892157833, 0.36610203401046015, 0.16689896372140534, 0.9624322297755521, 0.1478381603984824, 0.45423514524961806, 0.5564610482579242, 0.1844671322175948, 0.0823893170743808, 0.25409993152712307, 0.10475245199943273, 0.5863595004922162, 0.14733912690513562, 0.6504152422320895, 0.4339054908933866, 0.27614384058497166, 0.4019619228414033, 0.8631491210976487, 0.20159747073826084, 0.3273367915690062, 0.23866880928640288, 0.9557759456784265, 0.016351125161178537, 0.5320355909884844, 0.9010930260468242, 0.3780686420068593, 0.6375477164214856, 0.9850645305956225, 0.5366242762582065, 0.08835652070698918, 0.9877090693717305, 0.9775298646022268, 0.9759511830494418]
    )
    @test klu(A).F.nzval ≈ [0.20623152093700403, 0.09038754099133628, 1.0]
end

# Regression tests for JuliaSparse/KLU.jl#12: the default `deepcopy` copied the raw libklu
# pointers, so a copy shared the C objects (and their scratch space) with its source: a
# refactorization of one changed the other, and either could free them under the other.
@testset "copy and deepcopy of a KLUFactorization (issue #12)" begin
    Ap = increment!([0,4,1,1,2,2,0,1,2,3,4,4])
    Ai = increment!([0,4,0,2,1,2,1,4,3,2,1,2])
    A0 = sparse(Ap, Ai, [2.,1.,3.,4.,-1.,-3.,3.,6.,2.,1.,4.,2.], 5, 5)
    A1 = sparse(Ap, Ai, [2.,1.,3.,4.,-1.,-3.,3.,9.,2.,1.,4.,2.], 5, 5)
    independent(C, K) = C !== K && C.lock !== K.lock && C.common !== K.common &&
        C.colptr !== K.colptr && C.rowval !== K.rowval && C.nzval !== K.nzval &&
        (C._symbolic == C_NULL || C._symbolic != K._symbolic) &&
        (C._numeric == C_NULL || C._numeric != K._numeric)
    copy_of_dropped_source(dup, A) = dup(klu(A))

    @testset "$dup for $Tv elements and $Ti indices" for dup in (copy, deepcopy),
            Tv in (Float64, ComplexF64), Ti in Base.uniontypes(KLUITypes)
        A = convert(SparseMatrixCSC{Tv, Ti}, A0)
        B = convert(SparseMatrixCSC{Tv, Ti}, A1)
        b = Tv[8., 45., -3., 3., 19.]
        Ab, Bb = Matrix(A) \ b, Matrix(B) \ b

        @testset "Factored source" begin
            K = klu(A)
            C = dup(K)
            @test C isa KLUFactorization{Tv, Ti}
            @test independent(C, K)
            @test C._symbolic != C_NULL && C._numeric != C_NULL
            @test issuccess(C)
            @test (C.colptr, C.rowval, C.nzval) == (K.colptr, K.rowval, K.nzval)
            @test (C.L, C.U, C.F, C.p, C.q, C.Rs) == (K.L, K.U, K.F, K.p, K.q, K.Rs)
            @test C \ b ≈ Ab
            @test C' \ b ≈ Matrix(A)' \ b
            @test transpose(C) \ b ≈ transpose(Matrix(A)) \ b
            # the lock of the source is not carried over, even when held
            Base.@lock K begin
                @test !islocked(dup(K))
            end
        end

        @testset "Using one does not affect the other" begin
            K = klu(A)
            C = dup(K)
            klu!(C, B)
            @test C \ b ≈ Bb
            @test K \ b ≈ Ab
            @test nonzeros(K) == nonzeros(A)
            klu_factor!(C) # frees and replaces the numeric object of the copy
            @test C \ b ≈ Bb
            @test K \ b ≈ Ab
            C = dup(K)
            klu!(K, B)
            @test K \ b ≈ Bb
            @test C \ b ≈ Ab
            @test nonzeros(C) == nonzeros(A)
        end

        @testset "Freeing one does not affect the other" begin
            K = klu(A)
            C = dup(K)
            KLU._free_numeric(K)
            KLU._free_symbolic(K)
            @test C \ b ≈ Ab
            C = copy_of_dropped_source(dup, A)
            GC.gc()
            GC.gc()
            @test C \ b ≈ Ab
            @test klu!(C, B) \ b ≈ Bb
            K = klu(A)
            foreach(_ -> dup(K), 1:20)
            GC.gc()
            GC.gc()
            @test K \ b ≈ Ab
            C = dup(dup(dup(K)))
            @test independent(C, K)
            K = nothing
            GC.gc()
            @test C \ b ≈ Ab
        end

        @testset "Unanalyzed source" begin
            K = KLUFactorization(A)
            C = dup(K)
            @test independent(C, K)
            @test C._symbolic == C_NULL && C._numeric == C_NULL
            @test C \ b ≈ Ab
            @test K._symbolic == C_NULL && K._numeric == C_NULL
        end

        @testset "Symbolic-only source" begin
            K = klu(copy(A); full_factor = false)
            nonzeros(K) .= zero(Tv) # the values are not known yet
            C = dup(K)
            @test independent(C, K)
            @test C._symbolic != C_NULL && C._numeric == C_NULL
            @test C.symbolic.nz == K.symbolic.nz
            @test C.nblocks == K.nblocks
            @test unsafe_wrap(Array, C.symbolic.Q, 5) == unsafe_wrap(Array, K.symbolic.Q, 5)
            nonzeros(C) .= nonzeros(B)
            klu_factor!(C)
            @test C \ b ≈ Bb
            @test K._numeric == C_NULL
            @test all(iszero, nonzeros(K))
        end

        @testset "Settings in common" begin
            K = KLUFactorization(A)
            K.common.tol = 0.5
            K.common.btf = 0
            K.common.ordering = 1
            K.common.scale = 2
            for F in (K, klu_analyze!(dup(K)), klu_factor!(dup(K)))
                C = dup(F)
                @test C.common.tol == 0.5
                @test C.common.btf == 0
                @test C.common.ordering == 1
                @test C.common.scale == 2
                @test C \ b ≈ Ab
                @test C.symbolic.do_btf == 0
                @test C.symbolic.ordering == 1
                @test C.nblocks == 1
            end
        end

        @testset "User-provided permutations" begin
            P, Q = Ti[4, 3, 2, 1, 0], Ti[0, 2, 1, 4, 3]
            K = KLUFactorization(A)
            K.common.btf = 0
            klu_analyze!(K, P, Q)
            C = dup(K)
            @test C._numeric == C_NULL
            @test C.symbolic.ordering == 2
            @test unsafe_wrap(Array, C.symbolic.P, 5) == P
            @test unsafe_wrap(Array, C.symbolic.Q, 5) == Q
            @test C \ b ≈ Ab
            # BTF further permutes the given ordering
            K = KLUFactorization(A)
            klu_analyze!(K, P, Q)
            klu_factor!(K)
            C = dup(K)
            @test C.symbolic.ordering == 2
            @test C \ b ≈ Ab
        end

        @testset "nzval changed after factorization" begin
            for scale in (0, 1, 2), btf in (0, 1)
                K = KLUFactorization(A)
                K.common.scale = scale
                K.common.btf = btf
                klu_factor!(K)
                nonzeros(K) .= nonzeros(B)
                C = dup(K)
                # The copy factors the current values, while its source keeps its factors.
                @test nonzeros(C) == nonzeros(B)
                @test independent(C, K)
                @test C \ b ≈ Bb
                @test K \ b ≈ Ab
            end
            # klu(A) may share the values of A, which can then change under it
            A2 = copy(A)
            K = klu(A2)
            nonzeros(A2) .= nonzeros(B)
            @test dup(K) \ b ≈ Bb
            @test K \ b ≈ Ab
        end

        @testset "After a refactorization" begin
            K = klu(A)
            klu!(K, B)
            C = dup(K)
            @test nonzeros(C) == nonzeros(B)
            @test C \ b ≈ Bb
            @test klu!(C, A) \ b ≈ Ab
            @test K \ b ≈ Bb
        end

        @testset "Singular factorizations" begin
            S = convert(SparseMatrixCSC{Tv, Ti}, sparse([1, 2], [1, 2], [1., 0.]))
            R = convert(SparseMatrixCSC{Tv, Ti}, sparse([1, 2], [1, 2], [1., 2.]))
            c = Tv[1., 2.]
            for K in (klu(S; allowsingular = true), klu(S; check = false))
                C = dup(K)
                @test independent(C, K)
                @test C._numeric != C_NULL
                @test C.common.status == K.common.status == KLU.KLU_SINGULAR
                @test C.common.numerical_rank == K.common.numerical_rank
                @test C.common.singular_col == K.common.singular_col
                @test !issuccess(C)
                @test issuccess(C; allowsingular = true)
                @test_throws SingularException klu_factor!(dup(K))
                @test klu!(C, R) \ c ≈ Matrix(R) \ c
                @test issuccess(C)
                @test !issuccess(K)
            end
            K = klu(S; allowsingular = true)
            nonzeros(K) .= nonzeros(R)
            @test dup(K) \ c ≈ Matrix(R) \ c
        end

        @testset "Status of the source" begin
            K = klu(A)
            K.common.status = KLU.KLU_INVALID
            C = dup(K)
            @test C._numeric != C_NULL
            @test C.common.status == KLU.KLU_INVALID
            # a failed analysis with check = false leaves nothing to copy but the status
            K = klu(2, Ti[0, 2, 2], Ti[0, 0], Tv[1, 1]; check = false, full_factor = false)
            C = dup(K)
            @test K._symbolic == C_NULL == C._symbolic
            @test C.common.status == K.common.status == KLU.KLU_INVALID
        end
    end

    @testset "deepcopy keeps the structure of the copied object" begin
        A = copy(A0)
        K = klu(A)
        @test nonzeros(K) === nonzeros(A)
        K2, K3, A2 = deepcopy((K, K, A))
        @test K2 === K3
        @test nonzeros(K2) === nonzeros(A2)
        @test nonzeros(K2) !== nonzeros(A)
        b = [8., 45., -3., 3., 19.]
        @test deepcopy(K') \ b ≈ K' \ b
        @test deepcopy(transpose(K)) \ b ≈ transpose(K) \ b
        @test parent(deepcopy(K'))._numeric != K._numeric
    end

    # This is the situation from the issue: one copy of a factorization per task.
    @testset "Copies used concurrently for $Tv elements and $Ti indices" for
            Tv in (Float64, ComplexF64), Ti in Base.uniontypes(KLUITypes)
        n = 400
        A = convert(SparseMatrixCSC{Tv, Ti}, sparse(1:n, 1:n, 4.0 .+ rand(n)) + sprand(n, n, 4 / n))
        K = klu(A)
        ntasks = 4 * max(Threads.nthreads(), 2)
        copies = [deepcopy(K) for _ in 1:ntasks]
        residuals = zeros(ntasks)
        @sync for (t, C) in enumerate(copies)
            Threads.@spawn begin
                At = copy(A)
                nonzeros(At) .*= t
                b = rand(Tv, n)
                worst = 0.0
                for _ in 1:50
                    klu!(C, At)
                    worst = max(worst, norm(At * (C \ b) - b) / norm(b))
                    worst = max(worst, norm(At' * (C' \ b) - b) / norm(b))
                end
                residuals[t] = worst
            end
        end
        @test all(residuals .< 1e-10)
        # Copying while the source is being used by other tasks.
        @sync for t in 1:ntasks
            Threads.@spawn begin
                b = rand(Tv, n)
                worst = 0.0
                for _ in 1:10
                    worst = max(worst, norm(A * (K \ b) - b) / norm(b))
                    worst = max(worst, norm(A * (copy(K) \ b) - b) / norm(b))
                end
                residuals[t] = worst
            end
        end
        @test all(residuals .< 1e-10)
    end
end

@testset "check=false" begin
    for A in sparse.((Float64[1 2; 0 0], ComplexF64[1 2; 0 0]))
        @test_throws SingularException klu(A)
        @test !issuccess(klu(A; check = false))
        @test issuccess(klu(A; allowsingular=true); allowsingular=true)
    end
end

@testset "full_factor = false" begin
    for T in (Float64, ComplexF64, Float32, ComplexF32)
        A = sparse(Matrix{T}([1 0; 1 1]))
        nonzeros(A) .= zero(T)
        F = klu(A; full_factor = false)
        # we do need both here--for non 64 bit types, the nzval of F is not the same as A's.
        nonzeros(A) .= one(T)
        nonzeros(F) .= one(T)
        klu_factor!(F)
        b = ones(T, 2)
        @test F\b ≈ A\b
    end
end
# Regression test for JuliaSparse/KLU.jl#38: libklu keeps scratch workspace inside the
# Numeric object, so unsynchronized concurrent solves on one factorization return garbage.
@testset "Concurrent use of a shared factorization (issue #38)" begin
    n = 400
    A = sparse(1:n, 1:n, 4.0 .+ rand(n)) + sprand(n, n, 4 / n)
    for Tv in (Float64, ComplexF64), Ti in Base.uniontypes(KLUITypes)
        As = convert(SparseMatrixCSC{Tv, Ti}, A)
        F = klu(As)
        ntasks = 4 * max(Threads.nthreads(), 2)
        residuals = zeros(ntasks)
        @sync for t in 1:ntasks
            Threads.@spawn begin
                b = rand(Tv, n)
                worst = 0.0
                for i in 1:100
                    # mix plain, adjoint, and transposed solves plus multi-rhs
                    x = F \ b
                    worst = max(worst, norm(As * x - b) / norm(b))
                    x = F' \ b
                    worst = max(worst, norm(As' * x - b) / norm(b))
                    x = transpose(F) \ b
                    worst = max(worst, norm(transpose(As) * x - b) / norm(b))
                    X = solve!(F, repeat(b, 1, 3))
                    worst = max(worst, norm(As * X .- b) / norm(b))
                end
                residuals[t] = worst
            end
        end
        @test all(residuals .< 1e-10)

        # Concurrent refactorization (same values) interleaved with solves must stay consistent.
        @sync begin
            Threads.@spawn for _ in 1:50
                klu!(F, As)
            end
            Threads.@spawn begin
                b = rand(Tv, n)
                for _ in 1:200
                    x = F \ b
                    @test norm(As * x - b) / norm(b) < 1e-10
                end
            end
        end
        # a caller may hold the lock across several operations
        @test !islocked(F)
        Base.@lock F begin
            @test islocked(F)
            b = rand(Tv, n)
            klu!(F, As)
            @test norm(As * (F \ b) - b) / norm(b) < 1e-10
        end
        @test !islocked(F)
    end
    # Independent factorizations on different threads share nothing.
    @sync for _ in 1:8
        Threads.@spawn begin
            F = klu(A)
            b = rand(n)
            for _ in 1:20
                @test norm(A * (F \ b) - b) / norm(b) < 1e-10
            end
        end
    end
end

@testset "klu_factor! on an already factored object releases the old factor" begin
    A = sparse([1, 2], [1, 2], [1.0, 2.0])
    F = klu(A)
    p1 = F._numeric
    @test p1 != C_NULL
    klu_factor!(F)               # must not leak p1 (setproperty! frees it)
    @test F._numeric != C_NULL
    @test F \ [1.0, 2.0] ≈ [1.0, 1.0]
end
