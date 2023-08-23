#! /bin/bash julia --project generator.jl
using Pkg
using Pkg.Artifacts
using Clang
using Clang.Generators
using Clang.Generators.JLLEnvs
using SuiteSparse_jll

cd(@__DIR__)

# headers
SuiteSparse_toml = joinpath(dirname(pathof(SuiteSparse_jll)), "..", "StdlibArtifacts.toml")
SuiteSparse_dir = Pkg.Artifacts.ensure_artifact_installed("SuiteSparse", SuiteSparse_toml)

include_dir = joinpath(SuiteSparse_dir, "include") |> normpath

klu_h = joinpath(include_dir, "klu.h")
@assert isfile(klu_h)

# load common option
options = load_options(joinpath(@__DIR__, "generator.toml"))

# run generator for all platforms
for target in JLLEnvs.JLL_ENV_TRIPLES
    @info "processing $target"

    options["general"]["output_file_path"] = joinpath(@__DIR__, "..", "lib", "$target.jl")

    args = get_default_args(target)
    push!(args, "-I$include_dir")
    if startswith(target, "x86_64") || startswith(target, "powerpc64le") || startswith(target, "aarch64")
        push!(args, "-DSUN64 -DLONGBLAS='long long'")
    end

    header_files = [klu_h]

    ctx = create_context(header_files, args, options)

    build!(ctx, BUILDSTAGE_NO_PRINTING)

    # custom rewriter
    function rewrite!(dag::ExprDAG)
        replace!(get_nodes(dag)) do node
            filename = normpath(Clang.get_filename(node.cursor))
            if !contains(filename, "klu")
                return ExprNode(node.id, Generators.Skip(), node.cursor, Expr[], node.adj)
            end
            return node
        end
    end

    rewrite!(ctx.dag)

    # print
    build!(ctx, BUILDSTAGE_PRINTING_ONLY)
end
