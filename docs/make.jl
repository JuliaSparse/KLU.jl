using Documenter
using KLU

makedocs(
    sitename = "KLU",
    format = Documenter.HTML(),
    modules = [KLU];
    warnonly = [:missing_docs, :cross_references],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/JuliaSparse/KLU.jl.git",
    devbranch="main"
)
