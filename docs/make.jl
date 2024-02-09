using Documenter
using KLU

makedocs(
    modules = [KLU],
    sitename = "Interface to the KLU solver from SuiteSparse",
    authors = "Raye Kimmerer, Viral B. Shah",
    format = Documenter.HTML(),
    pages = ["KLU" => "index.md",
             "Reference" => "functions_list.md"
             ],
    warnonly = [:missing_docs, :cross_references],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/JuliaSparse/KLU.jl.git",
    devbranch="main"
)
