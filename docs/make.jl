using Documenter
using KLU

makedocs(
    sitename = "KLU",
    format = Documenter.HTML(),
    modules = [KLU]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/Wimmerer/KLU.jl.git"
)
