using Bigleaf
using Documenter, Latexify, DataFrames

# allow plot to work without display
# https://discourse.julialang.org/t/generation-of-documentation-fails-qt-qpa-xcb-could-not-connect-to-display/60988/2
ENV["GKSwstype"] = "100"

DocMeta.setdocmeta!(Bigleaf, :DocTestSetup, :(using Bigleaf, Latexify, DataFrames); 
    recursive=true, warn=false)
doctest(Bigleaf, manual = false)

makedocs(;
    # modules=[Bigleaf], # uncomment to show warnings on non-included docstrings
    authors="Thomas Wutzler <twutz@bgc-jena.mpg.de>, Jürgen Knauer <Juergen.Knauer@csiro.au> and contributors",
    #repo="https://github.com/earthyscience/Bigleaf.jl/blob/{commit}{path}#{line}",
    repo = Remotes.GitHub("earthyscience", "Bigleaf.jl"),
    sitename="Bigleaf.jl",
    doctestfilters=[r".*Info.*"],
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://earthyscience.github.io/Bigleaf.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        #"Meteorological variables" => "metorological_variables.md",
        #"Unit conversions" => "unit_conversions.md",
        "Walkthrough" => "walkthrough.md",
        hide("metorological_variables.md"),
        hide("stability_correction.md"),
        hide("surface_roughness.md"),
        hide("boundary_layer_conductance.md"),
        hide("surface_conductance.md"),
        hide("aerodynamic_conductance.md"),
        hide("evapotranspiration.md"),
        hide("surface_conductance.md"),
        hide("global_radiation.md"),
        hide("unit_conversions.md"),
        hide("BigleafConstants.md"),
        hide("filter_data.md"),
        #"Index" => "autodocs.md",
        ],
)

deploydocs(;
    repo="github.com/earthyscience/Bigleaf.jl",
    devbranch="main",
)
