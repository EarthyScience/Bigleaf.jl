using Test, SafeTestsets
const GROUP = get(ENV, "GROUP", "All") # defined in in CI.yml

@time begin
    if GROUP == "All" || GROUP == "Basic"
        #@safetestset "bigleaf_constants Tests" include("test/bigleaf_constants_test.jl")
        @time @safetestset "bigleaf_constants Tests" include("bigleaf_constants_test.jl")
        #@safetestset "Tests" include("test/util_test.jl")
        @time @safetestset "util Tests" include("util_test.jl")
        #@safetestset "Tests" include("test/unit_conversions_test.jl")
        @time @safetestset "unit_conversions Tests" include("unit_conversions_test.jl")
        #@safetestset "Tests" include("test/filter_data_test.jl")
        @time @safetestset "filter_data Tests" include("filter_data_test.jl")
        # @safetestset "Tests" include("test/meteorological_variables_test.jl")
        @time @safetestset "meteorological_variables Tests" include("meteorological_variables_test.jl")
        # @safetestset "Tests" include("test/sun_position_test.jl")
        @time @safetestset "sun_position Tests" include("sun_position_test.jl")
        #@safetestset "Tests" include("test/potential_radiation_test.jl")
        @time @safetestset "potential_radiation Tests" include("potential_radiation_test.jl")
        #@safetestset "Tests" include("test/stability_correction_test.jl")
        @time @safetestset "stability_correction Tests" include("stability_correction_test.jl")
        #@safetestset "Tests" include("test/surface_roughness_test.jl")
        @time @safetestset "surface_roughness Tests" include("surface_roughness_test.jl")
        #@safetestset "Tests" include("test/boundary_layer_conductance_test.jl")
        @time @safetestset "boundary_layer_conductance Tests" include("boundary_layer_conductance_test.jl")
        #@safetestset "Tests" include("test/aerodynamic_conductance_test.jl")
        @time @safetestset "aerodynamic_conductance Tests" include("aerodynamic_conductance_test.jl")
        #@safetestset "Tests" include("test/surface_conductance_test.jl")
        @time @safetestset "surface_conductance Tests" include("surface_conductance_test.jl")
        #@safetestset "Tests" include("test/evapotranspiration_test.jl")
        @time @safetestset "evapotranspiration Tests" include("evapotranspiration_test.jl")
        #@safetestset "Tests" include("test/undoc_test.jl")
        @time @safetestset "undocumented objects" include("undoc_test.jl")
    end
end

