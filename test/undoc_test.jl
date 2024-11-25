using Bigleaf
using Test

@testset "Docstrings" begin
    if VERSION >= v"1.11"
        @test_broken isempty(Docs.undocumented_names(Bigleaf))
    end
    #undoc = Docs.undocumented_names(Bigleaf)
    #@test_broken undoc == []
    #@test undoc = [:bar, :baz]
    #@test undoc = []
end
