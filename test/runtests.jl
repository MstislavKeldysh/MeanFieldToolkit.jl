using MeanFieldToolkit, TightBindingToolkit, FixedPointToolkit
using LinearAlgebra, JLD2, Logging
using Test

##### the solvers log every iteration; only show warnings and errors
global_logger(ConsoleLogger(stderr, Logging.Warn))

include("models.jl")

@testset "MeanFieldToolkit.jl" begin
    include("test_decompose.jl")
    include("test_bonds.jl")
    include("test_hubbard.jl")
    include("test_io.jl")
    include("test_bdg.jl")
end
