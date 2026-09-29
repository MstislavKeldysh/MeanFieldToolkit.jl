@testset "Saving, reading and resuming" begin
    dir             =   mktempdir()
    compact_keys    =   Set(["Iterations", "MFT_Energy", "Hopping_Order", "UC", "Gap", "mu", "Outputs", "Convergence"])

    @testset "debug = false writes only the compact output" begin
        file    =   joinpath(dir, "compact.jld2")
        SolveMFT!(SquareHubbard(kSize = 9), [0.1, 0.5, 0.5], file; max_iter = 30, tol = 1e-6, debug = false)
        @test Set(keys(load(file))) == compact_keys
        out     =   ReadMFT(file)
        @test out["Iterations"] isa Integer && length(out["Expectations"]) == 3 && out["MFT"] === nothing

        ##### random initial conditions; a checkpoint every iteration would show up if suppression failed
        file    =   joinpath(dir, "compact_random.jld2")
        SolveMFT!(SquareHubbard(kSize = 9), file; max_iter = 5, tol = 1e-12, checkpoint_interval = 1, debug = false)
        @test Set(keys(load(file))) == compact_keys
        @test ReadMFT(file)["Iterations"] == 5
        @test_throws ErrorException ResumeMFT!(file)
    end

    @testset "debug = true writes resumable checkpoints" begin
        file    =   joinpath(dir, "checkpoint.jld2")
        SolveMFT!(SquareHubbard(kSize = 9), [0.1, 0.5, 0.5], file; max_iter = 20, tol = 1e-12, checkpoint_interval = 5, debug = true)
        out     =   ReadMFT(file)
        @test out["MFT"] isa TBMFTModel && length(out["Expectations"]) == 3

        sc      =   ResumeMFT!(file; max_iter = 10, tol = 1e-12, checkpoint_interval = 5)
        @test sc isa SelfCons
        @test ReadMFT(file)["Iterations"] > out["Iterations"]
    end
end
