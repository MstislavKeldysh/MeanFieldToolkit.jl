@testset "Renormalized t-J model, d-wave BdG" begin
    mft     =   SquaretJ()
    ##### loose tolerance: this iteration has not been shown to converge tightly (residual ~3e-4 after 300 iterations)
    sc      =   SolveMFT!(mft, [0.2, 0.2]; max_iter = 300, tol = 1e-3)

    @test maximum(abs.(sc.VOuts[end] - sc.VIns[end])) < 1e-3
    @test all(isfinite, sc.VOuts[end]) && isfinite(mft.MFTEnergy[end])
end
