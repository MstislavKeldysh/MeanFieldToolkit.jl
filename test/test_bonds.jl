@testset "GetBondDictionary" begin
    A       =   ComplexF64[1.0 2.0im; 3.0 4.0]
    Z       =   zeros(ComplexF64, 2, 2)
    bond    =   (1, 2, [1, 0])
    reverse =   (2, 1, [-1, 0])
    onsite  =   (1, 1, [0, 0])

    ##### only the bond itself is in the lookup
    d   =   GetBondDictionary(Dict{Tuple, Matrix{ComplexF64}}(bond => A), bond, 2)
    @test d["ij"] == A
    @test d["ii"] == Z && d["jj"] == Z

    ##### only the reverse bond is in the lookup: recovered by Hermitian conjugation
    d   =   GetBondDictionary(Dict{Tuple, Matrix{ComplexF64}}(reverse => A), bond, 2)
    @test d["ij"] == adjoint(A)

    ##### on-site expectations are picked up for both ends of the bond
    d   =   GetBondDictionary(Dict{Tuple, Matrix{ComplexF64}}(bond => A, (1, 1, [0, 0]) => 2A, (2, 2, [0, 0]) => 3A), bond, 2)
    @test d["ii"] == 2A && d["jj"] == 3A

    ##### on-site bond: both lookups hit the same entry, so the /2 gives its Hermitian part.
    ##### This is the on-site exchange (Fock) amplitude and must not be zeroed.
    d   =   GetBondDictionary(Dict{Tuple, Matrix{ComplexF64}}(onsite => A), onsite, 2)
    @test d["ij"] ≈ (A + adjoint(A)) / 2
    @test d["ii"] == A && d["jj"] == A
end
