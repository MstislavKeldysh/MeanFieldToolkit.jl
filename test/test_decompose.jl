##### Analytic checks of the mean-field decompositions against textbook Hartree-Fock(-Bogoliubov)

@testset "Decompositions" begin
    U   =   2.0
    Chi =   ComplexF64[0.3 0.1+0.2im; 0.1-0.2im 0.6]     ##### generic Hermitian on-site <f^dag_a f_b>

    @testset "IntraQuarticToHopping: Hubbard Hartree + Fock" begin
        t   =   IntraQuarticToHopping(U .* Hubbard, Dict("ii" => Chi))["ii"]
        ##### U n_up n_down -> U <n_dn> n_up + U <n_up> n_dn - U <f^dag_dn f_up> f^dag_up f_dn - U <f^dag_up f_dn> f^dag_dn f_up
        @test t ≈ [U * Chi[2, 2]  -U * Chi[2, 1];
                  -U * Chi[1, 2]   U * Chi[1, 1]]
        @test t ≈ adjoint(t)
    end

    @testset "IntraQuarticToPairing: Hubbard" begin
        Delta   =   ComplexF64[0.0 0.4+0.1im; -0.4-0.1im 0.0]
        p       =   IntraQuarticToPairing(U .* Hubbard, Dict("ii" => Delta))["ii"]
        @test p ≈ [0.0  U * conj(Delta[1, 2]); 0.0 0.0]
    end

    @testset "InterQuarticToHopping: density-density V n_i n_j" begin
        V       =   0.7
        DD      =   DensityToPartonCoupling(Matrix{Float64}(I, 2, 2), Matrix{Float64}(I, 2, 2))
        Chi_jj  =   ComplexF64[0.2 0.0; 0.0 0.5]
        Chi_ij  =   ComplexF64[0.1+0.3im 0.05; -0.2im 0.4]
        t       =   InterQuarticToHopping(V .* DD, Dict("ii" => Chi, "jj" => Chi_jj, "ij" => Chi_ij))
        @test t["ii"] ≈ V * tr(Chi_jj) * I          ##### Hartree on i from the density on j
        @test t["jj"] ≈ V * tr(Chi) * I             ##### Hartree on j from the density on i
        @test t["ij"] ≈ -V .* conj.(Chi_ij)         ##### Fock on the bond
    end
end
