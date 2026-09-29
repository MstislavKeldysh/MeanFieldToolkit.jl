##### Small model builders shared by the tests (reduced versions of the scripts in Sample/)

const SpinVec   =   SpinMats(1 // 2)
const Hubbard   =   DensityToPartonCoupling([1.0 0.0; 0.0 0.0], [0.0 0.0; 0.0 1.0])   ##### n_up * n_down

##### Square lattice, two-site (Néel) unit cell, on-site Hubbard U, order parameters [s-hopping, Sz_1, Sz_2]
function SquareHubbard(; U::Float64 = 8.0, kSize::Int64 = 15, filling::Float64 = 0.5) :: TBMFTModel
    UC      =   UnitCell([[1.0, 1.0], [1.0, -1.0]], 2, 2)
    AddBasisSite!(UC, [0.0, 0.0])
    AddBasisSite!(UC, [1.0, 0.0])

    t1Param =   Param(-1.0, 2)
    AddIsotropicBonds!(t1Param, UC, 1.0, SpinVec[4], "t1")
    UParam  =   Param(U, 4)
    AddIsotropicBonds!(UParam, UC, 0.0, Hubbard, "Hubbard Interaction")
    CreateUnitCell!(UC, [t1Param])

    bz      =   BZ([kSize, kSize])
    FillBZ!(bz, UC)

    t_s     =   Param(-1.0, 2)
    AddIsotropicBonds!(t_s, UC, 1.0, SpinVec[4], "s Hopping")
    Neel_1  =   Param(1.0, 2)
    AddAnisotropicBond!(Neel_1, UC, 1, 1, [0, 0], SpinVec[3], 0.0, "Neel order_1")
    Neel_2  =   Param(1.0, 2)
    AddAnisotropicBond!(Neel_2, UC, 2, 2, [0, 0], SpinVec[3], 0.0, "Neel order_2")

    H       =   Hamiltonian(UC, bz)
    DiagonalizeHamiltonian!(H)
    M       =   Model(UC, bz, H; T = 1e-5, filling = filling, stat = -1)
    SolveModel!(M)

    return TBMFTModel(M, [t_s, Neel_1, Neel_2], [UParam], IntraQuarticToHopping)
end

##### 1D chain, two-site unit cell, on-site Hubbard U, order parameters [s-hopping, Sz_1, Sz_2]
function ChainHubbard(; U::Float64 = 6.0, kSize::Int64 = 31) :: TBMFTModel
    UC      =   UnitCell([[2.0]], 2, 2)
    AddBasisSite!(UC, [0.0])
    AddBasisSite!(UC, [1.0])

    t1Param =   Param(-1.0, 2)
    AddIsotropicBonds!(t1Param, UC, 1.0, SpinVec[4], "t1")
    UParam  =   Param(U, 4)
    AddIsotropicBonds!(UParam, UC, 0.0, Hubbard, "Hubbard Interaction")
    CreateUnitCell!(UC, [t1Param])

    bz      =   BZ([kSize])
    FillBZ!(bz, UC)

    t_s     =   Param(-1.0, 2)
    AddIsotropicBonds!(t_s, UC, 1.0, SpinVec[4], "s Hopping")
    Neel_1  =   Param(1.0, 2)
    AddAnisotropicBond!(Neel_1, UC, 1, 1, [0], SpinVec[3], 0.0, "Neel order_1")
    Neel_2  =   Param(1.0, 2)
    AddAnisotropicBond!(Neel_2, UC, 2, 2, [0], SpinVec[3], 0.0, "Neel order_2")

    H       =   Hamiltonian(UC, bz)
    DiagonalizeHamiltonian!(H)
    M       =   Model(UC, bz, H; T = 1e-5, filling = 0.5, stat = -1)
    SolveModel!(M)

    return TBMFTModel(M, [t_s, Neel_1, Neel_2], [UParam], IntraQuarticToHopping)
end

##### Renormalized square-lattice t-J model with d-wave pairing (reduced Sample/RenormalizedSquaretJ.jl), order parameters [s-hopping, d-wave pairing]
function SquaretJ(; filling::Float64 = 0.4, kSize::Int64 = 15) :: BdGMFTModel
    HoppingUC   =   UnitCell([[1.0, 1.0], [1.0, -1.0]], 2, 2)
    PairingUC   =   UnitCell([[1.0, 1.0], [1.0, -1.0]], 2, 2)
    for uc in (HoppingUC, PairingUC)
        AddBasisSite!(uc, [0.0, 0.0])
        AddBasisSite!(uc, [1.0, 0.0])
    end

    t, J        =   1.0, 0.2
    delta       =   1 - 2 * filling
    t1Param     =   Param(-t * 2 * delta / (1 + delta), 2)
    AddIsotropicBonds!(t1Param, HoppingUC, 1.0, SpinVec[4], "t1")
    JParam      =   Param(J * 4 / ((1 + delta)^2), 4)
    AddIsotropicBonds!(JParam, HoppingUC, 1.0, SpinToPartonCoupling(Matrix{Float64}(I, 3, 3), 1 // 2), "Heisenberg Interaction")
    CreateUnitCell!(HoppingUC, [t1Param])

    bz          =   BZ([kSize, kSize])
    FillBZ!(bz, HoppingUC)

    t_s         =   Param(1.0, 2)
    AddIsotropicBonds!(t_s, HoppingUC, 1.0, SpinVec[4], "s Hopping")
    p_d         =   Param(1.0, 2)
    AddAnisotropicBond!(p_d, PairingUC, 1, 2, [ 0,  0],  SpinVec[2], 1.0, "d_x^2-y^2 Pairing")
    AddAnisotropicBond!(p_d, PairingUC, 1, 2, [-1, -1],  SpinVec[2], 1.0, "d_x^2-y^2 Pairing")
    AddAnisotropicBond!(p_d, PairingUC, 1, 2, [ 0, -1], -SpinVec[2], 1.0, "d_x^2-y^2 Pairing")
    AddAnisotropicBond!(p_d, PairingUC, 1, 2, [-1,  0], -SpinVec[2], 1.0, "d_x^2-y^2 Pairing")

    H           =   Hamiltonian(HoppingUC, PairingUC, bz)
    DiagonalizeHamiltonian!(H)
    M           =   BdGModel(HoppingUC, PairingUC, bz, H; T = 0.001, filling = filling, stat = -1)
    SolveModel!(M)

    return BdGMFTModel(M, [t_s], [p_d], [JParam], InterQuarticToHopping, InterQuarticToPairing)
end
