function HF_ERPA(Params::Vector{Any},Orb::Vector{NOrb},Orb_NN::NNOrb,VNN::NNInt)
    # Calculation parameters review ...
    # A, Z, HbarOmega, N_max = Params[1], Params[2], Params[3], Params[4]
    # N_2max = 2*N_max
    # J_max = N_2max + 1
    # Orthogonalization = Params[5]
    # Input_File = Params[8]

    # Prepare Particle & Hole orbitals ...
    @time N_Particle, Particle, N_Hole, Hole = Make_ParticleHole_Orbitals(Params[4],Params[8],Orb)

    # Prepare 1p-1h phonon states ...
    @time N_Phonon, Phonon = Make_Phonon_Orbitals(N_Particle,Particle,N_Hole,Hole)

    # Initialize transition operators ...
    @time TrOp = Transition_Operators_Ini(Params[3],Params[4],Orb)

    @time TrOp = Transition_Operators_Transform(Params[4],Params[8],TrOp,Orb)

    # Count & pre-index all phonon states in JP subspaces ...
    @time N_nu, Orb_Phonon = HF_RPA_Phonon_Count(Params,N_Phonon,Phonon,Particle,Hole)

    # RPA standard iteration - ERPA(0)
    # Allocate matrices A & B ...
    @time A, B = HF_RPA_Allocate(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN,VNN)

    Params[5] = false

    # Solve TDA eigenvalue problem and RPA generalized-eigenvalue problem ...
    @time E_TDA, X_TDA, E_RPA, X_RPA, Y_RPA = HF_RPA_Diagonalize(Params,A,B,N_nu,Orb_Phonon,Phonon,Particle,Hole,TrOp)

    a_max = length(Orb)
    pRho_HF, nRho_HF = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)

    for a in 1:a_max
        pRho_HF[a,a] = Orb[a].pO
        nRho_HF[a,a] = Orb[a].nO
    end

    Params[5] = true

    # Construct OBDM ...
    @time pRho, nRho = HF_ERPA_OBDM(Params,pRho_HF,nRho_HF,N_nu,Orb_Phonon,Phonon,N_Particle,Particle,N_Hole,Hole,X_RPA,Y_RPA)

    # Diagonalize density operator ...

    # At present time it is diagonal ... ??? ... sth. went wrong?

    # Transform Hamiltonian - 1-body part & 2-body part to the new basis - where pRho, nRho are diagonal ...

    # Not need right now ...

    # Transform transition operators ...

    # Not need right now ...

    # Allocate A & B matrices ...

    @time A, B = HF_ERPA_Allocate(Params,pRho,nRho,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN,VNN)

    # Solve RPA eqs. ...
    @time E_RPA, X_RPA, Y_RPA = HF_ERPA_Diagonalize(Params,pRho,nRho,A,B,N_nu,Orb_Phonon,Phonon,Particle,Hole,TrOp)

    # ERPA solutions export ...
    @time HF_ERPA_Export(Params,N_nu,E_RPA)

    return
end