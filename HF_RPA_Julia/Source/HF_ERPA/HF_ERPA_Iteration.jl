function HF_ERPA_Iteration(Params::Vector{Any},N_nu::Matrix{Int64},Orb_Phonon::Matrix{Vector{Int64}},Phonon::Vector{PhState},Particle::pnVec,Hole::pnVec,Orb::Vector{NOrb},Orb_NN::NNOrb,VNN::NNInt,TrOp::TranOper)
    # Read parameters ...
    N_max = Params[4]
    a_max = div((N_max+1)*(N_max+2),2)

    # Iteration parameters ...
    Eta = 1.0
    Eps = 1e-6
    Iter = 0
    Iter_max = 100

    # Initialize density matrices ...
    pRho, nRho = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    pRho_0, nRho_0 = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    
    # Initialize iteration parameters ...
    Params_Iter = Any[Params[1],Params[2],Params[3],Params[4],false,Params[6],Params[7],Params[8]]
    
    # Standard RPA calculation - ERPA(0) ...

    # Allocate matrices A & B ...
    @time A, B = HF_ERPA_0_Allocate(Params_Iter,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN,VNN)

    # Solve RPA generalized-eigenvalue problem ...
    @time E_TDA, X_TDA, E_RPA, X_RPA, Y_RPA = HF_ERPA_0_Diagonalize(Params_Iter,A,B,N_nu)

    # ERPA iteration ...
    @time while (Eta > Eps) && (Iter < Iter_max)

        # Construct OBDM ...
        @time pRho, nRho = HF_ERPA_OBDM(Params_Iter,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,X_RPA,Y_RPA)

        # Allocate A & B matrices ...
        @time A, B = HF_ERPA_Allocate(Params_Iter,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN,VNN,pRho,nRho)

        # Solve RPA eqs. ...
        @time E_RPA, X_RPA, Y_RPA = HF_ERPA_I_Diagonalize(Params_Iter,A,B,N_nu)

        Sum = 0.0
        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                ME = (abs(pRho[a,b] - pRho_0[a,b]) + abs(nRho[a,b] - nRho_0[a,b])) / Float64(2*a_max^2)
                Sum += ME
            end
        end

        Eta = Sum
        Iter += 1

        println("\nERPA iteration:   Iteration Number = " * string(Iter) * ",\tEta = " * string(Eta))

        pRho_0, nRho_0 = deepcopy(pRho), deepcopy(nRho)

    end

    println("\nERPA iteration has finished ...")

    # Final calculation with/without Orthogonalization ...

    # Allocate A & B matrices ...
    @time A, B = HF_ERPA_Allocate(Params,N_nu,Orb_Phonon,Phonon,Particle,Hole,Orb,Orb_NN,VNN,pRho,nRho)

    # Solve RPA eqs. ...
    @time E_RPA, X_RPA, Y_RPA = HF_ERPA_Diagonalize(Params,A,B,N_nu,Orb_Phonon,Phonon,Particle,Hole,TrOp,pRho,nRho,X_RPA,Y_RPA)

    return E_RPA, X_RPA, Y_RPA, pRho, nRho
end