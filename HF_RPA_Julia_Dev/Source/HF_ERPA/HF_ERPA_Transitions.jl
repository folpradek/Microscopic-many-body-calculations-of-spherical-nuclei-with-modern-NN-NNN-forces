#=
struct ReducedMultipole
    E0::pnCVector
    E1::pnCVector
    E2::pnCVector
    E3::pnCVector
end

struct Transition
    ph::Vector{Float64}
    is::Vector{Float64}
    iv::Vector{Float64}
end

struct ReducedTransition
    E0::Transition
    E1::Transition
    E2::Transition
    E3::Transition
end
=#

function HF_ERPA_rM(Params::Vector{Any},N_nu::Matrix{Int64},Orb_Phonon::Matrix{Vector{Int64}},Phonon::Vector{PhState},Particle::pnVec,Hole::pnVec,X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}},TrOp::TranOper,pRho::Matrix{Float64},nRho::Matrix{Float64})
    # Read calculation parameters ...
    Orthogon = Params[5]

    # E0
    J = 0
    P = 1
    N_ph = N_nu[J+1,P]
    prME0 = zeros(ComplexF64,N_ph)
    nrME0 = zeros(ComplexF64,N_ph)
    for nu in 1:N_ph
        pME0Sum = ComplexF64(0.0)
        nME0Sum = ComplexF64(0.0)
        for Ind_ph in 1:N_ph
            ph = Orb_Phonon[J+1,P][Ind_ph]
            p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(ph,Phonon,Particle,Hole)

            if t_ph == -1
                ME0 = TrOp.E0.p[a_p,a_h] * (X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(pRho[a_h,a_h] - pRho[a_p,a_p])
                pME0Sum += ME0
            end

            if t_ph == 1
                ME0 = TrOp.E0.n[a_p,a_h] * (X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(nRho[a_h,a_h] - nRho[a_p,a_p])
                nME0Sum += ME0
            end

        end
        prME0[nu] = pME0Sum
        nrME0[nu] = nME0Sum
    end

    # E1
    J = 1
    P = 2
    N_ph = N_nu[J+1,P]
    prME1 = zeros(ComplexF64,N_ph)
    nrME1 = zeros(ComplexF64,N_ph)
    for nu in 1:N_ph
        if (Orthogon == true && nu != 1) || (Orthogon == false)
            pME1Sum = ComplexF64(0.0)
            nME1Sum = ComplexF64(0.0)
            for Ind_ph in 1:N_ph
                ph = Orb_Phonon[J+1,P][Ind_ph]
                p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(ph,Phonon,Particle,Hole)

                if t_ph == -1
                    ME1 = TrOp.E1.p[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(pRho[a_h,a_h] - pRho[a_p,a_p])
                    pME1Sum += ME1
                end

                if t_ph == 1
                    ME1 = TrOp.E1.n[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(nRho[a_h,a_h] - nRho[a_p,a_p])
                    nME1Sum += ME1
                end

            end
            prME1[nu] = pME1Sum
            nrME1[nu] = nME1Sum
        end
    end

    # E2
    J = 2
    P = 1
    N_ph = N_nu[J+1,P]
    prME2 = zeros(ComplexF64,N_ph)
    nrME2 = zeros(ComplexF64,N_ph)
    for nu in 1:N_ph
        pME2Sum = ComplexF64(0.0)
        nME2Sum = ComplexF64(0.0)
        for Ind_ph in 1:N_ph
            ph = Orb_Phonon[J+1,P][Ind_ph]
            p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(ph,Phonon,Particle,Hole)

            if t_ph == -1
                ME2 = TrOp.E2.p[a_p,a_h] * (X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(pRho[a_h,a_h] - pRho[a_p,a_p])
                pME2Sum += ME2
            end

            if t_ph == 1
                ME2 = TrOp.E2.n[a_p,a_h] * (X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(nRho[a_h,a_h] - nRho[a_p,a_p])
                nME2Sum += ME2
            end

        end
        prME2[nu] = pME2Sum
        nrME2[nu] = nME2Sum
    end

    # E3
    J = 3
    P = 2
    N_ph = N_nu[J+1,P]
    prME3 = zeros(ComplexF64,N_ph)
    nrME3 = zeros(ComplexF64,N_ph)
    for nu in 1:N_ph
        pME3Sum = ComplexF64(0.0)
        nME3Sum = ComplexF64(0.0)
        for Ind_ph in 1:N_ph
            ph = Orb_Phonon[J+1,P][Ind_ph]
            p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(ph,Phonon,Particle,Hole)

            if t_ph == -1
                ME3 = TrOp.E3.p[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(pRho[a_h,a_h] - pRho[a_p,a_p])
                pME3Sum += ME3
            end

            if t_ph == 1
                ME3 = TrOp.E3.n[a_p,a_h] * (-1.0 * X_RPA[J+1,P][Ind_ph,nu] + Y_RPA[J+1,P][Ind_ph,nu]) * sqrt(nRho[a_h,a_h] - nRho[a_p,a_p])
                nME3Sum += ME3
            end

        end
        prME3[nu] = pME3Sum
        nrME3[nu] = nME3Sum
    end

    rM0 = pnCVector(prME0, nrME0)
    rM1 = pnCVector(prME1, nrME1)
    rM2 = pnCVector(prME2, nrME2)
    rM3 = pnCVector(prME3, nrME3)

    rM = ReducedMultipole(rM0,rM1,rM2,rM3)


    return rM
end

function HF_ERPA_rB(Params::Vector{Any},N_nu::Matrix{Int64},rM::ReducedMultipole)
    Orthogon = Params[5]

    # E0
    J = 0
    P = 1
    N_ph = N_nu[J+1,P]

    rB_phE0 = zeros(Float64,N_ph)
    rB_isE0 = zeros(Float64,N_ph)
    rB_ivE0 = zeros(Float64,N_ph)

    @inbounds for nu in 1:N_ph
        rB_phE0[nu] = abs(rM.E0.p[nu])^2
        rB_isE0[nu] = 0.25 * abs(rM.E0.p[nu] + rM.E0.n[nu])^2
        rB_ivE0[nu] = 0.25 * abs(rM.E0.p[nu] - rM.E0.n[nu])^2
    end

    # E1
    J = 1
    P = 2
    N_ph = N_nu[J+1,P]

    rB_phE1 = zeros(Float64,N_ph)
    rB_isE1 = zeros(Float64,N_ph)
    rB_ivE1 = zeros(Float64,N_ph)

    if Orthogon == true
        @inbounds for nu in 1:N_ph
            rB_phE1[nu] = abs(rM.E1.p[nu])^2
            rB_isE1[nu] = 0.25 * abs(rM.E1.p[nu] + rM.E1.n[nu])^2
            rB_ivE1[nu] = 0.25 * abs(rM.E1.p[nu] - rM.E1.n[nu])^2
        end
    else
        A = Params[1]
        Z = Params[2]
        e_p = Float64(A - Z) / Float64(A)
        e_n = Float64(Z) / Float64(A)
        @inbounds for nu in 1:N_ph
            rB_phE1[nu] = abs(rM.E1.p[nu])^2
            rB_isE1[nu] = 0.25 * abs(rM.E1.p[nu] + rM.E1.n[nu])^2
            rB_ivE1[nu] = abs(e_p * rM.E1.p[nu] - e_n * rM.E1.n[nu])^2
        end
    end

    # E2
    J = 2
    P = 1
    N_ph = N_nu[J+1,P]

    rB_phE2 = zeros(Float64,N_ph)
    rB_isE2 = zeros(Float64,N_ph)
    rB_ivE2 = zeros(Float64,N_ph)

    @inbounds for nu in 1:N_ph
        rB_phE2[nu] = abs(rM.E2.p[nu])^2
        rB_isE2[nu] = 0.25 * abs(rM.E2.p[nu] + rM.E2.n[nu])^2
        rB_ivE2[nu] = 0.25 * abs(rM.E2.p[nu] - rM.E2.n[nu])^2
    end

    # E3
    J = 3
    P = 2
    N_ph = N_nu[J+1,P]

    rB_phE3 = zeros(Float64,N_ph)
    rB_isE3 = zeros(Float64,N_ph)
    rB_ivE3 = zeros(Float64,N_ph)

    @inbounds for nu in 1:N_ph
        rB_phE3[nu] = abs(rM.E3.p[nu])^2
        rB_isE3[nu] = 0.25 * abs(rM.E3.p[nu] + rM.E3.n[nu])^2
        rB_ivE3[nu] = 0.25 * abs(rM.E3.p[nu] - rM.E3.n[nu])^2
    end

    rB_E0 = Transition(rB_phE0,rB_isE0,rB_ivE0)
    rB_E1 = Transition(rB_phE1,rB_isE1,rB_ivE1)
    rB_E2 = Transition(rB_phE2,rB_isE2,rB_ivE2)
    rB_E3 = Transition(rB_phE3,rB_isE3,rB_ivE3)

    rB_ph = ReducedTransition(rB_E0,rB_E1,rB_E2,rB_E3)

    return rB_ph

end