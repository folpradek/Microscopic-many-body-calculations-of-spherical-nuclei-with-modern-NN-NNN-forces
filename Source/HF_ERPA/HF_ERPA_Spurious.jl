function HF_ERPA_Spur_Ini(N_ph::Int64,Orb_Phonon::Vector{Int64},Phonon::Vector{PhState},Particle::pnSVector,Hole::pnSVector,TrOp::TranOper,pRho::Matrix{Float64},nRho::Matrix{Float64},X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}})
    ph_CMS = Vector{ComplexF64}(undef,N_ph)
    @inbounds for ph in 1:N_ph
        Ind_ph = Orb_Phonon[ph]
        p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(Ind_ph,Phonon,Particle,Hole)
        if t_ph == 1
            Sum = 0.0
            @inbounds for nu in 1:N_ph
                ME_E1 = TrOp.E1.n[a_p,a_h] * sqrt(nRho[a_h,a_h] - nRho[a_p,a_p]) * conj(X_RPA[J_ph+1,P_ph][ph,nu]) * (X_RPA[J_ph+1,P_ph][ph,nu] + Y_RPA[J_ph+1,P_ph][ph,nu])
                Sum += ME_E1
            end
            ph_CMS[ph] = Sum
        elseif t_ph == -1
            Sum = 0.0
            @inbounds for nu in 1:N_ph
                ME_E1 = TrOp.E1.p[a_p,a_h] * sqrt(pRho[a_h,a_h] - pRho[a_p,a_p]) * conj(X_RPA[J_ph+1,P_ph][ph,nu]) * (X_RPA[J_ph+1,P_ph][ph,nu] + Y_RPA[J_ph+1,P_ph][ph,nu])
                Sum += ME_E1
            end
            ph_CMS[ph] = Sum
        end
    end

    @inbounds for ph in 1:N_ph
        ME = ph_CMS[ph]
        if abs(ME) < 1e-6
            ph_CMS[ph] = 0.0
        end
    end

    ph_CMS = ph_CMS ./ norm(ph_CMS)
    return ph_CMS
    
end

function HF_ERPA_Spur_Ortho(N_ph::Int64,Spur_State::Vector{ComplexF64})
    # Orthogonalization done by Householder reflection algorithm
    # For explanation see ...
    #
    # https://blogs.mathworks.com/cleve/2016/07/25/compare-gram-schmidt-and-householder-orthogonalization-algorithms/?s_tid=answers_rc2-1_p4_BOTH
    #
    # "More numerically stable than modified Gramm-Schmidt with same numerical cost"

    Spur_Ind = 0

    I = diagm(ones(ComplexF64,N_ph))

    @inbounds for a in 1:N_ph
        proj = 0.0
        @inbounds for b in 1:N_ph
            ME = Spur_State[b] * I[b,a]
            proj += ME
        end
        @inbounds for b in 1:N_ph
            ME = I[b,a] - proj * Spur_State[b]
            I[b,a] = ME
        end
    end

    @inbounds for a in 1:N_ph
        b = @views I[:,a] ./ norm(I[:,a]) 
        @views I[:,a] = b
    end

    U = zeros(ComplexF64, N_ph, N_ph)

    @views U[:,1] = I[:,1]

    @inbounds for nu in 2:N_ph
        u = @views I[:,nu]
        @inbounds for mu in 1:(nu - 1)
            u -= @views (U[:, mu]' * u) * U[:, mu]
        end
        if abs(norm(u)) < 1e-4
            Spur_Ind = nu
        end
        @views U[:, nu] = u ./ norm(u)
    end
    @views U[:,Spur_Ind] = Spur_State

    # Check orthonormality ... precision in order of 1e-15 or smaller is thanks to the Householder algorithm ... very precise
    #println("Orthogonality check (close to identity?) ...")
    #orthogonality_check = U' * U
    #display(orthogonality_check)

    return U, Spur_Ind
end