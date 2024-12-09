function HF_ERPA_Energy(Params::Vector{Any},N_nu::Matrix{Int64},Orb_Phonon::Matrix{Vector{Int64}},Phonon::Vector{PhState},Particle::pnVec,N_Hole::pnInt,Hole::pnVec,Orb::Vector{NOrb},Orb_NN::NNOrb,VNN::NNInt,X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}},pRho::Matrix{Float64},nRho::Matrix{Float64})
    N_max = Params[4]
    N_2max = 2*N_max
    J_max = N_2max + 1
    Orthogon = Params[5]

    E_0 = 0.0

    # Proton hole state energies ...
    @inbounds for h in 1:N_Hole.p
        E_h = Hole.p[h].E
        E_0 += E_h
    end
    # Neutron hole state energies ...
    @inbounds for h in 1:N_Hole.n
        E_h = Hole.n[h].E
        E_0 += E_h
    end

    E_0_1b = ComplexF64(0.0)
    E_0_2b_ph_pp = ComplexF64(0.0)
    E_0_2b_ph_nn = ComplexF64(0.0)
    E_0_2b_ph_pn = ComplexF64(0.0)
    E_0_2b_pphh_pp = ComplexF64(0.0)
    E_0_2b_pphh_nn = ComplexF64(0.0)


    # 1-body HF contributions ...
    @inbounds for J in 0:J_max
        @inbounds for P in 1:2
            N_ph = N_nu[J+1,P]

            @inbounds for nu in 1:N_ph

                @inbounds for ph in 1:N_ph
                    Ind_ph = Orb_Phonon[J+1,P][ph]
                    p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(Ind_ph,Phonon,Particle,Hole)
                    D = 0.0
                    if t_ph == -1
                        D = (pRho[a_h,a_h] - pRho[a_p,a_p])
                    elseif t_ph == 1
                        D = (nRho[a_h,a_h] - nRho[a_p,a_p])
                    end
                    ME = D * (E_p / Float64(j_p + 1) - E_h / Float64(j_h + 1) ) * abs(Y_RPA[J+1,P][ph,nu])^2 * Float64(2*J + 1)
                    E_0_1b += ME

                    @inbounds for mu in 1:N_ph

                        @inbounds for qg in 1:N_ph
                            Ind_qg = Orb_Phonon[J+1,P][qg]
                            q,g,t_qg,J_qg,P_qg,a_q,l_q,j_q,E_q,a_g,l_g,j_g,E_g = Phonon_Ind(Ind_qg,Phonon,Particle,Hole)
                            D = 0.0
                            if t_ph == -1
                                if t_qg == -1
                                    D = (pRho[a_h,a_h] - pRho[a_p,a_p]) * (pRho[a_g,a_g] - pRho[a_q,a_q])
                                elseif t_qg == 1
                                    D = (pRho[a_h,a_h] - pRho[a_p,a_p]) * (nRho[a_g,a_g] - nRho[a_q,a_q])
                                end
                            elseif t_ph == 1
                                if t_qg == -1
                                    D = (nRho[a_h,a_h] - nRho[a_p,a_p]) * (pRho[a_g,a_g] - pRho[a_q,a_q])
                                elseif t_qg == 1
                                    D = (nRho[a_h,a_h] - nRho[a_p,a_p]) * (nRho[a_g,a_g] - nRho[a_q,a_q])
                                end
                            end
                            ME = 0.5 * D * (E_h / Float64(j_h + 1) - E_p / Float64(j_p + 1)) * Y_RPA[J+1,P][ph,nu] * conj(Y_RPA[J+1,P][ph,mu]) * X_RPA[J+1,P][qg,mu] * conj(X_RPA[J+1,P][qg,nu]) * Float64(2*J + 1)^2
                            E_0_1b += ME
                        end

                    end
                end

            end

        end
    end
    println("\n E_0_1b = " * string(real(E_0_1b)) * " Mev")

    # 2-body Residual interaction ... particle-hole terms...
    @inbounds for J in 0:J_max
        @inbounds for P in 1:2
            N_ph = N_nu[J+1,P]

            @inbounds for ph in 1:N_ph
                Ind_ph = Orb_Phonon[J+1,P][ph]
                p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(Ind_ph,Phonon,Particle,Hole)

                @inbounds for qg in 1:N_ph
                    Ind_qg = Orb_Phonon[J+1,P][qg]
                    q,g,t_qg,J_qg,P_qg,a_q,l_q,j_q,E_q,a_g,l_g,j_g,E_g = Phonon_Ind(Ind_qg,Phonon,Particle,Hole)

                    @inbounds for nu in 1:N_ph

                        #if P_ph == P && P_qg == P && J_ph == J && J_qg == J
                            ME1 = 0.0
                            ME2 = 0.0

                            if t_ph == -1 && abs(j_p - j_q) <= 2*J && 2*J <= (j_p + j_q) && rem(l_p + l_q,2) == rem(l_h + l_g,2) && abs(j_h - j_g) <= 2*J && 2*J <= (j_h + j_g)
                                if t_qg == -1
                                    ME1 = - 0.25 * sqrt(pRho[a_h,a_h] - pRho[a_p,a_p]) * sqrt(pRho[a_g,a_g] - pRho[a_q,a_q]) * V2B(a_p,a_q,a_h,a_g,J,1,VNN.pp,Orb,Orb_NN) *
                                                (Y_RPA[J+1,P][ph,nu] * conj(X_RPA[J+1,P][qg,nu]) + X_RPA[J+1,P][ph,nu] * conj(Y_RPA[J+1,P][qg,nu])) * Float64(2*J + 1)
                                    E_0_2b_ph_pp += ME1
                                elseif t_qg == 1
                                    ME1 = sqrt(pRho[a_h,a_h] - pRho[a_p,a_p]) * sqrt(nRho[a_g,a_g] - nRho[a_q,a_q]) * V2B(a_p,a_q,a_h,a_g,J,0,VNN.pn,Orb,Orb_NN) *
                                                (Y_RPA[J+1,P][ph,nu] * conj(X_RPA[J+1,P][qg,nu]) + X_RPA[J+1,P][ph,nu] * conj(Y_RPA[J+1,P][qg,nu])) * Float64(2*J + 1)
                                    E_0_2b_ph_pn += ME1
                                end
                            elseif t_ph == 1 && abs(j_p - j_q) <= 2*J && 2*J <= (j_p + j_q) && rem(l_p + l_q,2) == rem(l_h + l_g,2) && abs(j_h - j_g) <= 2*J && 2*J <= (j_h + j_g)
                                if t_qg == 1
                                    ME1 =  - 0.25 * sqrt(nRho[a_h,a_h] - nRho[a_p,a_p]) * sqrt(nRho[a_g,a_g] - nRho[a_q,a_q]) * V2B(a_p,a_q,a_h,a_g,J,1,VNN.nn,Orb,Orb_NN) *
                                                (Y_RPA[J+1,P][ph,nu] * conj(X_RPA[J+1,P][qg,nu]) + X_RPA[J+1,P][ph,nu] * conj(Y_RPA[J+1,P][qg,nu])) * Float64(2*J + 1)
                                    E_0_2b_ph_nn += ME1
                                end
                            end

                            if t_ph == -1 && t_qg == 1 && abs(j_p - j_g) <= 2*J && 2*J <= (j_p + j_g) && rem(l_p + l_g,2) == rem(l_h + l_q,2) && abs(j_h - j_q) <= 2*J && 2*J <= (j_h + j_q)

                                ME2 = sqrt(pRho[a_h,a_h] - pRho[a_p,a_p]) * sqrt(nRho[a_g,a_g] - nRho[a_q,a_q]) * V2B(a_p,a_g,a_h,a_q,J,0,VNN.pn,Orb,Orb_NN) *
                                            (Y_RPA[J+1,P][ph,nu] * conj(Y_RPA[J+1,P][qg,nu]) + X_RPA[J+1,P][ph,nu] * conj(X_RPA[J+1,P][qg,nu])) * Float64(2*J + 1)
                                E_0_2b_ph_pn += ME2
                            end

                            #E_0 += ME1 + ME2

                        #end
                    end
                end
            end
        end
    end

    println("\n E_0_2b_ph_pp = " * string(real(E_0_2b_ph_pp)) * " Mev")
    println("\n E_0_2b_ph_nn = " * string(real(E_0_2b_ph_nn)) * " Mev")
    println("\n E_0_2b_ph_pn = " * string(real(E_0_2b_ph_pn)) * " Mev")

    # 2-body Residual interaction ... particle-particle & hole-hole terms...
    @inbounds for J in 0:J_max
        @inbounds for P in 1:2
            N_ph = N_nu[J+1,P]

            @inbounds for ph in 1:N_ph
                Ind_ph = Orb_Phonon[J+1,P][ph]
                p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(Ind_ph,Phonon,Particle,Hole)

                @inbounds for qg in 1:N_ph
                    Ind_qg = Orb_Phonon[J+1,P][qg]
                    q,g,t_qg,J_qg,P_qg,a_q,l_q,j_q,E_q,a_g,l_g,j_g,E_g = Phonon_Ind(Ind_qg,Phonon,Particle,Hole)

                    @inbounds for nu in 1:N_ph
                        @inbounds for mu in 1:N_ph

                            #if P_ph == P && P_qg == P && J_ph == J && J_qg == J
                                ME = 0.0

                                if abs(j_p - j_q) <= 2*J && 2*J <= (j_p + j_q)
                                    if t_ph == -1 && t_qg == -1

                                        ME = 0.25 * (pRho[a_h,a_h] - pRho[a_p,a_p]) * (pRho[a_g,a_g] - pRho[a_q,a_q]) *
                                            Y_RPA[J+1,P][ph,nu] * conj(Y_RPA[J+1,P][ph,mu]) * conj(X_RPA[J+1,P][qg,nu]) * X_RPA[J+1,P][qg,mu] *
                                            V2B(a_p,a_q,a_p,a_q,J,1,VNN.pp,Orb,Orb_NN) * Float64(2*J + 1)^2
                                        E_0_2b_pphh_pp += ME
                                    elseif t_ph == 1 && t_qg == 1

                                        ME = 0.25 * (nRho[a_h,a_h] - nRho[a_p,a_p]) * (nRho[a_g,a_g] - nRho[a_q,a_q]) *
                                            Y_RPA[J+1,P][ph,nu] * conj(Y_RPA[J+1,P][ph,mu]) * conj(X_RPA[J+1,P][qg,nu]) * X_RPA[J+1,P][qg,mu] *
                                            V2B(a_p,a_q,a_p,a_q,J,1,VNN.nn,Orb,Orb_NN) * Float64(2*J + 1)^2
                                        E_0_2b_pphh_nn += ME
                                    end
                                end
                                
                                if abs(j_h - j_g) <= 2*J && 2*J <= (j_h + j_g)
                                    if t_ph == -1 && t_qg == -1

                                        ME = 0.25 * (pRho[a_h,a_h] - pRho[a_p,a_p]) * (pRho[a_g,a_g] - pRho[a_q,a_q]) *
                                            Y_RPA[J+1,P][ph,nu] * conj(Y_RPA[J+1,P][ph,mu]) * conj(X_RPA[J+1,P][qg,nu]) * X_RPA[J+1,P][qg,mu] *
                                            V2B(a_h,a_g,a_h,a_g,J,1,VNN.pp,Orb,Orb_NN) * Float64(2*J + 1)^2
                                        E_0_2b_pphh_pp += ME
                                    elseif t_ph == 1 && t_qg == 1

                                        ME = 0.25 * (nRho[a_h,a_h] - nRho[a_p,a_p]) * (nRho[a_g,a_g] - nRho[a_q,a_q]) *
                                            Y_RPA[J+1,P][ph,nu] * conj(Y_RPA[J+1,P][ph,mu]) * conj(X_RPA[J+1,P][qg,nu]) * X_RPA[J+1,P][qg,mu] *
                                            V2B(a_h,a_g,a_h,a_g,J,1,VNN.nn,Orb,Orb_NN) * Float64(2*J + 1)^2
                                        E_0_2b_pphh_nn += ME
                                    end
                                end
                                
                                #E_0 += ME
                            #end
                        end
                    end
                end
            end
        end
    end

    println("\n E_0_2b_pphh_pp = " * string(real(E_0_2b_pphh_pp)) * " Mev")
    println("\n E_0_2b_pphh_nn = " * string(real(E_0_2b_pphh_nn)) * " Mev")

    E_0 += real(E_0_1b + E_0_2b_ph_pp + E_0_2b_ph_nn + E_0_2b_ph_pn + E_0_2b_pphh_pp + E_0_2b_pphh_nn)
    #E_0 = Float64(real(E_0))

    println("\nERPA ground-state correlation energy E_0^ERPA is ...     E_0 = " * string(round(E_0,digits=6)) * "\tMeV")
    println("\nTo get the total ground-state energy add the mean-field energy & subtract the total mean-field hole states energy (addition of V_0 term from Hamiltonian) ...")
    return E_0
end