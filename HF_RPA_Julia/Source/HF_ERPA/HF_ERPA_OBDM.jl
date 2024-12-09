function HF_ERPA_OBDM(Params::Vector{Any},N_nu::Matrix{Int64},Orb_Phonon::Matrix{Vector{Int64}},Phonon::Vector{PhState},Particle::pnVec,Hole::pnVec,Orb::Vector{NOrb},X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}})
    # Read parameters ...
    N_max = Params[4]
    N_2max = 2*N_max
    J_max = N_2max + 1
    a_max = div((N_max+1)*(N_max+2),2)
    
    # Make HF density matrices ...
    pRho_HF, nRho_HF = zeros(Float64,a_max,a_max), zeros(Float64,a_max,a_max)
    @inbounds for a in 1:a_max
        pRho_HF[a,a] = Orb[a].pO
        nRho_HF[a,a] = Orb[a].nO
    end

    # Initialize new density matrices ...
    pRho, nRho = deepcopy(pRho_HF), deepcopy(nRho_HF)
    pRho_old, nRho_old = deepcopy(pRho_HF), deepcopy(nRho_HF)

    Delta = 1.0
    Eps = 1e-6
    N_Iteration = 100
    Iteration = 0

    println("\nStarting HF ERPA OBDM calculation ...")

    @time while (Delta > Eps) && (N_Iteration > Iteration)

        pRho, nRho = deepcopy(pRho_HF), deepcopy(nRho_HF)
        pRho_t, nRho_t = deepcopy(ComplexF64.(pRho)), deepcopy(ComplexF64.(pRho))

        @inbounds for J in 0:J_max
            @inbounds for P in 1:2

                N_ph = N_nu[J+1,P]

                # Particles
                @inbounds Threads.@threads for ph in 1:N_ph
                    Ind_ph = Orb_Phonon[J+1,P][ph]
                    p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(Ind_ph,Phonon,Particle,Hole)

                    @inbounds for qh in 1:N_ph

                        Ind_qh = Orb_Phonon[J+1,P][qh]
                        q,h_q,t_qh,J_qh,P_qh,a_q,l_q,j_q,E_q,a_h_q,l_h_q,j_h_q,E_h = Phonon_Ind(Ind_qh,Phonon,Particle,Hole)

                        if a_h == a_h_q && j_p == j_q && l_p == l_q

                            pSum = ComplexF64(0.0)
                            nSum = ComplexF64(0.0)

                            @inbounds for nu in 1:N_ph

                                @inbounds for rf in 1:N_ph

                                    Ind_rf = Orb_Phonon[J+1,P][rf]
                                    r,f,t_rf,J_rf,P_rf,a_r,l_r,j_r,E_r,a_f,l_f,j_f,E_f = Phonon_Ind(Ind_rf,Phonon,Particle,Hole)

                                    @inbounds for mu in 1:N_ph

                                        if t_ph == -1 && t_qh == -1 && t_rf == -1

                                            ME = ComplexF64(-0.5 * Float64(2*J + 1) / Float64(j_p + 1) * (pRho_old[a_f,a_f] - pRho_old[a_r,a_r]) * 
                                                    sqrt((pRho_old[a_h,a_h] - pRho_old[a_p,a_p]) * (pRho_old[a_h,a_h] - pRho_old[a_q,a_q]))) *
                                                    X_RPA[J+1,P][mu,rf] * conj(X_RPA[J+1,P][nu,rf]) *  Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][mu,qh])
                                            pSum += ME

                                        end

                                        if t_ph == 1 && t_qh == 1 && t_rf == 1

                                            ME = ComplexF64(-0.5 * Float64(2*J + 1) / Float64(j_p + 1)  * (nRho_old[a_f,a_f] - nRho_old[a_r,a_r]) *
                                                    sqrt((nRho_old[a_h,a_h] - nRho_old[a_p,a_p]) * (nRho_old[a_h,a_h] - nRho_old[a_q,a_q]))) *
                                                    Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][mu,qh]) * X_RPA[J+1,P][mu,rf] * conj(X_RPA[J+1,P][nu,rf])
                                            nSum += ME

                                        end

                                    end
                                end

                                if t_ph == -1 && t_qh == -1
                                    ME = ComplexF64(sqrt((pRho_old[a_h,a_h] - pRho_old[a_p,a_p]) * (pRho_old[a_h,a_h] - pRho_old[a_q,a_q])) *
                                            Float64(2*J + 1) / Float64(j_p + 1)) * Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][nu,qh])
                                    pSum += ME
                                end

                                if t_ph == 1 && t_qh == 1
                                    ME = ComplexF64(sqrt((nRho_old[a_h,a_h] - nRho_old[a_p,a_p]) * (nRho_old[a_h,a_h] - nRho_old[a_q,a_q])) *
                                            Float64(2*J + 1) / Float64(j_p + 1)) * Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][nu,qh])
                                    nSum += ME
                                end

                            end
                            pRho_t[a_p,a_q] += pSum
                            nRho_t[a_p,a_q] += nSum
                        end
                    end
                end

                
                # Holes
                @inbounds Threads.@threads for ph in 1:N_ph
                    Ind_ph = Orb_Phonon[J+1,P][ph]
                    p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(Ind_ph,Phonon,Particle,Hole)

                    @inbounds for pg in 1:N_ph

                        Ind_pg = Orb_Phonon[J+1,P][pg]
                        p_g,g,t_pg,J_pg,P_pg,a_p_g,l_p_g,j_p_g,E_p_g,a_g,l_g,j_g,E_g = Phonon_Ind(Ind_pg,Phonon,Particle,Hole)

                        if a_p == a_p_g && j_h == j_g && l_h == l_g

                            @inbounds for nu in 1:N_ph

                                pSum = ComplexF64(0.0)
                                nSum = ComplexF64(0.0)

                                @inbounds for rf in 1:N_ph

                                    Ind_rf = Orb_Phonon[J+1,P][rf]
                                    r,f,t_rf,J_rf,P_rf,a_r,l_r,j_r,E_r,a_f,l_f,j_f,E_f = Phonon_Ind(Ind_rf,Phonon,Particle,Hole)

                                    @inbounds for mu in 1:N_ph

                                        if t_ph == -1 && t_pg == -1 && t_rf == -1

                                            ME = ComplexF64(0.5 * Float64(2*J + 1) / Float64(j_h + 1) * (pRho_old[a_f,a_f] - pRho_old[a_r,a_r]) *
                                                    sqrt((pRho_old[a_h,a_h] - pRho_old[a_p,a_p]) * (pRho_old[a_g,a_g] - pRho_old[a_p,a_p]))) *
                                                    Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][mu,pg]) * X_RPA[J+1,P][mu,rf] * conj(X_RPA[J+1,P][nu,rf]) 
                                            pSum += ME

                                        end

                                        if t_ph == 1 && t_pg == 1 && t_rf == 1

                                            ME = ComplexF64(0.5 * Float64(2*J + 1) / Float64(j_h + 1)  * (nRho_old[a_f,a_f] - nRho_old[a_r,a_r]) *
                                                    sqrt((nRho_old[a_h,a_h] - nRho_old[a_p,a_p]) * (nRho_old[a_g,a_g] - nRho_old[a_p,a_p]))) *
                                                    Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][mu,pg]) * X_RPA[J+1,P][mu,rf] * conj(X_RPA[J+1,P][nu,rf])
                                            nSum += ME

                                        end

                                    end
                                end

                                if t_ph == -1 && t_pg == -1
                                    ME = ComplexF64(- sqrt((pRho_old[a_h,a_h] - pRho_old[a_p,a_p]) * (pRho_old[a_g,a_g] - pRho_old[a_p,a_p])) *
                                            Float64(2*J + 1) / Float64(j_h + 1)) * Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][nu,pg])
                                    pSum += ME
                                end

                                if t_ph == 1 && t_pg == 1
                                    ME = ComplexF64(- sqrt((nRho_old[a_h,a_h] - nRho_old[a_p,a_p]) * (nRho_old[a_g,a_g] - nRho_old[a_p,a_p])) *
                                            Float64(2*J + 1) / Float64(j_h + 1)) * Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][nu,pg])
                                    nSum += ME
                                end

                                pRho_t[a_h,a_g] += pSum
                                nRho_t[a_h,a_g] += nSum
                            end
                        end
                    end
                end

            end
        end

        Sum = 0.0

        @inbounds for a in 1:a_max
            @inbounds for b in 1:a_max
                pRho[a,b], nRho[a,b] = Float64(real(pRho_t[a,b])), Float64(real(nRho_t[a,b]))
                ME = (abs(pRho[a,b] - pRho_old[a,b]) + abs(nRho[a,b] - nRho_old[a,b])) / Float64(2*a_max^2)
                Sum += ME
            end
        end

        Delta = Sum
        Iteration += 1

        println("\nOBDM iteration Number = " * string(Iteration) * ",\tDelta = " * string(Delta))

        pRho_old, nRho_old = deepcopy(pRho), deepcopy(nRho)

    end

    print("\nERPA oteration of OBDM has terminated ...")

    return pRho, nRho
end