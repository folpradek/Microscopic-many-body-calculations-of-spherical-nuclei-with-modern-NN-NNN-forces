function HF_ERPA_OBDM(Params::Vector{Any},pRho_old::Matrix{Float64},nRho_old::Matrix{Float64},N_nu::Matrix{Int64},Orb_Phonon::Matrix{Vector{Int64}},Phonon::Vector{PhState},N_Particle::pnInt,Particle::pnVec,N_Hole::pnInt,Hole::pnVec,X_RPA::Matrix{Matrix{ComplexF64}},Y_RPA::Matrix{Matrix{ComplexF64}})
    N_max = 3
    N_2max = 2*N_max
    J_max = N_2max + 1
    a_max = div((N_max+1)*(N_max+2),2)

    pRho, nRho = Matrix{Float64}(undef,a_max,a_max), Matrix{Float64}(undef,a_max,a_max)

    pRho_HF, nRho_HF = deepcopy(pRho_old), deepcopy(nRho_old)

    Delta = 1.0
    Eps = 1e-8
    N_Iteration = 100
    Iteration = 0

    println("\nStarting HF ERPA OBDM calculation ...")

    display(pRho_old)

    while (Delta > Eps) && (N_Iteration > Iteration)

        pRho, nRho = deepcopy(pRho_HF), deepcopy(nRho_HF)

        for J in 0:J_max
            for P in 1:2

                N_ph = N_nu[J+1,P]

                # Particles
                for nu in 1:N_ph
                    for ph in 1:N_ph

                        Ind_ph = Orb_Phonon[J+1,P][ph]
                        p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(Ind_ph,Phonon,Particle,Hole)

                        for qh in 1:N_ph

                            Ind_qh = Orb_Phonon[J+1,P][qh]
                            q,h_q,t_qh,J_qh,P_qh,a_q,l_q,j_q,E_q,a_h_q,l_h_q,j_h_q,E_h = Phonon_Ind(Ind_qh,Phonon,Particle,Hole)

                            if a_h == a_h_q && j_p == j_q && l_p == l_q

                                pSum = 0.0
                                nSum = 0.0

                                for mu in 1:N_ph

                                    for rf in 1:N_ph

                                        Ind_rf = Orb_Phonon[J+1,P][rf]
                                        r,f,t_rf,J_rf,P_rf,a_r,l_r,j_r,E_r,a_f,l_f,j_f,E_f = Phonon_Ind(Ind_rf,Phonon,Particle,Hole)

                                        if t_ph == -1 && t_qh == -1 && t_rf == -1

                                            ME = sqrt((pRho_old[a_h,a_h] - pRho_old[a_p,a_p]) * (pRho_old[a_h,a_h] - pRho_old[a_q,a_q])) *
                                                    Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][mu,qh]) * (pRho_old[a_f,a_f] - pRho_old[a_r,a_r]) *
                                                    X_RPA[J+1,P][mu,rf] * conj(X_RPA[J+1,P][nu,rf]) * (-0.5) * Float64(2*J + 1) / Float64(j_p + 1)
                                            pSum += ME

                                        end

                                        if t_ph == 1 && t_qh == 1 && t_rf == 1

                                            ME = sqrt((nRho_old[a_h,a_h] - nRho_old[a_p,a_p]) * (nRho_old[a_h,a_h] - nRho_old[a_q,a_q])) *
                                                    Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][mu,qh]) * (nRho_old[a_f,a_f] - nRho_old[a_r,a_r]) *
                                                    X_RPA[J+1,P][mu,rf] * conj(X_RPA[J+1,P][nu,rf]) * (-0.5) * Float64(2*J + 1) / Float64(j_p + 1)
                                            nSum += ME

                                        end

                                    end
                                end

                                if t_ph == -1 && t_qh == -1
                                    ME = sqrt((pRho_old[a_h,a_h] - pRho_old[a_p,a_p]) * (pRho_old[a_h,a_h] - pRho_old[a_q,a_q])) *
                                            Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][nu,qh]) * Float64(2*J + 1) / Float64(j_p + 1)
                                    pSum += ME
                                    pRho[a_p,a_q] += pSum
                                end

                                if t_ph == 1 && t_qh == 1
                                    ME = sqrt((nRho_old[a_h,a_h] - nRho_old[a_p,a_p]) * (nRho_old[a_h,a_h] - nRho_old[a_q,a_q])) *
                                            Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][nu,qh]) * Float64(2*J + 1) / Float64(j_p + 1)
                                    nSum += ME
                                    nRho[a_p,a_q] += nSum
                                end

                            end
                        end
                    end
                end

                
                # Holes
                for nu in 1:N_ph
                    for ph in 1:N_ph

                        Ind_ph = Orb_Phonon[J+1,P][ph]
                        p,h,t_ph,J_ph,P_ph,a_p,l_p,j_p,E_p,a_h,l_h,j_h,E_h = Phonon_Ind(Ind_ph,Phonon,Particle,Hole)

                        for pg in 1:N_ph

                            Ind_pg = Orb_Phonon[J+1,P][pg]
                            p_g,g,t_pg,J_pg,P_pg,a_p_g,l_p_g,j_p_g,E_p_g,a_g,l_g,j_g,E_g = Phonon_Ind(Ind_pg,Phonon,Particle,Hole)

                            if a_p == a_p_g && j_h == j_g && l_h == l_g

                                pSum = 0.0
                                nSum = 0.0

                                for mu in 1:N_ph

                                    for rf in 1:N_ph

                                        Ind_rf = Orb_Phonon[J+1,P][rf]
                                        r,f,t_rf,J_rf,P_rf,a_r,l_r,j_r,E_r,a_f,l_f,j_f,E_f = Phonon_Ind(Ind_rf,Phonon,Particle,Hole)

                                        if t_ph == -1 && t_pg == -1 && t_rf == -1

                                            ME = sqrt((pRho_old[a_h,a_h] - pRho_old[a_p,a_p]) * (pRho_old[a_g,a_g] - pRho_old[a_p,a_p])) *
                                                    Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][mu,pg]) * (pRho_old[a_f,a_f] - pRho_old[a_r,a_r]) *
                                                    X_RPA[J+1,P][mu,rf] * conj(X_RPA[J+1,P][nu,rf]) * 0.5 * Float64(2*J + 1) / Float64(j_h + 1)
                                            pSum += ME

                                        end

                                        if t_ph == 1 && t_pg == 1 && t_rf == 1

                                            ME = sqrt((nRho_old[a_h,a_h] - nRho_old[a_p,a_p]) * (nRho_old[a_g,a_g] - nRho_old[a_p,a_p])) *
                                                    Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][mu,pg]) * (nRho_old[a_f,a_f] - nRho_old[a_r,a_r]) *
                                                    X_RPA[J+1,P][mu,rf] * conj(X_RPA[J+1,P][nu,rf]) * 0.5 * Float64(2*J + 1) / Float64(j_h + 1)
                                            nSum += ME

                                        end

                                    end
                                end

                                if t_ph == -1 && t_pg == -1
                                    ME = - sqrt((pRho_old[a_h,a_h] - pRho_old[a_p,a_p]) * (pRho_old[a_g,a_g] - pRho_old[a_p,a_p])) *
                                            Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][nu,pg]) * Float64(2*J + 1) / Float64(j_h + 1)
                                    pSum += ME
                                    pRho[a_h,a_g] += pSum
                                end

                                if t_ph == 1 && t_pg == 1
                                    ME = - sqrt((nRho_old[a_h,a_h] - nRho_old[a_p,a_p]) * (nRho_old[a_g,a_g] - nRho_old[a_p,a_p])) *
                                            Y_RPA[J+1,P][nu,ph] * conj(Y_RPA[J+1,P][nu,pg]) * Float64(2*J + 1) / Float64(j_h + 1)
                                    nSum += ME
                                    nRho[a_h,a_g] += nSum
                                end

                            end
                        end
                    end
                end

            end
        end

        Sum = 0.0

        for a in 1:a_max
            for b in 1:a_max

                ME = (abs(pRho[a,b] - pRho_old[a,b]) + abs(nRho[a,b] - nRho_old[a,b])) / (2*a_max)
                Sum += ME
            
            end
        end
        Delta = Sum
        Iteration += 1

        println("\nIteration Number = " * string(Iteration) * ",\tDelta = " * string(Delta))

        pRho_old, nRho_old = deepcopy(pRho), deepcopy(nRho)

    end

    print("Iteration of ERPA OBDM has terminated ...")

    display(pRho)

    return pRho, nRho
end