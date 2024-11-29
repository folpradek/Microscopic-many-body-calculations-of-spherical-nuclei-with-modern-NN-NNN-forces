function HF_Radial_Potential(Params::Vector{Any},Orb::Vector{NOrb},pU::Matrix{Float64},nU::Matrix{Float64},Vp::Matrix{Float64},Vn::Matrix{Float64})
    # Read calculation parameters ...
    HbarOmega = Params[1]
    A = Params[5]
    N_max = Params[7]
    Output_File = Params[13]

    a_max = div((N_max+1)*(N_max+2),2)

    println("\nPreparing radial mean-field potentials...")
    HbarC = 197.326980

    r1 = 0.0
    r2 = 2.5 * 1.2 * A^(1/3)
    N_Sampling = 10000
    r_grid = range(r1, stop=r2, length=N_Sampling)
    r_grid = collect(r_grid)

    pRho_rad = @views readdlm("IO/" * Output_File * "/Densities/HF_Radial_Densities.dat", Float64)[:,2]
    nRho_rad = @views readdlm("IO/" * Output_File * "/Densities/HF_Radial_Densities.dat", Float64)[:,3]
    pNorm = Integrate_Trap(r_grid, r_grid.^2 .* pRho_rad)
    nNorm = Integrate_Trap(r_grid, r_grid.^2 .* nRho_rad)

    Vp_rad = zeros(Float64, N_Sampling)
    Vn_rad = zeros(Float64, N_Sampling)

    nu_proton = 0.5 * 939.565346 * HbarOmega / HbarC^2
    nu_neutron = 0.5 * 938.272013 * HbarOmega / HbarC^2

    @inbounds for i in 1:N_Sampling
        r = r_grid[i]
        pSum = 0.0
        nSum = 0.0
        @inbounds for a in 1:a_max
            Rad_proton = 0.0
            Rad_neutron = 0.0
            @inbounds for b in 1:a_max
                n = Orb[b].n
                l = Orb[b].l
                N_proton = sqrt(sqrt(2 * nu_proton^3 / π) * 2^(n + 2*l + 3) * factorial(n) * nu_proton^l / prod(2*n + 2*l + 1:-2:1))
                N_neutron = sqrt(sqrt(2 * nu_neutron^3 / π) * 2^(n + 2*l + 3) * factorial(n) * nu_neutron^l / prod(2*n + 2*l + 1:-2:1))

                Exp_proton = exp(-nu_proton * r^2)
                Exp_neutron = exp(-nu_neutron * r^2)
                
                Rad_proton = Rad_proton + N_proton * r^l * Exp_proton * GeneralizedLaguerre(n, l + 0.5, 2 * nu_proton * r^2) * pU[b,a]
                Rad_neutron = Rad_neutron + N_neutron * r^l * Exp_neutron * GeneralizedLaguerre(n, l + 0.5, 2 * nu_neutron * r^2) * nU[b,a]
            end
            pSum += Vp[a,a] * Orb[a].pO * (Rad_proton)^2 * (Orb[a].j + 1)
            nSum += Vn[a,a] * Orb[a].nO * (Rad_neutron)^2 * (Orb[a].j + 1)
        end
        Vp_rad[i] = pSum / pNorm
        Vn_rad[i] = nSum / nNorm
    end

    Density_File_name = "IO/" * Output_File * "/Densities/HF_RadialPotential.dat"

    open(Density_File_name, "w") do io
        writedlm(io, hcat(r_grid, Vp_rad, Vn_rad), "\t")
    end

    println("\nMean-field potential exported...")
    
    return
end