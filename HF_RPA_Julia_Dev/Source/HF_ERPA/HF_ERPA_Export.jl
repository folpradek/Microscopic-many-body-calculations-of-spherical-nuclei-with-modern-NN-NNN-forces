function HF_ERPA_Export(Params::Vector{Any},N_nu::Matrix{Int64},E_RPA::Matrix{Vector{ComplexF64}})
    
    # Export of RPA & TDA spectra ...
    @time HF_ERPA_Spectrum_Export(Params,N_nu,E_RPA)

    # Export of RPA & TDA plot-ready spectra ...
    @time HF_ERPA_Plot_Spectrum_Export(Params,N_nu,E_RPA)

    return
end


function HF_ERPA_Spectrum_Export(Params::Vector{Any},N_nu::Matrix{Int64},E_RPA::Matrix{Vector{ComplexF64}})
    N_max = Params[4]
    N_2max = 2*N_max
    J_max = N_2max + 1
    Orthogon = Params[5]
    Output_File = Params[8]

    if Orthogon == true
        Output_Path_RPA = Output_File * "/Spectra/ERPA_Ortho.dat"
    else
        Output_Path_RPA = Output_File * "/Spectra/ERPA_Spur.dat"
    end

    println("\nPreparing ERPA spectra export ...")

    # RPA spectrum export
    open(Output_Path_RPA, "w") do Write_File
        println(Write_File, "s\tJ\tP\tE")
        @inbounds for J in 0:J_max
            @inbounds for P in 1:2
                N_ph = N_nu[J+1,P]
                @inbounds for nu in 1:N_ph
                    if P == 1
                        println(Write_File, string(J) * "\t" * "+" * "\t" * string(E_RPA[J+1,P][nu]))
                    else
                        println(Write_File, string(J) * "\t" * "-" * "\t" * string(E_RPA[J+1,P][nu]))
                    end
                end
                println(Write_File,"\n")
            end
        end
    end

    println("\nRPA & TDA spectra export finished ...")

    return
end

function HF_ERPA_Plot_Spectrum_Export(Params::Vector{Any},N_nu::Matrix{Int64},E_RPA::Matrix{Vector{ComplexF64}})
    N_max = Params[4]
    N_2max = 2*N_max
    J_max = N_2max + 1
    Orthogon = Params[5]
    Output_File = Params[8]

    if Orthogon == true
        Output_Path_RPA = Output_File * "/Spectra/ERPA_Plot_Ortho.dat"
    else
        Output_Path_RPA = Output_File * "/Spectra/ERPA_Plot_Spur.dat"
    end

    # Spectra label gap parameter ...
    Delta = 1.0

    println("\nPreparing plot-ready export of RPAsolutions ...")

    N_ph = Int64(sum(N_nu))

    RPA_Solution = Matrix{Float64}(undef,4,N_ph+1)

    # 0+ ground state addition ...

    RPA_Solution[1,1] = 0.0
    RPA_Solution[2,1] = 1.0
    RPA_Solution[3,1] = 0.0
    RPA_Solution[4,1] = 0.0

    nu_count = 1
    @inbounds for J in 0:J_max
        @inbounds for P in 1:2
            N_qg = N_nu[J+1,P]
            @inbounds for nu in 1:N_qg
                nu_count += 1
                RPA_Solution[1,nu_count] = Float64(J)
                RPA_Solution[2,nu_count] = Float64(P)
                RPA_Solution[3,nu_count] = real(E_RPA[J+1,P][nu])
                RPA_Solution[4,nu_count] = real(E_RPA[J+1,P][nu])
            end
        end
    end

    E_RPA_full = @views RPA_Solution[3,:]

    Sort_RPA = sortperm(E_RPA_full, by = x -> real(x))

    RPA_Solution .= @views RPA_Solution[:,Sort_RPA]

    # Adjust the JP label positions ...
    @inbounds for nu in 2:(N_ph+1)

        e_RPA_1 = RPA_Solution[4,nu-1]
        e_RPA_2 = RPA_Solution[4,nu]
        if e_RPA_2 < e_RPA_1
            RPA_Solution[4,nu] = e_RPA_1 + 1.0
        elseif abs(e_RPA_2 - e_RPA_1) < Delta
            RPA_Solution[4,nu] = e_RPA_1 + 1.0
        end
    end

    # RPA spectrum plot data export ...
    open(Output_Path_RPA, "w") do Write_File
        println(Write_File, "J\tP\tE\tm")
        @inbounds for nu in 1:N_ph
            J = Int64(round(RPA_Solution[1,nu]))
            P = "P"
            if abs(RPA_Solution[2,nu] - 1.0) < 1e-3
                P = "+"
            else
                P = "-"
            end
            E = RPA_Solution[3,nu]
            E_m = RPA_Solution[4,nu]
            Row = string(string(J) * "\t" * P * "\t" * string(round(E,sigdigits = 5)) * "\t" * string(round(E_m,sigdigits = 5)))
            println(Write_File, Row)
        end
    end

    println("\nRPA & TDA plot-ready solutions exported ...")

    return
end