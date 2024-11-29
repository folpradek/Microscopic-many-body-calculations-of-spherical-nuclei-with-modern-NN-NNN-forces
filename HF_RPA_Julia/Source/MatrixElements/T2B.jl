function T2B(Orb::Vector{NOrb},a::Int64,b::Int64,c::Int64,d::Int64,J::Int64)
    Grad1 = ReducedGradient(Orb,a,c)
    Grad2 = ReducedGradient(Orb,b,d)
    j_a = Orb[a].j
    j_b = Orb[b].j
    j_c = Orb[c].j
    j_d = Orb[d].j
    TNN = (-1)^(div(j_b + j_c,2) + J) * Grad1 * Grad2 * f6j(j_a, j_b, 2*J, j_d, j_c, 2)
    return TNN
end

function ReducedGradient(Orb::Vector{NOrb},a::Int64,b::Int64)
    n_a = Orb[a].n
    l_a = Orb[a].l
    j_a = Orb[a].j
    n_b = Orb[b].n
    l_b = Orb[b].l
    j_b = Orb[b].j
    Amp = f6j(j_a,j_b,2,2*l_b,2*l_a,1) * sqrt((j_a + 1) * (j_b + 1)) * (-1)^(l_a + div(j_b + 1,2))
    Grad = Amp * (sqrt((l_b+1)*(n_b+l_b+3/2)) * KroneckerDelta(l_a,l_b+1)*KroneckerDelta(n_a,n_b) +
                  sqrt((l_b+1)*n_b) * KroneckerDelta(l_a,l_b+1)*KroneckerDelta(n_a,n_b-1) +
                  sqrt(l_b*(n_b+l_b+1/2)) * KroneckerDelta(l_a,l_b-1)*KroneckerDelta(n_a,n_b) +
                  sqrt(l_b*(n_b+1)) * KroneckerDelta(l_a,l_b-1)*KroneckerDelta(n_a,n_b+1))
    #=
    if n_a == n_b && l_a == (l_b + 1)
        Grad = (-1)^(l_a + div(j_b + 1,2)) * sqrt((j_a + 1) * (j_b + 1)) *
            f6j(j_a, j_b, 2, 2*l_b, 2*l_a, 1) * sqrt((l_b + 1) * (n_b + l_b + 3/2))
    elseif n_a == (n_b - 1) && l_a == (l_b + 1)
        Grad = (-1)^(l_a + div(j_b + 1,2)) * sqrt((j_a + 1) * (j_b + 1)) *
            f6j(j_a, j_b, 2, 2*l_b, 2*l_a, 1) * sqrt((l_b + 1) * n_b)
    elseif n_a == n_b && l_a == (l_b - 1)
        Grad = (-1)^(l_a + div(j_b + 1,2)) * sqrt((j_a + 1) * (j_b + 1)) *
            f6j(j_a, j_b, 2, 2*l_b, 2*l_a, 1) * sqrt(l_b * (n_b + l_b + 1/2))
    elseif n_a == (n_b + 1) && l_a == (l_b - 1)
        Grad = (-1)^(l_a + div(j_b + 1,2)) * sqrt((j_a + 1) * (j_b + 1)) *
            f6j(j_a, j_b, 2, 2*l_b, 2*l_a, 1) * sqrt(l_b * (n_b + 1))
    end
    =#
    return Grad
end