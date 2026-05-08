struct ThetaQuad
    nt::Int
    thet1::Vector{Float64}
    thet2::Vector{Float64}
    Sθ₁::Vector{Float64}
    Cθ₁::Vector{Float64}
    Sθ₂::Vector{Float64}
    Cθ₂::Vector{Float64}
    fwmt1::Vector{Float64}
    fwmt2::Vector{Float64}
    fwmt3::Vector{Float64}
    r₁::Matrix{Float64}
    r₂::Matrix{Float64}
    r₃::Matrix{Float64}
    rs₁::Matrix{Float64}
    rs₂::Matrix{Float64}
    rs₃::Matrix{Float64}
end

function makeThetaQuad(nt::Int, nr::Int, q::Int, s::Float64, wz1, wz2)
    zt  = Vector{Float64}(undef, nt)
    z2t = Vector{Float64}(undef, nt)

    # Fejér nodes and z2t
    @inbounds for i in 1:nt
        c₁ = π * (2 * i - 1) / (2 * nt)
        zt[i]  = cos(c₁)
        z2t[i] = sin(c₁ / 2)^2
    end

    # angular Fejér weights
    fwt = getF1W(nt)

    # polynomial change of variables
    thet1 = Vector{Float64}(undef, nt)
    thet2 = Vector{Float64}(undef, nt)
    wfunc!(thet1, q, z2t; α=-1.0, β=π/4)
    wfunc!(thet2, q, zt;  α=1.0,  β=π/8)

    Sθ₁ = similar(thet1); Cθ₁ = similar(thet1); Tθ₁ = similar(thet1)
    Sθ₂ = similar(thet1); Cθ₂ = similar(thet1); Tθ₂ = similar(thet1)

    @inbounds for i in 1:nt
        si₁, co₁ = sincos(thet1[i])
        si₂, co₂ = sincos(thet2[i])
        Sθ₁[i] = si₁;  Cθ₁[i] = co₁;  Tθ₁[i] = si₁ / co₁
        Sθ₂[i] = si₂;  Cθ₂[i] = co₂;  Tθ₂[i] = si₂ / co₂
    end

    # fwmt* (angular weights)
    fwmt1 = Vector{Float64}(undef, nt)
    fwmt2 = Vector{Float64}(undef, nt)
    fwmt3 = Vector{Float64}(undef, nt)

    tmp2  = Vector{Float64}(undef, nt)
    tmp3  = Vector{Float64}(undef, nt)
    dwfunc!(tmp2, q, z2t)
    dwfunc!(tmp3, q, zt)

    @inbounds for i in 1:nt
        fwmt1[i] = π * fwt[i] * tmp2[i] * (Cθ₁[i])^(2*(s - 1))   / 16
        fwmt2[i] = π * fwt[i] * tmp2[i] * (2*Cθ₁[i])^(2*(s - 1)) / 8
        fwmt3[i] = π * fwt[i] * tmp3[i] * (2*Cθ₂[i])^(2*(s - 1)) / 8
    end

    #Radial parts of singular integrals
    r₁ = Matrix{Float64}(undef, nt, nr) #sec(thet1) * w(-z2r)'
    r₂ = Matrix{Float64}(undef, nt, nr) #sec(thet1) * w(zr)' / 2
    r₃ = Matrix{Float64}(undef, nt, nr) #sec(thet2) * w(zr)' / 2
    rs₁= Matrix{Float64}(undef, nt, nr) #tan(thet1) * w(-z2r)'
    rs₂= Matrix{Float64}(undef, nt, nr) #tan(thet1) * w(zr)'
    rs₃= Matrix{Float64}(undef, nt, nr) #tan(thet2) * w(zr)'

    @inbounds for i in 1:nt
        c1 = Cθ₁[i]; t1 = Tθ₁[i]
        c2 = Cθ₂[i]; t2 = Tθ₂[i]
        @inbounds for j in 1:nr
            r₁[i,j]  = wz2[j] / c1
            r₂[i,j]  = wz1[j] / (2c1)
            r₃[i,j]  = wz1[j] / (2c2)
            rs₁[i,j] = t1 * wz2[j]
            rs₂[i,j] = t1 * wz1[j]
            rs₃[i,j] = t2 * wz1[j]
        end
    end

    return ThetaQuad(nt, thet1, thet2,
                     Sθ₁, Cθ₁, Sθ₂, Cθ₂,
                     fwmt1, fwmt2, fwmt3,
                     r₁, r₂, r₃, rs₁, rs₂, rs₃)
end

struct ThetaWork
    nt::Int
    d1::Vector{Float64}
    d2::Vector{Float64}
    y1::Matrix{Float64}
    y1tmp::Matrix{Float64}
    y2::Matrix{Float64}
    dfY::Matrix{Float64}
    TY::Array{Float64,3}
    J1::Matrix{Float64}
    Jt::Matrix{Float64}
    DIF::Matrix{Float64}
    Zx::Matrix{Float64}
    Zy::Matrix{Float64}
    DJ::Matrix{Float64} # Jacobian |detJ|
end

function makeThetaWork(nt::Int, nr::Int, N::Int)
    d1   = Vector{Float64}(undef, nt)
    d2   = similar(d1)
    y1   = Matrix{Float64}(undef, nt, nr)
    y1tmp= similar(y1)
    y2   = similar(y1)
    dfY  = Matrix{Float64}(undef, nt, nr)
    TY   = Array{Float64}(undef, nt, nr, N)
    J1   = Matrix{Float64}(undef, nt, nr)
    Jt   = Matrix{Float64}(undef, nt, nr)
    DIF  = Matrix{Float64}(undef, nt, nr)
    Zx   = Matrix{Float64}(undef, nt, nr)
    Zy   = Matrix{Float64}(undef, nt, nr)
    DJ   = Matrix{Float64}(undef, nt, nr)

    return ThetaWork(nt, d1, d2, y1, y1tmp, y2, dfY, TY, J1, Jt, DIF, Zx, Zy, DJ)
end

#For s>=0.5
function precompsH(d::abstractdomain, dp::domprop, s::Float64, p::Int
    ; n::Int=128)::Matrix{Float64}

    #Bookkeeping
    M = d.Npat  # number of patches
    N = dp.N    # number of nodes per patch per Axis
    Np = N * N  # number of nodes per patch

    # zp1 = zp+1, zp2 = xp-1 where zp = cos(pi*(2j-1)/(2N))
    zp  = Vector{Float64}(undef, N)
    zp1 = Vector{Float64}(undef, N)
    zp2 = Vector{Float64}(undef, N)

    for j in 1:N
        zp[j]  = cospi((2 * j - 1) / (2N))
        zp1[j] = 2 * cospi((2 * j - 1) / (4N))^2
        zp2[j] =-2 * sinpi((2 * j - 1) / (4N))^2
    end

    #Parameter q for another polynomial change of variables is theta. 
    # It is usually taken as 4
    q = 4

    #Define nodes in radial and weights
    nr = n     
    fwr = getF1W(nr)

    zr = Vector{Float64}(undef, nr)
    z2r = Vector{Float64}(undef, nr)

    @inbounds for i in 1:nr
        c₁ = π * (2 * i - 1) / (2 * nr)
        zr[i] = cos(c₁)
        z2r[i] = sin(c₁ / 2)^2
    end

    #Precomputing modified weights
    fwmr1 = Vector{Float64}(undef, nr)
    fwmr2 = Vector{Float64}(undef, nr)
    fwmr3 = Vector{Float64}(undef, nr)

    #To store w(zr) and w(-z2r)
    wz1 = Vector{Float64}(undef, nr)
    wz2 = similar(wz1)
    wz3 = similar(wz1)

    wfunc!(wz1, p, zr)
    wfunc!(wz2, p, z2r; α=-1.0)
    wfunc!(wz3, p, zr; α=-1.0)

    qw1func!(fwmr1, p, z2r, s)
    qw1func!(fwmr3, p, zr, s; α=-1.0)

    @. fwmr1 = fwr * fwmr1
    @. fwmr2 = fwr * fwmr3
    @. fwmr3 = fwr * fwmr3

    # Storage of temporary variables inside the for loop
    # These are independent on nt
    I = Matrix{Float64}(undef, N, N)

    yr  = Vector{Float64}(undef, nr)

    TNy = Matrix{Float64}(undef, nr, N)
    dfy = Vector{Float64}(undef, nr)

    #This is nt is the finest one
    nt = 2 * n
    dfYc= Matrix{Float64}(undef, nt, nr)

    QuadT = makeThetaQuad(nt, nr, q, s, wz1, wz2)
    WorkT = makeThetaWork(nt, nr, N)

    Sθ₁ = QuadT.Sθ₁
    Cθ₁ = QuadT.Cθ₁
    Sθ₂ = QuadT.Sθ₂
    Cθ₂ = QuadT.Cθ₂
    r₁  = QuadT.r₁
    r₂  = QuadT.r₂
    r₃  = QuadT.r₃
    rs₁ = QuadT.rs₁
    rs₂ = QuadT.rs₂
    rs₃ = QuadT.rs₃
    fwmt1 = QuadT.fwmt1
    fwmt2 = QuadT.fwmt2
    fwmt3 = QuadT.fwmt3

    d1 = WorkT.d1
    d2 = WorkT.d2
    y1 = WorkT.y1 
    y1tmp = WorkT.y1tmp 
    y2 = WorkT.y2 
    dfY= WorkT.dfY 
    TY = WorkT.TY 
    J1 = WorkT.J1 
    Jt = WorkT.Jt 
    DIF= WorkT.DIF 
    Zx = WorkT.Zx 
    Zy = WorkT.Zy 
    DJ = WorkT.DJ

    for i in 1:nt
        cθ = sqrt(2) * (cos(QuadT.thet2[i] + π / 4) / QuadT.Cθ₂[i])
        for j in 1:nr
            dfYc[i, j] = (wz3[j] + cθ * wz1[j])^(s - 1.0)
        end
    end

    J2  = Matrix{Float64}(undef, nr, N)
    G   = Matrix{Float64}(undef, nr, N)

    # This is the case when s>=0.5
    # Go over all the points except the boundary
    Lₚ = dp.pthgo[M+1] - 1

    Dhc = Vector{Float64}(undef,M)

    for k = 1:M
        Dhc[k] = d.pths[k].ck1 - d.pths[k].ck0
    end

    # knbd are patches which are not the boundary patches
    # and d.kd are patches which are touching the boundary
    knbd = setdiff(collect(1:M), d.kd)

    IntS = zeros(Float64, Np, Lₚ)

    #--------------------------------------------
    #------------Singular Integration------------
    #--Singular Integration of Boundary patches--
    @inbounds for j in 1:Np
        # j in the linear index of the Chebyshev
        # point. We will compute precomps of all 
        # points τₖ(x₁,x₂) at once by now varying
        # the patches k. This removes some 
        # repeated calculations and thereby saves 
        # time!
        qq, rr = divrem(j - 1, N)
        qq = qq + 1
        rr = rr + 1
        x1, x2 = zp[rr], zp[qq]
        x2p = zp1[qq]
        x1m = zp2[rr]
        x2m = zp2[qq]
        x̃₃ = -x1m * x2p
        x̃₄ = x1m * x2m

        # Singular and k is a boundary patch!
        #-----------------3rd part-----------------
        #------------1ˢᵗ half of part 3------------
        @. d1 = x1m * Cθ₁
        @. d2 = x2p * Sθ₁
        @. yr = x1 - x1m * wz1 / 2
        y1 .= yr'
        @. y2 = x2 - x2p * rs₂ / 2
        
        ChebyTN!(TNy, N, yr)
        qw3func!(dfy, p, zr, s; α=-1.0)

        @inbounds for i₂ in 1:N
            @views Ty = TY[:, :, i₂]
            ChebyT!(Ty, i₂ - 1, y2)
        end

        @inbounds for k in d.kd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₂, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            c₁ = (Dhc[k] * (-x1m / 4.0))^(s - 1.0)

            J2 .= (c₁ .* dfy .* fwr) .* TNy

            @inbounds for i₂ in 1:N
                @views Ty = TY[:, :, i₂]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₂), Jt', fwmt2)
            end

            mul!(I, J2', G)
            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₃ * I[m]
            end

        end

        # Singular and k is a boundary patch!
        #------------2ⁿᵈ half of part 3------------
        @. d1 = x1m * Sθ₂
        @. d2 = x2p * Cθ₂
        @. yr = x2 - x2p * wz1 / 2
        @. y1 = x1 - x1m * rs₃ / 2
        y2 .= yr'

        @. y1tmp = 1 - y1

        ChebyTN!(TNy, N, yr)

        @inbounds for i₁ in 1:N
            @views Ty = TY[:, :, i₁]
            ChebyT!(Ty, i₁ - 1, y1)
        end

        @inbounds for k in d.kd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₃, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            c₁ = (Dhc[k] * (-x1m / 4.0))^(s - 1.0)

            @. J1 = c₁ * abs(DIF)^(-2 * s) * DJ * dfYc

            J2 .= fwmr3 .* TNy

            @inbounds for i₁ in 1:N
                @views Ty = TY[:, :, i₁]
                @. Jt = Ty .* J1
                # G[:, i₁] = Jt' * fwmt1
                mul!(view(G, :, i₁), Jt', fwmt3)
            end

            # I = G' * J2   (N×nr * nr×N → N×N)
            mul!(I, G', J2)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₃ * I[m]
            end
        end

        # Singular and k is a boundary patch!
        #-----------------4th part-----------------
        #------------1ˢᵗ half of part 4------------
        @. d1 = x1m * Cθ₁
        @. d2 = x2m * Sθ₁
        @. yr = x1 - x1m * wz1 / 2
        y1 .= yr'
        @. y2 = x2 - x2m * rs₂ / 2

        ChebyTN!(TNy, N, yr)
        qw3func!(dfy, p, zr, s; α=-1.0)

        @inbounds for i₂ in 1:N
            @views Ty = TY[:, :, i₂]
            ChebyT!(Ty, i₂ - 1, y2)
        end

        @inbounds for k in d.kd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₂, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            c₁ = (Dhc[k] * (-x1m / 4.0))^(s - 1.0)

            J2 .= (c₁ .* dfy .* fwr) .* TNy

            @inbounds for i₂ in 1:N
                @views Ty = TY[:, :, i₂]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₂), Jt', fwmt2)
            end

            mul!(I, J2', G)
            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₄ * I[m]
            end

        end

        # Singular and k is a boundary patch!
        #------------2ⁿᵈ half of part 4------------
        @. d1 = x1m * Sθ₂
        @. d2 = x2m * Cθ₂
        @. yr = x2 - x2m * wz1 / 2
        @. y1 = x1 - x1m * rs₃ / 2
        y2 .= yr'

        @. y1tmp = 1 - y1

        ChebyTN!(TNy, N, yr)

        @inbounds for i₁ in 1:N
            @views Ty = TY[:, :, i₁]
            ChebyT!(Ty, i₁ - 1, y1)
        end

        @inbounds for k in d.kd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₃, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            c₁ = (Dhc[k] * (-x1m / 4.0))^(s - 1.0)

            @. J1 = c₁ * abs(DIF)^(-2 * s) * DJ * dfYc

            J2 .= fwmr3 .* TNy

            @inbounds for i₁ in 1:N
                @views Ty = TY[:, :, i₁]
                @. Jt = Ty .* J1
                # G[:, i₁] = Jt' * fwmt1
                mul!(view(G, :, i₁), Jt', fwmt3)
            end

            # I = G' * J2   (N×nr * nr×N → N×N)
            mul!(I, G', J2)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₄ * I[m]
            end
        end

    end

    #-Singular Integration of Interior patches-
    @inbounds for j in 1:Np
        # j in the linear index of the Chebyshev
        # point. We will compute precomps of all 
        # points τₖ(x₁,x₂) at once by now varying
        # the patches k. This removes some 
        # repeated calculations and thereby saves 
        # time!
        qq, rr = divrem(j - 1, N) .+ 1

        x1, x2 = zp[rr], zp[qq]
        x1p, x2p = zp1[rr], zp1[qq]
        x1m, x2m = zp2[rr], zp2[qq]
        x̃₁ = x1p * x2p
        x̃₂ = -x1p * x2m
        x̃₃ = -x1m * x2p
        x̃₄ = x1m * x2m
       
        #Singular integration is sum of four parts
        #I1,I2,I3,I4, they are all N*N matrices
        #-----------------1st part-----------------
        #------------1ˢᵗ half of part 1------------
        @. d1 = x1p * Cθ₁
        @. d2 = x2p * Sθ₁
        @. yr = x1 - x1p * wz2
        y1 .= yr'
        @. y2 = x2 - x2p * rs₁

        ChebyTN!(TNy, N, yr)
        @. yr = 1 - yr

        @inbounds for i₂ in 1:N
            @views Ty = TY[:, :, i₂]
            ChebyT!(Ty, i₂ - 1, y2)
        end

        @inbounds for k in 1:M

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            dfy .*= fwmr1
            J2 .= dfy .* TNy

            @inbounds for i₂ in 1:N
                @views Ty = TY[:, :, i₂]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₂), Jt', fwmt1)
            end

            mul!(I, J2', G)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₁ * I[m]
            end
        end

        #------------2ⁿᵈ half of part 1------------

        @. d1 = x1p * Sθ₁
        @. d2 = x2p * Cθ₁
        @. yr = x2 - x2p * wz2
        @. y1 = x1 - x1p * rs₁
        y2 .= yr'
        @. y1tmp = 1 - y1

        ChebyTN!(TNy, N, yr)

        @inbounds for i₁ in 1:N
            @views Ty = TY[:, :, i₁]
            ChebyT!(Ty, i₁ - 1, y1)
        end

        @inbounds for k in 1:M
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            J2 .= fwmr1 .* TNy

            @inbounds for i₁ in 1:N
                @views Ty = TY[:, :, i₁]
                @. Jt = Ty .* J1
                # G[:, i₁] = Jt' * fwmt1
                mul!(view(G, :, i₁), Jt', fwmt1)
            end

            # I = G' * J2   (N×nr * nr×N → N×N)
            mul!(I, G', J2)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₁ * I[m]
            end
        end

        #-----------------2nd part-----------------
        #------------1ˢᵗ half of part 2------------
        @. d1 = x1p * Cθ₁
        @. d2 = x2m * Sθ₁
        @. yr = x1 - x1p * wz2
        y1 .= yr'
        @. y2 = x2 - x2m * rs₁

        ChebyTN!(TNy, N, yr)
        @. yr = 1 - yr

        @inbounds for i₂ in 1:N
            @views Ty = TY[:, :, i₂]
            ChebyT!(Ty, i₂ - 1, y2)
        end

        @inbounds for k in 1:M
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            dfy .*= fwmr1
            J2 .= dfy .* TNy

            @inbounds for i₂ in 1:N
                @views Ty = TY[:, :, i₂]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₂), Jt', fwmt1)
            end

            mul!(I, J2', G)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₂ * I[m]
            end
        end

        #------------2ⁿᵈ half of part 2------------

        @. d1 = x1p * Sθ₁
        @. d2 = x2m * Cθ₁
        @. yr = x2 - x2m * wz2
        @. y1 = x1 - x1p * rs₁
        y2 .= yr'
        @. y1tmp = 1 - y1

        ChebyTN!(TNy, N, yr)

        @inbounds for i₁ in 1:N
            @views Ty = TY[:, :, i₁]
            ChebyT!(Ty, i₁ - 1, y1)
        end

        @inbounds for k in 1:M
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            J2 .= fwmr1 .* TNy

            @inbounds for i₁ in 1:N
                @views Ty = TY[:, :, i₁]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₁), Jt', fwmt1)
            end

            mul!(I, G', J2)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₂ * I[m]
            end
        end

        #-----------------3rd part-----------------
        #------------1ˢᵗ half of part 3------------
        @. d1 = x1m * Cθ₁
        @. d2 = x2p * Sθ₁
        @. yr = x1 - x1m * wz2
        y1 .= yr'
        @. y2 = x2 - x2p * rs₁

        ChebyTN!(TNy, N, yr)
        @. yr = 1 - yr

        @inbounds for i₂ in 1:N
            @views Ty = TY[:, :, i₂]
            ChebyT!(Ty, i₂ - 1, y2)
        end

        @inbounds for k in knbd

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            dfy .*= fwmr1
            J2 .= dfy .* TNy

            @inbounds for i₂ in 1:N
                @views Ty = TY[:, :, i₂]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₂), Jt', fwmt1)
            end

            mul!(I, J2', G)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₃ * I[m]
            end
        end

        #------------2ⁿᵈ half of part 3------------

        @. d1 = x1m * Sθ₁
        @. d2 = x2p * Cθ₁
        @. yr = x2 - x2p * wz2
        @. y1 = x1 - x1m * rs₁
        y2 .= yr'

        @. y1tmp = 1 - y1

        ChebyTN!(TNy, N, yr)

        @inbounds for i₁ in 1:N
            @views Ty = TY[:, :, i₁]
            ChebyT!(Ty, i₁ - 1, y1)
        end

        @inbounds for k in knbd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            J2 .= fwmr1 .* TNy

            @inbounds for i₁ in 1:N
                @views Ty = TY[:, :, i₁]
                @. Jt = Ty .* J1
                # G[:, i₁] = Jt' * fwmt1
                mul!(view(G, :, i₁), Jt', fwmt1)
            end

            # I = G' * J2   (N×nr * nr×N → N×N)
            mul!(I, G', J2)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₃ * I[m]
            end
        end

        #-----------------4th part-----------------
        #------------1ˢᵗ half of part 4------------
        @. d1 = x1m * Cθ₁
        @. d2 = x2m * Sθ₁
        @. yr = x1 - x1m * wz2
        y1 .= yr'
        @. y2 = x2 - x2m * rs₁

        ChebyTN!(TNy, N, yr)
        @. yr = 1 - yr

        @inbounds for i₂ in 1:N
            @views Ty = TY[:, :, i₂]
            ChebyT!(Ty, i₂ - 1, y2)
        end

        @inbounds for k in knbd

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            dfy .*= fwmr1
            J2 .= dfy .* TNy

            @inbounds for i₂ in 1:N
                @views Ty = TY[:, :, i₂]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₂), Jt', fwmt1)
            end

            mul!(I, J2', G)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₄ * I[m]
            end
        end
        #------------2ⁿᵈ half of part 4------------
        @. d1 = x1m * Sθ₁
        @. d2 = x2m * Cθ₁
        @. yr = x2 - x2m * wz2
        y2 .= yr'
        @. y1 = x1 - x1m * rs₁
        @. y1tmp = 1 - y1

        ChebyTN!(TNy, N, yr)

        @inbounds for i₁ in 1:N
            @views Ty = TY[:, :, i₁]
            ChebyT!(Ty, i₁ - 1, y1)
        end

        @inbounds for k in knbd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            J2 .= fwmr1 .* TNy

            @inbounds for i₁ in 1:N
                @views Ty = TY[:, :, i₁]
                @. Jt = Ty .* J1
                # G[:, i₁] = Jt' * fwmt1
                mul!(view(G, :, i₁), Jt', fwmt1)
            end

            # I = G' * J2   (N×nr * nr×N → N×N)
            mul!(I, G', J2)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₄ * I[m]
            end
        end

    end

    #Clearing Variables (mainly large matrices)
    QuadT = nothing
    WorkT = nothing
    #--------------------------------------------
    #----------Near Singular Integration---------
    # z1 = (1+z)/2, z2 = (1-z)/2 where z = cos(pi*(2j-1)/(2n)) 
    # are Chebyshev nodes in the interval [-1,1]. This is for 
    # near singular integrals 
    #n = n #Maybe n here can be fixed to 128 ? (only for small s)
    z = Vector{Float64}(undef, n)
    z1 = Vector{Float64}(undef, n)
    z2 = Vector{Float64}(undef, n)

    #Store Fejer 1st quadrature weights for near singular integrals
    fw = getF1W(n)  

    d1 = Vector{Float64}(undef, n)
    d2 = similar(d1)

    @inbounds for i in 1:n
        c₁ = π * (2 * i - 1) / (2 * n)
        z[i] = cos(c₁)
        z1[i] = cos(c₁ / 2)^2
        z2[i] = sin(c₁ / 2)^2
    end

    y1  = Vector{Float64}(undef, n)
    wdf = similar(y1)
    y1tmp=similar(y1)
    y2  = similar(y1)
    t1  = Matrix{Float64}(undef, n, n)  # meshgrid of y1/y2 
    t2  = Matrix{Float64}(undef, n, n)  # (column = y2[j], row = y1[i])
    DJ  = Matrix{Float64}(undef, n, n)  # To store the Jacobian
    DIF = similar(DJ)
    Zx  = similar(DJ)
    Zy  = similar(DJ)

    TN_y1 = Matrix{Float64}(undef, n, N)
    TN_y2 = Matrix{Float64}(undef, n, N)

    TNL = Matrix{Float64}(undef, N, n)
    TNR = Matrix{Float64}(undef, n, N)
    A   = Matrix{Float64}(undef, n, n)    
    Tmp = Matrix{Float64}(undef, N, n)
    dw  = Vector{Float64}(undef, n)
    dfy = Vector{Float64}(undef, n)

    I₁ = Matrix{Float64}(undef, N, N)
    I₂ = Matrix{Float64}(undef, N, N)
    qw = Vector{Float64}(undef, n)

    #-------------------------------------------
    # small helpers (no allocations)
    @inline function fill_meshgrid!(T1, T2, y1, y2)
        @inbounds for j in eachindex(y2)
            @inbounds for i in eachindex(y1)
                T1[i, j] = y1[i]
                T2[i, j] = y2[j]
            end
        end
        return nothing
    end
    #-------------------------------------------
    #A vector of Bool, initialized to true for all
    #indices from 1:Lₚ. They will be updated as 
    #False  for singular points. (and ofcourse the
    #points left are near singular, which are true)
    NSI = trues(Lₚ)

    @inbounds for j in 1:Np
        @inbounds for k in 1:M
            NSI[dp.pthgo[k]+j-1] = false
        end
    end

    #Now the big for loop for near singular points
    @inbounds for i in 1:Lₚ
        if NSI[i] == true

            k = dp.prepts[2,i]

            ll = i - k * Np

            α₁ = dp.invpts[1, ll]
            α₂ = dp.invpts[2, ll]

            if 1 ≤ α₁
                B1 = -winv(p, (α₁ - 1) / (α₁ + 1))
            elseif α₁ ≤ -1
                A1 = -winv(p, (-1 - α₁) / (1 - α₁))
            end

            if 1 ≤ α₂
                B2 = -winv(p, (α₂ - 1) / (α₂ + 1))
            elseif α₂ ≤ -1
                A2 = -winv(p, (-1 - α₂) / (1 - α₂))
            end

            #Not a boundary patch
            if k in knbd
                if 1 ≤ α₁
                    wfunc!(d1, p, z1; α=-B1, β=α₁ + 1.0)
                    @. y1 = α₁ - d1
                    @. y1tmp = 1 - y1
                    if 1 ≤ α₂
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* dw(B1*z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1; α=B1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= B1 * B2 * (α₁ + 1) * (α₂ + 1) / 4

                    elseif -1 < α₂ && α₂ < 1
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(B1*z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1; α=B1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(z1)'.*fw'
                        dwfunc!(dw, p, z1)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        #Left weights same as before
                        #Right weights: TN(y2).*dw(z2)'.*fw'
                        ChebyTN!(TN_y2, N, y2)
                        dwfunc!(dw, p, z2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = B1*(α₁ + 1)*( (α₂ + 1)*I₁ + (1 - α₂)*I₂ )/4

                    elseif α₂ ≤ -1
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(B1*z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1; α=B1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= B1 * A2 * (α₁ + 1) * (1 - α₂) / 4
                    end
                elseif -1 < α₁ && α₁ < 1
                    if 1 ≤ α₂
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d1, p, z2; α=-1.0, β=α₁ - 1.0)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z2) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = B2*(α₂ + 1)*( (α₁ + 1)*I₁ + (1 - α₁)*I₂ )/4

                    elseif α₂ ≤ -1
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d1, p, z2; α=-1.0, β=α₁ - 1.0)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z2) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = A2*(1 - α₂)*( (α₁ + 1)*I₁ + (1 - α₁)*I₂ )/4
                    else
                        error("Fatal error! mapinv not working correctly.")
                    end

                elseif α₁ ≤ -1
                    wfunc!(d1, p, z2; α=-A1, β=α₁ - 1.0)
                    @. y1 = α₁ - d1
                    @. y1tmp = 1 - y1
                    if 1 ≤ α₂
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* dw(A1*z2) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2; α=A1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= A1 * B2 * (1 - α₁) * (α₂ + 1) / 4

                    elseif -1 < α₂ && α₂ < 1
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(A1*z2) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2; α=A1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(z1)'.*fw'
                        dwfunc!(dw, p, z1)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        #Left weights same as before
                        #Right weights: TN(y2).*dw(z2)'.*fw'
                        ChebyTN!(TN_y2, N, y2)
                        dwfunc!(dw, p, z2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = A1*(1 - α₁)*( (α₂ + 1)*I₁ + (1 - α₂)*I₂ )/4
                    elseif α₂ ≤ -1
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2;

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(A1*z2) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2; α=A1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= A1 * A2 * (1 - α₁) * (1 - α₂) / 4
                    end
                end
            else
                hc = Dhc[k]
                #Near singular for boundary patches.
                if 1 ≤ α₁ 
                    #α₁ = 1
                    wfunc!(d1, p, z1; α=-1.0, β=2.0)
                    @. y1 = 1 - d1
                    if 1 ≤ α₂
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* qw .* TN(y1)'
                        qw2func!(qw, p, z1, s)
                        @. wdf = fw * qw
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= hc^(s-1) * B2 * (α₂ + 1) / 2

                    elseif -1 < α₂ && α₂ < 1
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* qw .* TN(y1)'
                        qw2func!(qw, p, z1, s)
                        @. wdf = fw * qw
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(z1)'.*fw'
                        dwfunc!(dw, p, z1)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        #Left weights same as before
                        #Right weights: TN(y2).*dw(z2)'.*fw'
                        ChebyTN!(TN_y2, N, y2)
                        dwfunc!(dw, p, z2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = hc^(s-1) * ( (α₂ + 1)*I₁ + (1 - α₂)*I₂ )/2

                    elseif α₂ ≤ -1
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(B1*z1) .* dfy .* TN(y1)'
                        qw2func!(qw, p, z1, s)
                        @. wdf = fw * qw
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= hc^(s-1) * A2 * (1 - α₂) / 2
                    end
                elseif -1 < α₁ && α₁ < 1
                    if 1 ≤ α₂
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d1, p, z; α=1.0, β=(α₁ - 1.0) / 2)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* qw .* TN(y1)'
                        qw2func!(qw, p, z, s)
                        @. wdf = fw * qw
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = B2*(α₂ + 1)*( (α₁ + 1)*I₁ + (1 - α₁) * (hc*(1-α₁)/4)^(s-1) * I₂ )/4

                    elseif α₂ ≤ -1
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d1, p, z; α=1.0, β=(α₁ - 1.0) / 2)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* qw .* TN(y1)'
                        qw2func!(qw, p, z, s)
                        @. wdf = fw * qw
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = A2*(1 - α₂)*( (α₁ + 1)*I₁ + (1 - α₁)* (hc*(1-α₁)/4)^(s-1) * I₂ )/4
                    else
                        error("Fatal error! mapinv not working correctly.")
                    end

                elseif α₁ ≤ -1
                    A1 = 1 - winv(p, 2 * (-1 - α₁) / (1 - α₁))
                    wfunc!(d1, p, z2; α=-A1, β=(α₁ - 1.0) / 2, γ=1.0)
                    @. y1 = α₁ - d1
                    if 1 ≤ α₂
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* qw .* TN(y1)'
                        qw2func!(qw, p, z2, s; α=-A1, β=1.0, γ=1.0)
                        @. wdf = fw * qw
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= (hc*(1-α₁)/4)^(s-1) * A1 * B2 * (1-α₁)*(α₂+1)/8

                    elseif -1 < α₂ && α₂ < 1
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* qw .* TN(y1)'
                        qw2func!(qw, p, z2, s; α=-A1, β=1.0, γ=1.0)
                        @. wdf = fw * qw
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(z1)'.*fw'
                        dwfunc!(dw, p, z1)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        #Left weights same as before
                        #Right weights: TN(y2).*dw(z2)'.*fw'
                        ChebyTN!(TN_y2, N, y2)
                        dwfunc!(dw, p, z2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = (hc*(1-α₁)/4)^(s-1) *A1*(1-α₁)*( (α₂+1)*I₁ + (1-α₂)*I₂ )/8
                        
                    elseif α₂ ≤ -1
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(B1*z1) .* dfy .* TN(y1)'
                        qw2func!(qw, p, z2, s; α=-A1, β=1.0, γ=1.0)
                        @. wdf = fw * qw
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= (hc*(1-α₁)/4)^(s-1) * A1*A2*(1-α₁)*(1-α₂)/8

                    end
                end

            end

            #Allocation free version of IntS[:,i] = reshape(I,(Np,1));
            @inbounds for m in 1:Np
                IntS[m, i] = I[m]
            end

        end
    end
  
    return IntS
end

#For s<0.5
function precompsL(d::abstractdomain, dp::domprop, s::Float64, p::Int
    ; n::Int=128)::Matrix{Float64}

    #Bookkeeping
    M = d.Npat  # number of patches
    N = dp.N    # number of nodes per patch per Axis
    Np = N * N  # number of nodes per patch

    Mbd= length(d.kd) #Number of bonundary patches
    nbd = Mbd * N
    # Col index where the bd near–singular block starts in `prepts`
    bdnsp = dp.pthgo[M+1] + nbd

    # zp1 = zp+1, zp2 = xp-1 where zp = cos(pi*(2j-1)/(2N))
    zp  = Vector{Float64}(undef, N)
    zp1 = Vector{Float64}(undef, N)
    zp2 = Vector{Float64}(undef, N)

    @inbounds for j in 1:N
        zp[j]  = cospi((2 * j - 1) / (2N))
        zp1[j] = 2 * cospi((2 * j - 1) / (4N))^2
        zp2[j] =-2 * sinpi((2 * j - 1) / (4N))^2
    end

    #Parameter q for another polynomial change of variables is theta. 
    # It is usually taken as 4
    q = 4

    #Define nodes in radial and weights
    nr = n     
    fwr = getF1W(nr)

    zr = Vector{Float64}(undef, nr)
    z2r = Vector{Float64}(undef, nr)

    @inbounds for i in 1:nr
        c₁ = π * (2 * i - 1) / (2 * nr)
        zr[i] = cos(c₁)
        z2r[i] = sin(c₁ / 2)^2
    end

    #Precomputing modified weights
    fwmr1 = Vector{Float64}(undef, nr)
    fwmr2 = Vector{Float64}(undef, nr)
    fwmr3 = Vector{Float64}(undef, nr)

    #To store w(zr) and w(-z2r)
    wz1 = Vector{Float64}(undef, nr)
    wz2 = similar(wz1)
    wz3 = similar(wz1)

    wfunc!(wz1, p, zr)
    wfunc!(wz2, p, z2r; α=-1.0)
    wfunc!(wz3, p, zr; α=-1.0)

    qw1func!(fwmr1, p, z2r, s)
    qw1func!(fwmr3, p, zr, s; α=-1.0)

    @. fwmr1 = fwr * fwmr1
    @. fwmr2 = fwr * fwmr3
    @. fwmr3 = fwr * fwmr3

    # Storage of temporary variables inside the for loop
    # These are independent on nt
    I = Matrix{Float64}(undef, N, N)

    yr  = Vector{Float64}(undef, nr)

    TNy = Matrix{Float64}(undef, nr, N)
    dfy = Vector{Float64}(undef, nr)

    #This is nt is the finest one
    nt = 2*n
    dfYc= Matrix{Float64}(undef, nt, nr)

    QuadT = makeThetaQuad(nt, nr, q, s, wz1, wz2)
    WorkT = makeThetaWork(nt, nr, N)

    Sθ₁ = QuadT.Sθ₁
    Cθ₁ = QuadT.Cθ₁
    Sθ₂ = QuadT.Sθ₂
    Cθ₂ = QuadT.Cθ₂
    r₁  = QuadT.r₁
    r₂  = QuadT.r₂
    r₃  = QuadT.r₃
    rs₁ = QuadT.rs₁
    rs₂ = QuadT.rs₂
    rs₃ = QuadT.rs₃
    fwmt1 = QuadT.fwmt1
    fwmt2 = QuadT.fwmt2
    fwmt3 = QuadT.fwmt3

    d1 = WorkT.d1
    d2 = WorkT.d2
    y1 = WorkT.y1 
    y1tmp = WorkT.y1tmp 
    y2 = WorkT.y2 
    dfY= WorkT.dfY 
    TY = WorkT.TY 
    J1 = WorkT.J1 
    Jt = WorkT.Jt 
    DIF= WorkT.DIF 
    Zx = WorkT.Zx 
    Zy = WorkT.Zy 
    DJ = WorkT.DJ

    @inbounds for i in 1:nt
        cθ = sqrt(2) * (cos(QuadT.thet2[i] + π / 4) / QuadT.Cθ₂[i])
        @inbounds for j in 1:nr
            dfYc[i, j] = (wz3[j] + cθ * wz1[j])^s
        end
    end

    J2  = Matrix{Float64}(undef, nr, N)
    G   = Matrix{Float64}(undef, nr, N)

    # This is the case when s<0.5
    Lₚ = size(dp.prepts,2) # > dp.pthgo[M+1] - 1

    # First we will go over all the points in the interioir

    Dhc = Vector{Float64}(undef,M)

    @inbounds for k = 1:M
        Dhc[k] = d.pths[k].ck1 - d.pths[k].ck0
    end

    # knbd are patches which are not the boundary patches
    # and d.kd are patches are touching the boundary
    knbd = setdiff(collect(1:M), d.kd)

    IntS = zeros(Float64, Np, Lₚ)

    #--------------------------------------------
    #------------Singular Integration------------
    #--Singular Integration of Boundary patches--
    @inbounds for j in 1:Np
        # j in the linear index of the Chebyshev
        # point. We will compute precomps of all 
        # points τₖ(x₁,x₂) at once by now varying
        # the patches k. This removes some 
        # repeated calculations and thereby saves 
        # time!
        qq, rr = divrem(j - 1, N)
        qq = qq + 1
        rr = rr + 1
        x1, x2 = zp[rr], zp[qq]
        x2p = zp1[qq]
        x1m = zp2[rr]
        x2m = zp2[qq]
        x̃₃ = -x1m * x2p
        x̃₄ = x1m * x2m

        # Singular and k is a boundary patch!
        #-----------------3rd part-----------------
        #------------1ˢᵗ half of part 3------------
        @. d1 = x1m * Cθ₁
        @. d2 = x2p * Sθ₁
        @. yr = x1 - x1m * wz1 / 2
        y1 .= yr'
        @. y2 = x2 - x2p * rs₂ / 2
        
        ChebyTN!(TNy, N, yr)

        @. dfy = wz3 ^ s

        @inbounds for i₂ in 1:N
            @views Ty = TY[:, :, i₂]
            ChebyT!(Ty, i₂ - 1, y2)
        end

        @inbounds for k in d.kd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₂, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            c₁ = (Dhc[k] * (-x1m / 4.0))^s

            J2 .= (c₁ .* dfy .* fwmr2) .* TNy

            @inbounds for i₂ in 1:N
                @views Ty = TY[:, :, i₂]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₂), Jt', fwmt2)
            end

            mul!(I, J2', G)
            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₃ * I[m]
            end

        end

        # Singular and k is a boundary patch!
        #------------2ⁿᵈ half of part 3------------
        @. d1 = x1m * Sθ₂
        @. d2 = x2p * Cθ₂
        @. yr = x2 - x2p * wz1 / 2
        @. y1 = x1 - x1m * rs₃ / 2
        y2 .= yr'

        ChebyTN!(TNy, N, yr)

        @inbounds for i₁ in 1:N
            @views Ty = TY[:, :, i₁]
            ChebyT!(Ty, i₁ - 1, y1)
        end

        @inbounds for k in d.kd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₃, d1, d2, k)

            c₁ = (Dhc[k] * (-x1m / 4.0))^s

            @. J1 = c₁ * abs(DIF)^(-2 * s) * DJ * dfYc

            J2 .= fwmr3 .* TNy

            @inbounds for i₁ in 1:N
                @views Ty = TY[:, :, i₁]
                @. Jt = Ty .* J1
                # G[:, i₁] = Jt' * fwmt1
                mul!(view(G, :, i₁), Jt', fwmt3)
            end

            # I = G' * J2   (N×nr * nr×N → N×N)
            mul!(I, G', J2)

            @inbounds  for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₃ * I[m]
            end
        end

        # Singular and k is a boundary patch!
        #-----------------4th part-----------------
        #------------1ˢᵗ half of part 4------------
        @. d1 = x1m * Cθ₁
        @. d2 = x2m * Sθ₁
        @. yr = x1 - x1m * wz1 / 2
        y1 .= yr'
        @. y2 = x2 - x2m * rs₂ / 2

        ChebyTN!(TNy, N, yr)

        @inbounds for i₂ in 1:N
            @views Ty = TY[:, :, i₂]
            ChebyT!(Ty, i₂ - 1, y2)
        end

        @inbounds for k in d.kd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₂, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            c₁ = (Dhc[k] * (-x1m / 4.0))^s

            J2 .= (c₁ .* dfy .* fwmr2) .* TNy

            @inbounds for i₂ in 1:N
                @views Ty = TY[:, :, i₂]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₂), Jt', fwmt2)
            end

            mul!(I, J2', G)
            @inbounds  for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₄ * I[m]
            end

        end

        # Singular and k is a boundary patch!
        #------------2ⁿᵈ half of part 4------------
        @. d1 = x1m * Sθ₂
        @. d2 = x2m * Cθ₂
        @. yr = x2 - x2m * wz1 / 2
        @. y1 = x1 - x1m * rs₃ / 2
        y2 .= yr'

        ChebyTN!(TNy, N, yr)

        @inbounds for i₁ in 1:N
            @views Ty = TY[:, :, i₁]
            ChebyT!(Ty, i₁ - 1, y1)
        end

        @inbounds for k in d.kd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₃, d1, d2, k)

            c₁ = (Dhc[k] * (-x1m / 4.0))^s

            @. J1 = c₁ * abs(DIF)^(-2 * s) * DJ * dfYc

            J2 .= fwmr3 .* TNy

            @inbounds for i₁ in 1:N
                @views Ty = TY[:, :, i₁]
                @. Jt = Ty .* J1
                # G[:, i₁] = Jt' * fwmt1
                mul!(view(G, :, i₁), Jt', fwmt3)
            end

            # I = G' * J2   (N×nr * nr×N → N×N)
            mul!(I, G', J2)

            @inbounds  for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₄ * I[m]
            end
        end

    end

    Lₚₛ = dp.pthgo[M+1] - 1
    #-Singular Integration For pts on bonudary-
    #These points are ofcourse always on bd patches
    @inbounds for j in 1:N
        # j in the linear index of the Chebyshev
        # point. We will compute precomps of all 
        # points τₖ(x₁,x₂) at once by now varying
        # the patches k.

        x2 = zp[j]
        x2p = zp1[j]
        x2m = zp2[j]
        x̃₁ = 2.0 * x2p
        x̃₂ = -2.0 * x2m

        #Singular integration is sum of four parts
        #I1,I2,I3,I4, they are all N*N matrices
        #-----------------1st part-----------------
        #------------1ˢᵗ half of part 1------------
        @. d1 = 2.0 * Cθ₁
        @. d2 = x2p * Sθ₁
        @. yr = 1.0 - 2.0 * wz2
        y1 .= yr'
        @. y2 = x2 - x2p * rs₁

        ChebyTN!(TNy, N, yr)
        @. yr = 2.0 * wz2

        @inbounds for i₂ in 1:N
            @views Ty = TY[:, :, i₂]
            ChebyT!(Ty, i₂ - 1, y2)
        end

        @inbounds for ℓ in 1:Mbd

            k = d.kd[ℓ]

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                1.0, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            dfy .*= fwmr1
            J2 .= dfy .* TNy

            @inbounds for i₂ in 1:N
                @views Ty = TY[:, :, i₂]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₂), Jt', fwmt1)
            end

            mul!(I, J2', G)

            @inbounds  for m in 1:Np
                IntS[m, Lₚₛ+(ℓ-1)*N+j] += x̃₁ * I[m]
            end
        end

        #------------2ⁿᵈ half of part 1------------

        @. d1 = 2.0 * Sθ₁
        @. d2 = x2p * Cθ₁
        @. yr = x2 - x2p * wz2
        @. y1 = 1.0 - 2.0 * rs₁
        y2 .= yr'
        @. y1tmp = 2.0 * rs₁

        ChebyTN!(TNy, N, yr)

        @inbounds for i₁ in 1:N
            @views Ty = TY[:, :, i₁]
            ChebyT!(Ty, i₁ - 1, y1)
        end

        @inbounds for ℓ in 1:Mbd

            k = d.kd[ℓ]

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                1.0, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            J2 .= fwmr1 .* TNy

            @inbounds for i₁ in 1:N
                @views Ty = TY[:, :, i₁]
                @. Jt = Ty .* J1
                # G[:, i₁] = Jt' * fwmt1
                mul!(view(G, :, i₁), Jt', fwmt1)
            end

            # I = G' * J2   (N×nr * nr×N → N×N)
            mul!(I, G', J2)

            @inbounds  for m in 1:Np
                IntS[m, Lₚₛ+(ℓ-1)*N+j] += x̃₁ * I[m]
            end
        end

        #-----------------2nd part-----------------
        #------------1ˢᵗ half of part 2------------
        @. d1 = 2.0 * Cθ₁
        @. d2 = x2m * Sθ₁
        @. yr = 1.0 - 2.0 * wz2
        y1 .= yr'
        @. y2 = x2 - x2m * rs₁

        ChebyTN!(TNy, N, yr)
        @. yr = 2.0 * wz2

        @inbounds for i₂ in 1:N
            @views Ty = TY[:, :, i₂]
            ChebyT!(Ty, i₂ - 1, y2)
        end

        @inbounds for ℓ in 1:Mbd

            k = d.kd[ℓ]

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                1.0, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            dfy .*= fwmr1
            J2 .= dfy .* TNy

            @inbounds for i₂ in 1:N
                @views Ty = TY[:, :, i₂]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₂), Jt', fwmt1)
            end

            mul!(I, J2', G)

            @inbounds for m in 1:Np
                IntS[m, Lₚₛ+(ℓ-1)*N+j] += x̃₂ * I[m]
            end
        end

        #------------2ⁿᵈ half of part 2------------

        @. d1 = 2.0 * Sθ₁
        @. d2 = x2m * Cθ₁
        @. yr = x2 - x2m * wz2
        @. y1 = 1.0 - 2.0 * rs₁
        y2 .= yr'
        @. y1tmp = 2.0 * rs₁

        ChebyTN!(TNy, N, yr)

        @inbounds for i₁ in 1:N
            @views Ty = TY[:, :, i₁]
            ChebyT!(Ty, i₁ - 1, y1)
        end

        @inbounds for ℓ in 1:Mbd

            k = d.kd[ℓ]

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                1.0, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            J2 .= fwmr1 .* TNy

            @inbounds for i₁ in 1:N
                @views Ty = TY[:, :, i₁]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₁), Jt', fwmt1)
            end

            mul!(I, G', J2)

            @inbounds for m in 1:Np
                IntS[m, Lₚₛ+(ℓ-1)*N+j] += x̃₂ * I[m]
            end
        end

    end

    #-Singular Integration of Interior patches-
    @inbounds for j in 1:Np
        # j in the linear index of the Chebyshev
        # point. We will compute precomps of all 
        # points τₖ(x₁,x₂) at once by now varying
        # the patches k. This removes some 
        # repeated calculations and thereby saves 
        # time!
        qq, rr = divrem(j - 1, N) .+ 1

        x1, x2 = zp[rr], zp[qq]
        x1p, x2p = zp1[rr], zp1[qq]
        x1m, x2m = zp2[rr], zp2[qq]
        x̃₁ = x1p * x2p
        x̃₂ = -x1p * x2m
        x̃₃ = -x1m * x2p
        x̃₄ = x1m * x2m
       
        #Singular integration is sum of four parts
        #I1,I2,I3,I4, they are all N*N matrices
        #-----------------1st part-----------------
        #------------1ˢᵗ half of part 1------------
        @. d1 = x1p * Cθ₁
        @. d2 = x2p * Sθ₁
        @. yr = x1 - x1p * wz2
        y1 .= yr'
        @. y2 = x2 - x2p * rs₁

        ChebyTN!(TNy, N, yr)
        @. yr = 1 - yr

        @inbounds for i₂ in 1:N
            @views Ty = TY[:, :, i₂]
            ChebyT!(Ty, i₂ - 1, y2)
        end

        @inbounds for k in 1:M

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            dfy .*= fwmr1
            J2 .= dfy .* TNy

            @inbounds for i₂ in 1:N
                @views Ty = TY[:, :, i₂]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₂), Jt', fwmt1)
            end

            mul!(I, J2', G)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₁ * I[m]
            end
        end

        #------------2ⁿᵈ half of part 1------------

        @. d1 = x1p * Sθ₁
        @. d2 = x2p * Cθ₁
        @. yr = x2 - x2p * wz2
        @. y1 = x1 - x1p * rs₁
        y2 .= yr'
        @. y1tmp = 1 - y1

        ChebyTN!(TNy, N, yr)

        @inbounds for i₁ in 1:N
            @views Ty = TY[:, :, i₁]
            ChebyT!(Ty, i₁ - 1, y1)
        end

        @inbounds for k in 1:M
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            J2 .= fwmr1 .* TNy

            @inbounds for i₁ in 1:N
                @views Ty = TY[:, :, i₁]
                @. Jt = Ty .* J1
                # G[:, i₁] = Jt' * fwmt1
                mul!(view(G, :, i₁), Jt', fwmt1)
            end

            # I = G' * J2   (N×nr * nr×N → N×N)
            mul!(I, G', J2)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₁ * I[m]
            end
        end

        #-----------------2nd part-----------------
        #------------1ˢᵗ half of part 2------------
        @. d1 = x1p * Cθ₁
        @. d2 = x2m * Sθ₁
        @. yr = x1 - x1p * wz2
        y1 .= yr'
        @. y2 = x2 - x2m * rs₁

        ChebyTN!(TNy, N, yr)
        @. yr = 1 - yr

        @inbounds for i₂ in 1:N
            @views Ty = TY[:, :, i₂]
            ChebyT!(Ty, i₂ - 1, y2)
        end

        @inbounds for k in 1:M
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            dfy .*= fwmr1
            J2 .= dfy .* TNy

            @inbounds for i₂ in 1:N
                @views Ty = TY[:, :, i₂]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₂), Jt', fwmt1)
            end

            mul!(I, J2', G)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₂ * I[m]
            end
        end

        #------------2ⁿᵈ half of part 2------------

        @. d1 = x1p * Sθ₁
        @. d2 = x2m * Cθ₁
        @. yr = x2 - x2m * wz2
        @. y1 = x1 - x1p * rs₁
        y2 .= yr'
        @. y1tmp = 1 - y1

        ChebyTN!(TNy, N, yr)

        @inbounds for i₁ in 1:N
            @views Ty = TY[:, :, i₁]
            ChebyT!(Ty, i₁ - 1, y1)
        end

        @inbounds for k in 1:M
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            J2 .= fwmr1 .* TNy

            @inbounds for i₁ in 1:N
                @views Ty = TY[:, :, i₁]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₁), Jt', fwmt1)
            end

            mul!(I, G', J2)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₂ * I[m]
            end
        end

        #-----------------3rd part-----------------
        #------------1ˢᵗ half of part 3------------
        @. d1 = x1m * Cθ₁
        @. d2 = x2p * Sθ₁
        @. yr = x1 - x1m * wz2
        y1 .= yr'
        @. y2 = x2 - x2p * rs₁

        ChebyTN!(TNy, N, yr)
        @. yr = 1 - yr

        @inbounds for i₂ in 1:N
            @views Ty = TY[:, :, i₂]
            ChebyT!(Ty, i₂ - 1, y2)
        end

        @inbounds for k in knbd

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            dfy .*= fwmr1
            J2 .= dfy .* TNy

            @inbounds for i₂ in 1:N
                @views Ty = TY[:, :, i₂]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₂), Jt', fwmt1)
            end

            mul!(I, J2', G)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₃ * I[m]
            end
        end

        #------------2ⁿᵈ half of part 3------------

        @. d1 = x1m * Sθ₁
        @. d2 = x2p * Cθ₁
        @. yr = x2 - x2p * wz2
        @. y1 = x1 - x1m * rs₁
        y2 .= yr'

        @. y1tmp = 1 - y1

        ChebyTN!(TNy, N, yr)

        @inbounds for i₁ in 1:N
            @views Ty = TY[:, :, i₁]
            ChebyT!(Ty, i₁ - 1, y1)
        end

        @inbounds for k in knbd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            J2 .= fwmr1 .* TNy

            @inbounds for i₁ in 1:N
                @views Ty = TY[:, :, i₁]
                @. Jt = Ty .* J1
                # G[:, i₁] = Jt' * fwmt1
                mul!(view(G, :, i₁), Jt', fwmt1)
            end

            # I = G' * J2   (N×nr * nr×N → N×N)
            mul!(I, G', J2)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₃ * I[m]
            end
        end

        #-----------------4th part-----------------
        #------------1ˢᵗ half of part 4------------
        @. d1 = x1m * Cθ₁
        @. d2 = x2m * Sθ₁
        @. yr = x1 - x1m * wz2
        y1 .= yr'
        @. y2 = x2 - x2m * rs₁

        ChebyTN!(TNy, N, yr)
        @. yr = 1 - yr

        @inbounds for i₂ in 1:N
            @views Ty = TY[:, :, i₂]
            ChebyT!(Ty, i₂ - 1, y2)
        end

        @inbounds for k in knbd

            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            @. J1 = abs(DIF)^(-2 * s) * DJ

            dfunc!(dfy, d, k, yr, s)

            dfy .*= fwmr1
            J2 .= dfy .* TNy

            @inbounds for i₂ in 1:N
                @views Ty = TY[:, :, i₂]
                @. Jt = Ty .* J1
                mul!(view(G, :, i₂), Jt', fwmt1)
            end

            mul!(I, J2', G)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₄ * I[m]
            end
        end
        #------------2ⁿᵈ half of part 4------------
        @. d1 = x1m * Sθ₁
        @. d2 = x2m * Cθ₁
        @. yr = x2 - x2m * wz2
        y2 .= yr'
        @. y1 = x1 - x1m * rs₁
        @. y1tmp = 1 - y1

        ChebyTN!(TNy, N, yr)

        @inbounds for i₁ in 1:N
            @views Ty = TY[:, :, i₁]
            ChebyT!(Ty, i₁ - 1, y1)
        end

        @inbounds for k in knbd
            diff_rmap!(DIF, Zx, Zy, DJ, d,
                x1, x2, y1, y2, r₁, d1, d2, k)

            dfunc!(dfY, d, k, y1tmp, s)

            @. J1 = abs(DIF)^(-2 * s) * DJ * dfY

            J2 .= fwmr1 .* TNy

            @inbounds for i₁ in 1:N
                @views Ty = TY[:, :, i₁]
                @. Jt = Ty .* J1
                # G[:, i₁] = Jt' * fwmt1
                mul!(view(G, :, i₁), Jt', fwmt1)
            end

            # I = G' * J2   (N×nr * nr×N → N×N)
            mul!(I, G', J2)

            @inbounds for m in 1:Np
                IntS[m, dp.pthgo[k]+j-1] += x̃₄ * I[m]
            end
        end

    end

    #Clearing Variables (mainly large structs)
    QuadT = nothing
    WorkT = nothing
    #--------------------------------------------
    #----------Near Singular Integration---------
    # z1 = (1+z)/2, z2 = (1-z)/2 where z = cos(pi*(2j-1)/(2n)) 
    # are Chebyshev nodes in the interval [-1,1]. This is for 
    # near singular integrals 
    #n = n #Maybe n here can be fixed to 128 ? (only for small s)
    z = Vector{Float64}(undef, n)
    z1 = Vector{Float64}(undef, n)
    z2 = Vector{Float64}(undef, n)

    #Store Fejer 1st quadrature weights for near singular integrals
    fw = getF1W(n)  

    d1 = Vector{Float64}(undef, n)
    d2 = similar(d1)

    @inbounds for i in 1:n
        c₁ = π * (2 * i - 1) / (2 * n)
        z[i] = cos(c₁)
        z1[i] = cos(c₁ / 2)^2
        z2[i] = sin(c₁ / 2)^2
    end

    y1  = Vector{Float64}(undef, n)
    wdf = similar(y1)
    y1tmp=similar(y1)
    y2  = similar(y1)
    t1  = Matrix{Float64}(undef, n, n)  # meshgrid of y1/y2 
    t2  = Matrix{Float64}(undef, n, n)  # (column = y2[j], row = y1[i])
    DJ  = Matrix{Float64}(undef, n, n)  # To store the Jacobian
    DIF = similar(DJ)
    Zx  = similar(DJ)
    Zy  = similar(DJ)

    TN_y1 = Matrix{Float64}(undef, n, N)
    TN_y2 = Matrix{Float64}(undef, n, N)

    TNL = Matrix{Float64}(undef, N, n)
    TNR = Matrix{Float64}(undef, n, N)
    A   = Matrix{Float64}(undef, n, n)    
    Tmp = Matrix{Float64}(undef, N, n)
    dw  = Vector{Float64}(undef, n)
    dfy = Vector{Float64}(undef, n)

    I₁ = Matrix{Float64}(undef, N, N)
    I₂ = Matrix{Float64}(undef, N, N)
    qw = Vector{Float64}(undef, n)

    #-------------------------------------------
    # small helpers (no allocations)
    @inline function fill_meshgrid!(T1, T2, y1, y2)
        @inbounds for j in eachindex(y2)
             for i in eachindex(y1)
                T1[i, j] = y1[i]
                T2[i, j] = y2[j]
            end
        end
        return nothing
    end
    #-------------------------------------------
    #A vector of Bool, initialized to true for all
    #indices from 1:Lₚ. They will be updated as 
    #False  for singular points. (and ofcourse the
    #points left are near singular, which are true)
    NSI = trues(Lₚ)

    @inbounds for j in 1:Np
        @inbounds for k in 1:M
            NSI[dp.pthgo[k]+j-1] = false
        end
    end

    @inbounds for j in dp.pthgo[M+1]:bdnsp-1
        NSI[j] = false
    end

    #Now the big for loop for near singular points
    @inbounds for i in 1:Lₚ
        if NSI[i] == true

            k = dp.prepts[2,i]

            # (prepts index) -> `ll` (column in invpts)
            ll = if i < bdnsp
                # interior near–singular to patch k
                i - k * Np
            else
                # boundary near–singulars
                i - M * Np - nbd
            end

            α₁ = dp.invpts[1, ll]
            α₂ = dp.invpts[2, ll]

            if 1 ≤ α₁
                B1 = -winv(p, (α₁ - 1) / (α₁ + 1))
            elseif α₁ ≤ -1
                A1 = -winv(p, (-1 - α₁) / (1 - α₁))
            end

            if 1 ≤ α₂
                B2 = -winv(p, (α₂ - 1) / (α₂ + 1))
            elseif α₂ ≤ -1
                A2 = -winv(p, (-1 - α₂) / (1 - α₂))
            end

            #Not a boundary patch
            if k in knbd
                if 1 ≤ α₁
                    wfunc!(d1, p, z1; α=-B1, β=α₁ + 1.0)
                    @. y1 = α₁ - d1
                    @. y1tmp = 1 - y1
                    if 1 ≤ α₂
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* dw(B1*z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1; α=B1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= B1 * B2 * (α₁ + 1) * (α₂ + 1) / 4

                    elseif -1 < α₂ && α₂ < 1
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(B1*z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1; α=B1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(z1)'.*fw'
                        dwfunc!(dw, p, z1)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        #Left weights same as before
                        #Right weights: TN(y2).*dw(z2)'.*fw'
                        ChebyTN!(TN_y2, N, y2)
                        dwfunc!(dw, p, z2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = B1*(α₁ + 1)*( (α₂ + 1)*I₁ + (1 - α₂)*I₂ )/4

                    elseif α₂ ≤ -1
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(B1*z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1; α=B1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= B1 * A2 * (α₁ + 1) * (1 - α₂) / 4
                    end
                elseif -1 < α₁ && α₁ < 1
                    if 1 ≤ α₂
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d1, p, z2; α=-1.0, β=α₁ - 1.0)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z2) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = B2*(α₂ + 1)*( (α₁ + 1)*I₁ + (1 - α₁)*I₂ )/4

                    elseif α₂ ≤ -1
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d1, p, z2; α=-1.0, β=α₁ - 1.0)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z2) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = A2*(1 - α₂)*( (α₁ + 1)*I₁ + (1 - α₁)*I₂ )/4
                    else
                        error("Fatal error! mapinv not working correctly.")
                    end

                elseif α₁ ≤ -1
                    wfunc!(d1, p, z2; α=-A1, β=α₁ - 1.0)
                    @. y1 = α₁ - d1
                    @. y1tmp = 1 - y1
                    if 1 ≤ α₂
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* dw(A1*z2) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2; α=A1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= A1 * B2 * (1 - α₁) * (α₂ + 1) / 4

                    elseif -1 < α₂ && α₂ < 1
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(A1*z2) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2; α=A1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(z1)'.*fw'
                        dwfunc!(dw, p, z1)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        #Left weights same as before
                        #Right weights: TN(y2).*dw(z2)'.*fw'
                        ChebyTN!(TN_y2, N, y2)
                        dwfunc!(dw, p, z2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = A1*(1 - α₁)*( (α₂ + 1)*I₁ + (1 - α₂)*I₂ )/4
                    elseif α₂ ≤ -1
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2;

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(A1*z2) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2; α=A1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= A1 * A2 * (1 - α₁) * (1 - α₂) / 4
                    end
                end
            else
                #Near singular for boundary patches.
                if 1 ≤ α₁ 
                    #α₁ = 1
                    wfunc!(d1, p, z1; α=-1.0, β=2.0)
                    @. y1 = 1 - d1
                    if 1 ≤ α₂
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, d1, s)   # dfy = dfac(k, 1-y1) = dfac(k, d1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* dw(z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1; α=1.0)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= B2 * (α₂ + 1) / 2

                    elseif -1 < α₂ && α₂ < 1
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, d1, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(z1)'.*fw'
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        #Left weights same as before
                        #Right weights: TN(y2).*dw(z2)'.*fw'
                        ChebyTN!(TN_y2, N, y2)
                        dwfunc!(dw, p, z2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = ((α₂ + 1) * I₁ + (1 - α₂) * I₂) / 2

                    elseif α₂ ≤ -1
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, 1.0, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, d1, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*=  A2 * (1 - α₂) / 2
                    end
                elseif -1 < α₁ && α₁ < 1
                    if 1 ≤ α₂
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d1, p, z; α=1.0, β=(α₁ - 1.0) / 2)
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        wfunc!(y1tmp, p, z; α=-1.0, β=(1.0 - α₁) / 2)

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = B2*(α₂ + 1)*( (α₁ + 1)*I₁ + (1 - α₁) * I₂ )/4

                    elseif α₂ ≤ -1
                        wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        @. y1tmp = 1 - y1

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw(z1) .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z1)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d1, p, z; α=1.0, β=(α₁ - 1.0) / 2)
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y1 = α₁ - d1
                        @. y2 = α₂ - d2 
                        wfunc!(y1tmp, p, z; α=-1.0, β=(1.0 - α₁) / 2)

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        # Left weights: fw .* dw .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = A2*(1 - α₂)*( (α₁ + 1) * I₁ + (1 - α₁) * I₂ )/4
                    else
                        error("Fatal error! mapinv not working correctly.")
                    end

                elseif α₁ ≤ -1
                    A1 = 1 - winv(p, 2 * (-1 - α₁) / (1 - α₁))
                    wfunc!(d1, p, z2; α=-A1, β=(α₁ - 1.0) / 2, γ=1.0)
                    @. y1 = α₁ - d1
                    if 1 ≤ α₂
                        wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        wfunc!(y1tmp, p, z2; α=A1, β=(1.0 - α₁) / 2, γ=-1.0)

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* dw .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2; α=-A1, β=1.0, γ=1.0)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(B2*z1)'.*fw'
                        dwfunc!(dw, p, z1; α=B2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*=  A1 * B2 * (1-α₁)*(α₂+1)/8

                    elseif -1 < α₂ && α₂ < 1
                        wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
                        @. y2 = α₂ - d2

                        wfunc!(y1tmp, p, z2; α=A1, β=(1.0 - α₁) / 2, γ=-1.0)

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* dw .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2; α=-A1, β=1.0, γ=1.0)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(z1)'.*fw'
                        dwfunc!(dw, p, z1)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₁, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        #------------2ⁿᵈ half ------------
                        wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        #Left weights same as before
                        #Right weights: TN(y2).*dw(z2)'.*fw'
                        ChebyTN!(TN_y2, N, y2)
                        dwfunc!(dw, p, z2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I₂, Tmp, TNR)  # I   = Tmp * TNR (N×N)

                        @. I = A1*(1-α₁)*( (α₂+1)*I₁ + (1-α₂)*I₂ )/8
                        
                    elseif α₂ ≤ -1
                        wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
                        @. y2 = α₂ - d2

                        wfunc!(y1tmp, p, z2; α=A1, β=(1.0 - α₁) / 2, γ=-1.0)

                        fill_meshgrid!(t1, t2, y1, y2)

                        diff_map!(DIF, Zx, Zy, DJ, d, α₁, α₂, t1, t2, d1, d2, k)

                        dfunc!(dfy, d, k, y1tmp, s)   # dfy = dfac(k, 1-y1)

                        ChebyTN!(TN_y1, N, y1)
                        ChebyTN!(TN_y2, N, y2)

                        #Left weights: fw .* dw .* dfy .* TN(y1)'
                        dwfunc!(dw, p, z2; α=-A1, β=1.0, γ=1.0)
                        @. wdf = fw * dw * dfy
                        @. TNL = TN_y1' * wdf'

                        # Right weights: TN(y2).*dw(A2*z2)'.*fw'
                        dwfunc!(dw, p, z2; α=A2)
                        @. TNR = TN_y2 * dw * fw   # n×N scaled row-wise
                        # Middle terms together
                        @. A = abs(DIF)^(-2s) * DJ   

                        mul!(Tmp, TNL, A)   # Tmp = TNL * A   (N×n)
                        mul!(I, Tmp, TNR)   # I   = Tmp * TNR (N×N)

                        I .*= A1*A2*(1-α₁)*(1-α₂)/8

                    end
                end
            end

            #Allocation free version of IntS[:,i] = reshape(I,(Np,1));
            @inbounds for m in 1:Np
                IntS[m, i] = I[m]
            end

        end
    end
  
    return IntS
end
