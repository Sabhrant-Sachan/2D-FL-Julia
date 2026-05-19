"""
Boundary-domain properties for SLP boundary equations,
SLP stands for Single layer potential. 

This is a boundary-only analogue of `domprop`. Patches
are curves in 2D space, which were formed due to patching
structure in the interioir. 

Target boundary points are indexed globally as

    i = 1, ..., Mbd*N

where

    kI = cld(i, N)              # boundary-patch block index
    k  = kd[kI]                 # actual domain patch index
    j  = i - (kI - 1)*N         # local Chebyshev index on that boundary patch

The target point is

    γ_k(cospi((2j - 1)/(2N))).

I store the real space coordinates in `tgtpts`. The patch/index data can
always be recovered from the column index `i` as shown above.

`prepts` has two rows:

    prepts[1, col] = boundary target index i
    prepts[2, col] = integration boundary patch k

For each boundary integration patch kd[kI], the precomputation layout is

    N regular boundary targets on that patch,
    followed by nspsz[kI] near-singular boundary targets.

`hmap` maps

    packkey(boundary target index i, integration patch k) -> column in prepts

for near-singular pairs.

`projpts` stores the projected parameter values for near-singular pairs only.
There is no `hmaproj`: if `col` is the column in `prepts` for a near-singular
point in boundary block `kI`, then its column in `projpts` is

    ll = col - kI*N

because each previous boundary block contributes exactly N regular points.
"""
mutable struct dompropbd
    N::Int
    kd::Vector{Int}
    del::Float64

    # lookups: packkey(boundary target index, integration patch k) -> column in prepts
    hmap::Dict{UInt64,Int}

    # 2 × (Mbd*N), physical boundary target coordinates
    tgtpts::Matrix{Float64}

    # 2 × (Mbd*N + Tnsp), row 1 = target index, row 2 = integration patch
    prepts::Matrix{Int}

    # length Tnsp, projected local parameter values on integration patch
    projpts::Vector{Float64}

    # length Mbd, start column of each boundary integration patch block in prepts
    pthgo::Vector{Int}

    # length Mbd, near-singular count for each boundary integration patch
    # nspsz[k₀] = number of boundary target points near-singular to 
    # integration boundary patch kd[k₀]
    nspsz::Vector{Int}

    function dompropbd(N::Int, del::Float64, dom::D) where {D<:abstractdomain}

        kd = dom.kd

        Mbd = length(kd)
        Nd = Mbd * N

        # ------------------------------------------------------------
        # 1. Build boundary target points.
        # ------------------------------------------------------------
        z = Vector{Float64}(undef, N)

        @inbounds for j in 1:N
            z[j] = cospi((2*j - 1) / (2N))
        end

        tgtpts = Matrix{Float64}(undef, 2, Nd)

        x = Vector{Float64}(undef, N)
        y = Vector{Float64}(undef, N)

        @inbounds for k₀ in 1:Mbd
            k = kd[k₀]
            ak = (k₀ - 1) * N
            gamx!(x, dom, z, k)
            gamy!(y, dom, z, k)

            for j in 1:N
                col = ak + j
                tgtpts[1, col] = x[j]
                tgtpts[2, col] = y[j]
            end
        end

        # ------------------------------------------------------------
        # 2. Count near-singular boundary targets for each integration
        #    boundary patch, so that I can preallocate exactly.
        # ------------------------------------------------------------
        nspsz = zeros(Int, Mbd)
        pthgo = zeros(Int, Mbd)

        Qbd = dom.Qptsbd
        Qe = similar(Qbd)

        @inbounds for k₀ in 1:Mbd
            @views V1 = Qe[:, k₀]
            @views V2 = Qbd[:, k₀]
            extendqua!(V1, V2, del)
        end

        insidebdy = Vector{Bool}(undef, Nd)

        @inbounds for k₀ in 1:Mbd
            k = kd[k₀]

            P1x, P1y = Qe[1, k₀], Qe[2, k₀]
            P2x, P2y = Qe[3, k₀], Qe[4, k₀]
            P3x, P3y = Qe[5, k₀], Qe[6, k₀]
            P4x, P4y = Qe[7, k₀], Qe[8, k₀]

            ptinqua!(insidebdy, tgtpts, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)

            cnt = 0

            for i in 1:Nd
                tgt_k₀ = cld(i, N)
                tgt_k = kd[tgt_k₀]

                if insidebdy[i] && tgt_k != k
                    if isnspbd(dom, tgtpts[1, i], tgtpts[2, i], k, del)
                        cnt += 1
                    end
                end
            end

            nspsz[k₀] = cnt
        end

        Tnsp = sum(nspsz)

        # ------------------------------------------------------------
        # 3. Allocate and fill prepts, projpts, hmap directly.
        # ------------------------------------------------------------
        prepts = Matrix{Int}(undef, 2, Nd + Tnsp)
        projpts = Vector{Float64}(undef, Tnsp)

        hmap = Dict{UInt64,Int}()
        sizehint!(hmap, Tnsp)

        col = 1
        projcol = 1

        @inbounds for k₀ in 1:Mbd
            k = kd[k₀]

            pthgo[k₀] = col

            # Regular boundary targets on this boundary patch.
            for j in 1:N
                tgt_i = (k₀ - 1) * N + j

                prepts[1, col] = tgt_i
                prepts[2, col] = k

                col += 1
            end

            # Now find and directly append the near-singular boundary targets
            # with respect to integration patch k.
            P1x, P1y = Qe[1, k₀], Qe[2, k₀]
            P2x, P2y = Qe[3, k₀], Qe[4, k₀]
            P3x, P3y = Qe[5, k₀], Qe[6, k₀]
            P4x, P4y = Qe[7, k₀], Qe[8, k₀]

            ptinqua!(insidebdy, tgtpts, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)

            for tgt_i in 1:Nd
                tgt_k₀ = cld(tgt_i, N)
                tgt_k = kd[tgt_k₀]

                if insidebdy[tgt_i] && tgt_k != k
                    if isnspbd(dom, tgtpts[1, tgt_i], tgtpts[2, tgt_i], k, del)

                        prepts[1, col] = tgt_i
                        prepts[2, col] = k

                        hmap[packkey(tgt_i, k)] = col

                        x1 = tgtpts[1, tgt_i]
                        x2 = tgtpts[2, tgt_i]

                        τ = projbd(dom, x1, x2, k)

                        # Snap endpoint values for corner-type roundoff.
                        τ = abs(τ - 1.0) < 1e-14 ? 1.0 :
                            abs(τ + 1.0) < 1e-14 ? -1.0 : τ

                        projpts[projcol] = τ

                        col += 1
                        projcol += 1
                    end
                end
            end
        end

        return new(N, kd, del, hmap, tgtpts, prepts, projpts, pthgo, nspsz)
    end
end

# These are preocmputations within the boundary, not for interioir points.
# They will be used for solving the densities βⱼ, 1 <= j <= nh on the boundary.  
function precompsSLP(d::D, dpbd::dompropbd, p::Int; 
    n::Int = 128)::Matrix{Float64} where {D<:abstractdomain}

    kd = dpbd.kd
    N = dpbd.N
    Mbd = length(kd)
    Nd = Mbd * N
    Tnsp = sum(dpbd.nspsz)

    Lₚ = Nd + Tnsp

    # ------------------------------------------------------------
    # Chebyshev boundary target nodes.
    # zp[j]  = cospi((2j - 1)/(2N))
    # zp1[j] = zp[j] + 1
    # zp2[j] = zp[j] - 1
    # So:
    #   zp1[j] = t + 1
    #  -zp2[j] = 1 - t
    # ------------------------------------------------------------
    zp  = Vector{Float64}(undef, N)
    zp1 = Vector{Float64}(undef, N)
    zp2 = Vector{Float64}(undef, N)

    @inbounds for j in 1:N
        c = π * (2 * j - 1) / (2 * N)
        zp[j]  = cos(c)
        zp1[j] = 2 * cos(c / 2)^2       # 1 + cos(c)
        zp2[j] = -2 * sin(c / 2)^2      # cos(c) - 1
    end

    # ------------------------------------------------------------
    # Quadrature nodes for singular / near-singular splitting.
    #   z1 = cospi((2*(0:n-1)+1)/(4*n)).^2
    #   z2 = sinpi((2*(0:n-1)+1)/(4*n)).^2
    # ------------------------------------------------------------
    z1 = Vector{Float64}(undef, n)
    z2 = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        c = π * (2 * i - 1) / (2 * n)
        z1[i] = cos(c / 2)^2
        z2[i] = sin(c / 2)^2
    end

    fw = getF1W(n)

    # ------------------------------------------------------------
    # Scratch arrays.
    # ------------------------------------------------------------
    d1 = Vector{Float64}(undef, n)
    d2 = Vector{Float64}(undef, n)

    y1 = Vector{Float64}(undef, n)
    y2 = Vector{Float64}(undef, n)

    Iv1 = Vector{Float64}(undef, n)
    Iv2 = Vector{Float64}(undef, n)

    I1 = Vector{Float64}(undef, N)
    I2 = Vector{Float64}(undef, N)

    # For diff_map! singular distance fixup.
    zx  = Vector{Float64}(undef, n)
    zy  = Vector{Float64}(undef, n)
    DJℓ = Vector{Float64}(undef, n)
    DIF = Vector{Float64}(undef, n)

    # Boundary coordinates and derivatives.
    TNy = Matrix{Float64}(undef, n, N)

    IntS = Matrix{Float64}(undef, N, Lₚ)

    gamk = Matrix{Float64}(undef, 2, n)
    dgamk= Matrix{Float64}(undef, 2, n)

    # ------------------------------------------------------------
    # Store w(-z1), w(-z2), w'(z1), w'(z2).
    # ------------------------------------------------------------
    wmz₁ = Vector{Float64}(undef, n)
    wmz₂ = Vector{Float64}(undef, n)
    dwz₁ = Vector{Float64}(undef, n)
    dwz₂ = Vector{Float64}(undef, n)

    wfunc!(wmz₁, p, z1; α = -1.0)
    wfunc!(wmz₂, p, z2; α = -1.0)

    dwfunc!(dwz₁, p, z1)
    dwfunc!(dwz₂, p, z2)

    inv2π = 1 / (2π)

    # Map actual boundary patch label k -> boundary block index kI (1 to Mbd).
    # Needed because dpbd.projpts is packed by boundary integration block.
    kd_to_kI = Dict{Int,Int}()
    sizehint!(kd_to_kI, Mbd)

    @inbounds for kI in 1:Mbd
        kd_to_kI[kd[kI]] = kI
    end

    # ------------------------------------------------------------
    # Main precomputation loop.
    # ------------------------------------------------------------
    @inbounds for i in 1:Lₚ

        # col is the boundary target index.
        col = dpbd.prepts[1, i]

        # k is the actual integration boundary patch.
        k = dpbd.prepts[2, i]

        # Recover the target boundary block and local Chebyshev index.
        tgt_kI = cld(col, N)
        ℓ = kd[tgt_kI]
        j = col - (tgt_kI - 1) * N

        if ℓ == k
            # Singular integration case.
            t = zp[j]

            @. d1 = zp1[j] * wmz₁
            @. d2 = zp2[j] * wmz₂

            @. y1 = t - d1
            @. y2 = t - d2

            # First side: y1.
            diff_map!(DIF, zx, zy, DJℓ, d, 1.0, t, 1.0, y1, 0.0, d1, k)

            ChebyTN!(TNy, N, y1)

            @. Iv1 = fw * log(DIF) * inv2π * dwz₁ * DJℓ

            # I1 = TNy' * Iv1
            mul!(I1, transpose(TNy), Iv1)

            # Second side: y2.
            diff_map!(DIF, zx, zy, DJℓ, d, 1.0, t, 1.0, y2, 0.0, d2, k)

            ChebyTN!(TNy, N, y2)

            @. Iv2 = fw * log(DIF) * inv2π * dwz₂ * DJℓ

            # I2 = TNy' * Iv2
            mul!(I2, transpose(TNy), Iv2)

            @inbounds for row in 1:N
                IntS[row, i] = (zp1[j] * I1[row] - zp2[j] * I2[row]) / 2
            end

        else
            # Near-singular integration case.
            x1 = dpbd.tgtpts[1, col]
            x2 = dpbd.tgtpts[2, col]

            # kI is the boundary block index of the integration patch k.
            kI = kd_to_kI[k]

            # Column in projpts because prepts is packed by boundary 
            # integration block. Each previous/current block contributes
            # N regular columns that are not represented in projpts.
            # Therefore for a near-singular prepts column i:
            ll = i - kI * N

            t₀ = dpbd.projpts[ll]
            if t₀ == 1.0
                # So only the y1-side contributes.
                @. y1 = 1.0 - 2.0 * wmz₁

                gam!(gamk, d, y1, k)
                dgam!(dgamk, d, y1, k)

                @inbounds for r in 1:n
                    dx = x1 - gamk[1, r]
                    dy = x2 - gamk[2, r]

                    DIF[r] = hypot(dx, dy)
                    DJℓ[r] = hypot(dgamk[1, r], dgamk[2, r])
                end

                ChebyTN!(TNy, N, y1)

                @. Iv1 = fw * log(DIF) * inv2π * dwz₁ * DJℓ

                mul!(I1, transpose(TNy), Iv1)

                @inbounds for row in 1:N
                    IntS[row, i] = I1[row]
                end

            elseif t₀ == -1.0
                # So only the y2-side contributes.
                @. y2 = -1.0 + 2.0 * wmz₂

                gam!(gamk, d, y2, k)
                dgam!(dgamk, d, y2, k)

                @inbounds for r in 1:n
                    dx = x1 - gamk[1, r]
                    dy = x2 - gamk[2, r]

                    DIF[r] = hypot(dx, dy)
                    DJℓ[r] = hypot(dgamk[1, r], dgamk[2, r])
                end

                ChebyTN!(TNy, N, y2)

                @. Iv2 = fw * log(DIF) * inv2π * dwz₂ * DJℓ

                mul!(I2, transpose(TNy), Iv2)

                @inbounds for row in 1:N
                    IntS[row, i] = I2[row]
                end

            else
                # Both sides contribute:
                @. y1 = t₀ - (t₀ + 1.0) * wmz₁
                @. y2 = t₀ + (1.0 - t₀) * wmz₂

                # First side: y1.
                gam!(gamk, d, y1, k)
                dgam!(dgamk, d, y1, k)

                @inbounds for r in 1:n
                    dx = x1 - gamk[1, r]
                    dy = x2 - gamk[2, r]

                    DIF[r] = hypot(dx, dy)
                    DJℓ[r] = hypot(dgamk[1, r], dgamk[2, r])
                end

                ChebyTN!(TNy, N, y1)

                @. Iv1 = fw * log(DIF) * inv2π * dwz₁ * DJℓ

                mul!(I1, transpose(TNy), Iv1)

                # Second side: y2.
                gam!(gamk, d, y2, k)
                dgam!(dgamk, d, y2, k)

                @inbounds for r in 1:n
                    dx = x1 - gamk[1, r]
                    dy = x2 - gamk[2, r]

                    DIF[r] = hypot(dx, dy)
                    DJℓ[r] = hypot(dgamk[1, r], dgamk[2, r])
                end

                ChebyTN!(TNy, N, y2)

                @. Iv2 = fw * log(DIF) * inv2π * dwz₂ * DJℓ

                mul!(I2, transpose(TNy), Iv2)

                @inbounds for row in 1:N
                    IntS[row, i] = ((t₀ + 1.0) * I1[row] + (1.0 - t₀) * I2[row]) / 2.0
                end
            end
        end
    end

    return IntS
end

function SLPbeta(d::D, dpbd::dompropbd)::Matrix{Float64} where {D<:abstractdomain}

    Mbd, p = length(d.kd), 6

    kd = dpbd.kd

    N  = dpbd.N 

    Nd = Mbd * N

    IntS = precompsSLP(d, dpbd, p)

    # Initialize the matrix A
    # For overdetermined systems, Julia backslash 
    # gives least-squares solution. 
    A = zeros(Float64, Nd + 1, Nd)

    gamvals = Vector{Float64}(undef, N)
    coeffs = Matrix{Float64}(undef, N, N)

    gamvals[1] = 1.0
    @inbounds for k in 2:N
        gamvals[k] = 2.0
    end

    @inbounds for j in 1:N
        c = π * (2j - 1) / (2N)
        @inbounds for k in 1:N
            coeffs[k, j] = gamvals[k] * cos(c * (k - 1)) / N
        end
    end

    # Define nr, number of nodes for regular integral, with weights
    nr = 32
    fwr = getF1W(nr)

    zr = Vector{Float64}(undef, nr)
    Ir = Vector{Float64}(undef, nr)

    @inbounds for i in 1:nr
        ph = π * (2 * i - 1) / (2 * nr)
        zr[i] = cos(ph)
    end

    # Store idct values
    t₁ = Vector{Float64}(undef, N)
    t₂ = Vector{Float64}(undef, nr)
    KER = Vector{Float64}(undef, nr)
    idctrg = Matrix{Float64}(undef, nr, N)

    @inbounds for i in 1:nr
        t₂[i] = π * (2 * i - 1) / (2 * nr)
    end

    @inline function gs(t::Float64, N::Int)::Float64
        th = 0.5 * t                 
        return sin(N * th) * cos((N - 1) * th) / sin(th)
    end

    @inbounds for j in 1:N
        t₁[j] = π * (2 * j - 1) / (2 * N)
        for i in 1:nr
            u = t₁[j] + t₂[i]
            v = t₁[j] - t₂[i]
            idctrg[i, j] = (gs(u, N) + gs(v, N) - 1.0) / N
        end
    end

    #For boundary
    gamk = Matrix{Float64}(undef, 2, nr)
    dgamk = Matrix{Float64}(undef, 2, nr)
    DJ   = Vector{Float64}(undef, nr) 

    inv_2pi = 1 / (2 * pi)
    
    #------------ PART I: compute matrix A ------------
    @inbounds for kk in 1:Mbd

        k = kd[kk]

        gam!(gamk, d, zr, k)
        gamp!(dgamk, d, zr, k)

        @views v = A[:, 1+N*(kk-1):N*kk]

        @inbounds for jj in 1:nr
            DJ[jj] = hypot(dgamk[1,jj],dgamk[2,jj])
        end 

        @inbounds for row in 1:Nd

            kI = cld(row, N)              
            ℓ  = kd[kI]                 
            j  = row - (kI - 1)*N
            if ℓ == k
                #Singular integration
                @views SI = IntS[:, dpbd.pthgo[kI]+j-1]
                @inbounds for jj in 1:N
                    @views c = coeffs[:, jj]
                    v[row, jj] = dot(c, SI)
                end
            else
                Ikey = packkey(row, k)
                col = get(dpbd.hmap, Ikey, 0)
                if col != 0

                    @views NSI = IntS[:, col]
                    @inbounds for jj in 1:N
                        @views c = coeffs[:, jj]
                        v[row, jj] = dot(c, NSI)
                    end

                else

                    x1 = dpbd.tgtpts[1, row]
                    x2 = dpbd.tgtpts[2, row]
                    @inbounds for jj in 1:nr
                        dx1 = x1 - gamk[1, jj]
                        dx2 = x2 - gamk[2, jj]
                        KER[jj] = inv_2pi * log(hypot(dx1, dx2))
                    end

                    @inbounds for jj in 1:N
                        @views Ur = idctrg[:, jj]
                        @. Ir = KER * Ur * DJ 
                        v[row, jj] = dot(fwr, Ir)
                    end
                    
                end

            end

        end

        @inbounds for jj in 1:N
            @views Ur = idctrg[:, jj] 
            @. Ir = Ur * DJ 
            v[Nd+1, jj] = dot(fwr, Ir)
        end
    end

    #------------ PART II: compute RHS b ------------
    beta = Matrix{Float64}(undef, Nd, d.nh)
    bv = Vector{Float64}(undef, Nd+1)

    @inbounds for H in 1:d.nh
        fill!(bv, 0.0)

        @inbounds for i in 1:Mbd
            if bdno(d, kd[i]) == H + 1
                aN = N*(i-1)
                @inbounds for jj in 1:N
                    bv[aN + jj] = 1
                end
            end
        end

        @views beta[:, H] .= A \ bv

    end

    return beta
end

function SLPeval(d::D, dp::domprop, β::Matrix{Float64}, IV;
    p::Int=6)::Matrix{Float64} where {D<:abstractdomain}

    # --------------- Sizes ---------------
    N = dp.N
    Np = N^2
    M = d.Npat
    Mbd = length(d.kd)

    Ni = M * Np
    Nb = Mbd * N
    Nt = Ni + Nb

    nh = size(β, 2)

    @assert size(β, 1) == Mbd * N

    # ------------------------------------------------------------
    # Unpack reusable arrays from IV
    # ------------------------------------------------------------
    (; IVbd, IVbdt1, IVbdt3) = IV

    (; nr_bdy, ns_near, fwr_bdy, fw_near, zr_bdy) = IVbd

    (; y₁, y₂, wmz₁, wmz₂, dwz₁, dwz₂,
        gamks₁, gamks₂, gamperpks₁, gamperpks₂) = IVbdt1

    (; Tbdt₀, Tbdt1, Tbdt2) = IVbdt3

    inv2π = 1.0 / (2π)

    BD_FAR = UInt8(0)
    BD_NEAR = UInt8(1)

    # Local DCT scratch for β.
    βv = Vector{Float64}(undef, N)
    βfv_bdy = Vector{Float64}(undef, nr_bdy)

    p_dct2_N = FFTW.plan_r2r!(βv, FFTW.REDFT10; flags=FFTW.MEASURE)
    p_dct3_bdy = FFTW.plan_r2r!(βfv_bdy, FFTW.REDFT01; flags=FFTW.MEASURE)

    # --------------- Output ---------------
    bZ = zeros(Float64, Nt, nh)

    # ------------------------------------------------------------
    # Boundary density storage
    #
    # βcoef    : Chebyshev coefficients on each boundary patch
    # β_bdy    : β values on regular boundary quadrature grid zr_bdy
    # β_neart1 : β values on fixed endpoint-smoothed grid near t = 1
    # β_neart2 : β values on fixed endpoint-smoothed grid near t = -1
    # ------------------------------------------------------------
    βcoef    = Matrix{Float64}(undef, Mbd * N, nh)
    β_bdy    = Matrix{Float64}(undef, Mbd * nr_bdy, nh)
    β_neart1 = Matrix{Float64}(undef, Mbd * ns_near, nh)
    β_neart2 = Matrix{Float64}(undef, Mbd * ns_near, nh)

    @inbounds for H in 1:nh
        @inbounds for ll in 1:Mbd

            aN    = (ll - 1) * N
            ibd   = (ll - 1) * nr_bdy
            inear = (ll - 1) * ns_near

            @inbounds for j in 1:N
                βv[j] = β[aN + j, H]
            end

            mul!(βv, p_dct2_N, βv)

            βv .*= 1.0 / N
            βv[1] /= 2.0

            @inbounds for j in 1:N
                βcoef[aN + j, H] = βv[j]
            end

            mul!(y₁, Tbdt1, βv)
            mul!(y₂, Tbdt2, βv)

            @inbounds for j in 1:ns_near
                β_neart1[inear + j, H] = y₁[j]
                β_neart2[inear + j, H] = y₂[j]
            end

            βv[1] *= 2.0

            fill!(βfv_bdy, 0.0)

            @inbounds for j in 1:N
                βfv_bdy[j] = βv[j]
            end

            mul!(βfv_bdy, p_dct3_bdy, βfv_bdy)

            βfv_bdy .*= 0.5

            @inbounds for j in 1:nr_bdy
                β_bdy[ibd + j, H] = βfv_bdy[j]
            end
        end
    end

    # ------------------------------------------------------------
    # Regular boundary geometry cache
    #
    # gamkbdy_all[:, rb] = γ_k(zr_bdy)
    # gampkn_all[rb]     = ‖γ'_k(zr_bdy)‖ = ‖γᵖᵉʳᵖ_k(zr_bdy)‖
    # ------------------------------------------------------------
    gamkbdy = Matrix{Float64}(undef, 2, nr_bdy)
    gampkbdy = Matrix{Float64}(undef, 2, nr_bdy)

    gamkbdy_all = Matrix{Float64}(undef, 2, nr_bdy * Mbd)
    gampkn_all = Vector{Float64}(undef, nr_bdy * Mbd)

    @inbounds for ll in 1:Mbd
        k = d.kd[ll]

        gam!(gamkbdy, d, zr_bdy, k)
        gamp!(gampkbdy, d, zr_bdy, k)

        @inbounds for j in 1:nr_bdy
            r = (ll - 1) * nr_bdy + j

            gamkbdy_all[1, r] = gamkbdy[1, j]
            gamkbdy_all[2, r] = gamkbdy[2, j]

            gampkn_all[r] = hypot(gampkbdy[1, j], gampkbdy[2, j])
        end
    end

    # ----------- Small scratch arrays ----------- 
    βdyn1 = Vector{Float64}(undef, ns_near)
    βdyn2 = Vector{Float64}(undef, ns_near)

    # ============================================================
    # Main SLP evaluation
    # ============================================================
    @inbounds for H in 1:nh

        jnear = 0
        jintp = 0

        @inbounds for row in 1:Ni

            x1 = dp.tgtpts[1, row]
            x2 = dp.tgtpts[2, row]

            SLP = 0.0

            # ====================================================
            # Mode 0: far from every boundary panel
            # ====================================================
            if dp.bdmode[row] == BD_FAR

                @inbounds for ll in 1:Mbd

                    rb = ((ll-1)*nr_bdy+1):(ll*nr_bdy)

                    @views gamkbdy = gamkbdy_all[:, rb]
                    @views gampkn = gampkn_all[rb]
                    @views βv_bdy = β_bdy[rb, H]

                    acc = 0.0

                    @inbounds for j in 1:nr_bdy
                        dx = x1 - gamkbdy[1, j]
                        dy = x2 - gamkbdy[2, j]

                        acc += (fwr_bdy[j] * log(hypot(dx, dy)) *
                               inv2π * gampkn[j] * βv_bdy[j])
                    end

                    SLP += acc
                end

            elseif dp.bdmode[row] == BD_NEAR
                # ====================================================
                # Mode 1: moderately near some boundary panels
                # ====================================================

                jnear += 1

                qnear1 = dp.bdnearptr[jnear]
                qnear2 = dp.bdnearptr[jnear+1] - 1

                @inbounds for ll in 1:Mbd

                    k = d.kd[ll]

                    near_flag = false
                    t₀ = 0.0

                    @inbounds for q in qnear1:qnear2
                        if dp.bdneark[q] == k
                            near_flag = true
                            t₀ = dp.bdneart[q]
                            break
                        end
                    end

                    if near_flag

                        aN = (ll - 1) * N
                        inear = (ll - 1) * ns_near

                        if t₀ == 1.0

                            @views βn = β_neart1[(inear+1):(inear+ns_near), H]

                            @. y₁ = 1.0 - 2.0 * wmz₁

                            gam!(gamks₁, d, y₁, k)
                            gamp!(gamperpks₁, d, y₁, k)

                            acc = 0.0

                            @inbounds for j in 1:ns_near
                                dx = x1 - gamks₁[1, j]
                                dy = x2 - gamks₁[2, j]

                                gampkn = hypot(gamperpks₁[1, j], gamperpks₁[2, j])

                                acc += (fw_near[j] * log(hypot(dx, dy)) *
                                        inv2π * dwz₁[j] * gampkn * βn[j])
                            end

                            SLP += acc

                        elseif t₀ == -1.0

                            @views βn = β_neart2[(inear+1):(inear+ns_near), H]

                            @. y₂ = -1.0 + 2.0 * wmz₂

                            gam!(gamks₂, d, y₂, k)
                            gamp!(gamperpks₂, d, y₂, k)

                            acc = 0.0

                            @inbounds for j in 1:ns_near
                                dx = x1 - gamks₂[1, j]
                                dy = x2 - gamks₂[2, j]

                                gampkn = hypot(gamperpks₂[1, j], gamperpks₂[2, j])

                                acc += (fw_near[j] * log(hypot(dx, dy)) *
                                       inv2π * dwz₂[j] * gampkn * βn[j])
                            end

                            SLP += acc

                        else

                            @views βc = βcoef[(aN+1):(aN+N), H]

                            # First side.
                            @. y₁ = t₀ - (t₀ + 1.0) * wmz₁

                            gam!(gamks₁, d, y₁, k)
                            gamp!(gamperpks₁, d, y₁, k)

                            ChebyTN!(Tbdt₀, N, y₁)
                            mul!(βdyn1, Tbdt₀, βc)

                            acc1 = 0.0

                            @inbounds for j in 1:ns_near
                                dx = x1 - gamks₁[1, j]
                                dy = x2 - gamks₁[2, j]

                                gampkn = hypot(gamperpks₁[1, j], gamperpks₁[2, j])

                                acc1 += (fw_near[j] * log(hypot(dx, dy)) *
                                        inv2π * dwz₁[j] * gampkn * βdyn1[j])
                            end

                            # Second side.
                            @. y₂ = t₀ + (1.0 - t₀) * wmz₂

                            gam!(gamks₂, d, y₂, k)
                            gamp!(gamperpks₂, d, y₂, k)

                            ChebyTN!(Tbdt₀, N, y₂)
                            mul!(βdyn2, Tbdt₀, βc)

                            acc2 = 0.0

                            @inbounds for j in 1:ns_near
                                dx = x1 - gamks₂[1, j]
                                dy = x2 - gamks₂[2, j]

                                gampkn = hypot(gamperpks₂[1, j], gamperpks₂[2, j])

                                acc2 += (fw_near[j] * log(hypot(dx, dy)) *
                                        inv2π * dwz₂[j] * gampkn * βdyn2[j])
                            end

                            SLP += ((t₀ + 1.0) * acc1 + (1.0 - t₀) * acc2) / 2.0
                        end

                    else

                        rb = ((ll-1)*nr_bdy+1):(ll*nr_bdy)

                        @views gamkbdy = gamkbdy_all[:, rb]
                        @views gampkn = gampkn_all[rb]
                        @views βv_bdy = β_bdy[rb, H]

                        acc = 0.0

                        @inbounds for j in 1:nr_bdy
                            dx = x1 - gamkbdy[1, j]
                            dy = x2 - gamkbdy[2, j]

                            acc += (fwr_bdy[j] * log(hypot(dx, dy)) *
                                   inv2π * gampkn[j] * βv_bdy[j])
                        end

                        SLP += acc
                    end
                end

            else
                # ====================================================
                # Mode 2: very near the boundary
                # For SLP, no normal interpolation is needed.
                # Use smoothed quadrature. Nearby panels are taken from
                # the first auxiliary interior point, m = 2.
                # ====================================================

                jintp += 1

                l = dp.bdclosest[jintp]
                tₛ = dp.bdt[jintp]

                # Use first auxiliary interior point for nearby panel data.
                a_aux = (jintp - 1) * dp.Lᵢₙ + 2

                qaux1 = dp.bdintpptr[a_aux]
                qaux2 = dp.bdintpptr[a_aux+1] - 1

                @inbounds for ll in 1:Mbd

                    k = d.kd[ll]

                    near_flag = false
                    t₀ = 0.0

                    if k == l
                        near_flag = true
                        t₀ = tₛ
                    else
                        @inbounds for q in qaux1:qaux2
                            if dp.bdintpk[q] == k
                                near_flag = true
                                t₀ = dp.bdintpt[q]
                                break
                            end
                        end
                    end

                    if near_flag

                        aN = (ll - 1) * N
                        inear = (ll - 1) * ns_near

                        if t₀ == 1.0

                            @views βn = β_neart1[(inear+1):(inear+ns_near), H]

                            @. y₁ = 1.0 - 2.0 * wmz₁

                            gam!(gamks₁, d, y₁, k)
                            gamp!(gamperpks₁, d, y₁, k)

                            acc = 0.0

                            @inbounds for j in 1:ns_near
                                dx = x1 - gamks₁[1, j]
                                dy = x2 - gamks₁[2, j]

                                gampkn = hypot(gamperpks₁[1, j], gamperpks₁[2, j])

                                acc += (fw_near[j] * log(hypot(dx, dy)) *
                                       inv2π * dwz₁[j] * gampkn * βn[j])
                            end

                            SLP += acc

                        elseif t₀ == -1.0

                            @views βn = β_neart2[(inear+1):(inear+ns_near), H]

                            @. y₂ = -1.0 + 2.0 * wmz₂

                            gam!(gamks₂, d, y₂, k)
                            gamp!(gamperpks₂, d, y₂, k)

                            acc = 0.0

                            @inbounds for j in 1:ns_near
                                dx = x1 - gamks₂[1, j]
                                dy = x2 - gamks₂[2, j]

                                gampkn = hypot(gamperpks₂[1, j], gamperpks₂[2, j])

                                acc += (fw_near[j] * log(hypot(dx, dy)) *
                                       inv2π * dwz₂[j] * gampkn * βn[j])
                            end

                            SLP += acc

                        else

                            @views βc = βcoef[(aN+1):(aN+N), H]

                            # First side.
                            @. y₁ = t₀ - (t₀ + 1.0) * wmz₁

                            gam!(gamks₁, d, y₁, k)
                            gamp!(gamperpks₁, d, y₁, k)

                            ChebyTN!(Tbdt₀, N, y₁)
                            mul!(βdyn1, Tbdt₀, βc)

                            acc1 = 0.0

                            @inbounds for j in 1:ns_near
                                dx = x1 - gamks₁[1, j]
                                dy = x2 - gamks₁[2, j]

                                gampkn = hypot(gamperpks₁[1, j], gamperpks₁[2, j])

                                acc1 += (fw_near[j] * log(hypot(dx, dy)) *
                                        inv2π * dwz₁[j] * gampkn * βdyn1[j])
                            end

                            # Second side.
                            @. y₂ = t₀ + (1.0 - t₀) * wmz₂

                            gam!(gamks₂, d, y₂, k)
                            gamp!(gamperpks₂, d, y₂, k)

                            ChebyTN!(Tbdt₀, N, y₂)
                            mul!(βdyn2, Tbdt₀, βc)

                            acc2 = 0.0

                            @inbounds for j in 1:ns_near
                                dx = x1 - gamks₂[1, j]
                                dy = x2 - gamks₂[2, j]

                                gampkn = hypot(gamperpks₂[1, j], gamperpks₂[2, j])

                                acc2 += (fw_near[j] * log(hypot(dx, dy)) *
                                        inv2π * dwz₂[j] * gampkn * βdyn2[j])
                            end

                            SLP += ((t₀ + 1.0) * acc1 + (1.0 - t₀) * acc2) / 2.0
                        end

                    else

                        rb = ((ll-1)*nr_bdy+1):(ll*nr_bdy)

                        @views gamkbdy = gamkbdy_all[:, rb]
                        @views gampkn = gampkn_all[rb]
                        @views βv_bdy = β_bdy[rb, H]

                        acc = 0.0

                        @inbounds for j in 1:nr_bdy
                            dx = x1 - gamkbdy[1, j]
                            dy = x2 - gamkbdy[2, j]

                            acc += (fwr_bdy[j] * log(hypot(dx, dy)) *
                                   inv2π * gampkn[j] * βv_bdy[j])
                        end

                        SLP += acc
                    end
                end
            end

            bZ[row, H] = SLP
        end

        # ------------------------------------------------------------
        # Boundary values are known.
        # Boundary component H+1 gets value 1, all others stay 0.
        # ------------------------------------------------------------
        @inbounds for ll in 1:Mbd
            if bdno(d, d.kd[ll]) == H + 1
                a = Ni + (ll - 1) * N

                @inbounds for j in 1:N
                    bZ[a+j, H] = 1.0
                end
            end
        end
    end

    return bZ
end

function SLPeval(d::D, dp::domprop, IV; δ::Float64 = 0.1,
     p::Int = 6)::Matrix{Float64} where {D<:abstractdomain}

    #This function evaluates the single layer potentials 
    #for interioir points
    N = dp.N

    dpbd = dompropbd(N, δ, d)

    # β: A rectangular matrix of size (Mbd*N) * nh where
    #    N is number of Chebyshev coefficients per patch and
    #    Mbd is number of boundary patches. The jth column
    #    of beta contains βⱼ function values over Chebyshev mesh.
    β = SLPbeta(d, dpbd)

    return SLPeval(d, dp, β, IV; p = p)
end

function SLPtest(d::annulus, dp::domprop, IV;
    δ::Float64 = 0.1,
    tolβ::Float64 = 5e-10,
    tolSLP::Float64 = 5e-10)

    # This test assumes a circular concentric annulus.
    # Outer boundary: reg 1:4, radius R2
    # Inner boundary: reg 5:8, radius R1

    @assert d.nh == 1 "SLPtest expects an annulus with one hole."

    @assert isapprox(d.R11, d.R12; rtol=0.0, atol=1e-14) "Outer boundary must be circular."
    @assert isapprox(d.R21, d.R22; rtol=0.0, atol=1e-14) "Inner boundary must be circular."

    R2 = d.R11
    R1 = d.R21

    @assert R1 < R2 "Expected inner radius R1 < outer radius R2."

    logratio = log(R1 / R2)

    βinner_exact =  1.0 / (R1 * logratio)
    βouter_exact = -1.0 / (R2 * logratio)

    N   = dp.N
    Np  = N^2
    M   = d.Npat
    Mbd = length(d.kd)

    Ni = M * Np
    Nb = Mbd * N
    Nt = Ni + Nb

    # ------------------------------------------------------------
    # 1. Compute β using the actual solver.
    # ------------------------------------------------------------
    dpbd = dompropbd(N, δ, d)

    β = SLPbeta(d, dpbd)

    @assert size(β, 1) == Mbd * N
    @assert size(β, 2) == d.nh

    # ------------------------------------------------------------
    # 2. Compare computed β against exact constant density.
    # ------------------------------------------------------------
    βexact = Matrix{Float64}(undef, Mbd * N, d.nh)

    @inbounds for ll in 1:Mbd
        k = d.kd[ll]
        reg = d.pths[k].reg

        βval = reg <= 4 ? βouter_exact : βinner_exact

        a = (ll - 1) * N
        for j in 1:N
            βexact[a + j, 1] = βval
        end
    end

    βerr = abs.(β .- βexact)

    βmaxerr, βworst = findmax(βerr)

    βrelerr = βmaxerr / maximum(abs.(βexact))

    # ------------------------------------------------------------
    # 3. Evaluate SLP using computed β.
    # ------------------------------------------------------------
    bZ = SLPeval(d, dp, β, IV)

    # ------------------------------------------------------------
    # 4. Exact SLP:
    #
    #     S[β](x) = log(r/R2) / log(R1/R2)
    #
    # centered at (d.A, d.B).
    # ------------------------------------------------------------
    exact = Vector{Float64}(undef, Nt)

    @inbounds for i in 1:Nt
        x = dp.tgtpts[1, i] - d.A
        y = dp.tgtpts[2, i] - d.B
        r = hypot(x, y)

        exact[i] = log(r / R2) / logratio
    end

    # Boundary exact values:
    # outer boundary = 0, inner boundary = 1.
    @inbounds for ll in 1:Mbd
        k = d.kd[ll]
        reg = d.pths[k].reg

        val = reg <= 4 ? 0.0 : 1.0

        a = Ni + (ll - 1) * N
        for j in 1:N
            exact[a + j] = val
        end
    end

    slperr = abs.(bZ[:, 1] .- exact)

    slpmaxerr, slpworst = findmax(slperr)

    slprelerr = slpmaxerr / maximum(abs.(exact))

    # ------------------------------------------------------------
    # 5. Report.
    # ------------------------------------------------------------
    println("SLP annulus full test")
    println("  R1 inner radius       : ", R1)
    println("  R2 outer radius       : ", R2)
    println("  exact β inner         : ", βinner_exact)
    println("  exact β outer         : ", βouter_exact)

    println()
    if βmaxerr < tolβ
        println("  β solve: Okay!")
    else
        println("  β solve: Bug!")
    end

    println("    max β error         : ", βmaxerr)
    println("    rel β error         : ", βrelerr)
    println("    worst β index       : ", βworst)
    println("    computed β          : ", β[βworst])
    println("    exact β             : ", βexact[βworst])
    println("    β difference        : ", β[βworst] - βexact[βworst])

    println()
    if slpmaxerr < tolSLP
        println("  SLP eval: Okay!")
    else
        println("  SLP eval: Bug!")
    end

    println("    max SLP error       : ", slpmaxerr)
    println("    rel SLP error       : ", slprelerr)
    println("    worst target index  : ", slpworst)
    println("    computed SLP        : ", bZ[slpworst, 1])
    println("    exact SLP           : ", exact[slpworst])
    println("    SLP difference      : ", bZ[slpworst, 1] - exact[slpworst])

    return β, βexact, bZ, exact, βerr, slperr
end