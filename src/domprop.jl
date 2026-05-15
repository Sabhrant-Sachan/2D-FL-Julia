@inline packkey(i::Int, k::Int) = (UInt64(i) << 32) | UInt64(k)

mutable struct domprop
    N::Int
    del::Float64

    # Boundary DLP classification cutoffs.
    # del_near:
    #     if dist(xᵢ, ∂Ωₖ) <= del_near for some boundary panel k,
    #     then xᵢ is considered boundary-near.
    # del_intp:
    #     if dist(xᵢ, ∂Ωₖ) < del_intp for some boundary panel k,
    #     then xᵢ is considered very-near and will be evaluated by
    #     full-boundary normal interpolation.
    del_near::Float64
    del_intp::Float64

    # hmap[packkey(i,k)] = colp
    # i = global target column in tgtpts
    # k = actual volume integration patch index
    # colp = column in prepts
    hmap::Dict{UInt64,Int}

    # Nb is number of boundary targets
    # 2 × (Ni + Nb)
    # First row is x coordinates, second y coordinate
    # columns 1:Ni are interior targets
    # columns Ni+1:Ni+Nb are boundary targets
    tgtpts::Matrix{Float64}
    #===
    For an interior target i <= Ni:
        kt = cld(i, Np)
        ip = i - (kt - 1)*Np
        q, r = divrem(ip - 1, N)
        local coordinates are:
            t1 = cospi((2*(r+1)-1)/(2N))
            t2 = cospi((2*(q+1)-1)/(2N))
    
    For a boundary target i > Ni:
        ib  = i - Ni
        ibp = cld(ib, N)
        kt  = dom.kd[ibp]
        ip  = ib - (ibp - 1)*N
        local boundary coordinate is:
            t = cospi((2*ip - 1)/(2N))
    ===#

    # 2 × (Ni + Nb + Tnsp)
    # row 1 = global target column i
    # row 2 = actual integration patch k
    prepts::Matrix{Int}

    # 2 × Tnsp, inverse local coordinates
    # for near-singular volume interactions
    invpts::Matrix{Float64}

    #========== Boundary-DLP target classification ==========
    bdmode[i], i = 1:Ni:
    0 = target i is far from every boundary panel:
    dist(xᵢ, ∂Ωₖ) > del_near for all k
    
    1 = target i is moderately near the boundary:
    del_intp <= dist(xᵢ, ∂Ωₖ) <= del_near for at least one k,
    and dist(xᵢ, ∂Ωₖ) >= del_intp for all k
    
    2 = target i is very near the boundary:
    dist(xᵢ, ∂Ωₖ) < del_intp for at least one k
    
    Boundary DLP evaluation logic:
    bdmode == 0:
        direct full-boundary quadrature
    bdmode == 1:
        direct quadrature on far panels,
        smoothed quadrature on panels listed in bdneark/bdneart
    bdmode == 2:
        full-boundary normal interpolation using bdclosest/bdt/bddist
    ===#
    bdmode::Vector{UInt8}

    #====
    For targets with bdmode == 1, stored in encounter order.
    If jnear is the running index over mode 1 targets, then
      q1 = bdnearptr[jnear]
      q2 = bdnearptr[jnear + 1] - 1
    and for q = q1:q2
      bdneark[q] = actual boundary patch index
      bdneart[q] = projected local parameter on that boundary patch 
    ===#
    bdnearptr::Vector{Int}
    bdneark::Vector{Int}
    bdneart::Vector{Float64}

    # For targets with bdmode == 2, stored in encounter order.
    # If jintp is the running index over mode 2 targets, then
    #   bdclosest[jintp] = closest actual boundary patch index
    #   bdt[jintp]       = closest projected local boundary parameter
    #   bddist[jintp]    = closest distance to the boundary
    bdclosest::Vector{Int}
    bdt::Vector{Float64}
    bddist::Vector{Float64}

    #========== Mode-2 full-DLP interpolation data ========== 
    Lᵢₙ is the number of normal interpolation nodes used for very-near
    boundary targets. The first node is the boundary point itself: bdxvals[1] = 0
    and the remaining nodes are safe interior offset distances, for example
    
        bdxvals[j] = del_intp * j / 2,   j = 2:Lᵢₙ.
    
    For a mode-2 target with closest boundary projection γ(t₀), we evaluate
    the full DLP at auxiliary points x_j = γ(t₀) - bdxvals[j] * ν(t₀),
    then interpolate these full-DLP values back to the actual distance bddist[jintp].
    
    We only need near-panel classification for the auxiliary interior points
    j = 2:Lᵢₙ. The boundary node j = 1 is handled separately by the boundary
    limiting value plus the jump term.
    ==========#
    Lᵢₙ::Int
    bdxvals::Vector{Float64}

    #==========
    For every mode-2 target and every auxiliary interior interpolation point,
    store the boundary panels that are within del_near of that auxiliary point.

    The storage is flattened. Let Naux = Lᵢₙ - 1 and let jintp be the running 
    index over mode-2 targets. For interpolation node j = 2:Lᵢₙ, define
    
        a = (jintp - 1) * Naux + (j - 1).
    
    Then the near-panel data for that auxiliary point is stored in
    
        q1 = bdintpptr[a]
        q2 = bdintpptr[a + 1] - 1
    
    and for q = q1:q2:
    
        bdintpk[q] = actual boundary patch index near this auxiliary point
        bdintpt[q] = projected local parameter of the auxiliary point
                     onto that boundary patch.
    
    These projections are stored for each auxiliary point separately because
    on curved geometries, such as the kite, the projected parameter can change
    noticeably as the auxiliary point moves along the normal line.
    ==========#
    bdintpptr::Vector{Int}       # length NP2 * (Lᵢₙ - 1) + 1
    bdintpk::Vector{Int}         # flattened near panel indices
    bdintpt::Vector{Float64}     # corresponding projection parameters

    # pthgo[k] = first column in prepts for actual volume patch k
    # pthgo[M+1] = first column of the regular boundary-target block in prepts
    pthgo::Vector{Int}

    # nspsz[k] = number of interior targets near volume patch k
    # nspsz[M+1] = total number of boundary targets near volume patches
    nspsz::Vector{Int}

    function domprop(N::Integer, del::Float64, del_near::Float64, 
        del_intp::Float64, dom::D; Lᵢₙ = 5) where {D<:abstractdomain}

        # ---------------------------------------------------------------------
        # Notation used inside this constructor only
        # ---------------------------------------------------------------------
        # M       number of volume patches
        # Mbd     number of boundary patches
        # Np      number of tensor-product nodes per volume patch, N^2
        # Ni      number of interior targets, M*Np
        # Nb      number of boundary targets, Mbd*N
        # Nt      total number of targets, Ni+Nb
        #
        # k       actual patch index
        # kb      boundary-local patch index, 1:Mbd
        # kt      actual target patch index
        # ibp     boundary-local patch index containing a boundary target
        #
        # i       global target column in tgtpts
        # ib      boundary target index, 1:Nb
        # ip      local node index inside a patch
        #
        # colp    column in prepts
        # coli    column in invpts
        #
        # nsp_i   number of interior targets near each volume patch, length M
        # nsp_b   number of boundary targets near each volume patch, length M
        #
        # Xi      view of interior target coordinates, tgtpts[:, 1:Ni]
        # Xb      view of boundary target coordinates, tgtpts[:, Ni+1:Nt]
        # in_i    inside mask for Xi, length Ni
        # in_b    inside mask for Xb, length Nb
        # ---------------------------------------------------------------------
        M = dom.Npat
        Mbd = length(dom.kd)
        Np = N^2
        Ni = M * Np
        Nb = Mbd * N
        Nt = Ni + Nb

        snap1(x) = abs(x - 1.0) < 1e-14 ? 1.0 : abs(x + 1.0) < 1e-14 ? -1.0 : x

        # ---------------------------------------------------------------------
        # 1. Chebyshev nodes and target points
        # ---------------------------------------------------------------------
        tgtpts = Matrix{Float64}(undef, 2, Nt)

        z = Vector{Float64}(undef, N)
        @inbounds for ip in 1:N
            z[ip] = cospi((2*ip - 1) / (2N))
        end

        zx = repeat(z, 1, N)
        zy = repeat(z', N, 1)

        @views Xi = tgtpts[:, 1:Ni]
        @views Xb = tgtpts[:, Ni+1:Nt]

        X = similar(zx)
        Y = similar(zy)
        x = similar(z)
        y = similar(z)

        i = 1
        @inbounds for k in 1:M
            mapxy!(X, Y, dom, zx, zy, k)

            for ip in 1:Np
                Xi[1, i] = X[ip]
                Xi[2, i] = Y[ip]
                i += 1
            end
        end

        col = 1
        @inbounds for k in dom.kd
            gamx!(x, dom, z, k)
            gamy!(y, dom, z, k)

            for ip in 1:N
                Xb[1, col] = x[ip]
                Xb[2, col] = y[ip]
                col += 1
            end
        end

        # ---------------------------------------------------------------------
        # 2. Extended quadrilaterals and near-singular counts
        # ---------------------------------------------------------------------
        Q = dom.Qpts
        Qe = similar(Q)

        @inbounds for k in 1:M
            @views extendqua!(Qe[:, k], Q[:, k], del)
        end

        in_i = Vector{Bool}(undef, Ni)
        in_b = Vector{Bool}(undef, Nb)

        nsp_i = zeros(Int, M)
        nsp_b = zeros(Int, M)

        # Count interior and boundary targets near each volume patch.
        @inbounds for k in 1:M
            P1x, P1y = Qe[1, k], Qe[2, k]
            P2x, P2y = Qe[3, k], Qe[4, k]
            P3x, P3y = Qe[5, k], Qe[6, k]
            P4x, P4y = Qe[7, k], Qe[8, k]

            ptinqua!(in_i, Xi, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)
            ptinqua!(in_b, Xb, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)

            cnt = 0
            for i in 1:Ni
                kt = cld(i, Np)
                if in_i[i] && kt != k
                    if isnsp(dom, tgtpts[1, i], tgtpts[2, i], k, del)
                        cnt += 1
                    end
                end
            end
            nsp_i[k] = cnt

            cnt = 0
            for ib in 1:Nb
                i = Ni + ib
                ibp = cld(ib, N)
                kt = dom.kd[ibp]

                if in_b[ib] && kt != k
                    if isnsp(dom, tgtpts[1, i], tgtpts[2, i], k, del)
                        cnt += 1
                    end
                end
            end
            nsp_b[k] = cnt
        end

        Tnsp_i = sum(nsp_i)
        Tnsp_b = sum(nsp_b)
        Tnsp = Tnsp_i + Tnsp_b

        # ---------------------------------------------------------------------
        # 3. Layout arrays and output allocation
        # ---------------------------------------------------------------------
        pthgo = zeros(Int, M + 1)
        nspsz = zeros(Int, M + 1)

        colp = 1
        @inbounds for k in 1:M
            pthgo[k] = colp
            nspsz[k] = nsp_i[k]
            colp += Np + nsp_i[k]
        end

        pthgo[M+1] = colp
        nspsz[M+1] = Tnsp_b

        hmap = Dict{UInt64,Int}()
        sizehint!(hmap, Tnsp)

        prepts = Matrix{Int}(undef, 2, Ni + Nb + Tnsp)
        invpts = Matrix{Float64}(undef, 2, Tnsp)

        # ---------------------------------------------------------------------
        # 4. Inversion tables
        # ---------------------------------------------------------------------
        Nr = dom.pths[M].reg
        ftbs = Vector{FTable}(undef, Nr)

        @inbounds for r in 1:Nr
            ftbs[r] = inFTable(10_001)
            fill_FTable!(ftbs[r], dom, r)
        end

        # ---------------------------------------------------------------------
        # 5. Fill prepts and invpts for volume integrations
        # ---------------------------------------------------------------------
        coli = 1

        @inbounds for k in 1:M
            colp = pthgo[k]

            # Regular self-patch interior targets.
            for ip in 1:Np
                i = (k - 1) * Np + ip

                prepts[1, colp] = i
                prepts[2, colp] = k

                colp += 1
            end

            # Near-singular interior targets for volume patch k.
            P1x, P1y = Qe[1, k], Qe[2, k]
            P2x, P2y = Qe[3, k], Qe[4, k]
            P3x, P3y = Qe[5, k], Qe[6, k]
            P4x, P4y = Qe[7, k], Qe[8, k]

            ptinqua!(in_i, Xi, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)

            for i in 1:Ni
                kt = cld(i, Np)

                if in_i[i] && kt != k
                    if isnsp(dom, tgtpts[1, i], tgtpts[2, i], k, del)
                        prepts[1, colp] = i
                        prepts[2, colp] = k
                        hmap[packkey(i, k)] = colp

                        src_reg = dom.pths[k].reg

                        if src_reg == dom.pths[kt].reg
                            ip = i - (kt - 1) * Np
                            q, r = divrem(ip - 1, N)
                            Zx, Zy = mapinv2(dom, z[r+1], z[q+1], kt, k)
                        else
                            ftb = ftbs[src_reg]
                            Zx, Zy = mapinv(ftb, dom, tgtpts[1, i], tgtpts[2, i], k)
                        end

                        invpts[1, coli] = snap1(Zx)
                        invpts[2, coli] = snap1(Zy)

                        colp += 1
                        coli += 1
                    end
                end
            end
        end

        # Regular boundary targets.
        colp = pthgo[M+1]

        @inbounds for ib in 1:Nb
            i = Ni + ib
            ibp = cld(ib, N)
            kt = dom.kd[ibp]

            prepts[1, colp] = i
            prepts[2, colp] = kt

            colp += 1
        end

        # Near-singular boundary targets for volume patch k.
        @inbounds for k in 1:M
            P1x, P1y = Qe[1, k], Qe[2, k]
            P2x, P2y = Qe[3, k], Qe[4, k]
            P3x, P3y = Qe[5, k], Qe[6, k]
            P4x, P4y = Qe[7, k], Qe[8, k]

            ptinqua!(in_b, Xb, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)

            for ib in 1:Nb
                i = Ni + ib
                ibp = cld(ib, N)
                kt = dom.kd[ibp]

                if in_b[ib] && kt != k
                    if isnsp(dom, tgtpts[1, i], tgtpts[2, i], k, del)
                        prepts[1, colp] = i
                        prepts[2, colp] = k
                        hmap[packkey(i, k)] = colp

                        src_reg = dom.pths[k].reg

                        if src_reg == dom.pths[kt].reg
                            ip = ib - (ibp - 1) * N
                            Zx, Zy = mapinv2(dom, 1.0, z[ip], kt, k)
                        else
                            ftb = ftbs[src_reg]
                            Zx, Zy = mapinv(ftb, dom, tgtpts[1, i], tgtpts[2, i], k)
                        end

                        invpts[1, coli] = snap1(Zx)
                        invpts[2, coli] = snap1(Zy)

                        colp += 1
                        coli += 1
                    end
                end
            end
        end

        # ---------------------------------------------------------------------
        # 6. Boundary-DLP target classification
        # ---------------------------------------------------------------------
        # Output:
        #   bdmode[i], i = 1:Ni
        # Mode 1 storage:
        #   bdnearptr, bdneark, bdneart
        # Mode 2 storage:
        #   bdclosest, bdt, bddist
        # ---------------------------------------------------------------------
        BD_FAR = UInt8(0)
        BD_NEAR = UInt8(1)
        BD_INTP = UInt8(2)

        bdmode = fill(BD_FAR, Ni)

        # Boundary-patch extended quads for the two DLP thresholds.
        Qbd = dom.Qptsbd

        Qebd_near = similar(Qbd)
        Qebd_intp = similar(Qbd)

        @inbounds for kb in 1:Mbd
            @views extendqua!(Qebd_near[:, kb], Qbd[:, kb], del_near)
            @views extendqua!(Qebd_intp[:, kb], Qbd[:, kb], del_intp)
        end

        # Flattened storage for mode-1 targets.
        # bdnearptr starts with 1. Whenever we encounter a mode-1 target,
        # we append its nearby panels to bdneark/bdneart and then push the
        # next starting position into bdnearptr.
        bdnearptr = Int[1]
        bdneark = Int[]
        bdneart = Float64[]

        # Storage for mode-2 targets.
        bdclosest = Int[]
        bdt = Float64[]
        bddist = Float64[]

        # Small local work arrays for one target only.
        near_k = Int[]
        near_t = Float64[]

        @inbounds for i in 1:Ni

            x1 = tgtpts[1, i]
            x2 = tgtpts[2, i]

            # -------------------------------------------------------------
            # First look for very-near panels: dist < del_intp.
            # If this happens for any panel, this target is mode 2.
            # We only store the closest one.
            # -------------------------------------------------------------
            bestdist = Inf
            bestk = 0

            for kb in 1:Mbd
                k = dom.kd[kb]

                P1x, P1y = Qebd_intp[1, kb], Qebd_intp[2, kb]
                P2x, P2y = Qebd_intp[3, kb], Qebd_intp[4, kb]
                P3x, P3y = Qebd_intp[5, kb], Qebd_intp[6, kb]
                P4x, P4y = Qebd_intp[7, kb], Qebd_intp[8, kb]

                if ptinqua(x1, x2, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)

                    if isnspbd(dom, x1, x2, k, del_intp)

                        τ = projbd(dom, x1, x2, k)
                        γx, γy = gam(dom, τ, k)

                        dx = γx - x1
                        dy = γy - x2
                        dist = hypot(dx, dy)

                        if dist < del_intp && dist < bestdist
                            bestdist = dist
                            bestk = k
                        end
                    end
                end
            end

            if bestdist < del_intp

                # Very close to the boundary: full-boundary interpolation.
                bdmode[i] = BD_INTP

                # Refine closest projection to high precision.
                ε, t0 = GSS(dist_gam, -1.0, 1.0, 1e-15, dom, x1, x2, bestk)
                t₀ = Bis(dist_gamBis, -1.0, 1.0, t0, 64, dom, x1, x2, bestk)

                push!(bdclosest, bestk)
                push!(bdt, snap1(t₀))
                push!(bddist, ε)

                continue
            end

            # -------------------------------------------------------------
            # Otherwise, look for moderate-near panels:
            #     del_intp <= dist <= del_near.
            # These are the panels that get smoothed quadrature.
            # -------------------------------------------------------------
            empty!(near_k)
            empty!(near_t)

            for kb in 1:Mbd
                k = dom.kd[kb]

                P1x, P1y = Qebd_near[1, kb], Qebd_near[2, kb]
                P2x, P2y = Qebd_near[3, kb], Qebd_near[4, kb]
                P3x, P3y = Qebd_near[5, kb], Qebd_near[6, kb]
                P4x, P4y = Qebd_near[7, kb], Qebd_near[8, kb]

                if ptinqua(x1, x2, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)

                    if isnspbd(dom, x1, x2, k, del_near)

                        τ = projbd(dom, x1, x2, k)
                        γx, γy = gam(dom, τ, k)

                        dx = γx - x1
                        dy = γy - x2
                        dist = hypot(dx, dy)

                        if del_intp <= dist <= del_near
                            push!(near_k, k)
                            push!(near_t, snap1(τ))
                        end
                    end
                end
            end

            if !isempty(near_k)

                # Moderately close: panel-wise smoothing on listed panels.
                bdmode[i] = BD_NEAR

                for q in eachindex(near_k)
                    push!(bdneark, near_k[q])
                    push!(bdneart, near_t[q])
                end

                push!(bdnearptr, length(bdneark) + 1)

            else
                # Far from all boundary panels.
                bdmode[i] = BD_FAR
            end
        end

        # ---------------------------------------------------------------------
        # 7. Mode-2 auxiliary-point near-panel classification
        # ---------------------------------------------------------------------
        # For bdmode == 2 targets, the original target is closer than del_intp
        # to the boundary. We avoid evaluating the DLP there directly.
        # Instead, for each such target, we use normal interpolation. If the
        # closest boundary projection is γ_l(t₀), and ν_l(t₀) is the outward
        # unit normal, then the auxiliary interpolation points are
        #
        #     z_j = γ_l(t₀) - bdxvals[j] * ν_l(t₀),     j = 1:Lᵢₙ.
        #
        # The first interpolation node is the boundary point bdxvals[1] = 0.
        # For every mode-2 target and every auxiliary interior point, we store
        # the boundary panels lying within del_near, along with the projected
        # parameter of that auxiliary point onto each such panel.

        # The storage is flattened. Let Naux = Lᵢₙ - 1. For mode-2 target index 
        # jintp and interpolation node j = 2:Lᵢₙ, define
        #
        #     a = (jintp - 1) * Naux + (j - 1).
        #
        # Then the near panels for this auxiliary point are stored in
        #
        #     q1 = bdintpptr[a]
        #     q2 = bdintpptr[a + 1] - 1
        #
        # and for q = q1:q2:
        #
        #     bdintpk[q]  = actual boundary patch index
        #     bdintpt[q] = projected local parameter onto that patch.
        # ---------------------------------------------------------------------

        @assert Lᵢₙ >= 2 "Lᵢₙ must be at least 2."

        bdxvals = zeros(Float64, Lᵢₙ)

        @inbounds for j in 2:Lᵢₙ
            bdxvals[j] = del_intp * j / 2
        end

        bdintpptr = Int[1]
        bdintpk = Int[]
        bdintpt = Float64[]

        γbd = Vector{Float64}(undef, 2)
        ν = Vector{Float64}(undef, 2)

        @inbounds for jintp in eachindex(bdclosest)

            l = bdclosest[jintp]
            t0 = bdt[jintp]

            # Boundary projection and outward unit normal.
            gam!(γbd, dom, t0, l)
            nu!(ν, dom, t0, l)

            # Only classify auxiliary interior points.
            # The boundary point j = 1 is handled by the boundary limiting value.
            for m in 2:Lᵢₙ

                δ = bdxvals[m]

                zδ1 = γbd[1] - δ * ν[1]
                zδ2 = γbd[2] - δ * ν[2]

                # Find all boundary panels within del_near of this auxiliary point.
                for kb in 1:Mbd
                    k = dom.kd[kb]

                    P1x, P1y = Qebd_near[1, kb], Qebd_near[2, kb]
                    P2x, P2y = Qebd_near[3, kb], Qebd_near[4, kb]
                    P3x, P3y = Qebd_near[5, kb], Qebd_near[6, kb]
                    P4x, P4y = Qebd_near[7, kb], Qebd_near[8, kb]

                    if ptinqua(zδ1, zδ2, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)

                        if isnspbd(dom, zδ1, zδ2, k, del_near)

                            τ = projbd(dom, zδ1, zδ2, k)

                            γx, γy = gam(dom, τ, k)
                            dist = hypot(γx - zδ1, γy - zδ2)

                            if dist <= del_near
                                push!(bdintpk, k)
                                push!(bdintpt, snap1(τ))
                            end
                        end
                    end
                end

                # End of one auxiliary-point near-panel list.
                push!(bdintpptr, length(bdintpk) + 1)
            end
        end

        return new(N, del, del_near, del_intp, hmap, tgtpts, prepts, invpts,
            bdmode, bdnearptr, bdneark, bdneart, bdclosest, bdt, bddist,
            Lᵢₙ, bdxvals, bdintpptr, bdintpk, bdintpt, pthgo, nspsz)
    end
end

function plotdp(dp::domprop, d::abstractdomain; label=:none)

    fig, ax = draw(d)  # your existing domain-drawing routine

    if label != :none
        ax.backgroundcolor = RGBf(0.5, 0.5, 0.5) 
    end
    # Pull x/y once from the vector of targpt structs
    @views xs = dp.tgtpts[1,:]
    @views ys = dp.tgtpts[2,:]

    nint = d.Npat * dp.N^2         # interior count
    nall = size(dp.tgtpts, 2)      # interior + boundary

    ms =  label != :none ? 12 : 9

    # Interior targets (black)
    scatter!(ax, @view(xs[1:nint]), @view(ys[1:nint]);
        marker=:circle, markersize=ms,
        strokewidth=0, color=:black)

    # Boundary targets (blue)
    scatter!(ax, @view(xs[nint+1:nall]), @view(ys[nint+1:nall]);
        marker=:circle, markersize=ms,
        strokewidth=0, color=:black)

    if label != :none
        labelstep = 1
        if label === :int
            sel = 1:labelstep:nint
            lbl = cld.(1:nint, dp.N^2)
            lbj = collect(1:nint) .- (lbl .- 1).*dp.N^2
        elseif label === :bd
            sel = (nint+1):labelstep:nall
            lbl = cld.(1:(nall-nint), dp.N)
            lbj = collect(1:(nall-nint)) .- (lbl .- 1).*dp.N
        end

        xsel = xs[sel] 
        ysel = ys[sel] 

        labs = latexstring.(string.(lbj))

        for (x, y, s) in zip(xsel, ysel, labs)
            text!(ax, x, y; text = s, color = :white,
              align = (:center, :center), fontsize = 2*ms, font = :bold)
        end

    end

    return (fig, ax)
end

"""
Plot the near–singular target points to patch `k` on top of the domain drawing.
- Interior near–singulars to `k` are shown in black.
- Boundary near–singulars to `k` are shown in blue.
Also draws the (black) quadrilateral for patch `k` and its δ–extended quad.
"""
function plotns(dp::domprop, d::abstractdomain, k::Integer)

    fig, ax = draw(d)

    # -- draw quadrilateral for patch k --
    @views P = d.Qpts[:, k]      # 8-vector: [x1,y1,x2,y2,x3,y3,x4,y4]

    Pext = Vector{Float64}(undef, 8)

    ix = [1, 3, 5, 7, 1]
    iy = [2, 4, 6, 8, 2]
    #plot!(plt, P[ix], P[iy], color=:black, lw=2, label="")

    # -- draw δ-extended quadrilateral --
    extendqua!(Pext, P, dp.del)
    lines!(ax, Pext[ix], Pext[iy], color=:black, linewidth=2)

    # convenient sizes
    M   = d.Npat
    Np  = dp.N^2
    Mbd = length(d.kd)

    # ---------------- Interior near–singulars to k ----------------
    # In prepts layout: for each patch k, block = Np interior pts (k=l),
    # followed by dp.nspsz[k] near–singulars (k given patch index).
    S = dp.pthgo[k] + Np               # start of near–singulars for patch k
    E = dp.pthgo[k+1] - 1              # end index (inclusive) for that block

    xin = Vector{Float64}(undef, E - S + 1)
    yin = Vector{Float64}(undef, E - S + 1)
    @inbounds for (j, i) in enumerate(S:E)
        ti = dp.prepts[1, i]
        xin[j] = dp.tgtpts[1, ti]
        yin[j] = dp.tgtpts[2, ti]
    end
    scatter!(ax, xin, yin; markersize=9, color=:black, strokewidth=0)
 

    # ---------------- Boundary near–singulars to k ----------------
    # After dp.pthgo[M+1] we have:
    #   - a block of all boundary targets (length Nb)
    #   - then all boundary near–singulars (various k).
    start_bd_block = dp.pthgo[M+1] + Mbd * dp.N

    xb = Float64[]
    yb = Float64[]

    @inbounds for i in start_bd_block:size(dp.prepts,2)
        if dp.prepts[2, i] == k
            ti = dp.prepts[1, i]
            push!(xb, dp.tgtpts[1, ti])
            push!(yb, dp.tgtpts[2, ti])
        end
    end

    !isempty(xb) && scatter!(ax, xb, yb;  markersize=9, color=:black, strokewidth=0)

    return (fig,ax)
end


"""
Plot boundary-DLP projections associated with boundary patch `k`.
Here `k` is the actual boundary patch index, i.e. `k ∈ d.kd`.
Mode-1 targets:
1) target point shown in black,
2) projected point shown in black,
3) arrow shown in black.
For Mode-2 targets same but in blue.

Also draws the del_near extended boundary quadrilateral in black,
and the del_intp extended boundary quadrilateral in blue.
"""
function plotprojbd(dp::domprop, d::abstractdomain, k::Integer)

    fig, ax = draw(d, 1; L=1024)

    # k is the actual boundary patch index.
    kb = findfirst(==(k), d.kd)
    kb === nothing && error("Patch k = $k is not a boundary patch.")

    # -------------------------------------------------------------------------
    # Draw extended boundary quadrilaterals.
    # -------------------------------------------------------------------------
    @views Pbd = d.Qptsbd[:, kb]

    Pext_near = Vector{Float64}(undef, 8)
    Pext_intp = Vector{Float64}(undef, 8)

    ix = [1, 3, 5, 7, 1]
    iy = [2, 4, 6, 8, 2]

    extendqua!(Pext_near, Pbd, dp.del_near)
    extendqua!(Pext_intp, Pbd, dp.del_intp)

    lines!(ax, Pext_near[ix], Pext_near[iy];
        color=:black, linewidth=2)

    lines!(ax, Pext_intp[ix], Pext_intp[iy];
        color=:blue, linewidth=2)

    # -------------------------------------------------------------------------
    # Collect mode-1 projections for patch k.
    # -------------------------------------------------------------------------
    xt1 = Float64[]
    yt1 = Float64[]
    xp1 = Float64[]
    yp1 = Float64[]

    xt2 = Float64[]
    yt2 = Float64[]
    xp2 = Float64[]
    yp2 = Float64[]

    jnear = 0
    jintp = 0

    Ni = length(dp.bdmode)

    @inbounds for i in 1:Ni

        mode = dp.bdmode[i]

        if mode == UInt8(1)
            jnear += 1

            q1 = dp.bdnearptr[jnear]
            q2 = dp.bdnearptr[jnear+1] - 1

            @inbounds for q in q1:q2
                if dp.bdneark[q] == k
                    τ = dp.bdneart[q]

                    x0 = dp.tgtpts[1, i]
                    y0 = dp.tgtpts[2, i]

                    x1 = gamx(d, τ, k)
                    y1 = gamy(d, τ, k)

                    push!(xt1, x0)
                    push!(yt1, y0)
                    push!(xp1, x1)
                    push!(yp1, y1)
                end
            end

        elseif mode == UInt8(2)
            jintp += 1

            if dp.bdclosest[jintp] == k
                τ = dp.bdt[jintp]

                x0 = dp.tgtpts[1, i]
                y0 = dp.tgtpts[2, i]

                x1 = gamx(d, τ, k)
                y1 = gamy(d, τ, k)

                push!(xt2, x0)
                push!(yt2, y0)
                push!(xp2, x1)
                push!(yp2, y1)
            end
        end
    end

    # -------------------------------------------------------------------------
    # Plot mode-1 arrows/projections.
    # -------------------------------------------------------------------------
    if !isempty(xt1)
        dmin1 = minimum(hypot.(xp1 .- xt1, yp1 .- yt1))
        @printf("Smallest mode-1 distance: %.2e\n", dmin1)

        arrows2d!(ax, xt1, yt1, xp1 .- xt1, yp1 .- yt1;
            tipwidth=8,
            tiplength=12,
            shaftwidth=1.5,
            color=:black)

        scatter!(ax, xt1, yt1;
            markersize=9,
            color=:black,
            strokewidth=0)

        scatter!(ax, xp1, yp1;
            markersize=9,
            color=:black,
            strokewidth=0)
    end

    # -------------------------------------------------------------------------
    # Plot mode-2 arrows/projections.
    # -------------------------------------------------------------------------
    if !isempty(xt2)
        dmin2 = minimum(hypot.(xp2 .- xt2, yp2 .- yt2))
        @printf("Smallest mode-2 distance: %.2e\n", dmin2)

        arrows2d!(ax, xt2, yt2, xp2 .- xt2, yp2 .- yt2;
            tipwidth=8,
            tiplength=12,
            shaftwidth=1.5,
            color=:blue)

        scatter!(ax, xt2, yt2;
            markersize=9,
            color=:blue,
            strokewidth=0)

        scatter!(ax, xp2, yp2;
            markersize=9,
            color=:blue,
            strokewidth=0)
    end

    return fig, ax
end

"""
Plot boundary-near interior target points associated with boundary patch `k`.
Here `k` is the actual boundary patch index, i.e. `k ∈ d.kd`.

i) Mode-1 targets are shown in black:
  del_intp <= dist(x, ∂Ω_k) <= del_near for this patch `k`.
ii) Mode-2 targets are shown in blue:
  dist(x, ∂Ω) < del_intp and the closest boundary patch is `k`.

Also draws the del_near-extended boundary quadrilateral in black,
and the del_intp-extended boundary quadrilateral in blue.
"""
function plotnsbd(dp::domprop, d::abstractdomain, k::Integer)

    fig, ax = draw(d, 1)

    # k is the actual boundary patch index.
    # kb is the boundary-local index such that d.kd[kb] == k.
    kb = findfirst(==(k), d.kd)
    kb === nothing && error("Patch k = $k is not a boundary patch. It must belong to d.kd.")

    # -------------------------------------------------------------------------
    # Draw extended boundary quadrilaterals.
    # -------------------------------------------------------------------------
    @views Pbd = d.Qptsbd[:, kb]

    Pext_near = Vector{Float64}(undef, 8)
    Pext_intp = Vector{Float64}(undef, 8)

    ix = [1, 3, 5, 7, 1]
    iy = [2, 4, 6, 8, 2]

    extendqua!(Pext_near, Pbd, dp.del_near)
    extendqua!(Pext_intp, Pbd, dp.del_intp)

    lines!(ax, Pext_near[ix], Pext_near[iy];
        color = :black, linewidth = 2)

    lines!(ax, Pext_intp[ix], Pext_intp[iy];
        color = :blue, linewidth = 2)

    # -------------------------------------------------------------------------
    # Collect mode-1 and mode-2 target points associated with patch k.
    # -------------------------------------------------------------------------
    x_mode1 = Float64[]
    y_mode1 = Float64[]

    x_mode2 = Float64[]
    y_mode2 = Float64[]

    jnear = 0
    jintp = 0

    Ni = length(dp.bdmode)

    @inbounds for i in 1:Ni

        mode = dp.bdmode[i]

        if mode == UInt8(1)
            # Mode 1: moderately near.
            # This target may be near several boundary panels.
            jnear += 1

            q1 = dp.bdnearptr[jnear]
            q2 = dp.bdnearptr[jnear + 1] - 1

            for q in q1:q2
                if dp.bdneark[q] == k
                    push!(x_mode1, dp.tgtpts[1, i])
                    push!(y_mode1, dp.tgtpts[2, i])
                    break
                end
            end

        elseif mode == UInt8(2)
            # Mode 2: very near.
            # Stored only by closest boundary patch.
            jintp += 1

            if dp.bdclosest[jintp] == k
                push!(x_mode2, dp.tgtpts[1, i])
                push!(y_mode2, dp.tgtpts[2, i])
            end
        end
    end

    # -------------------------------------------------------------------------
    # Plot target points.
    # -------------------------------------------------------------------------
    if !isempty(x_mode1)
        scatter!(ax, x_mode1, y_mode1;
            markersize = 9,
            color = :black,
            strokewidth = 0)
    end

    if !isempty(x_mode2)
        scatter!(ax, x_mode2, y_mode2;
            markersize = 9,
            color = :blue,
            strokewidth = 0)
    end

    return fig, ax
end

"""
Plot the mode-2 auxiliary interpolation geometry for a given jintp.

For the mode-2 target indexed by `jintp`, this plots:
- the original target point in red,
- the closest boundary projection in black,
- auxiliary interpolation points in blue,
- projections of each auxiliary point onto its nearby boundary panels in green,
- arrows from each auxiliary point to those projections,
- del_near-extended boundary quadrilaterals for the involved panels.
"""
function plotbdintp(dp::domprop, d::abstractdomain, jintp::Integer)

    NP2 = length(dp.bdclosest)
    1 <= jintp <= NP2 || error("jintp = $jintp is outside 1:$NP2")

    fig, ax = draw(d, 1; show=false)

    l  = dp.bdclosest[jintp]
    t0 = dp.bdt[jintp]
    ε  = dp.bddist[jintp]

    γbd = Vector{Float64}(undef, 2)
    ν   = Vector{Float64}(undef, 2)

    gam!(γbd, d, t0, l)
    nu!(ν, d, t0, l)

    # Original very-near target.
    xt = γbd[1] - ε * ν[1]
    yt = γbd[2] - ε * ν[2]

    # Plot original target and boundary projection.
    scatter!(ax, [xt], [yt];
        markersize=12,
        color=:red,
        strokewidth=0)

    scatter!(ax, [γbd[1]], [γbd[2]];
        markersize=12,
        color=:black,
        strokewidth=0)

    arrows2d!(ax, [xt], [yt], [γbd[1] - xt], [γbd[2] - yt];
        tipwidth=8,
        tiplength=12,
        shaftwidth=1.5,
        color=:red)

    # Draw del_near boxes only for panels appearing in bdintpk for this jintp.
    involved = Set{Int}()

    Naux = dp.Lᵢₙ - 1

    for m in 2:dp.Lᵢₙ
        a = (jintp - 1) * Naux + (m - 1)

        q1 = dp.bdintpptr[a]
        q2 = dp.bdintpptr[a + 1] - 1

        for q in q1:q2
            push!(involved, dp.bdintpk[q])
        end
    end

    ix = [1, 3, 5, 7, 1]
    iy = [2, 4, 6, 8, 2]
    Pext = Vector{Float64}(undef, 8)

    for k in involved
        kb = findfirst(==(k), d.kd)
        kb === nothing && continue

        @views Pbd = d.Qptsbd[:, kb]
        extendqua!(Pext, Pbd, dp.del_near)

        lines!(ax, Pext[ix], Pext[iy];
            color=:gray,
            linewidth=2)
    end

    # Plot auxiliary points and their projections.
    xaux = Float64[]
    yaux = Float64[]

    xproj = Float64[]
    yproj = Float64[]

    xarr = Float64[]
    yarr = Float64[]
    uarr = Float64[]
    varr = Float64[]

    @inbounds for m in 2:dp.Lᵢₙ

        δ = dp.bdxvals[m]

        zδ1 = γbd[1] - δ * ν[1]
        zδ2 = γbd[2] - δ * ν[2]

        push!(xaux, zδ1)
        push!(yaux, zδ2)

        a = (jintp - 1) * Naux + (m - 1)

        q1 = dp.bdintpptr[a]
        q2 = dp.bdintpptr[a + 1] - 1

        for q in q1:q2
            k = dp.bdintpk[q]
            τ = dp.bdintpt[q]

            xp = gamx(d, τ, k)
            yp = gamy(d, τ, k)

            push!(xproj, xp)
            push!(yproj, yp)

            push!(xarr, zδ1)
            push!(yarr, zδ2)
            push!(uarr, xp - zδ1)
            push!(varr, yp - zδ2)
        end
    end

    if !isempty(xaux)
        scatter!(ax, xaux, yaux;
            markersize=10,
            color=:blue,
            strokewidth=0)
    end

    if !isempty(xproj)
        arrows2d!(ax, xarr, yarr, uarr, varr;
            tipwidth=7,
            tiplength=10,
            shaftwidth=1.2,
            color=:green)

        scatter!(ax, xproj, yproj;
            markersize=8,
            color=:green,
            strokewidth=0)
    end

    @printf("jintp = %d\n", jintp)
    @printf("closest patch l = %d\n", l)
    @printf("t0 = %.16e\n", t0)
    @printf("epsilon = %.16e\n", ε)
    @printf("number of involved panels = %d\n", length(involved))
    @printf("involved panels = %s\n", collect(involved))

    display(GLMakie.Screen(), fig)

    return fig, ax
end

"""
Plot a closed-form function `f(x,y)` over the patched domain.

The function is evaluated on a dense closed Chebyshev grid on each patch,
mapped to physical space with `mapxy`. All patch surfaces share one reactive
color scale, so the colormap is consistent across the whole domain.
"""
function plotfunc(dp::domprop, d::abstractdomain, f!::Function)

    # -------------------------------------------------------------------------
    # Notation
    # -------------------------------------------------------------------------
    # n     number of Chebyshev nodes per coordinate direction used by dp
    # ℓ     number of plotting nodes per coordinate direction, 4n+1
    # M     number of patches
    # k     actual patch index
    # -------------------------------------------------------------------------
    n = dp.N
    ℓ = 4n + 1
    M = d.Npat

    z  = cospi.((0:4n) ./ (4n))
    zx = repeat(z, 1, ℓ)
    zy = repeat(z', ℓ, 1)
    X = similar(zx)
    Y = similar(zx)
    F = similar(zx)

    cmap = :jet1

    fig = Figure(size = (800, 650))
    ax  = Axis3(fig[1, 1];
        aspect = :equal,
        xlabel = latexstring("\$x_1\$"),
        ylabel = latexstring("\$x_2\$"),
        zlabel = latexstring("\$f(x_1,x_2)\$"),
        xlabelsize = 22,
        ylabelsize = 22,
        zlabelsize = 22,
        xticklabelsize = 18,
        yticklabelsize = 18,
        zticklabelsize = 18)

    widen_clims!(clims, Z) = begin
        zlo, zhi = extrema(Z)
        clo, chi = clims[]
        clims[] = (min(clo, zlo), max(chi, zhi))
        nothing
    end

    clims = Observable((Inf, -Inf))
    surf = nothing

    @inbounds for k in 1:M
        # compute into reusable buffers
        mapxy!(X, Y, d, zx, zy, k)
        f!(F, X, Y)

        # freeze this patch for plotting
        Zk = real.(F)
        Xk = copy(X)
        Yk = copy(Y)

        if k == 1
            clims[] = extrema(Zk)
        else
            widen_clims!(clims, Zk)
        end

        p = surface!(ax, Xk, Yk, Zk;
            color = Zk,
            colorrange = clims,
            colormap = cmap,
            shading = NoShading)

        if surf === nothing
            surf = p
        end
    end

    Colorbar(fig[1, 2], surf;
        label = latexstring("\$f(x_1,x_2)\$"),
        labelsize = 20,
        ticklabelsize = 16)

    display(GLMakie.Screen(), fig)
    return fig, ax
end

function plotfunc(dp::domprop, d::abstractdomain, f::AbstractVector; interpolate::Bool = false)
    # -------------------------------------------------------------------------
    # Notation
    # -------------------------------------------------------------------------
    # n     number of Chebyshev nodes per coordinate direction
    # M     number of patches
    # Np    number of tensor-product nodes per patch, n^2
    # k     actual patch index
    # ℓ     number of plotting nodes per coordinate direction in interpolated mode
    # -------------------------------------------------------------------------

    n  = dp.N
    M  = d.Npat
    Np = n^2

    cmap = :jet1

    fig = Figure(size = (800, 650))
    ax  = Axis3(fig[1, 1];
        aspect = :equal,
        xlabel = latexstring("\$x_1\$"),
        ylabel = latexstring("\$x_2\$"),
        zlabel = latexstring("\$f(x_1,x_2)\$"),
        xlabelsize = 22,
        ylabelsize = 22,
        zlabelsize = 22,
        xticklabelsize = 18,
        yticklabelsize = 18,
        zticklabelsize = 18)

    widen_clims!(clims, Z) = begin
        zlo, zhi = extrema(real.(Z))
        clo, chi = clims[]
        clims[] = (min(clo, zlo), max(chi, zhi))
        nothing
    end

    clims = Observable((Inf, -Inf))
    surf = nothing

    if !interpolate
        # ---------------------------------------------------------------------
        # Plot raw values on the open Chebyshev grid
        # ---------------------------------------------------------------------
        z  = cospi.((2 .* (1:n) .- 1) ./ (2n))
        zx = repeat(z, 1, n)
        zy = repeat(z', n, 1)

        @inbounds for k in 1:M
            idx = (k - 1) * Np + 1 : k * Np

            X, Y = mapxy(d, zx, zy, k)
            Z = real.(reshape(@view(f[idx]), n, n))

            if k == 1
                clims[] = extrema(Z)
            else
                widen_clims!(clims, Z)
            end

            p = surface!(ax, X, Y, Z;
                color = Z,
                colorrange = clims,
                colormap = cmap,
                shading = NoShading)

            if surf === nothing
                surf = p
            end

            scatter!(ax, X, Y, Z;
                marker = :circle,
                markersize = 9,
                color = :black)
        end

    else
        # ---------------------------------------------------------------------
        # Interpolate onto a denser closed Chebyshev grid, then plot
        # ---------------------------------------------------------------------
        ℓ  = 4n + 1
        z  = cospi.((0:4n) ./ (4n))
        zx = repeat(z, 1, ℓ)
        zy = repeat(z', ℓ, 1)

        T = ChebyTN(n, z)

        chebcoef = similar(f, M * Np)

        @inbounds for k in 1:M
            idx = (k - 1) * Np + 1 : k * Np
            F = reshape(@view(f[idx]), n, n)

            # DCT-II in second dimension
            C = (1 / n) .* FFTW.r2r(F, FFTW.REDFT10, 2)
            C[:, 1] ./= 2

            # DCT-II in first dimension
            C = (1 / n) .* FFTW.r2r(C, FFTW.REDFT10, 1)
            C[1, :] ./= 2

            chebcoef[idx] .= vec(C)
        end

        @inbounds for k in 1:M
            idx = (k - 1) * Np + 1 : k * Np

            C = reshape(@view(chebcoef[idx]), n, n)
            Z = real.(T * C * transpose(T))
            X, Y = mapxy(d, zx, zy, k)

            if k == 1
                clims[] = extrema(Z)
            else
                widen_clims!(clims, Z)
            end

            p = surface!(ax, X, Y, Z;
                color = Z,
                colorrange = clims,
                colormap = cmap,
                shading = NoShading)

            if surf === nothing
                surf = p
            end
        end
    end

    Colorbar(fig[1, 2], surf;
        label = latexstring("\$f(x_1,x_2)\$"),
        labelsize = 20,
        ticklabelsize = 16)

    display(GLMakie.Screen(), fig)
    return fig, ax
end

# ---------------------------------------------------------------
# u given in closed form; enforce u=0 on the boundary of domain.
# Closed Chebyshev grid includes the boundary, so we can zero edges.
# ---------------------------------------------------------------
function plotu(dp::domprop, d::abstractdomain, u::Function)

    # -------------------------------------------------------------------------
    # Notation
    # -------------------------------------------------------------------------
    # n     number of Chebyshev nodes per coordinate direction used by dp
    # ℓ     number of plotting nodes per coordinate direction, 4n+1
    # M     number of volume patches
    # k     actual patch index
    # -------------------------------------------------------------------------

    n = dp.N
    M = d.Npat
    ℓ = 4n + 1

    # Closed Chebyshev nodes on [-1,1]
    z  = cospi.((0:4n) ./ (4n))
    zx = repeat(z, 1, ℓ)
    zy = repeat(z', ℓ, 1)

    cmap = :jet1

    fig = Figure(size = (650, 650))
    ax  = Axis3(fig[1, 1];
        aspect = :equal,
        xlabel = latexstring("\$x_1\$"),
        ylabel = latexstring("\$x_2\$"),
        zlabel = latexstring("\$u(x_1,x_2)\$"),
        xlabelsize = 22,
        ylabelsize = 22,
        zlabelsize = 22,
        xticklabelsize = 18,
        yticklabelsize = 18,
        zticklabelsize = 18)

    widen_clims!(clims, Z) = begin
        zlo, zhi = extrema(Z)
        clo, chi = clims[]
        clims[] = (min(clo, zlo), max(chi, zhi))
        nothing
    end

    clims = Observable((Inf, -Inf))

    surf = nothing

    @inbounds for k in 1:M

        X, Y = mapxy(d, zx, zy, k)

        if isbdpatch(d, k)
            # Boundary patch: first row corresponds to the domain boundary,
            # so enforce homogeneous Dirichlet data there.
            Z = u.(X[2:ℓ, :], Y[2:ℓ, :])
            Z = vcat(zeros(Float64, ℓ)', Z)
        else
            Z = u.(X, Y)
        end

        if k == 1
            clims[] = extrema(Z)
        else
            widen_clims!(clims, Z)
        end

        p = surface!(ax, X, Y, Z;
            color=Z,
            colorrange=clims,
            colormap=cmap,
            shading=NoShading)

        surf === nothing && (surf = p)
    end

    Colorbar(fig[1, 2], surf;
        label = latexstring("\$u(x_1,x_2)\$"),
        labelsize = 20,
        ticklabelsize = 16)

    display(GLMakie.Screen(), fig)
    return fig, ax
end

function plotu(dp::domprop, d::abstractdomain, u::AbstractVector, s::Real; interpolate::Bool = true)

    # -------------------------------------------------------------------------
    # Notation
    # -------------------------------------------------------------------------
    # n     number of Chebyshev nodes per coordinate direction
    # M     number of volume patches
    # Np    number of tensor-product nodes per patch, n^2
    # i     global interior target index
    # ip    local node index inside patch
    # k     actual patch index
    # j1    first local tensor-product index
    # ℓ     number of plotting nodes per coordinate direction
    # v     plotted/scaled copy of u
    # -------------------------------------------------------------------------

    n = dp.N
    M = d.Npat
    Np = n^2

    cmap = :jet1

    # Work on a copy so the input vector u is not mutated.
    v = copy(u)

    fig = Figure(size = (850, 650))

    ax = Axis3(fig[1, 1];
        xlabel = latexstring("\$x_1\$"),
        ylabel = latexstring("\$x_2\$"),
        zlabel = latexstring("\$u(x_1,x_2)\$"),
        xlabelsize = 22,
        ylabelsize = 22,
        zlabelsize = 22,
        xticklabelsize = 18,
        yticklabelsize = 18,
        zticklabelsize = 18)

    # -------------------------------------------------------------------------
    # Small local helper: update color limits
    # -------------------------------------------------------------------------
    widen_clims!(clims, Z) = begin
        zlo, zhi = extrema(real.(Z))
        clo, chi = clims[]
        clims[] = (min(clo, zlo), max(chi, zhi))
        nothing
    end

    # -------------------------------------------------------------------------
    # Case 1: no interpolation, plot raw samples
    # -------------------------------------------------------------------------
    if !interpolate

        z = cospi.((2 .* (1:n) .- 1) ./ (2n))

        zx = repeat(z, 1, n)
        zy = repeat(z', n, 1)

        # For boundary patches, add the boundary row ξ = 1.
        zxb = vcat(ones(Float64, n)', zx)
        zyb = vcat(z', zy)

        # Apply distance factor to interior samples.
        @inbounds for i in 1:M*Np
            k = cld(i, Np)
            ip = i - (k - 1) * Np

            j1 = mod(ip - 1, n) + 1
            ρ = 2 * sinpi((2*j1 - 1) / (4n))^2

            v[i] *= dfunc(d, k, ρ, s)
        end

        clims = Observable((Inf, -Inf))

        @inbounds for k in 1:M
            idx = (k - 1) * Np + 1 : k * Np
            Z = reshape(@view(v[idx]), n, n)

            if isbdpatch(d, k)
                X, Y = mapxy(d, zxb, zyb, k)
                Z = vcat(zeros(Float64, n)', Z)
            else
                X, Y = mapxy(d, zx, zy, k)
            end

            if k == 1
                clims[] = extrema(real.(Z))
            else
                widen_clims!(clims, Z)
            end

            surface!(ax, X, Y, real.(Z);
                color = real.(Z),
                colorrange = clims,
                colormap = cmap,
                shading = NoShading)

            scatter!(ax, X, Y, real.(Z);
                marker = :circle,
                markersize = 9,
                color = :black)
        end

        display(GLMakie.Screen(), fig)
        return fig, ax
    end

    # -------------------------------------------------------------------------
    # Case 2: interpolation, plot denser Chebyshev samples
    # -------------------------------------------------------------------------

    ℓ = 4n + 1

    z = cospi.((0:4n) ./ (4n))

    zx = repeat(z, 1, ℓ)
    zy = repeat(z', ℓ, 1)

    T = ChebyTN(n, z)

    chebcoef = similar(v, M * Np)

    # Compute Chebyshev coefficients patch-by-patch.
    @inbounds for k in 1:M
        idx = (k - 1) * Np + 1 : k * Np
        F = reshape(@view(v[idx]), n, n)

        # DCT-II in second dimension.
        C = (1 / n) .* FFTW.r2r(F, FFTW.REDFT10, 2)
        C[:, 1] ./= 2

        # DCT-II in first dimension.
        C = (1 / n) .* FFTW.r2r(C, FFTW.REDFT10, 1)
        C[1, :] ./= 2

        chebcoef[idx] .= vec(C)
    end

    clims = Observable((Inf, -Inf))

    @inbounds for k in 1:M
        idx = (k - 1) * Np + 1 : k * Np

        C = reshape(@view(chebcoef[idx]), n, n)
        Z = T * C * transpose(T)

        X, Y = mapxy(d, zx, zy, k)

        # Apply distance factor row-by-row.
        for row in 2:ℓ
            ρ = 2 * sinpi((row - 1) / (8n))^2
            Z[row, :] .*= dfunc(d, k, ρ, s)
        end

        if isbdpatch(d, k)
            Z[1, :] .= 0
        else
            Z[1, :] .*= dfunc(d, k, 0.0, s)
        end

        if k == 1
            clims[] = extrema(real.(Z))
        else
            widen_clims!(clims, Z)
        end

        surface!(ax, X, Y, real.(Z);
            color = real.(Z),
            colorrange = clims,
            colormap = cmap,
            shading = NoShading)
    end

    display(GLMakie.Screen(), fig)
    return fig, ax
end

# outdir: folder where files will be written (default "paraview_out")
# basename: prefix for filenames (default "u_patch")
#
# Writes:
#   outdir/basename.vtm
# plus the per-block VTK files referenced by the .vtm container.
#
function u_to_paraview(dp::domprop, d::abstractdomain, u::AbstractVector, s::Real;
    outdir::AbstractString = "paraview_out", basename::AbstractString = "u_patch")

    n  = dp.N
    M  = d.Npat
    Np = n^2

    ℓ = 4n + 1
    z = cospi.((0:4n) ./ (4n))

    zx = repeat(z, 1, ℓ)
    zy = repeat(z', ℓ, 1)

    T = ChebyTN(n, z)

    chebcoef = similar(u, M * Np)
    for k in 1:M
        idx = (k - 1) * Np + 1 : k * Np
        fv = reshape(@views(u[idx]), n, n)

        c = (1 / n) .* FFTW.r2r(fv, FFTW.REDFT10, 2)
        c[:, 1] ./= 2

        c = (1 / n) .* FFTW.r2r(c, FFTW.REDFT10, 1)
        c[1, :] ./= 2

        chebcoef[idx] .= vec(c)
    end

    isdir(outdir) || mkpath(outdir)

    vtm = vtk_multiblock(joinpath(outdir, basename))

    for P in 1:M
        c = reshape(@view(chebcoef[(P - 1) * Np + 1 : P * Np]), n, n)
        Zc = T * c * transpose(T)
        X, Y = mapxy(d, zx, zy, P)

        for i in 2:ℓ
            Zc[i, :] .*= dfunc(d, P, 2 * sinpi((i - 1) / (8 * n))^2, s)
        end

        if isbdpatch(d, P)
            Zc[1, :] .= 0.0
        else
            Zc[1, :] .*= dfunc(d, P, 0.0, s)
        end

        nx, ny = size(X)
        X3 = reshape(X, nx, ny, 1)
        Y3 = reshape(Y, nx, ny, 1)
        Z3 = reshape(Zc, nx, ny, 1)

        vtk = vtk_grid(vtm, X3, Y3, Z3)
        vtk["u"] = Z3
        vtk_save(vtk)
    end

    vtk_save(vtm)
    return nothing
end

"""
    plotu(dp::domprop, d::abstractdomain; 
    file::AbstractString="uex_vals.bin", interpolate::Bool = false)

Read a raw binary file of real numbers (`Float64`) into a vector and plot it
using `plotu(dp, d, u::AbstractVector; interpolate=...)`.

- The file is assumed to be **raw**, with no headers, just function values in a vector form.
"""
function plotu(dp::domprop, d::abstractdomain, s::Real;
    file::AbstractString="uex_vals.bin", interpolate::Bool=true)

    # change into the folder only for the duration of the block
    u = cd(raw"C:\Users\sabhr\OneDrive\Desktop\Caltech\Research\Code Julia\2D FL") do
        bytes = stat(file).size
        sz = sizeof(Float64)
        bytes % sz == 0 || error("File size $bytes is not a multiple of sizeof(Float64) = $sz.")
        n = Int(bytes ÷ sz)

        open(file, "r") do io
            v = Vector{Float64}(undef, n)
            read!(io, v)              
            v
        end
    end

    return plotu(dp, d, u, s; interpolate=interpolate)
end


"""
    chkinvpts(dp::domprop, d::abstractdomain) -> Vector{Float64}

Check that inverse points stored in `dp.invpts` are correct.

For every near–singular precomputation entry `pt = dp.prepts[i]`, 
that is, those with `pt.k != pt.l`, we:
-- locate its running index `ll` inside the near–singular list (the columns
   of `dp.invpts`), and
-- map the stored inverse `(t,s) = dp.invpts[:,ll]` back to real space via
   `mapxy(d, t, s, pt.k)`, then
-- record the Euclidean error to the original target `(pt.x, pt.y)`.

Returns a vector `err` of length `Tnsp = sum(dp.nspsz)`, ordered exactly like
`dp.invpts` (first the interior near–singulars, then the boundary ones).
"""
function chkinvpts(dp::domprop, d::abstractdomain,flag=nothing)::Float64
    #bookkeeping
    N   = dp.N
    M   = d.Npat
    Np  = N^2
    Mbd = length(d.kd)
    nbd = Mbd * N
    
    # Col index where the bd near–singular block starts in `prepts`
    bdnsp = dp.pthgo[M+1] + nbd        

    Tnsp  = sum(dp.nspsz)

    err  = Vector{Float64}(undef, Tnsp)

    Lₚ = Np * M + nbd + Tnsp # Total cols in Prepts matrix

    z = Vector{Float64}(undef, N)

    @inbounds for j in 1:N
        z[j] = cos(π*(2j-1)/(2N))
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

    @inbounds for i in 1:Lₚ
        if NSI[i] == true
            #Target point with linear index ti near singuluar to patch k
            ti = dp.prepts[1, i]

            k = dp.prepts[2, i]

            # (prepts index) -> `ll` (column in invpts)
            if i < bdnsp
                # interior near–singular to patch k
                ll = i - k * Np
                # Building the target point
                ℓ = ceil(Int, ti / Np)

                jj = ti - (ℓ - 1) * Np

                q, r = divrem(jj - 1, N)

                x, y = z[r+1], z[q+1]

                tx, ty = mapxy(d, x, y, ℓ)
            else
                # boundary near–singulars
                ll = i - M * Np - nbd
                # Building the target point
                k₀ = ceil(Int, (ti - M * Np) / N)

                jj = ti - M * Np - (k₀ - 1) * N

                ℓ = d.kd[k₀]

                tx, ty = mapxy(d, 1.0, z[jj], ℓ)
            end

            # pull (t,s) and map back to real space on patch k
            t̂ = dp.invpts[1, ll]
            ŝ = dp.invpts[2, ll]

            zx, zy = mapxy(d, t̂, ŝ, k)

            E = hypot(tx - zx, ty - zy)
            # Euclidean error to the original target
            err[ll] = E < 1e-16 ? 1e-16 : E
        end
    end

    flag = isnothing(flag) ? 0 : 1

    if flag == 1
        fig = Figure(size=(650, 650))

        ax = Axis(fig[1, 1];
            xlabel=latexstring("\$x\$"),
            ylabel="Error",
            xlabelsize=22,
            ylabelsize=22,
            xticklabelsize=18,
            yticklabelsize=18,
            yscale=log10,                # log y-axis
            xminorgridvisible=true,
            yminorgridvisible=true,
        )

        xs = collect(1:length(err))
        lines!(ax, xs, err; linewidth=1.2)
        scatter!(ax, xs, err;
            marker=:circle,
            markersize=9,
            color=:transparent,      # hollow fill
            strokecolor=:black,      # outline color
            strokewidth=1.2)

        ax.title = "Inverse-map error"

        display(GLMakie.Screen(), fig)
    end

    return maximum(err)
end

function _preview(v::Vector{Int}; n::Int=3)
    L = length(v)
    if L == 0
        return "Int[ ]"
    elseif L <= 2n
        return "Int[" * join(v, ' ') * "]"
    else
        head = join(v[1:n], ' ')
        tail = join(v[end-n+1:end], ' ')
        return "Int[$head … $tail]"
    end
end

function _preview(v::Vector{Float64}; n::Int=3,
    digits::Int=2, scale::Union{Nothing,Float64}=nothing,
    suffix::String="")

    L = length(v)

    if L == 0
        return "$(Float64)[ ]"
    end

    vals = scale === nothing ? v : v ./ scale

    fmt(x) = string(round(x; digits=digits))

    if L <= 2n
        body = join((fmt(vals[i]) for i in 1:L), ' ')
    else
        head = join((fmt(vals[i]) for i in 1:n), ' ')
        tail = join((fmt(vals[i]) for i in L-n+1:L), ' ')
        body = "$head … $tail"
    end

    return "$Float64[$body]" * suffix
end

function Base.show(io::IO, ::MIME"text/plain", d::domprop)

    Ni = length(d.bdmode)

    nfar  = count(==(UInt8(0)), d.bdmode)
    nnear = count(==(UInt8(1)), d.bdmode)
    nintp = count(==(UInt8(2)), d.bdmode)

    total_near_pairs = length(d.bdneark)
    total_intp_pairs = length(d.bdintpk)

    Naux = d.Lᵢₙ - 1
    expected_intp_ptr_len = length(d.bdclosest) * Naux + 1

    println(io, "domain properties:")
    println(io, "  N:          ", d.N)
    println(io, "  del:        ", d.del)
    println(io, "  del_near:   ", d.del_near)
    println(io, "  del_intp:   ", d.del_intp)

    println(io)
    println(io, "  hmap       :  Dict  (", length(d.hmap), " entries)")
    println(io, "  tgtpts     :  Matrix{Float64} (", size(d.tgtpts, 1), "×", size(d.tgtpts, 2), ")")
    println(io, "  prepts     :  Matrix{Int}     (", size(d.prepts, 1), "×", size(d.prepts, 2), ")")
    println(io, "  invpts     :  Matrix{Float64} (", size(d.invpts, 1), "×", size(d.invpts, 2), ")")

    println(io)
    println(io, "  Bookkeeping:")
    println(io, "    pthgo:    ", _preview(d.pthgo))
    println(io, "    nspsz:    ", _preview(d.nspsz))

    println(io)
    println(io, "  Boundary-DLP target classification:")
    println(io, "    bdmode   :  Vector{UInt8} (", Ni, " interior targets)")
    println(io, "      mode 0 far/direct        : ", nfar)
    println(io, "      mode 1 near/smoothed     : ", nnear)
    println(io, "      mode 2 interp/full-DLP   : ", nintp)

    println(io)
    println(io, "  Mode-1 near-panel storage:")
    println(io, "    bdnearptr :  Vector{Int}     (", length(d.bdnearptr), ")")
    println(io, "    bdneark   :  Vector{Int}     (", length(d.bdneark), " panel entries)")
    println(io, "    bdneart   :  Vector{Float64} (", length(d.bdneart), " projections)")
    println(io, "    total near pairs: ", total_near_pairs)

    println(io)
    println(io, "  Mode-2 interpolation target storage:")
    println(io, "    bdclosest :  Vector{Int}     (", length(d.bdclosest), ")")
    println(io, "    bdt       :  Vector{Float64} (", length(d.bdt), " projections)")
    println(io, "    bddist    :  Vector{Float64} (", length(d.bddist), " distances)")

    println(io)
    println(io, "  Mode-2 auxiliary interpolation storage:")
    println(io, "    Lᵢₙ       :  ", d.Lᵢₙ)
    println(io, "    bdxvals   :  Vector{Float64} (", length(d.bdxvals), " interpolation distances)")
    println(io, "    bdintpptr :  Vector{Int}     (", length(d.bdintpptr), ")")
    println(io, "    bdintpk   :  Vector{Int}     (", length(d.bdintpk), " panel entries)")
    println(io, "    bdintpt   :  Vector{Float64} (", length(d.bdintpt), " projections)")
    println(io, "    total auxiliary near pairs: ", total_intp_pairs)
    println(io, "    expected bdintpptr length: ", expected_intp_ptr_len)

    if !isempty(d.bdclosest)
        println(io)
        println(io, "  Boundary interpolation preview:")
        println(io, "    bdclosest: ", _preview(d.bdclosest))
        println(io, "    bdt:       ", _preview(d.bdt))
        println(io, "    bddist:    ", _preview(d.bddist; scale = 1e-3, suffix = " × 10^{-3}"))
    end

    if !isempty(d.bdxvals)
        println(io)
        println(io, "  Auxiliary interpolation preview:")
        println(io, "    bdxvals:   ", _preview(d.bdxvals; scale = 1e-3, suffix = " × 10^{-3}"))
    end

    if !isempty(d.bdintpk)
        println(io)
        println(io, "  Auxiliary near-panel preview:")
        println(io, "    bdintpk:   ", _preview(d.bdintpk))
        println(io, "    bdintpt:   ", _preview(d.bdintpt))
    end

    return nothing
end

function memory_report(dp::domprop)
    fields = fieldnames(typeof(dp))

    rows = [(name, Base.summarysize(getfield(dp, name))) for name in fields]
    total = Base.summarysize(dp)

    @printf("domprop memory report\n")
    @printf("---------------------\n")
    @printf("%-20s %12s %12s\n", "field", "MiB", "percent")

    for (name, bytes) in sort(rows; by = x -> x[2], rev = true)
        @printf("%-20s %12.4f %11.2f%%\n",
            String(name),
            bytes / 1024^2,
            100 * bytes / total)
    end

    @printf("---------------------\n")
    @printf("%-20s %12.4f %11.2f%%\n",
        "TOTAL", total / 1024^2, 100.0)

    return nothing
end