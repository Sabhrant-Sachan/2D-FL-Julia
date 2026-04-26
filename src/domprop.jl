"""
Distance-to-boundary record for each target point `x`
[d,l,t] such that d = | x - γ_l(t) |
"""
struct distrec
    d::Float64   # distance to boundary of the domain
    l::Int       # boundary patch index achieving the min distance
    t::Float64   # local parameter on that boundary curve
end

@inline packkey(i::Int, j::Int) =(UInt64(i) << 32) | UInt64(j)

mutable struct domprop
    N::Int
    del::Float64
    
    delclsbd::Float64
    delfarbd::Float64

    # lookups: packkey(global target col, integration patch k) -> column in prepts
    hmap::Dict{UInt64,Int}

    # grouped data: M = d.Npat, Mbd = |dom.kd|
    # 2 × (M*N^2 + Mbd*N), representing number of tgt pts
    # 1st and 2nd row are x and y points in real space of tgt pt resp.
    # The patch in which the point is present is given by 
    #
    # ℓ = ceil(i/Np), if i<=M*Np (interioir target point)
    # k₀ = ceil((i-M*Np)/N), ℓ = d.kd(k₀), if i>M*Np (bd tgt point)
    #
    # local linear index is given by 
    #
    # j = i - (ℓ - 1)*Np, if i<=M*Np (j₂, j₁ = divrem(j-1, N) .+ 1)
    #       t = ( cospi((2j₁-1)/(2N)), cospi((2j₂-1)/(2N)) )
    # j = i - M*Np - (k₀ - 1)*N, if i>M*Np
    #                   t = cospi((2j-1)/(2N))
    #
    tgtpts::Matrix{Float64}  

    # 2 × (M*N^2 + Mbd*N + Tnsp), representing number of pre pts
    # 1st and 2nd row are column index of tgt pt and int. patch k resp.
    prepts::Matrix{Int}   

    # 2 × Tnsp (inverse images for near-singulars)
    invpts::Matrix{Float64}  

    # A Matrix with M*N^2 columns in which 1st row represents which interior 
    #target pt is close to the boundary. (target-tgt boundary-bd boolean-b matrix-m)
    # tgtbdbm[1, i] = true implies ith interior tgt point close to bd
    # tgtbdbm[2, i] = true implies ith interior tgt point far from bd
    tgtbdbm::Matrix{Bool}    

    # For all points which are close to the bd, store distrec struct. 
    distpts::Vector{distrec} 

    # bookkeeping
    pthgo::Vector{Int}   # length M+1 (starts in pre for each patch; last entry=boundary start)
    nspsz::Vector{Int}   # length M+1 (near-singular counts per patch; last entry=boundary block)

    function domprop(N::Integer, del::Float64, delclsbd::Float64, dom::D) where {D<:abstractdomain}
        #JNote that Julia is a column major language

        M = dom.Npat            # number of patches
        Np = N^2                # points per patch (interior)

        Mbd = length(dom.kd)    # number of patches touching the boundary

        nint = M * Np           # Total points in interior
        nbdy = Mbd * N          # Total points on boundary

        # --- 1) Build target points (interior first then boundary) ---
        tgtpts = Matrix{Float64}(undef, 2, nint + nbdy)

        #A Column vector or just vector in Julia
        z = Vector{Float64}(undef, N)

        for j in 1:N z[j] = cospi((2*j - 1)/(2N)) end

        #===
        For row vectors x,y in matlab with length nx and ny
        [Y,X]=meshgrid(y,x) is equivalent to the following:
        Y = ones(nx,1) .* y and X = x' .* ones(1,ny).
        But notice that here z is a column vector.
        For column vector x and y, in matlab we could also 
        write Y = repmat(y',nx,1) and X = repmat(x,1,ny)
        NOTE : Use reshape if data is complex!
        The line "[zy,zx] = meshgrid(z);" converts to
        ===#
        zx = repeat(z, 1, N)
        zy = repeat(z', N, 1)

        #Collection of all interioir points
        @views Xint = tgtpts[:, 1:nint]
        #Collection of all boundary points
        @views Xbdy = tgtpts[:, nint+1:nint + nbdy]

        # Preallocate X, Y and x,y
        X = similar(zx)
        Y = similar(zy)
        x = similar(z)
        y = similar(z)

        i, ℓ = 1, 1
        for k in 1:M
            mapxy!(X, Y, dom, zx, zy, k)      # X,Y are N×N
            # fill N^2 targpt records 
            @inbounds for j in 1:Np
                Xint[1, i] = X[j]
                Xint[2, i] = Y[j]
                i += 1
            end
        end

        # boundary Chebyshev nodes on each boundary curve
        for k in dom.kd
            gamx!(x, dom, z, k)  # length N
            gamy!(y, dom, z, k)  # length N
            @inbounds for j in 1:N
                Xbdy[1, ℓ] = x[j]
                Xbdy[2, ℓ] = y[j]
                i += 1
                ℓ += 1
            end
        end

        #------------ 2. Count near-singulars to preallocate ------------
        # For each integration patch k, we’ll gather which *target* indices are
        # near-singular w.r.t. k. 
        pthgo = zeros(Int, M + 1)
        nspsz = zeros(Int, M + 1)

        nspsz_int = zeros(Int, M)
        nspsz_bd  = zeros(Int, M)

        
        near_off_int = Vector{Int}(undef, M+1)
        near_off_bd  = Vector{Int}(undef, M+1)
        writepos_int = Vector{Int}(undef, M+1)
        writepos_bd  = Vector{Int}(undef, M+1)

        # Extended quad for each k
        Q = dom.Qpts    # 8×M (each column is [P1...P4]’)
        Qe = similar(Q)
        @inbounds for k in 1:M
            @views V1 = Qe[:, k]
            @views V2 = Q[:, k]
            extendqua!(V1, V2, del)
        end

        insideint = Vector{Bool}(undef, nint)
        insidebdy = Vector{Bool}(undef, nbdy)

        @inbounds for k in 1:M
            P1x, P1y = Qe[1, k], Qe[2, k]
            P2x, P2y = Qe[3, k], Qe[4, k]
            P3x, P3y = Qe[5, k], Qe[6, k]
            P4x, P4y = Qe[7, k], Qe[8, k]

            ptinqua!(insideint, Xint, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)
            ptinqua!(insidebdy, Xbdy, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)

            cnt = 0
            @inbounds for i in 1:nint
                if insideint[i] && cld(i, Np) != k
                    if isnsp(dom, tgtpts[1, i], tgtpts[2, i], k, del)
                        cnt += 1
                    end
                end
            end

            nspsz_int[k] = cnt

            cnt = 0
            @inbounds for i in 1:nbdy
                ti = nint + i
                kₒ = cld(i, N)
                if insidebdy[i] && dom.kd[kₒ] != k
                    if isnsp(dom, tgtpts[1, ti], tgtpts[2, ti], k, del)
                         cnt += 1
                    end
                end
            end
            nspsz_bd[k] = cnt

        end

        near_off_int[1] = 1
        near_off_bd[1] = 1

        @inbounds for k in 1:M
            near_off_int[k+1] = near_off_int[k] + nspsz_int[k]
            near_off_bd[k+1]  = near_off_bd[k]  +  nspsz_bd[k]
        end

        Tnsp_int = near_off_int[M+1] - 1
        Tnsp_bd = near_off_bd[M+1] - 1

        near_int_all = Vector{Int}(undef, Tnsp_int)
        near_bd_all = Vector{Int}(undef, Tnsp_bd)

        Tnsp = Tnsp_int + Tnsp_bd

        copy!(writepos_int, near_off_int)
        copy!( writepos_bd, near_off_bd )

        @inbounds for k in 1:M
            P1x, P1y = Qe[1, k], Qe[2, k]
            P2x, P2y = Qe[3, k], Qe[4, k]
            P3x, P3y = Qe[5, k], Qe[6, k]
            P4x, P4y = Qe[7, k], Qe[8, k]

            ptinqua!(insideint, Xint, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)
            ptinqua!(insidebdy, Xbdy, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)

            @inbounds for i in 1:nint
                if insideint[i] && cld(i, Np) != k
                    if isnsp(dom, tgtpts[1, i], tgtpts[2, i], k, del)
                        pos = writepos_int[k]
                        near_int_all[pos] = i
                        writepos_int[k] = pos + 1
                    end
                end
            end

            @inbounds for i in 1:nbdy
                ti = nint + i
                kₒ = cld(i, N)
                if insidebdy[i] && dom.kd[kₒ] != k
                    if isnsp(dom, tgtpts[1, ti], tgtpts[2, ti], k, del)
                        pos = writepos_bd[k]
                        near_bd_all[pos] = ti
                        writepos_bd[k] = pos + 1
                    end
                end
            end

        end

        #Now updating pthgo and nspsz
        @inbounds for k in 1:M
            nspsz[k] = nspsz_int[k]
        end

        nspsz[M+1] = sum(nspsz_bd)

        startpth = 1
        @inbounds for k in 1:M
            pthgo[k] = startpth
            startpth += Np + nspsz_int[k]   # Np interior + near-singulars for patch k
        end
        pthgo[M+1] = startpth

        #-------------------- 3) Preallocate prepts / invpts  --------------------
        # Layout of prepts:
        # For each k = 1..M, a block of Np interior points (with integration k),
        # followed by nspsz[k] near-singulars (with integration k).
        # Then one big block of boundary targets (nbd) with integration l=tgt.l,
        # followed by all boundary near-singulars (distributed over k).
        hmap = Dict{UInt64,Int}()
        sizehint!(hmap, Tnsp)  # number of pairs I expect to insert

        prepts = Matrix{Int}(undef, 2, nint + nbdy + Tnsp)
        invpts = Matrix{Float64}(undef, 2, Tnsp) 
        # # ------------- 4) Fill prepts/ invpts and hash maps -------------

        Nr = dom.pths[M].reg #Total number of regions in domain, 5 for disc 

        #Adding interioir pts + near singular
        ℓ, ℓi = 1, 1

        ftbs = Vector{FTable}(undef, Nr)

        for r in 1:Nr
            ftbs[r] = inFTable(10_001)     
            fill_FTable!(ftbs[r], dom, r)  # prefill for this region
        end

        for j in 1:M

            for i in 1:Np

                prepts[1, ℓ] = (j - 1) * Np + i
                prepts[2, ℓ] = j 

                ℓ += 1

            end

            for t in near_off_int[j]:(near_off_int[j+1]-1)

                i = near_int_all[t]

                tpl = cld(i, Np)

                prepts[1, ℓ] = i
                prepts[2, ℓ] = j

                hmap[packkey(i, j)] = ℓ
            
                #region of the jth patch 
                ξ = dom.pths[j].reg

                #We now compute the inverse of the near singular 
                #points with respect to patch j, if ith target point
                #is in the same "region" as the jth patch, then use
                #mapinv2, otherwise use mapinv
                if ξ == dom.pths[tpl].reg

                    tpj = i - (tpl - 1)*Np

                    q, r = divrem(tpj - 1, N)

                    #j1, j2 = r + 1, q + 1 

                    #t1 = cospi((2*r+1)/(2*N)); t2 = cospi((2*q+1)/(2*N));

                    Zx, Zy = mapinv2(dom, z[r+1], z[q+1], tpl, j)

                else
                    tpx = tgtpts[1, i]
                    tpy = tgtpts[2, i]
                    ftb = ftbs[ξ]

                    Zx, Zy = mapinv(ftb, dom, tpx, tpy, j)
                end

                Zx = abs(Zx - 1) < 1e-14 ? 1 :
                     abs(Zx + 1) < 1e-14 ? -1 : Zx

                Zy = abs(Zy - 1) < 1e-14 ? 1 :
                     abs(Zy + 1) < 1e-14 ? -1 : Zy

                invpts[1, ℓi] = Zx
                invpts[2, ℓi] = Zy

                ℓ, ℓi = ℓ + 1, ℓi + 1

            end

        end

        #adding boundary points + near singular    
        for i in 1:N*Mbd

            prepts[1, ℓ] = M*Np+i
            prepts[2, ℓ] = cld(i, N)

            ℓ += 1

        end

        for j = 1:M

            for t in near_off_bd[j]:(near_off_bd[j+1]-1)

                i = near_bd_all[t]

                k₀ = cld(i - M*Np, N) 

                tpl = dom.kd[k₀]

                prepts[1, ℓ] = i
                prepts[2, ℓ] = j 

                hmap[packkey(i, j)] = ℓ

                #region of the jth patch 
                ξ = dom.pths[j].reg

                if ξ == dom.pths[tpl].reg

                    tpj = i - M*Np - (k₀ - 1)*N

                    Zx, Zy = mapinv2(dom, 1.0, z[tpj], tpl, j)

                else
                    tpx = tgtpts[1, i]
                    tpy = tgtpts[2, i]
                    ftb = ftbs[ξ]

                    Zx, Zy = mapinv(ftb, dom, tpx, tpy, j)

                end

                Zx = abs(Zx - 1) < 1e-14 ? 1 :
                     abs(Zx + 1) < 1e-14 ? -1 : Zx

                Zy = abs(Zy - 1) < 1e-14 ? 1 :
                     abs(Zy + 1) < 1e-14 ? -1 : Zy

                invpts[1, ℓi] = Zx
                invpts[2, ℓi] = Zy

                ℓ, ℓi = ℓ + 1, ℓi + 1
            end

        end

        # ------------ 5) Distance-to-boundary for *interior* targets ------------
        delfarbd = 0.4          # Default. 

        # tgtbdbm has length M*Np, assume no point is close to bd 
        # and all points are far from the bd
        tgtbdbm = Matrix{Bool}(undef, 2, M * Np)

        #row 1: True = close to boundary
        #row 2: True = far from boundary
        @inbounds for j in 1:M*Np
            tgtbdbm[1, j] = false   # all pts are not close
            tgtbdbm[2, j] = true    # all pts are far
        end

        Qbd = dom.Qptsbd
        PE = similar(Qbd)
        W = Vector{Float64}(undef, 8)

        @inbounds for i in 1:Mbd
            @views V = PE[:, i]
            @views U = Qbd[:, i]
            extendqua!(V, U, delclsbd)
            extendqua!(W, U, delfarbd)

            ptinqua!(insideint, Xint, V[1], V[2], V[3], V[4], V[5], V[6], V[7], V[8])
            @inbounds for j in 1:nint
                if insideint[j]
                    tgtbdbm[1, j] = true
                end
            end

            ptinqua!(insideint, Xint, W[1], W[2], W[3], W[4], W[5], W[6], W[7], W[8])
            @inbounds for j in 1:nint
                if insideint[j]
                    tgtbdbm[2, j] = false
                end
            end
        end

        nclose = sum(@view tgtbdbm[1, :])
        distpts = Vector{distrec}(undef, nclose)
        bdx = Vector{Int}(undef, Mbd)

        ℓbd = 0
        #fill distpts
        for i in 1:nint

            # Point in inside δclsbd quadrilateral,
            # but not actually be close to bd
            if tgtbdbm[1, i]

                ℓbdx = 0

                x1 = tgtpts[1, i]
                x2 = tgtpts[2, i]

                for k in 1:Mbd

                    P1x, P1y = PE[1, k], PE[2, k]
                    P2x, P2y = PE[3, k], PE[4, k]
                    P3x, P3y = PE[5, k], PE[6, k]
                    P4x, P4y = PE[7, k], PE[8, k]

                    @views PQ = PE[:, k]
                    if ptinqua(x1, x2, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)
                        ℓbdx += 1
                        bdx[ℓbdx] = dom.kd[k]
                    end

                end

                # choose l with smallest distance by a quick 1D search
                ℓ = bdx[1]

                dmin, _ = GSS(dist_gam, -1.0, 1.0, 1e-8, dom, x1, x2, ℓ)

                for ℓℓ in 1:ℓbdx

                    k = bdx[ℓℓ]

                    dmin2, _ = GSS(dist_gam, -1.0, 1.0, 1e-8, dom, x1, x2, k)

                    if dmin2 < dmin
                        dmin = dmin2
                        ℓ = k
                    end

                end

                if dmin < delclsbd

                    ℓbd += 1

                    # refine location using your GSS+Bis
                    _, t0 = GSS(dist_gam, -1.0, 1.0, 1e-15, dom, x1, x2, ℓ)

                    tloc = Bis(dist_gamBis, -1.0, 1.0, t0, 64, dom, x1, x2, ℓ)

                    distpts[ℓbd] = distrec(dmin, ℓ, tloc)

                else

                    tgtbdbm[1, i] =  false

                    nclose = nclose - 1

                    resize!(distpts, nclose)

                end
                
           end

        end

        new(N, del, delclsbd, delfarbd, hmap, tgtpts, prepts,
            invpts, tgtbdbm, distpts, pthgo, nspsz)

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
Plot the target boundary boolean matrix on top of the domain drawing.
That is, plots points which are close to the (black) and far from the bd (blue)
Also draws the (black) quadrilateral for patch `k` and its δcls–extended quad.
and  draws the (blue) quadrilateral for patch `k` and its δfar–extended quad.
"""
function plotbm(dp::domprop, d::abstractdomain)
    fig, ax = draw(d)

    nint = d.Npat * dp.N^2
    ms = 9

    # Preallocate index buffers (max size nint)
    iblk = Vector{Int}(undef, nint)
    iblu = Vector{Int}(undef, nint)
    nblk = 0
    nblu = 0

    # Build index lists
    @inbounds for i in 1:nint
        if dp.tgtbdbm[1, i]
            nblk += 1
            iblk[nblk] = i
        end
        if dp.tgtbdbm[2, i]
            nblu += 1
            iblu[nblu] = i
        end
    end

    # Scatter black points in ONE call
    if nblk > 0
        xs = Vector{Float64}(undef, nblk)
        ys = Vector{Float64}(undef, nblk)
        @inbounds for j in 1:nblk
            i = iblk[j]
            xs[j] = dp.tgtpts[1, i]
            ys[j] = dp.tgtpts[2, i]
        end
        scatter!(ax, xs, ys; marker=:circle, markersize=ms, strokewidth=0, color=:black)
    end

    # Scatter blue points in ONE call
    if nblu > 0
        xs = Vector{Float64}(undef, nblu)
        ys = Vector{Float64}(undef, nblu)
        @inbounds for j in 1:nblu
            i = iblu[j]
            xs[j] = dp.tgtpts[1, i]
            ys[j] = dp.tgtpts[2, i]
        end
        scatter!(ax, xs, ys; marker=:circle, markersize=ms, strokewidth=0, color=:blue)
    end

    Mbd = length(d.kd)
    Pext = Vector{Float64}(undef, 8)

    ix = [1, 3, 5, 7, 1] 
    iy = [2, 4, 6, 8, 2]

    for k in 1:Mbd 
        # -- draw quadrilateral for patch k -- 
        @views P = d.Qptsbd[:, k] 
        # 8-vector: [x1,y1,x2,y2,x3,y3,x4,y4] 
        # -- draw δcls-extended quadrilateral --
        extendqua!(Pext, P, dp.delclsbd) 
        lines!(ax, Pext[ix], Pext[iy], color=:black, linewidth=2) 
        # -- draw δfar-extended quadrilateral -- 
        #extendqua!(Pext, P, dp.delfarbd) 
        #lines!(ax, Pext[ix], Pext[iy], color=:blue, linewidth=2) 
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
    @views P  = d.Qpts[:, k]      # 8-vector: [x1,y1,x2,y2,x3,y3,x4,y4]

    Pext = Vector{Float64}(undef, 8)

    ix = [1, 3, 5, 7, 1];  iy = [2, 4, 6, 8, 2]
    #plot!(plt, P[ix], P[iy], color=:black, lw=2, label="")

    # -- draw δ-extended quadrilateral --
    extendqua!(Pext, P, dp.del)
    lines!(ax, Pext[ix], Pext[iy], color=:black, linewidth=2)

    # convenient sizes
    M   = d.Npat
    Np  = dp.N^2
    Mbd = length(d.kd)
    nbd = Mbd * dp.N

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
    #   - a block of all boundary targets (length nbd)
    #   - then all boundary near–singulars (various k).
    start_bd_block = dp.pthgo[M+1] + nbd

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
 ========== Dispatch 1: f is given in closed-form ==========
The code below uses Makie’s reactive Observables to unify colors across many surfaces 
in a single pass. We create `clims = Observable((lo, hi))` and pass it directly to each
`surface!` via `colorrange = clims`. Plot attributes that receive an Observable keep a 
live reference to it: every surface “subscribes” to `clims`. Inside the loop, for each 
patch we compute Z once, widen 
          `clims[] = (min(current_lo, minimum(Z)), max(current_hi, maximum(Z)))`,
and Makie automatically re-colors all previously drawn surfaces—no re-plotting, no second 
pass, and no giant stitched matrices. This differs from reading snapshots like `clims[][1]`, 
which returns plain values that won’t update; only the Observable itself is reactive. We 
don’t need `lift`/`@lift` here because we bind the Observable directly to a plot attribute 
(not creating a derived Observable). If one prefer's non-reactive alternatives, use a two-pass 
approach (first gather global min/max, then plot with a fixed `colorrange`). Outside the plot
environment 
        lift(f, obs...) = create a new observable that updates whenever any input changes.
        on(obs) do v ... end = subscribe and run side-effects on changes.
For instance: clims = Observable((0.0, 1.0)); x = Observable(0.0); y = Observable(0.0);
                  S = Observable(0.0); 

on(clims) do (lo, hi)
    x[] = lo
    y[] = hi
    S[] = lo + hi
end

clims[] = (-1, 10); #x, y and S automatically update

x[]  # -> -1; y[]  # -> 10; S[]  # -> 9

Immutable payloads (numbers, tuples): changing them always means assign, so Observables 
auto-notify.Mutable payloads (arrays, dicts): mutations don’t notify—either reassign 
A[] = A[] (or a copy), or call notify(A) after you mutate. Use reassignment (A[] = new)
or notify(A) after in-place edits to trigger listeners.
"""
function plotfunc(dp::domprop, d::abstractdomain, f::Function)

    # Closed Chebyshev nodes on [-1,1] 
    z = cospi.((0:4*dp.N) ./ (4*dp.N))
    zx = repeat(z, 1, 4*dp.N+1)
    zy = repeat(z', 4*dp.N+1, 1)

    # figure + 3D axis with labels
    fig = Figure(size = (650, 650))

    ax  = Axis3(fig[1, 1];
        aspect = :equal,
        xlabel = latexstring("\$x\$"),
        ylabel = latexstring("\$y\$"),
        zlabel = latexstring("\$f(x,y)\$"),
        xlabelsize = 22, ylabelsize = 22, zlabelsize = 22,
        xticklabelsize = 18, yticklabelsize = 18, zticklabelsize = 18
    )

    M = d.Npat

    clims = Observable((0.0, 1.0))  # placeholder; will be updated from first patch

    X, Y = mapxy(d, zx, zy, 1)

    Z = real.(f.(X, Y))

    clims[] = extrema(Z)

    surface!(ax, X, Y, Z; color=Z, colorrange=clims)

    for P in 2:M
        X, Y = mapxy(d, zx, zy, P)

        Z = real.(f.(X, Y))

        # widen the shared colorrange and all already-plotted surfaces will update
        c0, c1 = clims[]
        z0, z1 = extrema(Z)

        clims[] = (min(c0, z0), max(c1, z1))

        # plot this patch using the shared colorrange (bound to the Observable)
        surface!(ax, X, Y, Z; color=Z, colorrange=clims)
    end

    display(GLMakie.Screen(), fig)

    return (fig, ax)
end

"""
This function plots a surface when the function values are already known on each patch’s 
n×n open Chebyshev grid. It computes the open nodes z = cospi((2j-1)/(2n)), maps them 
to each patch with mapxy(d, zx, zy, P), reshapes the corresponding slice of f to n×n,
and draws one surface! per patch. All surfaces share a single, reactive color scale 
via a common Observable clims, so the colormap is same across all the patches. Axes are 
styled with Axis3 aspect = :equal,LaTeX labels to make the plot beautiful.
If function is called with interpolate=true, it switches to interpolation on denser
closed Chebyshev mesh and then plotting it.
"""

function plotfunc(dp::domprop, d::abstractdomain, f::AbstractVector; interpolate::Bool = false)
    # grid sizes
    n  = dp.N
    M  = d.Npat
    Np = n^2

    # figure + axis
    fig = Figure(size = (650, 650))
    ax  = Axis3(fig[1, 1];
        aspect = :equal,
        xlabel = latexstring("\$x\$"),
        ylabel = latexstring("\$y\$"),
        zlabel = latexstring("\$f(x,y)\$"),
        xlabelsize = 22, ylabelsize = 22, zlabelsize = 22,
        xticklabelsize = 18, yticklabelsize = 18, zticklabelsize = 18
    )

    if !interpolate
        # ---------- NO INTERPOLATION: plot the raw samples ----------
        # Open Chebyshev (Gauss) nodes: z_j = cos( (2j-1)π / (2n) ), j=1..n
        z = cospi.((2 .* (1:n) .- 1) ./ (2n))
        zx = repeat(z, 1, n)         
        zy = repeat(z', n, 1)        

        # Single-pass, unified colormap via Observable
        X, Y = mapxy(d, zx, zy, 1)
        Z = reshape(@view(f[1:Np]), n, n) |> real
        clims = Observable(extrema(Z))
        surface!(ax, X, Y, Z; color=Z, colorrange=clims)
        scatter!(ax, X, Y, Z;
            marker=:circle,
            markersize=9,
            color=:black)

        for P in 2:M
            X, Y = mapxy(d, zx, zy, P)
            Z = reshape(@view(f[(P-1)*Np+1:P*Np]), n, n) |> real

            (lo, hi) = clims[]
            (zlo, zhi) = extrema(Z)
            clims[] = (min(lo, zlo), max(hi, zhi))

            surface!(ax, X, Y, Z; color=Z, colorrange=clims)
            scatter!(ax, X, Y, Z;
                marker=:circle,
                markersize=9,
                color=:black)
        end

        display(GLMakie.Screen(),fig)

        return (fig, ax)
    end

    # closed Chebyshev nodes on [-1,1] for plotting
    z  = cospi.((0:4*n) ./ (4*n))
    zx = repeat(z, 1, 4*n + 1)
    zy = repeat(z', 4*n + 1, 1)
    T  = ChebyTN(n, z)

    # Compute Chebyshev coefficients
    # https://www.fftw.org/fftw3_doc/1d-Real_002deven-DFTs-_0028DCTs_0029.html
    chebcoef = similar(f, M * Np)
    for k in 1:M
        idx = (k-1)*Np + 1 : k*Np
        fv  = reshape(@views(f[idx]), n, n)

        # columns first (dim = 2), DCT-II
        c = (1 / n) .* FFTW.r2r(fv, FFTW.REDFT10, 2)
        c[:, 1] ./= 2

        # rows next (dim = 1), DCT-II
        c = (1 / n) .* FFTW.r2r(c,  FFTW.REDFT10, 1)
        c[1, :] ./= 2

        chebcoef[idx] .= vec(c)
    end

    # first patch → initialize clims and plot
    c = reshape(@view(chebcoef[1:Np]), n, n)
    Z = T * c * transpose(T)
    X, Y = mapxy(d, zx, zy, 1)

    clims = Observable(extrema(real(Z)))
    surface!(ax, X, Y, real.(Z); color = real.(Z), colorrange = clims)

    # remaining patches (single pass, widen clims as we go)
    for P in 2:M
        c = reshape(@view(chebcoef[(P-1)*Np + 1 : P*Np]), n, n)
        Z  = T * c * transpose(T)
        X, Y = mapxy(d, zx, zy, P)

        z0, z1 = extrema(real(Z))
        c0, c1 = clims[]
        clims[] = (min(c0, z0), max(c1, z1))

        surface!(ax, X, Y, real.(Z); color = real.(Z), colorrange = clims)
    end

    display(GLMakie.Screen(), fig)
    return (fig, ax)
end

# ---------------------------------------------------------------
# u given in closed form; enforce u=0 on the boundary of domain.
# (Closed Chebyshev grid includes the boundary, so we can zero edges.)
# ---------------------------------------------------------------
function plotu(dp::domprop, d::abstractdomain, u::Function)
    n = dp.N
    ℓ = 4n+1
    z  = cospi.((0:4n) ./ (4n))        # closed Chebyshev nodes ∈ [-1,1]
    zx = repeat(z, 1, ℓ)
    zy = repeat(z', ℓ, 1)

    fig = Figure(size = (650, 650))
    ax  = Axis3(fig[1, 1];
        aspect = :equal,
        xlabel = latexstring("\$x\$"),
        ylabel = latexstring("\$y\$"),
        zlabel = latexstring("\$u(x,y)\$"),
        xlabelsize = 22, ylabelsize = 22, zlabelsize = 22,
        xticklabelsize = 18, yticklabelsize = 18, zticklabelsize = 18
    )

    # First patch → initialize clims after enforcing boundary zeros
    X, Y = mapxy(d, zx, zy, 1)
    # enforce homogeneous Dirichlet boundary
    # v varying u= +1 fixed
    if isbdpatch(d, 1)
        Z = u.(X[2:ℓ, :], Y[2:ℓ, :])
        Z = vcat(zeros(Float64, ℓ)', Z)
    else
        Z = u.(X, Y)
    end

    clims = Observable(extrema(Z))
    surface!(ax, X, Y, Z; color = Z, colorrange = clims)

    # Remaining patches (single pass; widen shared colorrange)
    for P in 2:d.Npat
        X, Y = mapxy(d, zx, zy, P)
        if isbdpatch(d, P)
            Z = u.(X[2:ℓ, :], Y[2:ℓ, :])
            Z = vcat(zeros(Float64, ℓ)', Z)
        else
            Z = u.(X, Y)
        end

        (lo, hi) = clims[]; (zlo, zhi) = extrema(Z)
        clims[] = (min(lo, zlo), max(hi, zhi))

        surface!(ax, X, Y, Z; color = Z, colorrange = clims)
    end

    display(GLMakie.Screen(), fig)
    return (fig, ax)
    #save("myplot.png",fig; px_per_unit=10)
end

#u is being mutated inside!
function plotu(dp::domprop, d::abstractdomain, u::AbstractVector, s::Real; interpolate::Bool=true)

    # grid sizes
    n = dp.N
    M = d.Npat
    Np = n^2
    #https://docs.makie.org/dev/explanations/colors
    cmap = :jet1 #:hsv

    # figure + axis
    fig = Figure(size= (850, 650))
    ax = Axis3(fig[1, 1];
        #aspect=:equal,
        xlabel=latexstring("\$x\$"),
        ylabel=latexstring("\$y\$"),
        zlabel=latexstring("\$u(x,y)\$"),
        xlabelsize=22, ylabelsize=22, zlabelsize=22,
        xticklabelsize=18, yticklabelsize=18, zticklabelsize=18
    )

    if !interpolate
        # ---------- NO INTERPOLATION: plot the raw samples ----------
        # Open Chebyshev (Gauss) nodes: z_j = cos( (2j-1)π / (2n) ), j=1..n
        z = cospi.((2 .* (1:n) .- 1) ./ (2n))
        zx = repeat(z, 1, n)
        zy = repeat(z', n, 1)

        zxb = vcat(ones(Float64, n)', zx)
        zyb = vcat(z', zy)

        for i in 1:M*Np

            k = cld(i, Np)
            j = i - (k - 1) * Np

            j1 = mod(j - 1, n) + 1
            
            u[i] = u[i] * dfunc(d, k, 2 * sinpi((2 * j1 - 1) / (4 * n))^2, s)

        end

        # Single-pass, unified colormap via Observable
        if isbdpatch(d, 1)
            #Put boundary points together 
            X, Y = mapxy(d, zxb, zyb, 1)
            Z = reshape(@view(u[1:Np]), n, n) 
            Z = vcat(zeros(Float64, n)', Z)
        else
            X, Y = mapxy(d, zx, zy, 1)
            Z = reshape(@view(u[1:Np]), n, n) 
        end

        clims = Observable(extrema(Z))
        surface!(ax, X, Y, Z; color=Z, colorrange=clims, colormap = cmap, shading = NoShading)
        scatter!(ax, X, Y, Z;
            marker=:circle,
            markersize=9,
            color=:black,
        )

        for P in 2:M
            if isbdpatch(d, P)
                X, Y = mapxy(d, zxb, zyb, P)
                Z = reshape(@view(u[(P-1)*Np+1:P*Np]), n, n)
                Z = vcat(zeros(Float64, n)', Z)
            else
                X, Y = mapxy(d, zx, zy, P)
                Z = reshape(@view(u[(P-1)*Np+1:P*Np]), n, n)
            end

            (lo, hi) = clims[]
            (zlo, zhi) = extrema(Z)
            clims[] = (min(lo, zlo), max(hi, zhi))

            surface!(ax, X, Y, Z; color=Z, colorrange=clims, colormap = cmap, shading = NoShading)
            scatter!(ax, X, Y, Z;
                marker=:circle,
                markersize=9,
                color=:black,
            )
        end

        display(GLMakie.Screen(), fig)

        return (fig, ax)
    end

    # ---------- INTERPOLATION: plot the denser samples ----------
    # closed Chebyshev nodes on [-1,1] for plotting
    ℓ = 4n+1
    z = cospi.((0:4*n) ./ (4*n))
    #First row of matrix zx is [1,1,...,1,1]
    zx = repeat(z, 1, ℓ)
    zy = repeat(z', ℓ, 1)

    T = ChebyTN(n, z)

    # Compute Chebyshev coefficients
    # https://www.fftw.org/fftw3_doc/1d-Real_002deven-DFTs-_0028DCTs_0029.html
    chebcoef = similar(u, M * Np)
    for k in 1:M
        idx = (k-1)*Np+1:k*Np
        fv = reshape(@views(u[idx]), n, n)

        # columns first (dim = 2), DCT-II
        c = (1 / n) .* FFTW.r2r(fv, FFTW.REDFT10, 2)
        c[:, 1] ./= 2

        # rows next (dim = 1), DCT-II
        c = (1 / n) .* FFTW.r2r(c, FFTW.REDFT10, 1)
        c[1, :] ./= 2

        chebcoef[idx] .= vec(c)
    end

    # first patch → initialize clims 
    c = reshape(@view(chebcoef[1:Np]), n, n)
    Z = T * c * transpose(T)
    X, Y = mapxy(d, zx, zy, 1)
    for i in 2:ℓ
        Z[i, :] = Z[i, :] .* dfunc(d, 1, 2 * sinpi((i - 1) / (8 * n))^2, s)
    end

    if isbdpatch(d, 1)
        #Multiply the distance factor
        Z[1,:] = zeros(Float64, ℓ)
    else
        Z[1, :] = Z[1, :] .* dfunc(d, 1, 0.0, s)
    end

    clims = Observable(extrema(real(Z)))

    surface!(ax, X, Y, real.(Z); color=real.(Z), colorrange=clims, colormap = cmap, shading = NoShading)

    # remaining patches (single pass, widen clims as we go)
    for P in 2:M

        c = reshape(@view(chebcoef[(P-1)*Np+1:P*Np]), n, n)
        Z = T * c * transpose(T)
        X, Y = mapxy(d, zx, zy, P)

        #Multiply the distance factor
        for i in 2:ℓ
            Z[i, :] = Z[i, :] .* dfunc(d, P, 2 * sinpi((i - 1) / (8 * n))^2, s)
        end

        if isbdpatch(d, P)
            Z[1,:] = zeros(Float64, ℓ)
        else
            Z[1, :] = Z[1, :] .* dfunc(d, P, 0.0, s)
        end

        z0, z1 = extrema(real(Z))
        c0, c1 = clims[]
        clims[] = (min(c0, z0), max(c1, z1))

        surface!(ax, X, Y, real.(Z); color=real.(Z), colorrange=clims, colormap = cmap, shading = NoShading)
    end

    display(GLMakie.Screen(), fig)
    return (fig, ax)

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

function Base.show(io::IO, ::MIME"text/plain", d::domprop)
    # headline
    println(io, "domain properties:")
    println(io, "  N:        ", d.N)
    println(io, "  del:      ", d.del)
    println(io, "  delclsbd: ", d.delclsbd)
    println(io, "  delfarbd: ", d.delfarbd)

    # quick counts
    Nint   = length(d.distpts)                  # interior targets
    maplen = length(d.hmap)

    println(io)
    println(io, "  hmap:     Dict  (", maplen,  " entries)")
    println(io, "  tgtpts:   Matrix{Float64}  (", size(d.tgtpts,1), "×", size(d.tgtpts,2), ")")
    println(io, "  prepts:   Matrix{Int}      (", size(d.prepts,1), "×", size(d.prepts,2), ")")
    println(io, "  invpts:   Matrix{Float64}  (", size(d.invpts,1), "×", size(d.invpts,2), ")")
    println(io, "  tgtbdbm:  Matrix{bool}     (", size(d.tgtbdbm,1), "×", size(d.tgtbdbm,2), ")")

    println(io, "  distpts:  Vector{distrec}  (length ", Nint , ")")

    # bookkeeping vectors (show size + a short preview)
    println(io)
    println(io, "  pthgo:    ", _preview(d.pthgo))
    println(io, "  nspsz:    ", _preview(d.nspsz))

    return
end
