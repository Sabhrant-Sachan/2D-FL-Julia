@inline packkey(i::Int, k::Int) = (UInt64(i) << 32) | UInt64(k)

mutable struct domprop
    N::Int
    del::Float64
    del_bd::Float64

    # hmap[packkey(i,k)] = colp
    # i = global target column in tgtpts
    # k = actual volume integration patch index
    # colp = column in prepts
    hmap::Dict{UInt64,Int}

    # Ni is number of interior targets
    # hmap_bd[packkey(i,k)] = colpbd
    # i = interior target column, 1:Ni
    # k = actual boundary integration patch index, i.e. one of dom.kd
    # colpbd = column in prepts_bd / projpts
    hmap_bd::Dict{UInt64,Int}

    # Nb is number of boundary targets
    # 2 × (Ni + Nb)
    # First row is x coordinates, second y coordinate
    # columns 1:Ni are interior targets
    # columns Ni+1:Ni+Nb are boundary targets
    tgtpts::Matrix{Float64}
    # For an interior target i <= Ni:
    #     kt = cld(i, Np)
    #     ip = i - (kt - 1)*Np
    #     q, r = divrem(ip - 1, N)
    #     local coordinates are:
    #         t1 = cospi((2*(r+1)-1)/(2N))
    #         t2 = cospi((2*(q+1)-1)/(2N))
    #
    # For a boundary target i > Ni:
    #     ib  = i - Ni
    #     ibp = cld(ib, N)
    #     kt  = dom.kd[ibp]
    #     ip  = ib - (ibp - 1)*N
    #     local boundary coordinate is:
    #         t = cospi((2*ip - 1)/(2N))

    # 2 × (Ni + Nb + Tnsp)
    # row 1 = global target column i
    # row 2 = actual integration patch k
    prepts::Matrix{Int}

    # 2 × Tnsp_bd
    # row 1 = interior target column i
    # row 2 = actual boundary integration patch k
    prepts_bd::Matrix{Int}

    # 2 × Tnsp
    # inverse local coordinates for near-singular volume interactions
    invpts::Matrix{Float64}

    # length Tnsp_bd
    # projected local boundary parameter for near-singular boundary interactions
    projpts::Vector{Float64}

    # pthgo[k] = first column in prepts for actual volume patch k
    # pthgo[M+1] = first column of the regular boundary-target block in prepts
    pthgo::Vector{Int}

    # nspsz[k] = number of interior targets near volume patch k
    # nspsz[M+1] = total number of boundary targets near volume patches
    nspsz::Vector{Int}

    function domprop(N::Integer, del::Float64, del_bd::Float64, dom::D) where {D<:abstractdomain}
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
        # colpbd  column in prepts_bd and projpts
        #
        # nsp_i   number of interior targets near each volume patch, length M
        # nsp_b   number of boundary targets near each volume patch, length M
        # nsp_bd  number of interior targets near each boundary patch, length Mbd
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

        Qbd = dom.Qptsbd
        Qebd = similar(Qbd)

        @inbounds for kb in 1:Mbd
            @views extendqua!(Qebd[:, kb], Qbd[:, kb], del_bd)
        end

        in_i = Vector{Bool}(undef, Ni)
        in_b = Vector{Bool}(undef, Nb)

        nsp_i = zeros(Int, M)
        nsp_b = zeros(Int, M)
        nsp_bd = zeros(Int, Mbd)

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

        # Count interior targets near each boundary patch.
        @inbounds for kb in 1:Mbd
            k = dom.kd[kb]

            P1x, P1y = Qebd[1, kb], Qebd[2, kb]
            P2x, P2y = Qebd[3, kb], Qebd[4, kb]
            P3x, P3y = Qebd[5, kb], Qebd[6, kb]
            P4x, P4y = Qebd[7, kb], Qebd[8, kb]

            ptinqua!(in_i, Xi, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)

            cnt = 0
            for i in 1:Ni
                if in_i[i]
                    if isnspbd(dom, tgtpts[1, i], tgtpts[2, i], k, del_bd)
                        cnt += 1
                    end
                end
            end
            nsp_bd[kb] = cnt
        end

        Tnsp_i = sum(nsp_i)
        Tnsp_b = sum(nsp_b)
        Tnsp = Tnsp_i + Tnsp_b
        Tnsp_bd = sum(nsp_bd)

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

        hmap_bd = Dict{UInt64,Int}()
        sizehint!(hmap_bd, Tnsp_bd)

        prepts = Matrix{Int}(undef, 2, Ni + Nb + Tnsp)
        invpts = Matrix{Float64}(undef, 2, Tnsp)

        prepts_bd = Matrix{Int}(undef, 2, Tnsp_bd)
        projpts = Vector{Float64}(undef, Tnsp_bd)

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
        # 6. Fill prepts_bd and projpts for boundary integrations
        # ---------------------------------------------------------------------
        colpbd = 1

        @inbounds for kb in 1:Mbd
            k = dom.kd[kb]

            P1x, P1y = Qebd[1, kb], Qebd[2, kb]
            P2x, P2y = Qebd[3, kb], Qebd[4, kb]
            P3x, P3y = Qebd[5, kb], Qebd[6, kb]
            P4x, P4y = Qebd[7, kb], Qebd[8, kb]

            ptinqua!(in_i, Xi, P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y)

            for i in 1:Ni
                if in_i[i]
                    if isnspbd(dom, tgtpts[1, i], tgtpts[2, i], k, del_bd)
                        prepts_bd[1, colpbd] = i
                        prepts_bd[2, colpbd] = k
                        hmap_bd[packkey(i, k)] = colpbd

                        τ = projbd(dom, tgtpts[1, i], tgtpts[2, i], k)
                        projpts[colpbd] = snap1(τ)

                        colpbd += 1
                    end
                end
            end
        end

        return new(N, del, del_bd, hmap, hmap_bd, tgtpts,
                   prepts, prepts_bd, invpts, projpts, pthgo, nspsz)
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
Plot the near-singular interior target points to boundary patch `k`
on top of the domain drawing. Here `k` is the actual boundary patch index, 
i.e. `k ∈ d.kd`.
- Near-singular interior targets to boundary patch `k` are shown in blue.
- Also draws the δ_bd-extended quadrilateral for boundary patch `k`.
"""
function plotnsbd(dp::domprop, d::abstractdomain, k::Integer)

    fig, ax = draw(d, 1)

    # k is the actual boundary patch index.
    # kb is the boundary-local index such that d.kd[kb] == k.
    kb = findfirst(==(k), d.kd)

    kb === nothing && error("Patch k = $k is not a boundary patch. It must belong to d.kd.")

    # -- draw δ_bd-extended boundary quadrilateral --
    @views Pbd = d.Qptsbd[:, kb]

    Pext = Vector{Float64}(undef, 8)

    ix = [1, 3, 5, 7, 1]
    iy = [2, 4, 6, 8, 2]

    extendqua!(Pext, Pbd, dp.del_bd)
    lines!(ax, Pext[ix], Pext[iy], color=:black, linewidth=2)

    # -- collect near-singular interior targets to boundary patch k --
    xbd = Float64[]
    ybd = Float64[]

    @inbounds for colpbd in 1:size(dp.prepts_bd, 2)
        if dp.prepts_bd[2, colpbd] == k
            i = dp.prepts_bd[1, colpbd]

            push!(xbd, dp.tgtpts[1, i])
            push!(ybd, dp.tgtpts[2, i])
        end
    end

    !isempty(xbd) && scatter!(ax, xbd, ybd; markersize=9, color=:blue, strokewidth=0)

    return (fig, ax)
end


"""
Plot the projection of near-singular interior targets onto boundary patch `k`.
Here `k` is the actual boundary patch index, i.e. `k ∈ d.kd`.
For every interior target point near-singular to boundary patch `k`:
- the original target point is shown in black,
- the projected boundary point is shown in blue,
- an arrow is drawn from the target point to its projection.
"""
function plotprojbd(dp::domprop, d::abstractdomain, k::Integer)

    fig, ax = draw(d, 1)

    # k is the actual boundary patch index.
    # kb is the boundary-local index such that d.kd[kb] == k.
    kb = findfirst(==(k), d.kd)
    kb === nothing && error("Patch k = $k is not a boundary patch. It must belong to d.kd.")

    # -- draw δ_bd-extended boundary quadrilateral --
    @views Pbd = d.Qptsbd[:, kb]

    Pext = Vector{Float64}(undef, 8)

    ix = [1, 3, 5, 7, 1]
    iy = [2, 4, 6, 8, 2]

    extendqua!(Pext, Pbd, dp.del_bd)
    lines!(ax, Pext[ix], Pext[iy], color=:black, linewidth=2)

    # -- collect target points and projected points --
    xt = Float64[]
    yt = Float64[]
    xp = Float64[]
    yp = Float64[]

    @inbounds for colpbd in 1:size(dp.prepts_bd, 2)

        if dp.prepts_bd[2, colpbd] == k

            i = dp.prepts_bd[1, colpbd]
            τ = dp.projpts[colpbd]

            # original near-singular target
            x0 = dp.tgtpts[1, i]
            y0 = dp.tgtpts[2, i]

            # projected point on boundary patch k
            x1 = gamx(d, τ, k)
            y1 = gamy(d, τ, k)

            push!(xt, x0)
            push!(yt, y0)
            push!(xp, x1)
            push!(yp, y1)

        end
    end

    if !isempty(xt)
        arrows2d!(ax, xt, yt, xp .- xt, yp .- yt;
            tipwidth=8, tiplength=12, shaftwidth=1.5, color=:gray)
        scatter!(ax, xt, yt; markersize=9, color=:black,strokewidth=0)
        scatter!(ax, xp, yp; markersize=9, color=:blue, strokewidth=0)
    end

    return fig, ax
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
    println(io, "  N:       ", d.N)
    println(io, "  del:     ", d.del)
    println(io, "  del_bd:  ", d.del_bd)

    # quick counts
    maplen = length(d.hmap)
    maplenbd = length(d.hmap_bd)

    println(io)
    println(io, "  hmap     :  Dict  (", maplen,  " entries)")
    println(io, "  hmap_bd  :  Dict  (", maplenbd,  " entries)")
    println(io, "  tgtpts   :  Matrix{Float64} (", size(d.tgtpts,1), "×", size(d.tgtpts,2), ")")
    println(io, "  prepts   :  Matrix{Int}     (", size(d.prepts,1), "×", size(d.prepts,2), ")")
    println(io, "  prepts_bd:  Matrix{Int}     (", size(d.prepts_bd,1), "×", size(d.prepts_bd,2), ")")
    println(io, "  invpts   :  Matrix{Float64} (", size(d.invpts,1), "×", size(d.invpts,2), ")")
    println(io, "  projpts  :  Vector{Float64} (", length(d.projpts), " projections)")

    # bookkeeping vectors (show size + a short preview)
    println(io)
    println(io, "  pthgo:    ", _preview(d.pthgo))
    println(io, "  nspsz:    ", _preview(d.nspsz))

    return
end
