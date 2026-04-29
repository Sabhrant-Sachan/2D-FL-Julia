"""
Boundary-domain properties for SLP boundary equations,
SLP stands for Single layer potential. 

This is a boundary-only analogue of `domprop`. Patches
are curves in 2D space, which were formed due to patching
structure in the interioir. 

Target boundary points are indexed globally as

    i = 1, ..., Mbd*N

where

    k₀ = cld(i, N)              # boundary-patch block index
    k  = kd[k₀]                 # actual domain patch index
    j  = i - (k₀ - 1)*N         # local Chebyshev index on that boundary patch

The target point is

    γ_k(cospi((2j - 1)/(2N))).

I store the real space coordinates in `tgtpts`. The patch/index data can
always be recovered from the column index `i` as shown above.

`prepts` has two rows:

    prepts[1, col] = boundary target index i
    prepts[2, col] = integration boundary patch k

For each boundary integration patch kd[k₀], the precomputation layout is

    N regular boundary targets on that patch,
    followed by nspsz[k₀] near-singular boundary targets.

`hmap` maps

    packkey(boundary target index i, integration patch k) -> column in prepts

for near-singular pairs.

`projpts` stores the projected parameter values for near-singular pairs only.
There is no `hmaproj`: if `col` is the column in `prepts` for a near-singular
point in boundary block `k₀`, then its column in `projpts` is

    ll = col - k₀*N

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

function SLPprecomps(d::abstractdomain, dpbd::dompropbd, p::Int; n::Int=128)::Matrix{Float64}

end

function SLPbeta(d::abstractdomain, N::Int)::Matrix{Float64}

Mbd = length(d.kd)

Nd = Mbd * N

δ, p = 0.1, 6

dpbd = dompropbd(N, δ, d)

IntS = SLPprecomps(d, dpbd, p)

end