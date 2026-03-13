function Axbdpth!(v::SubArray{Float64}, kM::Int, IntS::Matrix{Float64},
    d::abstractdomain, dp::domprop, s::Float64, IV)

    # --------- Unpack IV (NamedTuple) ---------
    (; IV1, IVr, IVbdth, IVt, IVbt1, IVbt2) = IV
    (; N, Np, Cs, M, Mbd) = IV1
    (; zx2, zy2, Tz1, Tz, coeffs) = IVbdth
    (; Zx₂, Zy₂, DJ₂, Ker₂, KIbd) = IVbt1
    (; Ubd, Ubd1, Ubd2, mfw, CT) = IVbt2
    fwr = IVr.fwr
    CN =  IVt.CN
    # v is M*Np + Mbd*N by Np matrix (preallocated given)
    # that is, size(v) == (M*Np + Mbd*N, Np)

    # --------- Compute regular patch (kM) ---------
    # Zx₂, Zy₂: images of Chebyshev grid (zx2,zy2) on patch kM
    # Determinant of map and distance function on regular patch
    mapxy_Dmap!(Zx₂, Zy₂, DJ₂, d, zx2, zy2, kM) # nr×nr Zx₂, Zy₂, DJ₂

    hc = d.pths[kM].ck1 - d.pths[kM].ck0

    Dhc = s>=0.5 ? hc^(s-1) : hc^s 

    # --------- interior rows  ---------
    Lₚ = M * Np

    @inbounds for row in 1:Lₚ

        x1 = dp.tgtpts[1, row]
        x2 = dp.tgtpts[2, row]

        l = cld(row, Np)
        j = row - (l - 1)*Np

        if l == kM
            # ---------- Singular patch ----------
            col = dp.pthgo[l] + j - 1

            @views SI = IntS[:, col]  # length Np vector

            @views IM = reshape(SI, N, N) # Integral matrix

            @inbounds for jj in 1:Np

                q, r = divrem(jj - 1, N)
                j1 = r + 1     # remainder
                j2 = q + 1     # quotient

                @views c1 = coeffs[:, j1]
                @views c2 = coeffs[:, j2]

                mul!(CN, c1, c2')

                v[row, jj] = Cs * dot(CN, IM)
            end

        else
            # --------- Check for near-singular patch ----------
            # dp.hmap is keyed by (row,kM) :: Tuple{Int,Int}
            if haskey(dp.hmap, packkey(row, kM))
                # ---------- Near-singular patch case ----------
                col = dp.hmap[packkey(row, kM)]

                @views NSI = IntS[:, col]

                @views IM = reshape(NSI, N, N)

                @inbounds for jj in 1:Np

                    q, r = divrem(jj - 1, N)
                    j1 = r + 1     # remainder
                    j2 = q + 1     # quotient

                    @views c1 = coeffs[:, j1]
                    @views c2 = coeffs[:, j2]

                    mul!(CN, c1, c2')

                    v[row, jj] = Cs * dot(CN, IM)
                end

            else
                # ---------- Regular patch  ----------
                @. Ker₂ = 1 / ((x1 - Zx₂)^2 + (x2 - Zy₂)^2)^s

                @inbounds for jj in 1:Np

                    q, r = divrem(jj - 1, N)
                    j1 = r + 1     # remainder
                    j2 = q + 1     # quotient

                    @views c1 = coeffs[:, j1]
                    @views c2 = coeffs[:, j2]

                    mul!(Ubd1, Tz1, c1)
                    mul!(Ubd2, Tz , c2)

                    mul!(Ubd, Ubd1, Ubd2')

                    @. KIbd = Ker₂ * Ubd * DJ₂ 

                    v[row, jj] = Cs * Dhc * dot(mfw, KIbd, fwr)
                end
            end
        end
    end

    #kM is a boundary patch, Cheby coefficiets of kM patch
    #will be non-zero and will therefore play a role 
    #in boundary equations
    # --------- Part 2: boundary rows ---------
    Lₚₘ = Lₚ + 1
    Lₚₙ = Lₚ + Mbd * N

    if s >= 0.5

        @inbounds for row in Lₚₘ:Lₚₙ

            k₀ = cld(row - M*Np, N)

            l, j = d.kd[k₀],  row - M*Np - (k₀ - 1)*N

            if l == kM

                @views CTj = CT[:, j]

                @inbounds for jj in 1:Np

                    q, r = divrem(jj - 1, N)
                    j1 = r + 1     # remainder
                    j2 = q + 1     # quotient

                    @views c1 = coeffs[:, j1]
                    @views c2 = coeffs[:, j2]

                    v[row, jj] = sum(c1) * dot(c2, CTj)

                end

            end
        end
    else

        @inbounds for row in Lₚₘ:Lₚₙ

            k₀ = cld(row - M*Np, N)

            l, j = d.kd[k₀],  row - M*Np - (k₀ - 1)*N

            if l == kM
                # ---------- Singular patch ----------
                col = dp.pthgo[M+1] + j - 1 + N*(ceil(Int, (row-M*Np)/N) - 1)

                @views SI = IntS[:, col]  # length Np vector

                @views IM = reshape(SI, N, N) # Integral matrix

                @inbounds for jj in 1:Np

                    q, r = divrem(jj - 1, N)
                    j1 = r + 1     # remainder
                    j2 = q + 1     # quotient

                    @views c1 = coeffs[:, j1]
                    @views c2 = coeffs[:, j2]

                    mul!(CN, c1, c2')

                    v[row, jj] = Cs * dot(CN, IM)
                end

            else
                if haskey(dp.hmap, packkey(row, kM))
                    # ---------- Near-singular patch case ----------
                    col = dp.hmap[packkey(row, kM)]

                    @views NSI = IntS[:, col]

                    @views IM = reshape(NSI, N, N)

                    @inbounds for jj in 1:Np

                        q, r = divrem(jj - 1, N)
                        j1 = r + 1     # remainder
                        j2 = q + 1     # quotient

                        @views c1 = coeffs[:, j1]
                        @views c2 = coeffs[:, j2]

                        mul!(CN, c1, c2')

                        v[row, jj] = Cs * dot(CN, IM)
                    end

                else
                    x1 = dp.tgtpts[1, row]
                    x2 = dp.tgtpts[2, row]
                    # ---------- Regular patch  ----------
                    @. Ker₂ = 1 / ((x1 - Zx₂)^2 + (x2 - Zy₂)^2)^s

                    @inbounds for jj in 1:Np

                        q, r = divrem(jj - 1, N)
                        j1 = r + 1     # remainder
                        j2 = q + 1     # quotient

                        @views c1 = coeffs[:, j1]
                        @views c2 = coeffs[:, j2]

                        mul!(Ubd1, Tz1, c1)
                        mul!(Ubd2, Tz, c2)

                        mul!(Ubd, Ubd1, Ubd2')

                        @. KIbd = Ker₂ * Ubd * DJ₂

                        v[row, jj] = Cs * Dhc * dot(mfw, KIbd, fwr)
                    end
                end
            end
        end

    end

    return nothing

end