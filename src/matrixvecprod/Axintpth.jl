function Axintpth!(v::SubArray{Float64}, kM::Int, IntS::Matrix{Float64},
    d::abstractdomain, dp::domprop, s::Float64, IV)

    # --------- Unpack IV (NamedTuple) ---------
    (; IV1, IVr, IVbdth, IVt) = IV
    (; N, Np, Cs, M, Mbd) = IV1
    (; fwr, zx, zt, zy, Df, idctrg) = IVr
    (; Zx, Zy, DJ, Ker, KIr, Ur, CN) = IVt
    coeffs = IVbdth.coeffs  # size N×N

    # v is M*Np + Mbd*N by Np matrix (preallocated given)
    # that is, size(v) == (M*Np + Mbd*N, Np)

    # --------- Compute regular patch (kM) ---------
    # Zx, Zy: images of Chebyshev grid (zx,zy) on patch kM
    # Determinant of map and distance function on regular patch
    mapxy_Dmap!(Zx, Zy, DJ, d, zx, zy, kM) # nr×nr Zx, Zy, DJ

    dfunc!(Df, d, kM, zt, s)  # zt = 1-zx

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
            # dp.hmap is keyed by (row,kM) :: Tuple{UInt,Int}
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
                @. Ker = 1 / ((x1 - Zx)^2 + (x2 - Zy)^2)^s

                @inbounds for jj in 1:Np

                    q, r = divrem(jj - 1, N)
                    j1 = r + 1     # remainder
                    j2 = q + 1     # quotient

                    @views g1 = idctrg[:, j1]
                    @views g2 = idctrg[:, j2]

                    mul!(Ur, g1, g2')

                    @. KIr = Ker * Ur * DJ * Df

                    v[row, jj] = Cs * dot(fwr, KIr, fwr)
                end
            end
        end
    end

    #Integrals involved for boundary points when s is small
    # --------- Part 2: boundary rows (if s < 0.5) ---------
    if s < 0.5

        Lₚₘ = Lₚ + 1
        Lₚₙ = Lₚ + Mbd * N

        @inbounds for row in Lₚₘ:Lₚₙ

            x1 = dp.tgtpts[1, row]
            x2 = dp.tgtpts[2, row]

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
                @. Ker = 1 / ((x1 - Zx)^2 + (x2 - Zy)^2)^s

                @inbounds for jj in 1:Np

                    q, r = divrem(jj - 1, N)
                    j1 = r + 1     # remainder
                    j2 = q + 1     # quotient

                    @views g1 = idctrg[:, j1]
                    @views g2 = idctrg[:, j2]

                    mul!(Ur, g1, g2')

                    @. KIr = Ker * Ur * DJ * Df

                    v[row, jj] = Cs * dot(fwr, KIr, fwr)
                end
            end
        end
    end

    return nothing
end
