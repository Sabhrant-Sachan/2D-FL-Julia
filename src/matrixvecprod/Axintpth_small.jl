#Function for small s case only
function Axintpth_small!(v::SubArray{Float64}, kM::Int, IntS::Matrix{Float64},
    d::D, dp::domprop, s::Float64, IV::IVT) where {D<:abstractdomain,IVT}
    # --------- Unpack IV (NamedTuple) ---------
    (; IV1, IVr, IVbdth, IVt) = IV
    (; N, Np, Cs, M, Mbd) = IV1
    (; nr, fwr, zx, zt, zy, Df, idctrg) = IVr
    (; Zx, Zy, DJ, KIr, Wrg, Grgtmp, Grg) = IVt
    (; coeffs) = IVbdth  # size N×N

    # v is M*Np + Mbd*N by Np matrix (preallocated given)
    # that is, size(v) == (M*Np + Mbd*N, Np)

    # --------- Compute regular patch (kM) ---------
    # Zx, Zy: images of Chebyshev grid (zx,zy) on patch kM
    # Determinant of map and distance function on regular patch
    mapxy_Dmap!(Zx, Zy, DJ, d, zx, zy, kM) # nr×nr Zx, Zy, DJ

    dfunc!(Df, d, kM, zt, s)  # zt = 1-zx

    @inbounds for j in 1:nr
        @inbounds for i in 1:nr
            Wrg[i, j] = fwr[i] * DJ[i, j] * Df[i, j] * fwr[j]
        end
    end

    # --------- interior rows  ---------
    Ni = M * Np

    @inbounds for row in 1:Ni

        x1 = dp.tgtpts[1, row]
        x2 = dp.tgtpts[2, row]

        l = cld(row, Np)
        j = row - (l - 1) * Np

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

                #mul!(CN, c1, c2')
                #v[row, jj] = Cs * dot(CN, IM)
                v[row, jj] = Cs * dot(c1, IM, c2)
            end

        else
            # --------- Check for near-singular patch ----------
            # dp.hmap is keyed by (row,kM) :: Tuple{UInt,Int}
            Ikey = packkey(row, kM)
            col = get(dp.hmap, Ikey, 0)

            if col != 0

                @views NSI = IntS[:, col]

                @views IM = reshape(NSI, N, N)

                @inbounds for jj in 1:Np

                    q, r = divrem(jj - 1, N)
                    j1 = r + 1     # remainder
                    j2 = q + 1     # quotient

                    @views c1 = coeffs[:, j1]
                    @views c2 = coeffs[:, j2]

                    #mul!(CN, c1, c2')
                    #v[row, jj] = Cs * dot(CN, IM)
                    v[row, jj] = Cs * dot(c1, IM, c2)
                end

            else
                @inbounds for j in 1:nr
                    @inbounds for i in 1:nr
                        dx = x1 - Zx[i, j]
                        dy = x2 - Zy[i, j]
                        r2 = dx * dx + dy * dy
                        KIr[i, j] = Wrg[i, j] * expm1(-s * log(r2))
                    end
                end

                mul!(Grgtmp, idctrg', KIr)  # N × nr
                mul!(Grg, Grgtmp, idctrg)   # N × N

                @inbounds for jj in 1:Np
                    v[row, jj] = Cs * Grg[jj]
                end
            end
        end
    end

    #Integrals involved for boundary points when s is small
    # --------- Part 2: boundary rows ---------

    Nb = Mbd * N
    Lₚₘ = Ni + 1
    Lₚₙ = Ni + Nb

    @inbounds for row in Lₚₘ:Lₚₙ

        x1 = dp.tgtpts[1, row]
        x2 = dp.tgtpts[2, row]

        Ikey = packkey(row, kM)
        col = get(dp.hmap, Ikey, 0)

        if col != 0

            @views NSI = IntS[:, col]

            @views IM = reshape(NSI, N, N)

            @inbounds for jj in 1:Np

                q, r = divrem(jj - 1, N)
                j1 = r + 1     # remainder
                j2 = q + 1     # quotient

                @views c1 = coeffs[:, j1]
                @views c2 = coeffs[:, j2]

                #mul!(CN, c1, c2')
                #v[row, jj] = Cs * dot(CN, IM)
                v[row, jj] = Cs * dot(c1, IM, c2)
            end

        else
            @inbounds for j in 1:nr
                @inbounds for i in 1:nr
                    dx = x1 - Zx[i, j]
                    dy = x2 - Zy[i, j]
                    r2 = dx * dx + dy * dy
                    KIr[i, j] = Wrg[i, j] * expm1(-s * log(r2))
                end
            end

            mul!(Grgtmp, idctrg', KIr)  # N × nr
            mul!(Grg, Grgtmp, idctrg)   # N × N

            @inbounds for jj in 1:Np
                v[row, jj] = Cs * Grg[jj]
            end
        end
    end


    return nothing
end