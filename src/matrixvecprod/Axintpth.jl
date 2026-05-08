function Axintpth!(v::SubArray{Float64}, kM::Int, IntS::Matrix{Float64},
    d::D, dp::domprop, s::Float64, IV::IVT) where {D<:abstractdomain, IVT}
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
                        KIr[i, j] = Wrg[i, j] / (r2^s)
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
    # --------- Part 2: boundary rows (if s < 0.5) ---------
    if s < 0.5

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
                        KIr[i, j] = Wrg[i, j] / (r2^s)
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

    return nothing
end

#Old way of computing regular patch integral
# ---------- Regular patch  ----------
# @. Ker = 1 / ((x1 - Zx)^2 + (x2 - Zy)^2)^s

# @inbounds for jj in 1:Np

#     q, r = divrem(jj - 1, N)
#     j1 = r + 1     # remainder
#     j2 = q + 1     # quotient

#     @views g1 = idctrg[:, j1]
#     @views g2 = idctrg[:, j2]

#     mul!(Ur, g1, g2')

#     @. KIr = Ker * Ur * DJ * Df

#     v[row, jj] = Cs * dot(fwr, KIr, fwr)
# end
# Optimized version of the old regular-patch loop.
#-------------------------------------------
# Old code, for each basis index jj, did:
#
#     g1 = idctrg[:, j1]
#     g2 = idctrg[:, j2]
#     Ur = g1 * g2'
#     KIr = Ker .* Ur .* DJ .* Df
#     v[row, jj] = Cs * dot(fwr, KIr, fwr)
#
# where jj corresponds to the tensor-product basis pair
#     jj = j1 + (j2 - 1) * N.
# Expanding the old formula:
#
# dot(fwr, KIr, fwr)= sum_i sum_j fwr[i] * KIr[i,j] * fwr[j]
#
# and since
#
#     KIr[i,j] = Ker[i,j] * Ur[i,j] * DJ[i,j] * Df[i,j]
#     Ur[i,j]  = g1[i] * g2[j]
#
# the old value was
#
# v[row,jj] = Cs * sum_i sum_j (
#   fwr[i] * Ker[i,j] * DJ[i,j] * Df[i,j] * fwr[j] * g1[i] * g2[j]).
#
# For a fixed target row, the only part that depends on the basis pair
# (j1,j2) is g1 = idctrg[:,j1] and g2 = idctrg[:,j2].
# Everything else can be grouped into one weighted kernel matrix:
#
#     W[i,j] = fwr[i] * Ker[i,j] * DJ[i,j] * Df[i,j] * fwr[j].
#
# Since
#     Ker[i,j] = 1 / ((x1 - Zx[i,j])^2 + (x2 - Zy[i,j])^2)^s,
# we build this weighted kernel as
#     KIr[i,j] = Wrg[i,j] / r2^s
# where
#     Wrg[i,j] = fwr[i] * DJ[i,j] * Df[i,j] * fwr[j]
#     r2       = (x1 - Zx[i,j])^2 + (x2 - Zy[i,j])^2.
# Thus after this loop, KIr stores the full weighted kernel matrix:
#
#     KIr = diag(fwr) * (Ker .* DJ .* Df) * diag(fwr).
#
# For one tensor-product basis pair, the old integral is therefore
#
#     Cs * idctrg[:,j1]' * KIr * idctrg[:,j2].
#
# Instead of computing this scalar separately for every pair (j1,j2),
# we compute all N^2 values at once:
#
#     Grg = idctrg' * KIr * idctrg.
#
# Then
#
#     Grg[j1,j2] = idctrg[:,j1]' * KIr * idctrg[:,j2],
#
# which is exactly the old value before multiplying by Cs.
# Julia matrices are stored column-major, so the linear index
#     jj = j1 + (j2 - 1) * N
# satisfies
#     Grg[jj] == Grg[j1,j2].
# Therefore
#     v[row,jj] = Cs * Grg[jj]
# gives the same result as the old loop, but avoids forming Ur and KIr
# separately for every basis pair. The expensive quadrature projection is
# now done with two BLAS matrix multiplications:
#     Grgtmp = idctrg' * KIr
#     Grg    = Grgtmp * idctrg
# This is mathematically equivalent to the old scalar loop, but much faster.