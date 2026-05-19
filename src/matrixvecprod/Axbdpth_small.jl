function Axbdpth_small!(v::SubArray{Float64}, kM::Int, IntS::Matrix{Float64},
    d::D, dp::domprop, s::Float64, IV::IVT) where {D<:abstractdomain, IVT}

    # --------- Unpack IV (NamedTuple) ---------
    (; IV1, IVr, IVbdth, IVbt1, IVbt2) = IV
    (; N, Np, Cs, M, Mbd, nbd) = IV1
    (; zx2, zy2, coeffs, coeffs_sum) = IVbdth
    (; Zx₂, Zy₂, DJ₂, KIbd) = IVbt1
    (; mfw, CT, B1bd, B2bd, Gbdtmp, Gbd, Wbd) = IVbt2

    (; nr, fwr) = IVr
    # v is M*Np + Mbd*N by Np matrix (preallocated given)
    # that is, size(v) == (M*Np + Mbd*N, Np)

    # --------- Compute regular patch (kM) ---------
    # Zx₂, Zy₂: images of Chebyshev grid (zx2,zy2) on patch kM
    # Determinant of map and distance function on regular patch
    mapxy_Dmap!(Zx₂, Zy₂, DJ₂, d, zx2, zy2, kM) # nr×nr Zx₂, Zy₂, DJ₂

    hc = d.pths[kM].ck1 - d.pths[kM].ck0

    Cs₂ = Cs * hc^s 

    @inbounds for j in 1:nr
        @inbounds for i in 1:nbd
            Wbd[i, j] = mfw[i] * DJ₂[i, j] * fwr[j]
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
            # dp.hmap is keyed by (row,kM) :: Tuple{Int,Int}
            Ikey = packkey(row, kM)
            col = get(dp.hmap, Ikey, 0)

            if col != 0
                # ---------- Near-singular patch case ----------
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
                x1 = dp.tgtpts[1, row]
                x2 = dp.tgtpts[2, row]

                # Build W = diag(mfw) * (Ker₂ .* DJ₂) * diag(fwr)
                @inbounds for jj in 1:nr
                    @inbounds for ii in 1:nbd
                        dx = x1 - Zx₂[ii, jj]
                        dy = x2 - Zy₂[ii, jj]
                        r2 = dx * dx + dy * dy
                        KIbd[ii, jj] = Wbd[ii, jj] * expm1(-s*log(r2))
                    end
                end
                # Gbd = B1bd' * W * B2bd
                mul!(Gbdtmp, B1bd', KIbd)  # N × nr
                mul!(Gbd, Gbdtmp, B2bd)    # N × N

                @inbounds for jj in 1:Np
                    v[row, jj] = Cs₂ * Gbd[jj]
                end
            end
        end
    end

    #kM is a boundary patch, Cheby coefficiets of kM patch
    #will be non-zero and will therefore play a role 
    #in boundary equations
    # --------- Part 2: boundary rows ---------
    Lₚₘ = Ni + 1
    Lₚₙ = Ni + Mbd * N


        @inbounds for row in Lₚₘ:Lₚₙ

            k₀ = cld(row - Ni, N)

            l, j = d.kd[k₀],  row - Ni - (k₀ - 1)*N

            if l == kM
                # ---------- Singular patch ----------
                col = dp.pthgo[M+1] + j - 1 + N * (k₀ - 1)

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
                    x1 = dp.tgtpts[1, row]
                    x2 = dp.tgtpts[2, row]

                    # Build W = diag(mfw) * (Ker₂ .* DJ₂) * diag(fwr)
                    @inbounds for jj in 1:nr
                        @inbounds for ii in 1:nbd
                            dx = x1 - Zx₂[ii, jj]
                            dy = x2 - Zy₂[ii, jj]
                            r2 = dx * dx + dy * dy
                            KIbd[ii, jj] = Wbd[ii, jj] * expm1(-s*log(r2))
                        end
                    end

                    # Gbd = B1bd' * W * B2bd
                    mul!(Gbdtmp, B1bd', KIbd)  # N × nr
                    mul!(Gbd, Gbdtmp, B2bd)    # N × N

                    @inbounds for jj in 1:Np
                        v[row, jj] = Cs₂ * Gbd[jj]
                    end
                end
            end
        end


    return nothing

end