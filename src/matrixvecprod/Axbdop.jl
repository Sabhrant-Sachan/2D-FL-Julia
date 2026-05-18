function Axbdop!(v::SubArray{Float64}, kI::Int, d::D, dp::domprop, 
    s::Float64, IV::IVT) where {D<:abstractdomain, IVT}
    # --------- Unpack IV (NamedTuple) ---------
    (; IV1, IVbd, IVbdt1, IVbdt2, IVbdt3, IVbdth, IVbt2) = IV

    (; N, Np, M, Mbd) = IV1

    (; nr_bdy, ns_near, fwr_bdy, fw_near, zr_bdy, idctbd,
       inv2π, BD_FAR, BD_NEAR, BD_INTP) = IVbd

    (; y₁, y₂, yt1, yt2,
       wmz₁, wmz₂, dwz₁, dwz₂,
       gamkbdy, gampkbdy,
       gamks₁, gamks₂,
       gamperpks₁, gamperpks₂,
       gamkt1, gamkt2,
       gamperpkt1, gamperpkt2) = IVbdt1

    (; kbd_bdy, kbd₁, kbd₂, μ₀, γt1, γt2) = IVbdt2

    (; fvals, FVals, γxbd, Tbd, Tbdt₀,
       Bbdt1, Bbdt2, Bbdt₀) = IVbdt3

    (; coeffs) = IVbdth

    (; CT) = IVbt2

    k = d.kd[kI]

    # v is M*Np + Mbd*N by N matrix (preallocated given)
    # that is, size(v) == (M*Np + Mbd*N, N)

    # --------- Compute boundary curve points ---------
    gam!(gamkbdy, d, zr_bdy, k)
    gamp!(gampkbdy, d, zr_bdy, k)

    gam!(gamkt1, d, yt1, k)
    gam!(gamkt2, d, yt2, k)
    gamp!(gamperpkt1, d, yt1, k)
    gamp!(gamperpkt2, d, yt2, k)

    # --------- interior rows  ---------
    Ni = M * Np

    jnear = 0 

    #Contribution of all points for the kth patch
    @inbounds for row in 1:Ni

        x1 = dp.tgtpts[1, row]
        x2 = dp.tgtpts[2, row]

        if dp.bdmode[row] == BD_FAR
            # Mode-0 point, far from all panels, direct computation
            @inbounds for j in 1:nr_bdy
                dx1 = gamkbdy[1, j] - x1
                dx2 = gamkbdy[2, j] - x2

                num = dx1 * gampkbdy[1, j] + dx2 * gampkbdy[2, j]
                den = dx1 * dx1 + dx2 * dx2

                kbd_bdy[j] = (num * fwr_bdy[j] * inv2π) / den
            end

            rview = @view v[row, :]

            mul!(rview, transpose(idctbd), kbd_bdy)

        elseif dp.bdmode[row] == BD_NEAR
            # Mode-1 point, near to some panels, but not close enough
            jnear += 1

            q1 = dp.bdnearptr[jnear]
            q2 = dp.bdnearptr[jnear+1] - 1

            near_flag = false
            @inbounds for q in q1:q2
                if dp.bdneark[q] == k
                    t₀ = dp.bdneart[q]
                    near_flag = true
                    if t₀ == 1.0

                        @inbounds for j in 1:ns_near
                            dx1 = gamkt1[1, j] - x1
                            dx2 = gamkt1[2, j] - x2

                            num = dx1 * gamperpkt1[1, j] + dx2 * gamperpkt1[2, j]
                            den = dx1 * dx1 + dx2 * dx2

                            kbd₁[j] = (num * fw_near[j] * dwz₁[j] * inv2π) / den
                        end

                        rview = @view v[row, :]

                        mul!(rview, transpose(Bbdt1), kbd₁)

                    elseif t₀ == -1.0
                        # So only the y₂-side contributes.

                        @inbounds for j in 1:ns_near
                            dx1 = gamkt2[1, j] - x1
                            dx2 = gamkt2[2, j] - x2

                            num = dx1 * gamperpkt2[1, j] + dx2 * gamperpkt2[2, j]
                            den = dx1 * dx1 + dx2 * dx2

                            kbd₂[j] = (num * fw_near[j] * dwz₂[j] * inv2π) / den
                        end

                        rview = @view v[row, :]

                        mul!(rview, transpose(Bbdt2), kbd₂)

                    else
                        # Both sides contribute:
                        @. y₁ = t₀ - (t₀ + 1.0) * wmz₁
                        @. y₂ = t₀ + (1.0 - t₀) * wmz₂

                        # First side: y₁.
                        gam!(gamks₁, d, y₁, k)
                        gamp!(gamperpks₁, d, y₁, k)

                        ChebyTN!(Tbdt₀, N, y₁)

                        mul!(Bbdt₀, Tbdt₀, coeffs)

                        @inbounds for j in 1:ns_near
                            dx1 = gamks₁[1, j] - x1
                            dx2 = gamks₁[2, j] - x2

                            num = dx1 * gamperpks₁[1, j] + dx2 * gamperpks₁[2, j]
                            den = dx1 * dx1 + dx2 * dx2

                            kbd₁[j] = (num * fw_near[j] * dwz₁[j] * inv2π) / den
                        end

                        rview = @view v[row, :]

                        mul!(rview, transpose(Bbdt₀), kbd₁, (t₀ + 1.0) / 2.0, 0.0)

                        # Second side: y₂.
                        gam!(gamks₂, d, y₂, k)
                        gamp!(gamperpks₂, d, y₂, k)

                        ChebyTN!(Tbdt₀, N, y₂)

                        mul!(Bbdt₀, Tbdt₀, coeffs)

                        @inbounds for j in 1:ns_near
                            dx1 = gamks₂[1, j] - x1
                            dx2 = gamks₂[2, j] - x2

                            num = dx1 * gamperpks₂[1, j] + dx2 * gamperpks₂[2, j]
                            den = dx1 * dx1 + dx2 * dx2

                            kbd₂[j] = (num * fw_near[j] * dwz₂[j] * inv2π) / den
                        end

                        mul!(rview, transpose(Bbdt₀), kbd₂, (1.0 - t₀) / 2.0, 1.0)

                    end

                    break
                end
            end

            if !near_flag
                @inbounds for j in 1:nr_bdy
                    dx1 = gamkbdy[1, j] - x1
                    dx2 = gamkbdy[2, j] - x2

                    num = dx1 * gampkbdy[1, j] + dx2 * gampkbdy[2, j]
                    den = dx1 * dx1 + dx2 * dx2

                    kbd_bdy[j] = (num * fwr_bdy[j] * inv2π) / den
                end

                rview = @view v[row, :]

                mul!(rview, transpose(idctbd), kbd_bdy)
            end

        end

    end

    jintp = 0

    @inbounds for row in 1:Ni
        dp.bdmode[row] != BD_INTP && continue

        jintp += 1

        l = dp.bdclosest[jintp]
        tₛ = dp.bdt[jintp]
        ϵ = dp.bddist[jintp]

        gam!(γxbd, d, tₛ, l)
        nu!(γt1, d, tₛ, l)

        # ------------------------------------------------------------
        # m = 1: boundary interpolation value.
        # ------------------------------------------------------------
        ChebyTN!(Tbd, N, tₛ)

        fv = @view FVals[:, 1]
        fill!(fv, 0.0)

        if k == l
            # Jump term for this source block only.
            mul!(fv, transpose(coeffs), Tbd, 0.5, 0.0)
        end

        a = (jintp - 1) * dp.Lᵢₙ + 1

        q1 = dp.bdintpptr[a]
        q2 = dp.bdintpptr[a+1] - 1

        if k == l
            DLP!(kbd_bdy, d, tₛ, zr_bdy, k, μ₀, γt1, γt2)

            @. kbd_bdy = kbd_bdy * fwr_bdy * inv2π

            mul!(fv, transpose(idctbd), kbd_bdy, 1.0, 1.0)

        else
            near_flag = false
            s_inv = 0.0

            for q in q1:q2
                if dp.bdintpk[q] == k
                    near_flag = true
                    s_inv = dp.bdintpt[q]
                    break
                end
            end

            if near_flag
                DLP!(kbd_bdy, d, s_inv, zr_bdy, k, μ₀, γt1, γt2)

                @. kbd_bdy = kbd_bdy * fwr_bdy * inv2π

                mul!(fv, transpose(idctbd), kbd_bdy, 1.0, 1.0)

            else
                for j in 1:nr_bdy
                    dx1 = gamkbdy[1, j] - γxbd[1]
                    dx2 = gamkbdy[2, j] - γxbd[2]

                    num = dx1 * gampkbdy[1, j] + dx2 * gampkbdy[2, j]
                    den = dx1 * dx1 + dx2 * dx2

                    kbd_bdy[j] = (num * fwr_bdy[j] * inv2π) / den
                end

                mul!(fv, transpose(idctbd), kbd_bdy, 1.0, 1.0)
            end
        end

        # ------------------------------------------------------------
        # m = 2:Lᵢₙ: auxiliary interior interpolation nodes.
        # ------------------------------------------------------------
        @inbounds for m in 2:dp.Lᵢₙ

            δ = dp.bdxvals[m]

            μ₀[1] = γxbd[1] - δ * γt1[1]
            μ₀[2] = γxbd[2] - δ * γt1[2]

            @inbounds for j in 1:N
                FVals[j, m] = 0.0
            end

            a = (jintp - 1) * dp.Lᵢₙ + m

            q1 = dp.bdintpptr[a]
            q2 = dp.bdintpptr[a+1] - 1

            near_flag = false
            tnear = 0.0

            for q in q1:q2
                if dp.bdintpk[q] == k
                    near_flag = true
                    tnear = dp.bdintpt[q]
                    break
                end
            end

            if near_flag

                t₀ = tnear

                if t₀ == 1.0

                    @inbounds for j in 1:ns_near
                        dx1 = gamkt1[1, j] - μ₀[1]
                        dx2 = gamkt1[2, j] - μ₀[2]

                        num = dx1 * gamperpkt1[1, j] + dx2 * gamperpkt1[2, j]
                        den = dx1 * dx1 + dx2 * dx2

                        kbd₁[j] = (num * fw_near[j] * dwz₁[j] * inv2π) / den
                    end

                    fv = @view FVals[:, m]

                    mul!(fv, transpose(Bbdt1), kbd₁, 1.0, 0.0)

                elseif t₀ == -1.0

                    @inbounds for j in 1:ns_near
                        dx1 = gamkt2[1, j] - μ₀[1]
                        dx2 = gamkt2[2, j] - μ₀[2]

                        num = dx1 * gamperpkt2[1, j] + dx2 * gamperpkt2[2, j]
                        den = dx1 * dx1 + dx2 * dx2

                        kbd₂[j] = (num * fw_near[j] * dwz₂[j] * inv2π) / den
                    end

                    fv = @view FVals[:, m]

                    mul!(fv, transpose(Bbdt2), kbd₂, 1.0, 0.0)

                else
                    # Both sides contribute:
                    @. y₁ = t₀ - (t₀ + 1.0) * wmz₁
                    @. y₂ = t₀ + (1.0 - t₀) * wmz₂

                    # First side: y₁.
                    gam!(gamks₁, d, y₁, k)
                    gamp!(gamperpks₁, d, y₁, k)

                    ChebyTN!(Tbdt₀, N, y₁)

                    mul!(Bbdt₀, Tbdt₀, coeffs)

                    @inbounds for j in 1:ns_near
                        dx1 = gamks₁[1, j] - μ₀[1]
                        dx2 = gamks₁[2, j] - μ₀[2]

                        num = dx1 * gamperpks₁[1, j] + dx2 * gamperpks₁[2, j]
                        den = dx1 * dx1 + dx2 * dx2

                        kbd₁[j] = (num * fw_near[j] * dwz₁[j] * inv2π) / den
                    end

                    fv = @view FVals[:, m]

                    mul!(fv, transpose(Bbdt₀), kbd₁, (t₀ + 1.0) / 2.0, 0.0)

                    # Second side: y₂.
                    gam!(gamks₂, d, y₂, k)
                    gamp!(gamperpks₂, d, y₂, k)

                    ChebyTN!(Tbdt₀, N, y₂)

                    mul!(Bbdt₀, Tbdt₀, coeffs)

                    @inbounds for j in 1:ns_near
                        dx1 = gamks₂[1, j] - μ₀[1]
                        dx2 = gamks₂[2, j] - μ₀[2]

                        num = dx1 * gamperpks₂[1, j] + dx2 * gamperpks₂[2, j]
                        den = dx1 * dx1 + dx2 * dx2

                        kbd₂[j] = (num * fw_near[j] * dwz₂[j] * inv2π) / den
                    end

                    mul!(fv, transpose(Bbdt₀), kbd₂, (1.0 - t₀) / 2.0, 1.0)
                end

            else

                for j in 1:nr_bdy
                    dx1 = gamkbdy[1, j] - μ₀[1]
                    dx2 = gamkbdy[2, j] - μ₀[2]

                    num = dx1 * gampkbdy[1, j] + dx2 * gampkbdy[2, j]
                    den = dx1 * dx1 + dx2 * dx2

                    kbd_bdy[j] = (num * fwr_bdy[j] * inv2π) / den
                end

                fv = @view FVals[:, m]

                mul!(fv, transpose(idctbd), kbd_bdy, 1.0, 0.0)
            end
        end

        @inbounds for j in 1:N
            @inbounds for i in 1:dp.Lᵢₙ
                fvals[i] =  FVals[j, i]
            end
            v[row, j] = nevill!(dp.bdxvals, fvals, ϵ)
        end
    end

    if s < 0.5

        Nb = Mbd * N
        Lₚₘ = Ni + 1
        Lₚₙ = Ni + Nb

        for row in Lₚₘ:Lₚₙ

            k₀ = cld(row - Ni, N)

            #Boundary patch number
            ptl = d.kd[k₀]

            #Linear index of point on the boundary patch
            ptj = row - Ni - (k₀ - 1) * N

            rview = @view v[row, :]

            gam!(γxbd, d, CT[2, ptj], ptl)

            if k == ptl
                # Self-panel boundary principal value.
                DLP!(kbd_bdy, d, CT[2, ptj], zr_bdy, k, μ₀, γt1, γt2)

                @. kbd_bdy = kbd_bdy * fwr_bdy * inv2π

                mul!(rview, transpose(idctbd), kbd_bdy)

            else
                near_flag = false
                s_inv = 0.0

                gam!(γt1, d, 1.0, k)
                gam!(γt2, d, -1.0, k)

                if (hypot(γt1[1] - γxbd[1], γt1[2] - γxbd[2]) < dp.del_intp ||
                    hypot(γt2[1] - γxbd[1], γt2[2] - γxbd[2]) < dp.del_intp)
                    #This should rarely occur 
                    near_flag = true
                    s_inv = bdinv(d, CT[2, ptj], ptl, k)
                end

                if near_flag
                    # Near off-diagonal boundary panel.
                    # Use the extended inverse parameter on patch k.
                    DLP!(kbd_bdy, d, s_inv, zr_bdy, k, μ₀, γt1, γt2)

                    @. kbd_bdy = kbd_bdy * fwr_bdy * inv2π

                    mul!(rview, transpose(idctbd), kbd_bdy)

                else
                    # Far off-diagonal panel: ordinary direct quadrature.

                    @inbounds for j in 1:nr_bdy
                        dx1 = gamkbdy[1, j] - γxbd[1]
                        dx2 = gamkbdy[2, j] - γxbd[2]

                        num = dx1 * gampkbdy[1, j] + dx2 * gampkbdy[2, j]
                        den = dx1 * dx1 + dx2 * dx2

                        kbd_bdy[j] = (num * fwr_bdy[j] * inv2π) / den
                    end

                    mul!(rview, transpose(idctbd), kbd_bdy)
                end
            end

        end
        #The limiting value from the interioir is 1/2,
        for j in 1:N
            v[Ni+N*(kI-1)+j, j] += 0.5
        end
    end

end