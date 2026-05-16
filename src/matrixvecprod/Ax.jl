function Ax!(v::AbstractVector{Float64}, u::AbstractVector{Float64}, IntS::Matrix{Float64},
    d::D, dp::domprop, s::Float64, IV::IVT) where {D<:abstractdomain,IVT}
    # Ax! :Matrix-vector product v = A*u (matrix-free), for GMRES and 
    #      other iterative solvers. This subroutine computes the Matrix vector 
    #      product Au without constructing the matrix A.
    #      This function is given as a function handle to GMRES
    #      to solve the linear system Au=B without storing the matrix A.

    # INPUTS :
    #     u - Input vector u of size M*Np+Mbd*N for which we compute Au
    #         We will not mutate u. Mutation in vector v is fine. 
    #     v - Input vector v of size M*Np+Mbd*N which is v=Au 
    #  IntS - Precomputation matrix of size Np*(M*Np+Dp.Tnsp+Mbd*N) given using
    #         "precomps.jl" subroutine. Np is number of points per patch.
    #     d - An instance of the domain d
    #    dp - This is an instance of the class domprop. For more
    #          information on dp, see comments in "domprop.jl" file.
    #     s - A parameter in (0,1) for the fractional laplacian operator.
    #    IV - A structure which contains packaged variables.

    # OUTPUTS: None
    #     
    # Author:
    #      Sabhrant Sachan
    #      Email : ssachan@caletch.edu

    # --------- Unpack IV ---------
    (; IV1, IVr, IVbdth, IVbd, IVt, IVbt1, IVbt2,
        IVbdt1, IVbdt2, IVbdt3, IVAf, IVAdct) = IV

    (; p_dct2_dim2, p_dct2_dim1,
        p_dct3_dim2, p_dct3_dim1,
        p_dct2_N, p_dct3_bdy) = IVAdct

    (; chebcoef, ufin,
        ζ_bdy, ζ_neart1, ζ_neart2,
        UV, UFV, ζv, ζfv_bdy,
        CNnr, TzT) = IVAf

    (; N, Np, Cs, M, Mbd) = IV1

    (; nr, fwr, nrp, zx, zt, zy, Df) = IVr

    (; zx2, zy2, Tz1) = IVbdth

    (; nr_bdy, ns_near, fwr_bdy, fw_near, zr_bdy,
        inv2π, BD_FAR, BD_NEAR, BD_INTP) = IVbd

    (; Zx, Zy, DJ, Ker, KIr, Ur, CN) = IVt

    (; Zx₂, Zy₂, DJ₂, Ker₂, KIbd) = IVbt1

    (; Ubd, mfw, CT) = IVbt2

    (; y₁, y₂, yt1, yt2,
        wmz₁, wmz₂, dwz₁, dwz₂,
        gamkbdy, gampkbdy,
        gamks₁, gamks₂,
        gamperpks₁, gamperpks₂,
        gamkt1, gamkt2,
        gamperpkt1, gamperpkt2) = IVbdt1

    (; kbd_bdy, kbd₁, kbd₂, μ₀, γt1, γt2) = IVbdt2

    (; fvals, γxbd, Tbd, Tbdt₀, Tbdt1, Tbdt2) = IVbdt3

    Ni = M * Np
    Nb = Mbd * N

    # ----- STEP 1: Get chebyshev coefficients and finer function values -----
    # ----- Compute Chebyshev coeffs for interior patches + refined ufin -----
    # Fill chebcoef for interior patches:
    @inbounds for k in 1:M
        ak₁ = (k - 1) * Np
        ak₂ = (k - 1) * nrp
        # load function values UV from u (fill UV directly)
        # UV is a matrix of size N * N of type Float64
        # capital letters matrix, lowercase vectors
        # UV, UFV (u values and u finer values). 
        @inbounds for i in 1:Np
            UV[i] = u[ak₁+i]
        end

        # DCT-II along dim2, then dim1
        mul!(UV, p_dct2_dim2, UV)

        UV .*= 1 / N

        # UV[:, 1] ./= 2
        @inbounds for i in 1:N
            UV[i, 1] /= 2
        end

        mul!(UV, p_dct2_dim1, UV)

        UV .*= 1 / N

        # UV[1, :] ./= 2
        @inbounds for i in 1:N
            UV[1, i] /= 2
        end

        # store into chebcoef
        @inbounds for i in 1:Np
            chebcoef[ak₁+i] = UV[i]
        end

        #UFV is a nr * nr zero float64 matrix 
        UFV .= 0
        # UFV[1:N, 1:N] .= UV
        @inbounds for j in 1:N
            @inbounds for i in 1:N
                UFV[i, j] = UV[i, j]
            end
        end

        # Undo the first-mode scaling for dim = 2
        #UFV[:, 1] .*= 2
        @inbounds for i in 1:nr
            UFV[i, 1] *= 2
        end

        # Invert along dim = 2
        # UFV = 0.5 .* FFTW.r2r(UFV, FFTW.REDFT01, 2)
        mul!(UFV, p_dct3_dim2, UFV)

        UFV .*= 0.5

        # Undo the first-mode scaling for dim = 1
        # UFV[1, :] .*= 2
        @inbounds for i in 1:nr
            UFV[1, i] *= 2
        end

        # Invert along dim = 1
        # UFV = 0.5 .* FFTW.r2r(UFV, FFTW.REDFT01, 1)
        mul!(UFV, p_dct3_dim1, UFV)

        UFV .*= 0.5

        # store into chebcoef
        @inbounds for i in 1:nrp
            ufin[ak₂+i] = UFV[i]
        end

    end

    # ----- Boundary ζ: coefficients + values on nr_bdy grid -----
    # For k=1:Mbd boundary patches, u has ζ values on N Cheby nodes for each bd patch.
    # We will also compute and store ζ_neart1 and ζ_neart2 once over all bd patches.
    # ----- Boundary ζ: coefficients + values on fixed boundary grids -----
    @inbounds for kb in 1:Mbd
        ak = Ni + (kb - 1) * N
        ibd = (kb - 1) * nr_bdy
        inear = (kb - 1) * ns_near

        # Load ζ values on N boundary Chebyshev nodes.
        for i in 1:N
            ζv[i] = u[ak+i]
        end

        # Convert to Chebyshev coefficients.
        mul!(ζv, p_dct2_N, ζv)

        ζv .*= 1 / N
        ζv[1] /= 2

        # Store coefficients in the global coefficient vector.
        for i in 1:N
            chebcoef[ak+i] = ζv[i]
        end

        # Fixed endpoint-smoothed grids yt1/yt2:
        # these are not Chebyshev grids, so use Chebyshev interpolation.
        mul!(y₁, Tbdt1, ζv)
        mul!(y₂, Tbdt2, ζv)

        for i in 1:ns_near
            ζ_neart1[inear+i] = y₁[i]
            ζ_neart2[inear+i] = y₂[i]
        end

        # Regular boundary quadrature grid zr_bdy:
        # this is a Chebyshev grid of size nr_bdy, so use DCT-III.
        ζv[1] *= 2

        fill!(ζfv_bdy, 0.0)

        for i in 1:N
            ζfv_bdy[i] = ζv[i]
        end

        mul!(ζfv_bdy, p_dct3_bdy, ζfv_bdy)

        ζfv_bdy .*= 0.5

        for i in 1:nr_bdy
            ζ_bdy[ibd+i] = ζfv_bdy[i]
        end
    end

    # ----- STEP 2: Singular Integrals over all rows and DLP contributions -----
    @inbounds for row in 1:Ni

        ℓ = cld(row, Np)
        j = row - (ℓ - 1) * Np

        col = dp.pthgo[ℓ] + j - 1

        @views cf = chebcoef[((ℓ-1)*Np+1):(ℓ*Np)]

        @views SI = IntS[:, col]

        v[row] = Cs * dot(SI, cf)
    end

    if s < 0.5

        Lₚₘ = Ni + 1
        Lₚₙ = Ni + Nb

        @inbounds for row in Lₚₘ:Lₚₙ

            k₀ = cld(row - Ni, N)

            ℓ = d.kd[k₀]

            j = row - M * Np - (k₀ - 1) * N

            #---- Add singular contributions first ----
            col = dp.pthgo[M+1] + j - 1 + N * (k₀ - 1)

            @views cf = chebcoef[((ℓ-1)*Np+1):(ℓ*Np)]

            @views SI = IntS[:, col]

            v[row] = Cs * dot(SI, cf)

        end

    end

    #If s<0.5, all values of vector v are now initlized.
    #If s>=0.5, all rows from 1:Ni of v are now initlized.

    # ----- Contribution from DLP starts here -----

    for ll in 1:Mbd
        a_ll = Ni + (ll - 1) * N

        k = d.kd[ll]

        jnear = 0

        gam!(gamkbdy, d, zr_bdy, k)
        gamp!(gampkbdy, d, zr_bdy, k)

        gam!(gamkt1, d, yt1, k)
        gam!(gamkt2, d, yt2, k)
        gamp!(gamperpkt1, d, yt1, k)
        gamp!(gamperpkt2, d, yt2, k)

        #Contribution of all points for the ℓth patch
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

                @views ζv = ζ_bdy[(1+(ll-1)*nr_bdy):(ll*nr_bdy)]

                v[row] += dot(ζv, kbd_bdy)

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

                            @views ζv = ζ_neart1[(1+(ll-1)*ns_near):(ll*ns_near)]

                            v[row] = v[row] + dot(ζv, kbd₁)

                        elseif t₀ == -1.0
                            # So only the y₂-side contributes.

                            @inbounds for j in 1:ns_near
                                dx1 = gamkt2[1, j] - x1
                                dx2 = gamkt2[2, j] - x2

                                num = dx1 * gamperpkt2[1, j] + dx2 * gamperpkt2[2, j]
                                den = dx1 * dx1 + dx2 * dx2

                                kbd₂[j] = (num * fw_near[j] * dwz₂[j] * inv2π) / den
                            end

                            @views ζv = ζ_neart2[(1+(ll-1)*ns_near):(ll*ns_near)]

                            v[row] = v[row] + dot(ζv, kbd₂)
                        else
                            # Both sides contribute:
                            @. y₁ = t₀ - (t₀ + 1.0) * wmz₁
                            @. y₂ = t₀ + (1.0 - t₀) * wmz₂

                            # First side: y₁.
                            gam!(gamks₁, d, y₁, k)
                            gamp!(gamperpks₁, d, y₁, k)

                            ChebyTN!(Tbdt₀, N, y₁)

                            @views ζv = chebcoef[a_ll+1:a_ll+N]
                            #Re-using y₁ here, which is fine
                            mul!(y₁, Tbdt₀, ζv)

                            @inbounds for j in 1:ns_near
                                dx1 = gamks₁[1, j] - x1
                                dx2 = gamks₁[2, j] - x2

                                num = dx1 * gamperpks₁[1, j] + dx2 * gamperpks₁[2, j]
                                den = dx1 * dx1 + dx2 * dx2

                                kbd₁[j] = (num * fw_near[j] * dwz₁[j] * inv2π) / den
                            end

                            # Second side: y₂.
                            gam!(gamks₂, d, y₂, k)
                            gamp!(gamperpks₂, d, y₂, k)

                            ChebyTN!(Tbdt₀, N, y₂)
                            #Re-using y₂ here, which is fine
                            mul!(y₂, Tbdt₀, ζv)

                            @inbounds for j in 1:ns_near
                                dx1 = gamks₂[1, j] - x1
                                dx2 = gamks₂[2, j] - x2

                                num = dx1 * gamperpks₂[1, j] + dx2 * gamperpks₂[2, j]
                                den = dx1 * dx1 + dx2 * dx2

                                kbd₂[j] = (num * fw_near[j] * dwz₂[j] * inv2π) / den
                            end

                            v[row] = v[row] + ((t₀ + 1.0) * dot(kbd₁, y₁)
                                               +
                                               (1.0 - t₀) * dot(kbd₂, y₂)) / 2.0

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

                    @views ζv = ζ_bdy[(1+(ll-1)*nr_bdy):(ll*nr_bdy)]

                    v[row] = v[row] + dot(ζv, kbd_bdy)
                end

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

        # Find boundary-local index of closest patch l.
        ll_l = findfirst(==(l), d.kd)
        ll_l === nothing && error("closest boundary patch l=$l not found in d.kd")

        a_l = Ni + (ll_l - 1) * N

        # ------------------------------------------------------------
        # m = 1: boundary interpolation value.
        # ------------------------------------------------------------
        @views Cbd = chebcoef[(a_l+1):(a_l+N)]

        ChebyTN!(Tbd, N, tₛ)

        ζₛ = dot(Tbd, Cbd)

        fvals[1] = 0.5 * ζₛ

        a = (jintp - 1) * dp.Lᵢₙ + 1

        q1 = dp.bdintpptr[a]
        q2 = dp.bdintpptr[a+1] - 1

        for ll in 1:Mbd
            k = d.kd[ll]

            @views ζb = ζ_bdy[(1+(ll-1)*nr_bdy):(ll*nr_bdy)]

            if k == l
                DLP!(kbd_bdy, d, tₛ, zr_bdy, k, μ₀, γt1, γt2)

                @. kbd_bdy = kbd_bdy * fwr_bdy * inv2π

                fvals[1] += dot(ζb, kbd_bdy)

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
                    # Near off-diagonal boundary panel.
                    # bdintpt stores extended inverse s on source patch k.
                    DLP!(kbd_bdy, d, s_inv, zr_bdy, k, μ₀, γt1, γt2)

                    @. kbd_bdy = kbd_bdy * fwr_bdy * inv2π

                    fvals[1] += dot(ζb, kbd_bdy)

                else
                    gam!(gamkbdy, d, zr_bdy, k)
                    gamp!(gampkbdy, d, zr_bdy, k)

                    for j in 1:nr_bdy
                        dx1 = gamkbdy[1, j] - γxbd[1]
                        dx2 = gamkbdy[2, j] - γxbd[2]

                        num = dx1 * gampkbdy[1, j] + dx2 * gampkbdy[2, j]
                        den = dx1 * dx1 + dx2 * dx2

                        kbd_bdy[j] = (num * fwr_bdy[j] * inv2π) / den
                    end

                    fvals[1] += dot(ζb, kbd_bdy)
                end
            end
        end

        # ------------------------------------------------------------
        # m = 2:Lᵢₙ: auxiliary interior interpolation nodes.
        # ------------------------------------------------------------
        nu!(γt1, d, tₛ, l)

        for m in 2:dp.Lᵢₙ

            δ = dp.bdxvals[m]

            μ₀[1] = γxbd[1] - δ * γt1[1]
            μ₀[2] = γxbd[2] - δ * γt1[2]

            fvals[m] = 0.0

            a = (jintp - 1) * dp.Lᵢₙ + m

            q1 = dp.bdintpptr[a]
            q2 = dp.bdintpptr[a+1] - 1

            for ll in 1:Mbd
                k = d.kd[ll]

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
                        gam!(gamks₁, d, yt1, k)
                        gamp!(gamperpks₁, d, yt1, k)

                        for j in 1:ns_near
                            dx1 = gamks₁[1, j] - μ₀[1]
                            dx2 = gamks₁[2, j] - μ₀[2]

                            num = dx1 * gamperpks₁[1, j] + dx2 * gamperpks₁[2, j]
                            den = dx1 * dx1 + dx2 * dx2

                            kbd₁[j] = (num * fw_near[j] * dwz₁[j] * inv2π) / den
                        end

                        @views ζn = ζ_neart1[(1+(ll-1)*ns_near):(ll*ns_near)]

                        fvals[m] += dot(ζn, kbd₁)

                    elseif t₀ == -1.0
                        gam!(gamks₂, d, yt2, k)
                        gamp!(gamperpks₂, d, yt2, k)

                        for j in 1:ns_near
                            dx1 = gamks₂[1, j] - μ₀[1]
                            dx2 = gamks₂[2, j] - μ₀[2]

                            num = dx1 * gamperpks₂[1, j] + dx2 * gamperpks₂[2, j]
                            den = dx1 * dx1 + dx2 * dx2

                            kbd₂[j] = (num * fw_near[j] * dwz₂[j] * inv2π) / den
                        end

                        @views ζn = ζ_neart2[(1+(ll-1)*ns_near):(ll*ns_near)]

                        fvals[m] += dot(ζn, kbd₂)

                    else
                        # Both sides contribute.
                        @. y₁ = t₀ - (t₀ + 1.0) * wmz₁
                        @. y₂ = t₀ + (1.0 - t₀) * wmz₂

                        # Source-panel boundary coefficients.
                        a_k = Ni + (ll - 1) * N
                        @views ζcoefk = chebcoef[(a_k+1):(a_k+N)]

                        # First side.
                        gam!(gamks₁, d, y₁, k)
                        gamp!(gamperpks₁, d, y₁, k)

                        ChebyTN!(Tbdt₀, N, y₁)
                        mul!(y₁, Tbdt₀, ζcoefk)

                        for j in 1:ns_near
                            dx1 = gamks₁[1, j] - μ₀[1]
                            dx2 = gamks₁[2, j] - μ₀[2]

                            num = dx1 * gamperpks₁[1, j] + dx2 * gamperpks₁[2, j]
                            den = dx1 * dx1 + dx2 * dx2

                            kbd₁[j] = (num * fw_near[j] * dwz₁[j] * inv2π) / den
                        end

                        # Second side.
                        gam!(gamks₂, d, y₂, k)
                        gamp!(gamperpks₂, d, y₂, k)

                        ChebyTN!(Tbdt₀, N, y₂)
                        mul!(y₂, Tbdt₀, ζcoefk)

                        for j in 1:ns_near
                            dx1 = gamks₂[1, j] - μ₀[1]
                            dx2 = gamks₂[2, j] - μ₀[2]

                            num = dx1 * gamperpks₂[1, j] + dx2 * gamperpks₂[2, j]
                            den = dx1 * dx1 + dx2 * dx2

                            kbd₂[j] = (num * fw_near[j] * dwz₂[j] * inv2π) / den
                        end

                        fvals[m] += ((t₀ + 1.0) * dot(kbd₁, y₁) +
                                     (1.0 - t₀) * dot(kbd₂, y₂)) / 2.0
                    end

                else
                    gam!(gamkbdy, d, zr_bdy, k)
                    gamp!(gampkbdy, d, zr_bdy, k)

                    for j in 1:nr_bdy
                        dx1 = gamkbdy[1, j] - μ₀[1]
                        dx2 = gamkbdy[2, j] - μ₀[2]

                        num = dx1 * gampkbdy[1, j] + dx2 * gampkbdy[2, j]
                        den = dx1 * dx1 + dx2 * dx2

                        kbd_bdy[j] = (num * fwr_bdy[j] * inv2π) / den
                    end

                    @views ζb = ζ_bdy[(1+(ll-1)*nr_bdy):(ll*nr_bdy)]

                    fvals[m] += dot(ζb, kbd_bdy)
                end
            end
        end

        v[row] += nevill!(dp.bdxvals, fvals, ϵ)
    end

    # ----- STEP 3: Volumetruc Integrals over all rows -----
    for k in 1:M
        isbdflag = (k in d.kd)

        if isbdflag
            # Zx₂, Zy₂: images of Chebyshev grid (zx2,zy2) on patch k
            mapxy_Dmap!(Zx₂, Zy₂, DJ₂, d, zx2, zy2, k) # nr×nr Zx₂, Zy₂, DJ₂

            hc = d.pths[k].ck1 - d.pths[k].ck0

            Dhc = s >= 0.5 ? hc^(s - 1) : hc^s

            #CN .= reshape(chebcoef[(k - 1) * Np + 1 : k * Np], N, N)
            ak = (k - 1) * Np

            @inbounds for i in 1:Np
                CN[i] = chebcoef[ak+i]
            end

            mul!(CNnr, CN, TzT)   # N × nr
            mul!(Ubd, Tz1, CNnr)            # nbd × nr

            @. Ubd = Ubd * DJ₂
        else
            # Zx, Zy: images of Chebyshev grid (zx,zy) on patch k
            mapxy_Dmap!(Zx, Zy, DJ, d, zx, zy, k) # nr×nr Zx, Zy, DJ

            dfunc!(Df, d, k, zt, s)

            #UFV .= reshape(ufin[(k - 1) * nrp + 1 : k * nrp], nr, nr)
            ak = (k - 1) * nrp

            @inbounds for i in 1:nrp
                UFV[i] = ufin[ak+i]
            end

            @. Ur = UFV * DJ * Df
        end

        for row in 1:Ni

            ℓ = cld(row, Np)
            k == ℓ && continue

            Ikey = packkey(row, k)
            col = get(dp.hmap, Ikey, 0)

            if col != 0
                # ---------- Near-singular patch case ----------
                @views ck = chebcoef[(k-1)*Np+1:k*Np]

                @views NSI = IntS[:, col]

                v[row] += Cs * dot(ck, NSI)
            else
                # ---------- Regular patch case ----------
                x1 = dp.tgtpts[1, row]
                x2 = dp.tgtpts[2, row]

                if isbdflag
                    @. Ker₂ = ((x1 - Zx₂)^2 + (x2 - Zy₂)^2)^(-s)
                    @. KIbd = Ker₂ * Ubd
                    v[row] += Cs * Dhc * dot(mfw, KIbd, fwr)

                else
                    @. Ker = ((x1 - Zx)^2 + (x2 - Zy)^2)^(-s)
                    @. KIr = Ker * Ur
                    # Computes fwr' * KIr * fwr
                    v[row] += Cs * dot(fwr, KIr, fwr)
                end

            end
        end

        if s < 0.5

            Lₚₘ = Ni + 1
            Lₚₙ = Ni + Nb

            for row in Lₚₘ:Lₚₙ

                k₀ = cld(row - Ni, N)

                ℓ = d.kd[k₀]

                k == ℓ && continue

                Ikey = packkey(row, k)
                col = get(dp.hmap, Ikey, 0)

                if col != 0
                    # ---------- Near-singular patch case ----------
                    @views ck = chebcoef[(k-1)*Np+1:k*Np]

                    @views NSI = IntS[:, col]

                    v[row] += Cs * dot(ck, NSI)
                else
                    # ---------- Regular patch case ----------
                    x1 = dp.tgtpts[1, row]
                    x2 = dp.tgtpts[2, row]

                    if isbdflag
                        @. Ker₂ = ((x1 - Zx₂)^2 + (x2 - Zy₂)^2)^(-s)
                        @. KIbd = Ker₂ * Ubd
                        v[row] += Cs * Dhc * dot(mfw, KIbd, fwr)

                    else
                        @. Ker = ((x1 - Zx)^2 + (x2 - Zy)^2)^(-s)
                        @. KIr = Ker * Ur
                        # Computes fwr' * KIr * fwr
                        v[row] += Cs * dot(fwr, KIr, fwr)
                    end

                end
            end

        end
    end

    # ----- STEP 4: Remaining boundary rows -----
    if s >= 0.5

        @inbounds for j in 1:N
            @views Ct = CT[:, j]
            @inbounds for k in 1:Mbd
                row = Ni + (k - 1) * N + j
                ℓ = d.kd[k]
                #CN .= reshape(chebcoef[(ℓ - 1) * Np + 1 : ℓ * Np], N, N)
                aℓ = (ℓ - 1) * Np
                @inbounds for i in 1:Np
                    CN[i] = chebcoef[aℓ+i]
                end
                #ζv is a vector if size N, temporarily being used!
                mul!(ζv, CN, Ct)
                v[row] = sum(ζv)
            end
        end

    else #s<0.5 case

        Lₚₘ = Ni + 1
        Lₚₙ = Ni + Nb

        for row in Lₚₘ:Lₚₙ

            k₀ = cld(row - Ni, N)

            #Boundary patch number
            ptl = d.kd[k₀]

            #Linear index of point on the boundary patch
            ptj = row - Ni - (k₀ - 1) * N

            gam!(γxbd, d, CT[2, ptj], ptl)

            # ---- boundary base term ----
            v[row] += u[row] / 2

            for ll in 1:Mbd
                k = d.kd[ll]
                @views ζv = ζ_bdy[(1+(ll-1)*nr_bdy):(ll*nr_bdy)]

                if k == ptl
                    # Self-panel boundary principal value.
                    DLP!(kbd_bdy, d, CT[2, ptj], zr_bdy, k, μ₀, γt1, γt2)

                    @. kbd_bdy = kbd_bdy * fwr_bdy * inv2π

                    v[row] += dot(ζv, kbd_bdy)

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

                        v[row] += dot(ζv, kbd_bdy)

                    else
                        # Far off-diagonal panel: ordinary direct quadrature.
                        gam!(gamkbdy, d, zr_bdy, k)
                        gamp!(gampkbdy, d, zr_bdy, k)

                        @inbounds for j in 1:nr_bdy
                            dx1 = gamkbdy[1, j] - γxbd[1]
                            dx2 = gamkbdy[2, j] - γxbd[2]

                            num = dx1 * gampkbdy[1, j] + dx2 * gampkbdy[2, j]
                            den = dx1 * dx1 + dx2 * dx2

                            kbd_bdy[j] = (num * fwr_bdy[j] * inv2π) / den
                        end

                        v[row] += dot(ζv, kbd_bdy)
                    end
                end
            end

        end

    end

    return nothing
end

