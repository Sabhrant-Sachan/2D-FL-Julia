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
    (; IV1, IVr, IVbdth, IVbd, IVt, IVbt1, IVbt2, IVbdt1, IVbdt2, IVbdt3, IVAf, IVAdct) = IV

    (; p_dct2_dim2, p_dct2_dim1, p_dct3_dim2, p_dct3_dim1, p_dct2_N, p_dct3_N1, p_dct3_N2)= IVAdct

    (; chebcoef, ufin, ζ₁, ζ₂, ζ₂coeff, UV, UFV, ζv, ζfv₁, ζfv₂, CNnr,  TzT) = IVAf

    (; N, Np, Cs, M, Mbd) = IV1

    (; nr, fwr, nrp, zx, zt, zy, Df) = IVr

    (; zx2, zy2, Tz1) = IVbdth
    
    (; N₁, N₂, fw₁, fw₂, dₙₕ) = IVbd

    (; Zx, Zy, DJ, Ker, KIr, Ur, CN) = IVt

    (; Zx₂, Zy₂, DJ₂, Ker₂, KIbd) = IVbt1

    (; Ubd, mfw, CT) = IVbt2

    (; y₁, y₂, gamk₁, gamp₁, gamk₂, gamp₂) = IVbdt1

    (; kbd₁, kbd₂, γx, μ₀, γt1, γt2) = IVbdt2

    (; Lᵢₙ, xvals, fvals, γxbd, kdx, Tbd₂) = IVbdt3
     
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
            UV[i] = u[ak₁ + i]
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
            chebcoef[ak₁ + i] = UV[i]
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
            ufin[ak₂ + i] = UFV[i]
        end

    end
    
    ζ₂coeff .= 0 

    # ----- Boundary ζ: coefficients + values on N₁ and N₂ grids -----
    # For k=1:Mbd boundary patches, u has ζ values on N Cheby nodes for each bd patch.
    @inbounds for k in 1:Mbd
        ak = Ni + (k - 1) * N
        i₁ = (k - 1) * N₁
        i₂ = (k - 1) * N₂
        # load ζv from u boundary values
        @inbounds for i in 1:N
            ζv[i] = u[ak + i]
        end

        mul!(ζv, p_dct2_N, ζv)

        ζv .*= 1 / N

        ζv[1] /= 2

        # store these N coeffs
        @inbounds for i in 1:N
            chebcoef[ak + i]  = ζv[i]
            ζ₂coeff[i₂ + i] = ζv[i]
        end

        # ζv are now the Chebyshev coefficients of ζ 
        # Evaluate ζ finer values over N₁ and N₂ 
        ζv[1] *= 2

        ζfv₁ .= 0 
        ζfv₂ .= 0

        # ζfv₁[1:N] .= ζv
        # ζfv₂[1:N] .= ζv
        @inbounds for i in 1:N
            ζfv₁[i] = ζv[i]
            ζfv₂[i] = ζv[i]
        end
        
        mul!(ζfv₁, p_dct3_N1, ζfv₁) 
        mul!(ζfv₂, p_dct3_N2, ζfv₂) 

        ζfv₁ .*= 0.5
        ζfv₂ .*= 0.5

        @inbounds for i in 1:N₁
            ζ₁[i₁ + i] = ζfv₁[i]
        end
        
        @inbounds for i in 1:N₂
            ζ₂[i₂ + i] = ζfv₂[i]
        end

    end

    # ----- STEP 2: Singular Integrals over all rows and DLP contributions -----
    @inbounds for row in 1:Ni

        ℓ = cld(row, Np)
        j = row - (ℓ - 1) * Np

        col = dp.pthgo[ℓ] + j - 1

        @views cf = chebcoef[((ℓ - 1) * Np + 1):(ℓ * Np)] 

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
        ℓ = d.kd[ll]
        ℓbd = 0

        gam!(gamk₁, d, y₁, ℓ)
        gamp!(gamp₁, d, y₁, ℓ)
        gam!(gamk₂, d, y₂, ℓ)
        gamp!(gamp₂, d, y₂, ℓ)

        #Contribution of all points for the ℓth patch
        @inbounds for row in 1:Ni

            γx[1] = dp.tgtpts[1, row]
            γx[2] = dp.tgtpts[2, row]

            if !dp.tgtbdbm[1, row]
                # Point far away from the boundary, direct computation
                if dp.tgtbdbm[2, row]

                    @inbounds for j in 1:N₁
                        dx1 = gamk₁[1, j] - γx[1]
                        dx2 = gamk₁[2, j] - γx[2]

                        num = dx1 * gamp₁[1, j] + dx2 * gamp₁[2, j]
                        den = dx1 * dx1 + dx2 * dx2

                        kbd₁[j] = (num * fw₁[j]) / den
                    end

                    @views ζv₁ = ζ₁[(1 + (ll - 1) * N₁) : (ll * N₁)]

                    v[row] = v[row] + dot(ζv₁, kbd₁) / (2 * pi)

                else

                    @inbounds for j in 1:N₂
                        dx1 = gamk₂[1, j] - γx[1]
                        dx2 = gamk₂[2, j] - γx[2]

                        num = dx1 * gamp₂[1, j] + dx2 * gamp₂[2, j]
                        den = dx1 * dx1 + dx2 * dx2

                        kbd₂[j] = (num * fw₂[j]) / den
                    end

                    @views ζv₂ = ζ₂[(1 + (ll - 1) * N₂) : (ll * N₂)]

                    v[row] = v[row] + dot(ζv₂, kbd₂) / (2 * pi)

                end
            else

                ℓbd += 1
                # Point close to the boundary, Interpolation!
                l = dp.distpts[ℓbd].l
                ϵ = dp.distpts[ℓbd].d

                knbh!(kdx, d, l, dₙₕ, γt1, γt2)

                if !(ℓ in kdx)

                    @inbounds for j in 1:N₁
                        dx1 = gamk₁[1, j] - γx[1]
                        dx2 = gamk₁[2, j] - γx[2]

                        num = dx1 * gamp₁[1, j] + dx2 * gamp₁[2, j]
                        den = dx1 * dx1 + dx2 * dx2

                        kbd₁[j] = (num * fw₁[j]) / den
                    end

                    @views ζv₁ = ζ₁[(1 + (ll - 1) * N₁):(ll * N₁)]

                    v[row] = v[row] + dot(ζv₁, kbd₁) / (2 * pi)

                else
                    tₛ = dp.distpts[ℓbd].t

                    ChebyTN!(Tbd₂, N₂, tₛ)

                    @views Cbd = ζ₂coeff[(1 + (ll - 1) * N₂):(ll * N₂)]

                    gam!(γxbd, d, tₛ, l)

                    DLP!(kbd₁, d, tₛ, l, y₁, ℓ, μ₀, γt1, γt2)

                    nu!(γt1, d, tₛ, l)

                    @. kbd₁ = kbd₁ * fw₁

                    @views ζv₁ = ζ₁[(1+(ll-1)*N₁):(ll*N₁)]

                    @views ζv₂ = ζ₂[(1+(ll-1)*N₂):(ll*N₂)]

                    fvals[1] = dot(ζv₁, kbd₁) / (2π)

                    if ℓ == l

                        ζₛ = dot(Tbd₂, Cbd)

                        fvals[1] += ζₛ / 2
                    end

                    for i in 2:Lᵢₙ

                        @. μ₀ = γxbd - xvals[i] * γt1

                        @inbounds for jj in 1:N₂
                            dx1 = gamk₂[1, jj] - μ₀[1]
                            dx2 = gamk₂[2, jj] - μ₀[2]

                            num = dx1 * gamp₂[1, jj] + dx2 * gamp₂[2, jj]
                            den = dx1 * dx1 + dx2 * dx2

                            kbd₂[jj] = (num * fw₂[jj]) / den
                        end

                        fvals[i] = dot(ζv₂, kbd₂) / (2π)

                    end

                    v[row] = v[row] + nevill!(xvals, fvals, ϵ)

                end

            end

        end
    end

    # ----- STEP 3: Volumetruc Integrals over all rows
    #               and DLP contributions in interioir  -----
    for k in 1:M
        isbdflag = (k in d.kd)

        if isbdflag
            # Zx₂, Zy₂: images of Chebyshev grid (zx2,zy2) on patch k
            mapxy_Dmap!(Zx₂, Zy₂, DJ₂, d, zx2, zy2, k) # nr×nr Zx₂, Zy₂, DJ₂

            hc = d.pths[k].ck1 - d.pths[k].ck0

            Dhc = s>=0.5 ? hc^(s-1) : hc^s 

            #CN .= reshape(chebcoef[(k - 1) * Np + 1 : k * Np], N, N)
            ak = (k - 1) * Np

            @inbounds for i in 1:Np
                CN[i] = chebcoef[ak + i]
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
                UFV[i] = ufin[ak + i]
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
                @views ck = chebcoef[(k - 1) * Np + 1 : k * Np]

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
                row = Ni + (k - 1)*N + j
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

        for row in Lₚₘ : Lₚₙ

            k₀ = cld(row - Ni, N)
            
            #Boundary patch number
            ptl = d.kd[k₀]

            #Linear index of point on the boundary patch
            ptj = row - Ni - (k₀ - 1) * N

            # ---- boundary base term ----
            for ll = 1:Mbd
                ℓ = d.kd[ll]

                @views ζvll = ζ₁[1 + (ll - 1) * N₁ : ll * N₁]

                DLP!(kbd₁, d, CT[2,ptj], ptl, y₁, ℓ, μ₀, γt1, γt2)

                @. kbd₁ = kbd₁ * fw₁

                v[row] += dot(ζvll, kbd₁)/(2π)
            end

            v[row] += u[row] / 2
            
        end

    end

    return nothing
end

