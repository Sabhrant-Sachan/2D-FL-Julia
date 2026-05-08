function Axbdop!(v::SubArray{Float64}, kI::Int, d::D, dp::domprop, 
    s::Float64, IV::IVT) where {D<:abstractdomain, IVT}
    # --------- Unpack IV (NamedTuple) ---------
    (; IV1, IVbdth, IVbd, IVbt2, IVbdt1, IVbdt2, IVbdt3) = IV
    (; N, Np, M, Mbd) = IV1
    (; N₁, N₂, fw₁, fw₂, idctbd₁, idctbd₂, dₙₕ) = IVbd 
    (; y₁, y₂, gamk₁, gamp₁, gamk₂, gamp₂) = IVbdt1
    (; kbd₁, kbd₂, γx, μ₀, γt1, γt2) = IVbdt2
    (; Lᵢₙ, xvals, fvals, γxbd, kdx, Tbd) = IVbdt3
    (; coeffs) = IVbdth

    # v is M*Np + Mbd*N by N matrix (preallocated given)
    # that is, size(v) == (M*Np + Mbd*N, N)

    k = d.kd[kI];
    # --------- Compute boundary curve points ---------
    gam!(gamk₁, d, y₁, k)
    gam!(gamk₂, d, y₂, k)
    gamp!(gamp₁, d, y₁, k)
    gamp!(gamp₂, d, y₂, k)

    # --------- interior rows  ---------
    Ni = M * Np

    ℓbd = 0

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

                rview = @view v[row, :]

                mul!(rview, transpose(idctbd₁), kbd₁, 1/(2π), 0.0)

            else

                @inbounds for j in 1:N₂
                    dx1 = gamk₂[1, j] - γx[1]
                    dx2 = gamk₂[2, j] - γx[2]

                    num = dx1 * gamp₂[1, j] + dx2 * gamp₂[2, j]
                    den = dx1 * dx1 + dx2 * dx2

                    kbd₂[j] = (num * fw₂[j]) / den
                end

                rview = @view v[row, :]

                mul!(rview, transpose(idctbd₂), kbd₂, 1/(2π), 0.0)

            end
        else

            ℓbd += 1 
            # Point close to the boundary, Interpolation!
            l = dp.distpts[ℓbd].l

            ϵ = dp.distpts[ℓbd].d

            knbh!(kdx, d, l, dₙₕ, γt1, γt2)

            if !(k in kdx) 

                @inbounds for j in 1:N₁
                    dx1 = gamk₁[1, j] - γx[1]
                    dx2 = gamk₁[2, j] - γx[2]

                    num = dx1 * gamp₁[1, j] + dx2 * gamp₁[2, j]
                    den = dx1 * dx1 + dx2 * dx2

                    kbd₁[j] = (num * fw₁[j]) / den
                end

                rview = @view v[row, :]

                mul!(rview, transpose(idctbd₁), kbd₁, 1/(2π), 0.0)

            else 
                tₛ = dp.distpts[ℓbd].t

                ChebyTN!(Tbd, N, tₛ)

                gam!(γxbd, d, tₛ, l)

                DLP!(kbd₁, d, tₛ, l, y₁, k, μ₀, γt1, γt2)

                nu!(γt1, d, tₛ, l)

                @. kbd₁ = kbd₁ * fw₁

                @inbounds for j in 1:N

                    @views ζ1 = idctbd₁[:,j]

                    @views ζ2 = idctbd₂[:,j]

                    fvals[1] = dot(ζ1, kbd₁) / (2π)

                    if k==l
                        # Cbd .= ζ2    
                        # mul!(Cbd, pCbd, Cbd)
                        # @. Cbd = Cbd / N₂
                        # Cbd[1] = Cbd[1] / 2
                        # ChebyTN!(Tbd₂, N₂, tₛ)
                        # ζₛ = dot(Tbd₂, Cbd)
                        # fvals[1] += ζₛ/2 

                        @views Cbdj = coeffs[:, j]

                        ζₛ = dot(Tbd, Cbdj)

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

                        fvals[i] = dot(ζ2, kbd₂) / (2π)

                    end

                    v[row, j] = nevill!(xvals, fvals, ϵ)
                end

            end

        end

    end

    if s<0.5

        Nb = Mbd * N
        Lₚₘ = Ni + 1
        Lₚₙ = Ni + Nb
        CT = IVbt2.CT

        for row in Lₚₘ : Lₚₙ

            k₀ = cld(row - Ni, N)
            
            ptl = d.kd[k₀]

            ptj = row - Ni - (k₀ - 1) * N

            rview = @view v[row, :]

            DLP!(kbd₁, d, CT[2,ptj], ptl, y₁, k, μ₀, γt1, γt2)

            @. kbd₁ = kbd₁ * fw₁

            mul!(rview, transpose(idctbd₁), kbd₁, 1/(2π), 0.0)
        end
        #The limiting value from the interioir is 1/2,
        for j in 1:N
            v[Ni+N*(kI-1)+j, j] += 0.5
        end
    end

end