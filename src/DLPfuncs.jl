function testDLP(d::D, dp::domprop; p::Int=2, nr::Int=64, 
    plot_err::Bool=true)::Float64 where {D<:abstractdomain}

    #This is just for testing the DLP evaluation, and not the solver.
    #Zeta density is taken as 1. 
    M = d.Npat
    N = dp.N
    Mbd = length(d.kd)
    Ni = M * N^2

    v = zeros(Float64, Ni)
    # ------------------------------------------------------------
    # Quadrature nodes for regular integration.
    #   zr= cospi((2*(0:n-1)+1)/(2*n))
    # ------------------------------------------------------------
    fwr = getF1W(nr)

    zr = Vector{Float64}(undef, nr)

    @inbounds for i in 1:nr
        zr[i] = cos(π * (2 * i - 1) / (2 * nr))
    end

    gamkr = Matrix{Float64}(undef, 2, nr)
    gamperpkr = Matrix{Float64}(undef, 2, nr)
    ker_reg = Vector{Float64}(undef, nr)
    # ------------------------------------------------------------
    # Quadrature nodes for singular / near-singular splitting.
    #   z1 = cospi((2*(0:n-1)+1)/(4*n)).^2
    #   z2 = sinpi((2*(0:n-1)+1)/(4*n)).^2
    # ------------------------------------------------------------
    n = 2*nr
    fw = getF1W(n)

    z1 = Vector{Float64}(undef, n)
    z2 = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        c = π * (2 * i - 1) / (2 * n)
        z1[i] = cos(c / 2)^2
        z2[i] = sin(c / 2)^2
    end

    # ------------------------------------------------------------
    # Scratch arrays.
    # ------------------------------------------------------------
    y1 = Vector{Float64}(undef, n)
    y2 = Vector{Float64}(undef, n)

    yt1 = Vector{Float64}(undef, n)
    yt2 = Vector{Float64}(undef, n)

    Iv1 = Vector{Float64}(undef, n)
    Iv2 = Vector{Float64}(undef, n)

    gamks₁ = Matrix{Float64}(undef, 2, n)
    gamks₂ = Matrix{Float64}(undef, 2, n)
    gamkt1 = Matrix{Float64}(undef, 2, n)
    gamkt2 = Matrix{Float64}(undef, 2, n)
    # gamp! returns γᵖᵉʳᵖ = (γ₂′, -γ₁′), so this is (γ - x) ⋅ γᵖᵉʳᵖ.
    gamperpks₁ = Matrix{Float64}(undef, 2, n)
    gamperpks₂ = Matrix{Float64}(undef, 2, n)
    gamperpkt1 = Matrix{Float64}(undef, 2, n)
    gamperpkt2 = Matrix{Float64}(undef, 2, n)

    μ₀    = Vector{Float64}(undef, 2)
    γt1   = Vector{Float64}(undef, 2)
    γt2   = Vector{Float64}(undef, 2)
    
    #For interpolating 
    fvals = Vector{Float64}(undef, dp.Lᵢₙ)
    γxbd  = Vector{Float64}(undef, 2)

    # ------------------------------------------------------------
    # Store w(-z1), w(-z2), w'(z1), w'(z2).
    # ------------------------------------------------------------
    wmz₁ = Vector{Float64}(undef, n)
    wmz₂ = Vector{Float64}(undef, n)
    dwz₁ = Vector{Float64}(undef, n)
    dwz₂ = Vector{Float64}(undef, n)

    wfunc!(wmz₁, p, z1; α=-1.0)
    wfunc!(wmz₂, p, z2; α=-1.0)

    @. yt1 = 1.0 - 2.0 * wmz₁
    @. yt2 = -1.0 + 2.0 * wmz₂

    dwfunc!(dwz₁, p, z1)
    dwfunc!(dwz₂, p, z2)

    inv2π = 1 / (2π)

    BD_FAR = UInt8(0)
    BD_NEAR = UInt8(1)
    BD_INTP = UInt8(2)

    for ll in 1:Mbd
        k = d.kd[ll]
        jnear = 0

        gam!(gamkr, d, zr, k)
        gamp!(gamperpkr, d, zr, k)

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
                @inbounds for j in 1:nr
                    dx1 = gamkr[1, j] - x1
                    dx2 = gamkr[2, j] - x2

                    num = dx1 * gamperpkr[1, j] + dx2 * gamperpkr[2, j]
                    den = dx1 * dx1 + dx2 * dx2

                    ker_reg[j] = (num * fwr[j] * inv2π) / den
                end

                v[row] = v[row] + sum(ker_reg)

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
                            
                            @inbounds for j in 1:n
                                dx1 = gamkt1[1, j] - x1
                                dx2 = gamkt1[2, j] - x2

                                num = dx1 * gamperpkt1[1, j] + dx2 * gamperpkt1[2, j]
                                den = dx1 * dx1 + dx2 * dx2

                                Iv1[j] = (num * fw[j] * dwz₁[j] * inv2π) / den
                            end

                            v[row] = v[row] + sum(Iv1)

                        elseif t₀ == -1.0
                            # So only the y2-side contributes.

                            @inbounds for j in 1:n
                                dx1 = gamkt2[1, j] - x1
                                dx2 = gamkt2[2, j] - x2

                                num = dx1 * gamperpkt2[1, j] + dx2 * gamperpkt2[2, j]
                                den = dx1 * dx1 + dx2 * dx2

                                Iv2[j] = (num * fw[j] * dwz₂[j] * inv2π) / den
                            end

                            v[row] = v[row] + sum(Iv2)
                        else
                            # Both sides contribute:
                            @. y1 = t₀ - (t₀ + 1.0) * wmz₁
                            @. y2 = t₀ + (1.0 - t₀) * wmz₂

                            # First side: y1.
                            gam!(gamks₁, d, y1, k)
                            gamp!(gamperpks₁, d, y1, k)

                            @inbounds for j in 1:n
                                dx1 = gamks₁[1, j] - x1
                                dx2 = gamks₁[2, j] - x2

                                num = dx1 * gamperpks₁[1, j] + dx2 * gamperpks₁[2, j]
                                den = dx1 * dx1 + dx2 * dx2

                                Iv1[j] = (num * fw[j] * dwz₁[j] * inv2π) / den
                            end

                            # Second side: y2.
                            gam!(gamks₂, d, y2, k)
                            gamp!(gamperpks₂, d, y2, k)

                            @inbounds for j in 1:n
                                dx1 = gamks₂[1, j] - x1
                                dx2 = gamks₂[2, j] - x2

                                num = dx1 * gamperpks₂[1, j] + dx2 * gamperpks₂[2, j]
                                den = dx1 * dx1 + dx2 * dx2

                                Iv2[j] = (num * fw[j] * dwz₂[j] * inv2π) / den
                            end

                            v[row] = v[row] + ((t₀ + 1.0) * sum(Iv1) + (1.0 - t₀) * sum(Iv2)) / 2.0

                        end
                        
                        break
                    end
                end

                if !near_flag
                    @inbounds for j in 1:nr
                        dx1 = gamkr[1, j] - x1
                        dx2 = gamkr[2, j] - x2

                        num = dx1 * gamperpkr[1, j] + dx2 * gamperpkr[2, j]
                        den = dx1 * dx1 + dx2 * dx2

                        ker_reg[j] = (num * fwr[j] * inv2π) / den
                    end

                    v[row] = v[row] + sum(ker_reg)
                end

            end
                
        end

    end

    jintp = 0

    @inbounds for row in 1:Ni
        dp.bdmode[row] != BD_INTP && continue

        # Mode-2 point: very close to the boundary.
        # We evaluate full DLP at interpolation nodes, then interpolate
        # back to the true distance ε.
        jintp += 1

        l = dp.bdclosest[jintp]
        tₛ = dp.bdt[jintp]
        ϵ = dp.bddist[jintp]

        gam!(γxbd, d, tₛ, l)

        # ------------------------------------------------------------
        # m = 1: boundary interpolation value.
        #
        # For ζ ≡ 1:
        #     interior boundary limit = principal value + 1/2.
        #
        # Self panel k == l:
        #     use diagonal DLP formula with target parameter tₛ on patch l.
        #
        # Near off-diagonal panels stored in bdintp for m = 1:
        #     bdintpt stores extended inverse s satisfying γ_k(s)=γ_l(tₛ).
        #     Use diagonal/continued formula on patch k.
        #
        # Far off-diagonal panels:
        #     use ordinary off-diagonal formula.
        # ------------------------------------------------------------
        fvals[1] = 0.5

        a = (jintp - 1) * dp.Lᵢₙ + 1

        q1 = dp.bdintpptr[a]
        q2 = dp.bdintpptr[a+1] - 1

        for ll in 1:Mbd
            k = d.kd[ll]

            if k == l
                # Self-panel boundary principal value.
                DLP!(ker_reg, d, tₛ, zr, k, μ₀, γt1, γt2)

                @. ker_reg = ker_reg * fwr

                fvals[1] += sum(ker_reg) * inv2π

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
                    # Use the extended inverse parameter on patch k.
                    DLP!(ker_reg, d, s_inv, zr, k, μ₀, γt1, γt2)

                    @. ker_reg = ker_reg * fwr

                    fvals[1] += sum(ker_reg) * inv2π

                else
                    # Far off-diagonal panel: ordinary direct quadrature.
                    gam!(gamkr, d, zr, k)
                    gamp!(gamperpkr, d, zr, k)

                    @inbounds for j in 1:nr
                        dx1 = gamkr[1, j] - γxbd[1]
                        dx2 = gamkr[2, j] - γxbd[2]

                        num = dx1 * gamperpkr[1, j] + dx2 * gamperpkr[2, j]
                        den = dx1 * dx1 + dx2 * dx2

                        ker_reg[j] = (num * fwr[j] * inv2π) / den
                    end

                    fvals[1] += sum(ker_reg)
                end
            end
        end

        # ------------------------------------------------------------
        # m = 2:Lᵢₙ: auxiliary interior interpolation nodes.
        #
        # Here bdintpt stores the ordinary closest projected parameter
        # onto each near panel.
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

                        @inbounds for j in 1:n
                            dx1 = gamks₁[1, j] - μ₀[1]
                            dx2 = gamks₁[2, j] - μ₀[2]

                            num = dx1 * gamperpks₁[1, j] + dx2 * gamperpks₁[2, j]
                            den = dx1 * dx1 + dx2 * dx2

                            Iv1[j] = (num * fw[j] * dwz₁[j] * inv2π) / den
                        end

                        fvals[m] += sum(Iv1)

                    elseif t₀ == -1.0
                        gam!(gamks₂, d, yt2, k)
                        gamp!(gamperpks₂, d, yt2, k)

                        @inbounds for j in 1:n
                            dx1 = gamks₂[1, j] - μ₀[1]
                            dx2 = gamks₂[2, j] - μ₀[2]

                            num = dx1 * gamperpks₂[1, j] + dx2 * gamperpks₂[2, j]
                            den = dx1 * dx1 + dx2 * dx2

                            Iv2[j] = (num * fw[j] * dwz₂[j] * inv2π) / den
                        end

                        fvals[m] += sum(Iv2)

                    else
                        @. y1 = t₀ - (t₀ + 1.0) * wmz₁
                        @. y2 = t₀ + (1.0 - t₀) * wmz₂

                        gam!(gamks₁, d, y1, k)
                        gamp!(gamperpks₁, d, y1, k)

                        @inbounds for j in 1:n
                            dx1 = gamks₁[1, j] - μ₀[1]
                            dx2 = gamks₁[2, j] - μ₀[2]

                            num = dx1 * gamperpks₁[1, j] + dx2 * gamperpks₁[2, j]
                            den = dx1 * dx1 + dx2 * dx2

                            Iv1[j] = (num * fw[j] * dwz₁[j] * inv2π) / den
                        end

                        gam!(gamks₂, d, y2, k)
                        gamp!(gamperpks₂, d, y2, k)

                        @inbounds for j in 1:n
                            dx1 = gamks₂[1, j] - μ₀[1]
                            dx2 = gamks₂[2, j] - μ₀[2]

                            num = dx1 * gamperpks₂[1, j] + dx2 * gamperpks₂[2, j]
                            den = dx1 * dx1 + dx2 * dx2

                            Iv2[j] = (num * fw[j] * dwz₂[j] * inv2π) / den
                        end

                        fvals[m] += ((t₀ + 1.0) * sum(Iv1) + (1.0 - t₀) * sum(Iv2)) / 2.0
                    end

                else
                    gam!(gamkr, d, zr, k)
                    gamp!(gamperpkr, d, zr, k)

                    @inbounds for j in 1:nr
                        dx1 = gamkr[1, j] - μ₀[1]
                        dx2 = gamkr[2, j] - μ₀[2]

                        num = dx1 * gamperpkr[1, j] + dx2 * gamperpkr[2, j]
                        den = dx1 * dx1 + dx2 * dx2

                        ker_reg[j] = (num * fwr[j] * inv2π) / den
                    end

                    fvals[m] += sum(ker_reg)
                end
            end
        end

        v[row] = nevill!(dp.bdxvals, fvals, ϵ)
    end

    Errv = abs.(v .- 1.0)
    Err, imax = findmax(Errv)

    @printf("Error in DLP forward : %.2e \n", Err)
    @printf("Worst target index: %d, computed value: %.16e \n", imax, v[imax])

    if plot_err
        plotfunc(dp, d, Errv)
    end

    return Err

end
