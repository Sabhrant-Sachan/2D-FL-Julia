function precompsDLP(d::D, dp::domprop, p::Int; 
    n::Int=128)::Matrix{Float64} where {D<:abstractdomain}

    N = dp.N

    Lₚ = size(dp.prepts_bd, 2)

    # ------------------------------------------------------------
    # Quadrature nodes for singular / near-singular splitting.
    #   z1 = cospi((2*(0:n-1)+1)/(4*n)).^2
    #   z2 = sinpi((2*(0:n-1)+1)/(4*n)).^2
    # ------------------------------------------------------------
    z1 = Vector{Float64}(undef, n)
    z2 = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        c = π * (2 * i - 1) / (2 * n)
        z1[i] = cos(c / 2)^2
        z2[i] = sin(c / 2)^2
    end

    fw = getF1W(n)

    # ------------------------------------------------------------
    # Scratch arrays.
    # ------------------------------------------------------------
    y1 = Vector{Float64}(undef, n)
    y2 = Vector{Float64}(undef, n)

    Iv1 = Vector{Float64}(undef, n)
    Iv2 = Vector{Float64}(undef, n)

    I1 = Vector{Float64}(undef, N)
    I2 = Vector{Float64}(undef, N)

    # Boundary coordinates and derivatives.
    TNy = Matrix{Float64}(undef, n, N)

    IntS = Matrix{Float64}(undef, N, Lₚ)

    gamk = Matrix{Float64}(undef, 2, n)
    # gamp! returns γᵖᵉʳᵖ = (γ₂′, -γ₁′), so this is (γ - x) ⋅ γᵖᵉʳᵖ.
    gamperpk = Matrix{Float64}(undef, 2, n)

    # ------------------------------------------------------------
    # Store w(-z1), w(-z2), w'(z1), w'(z2).
    # ------------------------------------------------------------
    wmz₁ = Vector{Float64}(undef, n)
    wmz₂ = Vector{Float64}(undef, n)
    dwz₁ = Vector{Float64}(undef, n)
    dwz₂ = Vector{Float64}(undef, n)

    wfunc!(wmz₁, p, z1; α=-1.0)
    wfunc!(wmz₂, p, z2; α=-1.0)

    dwfunc!(dwz₁, p, z1)
    dwfunc!(dwz₂, p, z2)

    inv2π = 1 / (2π)

    # ------------------------------------------------------------
    # Main precomputation loop.
    # ------------------------------------------------------------
    @inbounds for i in 1:Lₚ

        # col is the interioir point target index.
        col = dp.prepts_bd[1, i]

        # k is the actual integration boundary patch.
        k = dp.prepts_bd[2, i]

        x1 = dp.tgtpts[1, col]
        x2 = dp.tgtpts[2, col]

        t₀ = dp.projpts[i]

        t₀x = gamx(d, t₀, k)
        t₀y = gamy(d, t₀, k)

        dist_pt_to_bd = hypot(t₀x - x1, t₀y - x2)

        if dist_pt_to_bd < 5e-3
            #Interpolation needed!

        else

            if t₀ == 1.0
                # So only the y1-side contributes.
                @. y1 = 1.0 - 2.0 * wmz₁

                gam!(gamk, d, y1, k)
                gamp!(gamperpk, d, y1, k)

                @inbounds for j in 1:n
                    dx1 = gamk[1, j] - x1
                    dx2 = gamk[2, j] - x2

                    num = dx1 * gamperpk[1, j] + dx2 * gamperpk[2, j]
                    den = dx1 * dx1 + dx2 * dx2

                    Iv1[j] = (num * fw[j] * dwz₁[j] * inv2π) / den
                end

                ChebyTN!(TNy, N, y1)

                mul!(I1, transpose(TNy), Iv1)

                @inbounds for row in 1:N
                    IntS[row, i] = I1[row]
                end

            elseif t₀ == -1.0
                # So only the y2-side contributes.
                @. y2 = -1.0 + 2.0 * wmz₂

                gam!(gamk, d, y2, k)
                gamp!(gamperpk, d, y2, k)

                @inbounds for j in 1:n
                    dx1 = gamk[1, j] - x1
                    dx2 = gamk[2, j] - x2

                    num = dx1 * gamperpk[1, j] + dx2 * gamperpk[2, j]
                    den = dx1 * dx1 + dx2 * dx2

                    Iv2[j] = (num * fw[j] * dwz₂[j] * inv2π) / den
                end

                ChebyTN!(TNy, N, y2)

                mul!(I2, transpose(TNy), Iv2)

                @inbounds for row in 1:N
                    IntS[row, i] = I2[row]
                end

            else
                # Both sides contribute:
                @. y1 = t₀ - (t₀ + 1.0) * wmz₁
                @. y2 = t₀ + (1.0 - t₀) * wmz₂

                # First side: y1.
                gam!(gamk, d, y1, k)
                gamp!(gamperpk, d, y1, k)

                @inbounds for j in 1:n
                    dx1 = gamk[1, j] - x1
                    dx2 = gamk[2, j] - x2

                    num = dx1 * gamperpk[1, j] + dx2 * gamperpk[2, j]
                    den = dx1 * dx1 + dx2 * dx2

                    Iv1[j] = (num * fw[j] * dwz₁[j] * inv2π) / den
                end

                ChebyTN!(TNy, N, y1)

                mul!(I1, transpose(TNy), Iv1)

                # Second side: y2.
                gam!(gamk, d, y2, k)
                gamp!(gamperpk, d, y2, k)

                @inbounds for j in 1:n
                    dx1 = gamk[1, j] - x1
                    dx2 = gamk[2, j] - x2

                    num = dx1 * gamperpk[1, j] + dx2 * gamperpk[2, j]
                    den = dx1 * dx1 + dx2 * dx2

                    Iv2[j] = (num * fw[j] * dwz₂[j] * inv2π) / den
                end

                ChebyTN!(TNy, N, y2)

                mul!(I2, transpose(TNy), Iv2)

                @inbounds for row in 1:N
                    IntS[row, i] = ((t₀ + 1.0) * I1[row] + (1.0 - t₀) * I2[row]) / 2.0
                end
            end

        end

    end

    return IntS

end

function testDLPint(d::D, dp::domprop, IntS::Matrix{Float64}; 
    nreg::Int = 64, plot_err::Bool=false)::Float64 where {D<:abstractdomain}

    M = d.Npat
    N = dp.N
    Mbd = length(d.kd)
    Ni = M * N^2

    # ------------------------------------------------------------
    # Quadrature nodes for regular integration.
    #   z= cospi((2*(0:n-1)+1)/(2*n))
    # ------------------------------------------------------------
    n = nreg
    fw = getF1W(n)

    z = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        z[i] = cos(π * (2 * i - 1) / (2 * n))
    end
    # ------------------------------------------------------------
    # Scratch arrays.
    # ------------------------------------------------------------
    Iv = Vector{Float64}(undef, n)

    gamk = Matrix{Float64}(undef, 2, n)
    # gamp! returns γᵖᵉʳᵖ = (γ₂′, -γ₁′), so this is (γ - x) ⋅ γᵖᵉʳᵖ.
    gamperpk = Matrix{Float64}(undef, 2, n)

    # v contains the output of the DLP evaluation 
    # for each interioir point
    v = Vector{Float64}(undef, Ni)
    fill!(v, 0.0)

    inv2π = 1 / (2π)

    for kb in 1:Mbd
        k = d.kd[kb]

        gam!(gamk, d, z, k)
        gamp!(gamperpk, d, z, k)

        for i in 1:Ni
            x1 = dp.tgtpts[1, i]
            x2 = dp.tgtpts[2, i]

            Ikey = packkey(i, k)
            col = get(dp.hmap_bd, Ikey, 0)

            if col != 0
                # ---------- Near-singular patch case ----------
                v[i] += IntS[1, col]
            else
                x1 = dp.tgtpts[1, i]
                x2 = dp.tgtpts[2, i]

                @inbounds for j in 1:n
                    dx1 = gamk[1, j] - x1
                    dx2 = gamk[2, j] - x2

                    num = dx1 * gamperpk[1, j] + dx2 * gamperpk[2, j]
                    den = dx1 * dx1 + dx2 * dx2

                    Iv[j] = (num * inv2π) / den
                end
                v[i] += dot(fw, Iv)
            end
        end
    end

    #Now compare it with the exact value (which is 1)

    Errv = abs.(v .- 1.0)

    Err, imax = findmax(Errv)

    @printf("Error in DLP forward : %.2e \n", Err)

    @printf("Worst target index: %d, computed value: %.16e \n", imax, v[imax])

    if plot_err
        plotfunc(dp, d, Errv)
    end
    
    return Err

end