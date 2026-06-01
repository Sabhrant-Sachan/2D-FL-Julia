@inline function fill_meshgrid!(T1, T2, y1, y2)
   @inbounds for j in eachindex(y2)
      @inbounds for i in eachindex(y1)
         T1[i, j] = y1[i]
         T2[i, j] = y2[j]
      end
   end
   return nothing
end

function bvec(d::abstractdomain, dp::domprop, s::Float64, f!::Function
   ; n::Int=64)::Vector{Float64}

   # n₂ is defined as the number of nodes for integration of regular integrals
   # fw₂is defined as the Fejer 1st weights for corresponding nodes
   n₂ = n
   fw₂ = getF1W(n₂)

   # Increasing n by twice, which denotes nodes for Singular integration
   # corresponding to which we take fejer 1st weights fw
   n = 2n
   fw = getF1W(n)

   # The polynomial change-of-variables with parameter p
   p = 4

   #bookkeeping
   M = d.Npat
   Mbd = length(d.kd)
   N = dp.N
   Np = N^2
   Nd = Mbd * N

   # Chebyshev nodes (These are column vectors)
   z = cospi.((2 .* (0:n₂-1) .+ 1) ./ (2n₂))        # length n2 
   z1 = cospi.((2 .* (0:n-1) .+ 1) ./ (4n)) .^ 2      # length n
   z2 = sinpi.((2 .* (0:n-1) .+ 1) ./ (4n)) .^ 2      # length n
   zp = cospi.((2 .* (0:N-1) .+ 1) ./ (2N))         # length N

   zx = repeat(z, 1, n₂)
   zy = repeat(z', n₂, 1)

   # Initializing b vector
   Ni = M * Np
   
   b = zeros(Float64, Ni + Nd)
   # modified weights
   fwm = similar(fw)
   dwfunc!(fwm, p, z2)   # fwm := dw(z2)
   @. fwm = fw * fwm

   # whether we also fill RHS on boundary points
   is_s_small = s < 0.5 ? 1 : 0

   # ----------- reusable work buffers (singular, size n or n×n) ----------
   d1 = similar(z2)              # n length
   d2 = similar(z2)              # n length
   dz = similar(z2)              # n length
   y1 = similar(z2)              # n length
   y2 = similar(z2)              # n length
   t1 = Matrix{Float64}(undef, n, n)  # meshgrid of y1/y2 (column = y2[j], row = y1[i])
   t2 = Matrix{Float64}(undef, n, n)

   wfunc!(dz, p, z2; α=-1.0) #dz = w(-z2)

   Zx = similar(t1)              # mapped x on singular grid
   Zy = similar(t1)              # mapped y on singular grid
   DJ = similar(t1)              # Jacobian |detJ|
   DIF = similar(t1)              # ‖τ(u,v)−τ(u2,v2)‖ (via diff_map!)
   LOG = similar(t1)              # log(DIF)
   FXY = similar(t1)              # f(Zx, Zy)

   # weights that change per case
   wL = similar(fw)               # n length
   wR = similar(fw)               # n length

   # ---------- reusable work buffers (regular, size n₂×n₂) ----------
   ZxR = Matrix{Float64}(undef, n₂, n₂)
   ZyR = Matrix{Float64}(undef, n₂, n₂)
   DJR = Matrix{Float64}(undef, n₂, n₂)
   FXR = Matrix{Float64}(undef, n₂, n₂)
   H = similar(DJR)

   tolbd = 1e-13

   # -----------------------------------------------------------------------
   # main loop over target points (interior first, then optional boundary)
   # -----------------------------------------------------------------------
   last = Ni + is_s_small * Nd

   VOL_NEAR = UInt8(1)
   VOL_CLOSE = UInt8(2)

   #Patch perspective
   @inbounds for k in 1:M
      # Precompute regular 2D Chebyshev quad on (zx,zy)
      mapxy_Dmap!(ZxR, ZyR, DJR, d, zx, zy, k)
      f!(FXR, ZxR, ZyR)
      @. FXR = FXR * DJR

      qin = dp.invgo[k]
      qbd = dp.invgo[M+k]

      # ----- Weightage of kth patch for ith point -----
      @inbounds for i in 1:last

         isbd = i > Ni

         # tgtpts = [x, y, ℓ, j]  
         # x1, x2: Gauss nodes for the row/col of the (interior) collocation grid,
         # or boundary Gauss node if i > M*Np
         if isbd
            # boundary point
            k₀ = cld(i - Ni, N)
            l = d.kd[k₀]
            j = i - Ni - (k₀ - 1) * N
            x1 = 1.0
            x2 = zp[j]
         else
            l = cld(i, Np)
            j = i - (l - 1) * Np
            q, r = divrem(j - 1, N)
            x1 = zp[r+1]
            x2 = zp[q+1]
         end

         if l == k

            # -------- singular part on patch l: 4 contributions ----------
            # Part 1
            @. d1 = (x1 + 1.0) * dz
            @. d2 = (x2 + 1.0) * dz
            @. y1 = x1 - d1
            @. y2 = x2 - d2
            fill_meshgrid!(t1, t2, y1, y2)

            diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, l)
            @. LOG = log(DIF)
            f!(FXY, Zx, Zy)

            @. FXY = LOG * FXY * DJ
            I1 = dot(fwm, FXY, fwm)

            # Part 2 (reuse d1,y1; update d2,y2,t1,t2)
            @. d2 = (x2 - 1.0) * dz
            @. y2 = x2 - d2
            fill_meshgrid!(t1, t2, y1, y2)

            diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, l)
            @. LOG = log(DIF)
            f!(FXY, Zx, Zy)

            @. FXY = LOG * FXY * DJ
            I2 = dot(fwm, FXY, fwm)

            # Part 3 (update d1,y1 and d2,y2)
            @. d1 = (x1 - 1.0) * dz
            @. d2 = (x2 + 1.0) * dz
            @. y1 = x1 - d1
            @. y2 = x2 - d2
            fill_meshgrid!(t1, t2, y1, y2)

            diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, l)
            @. LOG = log(DIF)
            f!(FXY, Zx, Zy)

            @. FXY = LOG * FXY * DJ
            I3 = dot(fwm, FXY, fwm)

            # Part 4 (update d2,y2)
            @. d2 = (x2 - 1.0) * dz
            @. y2 = x2 - d2
            fill_meshgrid!(t1, t2, y1, y2)
            diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, l)
            @. LOG = log(DIF)
            f!(FXY, Zx, Zy)

            @. FXY = LOG * FXY * DJ
            I4 = dot(fwm, FXY, fwm)

            # combine four parts
            b[i] += ((x1 + 1) * (x2 + 1) * I1 + (x1 + 1) * (1 - x2) * I2 +
                     (1 - x1) * (x2 + 1) * I3 + (1 - x1) * (1 - x2) * I4) / 4

         elseif dp.volmode[i, k] == VOL_CLOSE
            # preimage (alp1, alp2)

            if isbd
               alp1 = dp.invpts[1, qbd]
               alp2 = dp.invpts[2, qbd]
               qbd += 1

            else
               alp1 = dp.invpts[1, qin]
               alp2 = dp.invpts[2, qin]
               qin += 1

            end

            # w-inverses (only defined for |α| ≥ 1)
            A1 = alp1 <= -1 ? -winv(p, (-1 - alp1) / (1 - alp1)) : 0.0
            B1 = alp1 >= 1 ? -winv(p, (alp1 - 1) / (alp1 + 1)) : 0.0
            A2 = alp2 <= -1 ? -winv(p, (-1 - alp2) / (1 - alp2)) : 0.0
            B2 = alp2 >= 1 ? -winv(p, (alp2 - 1) / (alp2 + 1)) : 0.0

            if alp1 >= 1
               # d1 from alp1 + B1*z1  ; y1 = alp1 - d1
               #@. d1 = (alp1 + 1) * w(-B1 * z1)
               wfunc!(d1, p, z1; α=-B1, β=alp1 + 1)
               @. y1 = alp1 - d1
               if alp2 >= 1
                  #@. d2 = (alp2 + 1) * w(-B2 * z1)
                  wfunc!(d2, p, z1; α=-B2, β=alp2 + 1)
                  @. y2 = alp2 - d2
                  #@. wL = fw * dw(B1 * z1)
                  #@. wR = fw * dw(B2 * z1)
                  dwfunc!(wL, p, z1; α=B1)
                  dwfunc!(wR, p, z1; α=B2)
                  @. wL = fw * wL
                  @. wR = fw * wR

                  fill_meshgrid!(t1, t2, y1, y2)
                  diff_map!(DIF, Zx, Zy, DJ, d, alp1, alp2, t1, t2, d1, d2, k)
                  @. LOG = log(DIF)
                  f!(FXY, Zx, Zy)

                  @. FXY = LOG * FXY * DJ
                  I = dot(wL, FXY, wR)
                  b[i] += B1 * B2 * (alp1 + 1) * (alp2 + 1) * I / 4

               elseif -1 < alp2 && alp2 < 1
                  # I1
                  #@. d2 = (alp2 + 1) * w(-z1)
                  wfunc!(d2, p, z1; α=-1.0, β=alp2 + 1.0)
                  @. y2 = alp2 - d2
                  fill_meshgrid!(t1, t2, y1, y2)
                  diff_map!(DIF, Zx, Zy, DJ, d, alp1, alp2, t1, t2, d1, d2, k)
                  @. LOG = log(DIF)
                  f!(FXY, Zx, Zy)
                  #@. wL = fw * dw(B1 * z1)
                  #@. wR = fw * dw(z1)
                  dwfunc!(wL, p, z1; α=B1)
                  dwfunc!(wR, p, z1; α=1.0)
                  @. wL = fw * wL
                  @. wR = fw * wR

                  @. FXY = LOG * FXY * DJ
                  I1 = dot(wL, FXY, wR)
                  # I2
                  @. d2 = (alp2 - 1.0) * dz
                  @. y2 = alp2 - d2
                  fill_meshgrid!(t1, t2, y1, y2)
                  diff_map!(DIF, Zx, Zy, DJ, d, alp1, alp2, t1, t2, d1, d2, k)
                  @. LOG = log(DIF)
                  f!(FXY, Zx, Zy)
                  #@. wR = fw * dw(z2)
                  dwfunc!(wR, p, z2; α=1.0)
                  @. wR = fw * wR

                  @. FXY = LOG * FXY * DJ
                  I2 = dot(wL, FXY, wR)
                  b[i] += B1 * (alp1 + 1) * ((alp2 + 1) * I1 + (1 - alp2) * I2) / 4

               else # alp2 <= -1
                  #@. d2 = (alp2 - 1) * w(-A2 * z2)
                  wfunc!(d2, p, z2; α=-A2, β=alp2 - 1.0)
                  @. y2 = alp2 - d2
                  fill_meshgrid!(t1, t2, y1, y2)
                  diff_map!(DIF, Zx, Zy, DJ, d, alp1, alp2, t1, t2, d1, d2, k)
                  @. LOG = log(DIF)
                  f!(FXY, Zx, Zy)
                  #@. wL = fw * dw(B1 * z1)
                  #@. wR = fw * dw(A2 * z2)
                  dwfunc!(wL, p, z1; α=B1)
                  dwfunc!(wR, p, z2; α=A2)
                  @. wL = fw * wL
                  @. wR = fw * wR

                  @. FXY = LOG * FXY * DJ
                  I = dot(wL, FXY, wR)
                  b[i] += B1 * A2 * (alp1 + 1) * (1 - alp2) * I / 4
               end

            elseif -1 < alp1 && alp1 < 1
               # (vi) and (viii) cases — mirror of above with roles swapped
               # vi: alp2 ≥ 1
               if alp2 >= 1
                  #@. d1 = (alp1 + 1) * w(-z1)
                  wfunc!(d1, p, z1; α=-1.0, β=alp1 + 1.0)
                  @. y1 = alp1 - d1
                  #@. d2 = (alp2 + 1) * w(-B2 * z1)
                  wfunc!(d2, p, z1; α=-B2, β=alp2 + 1.0)
                  @. y2 = alp2 - d2
                  fill_meshgrid!(t1, t2, y1, y2)
                  diff_map!(DIF, Zx, Zy, DJ, d, alp1, alp2, t1, t2, d1, d2, k)
                  @. LOG = log(DIF)
                  f!(FXY, Zx, Zy)
                  #@. wL = fw * dw(z1)
                  #@. wR = fw * dw(B2 * z1)
                  dwfunc!(wL, p, z1; α=1.0)
                  dwfunc!(wR, p, z1; α=B2)
                  @. wL = fw * wL
                  @. wR = fw * wR

                  @. FXY = LOG * FXY * DJ
                  I1 = dot(wL, FXY, wR)

                  @. d1 = (alp1 - 1.0) * dz
                  @. y1 = alp1 - d1
                  #@. d2 = (alp2 + 1) * w(-B2 * z1)
                  wfunc!(d2, p, z1; α=-B2, β=alp2 + 1.0)
                  @. y2 = alp2 - d2
                  fill_meshgrid!(t1, t2, y1, y2)
                  diff_map!(DIF, Zx, Zy, DJ, d, alp1, alp2, t1, t2, d1, d2, k)
                  @. LOG = log(DIF)
                  f!(FXY, Zx, Zy)
                  #@. wL = fw * dw(z2)
                  dwfunc!(wL, p, z2; α=1.0)
                  @. wL = fw * wL

                  @. FXY = LOG * FXY * DJ
                  I2 = dot(wL, FXY, wR)
                  b[i] += B2 * (alp2 + 1) * ((alp1 + 1) * I1 + (1 - alp1) * I2) / 4

               elseif alp2 <= -1
                  #@. d1 = (alp1 + 1) * w(-z1)
                  wfunc!(d1, p, z1; α=-1.0, β=alp1 + 1.0)
                  @. y1 = alp1 - d1
                  #@. d2 = (alp2 - 1) * w(-A2 * z2)
                  wfunc!(d2, p, z2; α=-A2, β=alp2 - 1.0)
                  @. y2 = alp2 - d2
                  fill_meshgrid!(t1, t2, y1, y2)
                  diff_map!(DIF, Zx, Zy, DJ, d, alp1, alp2, t1, t2, d1, d2, k)
                  @. LOG = log(DIF)
                  f!(FXY, Zx, Zy)
                  #@. wL = fw * dw(z1)
                  #@. wR = fw * dw(A2 * z2)
                  dwfunc!(wL, p, z1; α=1.0)
                  dwfunc!(wR, p, z2; α=A2)
                  @. wL = fw * wL
                  @. wR = fw * wR

                  @. FXY = LOG * FXY * DJ
                  I1 = dot(wL, FXY, wR)

                  @. d1 = (alp1 - 1.0) * dz
                  @. y1 = alp1 - d1
                  #@. d2 = (alp2 - 1) * w(-A2 * z2)
                  wfunc!(d2, p, z2; α=-A2, β=alp2 - 1.0)
                  @. y2 = alp2 - d2
                  fill_meshgrid!(t1, t2, y1, y2)
                  diff_map!(DIF, Zx, Zy, DJ, d, alp1, alp2, t1, t2, d1, d2, k)
                  @. LOG = log(DIF)
                  f!(FXY, Zx, Zy)
                  #@. wL = fw * dw(z2)
                  dwfunc!(wL, p, z2; α=1.0)
                  @. wL = fw * wL

                  @. FXY = LOG * FXY * DJ
                  I2 = dot(wL, FXY, wR)
                  b[i] += A2 * (1 - alp2) * ((alp1 + 1) * I1 + (1 - alp1) * I2) / 4

               else
                  error("mapinv near-sing: unexpected alp2 in (-1,1)")
               end

            else  # alp1 <= -1
               #@. d1 = (alp1 - 1) * w(-A1 * z2)
               wfunc!(d1, p, z2; α=-A1, β=alp1 - 1.0)
               @. y1 = alp1 - d1
               if alp2 >= 1
                  #@. d2 = (alp2 + 1) * w(-B2 * z1)
                  wfunc!(d2, p, z1; α=-B2, β=alp2 + 1.0)
                  @. y2 = alp2 - d2
                  fill_meshgrid!(t1, t2, y1, y2)
                  diff_map!(DIF, Zx, Zy, DJ, d, alp1, alp2, t1, t2, d1, d2, k)
                  @. LOG = log(DIF)
                  f!(FXY, Zx, Zy)
                  #@. wL = fw * dw(A1 * z2)
                  #@. wR = fw * dw(B2 * z1)
                  dwfunc!(wL, p, z2; α=A1)
                  dwfunc!(wR, p, z1; α=B2)
                  @. wL = fw * wL
                  @. wR = fw * wR
                  @. FXY = LOG * FXY * DJ
                  I = dot(wL, FXY, wR)
                  b[i] += A1 * B2 * (1 - alp1) * (alp2 + 1) * I / 4

               elseif -1 < alp2 && alp2 < 1
                  # I1
                  #@. d2 = (alp2 + 1) * w(-z1)
                  wfunc!(d2, p, z1; α=-1.0, β=alp2 + 1.0)
                  @. y2 = alp2 - d2
                  fill_meshgrid!(t1, t2, y1, y2)
                  diff_map!(DIF, Zx, Zy, DJ, d, alp1, alp2, t1, t2, d1, d2, k)
                  @. LOG = log(DIF)
                  f!(FXY, Zx, Zy)
                  #@. wL = fw * dw(A1 * z2)
                  #@. wR = fw * dw(z1)
                  dwfunc!(wL, p, z2; α=A1)
                  dwfunc!(wR, p, z1; α=1.0)
                  @. wL = fw * wL
                  @. wR = fw * wR

                  @. FXY = LOG * FXY * DJ
                  I1 = dot(wL, FXY, wR)

                  # I2
                  @. d2 = (alp2 - 1.0) * dz
                  @. y2 = alp2 - d2
                  fill_meshgrid!(t1, t2, y1, y2)
                  diff_map!(DIF, Zx, Zy, DJ, d, alp1, alp2, t1, t2, d1, d2, k)
                  @. LOG = log(DIF)
                  f!(FXY, Zx, Zy)
                  #@. wR = fw * dw(z2)
                  dwfunc!(wR, p, z2; α=1.0)
                  @. wR = fw * wR

                  @. FXY = LOG * FXY * DJ
                  I2 = dot(wL, FXY, wR)
                  b[i] += A1 * (1 - alp1) * ((alp2 + 1) * I1 + (1 - alp2) * I2) / 4

               else # alp2 <= -1
                  #@. d2 = (alp2 - 1) * w(-A2 * z2)
                  wfunc!(d2, p, z2; α=-A2, β=alp2 - 1.0)
                  @. y2 = alp2 - d2
                  fill_meshgrid!(t1, t2, y1, y2)
                  diff_map!(DIF, Zx, Zy, DJ, d, alp1, alp2, t1, t2, d1, d2, k)
                  @. LOG = log(DIF)
                  f!(FXY, Zx, Zy)
                  #@. wL = fw * dw(A1 * z2)
                  #@. wR = fw * dw(A2 * z2)
                  dwfunc!(wL, p, z2; α=A1)
                  dwfunc!(wR, p, z2; α=A2)
                  @. wL = fw * wL
                  @. wR = fw * wR

                  @. FXY = LOG * FXY * DJ
                  I = dot(wL, FXY, wR)
                  b[i] += A1 * A2 * (1 - alp1) * (1 - alp2) * I / 4
               end
            end

         elseif dp.volmode[i, k] == VOL_NEAR

            xx = dp.tgtpts[1, i]
            yy = dp.tgtpts[2, i]

            if isbd
               x1 = dp.invpts[1, qbd]
               x2 = dp.invpts[2, qbd]
               qbd += 1

            else
               x1 = dp.invpts[1, qin]
               x2 = dp.invpts[2, qin]
               qin += 1

            end

            #Eight possible cases, 4 corners and 4 edges:
            if abs(x1 - 1.0) < tolbd

               if abs(x2 - 1.0) < tolbd
                  #Corner (x1,x2) = (1.0,1.0)

                  @. y1 = 1.0 - 2.0 * dz
                  @. y2 = 1.0 - 2.0 * dz
                  fill_meshgrid!(t1, t2, y1, y2)

                  mapxy_Dmap!(Zx, Zy, DJ, d, t1, t2, k)
                  f!(FXY, Zx, Zy)
                  @. LOG = 0.5 * log((xx - Zx)^2 + (yy - Zy)^2)
                  @. FXY = LOG * FXY * DJ
                  b[i] += dot(fwm, FXY, fwm)

               elseif abs(x2 + 1.0) < tolbd
                  #Corner (x1,x2) = (1.0,-1.0)

                  @. y1 = 1.0 - 2.0 * dz
                  @. y2 = -1.0 + 2.0 * dz
                  fill_meshgrid!(t1, t2, y1, y2)

                  mapxy_Dmap!(Zx, Zy, DJ, d, t1, t2, k)
                  f!(FXY, Zx, Zy)
                  @. LOG = 0.5 * log((xx - Zx)^2 + (yy - Zy)^2)
                  @. FXY = LOG * FXY * DJ
                  b[i] += dot(fwm, FXY, fwm)

               else
                  # -1.0 < x2 < 1.0, (1.0, x2)
                  # Part 1
                  @. y1 = 1.0 - 2.0 * dz
                  @. y2 = x2 - (x2 + 1.0) * dz
                  fill_meshgrid!(t1, t2, y1, y2)

                  mapxy_Dmap!(Zx, Zy, DJ, d, t1, t2, k)
                  f!(FXY, Zx, Zy)
                  @. LOG = 0.5 * log((xx - Zx)^2 + (yy - Zy)^2)
                  @. FXY = LOG * FXY * DJ
                  I1 = dot(fwm, FXY, fwm)

                  # Part 2 (reuse y1; update y2,t1,t2)
                  @. y2 = x2 - (x2 - 1.0) * dz
                  fill_meshgrid!(t1, t2, y1, y2)

                  mapxy_Dmap!(Zx, Zy, DJ, d, t1, t2, k)
                  f!(FXY, Zx, Zy)
                  @. LOG = 0.5 * log((xx - Zx)^2 + (yy - Zy)^2)
                  @. FXY = LOG * FXY * DJ
                  I2 = dot(fwm, FXY, fwm)

                  b[i] += ((x2 + 1) * I1 + (1 - x2) * I2) / 2

               end

            elseif abs(x1 + 1.0) < tolbd

               if abs(x2 - 1.0) < tolbd
                  #Corner (x1,x2) = (-1.0,1.0)

                  @. y1 = -1.0 + 2.0 * dz
                  @. y2 = 1.0 - 2.0 * dz
                  fill_meshgrid!(t1, t2, y1, y2)

                  mapxy_Dmap!(Zx, Zy, DJ, d, t1, t2, k)
                  f!(FXY, Zx, Zy)
                  @. LOG = 0.5 * log((xx - Zx)^2 + (yy - Zy)^2)
                  @. FXY = LOG * FXY * DJ
                  b[i] += dot(fwm, FXY, fwm)

               elseif abs(x2 + 1.0) < tolbd
                  #Corner (x1,x2) = (-1.0,-1.0)

                  @. y1 = -1.0 + 2.0 * dz
                  @. y2 = -1.0 + 2.0 * dz
                  fill_meshgrid!(t1, t2, y1, y2)

                  mapxy_Dmap!(Zx, Zy, DJ, d, t1, t2, k)
                  f!(FXY, Zx, Zy)
                  @. LOG = 0.5 * log((xx - Zx)^2 + (yy - Zy)^2)
                  @. FXY = LOG * FXY * DJ
                  b[i] += dot(fwm, FXY, fwm)

               else
                  # -1.0 < x2 < 1.0, (-1.0, x2)
                  # Part 3
                  @. y1 = -1.0 + 2.0 .* dz
                  @. y2 = x2 - (x2 + 1.0) .* dz
                  fill_meshgrid!(t1, t2, y1, y2)

                  mapxy_Dmap!(Zx, Zy, DJ, d, t1, t2, k)
                  f!(FXY, Zx, Zy)
                  @. LOG = 0.5 * log((xx - Zx)^2 + (yy - Zy)^2)
                  @. FXY = LOG * FXY * DJ
                  I3 = dot(fwm, FXY, fwm)

                  # Part 4 (reuse y1; update y2,t1,t2)
                  @. y2 = x2 - (x2 - 1.0) * dz
                  fill_meshgrid!(t1, t2, y1, y2)

                  mapxy_Dmap!(Zx, Zy, DJ, d, t1, t2, k)
                  f!(FXY, Zx, Zy)
                  @. LOG = 0.5 * log((xx - Zx)^2 + (yy - Zy)^2)
                  @. FXY = LOG * FXY * DJ
                  I4 = dot(fwm, FXY, fwm)

                  b[i] += ((x2 + 1) * I3 + (1 - x2) * I4) / 2

               end

            elseif abs(x2 - 1.0) < tolbd
               #Edge case -1.0 < x1 < 1.0, (x1,  1.0)

               # Part 1
               @. y1 = x1 - (x1 + 1.0) * dz
               @. y2 = 1.0 - 2.0 * dz
               fill_meshgrid!(t1, t2, y1, y2)

               mapxy_Dmap!(Zx, Zy, DJ, d, t1, t2, k)
               f!(FXY, Zx, Zy)
               @. LOG = 0.5 * log((xx - Zx)^2 + (yy - Zy)^2)
               @. FXY = LOG * FXY * DJ
               I1 = dot(fwm, FXY, fwm)

               # Part 3
               @. y1 = x1 - (x1 - 1.0) * dz
               fill_meshgrid!(t1, t2, y1, y2)

               mapxy_Dmap!(Zx, Zy, DJ, d, t1, t2, k)
               f!(FXY, Zx, Zy)
               @. LOG = 0.5 * log((xx - Zx)^2 + (yy - Zy)^2)
               @. FXY = LOG * FXY * DJ
               I3 = dot(fwm, FXY, fwm)

               # combine four parts
               b[i] += ((x1 + 1) * I1 + (1 - x1) * I3 ) / 2

            elseif abs(x2 + 1.0) < tolbd
               #Edge case -1.0 < x1 < 1.0, (x1, -1.0)

               # Part 2 
               @. y1 = x1 - (x1 + 1.0) * dz
               @. y2 = -1.0 + 2.0 * dz
               fill_meshgrid!(t1, t2, y1, y2)

               mapxy_Dmap!(Zx, Zy, DJ, d, t1, t2, k)
               f!(FXY, Zx, Zy)
               @. LOG = 0.5 * log((xx - Zx)^2 + (yy - Zy)^2)
               @. FXY = LOG * FXY * DJ
               I2 = dot(fwm, FXY, fwm)

               # Part 4 
               @. y1 = x1 - (x1 - 1.0) * dz
               fill_meshgrid!(t1, t2, y1, y2)

               mapxy_Dmap!(Zx, Zy, DJ, d, t1, t2, k)
               f!(FXY, Zx, Zy)
               @. LOG = 0.5 * log((xx - Zx)^2 + (yy - Zy)^2)
               @. FXY = LOG * FXY * DJ
               I4 = dot(fwm, FXY, fwm)

               # combine four parts
               b[i] += ((x1 + 1) * I2 + (1 - x1) * I4) / 2

            else
               error("VOL_NEAR projection is not on patch boundary: row=$i, patch=$k, x1=$x1, x2=$x2")

            end

         else

            xx = dp.tgtpts[1, i]
            yy = dp.tgtpts[2, i]
            # diff = sqrt((x-ϕx)^2+(y-ϕy)^2);  pt.x/pt.y are physical coords
            @. H = 0.5 * log((xx - ZxR)^2 + (yy - ZyR)^2) * FXR
            b[i] += dot(fw₂, H, fw₂)
         end
      end
   end

   b ./= (2π)

   return b

end

#------------------------------------------------------------------------
struct NSRWork
   z::Vector{Float64}
   z1::Vector{Float64}
   z2::Vector{Float64}
   dz::Vector{Float64}
   fw::Vector{Float64}
   fwm::Vector{Float64}

   d1::Vector{Float64}
   d2::Vector{Float64}
   y1::Vector{Float64}
   y2::Vector{Float64}
   y1tmp::Vector{Float64}
   wdf::Vector{Float64}
   dw::Vector{Float64}
   dfy::Vector{Float64}
   dfys::Vector{Float64}
   wr::Vector{Float64}

   t1::Matrix{Float64}
   t2::Matrix{Float64}
   DJ::Matrix{Float64}
   DIF::Matrix{Float64}
   Zx::Matrix{Float64}
   Zy::Matrix{Float64}

   FX::Matrix{Float64}
   KS::Matrix{Float64}
   βS::Matrix{Float64}
   HS::Matrix{Float64}

   αₛ::Float64
   f!::Function
end

function NSRWork(n::Int, p::Int, αₛ::Float64, f!::Function)
   z = Vector{Float64}(undef, n)
   z1 = similar(z)
   z2 = similar(z)
   dz = similar(z)

   fw = getF1W(n)

   @inbounds for i in 1:n
      c = π * (2i - 1) / (2n)
      z[i] = cos(c)
      z1[i] = cos(c / 2)^2
      z2[i] = sin(c / 2)^2
   end

   wfunc!(dz, p, z2; α = -1.0)

   fwm = similar(fw)
   dwfunc!(fwm, p, z2)
   @. fwm = fw * fwm

   d1 = similar(z)
   d2 = similar(z)
   y1 = similar(z)
   y2 = similar(z)
   y1tmp = similar(z)
   wdf = similar(z)
   dw = similar(z)
   dfy = similar(z)
   dfys = similar(z)
   wr = similar(z)

   t1 = Matrix{Float64}(undef, n, n)
   t2 = Matrix{Float64}(undef, n, n)
   DJ = Matrix{Float64}(undef, n, n)
   DIF = Matrix{Float64}(undef, n, n)
   Zx = Matrix{Float64}(undef, n, n)
   Zy = Matrix{Float64}(undef, n, n)

   FX = Matrix{Float64}(undef, n, n)
   KS = Matrix{Float64}(undef, n, n)
   βS = Matrix{Float64}(undef, n, n)
   HS = Matrix{Float64}(undef, n, n)

   return NSRWork(
      z, z1, z2, dz, fw, fwm,
      d1, d2, y1, y2, y1tmp, wdf, dw, dfy, 
      dfys, wr, t1, t2, DJ, DIF, Zx, Zy,
      FX, KS, βS, HS, αₛ, f!
   )
end

@inline function eval_close_piece_Rs(
   work, d::abstractdomain,
   k::Int, α₁::Float64, α₂::Float64,
   y1::AbstractVector{Float64},
   y2::AbstractVector{Float64},
   d1::AbstractVector{Float64},
   d2::AbstractVector{Float64},
   dfy::AbstractVector{Float64},
   dfys::AbstractVector{Float64},
   wleft::AbstractVector{Float64},
   wright::AbstractVector{Float64},
   n::Int, s::Float64)::Float64

   fill_meshgrid!(work.t1, work.t2, y1, y2)

   diff_map!(work.DIF, work.Zx, work.Zy, work.DJ,
      d, α₁, α₂, work.t1, work.t2, d1, d2, k)

   work.f!(work.FX, work.Zx, work.Zy)

   @inbounds for j in 1:n
      @inbounds for i in 1:n
         r2 = work.DIF[i, j]^2

         dval = dfy[i]
         ds = dfys[i]

         logd = log(dval)
         logr2 = log(r2)

         βt1 = logd - logr2
         βt2 = logd

         βt3 = phi2_slr(βt1; α=s)
         βt4 = phi2_slr(βt2; α=s)

         work.βS[i, j] = βt1^2 * βt3 - βt2^2 * βt4
         work.HS[i, j] = ds * expm1(-s * logr2) / s
      end
   end

   @. work.KS = (work.βS / (4.0 * pi) - work.αₛ * work.HS) * work.FX * work.DJ

   return dot(wleft, work.KS, wright)
end

function NSclose_Rs!(work, d::abstractdomain, k::Int,
   α₁::Float64, α₂::Float64, s::Float64,
   p::Int, n::Int, knbd)::Float64

   (; z, z1, z2, fw, d1, d2, y1, y2,
      y1tmp, wdf, dw, dfy, dfys, wr) = work

   if 1.0 <= α₁
      B1 = -winv(p, (α₁ - 1.0) / (α₁ + 1.0))
   elseif α₁ <= -1.0
      A1 = -winv(p, (-1.0 - α₁) / (1.0 - α₁))
   end

   if 1.0 <= α₂
      B2 = -winv(p, (α₂ - 1.0) / (α₂ + 1.0))
   elseif α₂ <= -1.0
      A2 = -winv(p, (-1.0 - α₂) / (1.0 - α₂))
   end

   if k in knbd

      if 1.0 <= α₁

         wfunc!(d1, p, z1; α=-B1, β=α₁ + 1.0)
         @. y1 = α₁ - d1
         @. y1tmp = 1.0 - y1

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)
         #Distance factor is not included in weight
         #It is purely part of the knerel now!
         dwfunc!(dw, p, z1; α=B1)
         @. wdf = fw * dw

         if 1.0 <= α₂

            wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z1; α=B2)
            @. wr = fw * dw

            I = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I *= B1 * B2 * (α₁ + 1.0) * (α₂ + 1.0) / 4.0

         elseif -1.0 < α₂ && α₂ < 1.0

            # Lower part in α₂.
            wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z1)
            @. wr = fw * dw

            I₁ = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            # Upper part in α₂.
            wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z2)
            @. wr = fw * dw

            I₂ = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I = B1 * (α₁ + 1.0) *
                   ((α₂ + 1.0) * I₁ + (1.0 - α₂) * I₂) / 4.0

         elseif α₂ <= -1.0

            wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z2; α=A2)
            @. wr = fw * dw

            I = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I *= B1 * A2 * (α₁ + 1.0) * (1.0 - α₂) / 4.0
         end

      elseif -1.0 < α₁ && α₁ < 1.0

         if 1.0 <= α₂

            # Left part in α₁.
            wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
            @. y1 = α₁ - d1
            @. y1tmp = 1.0 - y1

            dfunc!(dfy, d, k, y1tmp)

            dfunc!(dfys, d, k, y1tmp, s)

            dwfunc!(dw, p, z1)
            @. wdf = fw * dw

            wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z1; α=B2)
            @. wr = fw * dw

            I₁ = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            # Right part in α₁.
            wfunc!(d1, p, z2; α=-1.0, β=α₁ - 1.0)
            @. y1 = α₁ - d1
            @. y1tmp = 1.0 - y1

            dfunc!(dfy, d, k, y1tmp)

            dfunc!(dfys, d, k, y1tmp, s)

            dwfunc!(dw, p, z2)
            @. wdf = fw * dw

            I₂ = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I = B2 * (α₂ + 1.0) *
                   ((α₁ + 1.0) * I₁ + (1.0 - α₁) * I₂) / 4.0

         elseif α₂ <= -1.0

            # Left part in α₁.
            wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
            @. y1 = α₁ - d1
            @. y1tmp = 1.0 - y1

            dfunc!(dfy, d, k, y1tmp)

            dfunc!(dfys, d, k, y1tmp, s)

            dwfunc!(dw, p, z1)
            @. wdf = fw * dw

            wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z2; α=A2)
            @. wr = fw * dw

            I₁ = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            # Right part in α₁.
            wfunc!(d1, p, z2; α=-1.0, β=α₁ - 1.0)
            @. y1 = α₁ - d1
            @. y1tmp = 1.0 - y1

            dfunc!(dfy, d, k, y1tmp)

            dfunc!(dfys, d, k, y1tmp, s)

            dwfunc!(dw, p, z2)
            @. wdf = fw * dw

            I₂ = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I = A2 * (1.0 - α₂) *
                   ((α₁ + 1.0) * I₁ + (1.0 - α₁) * I₂) / 4.0

         else
            error("Fatal error! mapinv not working correctly.")
         end

      elseif α₁ <= -1.0

         wfunc!(d1, p, z2; α=-A1, β=α₁ - 1.0)
         @. y1 = α₁ - d1
         @. y1tmp = 1.0 - y1

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)

         dwfunc!(dw, p, z2; α=A1)
         @. wdf = fw * dw

         if 1.0 <= α₂

            wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z1; α=B2)
            @. wr = fw * dw

            I = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I *= A1 * B2 * (1.0 - α₁) * (α₂ + 1.0) / 4.0

         elseif -1.0 < α₂ && α₂ < 1.0

            # Lower part in α₂.
            wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z1)
            @. wr = fw * dw

            I₁ = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            # Upper part in α₂.
            wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z2)
            @. wr = fw * dw

            I₂ = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I = A1 * (1.0 - α₁) *
                   ((α₂ + 1.0) * I₁ + (1.0 - α₂) * I₂) / 4.0

         elseif α₂ <= -1.0

            wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z2; α=A2)
            @. wr = fw * dw

            I = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I *= A1 * A2 * (1.0 - α₁) * (1.0 - α₂) / 4.0
         end
      end

   else
      # Boundary-touching volume patch.

      if 1.0 <= α₁

         # α₁ = 1 boundary-side branch.
         wfunc!(d1, p, z1; α=-1.0, β=2.0)
         @. y1 = 1.0 - d1

         dfunc!(dfy, d, k, d1)

         dfunc!(dfys, d, k, d1, s)

         dwfunc!(dw, p, z1)
         @. wdf = fw * dw

         if 1.0 <= α₂

            wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z1; α=B2)
            @. wr = fw * dw

            I = eval_close_piece_Rs(work, d, k, 1.0, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I *= B2 * (α₂ + 1.0) / 2.0

         elseif -1.0 < α₂ && α₂ < 1.0

            # Lower part in α₂.
            wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z1)
            @. wr = fw * dw

            I₁ = eval_close_piece_Rs(work, d, k, 1.0, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            # Upper part in α₂.
            wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z2)
            @. wr = fw * dw

            I₂ = eval_close_piece_Rs(work, d, k, 1.0, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I = ((α₂ + 1.0) * I₁ + (1.0 - α₂) * I₂) / 2.0

         elseif α₂ <= -1.0

            wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z2; α=A2)
            @. wr = fw * dw

            I = eval_close_piece_Rs(work, d, k, 1.0, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I *= A2 * (1.0 - α₂) / 2.0
         end

      elseif -1.0 < α₁ && α₁ < 1.0

         if 1.0 <= α₂

            # Left part in α₁: regular dfunc treatment.
            wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
            @. y1 = α₁ - d1
            @. y1tmp = 1.0 - y1

            dfunc!(dfy, d, k, y1tmp)

            dfunc!(dfys, d, k, y1tmp, s)

            dwfunc!(dw, p, z1)
            @. wdf = fw * dw

            wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z1; α=B2)
            @. wr = fw * dw

            I₁ = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            # Right part in α₁: boundary side, still d^s.
            wfunc!(d1, p, z; α=1.0, β=(α₁ - 1.0) / 2.0)
            @. y1 = α₁ - d1

            wfunc!(y1tmp, p, z; α=-1.0, β=(1.0 - α₁) / 2.0)

            dfunc!(dfy, d, k, y1tmp)

            dfunc!(dfys, d, k, y1tmp, s)

            dwfunc!(dw, p, z)
            @. wdf = fw * dw

            I₂ = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I = B2 * (α₂ + 1.0) * ((α₁ + 1.0) * I₁ + (1.0 - α₁) * I₂) / 4.0

         elseif α₂ <= -1.0

            # Left part in α₁: regular dfunc treatment.
            wfunc!(d1, p, z1; α=-1.0, β=α₁ + 1.0)
            @. y1 = α₁ - d1
            @. y1tmp = 1.0 - y1

            dfunc!(dfy, d, k, y1tmp)

            dfunc!(dfys, d, k, y1tmp, s)

            dwfunc!(dw, p, z1)
            @. wdf = fw * dw

            wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z2; α=A2)
            @. wr = fw * dw

            I₁ = eval_close_piece_Rs( work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            # Right part in α₁: boundary side, still d^s.
            wfunc!(d1, p, z; α=1.0, β=(α₁ - 1.0) / 2.0)
            @. y1 = α₁ - d1

            wfunc!(y1tmp, p, z; α=-1.0, β=(1.0 - α₁) / 2.0)

            dfunc!(dfy, d, k, y1tmp)

            dfunc!(dfys, d, k, y1tmp, s)

            dwfunc!(dw, p, z)
            @. wdf = fw * dw

            I₂ = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I = A2 * (1.0 - α₂) * ((α₁ + 1.0) * I₁ + (1.0 - α₁) * I₂) / 4.0

         else
            error("Fatal error! mapinv not working correctly.")
         end

      elseif α₁ <= -1.0

         # Special boundary-patch left-side formula from precompsL.
         A1 = 1.0 - winv(p, 2.0 * (-1.0 - α₁) / (1.0 - α₁))

         wfunc!(d1, p, z2; α=-A1, β=(α₁ - 1.0) / 2.0, γ=1.0)
         @. y1 = α₁ - d1

         wfunc!(y1tmp, p, z2; α=A1, β=(1.0 - α₁) / 2.0, γ=-1.0)

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)

         dwfunc!(dw, p, z2; α=-A1, β=1.0, γ=1.0)
         @. wdf = fw * dw

         if 1.0 <= α₂

            wfunc!(d2, p, z1; α=-B2, β=α₂ + 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z1; α=B2)
            @. wr = fw * dw

            I = eval_close_piece_Rs( work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I *= A1 * B2 * (1.0 - α₁) * (α₂ + 1.0) / 8.0

         elseif -1.0 < α₂ && α₂ < 1.0

            # Lower part in α₂.
            wfunc!(d2, p, z1; α=-1.0, β=α₂ + 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z1)
            @. wr = fw * dw

            I₁ = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            # Upper part in α₂.
            wfunc!(d2, p, z2; α=-1.0, β=α₂ - 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z2)
            @. wr = fw * dw

            I₂ = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I = A1 * (1.0 - α₁) * ((α₂ + 1.0) * I₁ + (1.0 - α₂) * I₂) / 8.0

         elseif α₂ <= -1.0

            wfunc!(d2, p, z2; α=-A2, β=α₂ - 1.0)
            @. y2 = α₂ - d2

            dwfunc!(dw, p, z2; α=A2)
            @. wr = fw * dw

            I = eval_close_piece_Rs(work, d, k, α₁, α₂,
               y1, y2, d1, d2, dfy, dfys, wdf, wr, n, s)

            I *= A1 * A2 * (1.0 - α₁) * (1.0 - α₂) / 8.0
         end
      end
   end

   return I
end

@inline function eval_piece_Rs(
   work, d::abstractdomain, 
   k::Int, xt::Float64, yt::Float64,
   y1::AbstractVector{Float64},
   y2::AbstractVector{Float64},
   dfy::AbstractVector{Float64},
   dfys::AbstractVector{Float64},
   wleft::AbstractVector{Float64},
   wright::AbstractVector{Float64},
   n::Int, s::Float64)::Float64

   fill_meshgrid!(work.t1, work.t2, y1, y2)

   mapxy_Dmap!(work.Zx, work.Zy, work.DJ, d, work.t1, work.t2, k)

   work.f!(work.FX, work.Zx, work.Zy)

   @inbounds for j in 1:n
      @inbounds for i in 1:n
         r2 = (xt - work.Zx[i, j])^2 + (yt - work.Zy[i, j])^2

         dval = dfy[i]       # raw d(y)
         ds = dfys[i]      # d(y)^s

         logd = log(dval)
         logr2 = log(r2)

         βt1 = logd - logr2
         βt2 = logd

         βt3 = phi2_slr(βt1; α=s)
         βt4 = phi2_slr(βt2; α=s)

         work.βS[i, j] = βt1^2 * βt3 - βt2^2 * βt4
         work.HS[i, j] = ds * expm1(-s * logr2) / s
      end
   end

   @. work.KS = ( work.βS / (4 * pi) - work.αₛ * work.HS ) * work.FX * work.DJ

   return dot(wleft, work.KS ,wright)
end

function NSnear_Rs!(work, d::abstractdomain,
   k::Int, xt::Float64, yt::Float64,
   u0::Float64, v0::Float64, s::Float64,
   pn::Int, n::Int, knbd;
   tolbd::Float64=1e-12)::Float64

   (; z, fw, d1, y1, y2, y1tmp, wdf, dw, dfy, dfys, dz, fwm) = work

   if k in knbd
      # Non-boundary volume patch.
      if abs(u0 - 1.0) < tolbd
         # Projection on right wall: (1, v0).

         @. y1 = 1.0 - 2.0 * dz
         @. y1tmp = 2.0 * dz

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)


         if abs(v0 - 1.0) < tolbd

            @. y2 = 1.0 - 2.0 * dz

            I = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

         elseif abs(v0 + 1.0) < tolbd

            @. y2 = -1.0 + 2.0 * dz

            I = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

         else

            # Lower part in v.
            @. y2 = v0 - (v0 + 1.0) * dz

            I₁ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

            # Upper part in v.
            @. y2 = v0 - (v0 - 1.0) * dz

            I₂ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

            I = ((v0 + 1.0) * I₁ + (1.0 - v0) * I₂) / 2.0
         end

      elseif abs(u0 + 1.0) < tolbd
         # Projection on left wall: (-1, v0).

         @. y1 = -1.0 + 2.0 * dz
         @. y1tmp = 2.0 * (1.0 - dz)

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)

         if abs(v0 - 1.0) < tolbd

            @. y2 = 1.0 - 2.0 * dz

            I = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

         elseif abs(v0 + 1.0) < tolbd

            @. y2 = -1.0 + 2.0 * dz

            I = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

         else

            # Lower part in v.
            @. y2 = v0 - (v0 + 1.0) * dz

            I₁ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

            # Upper part in v.
            @. y2 = v0 - (v0 - 1.0) * dz

            I₂ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

            I = ((v0 + 1.0) * I₁ + (1.0 - v0) * I₂) / 2.0
         end

      elseif abs(v0 - 1.0) < tolbd
         # Projection on top wall: (u0, 1).

         @. y2 = 1.0 - 2.0 * dz

         # Left part in u.
         @. y1 = u0 - (u0 + 1.0) * dz
         @. y1tmp = 1.0 - y1

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)

         I₁ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

         # Right part in u.
         @. y1 = u0 - (u0 - 1.0) * dz
         @. y1tmp = 1.0 - y1

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)

         I₂ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

         I = ((u0 + 1.0) * I₁ + (1.0 - u0) * I₂) / 2.0

      elseif abs(v0 + 1.0) < tolbd
         # Projection on bottom wall: (u0, -1).

         @. y2 = -1.0 + 2.0 * dz

         # Left part in u.
         @. y1 = u0 - (u0 + 1.0) * dz
         @. y1tmp = 1.0 - y1

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)

         I₁ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

         # Right part in u.
         @. y1 = u0 - (u0 - 1.0) * dz
         @. y1tmp = 1.0 - y1

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)

         I₂ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

         I = ((u0 + 1.0) * I₁ + (1.0 - u0) * I₂) / 2.0

      else
         error("VOL_NEAR projection is not on patch boundary: patch=$k, u0=$u0, v0=$v0")
      end

   else
      # Boundary-touching patch.
      if abs(u0 - 1.0) < tolbd
         # Projection on right wall: (1, v0).
         # Projection side is also the physical boundary side.

         @. y1 = 1.0 - 2.0 * dz
         @. y1tmp = 2.0 * dz

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)

         if abs(v0 - 1.0) < tolbd

            @. y2 = 1.0 - 2.0 * dz

            I = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

         elseif abs(v0 + 1.0) < tolbd

            @. y2 = -1.0 + 2.0 * dz

            I = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

         else

            # Lower part in v.
            @. y2 = v0 - (v0 + 1.0) * dz

            I₁ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

            # Upper part in v.
            @. y2 = v0 - (v0 - 1.0) * dz

            I₂ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

            I = ((v0 + 1.0) * I₁ + (1.0 - v0) * I₂) / 2.0
         end

      elseif abs(u0 + 1.0) < tolbd
         # Projection on left wall: (-1, v0).
         # Need smoothing for left projection side and right boundary side.

         wfunc!(y1, pn, z)
         @. y1 = y1 - 1.0

         wfunc!(y1tmp, pn, z; α=-1.0)   # y1tmp = w(-z) = 1 - y1

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)

         dwfunc!(dw, pn, z)
         @. wdf = fw * dw

         if abs(v0 - 1.0) < tolbd

            @. y2 = 1.0 - 2.0 * dz

            I = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, wdf, fwm, n, s)

         elseif abs(v0 + 1.0) < tolbd

            @. y2 = -1.0 + 2.0 * dz

            I = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, wdf, fwm, n, s)

         else

            # Lower part in v.
            @. y2 = v0 - (v0 + 1.0) * dz

            I₁ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, wdf, fwm, n, s)

            # Upper part in v.
            @. y2 = v0 - (v0 - 1.0) * dz

            I₂ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, wdf, fwm, n, s)

            I = ((v0 + 1.0) * I₁ + (1.0 - v0) * I₂) / 2.0
         end

      elseif abs(v0 - 1.0) < tolbd
         # Projection on top wall: (u0, 1).
         # Split u into regular-left and boundary-right pieces.

         @. y2 = 1.0 - 2.0 * dz

         # Left part in u: regular dfunc.
         @. y1 = u0 - (u0 + 1.0) * dz
         @. y1tmp = 1.0 - y1
         
         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)

         I₁ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

         # Right part in u: boundary compensated rule.
         wfunc!(d1, pn, z)
         @. y1 = (1.0 - u0) * d1 / 2.0 + u0

         wfunc!(y1tmp, pn, z; α=-1.0, β=(1.0 - u0) / 2.0)
         
         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)

         dwfunc!(dw, pn, z)
         @. wdf = fw * dw

         I₂ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, wdf, fwm, n, s)

         I = ((u0 + 1.0) * I₁ + (1.0 - u0) * I₂) / 2.0

      elseif abs(v0 + 1.0) < tolbd
         # Projection on bottom wall: (u0, -1).
         # Split u into regular-left and boundary-right pieces.

         @. y2 = -1.0 + 2.0 * dz

         # Left part in u: regular dfunc.
         @. y1 = u0 - (u0 + 1.0) * dz
         @. y1tmp = 1.0 - y1

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)

         I₁ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, fwm, fwm, n, s)

         # Right part in u: boundary compensated rule.
         wfunc!(d1, pn, z)
         @. y1 = (1.0 - u0) * d1 / 2.0 + u0

         wfunc!(y1tmp, pn, z; α=-1.0, β=(1.0 - u0) / 2.0)

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)

         dwfunc!(dw, pn, z)
         @. wdf = fw * dw

         I₂ = eval_piece_Rs(work, d, k, xt, yt,
               y1, y2, dfy, dfys, wdf, fwm, n, s)

         I = ((u0 + 1.0) * I₁ + (1.0 - u0) * I₂) / 2.0

      else
         error("VOL_NEAR projection is not on patch boundary: patch=$k, u0=$u0, v0=$v0")
      end
   end

   return I
end

function phi2_slr(z::Float64; α::Float64=1.0)::Float64

   z = α * z
   az = abs(z)
   out = 0.0

   if az <= 5e-3

      term = 0.5
      out = term

      @inbounds for k in 0:6
         term *= z / (k + 3)
         out += term
      end

   elseif az <= 0.5
      term = 0.5
      out = term

      @inbounds for k in 0:14
         term *= z / (k + 3)
         out += term
      end

   else
      out = (expm1(z) - z) / (z * z)
   end

   return out
end

#For s<0.5
function bvec_small(d::abstractdomain, dp::domprop, s::Float64, f!::Function
   ; n::Int=64)::Vector{Float64}

   #------------------- Computing the αₛ constant -------------------
   EULER_GAMMA = 0.5772156649015328606
   ZETA3 = 1.2020569031595942854
   ZETA5 = 1.0369277551433699263
   ZETA7 = 1.0083492773819228268

   Eₛ = if s <= 0.02
      2 * s * ((log(2.0) - EULER_GAMMA) - ZETA3 * s^2 / 3 -
               ZETA5 * s^4 / 5 - ZETA7 * s^6 / 7)
   else
      2 * s * log(2.0) + loggamma(1 + s) - loggamma(1 - s)
   end

   αₛ = -expm1(Eₛ) / (4 * π * s)

   # The polynomial change-of-variables with parameter p
   p = 4

   # n₂ is defined as the number of nodes for integration of regular integrals
   # fw₂is defined as the Fejer 1st weights for corresponding nodes
   nr = n
   fwr = getF1W(nr)

   #bookkeeping
   M = d.Npat
   Mbd = length(d.kd)
   N = dp.N
   Np = N^2
   Ni = M * Np
   Nd = Mbd * N

   Nt = Ni + Nd

   b = zeros(Float64, Ni + Nd)

   bdloc = Dict{Int,Int}()
   sizehint!(bdloc, Mbd)

   @inbounds for ll in 1:Mbd
      bdloc[d.kd[ll]] = ll
   end

   # ------------ Reg Integration ------------
   # Store regular points
   zr = Vector{Float64}(undef, nr)
   zx = Matrix{Float64}(undef, nr, nr)
   zy = Matrix{Float64}(undef, nr, nr)
   zt = Vector{Float64}(undef, nr)

   @inbounds for i in 1:nr
      ph = π * (2 * i - 1) / (2 * nr)
      cs = cos(ph)
      ss = 2.0 * sin(ph / 2.0)^2 #1-cs
      zr[i] = cs
      zt[i] = ss   # 1-zx
      @inbounds for j in 1:nr
         # "meshgrid(zr)" :-> zx, zy of size (nr, nr)
         zy[j, i] = cs   # rows are zr
         zx[i, j] = cs   # columns are zr
      end
   end

   Zx = Matrix{Float64}(undef, nr, nr)
   Zy = Matrix{Float64}(undef, nr, nr)
   DJ = Matrix{Float64}(undef, nr, nr)
   FX = Matrix{Float64}(undef, nr, nr)
   HS = Matrix{Float64}(undef, nr, nr)
   βS = Matrix{Float64}(undef, nr, nr)

   KS = Matrix{Float64}(undef, nr, nr)
   Df = Vector{Float64}(undef, nr)
   Dfs= Vector{Float64}(undef, nr)

   # ------------ Bd Integration ------------
   # For boundary patches, which have singularity on boundary
   nbd = 2*n
   fwbd = getF1W(nbd)

   # z1 = cospi((2*(1:nbd)-1)/(4*nbd)).^2
   z1 = Vector{Float64}(undef, nbd)
   # We need yvec = 1 - 2*w(-z1).
   wz1 = similar(z1)
   # [zy2, zx2] = meshgrid(zr, 1 - 2*w(-z1));
   zy2 = Matrix{Float64}(undef, nbd, nr)
   zx2 = Matrix{Float64}(undef, nbd, nr)
   zt2 = Vector{Float64}(undef, nbd)
   mfw = Vector{Float64}(undef, nbd)

   @inbounds for i in 1:nbd
      z1[i] = cos(π * (2 * i - 1) / (4 * nbd))^2
   end

   dwfunc!(mfw, p, z1)

   @. mfw = fwbd * mfw

   wfunc!(wz1, p, z1; α=-1.0)

   @inbounds for i in 1:nbd 
      yv = 2 * wz1[i]
      zt2[i] = yv
      @inbounds for j in 1:nr
         zy2[i, j] = zr[j]     # each row of zy2 is zr
         zx2[i, j] = 1 - yv   # each column of zx2 is 1-yv
      end
   end

   Zx₂ = Matrix{Float64}(undef, nbd, nr)
   Zy₂ = Matrix{Float64}(undef, nbd, nr)
   DJ₂ = Matrix{Float64}(undef, nbd, nr)
   FX₂ = Matrix{Float64}(undef, nbd, nr)
   HS₂ = Matrix{Float64}(undef, nbd, nr)
   βS₂ = Matrix{Float64}(undef, nbd, nr)

   KS₂ = Matrix{Float64}(undef, nbd, nr)
   Df₂ = Vector{Float64}(undef, nbd)
   Df₂s= Vector{Float64}(undef, nbd)

   #Singular or near singular integration 
   ns= 2*n

   # Near/close work.
   work = NSRWork(ns, p, αₛ, f!)

   # -----------------------------------------------------------------------
   # main loop over target points (interior first, then optional boundary)
   # -----------------------------------------------------------------------

   VOL_NEAR = UInt8(1)
   VOL_CLOSE = UInt8(2)

   # knbd are patches which are not the boundary patches
   # and d.kd are patches are touching the boundary
   knbd = setdiff(collect(1:M), d.kd)

   #Just do NS and regular integrals first
   #Patch perspective
   @inbounds for k in 1:M
      
      qin = dp.invgo[k]
      qbd = dp.invgo[M+k]

      kb = get(bdloc, k, 0)
      isbdflag = kb != 0

      if isbdflag
         mapxy_Dmap!(Zx₂, Zy₂, DJ₂, d, zx2, zy2, k)

         f!(FX₂, Zx₂, Zy₂)

         dfunc!(Df₂, d, k, zt2)

         dfunc!(Df₂s, d, k, zt2, s)

      else
         mapxy_Dmap!(Zx, Zy, DJ, d, zx, zy, k)

         f!(FX, Zx, Zy)

         dfunc!(Df, d, k, zt)

         dfunc!(Dfs, d, k, zt, s)
      end

      # ----- Weightage of kth patch for ith point -----
      @inbounds for i in 1:Nt

         isbd = i > Ni

         # tgtpts = [x, y, l, j]  
         # x1, x2: nodes for the row/col of the (interior) collocation grid,
         # or boundary node if i > M*Np
         if isbd
            # boundary point
            k₀ = cld(i - Ni, N)
            l = d.kd[k₀]
         else
            l = cld(i, Np)
         end

         l == k && continue

         if dp.volmode[i, k] == VOL_CLOSE
            # preimage (α₁, α₂)

            if isbd
               α₁ = dp.invpts[1, qbd]
               α₂ = dp.invpts[2, qbd]
               qbd += 1

            else
               α₁ = dp.invpts[1, qin]
               α₂ = dp.invpts[2, qin]
               qin += 1

            end

            b[i] += NSclose_Rs!(work, d, k, α₁, α₂, s, p, ns, knbd)

         elseif dp.volmode[i, k] == VOL_NEAR

            xx = dp.tgtpts[1, i]
            yy = dp.tgtpts[2, i]

            if isbd
               u0 = dp.invpts[1, qbd]
               v0 = dp.invpts[2, qbd]
               qbd += 1

            else
               u0 = dp.invpts[1, qin]
               v0 = dp.invpts[2, qin]
               qin += 1

            end

            b[i] += NSnear_Rs!(work, d, k, xx, yy, u0, v0, s, p, ns, knbd)

         else

            xx = dp.tgtpts[1, i]
            yy = dp.tgtpts[2, i]

            if isbdflag

               @. KS₂ = (xx - Zx₂)^2 + (yy - Zy₂)^2

               @inbounds for jj in 1:nr
                  @inbounds for ii in 1:nbd
                     r2 = KS₂[ii, jj]

                     dval = Df₂[ii]          # raw distance d(y)
                     logd = log(dval)
                     logr2 = log(r2)

                     βt1 = logd - logr2     # log(d(y) / |x-y|^2)
                     βt2 = logd             # log(d(y))

                     βt3 = phi2_slr(βt1; α=s)
                     βt4 = phi2_slr(βt2; α=s)

                     βS₂[ii, jj] = βt1^2 * βt3 - βt2^2 * βt4

                     HS₂[ii, jj] = Df₂s[ii] * expm1(-s * logr2) / s
                  end
               end

               @. KS₂ = (βS₂ / (4.0 * pi) - αₛ * HS₂) * FX₂ * DJ₂

               b[i] += dot(mfw, KS₂, fwr)

            else

               @. KS = (xx - Zx)^2 + (yy - Zy)^2

               @inbounds for jj in 1:nr
                  @inbounds for ii in 1:nr
                     r2 = KS[ii, jj]

                     dval = Df[ii]           # raw distance d(y)
                     logd = log(dval)
                     logr2 = log(r2)

                     βt1 = logd - logr2
                     βt2 = logd

                     βt3 = phi2_slr(βt1; α=s)
                     βt4 = phi2_slr(βt2; α=s)

                     βS[ii, jj] = βt1^2 * βt3 - βt2^2 * βt4

                     HS[ii, jj] = Dfs[ii] * expm1(-s * logr2) / s
                  end
               end

               @. KS = (βS / (4.0 * pi) - αₛ * HS) * FX * DJ

               b[i] += dot(fwr, KS, fwr)

            end
         end
      end
   end

   # -----------------------------------------------------------------------
   #Do Singular integraltion now! 
   # zp1 = zp+1, zp2 = xp-1 where zp = cos(pi*(2j-1)/(2N))
   zp = Vector{Float64}(undef, N)
   zp1 = Vector{Float64}(undef, N)
   zp2 = Vector{Float64}(undef, N)

   @inbounds for j in 1:N
      zp[j] = cospi((2 * j - 1) / (2N))
      zp1[j] = 2 * cospi((2 * j - 1) / (4N))^2
      zp2[j] = -2 * sinpi((2 * j - 1) / (4N))^2
   end

   # z1 = (1+z)/2, z2 = (1-z)/2 where z = cos(pi*(2j-1)/(2n)) 
   # are Chebyshev nodes in the interval [-1,1]. This is for 
   # near singular integrals 
   n = ns
   fw = getF1W(ns)

   z = Vector{Float64}(undef, n)
   z1 = Vector{Float64}(undef, n)
   z2 = Vector{Float64}(undef, n)

   d1 = Vector{Float64}(undef, n)
   d2 = Vector{Float64}(undef, n)

   @inbounds for i in 1:n
      c₁ = π * (2 * i - 1) / (2 * n)
      z[i] = cos(c₁)
      z1[i] = cos(c₁ / 2)^2
      z2[i] = sin(c₁ / 2)^2
   end

   wmz₂ = Vector{Float64}(undef, n)
   wz = Vector{Float64}(undef, n)
   wmz = Vector{Float64}(undef, n)

   wfunc!(wz, p, z) #wz = w(z)
   wfunc!(wmz, p, z; α=-1.0) #wmz = w(z)
   wfunc!(wmz₂, p, z2; α=-1.0) #wmz₂ = w(-z2)

   fwm = similar(fw)
   fwl = similar(fw)
   dwfunc!(fwl, p, z)
   dwfunc!(fwm, p, z2)   # fwm := dw(z2)
   @. fwm = fw * fwm
   @. fwl = fw * fwl

   y1 = Vector{Float64}(undef, n)
   wdf = similar(y1)
   y1tmp = similar(y1)
   y2 = similar(y1)
   t1 = Matrix{Float64}(undef, n, n)  # meshgrid of y1/y2 
   t2 = Matrix{Float64}(undef, n, n)  # (column = y2[j], row = y1[i])
   DJ = Matrix{Float64}(undef, n, n)  # To store the Jacobian
   DIF = similar(DJ)
   Zx = similar(DJ)
   Zy = similar(DJ)


   FX = Matrix{Float64}(undef, n, n)
   HS = Matrix{Float64}(undef, n, n)
   βS = Matrix{Float64}(undef, n, n)

   KS = Matrix{Float64}(undef, n, n)
   Df = Vector{Float64}(undef, n)


   A = Matrix{Float64}(undef, n, n)
   dfy = Vector{Float64}(undef, n)
   dfys = Vector{Float64}(undef, n)

   
   #-Singular Integration of Interior patches-
   @inbounds for j in 1:Np
      # j in the linear index of the Chebyshev
      # point. We will compute precomps of all 
      # points τₖ(x₁,x₂) at once by now varying
      # the patches k. This removes some 
      # repeated calculations and thereby saves 
      # time!
      qq, rr = divrem(j - 1, N)
      qq = qq + 1
      rr = rr + 1

      x1, x2 = zp[rr], zp[qq]
      x1p, x2p = zp1[rr], zp1[qq]
      x1m, x2m = zp2[rr], zp2[qq]
      x̃₁ = x1p * x2p / 4.0
      x̃₂ = -x1p * x2m / 4.0
      x̃₃ = -x1m * x2p / 4.0
      x̃₄ = x1m * x2m / 4.0

      #Singular integration is sum of four parts
      #-----------------1st part-----------------

      @. d1 = x1p * wmz₂
      @. d2 = x2p * wmz₂
      @. y1 = x1 - d1
      @. y2 = x2 - d2
      @. y1tmp = 1 - y1

      fill_meshgrid!(t1, t2, y1, y2)

      @inbounds for k in 1:M

         diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, k)

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)   # dfys = dfac(k, 1-y1)

         f!(FX, Zx, Zy)

         @inbounds for jj in 1:n
            @inbounds for ii in 1:n
               r2 = DIF[ii, jj]^2

               dval = dfy[ii]
               ds = dfys[ii]

               logd = log(dval)
               logr2 = log(r2)

               βt1 = logd - logr2
               βt2 = logd

               βt3 = phi2_slr(βt1; α=s)
               βt4 = phi2_slr(βt2; α=s)

               βS[ii, jj] = βt1^2 * βt3 - βt2^2 * βt4
               HS[ii, jj] = ds * expm1(-s * logr2) / s
            end
         end

         @. KS = (βS / (4.0 * pi) - αₛ * HS) * FX * DJ

         b[j+(k-1)*Np] += x̃₁ * dot(fwm, KS, fwm)

      end

      #-----------------2nd part-----------------
      #(reuse d1, y1, TN_y1, y1tmp; update d2,y2,t1,t2)
      @. d2 = x2m * wmz₂
      @. y2 = x2 - d2

      fill_meshgrid!(t1, t2, y1, y2)

      @inbounds for k in 1:M

         diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, k)

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)   # dfys = dfac(k, 1-y1)

         f!(FX, Zx, Zy)

         @inbounds for jj in 1:n
            @inbounds for ii in 1:n
               r2 = DIF[ii, jj]^2

               dval = dfy[ii]
               ds = dfys[ii]

               logd = log(dval)
               logr2 = log(r2)

               βt1 = logd - logr2
               βt2 = logd

               βt3 = phi2_slr(βt1; α=s)
               βt4 = phi2_slr(βt2; α=s)

               βS[ii, jj] = βt1^2 * βt3 - βt2^2 * βt4
               HS[ii, jj] = ds * expm1(-s * logr2) / s
            end
         end

         @. KS = (βS / (4.0 * pi) - αₛ * HS) * FX * DJ

         b[j+(k-1)*Np] += x̃₂ * dot(fwm, KS, fwm)

      end

      #-----------------3rd part-----------------

      @. d1 = x1m * wmz₂
      @. d2 = x2p * wmz₂
      @. y1 = x1 - d1
      @. y2 = x2 - d2
      @. y1tmp = 1 - y1

      fill_meshgrid!(t1, t2, y1, y2)

      @inbounds for k in knbd

         diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, k)

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)   # dfys = dfac(k, 1-y1)

         f!(FX, Zx, Zy)

         @inbounds for jj in 1:n
            @inbounds for ii in 1:n
               r2 = DIF[ii, jj]^2

               dval = dfy[ii]
               ds = dfys[ii]

               logd = log(dval)
               logr2 = log(r2)

               βt1 = logd - logr2
               βt2 = logd

               βt3 = phi2_slr(βt1; α=s)
               βt4 = phi2_slr(βt2; α=s)

               βS[ii, jj] = βt1^2 * βt3 - βt2^2 * βt4
               HS[ii, jj] = ds * expm1(-s * logr2) / s
            end
         end

         @. KS = (βS / (4.0 * pi) - αₛ * HS) * FX * DJ

         b[j+(k-1)*Np] += x̃₃ * dot(fwm, KS, fwm)

      end

      #-----------------4th part-----------------
      #(reuse d1, y1, TN_y1, y1tmp; update d2,y2,t1,t2)
      @. d2 = x2m * wmz₂
      @. y2 = x2 - d2

      fill_meshgrid!(t1, t2, y1, y2)

      @inbounds for k in knbd

         diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, k)

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)   # dfys = dfac(k, 1-y1)

         f!(FX, Zx, Zy)

         @inbounds for jj in 1:n
            @inbounds for ii in 1:n
               r2 = DIF[ii, jj]^2

               dval = dfy[ii]
               ds = dfys[ii]

               logd = log(dval)
               logr2 = log(r2)

               βt1 = logd - logr2
               βt2 = logd

               βt3 = phi2_slr(βt1; α=s)
               βt4 = phi2_slr(βt2; α=s)

               βS[ii, jj] = βt1^2 * βt3 - βt2^2 * βt4
               HS[ii, jj] = ds * expm1(-s * logr2) / s
            end
         end

         @. KS = (βS / (4.0 * pi) - αₛ * HS) * FX * DJ

         b[j+(k-1)*Np] += x̃₄ * dot(fwm, KS, fwm)

      end

   end

   #--Singular Integration of patches touching the boundary--
   @inbounds for j in 1:Np
      # j in the linear index of the Chebyshev
      # point. We will compute precomps of all 
      # points τₖ(x₁,x₂) at once by now varying
      # the patches k. This removes some 
      # repeated calculations and thereby saves 
      # time!
      qq, rr = divrem(j - 1, N)
      qq = qq + 1
      rr = rr + 1
      x1, x2 = zp[rr], zp[qq]
      x2p = zp1[qq]
      x1m, x2m = zp2[rr], zp2[qq]
      x̃₃ = -x1m * x2p / 4.0
      x̃₄ = x1m * x2m / 4.0

      # Singular and k is a boundary patch!
      #-----------------3rd part-----------------

      @. d1 = x1m * wz / 2.0
      @. d2 = x2p * wmz₂
      @. y1 = x1 - d1
      @. y2 = x2 - d2
      @. y1tmp = -x1m * wmz / 2.0

      fill_meshgrid!(t1, t2, y1, y2)

      @inbounds for k in d.kd

         diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, k)

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)   # dfys = dfac(k, 1-y1)

         f!(FX, Zx, Zy)

         @inbounds for jj in 1:n
            @inbounds for ii in 1:n
               r2 = DIF[ii, jj]^2

               dval = dfy[ii]
               ds = dfys[ii]

               logd = log(dval)
               logr2 = log(r2)

               βt1 = logd - logr2
               βt2 = logd

               βt3 = phi2_slr(βt1; α=s)
               βt4 = phi2_slr(βt2; α=s)

               βS[ii, jj] = βt1^2 * βt3 - βt2^2 * βt4
               HS[ii, jj] = ds * expm1(-s * logr2) / s
            end
         end

         @. KS = (βS / (4.0 * pi) - αₛ * HS) * FX * DJ

         b[j+(k-1)*Np] += x̃₃ * dot(fwl, KS, fwm)
 
      end

      #-----------------4th part-----------------
      #(reuse d1, y1, TN_y1, y1tmp; update d2,y2,t1,t2)
      @. d2 = x2m * wmz₂
      @. y2 = x2 - d2

      fill_meshgrid!(t1, t2, y1, y2)

      @inbounds for k in d.kd

         diff_map!(DIF, Zx, Zy, DJ, d, x1, x2, t1, t2, d1, d2, k)

         dfunc!(dfy, d, k, y1tmp)

         dfunc!(dfys, d, k, y1tmp, s)   # dfys = dfac(k, 1-y1)

         f!(FX, Zx, Zy)

         @inbounds for jj in 1:n
            @inbounds for ii in 1:n
               r2 = DIF[ii, jj]^2

               dval = dfy[ii]
               ds = dfys[ii]

               logd = log(dval)
               logr2 = log(r2)

               βt1 = logd - logr2
               βt2 = logd

               βt3 = phi2_slr(βt1; α=s)
               βt4 = phi2_slr(βt2; α=s)

               βS[ii, jj] = βt1^2 * βt3 - βt2^2 * βt4
               HS[ii, jj] = ds * expm1(-s * logr2) / s
            end
         end

         @. KS = (βS / (4.0 * pi) - αₛ * HS) * FX * DJ

         b[j+(k-1)*Np] += x̃₄ * dot(fwl, KS, fwm)

      end

   end

   #-Singular Integration For pts on bonudary-
   #These points are ofcourse always on bd patches
   @inbounds for j in 1:N
      # j in the linear index of the Chebyshev
      # point. We will compute precomps of all 
      # points τₖ(x₁,x₂) at once by now varying
      # the patches k.

      x2 = zp[j]
      x2p = zp1[j]
      x2m = zp2[j]
      x̃₁ = x2p / 2.0
      x̃₂ = -x2m / 2.0

      #Singular integration is sum of two parts here
      #-----------------1st part-----------------

      @. d1 = 2.0 * wmz₂
      @. d2 = x2p * wmz₂
      @. y1 = 1.0 - d1
      @. y2 = x2 - d2

      fill_meshgrid!(t1, t2, y1, y2)

      @inbounds for ℓ in 1:Mbd

         k = d.kd[ℓ]

         diff_map!(DIF, Zx, Zy, DJ, d, 1.0, x2, t1, t2, d1, d2, k)

         dfunc!(dfy, d, k, d1)

         dfunc!(dfys, d, k, d1, s)   # dfy = dfac(k, 1-y1=d1)

         f!(FX, Zx, Zy)

         @inbounds for jj in 1:n
            @inbounds for ii in 1:n
               r2 = DIF[ii, jj]^2

               dval = dfy[ii]
               ds = dfys[ii]

               logd = log(dval)
               logr2 = log(r2)

               βt1 = logd - logr2
               βt2 = logd

               βt3 = phi2_slr(βt1; α=s)
               βt4 = phi2_slr(βt2; α=s)

               βS[ii, jj] = βt1^2 * βt3 - βt2^2 * βt4
               HS[ii, jj] = ds * expm1(-s * logr2) / s
            end
         end

         @. KS = (βS / (4.0 * pi) - αₛ * HS) * FX * DJ

         b[Ni+j+(ℓ-1)*N] += x̃₁ * dot(fwm, KS, fwm)

      end

      #-----------------2nd part-----------------
      #(reuse d1, y1, TN_y1; update d2,y2,t1,t2)
      @. d2 = x2m * wmz₂
      @. y2 = x2 - d2

      fill_meshgrid!(t1, t2, y1, y2)

      @inbounds for ℓ in 1:Mbd

         k = d.kd[ℓ]

         diff_map!(DIF, Zx, Zy, DJ, d, 1.0, x2, t1, t2, d1, d2, k)

         dfunc!(dfy, d, k, d1) 

         dfunc!(dfys, d, k, d1, s)   # dfy = dfac(k, 1-y1)

         f!(FX, Zx, Zy)

         @inbounds for jj in 1:n
            @inbounds for ii in 1:n
               r2 = DIF[ii, jj]^2

               dval = dfy[ii]
               ds = dfys[ii]

               logd = log(dval)
               logr2 = log(r2)

               βt1 = logd - logr2
               βt2 = logd

               βt3 = phi2_slr(βt1; α=s)
               βt4 = phi2_slr(βt2; α=s)

               βS[ii, jj] = βt1^2 * βt3 - βt2^2 * βt4
               HS[ii, jj] = ds * expm1(-s * logr2) / s
            end
         end

         @. KS = (βS / (4.0 * pi) - αₛ * HS) * FX * DJ

         b[Ni+j+(ℓ-1)*N] += x̃₂ * dot(fwm, KS, fwm)

      end

   end

   return b

end

