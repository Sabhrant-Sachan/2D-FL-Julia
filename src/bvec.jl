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
