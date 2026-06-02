# ------------------------------------------------------------
# ellipseNh: This domain represents an ellipse with nh ≥ 2 elliptical holes.
#
# If nh == 0, use ellipse. If nh == 1, use annulus.If nh ≥ 2, use ellipseNh.
#
# Parameters:
#
#   A, B: Vectors of length nh + 1.
#       (A[1], B[1]) is the center of the outer ellipse.
#       (A[j+1], B[j+1]) is the center of the j-th hole.
#
#   R1, R2: Vectors of length nh + 1.
#       R1[j], R2[j] are the semi-axis lengths of ellipse j.
#       j = 1 is the outer ellipse. j = 2:nh+1 are the holes.
#
#   tht:Vector of length nh + 1.
#       tht[j] is the rotation angle of ellipse j.
#
#   nh: Number of holes.
#
#   T:  A 13 × nh matrix storing the angular and point data needed
#       to define the 8 regions around each hole.
#
#       For hole j:
#          T[1:4, j]   : outer ellipse angle breakpoints
#          T[5:9, j]   : hole ellipse angle breakpoints
#          T[10:11, j] : midpoint M_j
#          T[12:13, j] : previous midpoint M_{j-1}
#
#   pths:Vector of Patch.
#        Each patch stores:
#          reg, ck0, ck1, tk0, tk1
#
#        Actual region numbers are:
#          reg = rloc + 8*(hole_index - 1),
#          where rloc ∈ 1:8.
#   kd:
#       Boundary-touching patch indices.
#       For ellipseNh, this follows the MATLAB constructor:
#       every patch with ck1 == 1.0 is added to kd.
#   Qpts:
#       8 × Npat matrix of bounding quadrilaterals for all patches.
#
#   Qptsbd:
#       8 × length(kd) matrix of boundary quadrilaterals.
# ------------------------------------------------------------

mutable struct ellipseNh <: abstractdomain

   A::Vector{Float64}
   B::Vector{Float64}

   R1::Vector{Float64}
   R2::Vector{Float64}

   tht::Vector{Float64}

   nh::Int

   T::Matrix{Float64}

   kd::Vector{Int}
   Npat::Int
   pths::Vector{Patch}

   Qpts::Matrix{Float64}
   Qptsbd::Matrix{Float64}

   function ellipseNh(; b,
      nh::Int=3,
      a=nothing,
      A1=nothing, B1=nothing,
      Lm=nothing, tht=nothing,
      R1=nothing, R2=nothing,
      ck=nothing, tk=nothing)

      @assert isa(b, AbstractVector{Int}) && length(b) == 8
      @assert all(b .>= 1)
      @assert nh >= 2 "ellipseNh is intended for nh ≥ 2; use annulus for nh = 1."

      A1 = something(A1, 0.0)
      B1 = something(B1, 0.0)
      Lm = something(Lm, 0.5)

      L = 0.45

      # Rotation angles.
      if isnothing(tht)

         tht_vec = zeros(Float64, nh + 1)

      elseif isa(tht, Number)

         tht_vec = fill(Float64(tht), nh + 1)

      else

         tht_vec = collect(Float64, tht)
         @assert length(tht_vec) == nh + 1 "tht must have length nh + 1."
      end

      # If scalar R1/R2 are supplied 
      #   R1 = R1 * [1, ones(1,nh)/(nh+2)]
      #   R2 = R2 * [1, ones(1,nh)/(nh+2)]
      R1_vec = Vector{Float64}(undef, nh + 1)
      R2_vec = Vector{Float64}(undef, nh + 1)

      if isnothing(R1)

         R1_vec[1] = 1.0
         for j in 1:nh
            R1_vec[j+1] = 1 / (nh + 2)
         end

      elseif isa(R1, Number)

         R1_vec[1] = Float64(R1)
         for j in 1:nh
            R1_vec[j+1] = R1_vec[1] / (nh + 2)
         end
         @assert R1 > 0

      else

         R1_vec = collect(Float64, R1)
         @assert all(R1_vec .> 0.0)
         @assert length(R1_vec) == nh + 1 "R1 must have length nh + 1."

      end

      if isnothing(R2)

         R2_vec[1] = 1.0
         for j in 1:nh
            R2_vec[j+1] = 1 / (nh + 2)
         end

      elseif isa(R2, Number)

         R2_vec[1] = Float64(R2)
         for j in 1:nh
            R2_vec[j+1] = R2_vec[1] / (nh + 2)
         end
         @assert R2 > 0
      else

         R2_vec = collect(Float64, R2)
         @assert all(R2_vec .> 0.0)
         @assert length(R2_vec) == nh + 1 "R2 must have length nh + 1."

      end

      # Compute centers. The hole centers are placed on the smaller ellipse
      Avec = zeros(Float64, nh + 1)
      Bvec = zeros(Float64, nh + 1)

      Avec[1], Bvec[1] = A1, B1

      s, c = sincos(tht_vec[1])

      @inbounds for i in 1:nh

         tt = 2.0 * (i - 1) / nh

         Avec[i+1], Bvec[i+1] = A1, B1

         st, ct = sincos(π * tt)

         Avec[i+1] += L * (R1_vec[1] * ct * c - R2_vec[1] * st * s)
         Bvec[i+1] += L * (R1_vec[1] * ct * s + R2_vec[1] * st * c)

      end

      # Build T.
      T = zeros(Float64, 13, nh)

      # Outer ellipse frame.
      s1, c1 = sincos(tht_vec[1])
      e1x1 = c1
      e1y1 = s1
      e2x1 = -s1
      e2y1 = c1

      # ------------------------------------------------------------
      # STEP 1.
      # For each pair of consecutive holes, compute midpoint M_j and
      # the corresponding outer-ellipse angle T[4,j].
      # ------------------------------------------------------------
      @inbounds for j in 1:nh

         hcur = j + 1
         hnext = (j == nh) ? 2 : (j + 2)

         Ccurx = Avec[hcur]
         Ccury = Bvec[hcur]

         Cnextx = Avec[hnext]
         Cnexty = Bvec[hnext]

         # Midpoint between current and next hole centers.
         Mx = 0.5 * (Cnextx + Ccurx)
         My = 0.5 * (Cnexty + Ccury)

         # Cp = [B, -A], so Cp(next) - Cp(cur).
         Vx = Cnexty - Ccury
         Vy = -Cnextx + Ccurx

         Et, Lt = ellipserinter(
            A1, B1, R1_vec[1], R2_vec[1], e1x1,
            e1y1, e2x1, e2y1, Mx, My, Vx, Vy)

         @assert isfinite(Et) && isfinite(Lt) "ellipseNh constructor: ellipserinter failed in STEP 1, j=$j"

         T[10, j] = Mx + Lm * Lt * Vx
         T[11, j] = My + Lm * Lt * Vy

         T[4, j] = Et

         if j < nh
            T[12, j+1] = T[10, j]
            T[13, j+1] = T[11, j]
            T[1, j+1] = T[4, j]
         else
            T[12, 1] = T[10, j]
            T[13, 1] = T[11, j]
            T[1, 1] = T[4, j] - 2.0 * π
         end

      end

      # ------------------------------------------------------------
      # STEP 2.
      #
      # Compute T[7:9,j], the three angle breakpoints on the j-th hole.
      # ------------------------------------------------------------
      @inbounds for j in 1:nh

         hcur = j + 1

         sj, cj = sincos(tht_vec[hcur])

         e1xj = cj
         e1yj = sj
         e2xj = -sj
         e2yj = cj

         Cx = Avec[hcur]
         Cy = Bvec[hcur]

         # Current midpoint M_j to current hole center.
         Mx = T[10, j]
         My = T[11, j]

         Et1, _ = ellipserinter(Cx, Cy, R1_vec[hcur], R2_vec[hcur],
            e1xj, e1yj, e2xj, e2yj, Mx, My, Cx - Mx, Cy - My)

         @assert isfinite(Et1) "ellipserinter failed in STEP 2a, j=$j"

         # Outer center to current hole center.
         Et2, _ = ellipserinter(Cx, Cy, R1_vec[hcur], R2_vec[hcur],
            e1xj, e1yj, e2xj, e2yj, A1, B1, Cx - A1, Cy - B1)

         @assert isfinite(Et2) "ellipserinter failed in STEP 2b, j=$j"

         # Previous midpoint M_{j-1} to current hole center.
         Mx = T[12, j]
         My = T[13, j]

         Et3, _ = ellipserinter(Cx, Cy, R1_vec[hcur], R2_vec[hcur],
            e1xj, e1yj, e2xj, e2yj, Mx, My, Cx - Mx, Cy - My)

         @assert isfinite(Et3) "ellipserinter failed in STEP 2c, j=$j"

         T[7, j] = Et1
         T[8, j] = Et2
         T[9, j] = Et3

      end

      # ------------------------------------------------------------
      # STEP 3.
      #
      # Compute T[5:6,j] on each hole and T[2:3,j] on the outer ellipse.
      # ------------------------------------------------------------
      alp = π / (nh + 1)

      @inbounds for hcur in 2:(nh+1)

         j = hcur - 1

         sh, ch = sincos(tht_vec[hcur])

         t1 = if T[7, j] > T[9, j]
            0.5 * (T[7, j] + T[9, j])
         else
            0.5 * (T[7, j] + T[9, j]) + π
         end

         tm = t1 - alp
         tp = t1 + alp

         T[5, j] = tm
         T[6, j] = tp

         stm, ctm = sincos(tm)

         gx = Avec[hcur] + R1_vec[hcur] * ctm * ch - R2_vec[hcur] * stm * sh
         gy = Bvec[hcur] + R1_vec[hcur] * ctm * sh + R2_vec[hcur] * stm * ch

         dgx = -R1_vec[hcur] * stm * ch - R2_vec[hcur] * ctm * sh
         dgy = -R1_vec[hcur] * stm * sh + R2_vec[hcur] * ctm * ch

         Et1, _ = ellipserinter(A1, B1, R1_vec[1], R2_vec[1],
            e1x1, e1y1, e2x1, e2y1, gx, gy, dgy, -dgx)

         @assert isfinite(Et1) "ellipserinter failed in STEP 3a, j=$j"

         stp, ctp = sincos(tp)

         gx = Avec[hcur] + R1_vec[hcur] * ctp * ch - R2_vec[hcur] * stp * sh
         gy = Bvec[hcur] + R1_vec[hcur] * ctp * sh + R2_vec[hcur] * stp * ch

         dgx = -R1_vec[hcur] * stp * ch - R2_vec[hcur] * ctp * sh
         dgy = -R1_vec[hcur] * stp * sh + R2_vec[hcur] * ctp * ch

         Et2, _ = ellipserinter(A1, B1, R1_vec[1], R2_vec[1],
            e1x1, e1y1, e2x1, e2y1, gx, gy, dgy, -dgx)

         @assert isfinite(Et2) "ellipserinter failed in STEP 3b, j=$j"

         T[2, j] = Et1
         T[3, j] = Et2

      end

      # Fix angular ordering.
      @inbounds for j in 1:nh

         for kk in 1:3
            if T[kk, j] > T[kk+1, j]
               T[kk, j] -= 2.0π
            end
         end

         if T[5, j] > π / 2.0
            for kk in 5:9
               T[kk, j] -= 2.0π
            end
         elseif T[5, j] < -π / 2.0
            for kk in 5:9
               T[kk, j] += 2.0π
            end
         end

         for kk in 5:8
            if T[kk, j] > T[kk+1, j]
               T[kk+1, j] += 2.0π
            end
         end

      end

      # Patch partitions
      if isnothing(a)
         avec = ceil.(Int, 2 .* b ./ 3)
      else
         avec = collect(Int, a)
         @assert length(avec) == 8
         @assert all(avec .>= 1)
      end

      ck_vec = Vector{Vector{Float64}}(undef, 8)
      tk_vec = Vector{Vector{Float64}}(undef, 8)

      if isnothing(ck)

         @inbounds for kk in 1:8
            ck_vec[kk] = collect((0:avec[kk]) ./ avec[kk])
         end

      elseif isa(ck, AbstractMatrix)

         @inbounds for kk in 1:8
            ck_vec[kk] = collect(Float64, @view ck[1:(avec[kk]+1), kk])
         end

      else

         @assert length(ck) == 8
         @inbounds for kk in 1:8
            ck_vec[kk] = collect(Float64, ck[kk])
         end

      end

      if isnothing(tk)

         @inbounds for kk in 1:8
            tk_vec[kk] = collect((0:b[kk]) ./ b[kk])
         end

      elseif isa(tk, AbstractMatrix)

         @inbounds for kk in 1:8
            tk_vec[kk] = collect(Float64, @view tk[1:(b[kk]+1), kk])
         end

      else

         @assert length(tk) == 8
         @inbounds for kk in 1:8
            tk_vec[kk] = collect(Float64, tk[kk])
         end

      end

      @inbounds for kk in 1:8
         @assert length(ck_vec[kk]) == avec[kk] + 1 "bad ck partition in local region $kk"
         @assert length(tk_vec[kk]) == b[kk] + 1 "bad tk partition in local region $kk"
      end

      # Build patches.
      Npat = nh * sum(avec .* b)
      pths = Vector{Patch}(undef, Npat)

      ip = 1

      @inbounds for l in 1:nh

         for rloc in 1:8

            reg = rloc + 8 * (l - 1)

            for ic in 1:avec[rloc]

               ck0 = ck_vec[rloc][ic]
               ck1 = ck_vec[rloc][ic+1]

               for it in 1:b[rloc]

                  tk0 = tk_vec[rloc][it]
                  tk1 = tk_vec[rloc][it+1]

                  pths[ip] = Patch(reg, ck0, ck1, tk0, tk1)

                  ip += 1

               end

            end

         end

      end

      # Boundary patch indices.
      kd = Int[]

      @inbounds for k in 1:Npat
         if pths[k].ck1 == 1.0
            push!(kd, k)
         end
      end

      # Allocate Qpts and Qptsbd
      Qpts = Matrix{Float64}(undef, 8, Npat)
      Qptsbd = Matrix{Float64}(undef, 8, length(kd))

      d = new(Avec, Bvec, R1_vec, R2_vec, tht_vec,
         nh, T, kd, Npat, pths, Qpts, Qptsbd)

      @inbounds for k in 1:Npat
         @views boundquad!(d.Qpts[:, k], d, k)
      end

      @inbounds for (ℓ, k) in enumerate(d.kd)
         @views boundquadbd!(d.Qptsbd[:, ℓ], d, k)
      end

      return d
   end

end

function mapx(d::ellipseNh, u::Float64, v::Float64, k::Int)::Float64

   Zx, _ = mapxy(d, u, v, k)

   return Zx

end

function mapy(d::ellipseNh, u::Float64, v::Float64, k::Int)::Float64

   _, Zy = mapxy(d, u, v, k)

   return Zy

end

function mapxy(d::ellipseNh, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   ξ1 = muladd(0.5 * hc, u, 0.5 * (p.ck0 + p.ck1))
   ξ2 = muladd(0.5 * ht, v, 0.5 * (p.tk0 + p.tk1))

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   Tj = @view d.T[:, H]

   s1, c1 = sincos(d.tht[1])
   sj, cj = sincos(d.tht[Hj])

   A1 = d.A[1]
   B1 = d.B[1]
   R11 = d.R1[1]
   R12 = d.R2[1]

   Aj = d.A[Hj]
   Bj = d.B[Hj]
   Rj1 = d.R1[Hj]
   Rj2 = d.R2[Hj]

   if r == 1

      st5, ct5 = sincos(Tj[5])
      st2, ct2 = sincos(Tj[2])

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      Xx = (1.0 - ξ2) * Tj[12] + ξ2 * (gjx5 + g1x2) / 2.0
      Xy = (1.0 - ξ2) * Tj[13] + ξ2 * (gjy5 + g1y2) / 2.0

      θ = muladd(ξ2, Tj[2] - Tj[1], Tj[1])
      st, ct = sincos(θ)

      Yx = A1 + R11 * ct * c1 - R12 * st * s1
      Yy = B1 + R11 * ct * s1 + R12 * st * c1

   elseif r == 2

      θ = muladd(ξ2, Tj[3] - Tj[2], Tj[2])
      st, ct = sincos(θ)

      Yx = A1 + R11 * ct * c1 - R12 * st * s1
      Yy = B1 + R11 * ct * s1 + R12 * st * c1

      θj = muladd(ξ2, Tj[6] - Tj[5], Tj[5])
      stj, ctj = sincos(θj)

      Yjx = Aj + Rj1 * ctj * cj - Rj2 * stj * sj
      Yjy = Bj + Rj1 * ctj * sj + Rj2 * stj * cj

      Xx = (Yx + Yjx) / 2.0
      Xy = (Yy + Yjy) / 2.0

   elseif r == 3

      st6, ct6 = sincos(Tj[6])
      st3, ct3 = sincos(Tj[3])

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      Xx = (1.0 - ξ2) * (gjx6 + g1x3) / 2.0 + ξ2 * Tj[10]
      Xy = (1.0 - ξ2) * (gjy6 + g1y3) / 2.0 + ξ2 * Tj[11]

      θ = muladd(ξ2, Tj[4] - Tj[3], Tj[3])
      st, ct = sincos(θ)

      Yx = A1 + R11 * ct * c1 - R12 * st * s1
      Yy = B1 + R11 * ct * s1 + R12 * st * c1

   elseif r == 4

      θj = muladd(ξ2, Tj[6] - Tj[5], Tj[5])
      stj, ctj = sincos(θj)

      Yx = Aj + Rj1 * ctj * cj - Rj2 * stj * sj
      Yy = Bj + Rj1 * ctj * sj + Rj2 * stj * cj

      θ = muladd(ξ2, Tj[3] - Tj[2], Tj[2])
      st, ct = sincos(θ)

      Yjx = A1 + R11 * ct * c1 - R12 * st * s1
      Yjy = B1 + R11 * ct * s1 + R12 * st * c1

      Xx = (Yjx + Yx) / 2.0
      Xy = (Yjy + Yy) / 2.0

   elseif r == 5

      st6, ct6 = sincos(Tj[6])
      st3, ct3 = sincos(Tj[3])

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      Xx = (1.0 - ξ2) * (gjx6 + g1x3) / 2.0 + ξ2 * Tj[10]
      Xy = (1.0 - ξ2) * (gjy6 + g1y3) / 2.0 + ξ2 * Tj[11]

      θ = muladd(ξ2, Tj[7] - Tj[6], Tj[6])
      stj, ctj = sincos(θ)

      Yx = Aj + Rj1 * ctj * cj - Rj2 * stj * sj
      Yy = Bj + Rj1 * ctj * sj + Rj2 * stj * cj

   elseif r == 6

      Xx = (1.0 - ξ2) * Tj[10] + ξ2 * A1
      Xy = (1.0 - ξ2) * Tj[11] + ξ2 * B1

      θ = muladd(ξ2, Tj[8] - Tj[7], Tj[7])
      stj, ctj = sincos(θ)

      Yx = Aj + Rj1 * ctj * cj - Rj2 * stj * sj
      Yy = Bj + Rj1 * ctj * sj + Rj2 * stj * cj

   elseif r == 7

      Xx = (1.0 - ξ2) * A1 + ξ2 * Tj[12]
      Xy = (1.0 - ξ2) * B1 + ξ2 * Tj[13]

      θ = muladd(ξ2, Tj[9] - Tj[8], Tj[8])
      stj, ctj = sincos(θ)

      Yx = Aj + Rj1 * ctj * cj - Rj2 * stj * sj
      Yy = Bj + Rj1 * ctj * sj + Rj2 * stj * cj

   elseif r == 8

      st5, ct5 = sincos(Tj[5])
      st2, ct2 = sincos(Tj[2])

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      Xx = (1.0 - ξ2) * Tj[12] + ξ2 * (gjx5 + g1x2) / 2.0
      Xy = (1.0 - ξ2) * Tj[13] + ξ2 * (gjy5 + g1y2) / 2.0

      θ = muladd(ξ2, Tj[5] - Tj[9] + 2.0π, Tj[9])
      stj, ctj = sincos(θ)

      Yx = Aj + Rj1 * ctj * cj - Rj2 * stj * sj
      Yy = Bj + Rj1 * ctj * sj + Rj2 * stj * cj

   else

      throw(ArgumentError("mapxy local region must be 1:8; got r=$r, reg=$reg"))

   end

   Zx = (1.0 - ξ1) * Xx + ξ1 * Yx
   Zy = (1.0 - ξ1) * Xy + ξ1 * Yy

   return Zx, Zy

end

function mapxy!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
   d::ellipseNh, u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   αc = 0.5 * (p.ck1 - p.ck0)
   βc = 0.5 * (p.ck0 + p.ck1)

   αt = 0.5 * (p.tk1 - p.tk0)
   βt = 0.5 * (p.tk0 + p.tk1)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   # T values for this hole block.
   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]
   T10 = d.T[10, H]
   T11 = d.T[11, H]
   T12 = d.T[12, H]
   T13 = d.T[13, H]

   # Outer ellipse constants.
   s1, c1 = sincos(d.tht[1])
   A1 = d.A[1]
   B1 = d.B[1]
   R11 = d.R1[1]
   R12 = d.R2[1]

   # Hole ellipse constants.
   sj, cj = sincos(d.tht[Hj])
   Aj = d.A[Hj]
   Bj = d.B[Hj]
   Rj1 = d.R1[Hj]
   Rj2 = d.R2[Hj]

   if r == 1

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      XxR = 0.5 * (gjx5 + g1x2)
      XyR = 0.5 * (gjy5 + g1y2)

      dθ = T2 - T1

      @inbounds for i in eachindex(Zx, Zy, u, v)

         ξ1 = muladd(αc, u[i], βc)
         ξ2 = muladd(αt, v[i], βt)

         Xx = muladd(ξ2, XxR - T12, T12)
         Xy = muladd(ξ2, XyR - T13, T13)

         θ = muladd(ξ2, dθ, T1)
         st, ct = sincos(θ)

         Yx = A1 + R11 * ct * c1 - R12 * st * s1
         Yy = B1 + R11 * ct * s1 + R12 * st * c1

         Zx[i] = muladd(ξ1, Yx - Xx, Xx)
         Zy[i] = muladd(ξ1, Yy - Xy, Xy)

      end

   elseif r == 2

      dθo = T3 - T2
      dθh = T6 - T5

      @inbounds for i in eachindex(Zx, Zy, u, v)

         ξ1 = muladd(αc, u[i], βc)
         ξ2 = muladd(αt, v[i], βt)

         θ = muladd(ξ2, dθo, T2)
         st, ct = sincos(θ)

         Yx = A1 + R11 * ct * c1 - R12 * st * s1
         Yy = B1 + R11 * ct * s1 + R12 * st * c1

         θj = muladd(ξ2, dθh, T5)
         stj, ctj = sincos(θj)

         Yjx = Aj + Rj1 * ctj * cj - Rj2 * stj * sj
         Yjy = Bj + Rj1 * ctj * sj + Rj2 * stj * cj

         Xx = 0.5 * (Yx + Yjx)
         Xy = 0.5 * (Yy + Yjy)

         Zx[i] = muladd(ξ1, Yx - Xx, Xx)
         Zy[i] = muladd(ξ1, Yy - Xy, Xy)

      end

   elseif r == 3

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      XxL = 0.5 * (gjx6 + g1x3)
      XyL = 0.5 * (gjy6 + g1y3)

      dθ = T4 - T3

      @inbounds for i in eachindex(Zx, Zy, u, v)

         ξ1 = muladd(αc, u[i], βc)
         ξ2 = muladd(αt, v[i], βt)

         Xx = muladd(ξ2, T10 - XxL, XxL)
         Xy = muladd(ξ2, T11 - XyL, XyL)

         θ = muladd(ξ2, dθ, T3)
         st, ct = sincos(θ)

         Yx = A1 + R11 * ct * c1 - R12 * st * s1
         Yy = B1 + R11 * ct * s1 + R12 * st * c1

         Zx[i] = muladd(ξ1, Yx - Xx, Xx)
         Zy[i] = muladd(ξ1, Yy - Xy, Xy)

      end

   elseif r == 4

      dθh = T6 - T5
      dθo = T3 - T2

      @inbounds for i in eachindex(Zx, Zy, u, v)

         ξ1 = muladd(αc, u[i], βc)
         ξ2 = muladd(αt, v[i], βt)

         θj = muladd(ξ2, dθh, T5)
         stj, ctj = sincos(θj)

         Yx = Aj + Rj1 * ctj * cj - Rj2 * stj * sj
         Yy = Bj + Rj1 * ctj * sj + Rj2 * stj * cj

         θ = muladd(ξ2, dθo, T2)
         st, ct = sincos(θ)

         Yjx = A1 + R11 * ct * c1 - R12 * st * s1
         Yjy = B1 + R11 * ct * s1 + R12 * st * c1

         Xx = 0.5 * (Yjx + Yx)
         Xy = 0.5 * (Yjy + Yy)

         Zx[i] = muladd(ξ1, Yx - Xx, Xx)
         Zy[i] = muladd(ξ1, Yy - Xy, Xy)

      end

   elseif r == 5

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      XxL = 0.5 * (gjx6 + g1x3)
      XyL = 0.5 * (gjy6 + g1y3)

      dθ = T7 - T6

      @inbounds for i in eachindex(Zx, Zy, u, v)

         ξ1 = muladd(αc, u[i], βc)
         ξ2 = muladd(αt, v[i], βt)

         Xx = muladd(ξ2, T10 - XxL, XxL)
         Xy = muladd(ξ2, T11 - XyL, XyL)

         θ = muladd(ξ2, dθ, T6)
         stj, ctj = sincos(θ)

         Yx = Aj + Rj1 * ctj * cj - Rj2 * stj * sj
         Yy = Bj + Rj1 * ctj * sj + Rj2 * stj * cj

         Zx[i] = muladd(ξ1, Yx - Xx, Xx)
         Zy[i] = muladd(ξ1, Yy - Xy, Xy)

      end

   elseif r == 6

      dθ = T8 - T7

      @inbounds for i in eachindex(Zx, Zy, u, v)

         ξ1 = muladd(αc, u[i], βc)
         ξ2 = muladd(αt, v[i], βt)

         Xx = muladd(ξ2, A1 - T10, T10)
         Xy = muladd(ξ2, B1 - T11, T11)

         θ = muladd(ξ2, dθ, T7)
         stj, ctj = sincos(θ)

         Yx = Aj + Rj1 * ctj * cj - Rj2 * stj * sj
         Yy = Bj + Rj1 * ctj * sj + Rj2 * stj * cj

         Zx[i] = muladd(ξ1, Yx - Xx, Xx)
         Zy[i] = muladd(ξ1, Yy - Xy, Xy)

      end

   elseif r == 7

      dθ = T9 - T8

      @inbounds for i in eachindex(Zx, Zy, u, v)

         ξ1 = muladd(αc, u[i], βc)
         ξ2 = muladd(αt, v[i], βt)

         Xx = muladd(ξ2, T12 - A1, A1)
         Xy = muladd(ξ2, T13 - B1, B1)

         θ = muladd(ξ2, dθ, T8)
         stj, ctj = sincos(θ)

         Yx = Aj + Rj1 * ctj * cj - Rj2 * stj * sj
         Yy = Bj + Rj1 * ctj * sj + Rj2 * stj * cj

         Zx[i] = muladd(ξ1, Yx - Xx, Xx)
         Zy[i] = muladd(ξ1, Yy - Xy, Xy)

      end

   elseif r == 8

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      XxR = 0.5 * (gjx5 + g1x2)
      XyR = 0.5 * (gjy5 + g1y2)

      dθ = T5 - T9 + 2.0π

      @inbounds for i in eachindex(Zx, Zy, u, v)

         ξ1 = muladd(αc, u[i], βc)
         ξ2 = muladd(αt, v[i], βt)

         Xx = muladd(ξ2, XxR - T12, T12)
         Xy = muladd(ξ2, XyR - T13, T13)

         θ = muladd(ξ2, dθ, T9)
         stj, ctj = sincos(θ)

         Yx = Aj + Rj1 * ctj * cj - Rj2 * stj * sj
         Yy = Bj + Rj1 * ctj * sj + Rj2 * stj * cj

         Zx[i] = muladd(ξ1, Yx - Xx, Xx)
         Zy[i] = muladd(ξ1, Yy - Xy, Xy)

      end

   else

      throw(ArgumentError("mapxy! local region must be 1:8; got r=$r, reg=$reg"))

   end

   return nothing

end

function draw(d::ellipseNh, flag=nothing; L::Int=33, show::Bool=true)

   clrs = (
      RGBf(0.8500, 0.3250, 0.0980),  # dark orange
      RGBf(0.4940, 0.1840, 0.5560),  # deep purple
      RGBf(0.4660, 0.6740, 0.1880),  # sea green
      RGBf(0.9290, 0.6270, 0.8900),  # light pink
      RGBf(0.3010, 0.7450, 0.9330),  # sky blue
      RGBf(1.0000, 0.7320, 0.4000),  # coral
      RGBf(0.6350, 0.0780, 0.1840),  # maroon
      RGBf(0.9880, 0.4350, 0.7880),  # light purple
   )

   colors = Vector{RGBf}(undef, 8 * d.nh)

   @inbounds for jh in 1:d.nh
      jj = 8 * (jh - 1)

      colors[jj+1] = clrs[1]
      colors[jj+2] = clrs[2]
      colors[jj+3] = clrs[3]
      colors[jj+4] = clrs[4]
      colors[jj+5] = clrs[5]
      colors[jj+6] = clrs[6]
      colors[jj+7] = clrs[7]
      colors[jj+8] = clrs[8]
   end

   return draw_geom(d, colors; flag=flag, L=L, show=show)

end

#-----------------------
"""
  gamx(d::ellipseNh, t::Float64, k::Int)
  gamx!(out::StridedArray{Float64}, d::ellipseNh, t::StridedArray{Float64}, k)
First (x) coordinate of the boundary parametrization γ(t).

- `gamx(d, t, k)` computes the **right boundary** of patch `k`, where `t ∈ [-1,1]`.

`t` may be a scalar or an array; the return has the same shape.
"""
function gamx(d::ellipseNh, t::Float64, k::Int)::Float64

   p = d.pths[k]

   αt = 0.5 * (p.tk1 - p.tk0)
   βt = 0.5 * (p.tk0 + p.tk1)

   ξ2 = muladd(αt, t, βt)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if r == 1

      θ = muladd(ξ2, T2 - T1, T1)

      st, ct = sincos(θ)
      s1, c1 = sincos(d.tht[1])

      return d.A[1] + d.R1[1] * ct * c1 - d.R2[1] * st * s1

   elseif r == 2

      θ = muladd(ξ2, T3 - T2, T2)

      st, ct = sincos(θ)
      s1, c1 = sincos(d.tht[1])

      return d.A[1] + d.R1[1] * ct * c1 - d.R2[1] * st * s1

   elseif r == 3

      θ = muladd(ξ2, T4 - T3, T3)

      st, ct = sincos(θ)
      s1, c1 = sincos(d.tht[1])

      return d.A[1] + d.R1[1] * ct * c1 - d.R2[1] * st * s1

   elseif r == 4

      θ = muladd(ξ2, T6 - T5, T5)

      st, ct = sincos(θ)
      sj, cj = sincos(d.tht[Hj])

      return d.A[Hj] + d.R1[Hj] * ct * cj - d.R2[Hj] * st * sj

   elseif r == 5

      θ = muladd(ξ2, T7 - T6, T6)

      st, ct = sincos(θ)
      sj, cj = sincos(d.tht[Hj])

      return d.A[Hj] + d.R1[Hj] * ct * cj - d.R2[Hj] * st * sj

   elseif r == 6

      θ = muladd(ξ2, T8 - T7, T7)

      st, ct = sincos(θ)
      sj, cj = sincos(d.tht[Hj])

      return d.A[Hj] + d.R1[Hj] * ct * cj - d.R2[Hj] * st * sj

   elseif r == 7

      θ = muladd(ξ2, T9 - T8, T8)

      st, ct = sincos(θ)
      sj, cj = sincos(d.tht[Hj])

      return d.A[Hj] + d.R1[Hj] * ct * cj - d.R2[Hj] * st * sj

   elseif r == 8

      θ = muladd(ξ2, T5 - T9 + 2.0π, T9)

      st, ct = sincos(θ)
      sj, cj = sincos(d.tht[Hj])

      return d.A[Hj] + d.R1[Hj] * ct * cj - d.R2[Hj] * st * sj

   else

      throw(ArgumentError("gamx local region must be 1:8; got r=$r, reg=$reg"))

   end

end

function gamx!(out::StridedArray{Float64}, d::ellipseNh,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   αt = 0.5 * (p.tk1 - p.tk0)
   βt = 0.5 * (p.tk0 + p.tk1)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if r <= 3

      s1, c1 = sincos(d.tht[1])

      A = d.A[1]
      R1 = d.R1[1]
      R2 = d.R2[1]

      θ0, dθ = if r == 1
         T1, T2 - T1
      elseif r == 2
         T2, T3 - T2
      else
         T3, T4 - T3
      end

      @inbounds for i in eachindex(out, t)
         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         out[i] = A + R1 * ct * c1 - R2 * st * s1
      end

   elseif 4 <= r <= 8

      sj, cj = sincos(d.tht[Hj])

      A = d.A[Hj]
      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      θ0, dθ = if r == 4
         T5, T6 - T5
      elseif r == 5
         T6, T7 - T6
      elseif r == 6
         T7, T8 - T7
      elseif r == 7
         T8, T9 - T8
      else
         T9, T5 - T9 + 2.0π
      end

      @inbounds for i in eachindex(out, t)
         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         out[i] = A + R1 * ct * cj - R2 * st * sj
      end

   else

      throw(ArgumentError("gamx! local region must be 1:8; got r=$r, reg=$reg"))

   end

   return nothing

end

"""
    gamy(d::ellipseNh, t::Float64, k::Int)
    gamy!(out::StridedArray{Float64}, d::ellipseNh, t::StridedArray{Float64}, k)

Second (y) coordinate of the boundary parametrization γ(t).

- `gamy(d, t, k)` computes the right boundary of patch `k`, where `t ∈ [-1,1]`.

`t` can be a scalar `Float64` or an array; the return has the same shape.
"""
function gamy(d::ellipseNh, t::Float64, k::Int)::Float64

   p = d.pths[k]

   αt = 0.5 * (p.tk1 - p.tk0)
   βt = 0.5 * (p.tk0 + p.tk1)

   ξ2 = muladd(αt, t, βt)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if r == 1

      θ = muladd(ξ2, T2 - T1, T1)

      st, ct = sincos(θ)
      s1, c1 = sincos(d.tht[1])

      return d.B[1] + d.R1[1] * ct * s1 + d.R2[1] * st * c1

   elseif r == 2

      θ = muladd(ξ2, T3 - T2, T2)

      st, ct = sincos(θ)
      s1, c1 = sincos(d.tht[1])

      return d.B[1] + d.R1[1] * ct * s1 + d.R2[1] * st * c1

   elseif r == 3

      θ = muladd(ξ2, T4 - T3, T3)

      st, ct = sincos(θ)
      s1, c1 = sincos(d.tht[1])

      return d.B[1] + d.R1[1] * ct * s1 + d.R2[1] * st * c1

   elseif r == 4

      θ = muladd(ξ2, T6 - T5, T5)

      st, ct = sincos(θ)
      sj, cj = sincos(d.tht[Hj])

      return d.B[Hj] + d.R1[Hj] * ct * sj + d.R2[Hj] * st * cj

   elseif r == 5

      θ = muladd(ξ2, T7 - T6, T6)

      st, ct = sincos(θ)
      sj, cj = sincos(d.tht[Hj])

      return d.B[Hj] + d.R1[Hj] * ct * sj + d.R2[Hj] * st * cj

   elseif r == 6

      θ = muladd(ξ2, T8 - T7, T7)

      st, ct = sincos(θ)
      sj, cj = sincos(d.tht[Hj])

      return d.B[Hj] + d.R1[Hj] * ct * sj + d.R2[Hj] * st * cj

   elseif r == 7

      θ = muladd(ξ2, T9 - T8, T8)

      st, ct = sincos(θ)
      sj, cj = sincos(d.tht[Hj])

      return d.B[Hj] + d.R1[Hj] * ct * sj + d.R2[Hj] * st * cj

   elseif r == 8

      θ = muladd(ξ2, T5 - T9 + 2.0π, T9)

      st, ct = sincos(θ)
      sj, cj = sincos(d.tht[Hj])

      return d.B[Hj] + d.R1[Hj] * ct * sj + d.R2[Hj] * st * cj

   else

      throw(ArgumentError("gamy(ellipseNh): local region must be 1:8; got r=$r, reg=$reg"))

   end

end

function gamy!(out::StridedArray{Float64}, d::ellipseNh,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   αt = 0.5 * (p.tk1 - p.tk0)
   βt = 0.5 * (p.tk0 + p.tk1)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if r <= 3

      s1, c1 = sincos(d.tht[1])

      B = d.B[1]
      R1 = d.R1[1]
      R2 = d.R2[1]

      θ0, dθ = if r == 1
         T1, T2 - T1
      elseif r == 2
         T2, T3 - T2
      else
         T3, T4 - T3
      end

      @inbounds for i in eachindex(out, t)
         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         out[i] = B + R1 * ct * s1 + R2 * st * c1
      end

   elseif 4 <= r <= 8

      sj, cj = sincos(d.tht[Hj])

      B = d.B[Hj]
      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      θ0, dθ = if r == 4
         T5, T6 - T5
      elseif r == 5
         T6, T7 - T6
      elseif r == 6
         T7, T8 - T7
      elseif r == 7
         T8, T9 - T8
      else
         T9, T5 - T9 + 2.0π
      end

      @inbounds for i in eachindex(out, t)
         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         out[i] = B + R1 * ct * sj + R2 * st * cj
      end

   else

      throw(ArgumentError("gamy!(ellipseNh): local region must be 1:8; got r=$r, reg=$reg"))

   end

   return nothing

end

function gam(d::ellipseNh, t::Float64, k::Int)::Tuple{Float64,Float64}

   return gamx(d, t, k), gamy(d, t, k)

end

function gam!(out::Vector{Float64}, d::ellipseNh, t::Float64, k::Int)

   p = d.pths[k]

   αt = 0.5 * (p.tk1 - p.tk0)
   βt = 0.5 * (p.tk0 + p.tk1)

   ξ2 = muladd(αt, t, βt)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if r <= 3

      sth, cth = sincos(d.tht[1])

      A = d.A[1]
      B = d.B[1]
      R1 = d.R1[1]
      R2 = d.R2[1]

      θ0, dθ = if r == 1
         T1, T2 - T1
      elseif r == 2
         T2, T3 - T2
      else
         T3, T4 - T3
      end

      θ = muladd(ξ2, dθ, θ0)

      st, ct = sincos(θ)

      out[1] = A + R1 * ct * cth - R2 * st * sth
      out[2] = B + R1 * ct * sth + R2 * st * cth

   elseif 4 <= r <= 8

      sth, cth = sincos(d.tht[Hj])

      A = d.A[Hj]
      B = d.B[Hj]
      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      θ0, dθ = if r == 4
         T5, T6 - T5
      elseif r == 5
         T6, T7 - T6
      elseif r == 6
         T7, T8 - T7
      elseif r == 7
         T8, T9 - T8
      else
         T9, T5 - T9 + 2.0π
      end

      θ = muladd(ξ2, dθ, θ0)

      st, ct = sincos(θ)

      out[1] = A + R1 * ct * cth - R2 * st * sth
      out[2] = B + R1 * ct * sth + R2 * st * cth

   else

      throw(ArgumentError("gam!(ellipseNh): local region must be 1:8; got r=$r, reg=$reg"))

   end

   return nothing

end

function gam!(out::Matrix{Float64}, d::ellipseNh, t::Vector{Float64}, k::Int)

   p = d.pths[k]

   αt = 0.5 * (p.tk1 - p.tk0)
   βt = 0.5 * (p.tk0 + p.tk1)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if r <= 3

      sth, cth = sincos(d.tht[1])

      A = d.A[1]
      B = d.B[1]
      R1 = d.R1[1]
      R2 = d.R2[1]

      θ0, dθ = if r == 1
         T1, T2 - T1
      elseif r == 2
         T2, T3 - T2
      else
         T3, T4 - T3
      end

      @inbounds for i in eachindex(t)

         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         out[1, i] = A + R1 * ct * cth - R2 * st * sth
         out[2, i] = B + R1 * ct * sth + R2 * st * cth

      end

   elseif 4 <= r <= 8

      sth, cth = sincos(d.tht[Hj])

      A = d.A[Hj]
      B = d.B[Hj]
      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      θ0, dθ = if r == 4
         T5, T6 - T5
      elseif r == 5
         T6, T7 - T6
      elseif r == 6
         T7, T8 - T7
      elseif r == 7
         T8, T9 - T8
      else
         T9, T5 - T9 + 2.0π
      end

      @inbounds for i in eachindex(t)

         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         out[1, i] = A + R1 * ct * cth - R2 * st * sth
         out[2, i] = B + R1 * ct * sth + R2 * st * cth

      end

   else

      throw(ArgumentError("gam!(ellipseNh): local region must be 1:8; got r=$r, reg=$reg"))

   end

   return nothing

end

function drawbd(d::ellipseNh, flag=nothing; L::Int=33, show::Bool=true)

   clrs = (
      RGBf(0.8500, 0.3250, 0.0980),  # dark orange
      RGBf(0.4940, 0.1840, 0.5560),  # deep purple
      RGBf(0.4660, 0.6740, 0.1880),  # sea green
      RGBf(0.9290, 0.6270, 0.8900),  # light pink
      RGBf(0.3010, 0.7450, 0.9330),  # sky blue
      RGBf(1.0000, 0.7320, 0.4000),  # coral
      RGBf(0.6350, 0.0780, 0.1840),  # maroon
      RGBf(0.9880, 0.4350, 0.7880),  # light purple
   )

   colors = Vector{RGBf}(undef, 8 * d.nh)

   @inbounds for jh in 1:d.nh
      jj = 8 * (jh - 1)

      colors[jj+1] = clrs[1]
      colors[jj+2] = clrs[2]
      colors[jj+3] = clrs[3]
      colors[jj+4] = clrs[4]
      colors[jj+5] = clrs[5]
      colors[jj+6] = clrs[6]
      colors[jj+7] = clrs[7]
      colors[jj+8] = clrs[8]
   end

   return drawbd_geom(d, colors; flag=flag, L=L, show=show)

end

"""
    dgamx(d::ellipseNh, t::Float64, k::Int) -> Float64

Derivative of the first coordinate of the boundary parametrization.

- With `k`: derivative of the **right boundary** of patch `k`
  (i.e., `∂/∂t mapx(d, 1, t, k)` for `t ∈ [-1,1]`).
"""
function dgamx(d::ellipseNh, t::Float64, k::Int)::Float64

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = 0.5 * (p.tk0 + p.tk1)

   ξ2 = muladd(αt, t, βt)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if r <= 3

      sth, cth = sincos(d.tht[1])

      R1 = d.R1[1]
      R2 = d.R2[1]

      θ0, dθ = if r == 1
         T1, T2 - T1
      elseif r == 2
         T2, T3 - T2
      else
         T3, T4 - T3
      end

      θ = muladd(ξ2, dθ, θ0)
      st, ct = sincos(θ)

      return -αt * dθ * (R1 * st * cth + R2 * ct * sth)

   elseif 4 <= r <= 8

      sth, cth = sincos(d.tht[Hj])

      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      θ0, dθ = if r == 4
         T5, T6 - T5
      elseif r == 5
         T6, T7 - T6
      elseif r == 6
         T7, T8 - T7
      elseif r == 7
         T8, T9 - T8
      else
         T9, T5 - T9 + 2.0π
      end

      θ = muladd(ξ2, dθ, θ0)
      st, ct = sincos(θ)

      return -αt * dθ * (R1 * st * cth + R2 * ct * sth)

   else

      throw(ArgumentError("dgamx local region must be 1:8; got r=$r, reg=$reg"))

   end

end

function dgamx!(out::StridedArray{Float64}, d::ellipseNh,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = 0.5 * (p.tk0 + p.tk1)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if r <= 3

      sth, cth = sincos(d.tht[1])

      R1 = d.R1[1]
      R2 = d.R2[1]

      θ0, dθ = if r == 1
         T1, T2 - T1
      elseif r == 2
         T2, T3 - T2
      else
         T3, T4 - T3
      end

      cfac = -αt * dθ

      @inbounds for i in eachindex(out, t)
         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         out[i] = cfac * (R1 * st * cth + R2 * ct * sth)
      end

   elseif 4 <= r <= 8

      sth, cth = sincos(d.tht[Hj])

      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      θ0, dθ = if r == 4
         T5, T6 - T5
      elseif r == 5
         T6, T7 - T6
      elseif r == 6
         T7, T8 - T7
      elseif r == 7
         T8, T9 - T8
      else
         T9, T5 - T9 + 2.0π
      end

      cfac = -αt * dθ

      @inbounds for i in eachindex(out, t)
         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         out[i] = cfac * (R1 * st * cth + R2 * ct * sth)
      end

   else

      throw(ArgumentError("dgamx! local region must be 1:8; got r=$r, reg=$reg"))

   end

   return nothing

end

"""
    dgamy(d::ellipseNh, t::Float64, k::Int) -> Float64

Derivative of the second coordinate of the boundary parametrization γ(t).

- With `k`: derivative of the **right boundary** of patch `k`
  (`∂/∂t mapy(d, 1, t, k)` for `t ∈ [-1,1]`).
"""
function dgamy(d::ellipseNh, t::Float64, k::Int)::Float64

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = 0.5 * (p.tk0 + p.tk1)

   ξ2 = muladd(αt, t, βt)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if r <= 3

      sth, cth = sincos(d.tht[1])

      R1 = d.R1[1]
      R2 = d.R2[1]

      θ0, dθ = if r == 1
         T1, T2 - T1
      elseif r == 2
         T2, T3 - T2
      else
         T3, T4 - T3
      end

      θ = muladd(ξ2, dθ, θ0)
      st, ct = sincos(θ)

      return αt * dθ * (-R1 * st * sth + R2 * ct * cth)

   elseif 4 <= r <= 8

      sth, cth = sincos(d.tht[Hj])

      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      θ0, dθ = if r == 4
         T5, T6 - T5
      elseif r == 5
         T6, T7 - T6
      elseif r == 6
         T7, T8 - T7
      elseif r == 7
         T8, T9 - T8
      else
         T9, T5 - T9 + 2.0π
      end

      θ = muladd(ξ2, dθ, θ0)
      st, ct = sincos(θ)

      return αt * dθ * (-R1 * st * sth + R2 * ct * cth)

   else

      throw(ArgumentError("dgamy local region must be 1:8; got r=$r, reg=$reg"))

   end

end

function dgamy!(out::StridedArray{Float64}, d::ellipseNh,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = 0.5 * (p.tk0 + p.tk1)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if r <= 3

      sth, cth = sincos(d.tht[1])

      R1 = d.R1[1]
      R2 = d.R2[1]

      θ0, dθ = if r == 1
         T1, T2 - T1
      elseif r == 2
         T2, T3 - T2
      else
         T3, T4 - T3
      end

      cfac = αt * dθ

      @inbounds for i in eachindex(out, t)
         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         out[i] = cfac * (-R1 * st * sth + R2 * ct * cth)
      end

   elseif 4 <= r <= 8

      sth, cth = sincos(d.tht[Hj])

      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      θ0, dθ = if r == 4
         T5, T6 - T5
      elseif r == 5
         T6, T7 - T6
      elseif r == 6
         T7, T8 - T7
      elseif r == 7
         T8, T9 - T8
      else
         T9, T5 - T9 + 2.0π
      end

      cfac = αt * dθ

      @inbounds for i in eachindex(out, t)
         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         out[i] = cfac * (-R1 * st * sth + R2 * ct * cth)
      end

   else

      throw(ArgumentError("dgamy! local region must be 1:8; got r=$r, reg=$reg"))

   end

   return nothing

end

function dgam(d::ellipseNh, t::Float64, k::Int)::Tuple{Float64,Float64}

   return dgamx(d, t, k), dgamy(d, t, k)

end

function dgam!(out::Vector{Float64}, d::ellipseNh, t::Float64, k::Int)

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = 0.5 * (p.tk0 + p.tk1)

   ξ2 = muladd(αt, t, βt)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if r <= 3

      sth, cth = sincos(d.tht[1])

      R1 = d.R1[1]
      R2 = d.R2[1]

      θ0, dθ = if r == 1
         T1, T2 - T1
      elseif r == 2
         T2, T3 - T2
      else
         T3, T4 - T3
      end

      θ = muladd(ξ2, dθ, θ0)
      st, ct = sincos(θ)

      cfac = αt * dθ

      out[1] = -cfac * (R1 * st * cth + R2 * ct * sth)
      out[2] = cfac * (-R1 * st * sth + R2 * ct * cth)

   elseif 4 <= r <= 8

      sth, cth = sincos(d.tht[Hj])

      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      θ0, dθ = if r == 4
         T5, T6 - T5
      elseif r == 5
         T6, T7 - T6
      elseif r == 6
         T7, T8 - T7
      elseif r == 7
         T8, T9 - T8
      else
         T9, T5 - T9 + 2.0π
      end

      θ = muladd(ξ2, dθ, θ0)
      st, ct = sincos(θ)

      cfac = αt * dθ

      out[1] = -cfac * (R1 * st * cth + R2 * ct * sth)
      out[2] = cfac * (-R1 * st * sth + R2 * ct * cth)

   else

      throw(ArgumentError("dgam! local region must be 1:8; got r=$r, reg=$reg"))

   end

   return nothing

end

function dgam!(out::Matrix{Float64}, d::ellipseNh, t::Vector{Float64}, k::Int)

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = 0.5 * (p.tk0 + p.tk1)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if r <= 3

      sth, cth = sincos(d.tht[1])

      R1 = d.R1[1]
      R2 = d.R2[1]

      θ0, dθ = if r == 1
         T1, T2 - T1
      elseif r == 2
         T2, T3 - T2
      else
         T3, T4 - T3
      end

      cfac = αt * dθ

      @inbounds for i in eachindex(t)

         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         out[1, i] = -cfac * (R1 * st * cth + R2 * ct * sth)
         out[2, i] = cfac * (-R1 * st * sth + R2 * ct * cth)

      end

   elseif 4 <= r <= 8

      sth, cth = sincos(d.tht[Hj])

      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      θ0, dθ = if r == 4
         T5, T6 - T5
      elseif r == 5
         T6, T7 - T6
      elseif r == 6
         T7, T8 - T7
      elseif r == 7
         T8, T9 - T8
      else
         T9, T5 - T9 + 2.0π
      end

      cfac = αt * dθ

      @inbounds for i in eachindex(t)

         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         out[1, i] = -cfac * (R1 * st * cth + R2 * ct * sth)
         out[2, i] = cfac * (-R1 * st * sth + R2 * ct * cth)

      end

   else

      throw(ArgumentError("dgam!(ellipseNh): local region must be 1:8; got r=$r, reg=$reg"))

   end

   return nothing

end

function gamp(d::ellipseNh, t::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = 0.5 * (p.tk0 + p.tk1)

   ξ2 = muladd(αt, t, βt)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if r <= 3

      sth, cth = sincos(d.tht[1])

      R1 = d.R1[1]
      R2 = d.R2[1]

      θ0, dθ = if r == 1
         T1, T2 - T1
      elseif r == 2
         T2, T3 - T2
      else
         T3, T4 - T3
      end

      θ = muladd(ξ2, dθ, θ0)
      st, ct = sincos(θ)

      cfac = αt * dθ

      dx = -cfac * (R1 * st * cth + R2 * ct * sth)
      dy =  cfac * (-R1 * st * sth + R2 * ct * cth)

      return dy, -dx

   elseif 4 <= r <= 8

      sth, cth = sincos(d.tht[Hj])

      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      θ0, dθ = if r == 4
         T5, T6 - T5
      elseif r == 5
         T6, T7 - T6
      elseif r == 6
         T7, T8 - T7
      elseif r == 7
         T8, T9 - T8
      else
         T9, T5 - T9 + 2.0 * pi
      end

      θ = muladd(ξ2, dθ, θ0)
      st, ct = sincos(θ)

      cfac = αt * dθ

      dx = -cfac * (R1 * st * cth + R2 * ct * sth)
      dy =  cfac * (-R1 * st * sth + R2 * ct * cth)

      return -dy, dx

   else

      throw(ArgumentError("gamp local region must be 1:8; got r=$r, reg=$reg"))

   end

end

function gamp!(out::Matrix{Float64}, d::ellipseNh, t::Vector{Float64}, k::Int)

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = 0.5 * (p.tk0 + p.tk1)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if r <= 3

      σ = 1.0

      sth, cth = sincos(d.tht[1])

      R1 = d.R1[1]
      R2 = d.R2[1]

      θ0, dθ = if r == 1
         T1, T2 - T1
      elseif r == 2
         T2, T3 - T2
      else
         T3, T4 - T3
      end

      cfac = αt * dθ

      @inbounds for i in eachindex(t)

         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         dx = -cfac * (R1 * st * cth + R2 * ct * sth)
         dy =  cfac * (-R1 * st * sth + R2 * ct * cth)

         out[1, i] = σ * dy
         out[2, i] = σ * (-dx)

      end

   elseif 4 <= r <= 8

      σ = -1.0

      sth, cth = sincos(d.tht[Hj])

      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      θ0, dθ = if r == 4
         T5, T6 - T5
      elseif r == 5
         T6, T7 - T6
      elseif r == 6
         T7, T8 - T7
      elseif r == 7
         T8, T9 - T8
      else
         T9, T5 - T9 + 2.0 * pi
      end

      cfac = αt * dθ

      @inbounds for i in eachindex(t)

         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         dx = -cfac * (R1 * st * cth + R2 * ct * sth)
         dy =  cfac * (-R1 * st * sth + R2 * ct * cth)

         out[1, i] = σ * dy
         out[2, i] = σ * (-dx)

      end

   else

      throw(ArgumentError("gamp! local region must be 1:8; got r=$r, reg=$reg"))

   end

   return nothing

end

function nu!(out::Matrix{Float64}, d::ellipseNh, t::Vector{Float64}, k::Int)

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = 0.5 * (p.tk0 + p.tk1)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if r <= 3

      σ = 1.0

      sth, cth = sincos(d.tht[1])

      R1 = d.R1[1]
      R2 = d.R2[1]

      θ0, dθ = if r == 1
         T1, T2 - T1
      elseif r == 2
         T2, T3 - T2
      else
         T3, T4 - T3
      end

      cfac = αt * dθ

      @inbounds for i in eachindex(t)

         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         dx = -cfac * (R1 * st * cth + R2 * ct * sth)
         dy =  cfac * (-R1 * st * sth + R2 * ct * cth)

         gp1 = σ * dy
         gp2 = σ * (-dx)

         S = hypot(gp1, gp2)

         out[1, i] = gp1 / S
         out[2, i] = gp2 / S

      end

   elseif 4 <= r <= 8

      σ = -1.0

      sth, cth = sincos(d.tht[Hj])

      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      θ0, dθ = if r == 4
         T5, T6 - T5
      elseif r == 5
         T6, T7 - T6
      elseif r == 6
         T7, T8 - T7
      elseif r == 7
         T8, T9 - T8
      else
         T9, T5 - T9 + 2.0 * pi
      end

      cfac = αt * dθ

      @inbounds for i in eachindex(t)

         ξ2 = muladd(αt, t[i], βt)
         θ = muladd(ξ2, dθ, θ0)

         st, ct = sincos(θ)

         dx = -cfac * (R1 * st * cth + R2 * ct * sth)
         dy =  cfac * (-R1 * st * sth + R2 * ct * cth)

         gp1 = σ * dy
         gp2 = σ * (-dx)

         S = hypot(gp1, gp2)

         out[1, i] = gp1 / S
         out[2, i] = gp2 / S

      end

   else

      throw(ArgumentError("nu! local region must be 1:8; got r=$r, reg=$reg"))

   end

   return nothing

end

# The following function is specific to domains with nₕ >= 1.
function bdno(d::ellipseNh, k::Int)::Int

   reg = d.pths[k].reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   return r <= 3 ? 1 : H + 1

end

# This function finds s such that γ_l(t) = γ_k(s),
# allowing s to lie outside [-1, 1].
#
# For ellipseNh, this is only geometrically meaningful when patches l and k
# lie on the same boundary component:
#   outer component: local regions 1:3, for any hole block H
#   hole H component: local regions 4:8, within the same hole block H
function bdinv(d::ellipseNh, t::Float64, l::Int, k::Int)::Float64

   @inbounds begin

      pl = d.pths[l]
      pk = d.pths[k]

      regl = pl.reg
      regk = pk.reg

      Hl = cld(regl, 8)
      Hk = cld(regk, 8)

      rl = regl - 8 * (Hl - 1)
      rk = regk - 8 * (Hk - 1)

      outer_l = rl <= 3
      outer_k = rk <= 3

      if outer_l != outer_k
         throw(ArgumentError(
            "bdinv ellipseNh requires l and k on the same boundary type; got reg(l)=$regl, reg(k)=$regk"
         ))
      end

      if !outer_l && Hl != Hk
         throw(ArgumentError(
            "bdinv ellipseNh requires hole patches on the same hole; got H(l)=$Hl, H(k)=$Hk"
         ))
      end

      # Map t on patch l from [-1,1] to ξ_l ∈ [pl.tk0, pl.tk1].
      ξl = muladd(0.5 * (pl.tk1 - pl.tk0), t, 0.5 * (pl.tk0 + pl.tk1))

      if regl == regk
         return xi_inv(ξl, pk.tk0, pk.tk1)
      end

      # ------------------------------------------------------------
      # Source angle θ = Δθ_l * ξ_l + θ0_l.
      # ------------------------------------------------------------
      θ0_l, Δθ_l = if rl == 1

         d.T[1, Hl], d.T[2, Hl] - d.T[1, Hl]

      elseif rl == 2

         d.T[2, Hl], d.T[3, Hl] - d.T[2, Hl]

      elseif rl == 3

         d.T[3, Hl], d.T[4, Hl] - d.T[3, Hl]

      elseif rl == 4

         d.T[5, Hl], d.T[6, Hl] - d.T[5, Hl]

      elseif rl == 5

         d.T[6, Hl], d.T[7, Hl] - d.T[6, Hl]

      elseif rl == 6

         d.T[7, Hl], d.T[8, Hl] - d.T[7, Hl]

      elseif rl == 7

         d.T[8, Hl], d.T[9, Hl] - d.T[8, Hl]

      elseif rl == 8

         d.T[9, Hl], d.T[5, Hl] - d.T[9, Hl] + 2.0π

      else

         throw(ArgumentError("bdinv ellipseNh expects source local region ∈ 1:8; got reg=$regl"))

      end

      θ = muladd(Δθ_l, ξl, θ0_l)

      # ------------------------------------------------------------
      # Target angle map θ = Δθ_k * ξ_k + θ0_k.
      # ------------------------------------------------------------
      θ0_k, Δθ_k = if rk == 1

         d.T[1, Hk], d.T[2, Hk] - d.T[1, Hk]

      elseif rk == 2

         d.T[2, Hk], d.T[3, Hk] - d.T[2, Hk]

      elseif rk == 3

         d.T[3, Hk], d.T[4, Hk] - d.T[3, Hk]

      elseif rk == 4

         d.T[5, Hk], d.T[6, Hk] - d.T[5, Hk]

      elseif rk == 5

         d.T[6, Hk], d.T[7, Hk] - d.T[6, Hk]

      elseif rk == 6

         d.T[7, Hk], d.T[8, Hk] - d.T[7, Hk]

      elseif rk == 7

         d.T[8, Hk], d.T[9, Hk] - d.T[8, Hk]

      elseif rk == 8

         d.T[9, Hk], d.T[5, Hk] - d.T[9, Hk] + 2.0π

      else

         throw(ArgumentError("bdinv ellipseNh expects target local region ∈ 1:8; got reg=$regk"))

      end

      # Shift θ by multiples of 2π so it is closest to the target interval.
      center_k = θ0_k + 0.5 * Δθ_k
      m = round((center_k - θ) / (2.0π))
      θs = θ + 2.0π * m

      # Convert angle to target ξ_k, then ξ_k to local coordinate s.
      ξk = (θs - θ0_k) / Δθ_k

      return xi_inv(ξk, pk.tk0, pk.tk1)

   end

end

"""
    DLP!(out, d::ellipseNh, t, tau, k, x, G, GP)

Double-layer kernel on a boundary patch of ellipseNh.

For the diagonal/self-boundary case, each boundary component is an ellipse,
so the expression simplifies to

outer boundary pieces, local regions 1:3:

    K =  (ht*ak*R1*R2/4) / ((R1*sin(tt))^2 + (R2*cos(tt))^2)

hole boundary pieces, local regions 4:8:

    K = -(ht*ak*R1*R2/4) / ((R1*sin(tt))^2 + (R2*cos(tt))^2)

where

    tt = ak * (t̂ + τ̂)/2 + bk.
"""
function DLP!(out::StridedArray{Float64}, d::ellipseNh, t::Float64,
   tau::StridedArray{Float64}, k::Int, x::Vector{Float64},
   G::Vector{Float64}, GP::Vector{Float64})

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = 0.5 * (p.tk0 + p.tk1)

   tm = muladd(αt, t, βt)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   ak, bk = if r == 1

      d.T[2, H] - d.T[1, H], d.T[1, H]

   elseif r == 2

      d.T[3, H] - d.T[2, H], d.T[2, H]

   elseif r == 3

      d.T[4, H] - d.T[3, H], d.T[3, H]

   elseif r == 4

      d.T[6, H] - d.T[5, H], d.T[5, H]

   elseif r == 5

      d.T[7, H] - d.T[6, H], d.T[6, H]

   elseif r == 6

      d.T[8, H] - d.T[7, H], d.T[7, H]

   elseif r == 7

      d.T[9, H] - d.T[8, H], d.T[8, H]

   elseif r == 8

      d.T[5, H] - d.T[9, H] + 2.0π, d.T[9, H]

   else

      throw(ArgumentError("DLP!(ellipseNh): local region must be 1:8; got r=$r, reg=$reg"))

   end

   if r <= 3

      R1 = d.R1[1]
      R2 = d.R2[1]

      pref = 0.25 * ht * ak * R1 * R2
      halfak = 0.5 * ak

      @inbounds for i in eachindex(out, tau)

         τm = muladd(αt, tau[i], βt)

         tt = muladd(halfak, tm + τm, bk)

         st, ct = sincos(tt)

         den = (R1 * st) * (R1 * st) + (R2 * ct) * (R2 * ct)

         out[i] = pref / den

      end

   else

      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      pref = -0.25 * ht * ak * R1 * R2
      halfak = 0.5 * ak

      @inbounds for i in eachindex(out, tau)

         τm = muladd(αt, tau[i], βt)

         tt = muladd(halfak, tm + τm, bk)

         st, ct = sincos(tt)

         den = (R1 * st) * (R1 * st) + (R2 * ct) * (R2 * ct)

         out[i] = pref / den

      end

   end

   return nothing

end

#-----------------------

function dergam(d::ellipseNh, t::Float64, j::Int)::Tuple{Float64,Float64,Float64,Float64}

   st, ct = sincos(t)
   sj, cj = sincos(d.tht[j])

   gx = d.A[j] + d.R1[j] * ct * cj - d.R2[j] * st * sj
   gy = d.B[j] + d.R1[j] * ct * sj + d.R2[j] * st * cj

   dgx = -d.R1[j] * st * cj - d.R2[j] * ct * sj
   dgy = -d.R1[j] * st * sj + d.R2[j] * ct * cj

   return gx, gy, dgx, dgy

end

function Dmap!(out::StridedArray{Float64}, d::ellipseNh,
   u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   αc = 0.5 * (p.ck1 - p.ck0)
   βc = 0.5 * (p.ck0 + p.ck1)

   αt = 0.5 * (p.tk1 - p.tk0)
   βt = 0.5 * (p.tk0 + p.tk1)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1  = d.T[1, H]
   T2  = d.T[2, H]
   T3  = d.T[3, H]
   T4  = d.T[4, H]
   T5  = d.T[5, H]
   T6  = d.T[6, H]
   T7  = d.T[7, H]
   T8  = d.T[8, H]
   T9  = d.T[9, H]
   T10 = d.T[10, H]
   T11 = d.T[11, H]
   T12 = d.T[12, H]
   T13 = d.T[13, H]

   # Outer ellipse constants.
   s1, c1 = sincos(d.tht[1])
   A1 = d.A[1]
   B1 = d.B[1]
   R11 = d.R1[1]
   R12 = d.R2[1]

   # Hole ellipse constants.
   sj, cj = sincos(d.tht[Hj])
   Aj = d.A[Hj]
   Bj = d.B[Hj]
   Rj1 = d.R1[Hj]
   Rj2 = d.R2[Hj]

   if r == 1

      dθ = T2 - T1

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      XRx = 0.5 * (gjx5 + g1x2)
      XRy = 0.5 * (gjy5 + g1y2)

      Xpx = XRx - T12
      Xpy = XRy - T13

      @inbounds for i in eachindex(out, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         Xx = muladd(xi2, Xpx, T12)
         Xy = muladd(xi2, Xpy, T13)

         θ = muladd(xi2, dθ, T1)
         st, ct = sincos(θ)

         Yx = A1 + R11 * ct * c1 - R12 * st * s1
         Yy = B1 + R11 * ct * s1 + R12 * st * c1

         dYx = dθ * (-R11 * st * c1 - R12 * ct * s1)
         dYy = dθ * (-R11 * st * s1 + R12 * ct * c1)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

         out[i] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 2

      dθo = T3 - T2
      dθh = T6 - T5

      @inbounds for i in eachindex(out, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         θo = muladd(xi2, dθo, T2)
         sto, cto = sincos(θo)

         gx = A1 + R11 * cto * c1 - R12 * sto * s1
         gy = B1 + R11 * cto * s1 + R12 * sto * c1

         dgx = -R11 * sto * c1 - R12 * cto * s1
         dgy = -R11 * sto * s1 + R12 * cto * c1

         θh = muladd(xi2, dθh, T5)
         sth, cth = sincos(θh)

         gjx = Aj + Rj1 * cth * cj - Rj2 * sth * sj
         gjy = Bj + Rj1 * cth * sj + Rj2 * sth * cj

         dgjx = -Rj1 * sth * cj - Rj2 * cth * sj
         dgjy = -Rj1 * sth * sj + Rj2 * cth * cj

         dux = 0.5 * αc * (gx - gjx)
         duy = 0.5 * αc * (gy - gjy)

         dvx = 0.5 * αt * ((1.0 - xi1) * dθh * dgjx + (1.0 + xi1) * dθo * dgx)
         dvy = 0.5 * αt * ((1.0 - xi1) * dθh * dgjy + (1.0 + xi1) * dθo * dgy)

         out[i] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 3

      dθ = T4 - T3

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      XLx = 0.5 * (gjx6 + g1x3)
      XLy = 0.5 * (gjy6 + g1y3)

      Xpx = T10 - XLx
      Xpy = T11 - XLy

      @inbounds for i in eachindex(out, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         Xx = muladd(xi2, Xpx, XLx)
         Xy = muladd(xi2, Xpy, XLy)

         θ = muladd(xi2, dθ, T3)
         st, ct = sincos(θ)

         Yx = A1 + R11 * ct * c1 - R12 * st * s1
         Yy = B1 + R11 * ct * s1 + R12 * st * c1

         dYx = dθ * (-R11 * st * c1 - R12 * ct * s1)
         dYy = dθ * (-R11 * st * s1 + R12 * ct * c1)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

         out[i] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 4

      dθo = T3 - T2
      dθh = T6 - T5

      @inbounds for i in eachindex(out, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         θo = muladd(xi2, dθo, T2)
         sto, cto = sincos(θo)

         gx = A1 + R11 * cto * c1 - R12 * sto * s1
         gy = B1 + R11 * cto * s1 + R12 * sto * c1

         dgx = -R11 * sto * c1 - R12 * cto * s1
         dgy = -R11 * sto * s1 + R12 * cto * c1

         θh = muladd(xi2, dθh, T5)
         sth, cth = sincos(θh)

         gjx = Aj + Rj1 * cth * cj - Rj2 * sth * sj
         gjy = Bj + Rj1 * cth * sj + Rj2 * sth * cj

         dgjx = -Rj1 * sth * cj - Rj2 * cth * sj
         dgjy = -Rj1 * sth * sj + Rj2 * cth * cj

         dux = 0.5 * αc * (gjx - gx)
         duy = 0.5 * αc * (gjy - gy)

         dvx = 0.5 * αt * ((1.0 - xi1) * dθo * dgx + (1.0 + xi1) * dθh * dgjx)
         dvy = 0.5 * αt * ((1.0 - xi1) * dθo * dgy + (1.0 + xi1) * dθh * dgjy)

         out[i] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 5

      dθ = T7 - T6

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      XLx = 0.5 * (gjx6 + g1x3)
      XLy = 0.5 * (gjy6 + g1y3)

      Xpx = T10 - XLx
      Xpy = T11 - XLy

      @inbounds for i in eachindex(out, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         Xx = muladd(xi2, Xpx, XLx)
         Xy = muladd(xi2, Xpy, XLy)

         θ = muladd(xi2, dθ, T6)
         st, ct = sincos(θ)

         Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

         out[i] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 6

      dθ = T8 - T7

      Xpx = A1 - T10
      Xpy = B1 - T11

      @inbounds for i in eachindex(out, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         Xx = muladd(xi2, Xpx, T10)
         Xy = muladd(xi2, Xpy, T11)

         θ = muladd(xi2, dθ, T7)
         st, ct = sincos(θ)

         Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

         out[i] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 7

      dθ = T9 - T8

      Xpx = T12 - A1
      Xpy = T13 - B1

      @inbounds for i in eachindex(out, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         Xx = muladd(xi2, Xpx, A1)
         Xy = muladd(xi2, Xpy, B1)

         θ = muladd(xi2, dθ, T8)
         st, ct = sincos(θ)

         Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

         out[i] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 8

      dθ = T5 - T9 + 2.0π

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      XRx = 0.5 * (gjx5 + g1x2)
      XRy = 0.5 * (gjy5 + g1y2)

      Xpx = XRx - T12
      Xpy = XRy - T13

      @inbounds for i in eachindex(out, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         Xx = muladd(xi2, Xpx, T12)
         Xy = muladd(xi2, Xpy, T13)

         θ = muladd(xi2, dθ, T9)
         st, ct = sincos(θ)

         Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

         out[i] = abs(dux * dvy - dvx * duy)

      end

   else

      throw(ArgumentError("Dmap!(ellipseNh): local region must be 1:8; got r=$r, reg=$reg"))

   end

   return nothing

end

"""
A combination of mapxy! and Dmap! function (No allocations!)
The purpose of this function is to reduce computations  
related to cosine and sine's.  
"""
function mapxy_Dmap!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
   DJ::StridedArray{Float64}, d::ellipseNh, u::StridedArray{Float64},
   v::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   αc = 0.5 * (p.ck1 - p.ck0)
   βc = 0.5 * (p.ck0 + p.ck1)

   αt = 0.5 * (p.tk1 - p.tk0)
   βt = 0.5 * (p.tk0 + p.tk1)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1  = d.T[1, H]
   T2  = d.T[2, H]
   T3  = d.T[3, H]
   T4  = d.T[4, H]
   T5  = d.T[5, H]
   T6  = d.T[6, H]
   T7  = d.T[7, H]
   T8  = d.T[8, H]
   T9  = d.T[9, H]
   T10 = d.T[10, H]
   T11 = d.T[11, H]
   T12 = d.T[12, H]
   T13 = d.T[13, H]

   # Outer ellipse constants.
   s1, c1 = sincos(d.tht[1])
   A1 = d.A[1]
   B1 = d.B[1]
   R11 = d.R1[1]
   R12 = d.R2[1]

   # Hole ellipse constants.
   sj, cj = sincos(d.tht[Hj])
   Aj = d.A[Hj]
   Bj = d.B[Hj]
   Rj1 = d.R1[Hj]
   Rj2 = d.R2[Hj]

   if r == 1

      dθ = T2 - T1

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      XRx = 0.5 * (gjx5 + g1x2)
      XRy = 0.5 * (gjy5 + g1y2)

      Xpx = XRx - T12
      Xpy = XRy - T13

      @inbounds for i in eachindex(Zx, Zy, DJ, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         Xx = muladd(xi2, Xpx, T12)
         Xy = muladd(xi2, Xpy, T13)

         θ = muladd(xi2, dθ, T1)
         st, ct = sincos(θ)

         Yx = A1 + R11 * ct * c1 - R12 * st * s1
         Yy = B1 + R11 * ct * s1 + R12 * st * c1

         dYx = dθ * (-R11 * st * c1 - R12 * ct * s1)
         dYy = dθ * (-R11 * st * s1 + R12 * ct * c1)

         Zx[i] = muladd(xi1, Yx - Xx, Xx)
         Zy[i] = muladd(xi1, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

         DJ[i] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 2

      dθo = T3 - T2
      dθh = T6 - T5

      @inbounds for i in eachindex(Zx, Zy, DJ, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         θo = muladd(xi2, dθo, T2)
         sto, cto = sincos(θo)

         gx = A1 + R11 * cto * c1 - R12 * sto * s1
         gy = B1 + R11 * cto * s1 + R12 * sto * c1

         dgx = -R11 * sto * c1 - R12 * cto * s1
         dgy = -R11 * sto * s1 + R12 * cto * c1

         θh = muladd(xi2, dθh, T5)
         sth, cth = sincos(θh)

         gjx = Aj + Rj1 * cth * cj - Rj2 * sth * sj
         gjy = Bj + Rj1 * cth * sj + Rj2 * sth * cj

         dgjx = -Rj1 * sth * cj - Rj2 * cth * sj
         dgjy = -Rj1 * sth * sj + Rj2 * cth * cj

         Xx = 0.5 * (gx + gjx)
         Xy = 0.5 * (gy + gjy)

         Zx[i] = muladd(xi1, gx - Xx, Xx)
         Zy[i] = muladd(xi1, gy - Xy, Xy)

         dux = 0.5 * αc * (gx - gjx)
         duy = 0.5 * αc * (gy - gjy)

         dvx = 0.5 * αt * ((1.0 - xi1) * dθh * dgjx + (1.0 + xi1) * dθo * dgx)
         dvy = 0.5 * αt * ((1.0 - xi1) * dθh * dgjy + (1.0 + xi1) * dθo * dgy)

         DJ[i] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 3

      dθ = T4 - T3

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      XLx = 0.5 * (gjx6 + g1x3)
      XLy = 0.5 * (gjy6 + g1y3)

      Xpx = T10 - XLx
      Xpy = T11 - XLy

      @inbounds for i in eachindex(Zx, Zy, DJ, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         Xx = muladd(xi2, Xpx, XLx)
         Xy = muladd(xi2, Xpy, XLy)

         θ = muladd(xi2, dθ, T3)
         st, ct = sincos(θ)

         Yx = A1 + R11 * ct * c1 - R12 * st * s1
         Yy = B1 + R11 * ct * s1 + R12 * st * c1

         dYx = dθ * (-R11 * st * c1 - R12 * ct * s1)
         dYy = dθ * (-R11 * st * s1 + R12 * ct * c1)

         Zx[i] = muladd(xi1, Yx - Xx, Xx)
         Zy[i] = muladd(xi1, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

         DJ[i] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 4

      dθo = T3 - T2
      dθh = T6 - T5

      @inbounds for i in eachindex(Zx, Zy, DJ, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         θo = muladd(xi2, dθo, T2)
         sto, cto = sincos(θo)

         gx = A1 + R11 * cto * c1 - R12 * sto * s1
         gy = B1 + R11 * cto * s1 + R12 * sto * c1

         dgx = -R11 * sto * c1 - R12 * cto * s1
         dgy = -R11 * sto * s1 + R12 * cto * c1

         θh = muladd(xi2, dθh, T5)
         sth, cth = sincos(θh)

         gjx = Aj + Rj1 * cth * cj - Rj2 * sth * sj
         gjy = Bj + Rj1 * cth * sj + Rj2 * sth * cj

         dgjx = -Rj1 * sth * cj - Rj2 * cth * sj
         dgjy = -Rj1 * sth * sj + Rj2 * cth * cj

         Xx = 0.5 * (gx + gjx)
         Xy = 0.5 * (gy + gjy)

         Zx[i] = muladd(xi1, gjx - Xx, Xx)
         Zy[i] = muladd(xi1, gjy - Xy, Xy)

         dux = 0.5 * αc * (gjx - gx)
         duy = 0.5 * αc * (gjy - gy)

         dvx = 0.5 * αt * ((1.0 - xi1) * dθo * dgx + (1.0 + xi1) * dθh * dgjx)
         dvy = 0.5 * αt * ((1.0 - xi1) * dθo * dgy + (1.0 + xi1) * dθh * dgjy)

         DJ[i] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 5

      dθ = T7 - T6

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      XLx = 0.5 * (gjx6 + g1x3)
      XLy = 0.5 * (gjy6 + g1y3)

      Xpx = T10 - XLx
      Xpy = T11 - XLy

      @inbounds for i in eachindex(Zx, Zy, DJ, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         Xx = muladd(xi2, Xpx, XLx)
         Xy = muladd(xi2, Xpy, XLy)

         θ = muladd(xi2, dθ, T6)
         st, ct = sincos(θ)

         Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         Zx[i] = muladd(xi1, Yx - Xx, Xx)
         Zy[i] = muladd(xi1, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

         DJ[i] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 6

      dθ = T8 - T7

      Xpx = A1 - T10
      Xpy = B1 - T11

      @inbounds for i in eachindex(Zx, Zy, DJ, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         Xx = muladd(xi2, Xpx, T10)
         Xy = muladd(xi2, Xpy, T11)

         θ = muladd(xi2, dθ, T7)
         st, ct = sincos(θ)

         Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         Zx[i] = muladd(xi1, Yx - Xx, Xx)
         Zy[i] = muladd(xi1, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

         DJ[i] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 7

      dθ = T9 - T8

      Xpx = T12 - A1
      Xpy = T13 - B1

      @inbounds for i in eachindex(Zx, Zy, DJ, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         Xx = muladd(xi2, Xpx, A1)
         Xy = muladd(xi2, Xpy, B1)

         θ = muladd(xi2, dθ, T8)
         st, ct = sincos(θ)

         Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         Zx[i] = muladd(xi1, Yx - Xx, Xx)
         Zy[i] = muladd(xi1, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

         DJ[i] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 8

      dθ = T5 - T9 + 2.0π

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      XRx = 0.5 * (gjx5 + g1x2)
      XRy = 0.5 * (gjy5 + g1y2)

      Xpx = XRx - T12
      Xpy = XRy - T13

      @inbounds for i in eachindex(Zx, Zy, DJ, u, v)

         xi1 = muladd(αc, u[i], βc)
         xi2 = muladd(αt, v[i], βt)

         Xx = muladd(xi2, Xpx, T12)
         Xy = muladd(xi2, Xpy, T13)

         θ = muladd(xi2, dθ, T9)
         st, ct = sincos(θ)

         Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         Zx[i] = muladd(xi1, Yx - Xx, Xx)
         Zy[i] = muladd(xi1, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

         DJ[i] = abs(dux * dvy - dvx * duy)

      end

   else

      throw(ArgumentError("mapxy_Dmap!(ellipseNh): local region must be 1:8; got r=$r, reg=$reg"))

   end

   return nothing

end

function chk_map(d::ellipseNh; n::Int=32, tol::Float64=5e-14)

   Iex = 0.0

   Iex += d.R1[1] * d.R2[1] * ((d.R1[1]^2 + d.R2[1]^2) / 4 + d.A[1]^2 + d.B[1]^2)

   @inbounds for i in 2:(d.nh+1)

      Iex -= d.R1[i] * d.R2[i] * ((d.R1[i]^2 + d.R2[i]^2) / 4 + d.A[i]^2 + d.B[i]^2)

   end

   Iex *= π

   return chkmap_geom(d, Iex; n=n, tol=tol)

end

"""
    jinvmap(d::ellipseNh, u::Float64, v::Float64, r::Int) -> tuple

Inverse Jacobian of the ellipseNh regional mapping at `(u,v) ∈ [-1,1]^2`.

Here `r` is the actual region number, not just the local region.
For ellipseNh,

    H  = cld(r, 8)
    rl = r - 8*(H - 1)

where `H` is the hole block and `rl ∈ 1:8` is the local region.

Returns `(J11, J12, J21, J22)` for the inverse of

    [dux  dvx
     duy  dvy]

with `hc = ht = 1`.
"""
function jinvmap(d::ellipseNh, u::Float64, v::Float64, r::Int)

   @inbounds begin

      αc = 0.5
      βc = 0.5
      αt = 0.5
      βt = 0.5

      xi1 = muladd(αc, u, βc)
      xi2 = muladd(αt, v, βt)

      H = cld(r, 8)
      rl = r - 8 * (H - 1)

      Hj = H + 1

      T1  = d.T[1, H]
      T2  = d.T[2, H]
      T3  = d.T[3, H]
      T4  = d.T[4, H]
      T5  = d.T[5, H]
      T6  = d.T[6, H]
      T7  = d.T[7, H]
      T8  = d.T[8, H]
      T9  = d.T[9, H]
      T10 = d.T[10, H]
      T11 = d.T[11, H]
      T12 = d.T[12, H]
      T13 = d.T[13, H]

      # Outer ellipse constants.
      s1, c1 = sincos(d.tht[1])
      A1 = d.A[1]
      B1 = d.B[1]
      R11 = d.R1[1]
      R12 = d.R2[1]

      # Hole ellipse constants.
      sj, cj = sincos(d.tht[Hj])
      Aj = d.A[Hj]
      Bj = d.B[Hj]
      Rj1 = d.R1[Hj]
      Rj2 = d.R2[Hj]

      if rl == 1

         dθ = T2 - T1

         st5, ct5 = sincos(T5)
         st2, ct2 = sincos(T2)

         gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
         gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

         g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
         g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

         XRx = 0.5 * (gjx5 + g1x2)
         XRy = 0.5 * (gjy5 + g1y2)

         Xpx = XRx - T12
         Xpy = XRy - T13

         Xx = muladd(xi2, Xpx, T12)
         Xy = muladd(xi2, Xpy, T13)

         θ = muladd(xi2, dθ, T1)
         st, ct = sincos(θ)

         Yx = A1 + R11 * ct * c1 - R12 * st * s1
         Yy = B1 + R11 * ct * s1 + R12 * st * c1

         dYx = dθ * (-R11 * st * c1 - R12 * ct * s1)
         dYy = dθ * (-R11 * st * s1 + R12 * ct * c1)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

      elseif rl == 2

         dθo = T3 - T2
         dθh = T6 - T5

         θo = muladd(xi2, dθo, T2)
         sto, cto = sincos(θo)

         gx = A1 + R11 * cto * c1 - R12 * sto * s1
         gy = B1 + R11 * cto * s1 + R12 * sto * c1

         dgx = -R11 * sto * c1 - R12 * cto * s1
         dgy = -R11 * sto * s1 + R12 * cto * c1

         θh = muladd(xi2, dθh, T5)
         sth, cth = sincos(θh)

         gjx = Aj + Rj1 * cth * cj - Rj2 * sth * sj
         gjy = Bj + Rj1 * cth * sj + Rj2 * sth * cj

         dgjx = -Rj1 * sth * cj - Rj2 * cth * sj
         dgjy = -Rj1 * sth * sj + Rj2 * cth * cj

         dux = 0.5 * αc * (gx - gjx)
         duy = 0.5 * αc * (gy - gjy)

         dvx = 0.5 * αt * ((1.0 - xi1) * dθh * dgjx + (1.0 + xi1) * dθo * dgx)
         dvy = 0.5 * αt * ((1.0 - xi1) * dθh * dgjy + (1.0 + xi1) * dθo * dgy)

      elseif rl == 3

         dθ = T4 - T3

         st6, ct6 = sincos(T6)
         st3, ct3 = sincos(T3)

         gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
         gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

         g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
         g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

         XLx = 0.5 * (gjx6 + g1x3)
         XLy = 0.5 * (gjy6 + g1y3)

         Xpx = T10 - XLx
         Xpy = T11 - XLy

         Xx = muladd(xi2, Xpx, XLx)
         Xy = muladd(xi2, Xpy, XLy)

         θ = muladd(xi2, dθ, T3)
         st, ct = sincos(θ)

         Yx = A1 + R11 * ct * c1 - R12 * st * s1
         Yy = B1 + R11 * ct * s1 + R12 * st * c1

         dYx = dθ * (-R11 * st * c1 - R12 * ct * s1)
         dYy = dθ * (-R11 * st * s1 + R12 * ct * c1)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

      elseif rl == 4

         dθo = T3 - T2
         dθh = T6 - T5

         θo = muladd(xi2, dθo, T2)
         sto, cto = sincos(θo)

         gx = A1 + R11 * cto * c1 - R12 * sto * s1
         gy = B1 + R11 * cto * s1 + R12 * sto * c1

         dgx = -R11 * sto * c1 - R12 * cto * s1
         dgy = -R11 * sto * s1 + R12 * cto * c1

         θh = muladd(xi2, dθh, T5)
         sth, cth = sincos(θh)

         gjx = Aj + Rj1 * cth * cj - Rj2 * sth * sj
         gjy = Bj + Rj1 * cth * sj + Rj2 * sth * cj

         dgjx = -Rj1 * sth * cj - Rj2 * cth * sj
         dgjy = -Rj1 * sth * sj + Rj2 * cth * cj

         dux = 0.5 * αc * (gjx - gx)
         duy = 0.5 * αc * (gjy - gy)

         dvx = 0.5 * αt * ((1.0 - xi1) * dθo * dgx + (1.0 + xi1) * dθh * dgjx)
         dvy = 0.5 * αt * ((1.0 - xi1) * dθo * dgy + (1.0 + xi1) * dθh * dgjy)

      elseif rl == 5

         dθ = T7 - T6

         st6, ct6 = sincos(T6)
         st3, ct3 = sincos(T3)

         gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
         gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

         g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
         g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

         XLx = 0.5 * (gjx6 + g1x3)
         XLy = 0.5 * (gjy6 + g1y3)

         Xpx = T10 - XLx
         Xpy = T11 - XLy

         Xx = muladd(xi2, Xpx, XLx)
         Xy = muladd(xi2, Xpy, XLy)

         θ = muladd(xi2, dθ, T6)
         st, ct = sincos(θ)

         Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

      elseif rl == 6

         dθ = T8 - T7

         Xpx = A1 - T10
         Xpy = B1 - T11

         Xx = muladd(xi2, Xpx, T10)
         Xy = muladd(xi2, Xpy, T11)

         θ = muladd(xi2, dθ, T7)
         st, ct = sincos(θ)

         Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

      elseif rl == 7

         dθ = T9 - T8

         Xpx = T12 - A1
         Xpy = T13 - B1

         Xx = muladd(xi2, Xpx, A1)
         Xy = muladd(xi2, Xpy, B1)

         θ = muladd(xi2, dθ, T8)
         st, ct = sincos(θ)

         Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

      elseif rl == 8

         dθ = T5 - T9 + 2.0π

         st5, ct5 = sincos(T5)
         st2, ct2 = sincos(T2)

         gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
         gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

         g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
         g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

         XRx = 0.5 * (gjx5 + g1x2)
         XRy = 0.5 * (gjy5 + g1y2)

         Xpx = XRx - T12
         Xpy = XRy - T13

         Xx = muladd(xi2, Xpx, T12)
         Xy = muladd(xi2, Xpy, T13)

         θ = muladd(xi2, dθ, T9)
         st, ct = sincos(θ)

         Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = αt * ((1.0 - xi1) * Xpx + xi1 * dYx)
         dvy = αt * ((1.0 - xi1) * Xpy + xi1 * dYy)

      else

         throw(ArgumentError("jinvmap(ellipseNh): local region must be 1:8; got r=$r"))

      end

   end

   detJ = dux * dvy - dvx * duy
   invdet = 1.0 / detJ

   J11 =  invdet * dvy
   J12 = -invdet * dvx
   J21 = -invdet * duy
   J22 =  invdet * dux

   return J11, J12, J21, J22

end

"""
  mapinv(d::ellipseNh, u::Float64, v::Float64, k::Int) -> Tuple{Float64,Float64}

Inverse mapping: given a physical point `(u,v)` on patch `k`, return
the reference coordinates `[t; s]` in `[-1,1]^2` such that `mapxy(d, t, s, k) = (u,v)`.
"""
function Xx(s::Float64, d::ellipseNh, r::Int)

   ξ = muladd(0.5, s, 0.5)

   H = cld(r, 8)
   rl = r - 8 * (H - 1)

   Hj = H + 1

   T2  = d.T[2, H]
   T3  = d.T[3, H]
   T5  = d.T[5, H]
   T6  = d.T[6, H]
   T10 = d.T[10, H]
   T12 = d.T[12, H]

   s1, c1 = sincos(d.tht[1])
   sj, cj = sincos(d.tht[Hj])

   A1 = d.A[1]
   R11 = d.R1[1]
   R12 = d.R2[1]

   Aj = d.A[Hj]
   Rj1 = d.R1[Hj]
   Rj2 = d.R2[Hj]

   if rl == 1 || rl == 8

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1

      XR = 0.5 * (gjx5 + g1x2)

      return muladd(ξ, XR - T12, T12)

   elseif rl == 2 || rl == 4

      θo = muladd(ξ, T3 - T2, T2)
      sto, cto = sincos(θo)

      θh = muladd(ξ, T6 - T5, T5)
      sth, cth = sincos(θh)

      gx = A1 + R11 * cto * c1 - R12 * sto * s1
      gjx = Aj + Rj1 * cth * cj - Rj2 * sth * sj

      return 0.5 * (gx + gjx)

   elseif rl == 3 || rl == 5

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1

      XL = 0.5 * (gjx6 + g1x3)

      return muladd(ξ, T10 - XL, XL)

   elseif rl == 6

      return muladd(ξ, A1 - T10, T10)

   elseif rl == 7

      return muladd(ξ, T12 - A1, A1)

   else

      throw(ArgumentError("Xx local region must be 1:8; got rl=$rl, r=$r"))

   end

end

function Xy(s::Float64, d::ellipseNh, r::Int)

   ξ = muladd(0.5, s, 0.5)

   H = cld(r, 8)
   rl = r - 8 * (H - 1)

   Hj = H + 1

   T2  = d.T[2, H]
   T3  = d.T[3, H]
   T5  = d.T[5, H]
   T6  = d.T[6, H]
   T11 = d.T[11, H]
   T13 = d.T[13, H]

   s1, c1 = sincos(d.tht[1])
   sj, cj = sincos(d.tht[Hj])

   B1 = d.B[1]
   R11 = d.R1[1]
   R12 = d.R2[1]

   Bj = d.B[Hj]
   Rj1 = d.R1[Hj]
   Rj2 = d.R2[Hj]

   if rl == 1 || rl == 8

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      YR = 0.5 * (gjy5 + g1y2)

      return muladd(ξ, YR - T13, T13)

   elseif rl == 2 || rl == 4

      θo = muladd(ξ, T3 - T2, T2)
      sto, cto = sincos(θo)

      θh = muladd(ξ, T6 - T5, T5)
      sth, cth = sincos(θh)

      gy = B1 + R11 * cto * s1 + R12 * sto * c1
      gjy = Bj + Rj1 * cth * sj + Rj2 * sth * cj

      return 0.5 * (gy + gjy)

   elseif rl == 3 || rl == 5

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      YL = 0.5 * (gjy6 + g1y3)

      return muladd(ξ, T11 - YL, YL)

   elseif rl == 6

      return muladd(ξ, B1 - T11, T11)

   elseif rl == 7

      return muladd(ξ, T13 - B1, B1)

   else

      throw(ArgumentError("Xy local region must be 1:8; got rl=$rl, r=$r"))

   end

end

function Yx(s::Float64, d::ellipseNh, r::Int)

   ξ = muladd(0.5, s, 0.5)

   H = cld(r, 8)
   rl = r - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if rl <= 3

      s1, c1 = sincos(d.tht[1])

      A = d.A[1]
      R1 = d.R1[1]
      R2 = d.R2[1]

      θ0, dθ = if rl == 1
         T1, T2 - T1
      elseif rl == 2
         T2, T3 - T2
      else
         T3, T4 - T3
      end

      θ = muladd(ξ, dθ, θ0)
      st, ct = sincos(θ)

      return A + R1 * ct * c1 - R2 * st * s1

   elseif 4 <= rl <= 8

      sj, cj = sincos(d.tht[Hj])

      A = d.A[Hj]
      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      θ0, dθ = if rl == 4
         T5, T6 - T5
      elseif rl == 5
         T6, T7 - T6
      elseif rl == 6
         T7, T8 - T7
      elseif rl == 7
         T8, T9 - T8
      else
         T9, T5 - T9 + 2.0π
      end

      θ = muladd(ξ, dθ, θ0)
      st, ct = sincos(θ)

      return A + R1 * ct * cj - R2 * st * sj

   else

      throw(ArgumentError("Yx local region must be 1:8; got rl=$rl, r=$r"))

   end

end

function Yy(s::Float64, d::ellipseNh, r::Int)

   ξ = muladd(0.5, s, 0.5)

   H = cld(r, 8)
   rl = r - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]

   if rl <= 3

      s1, c1 = sincos(d.tht[1])

      B = d.B[1]
      R1 = d.R1[1]
      R2 = d.R2[1]

      θ0, dθ = if rl == 1
         T1, T2 - T1
      elseif rl == 2
         T2, T3 - T2
      else
         T3, T4 - T3
      end

      θ = muladd(ξ, dθ, θ0)
      st, ct = sincos(θ)

      return B + R1 * ct * s1 + R2 * st * c1

   elseif 4 <= rl <= 8

      sj, cj = sincos(d.tht[Hj])

      B = d.B[Hj]
      R1 = d.R1[Hj]
      R2 = d.R2[Hj]

      θ0, dθ = if rl == 4
         T5, T6 - T5
      elseif rl == 5
         T6, T7 - T6
      elseif rl == 6
         T7, T8 - T7
      elseif rl == 7
         T8, T9 - T8
      else
         T9, T5 - T9 + 2.0π
      end

      θ = muladd(ξ, dθ, θ0)
      st, ct = sincos(θ)

      return B + R1 * ct * sj + R2 * st * cj

   else

      throw(ArgumentError("Yy local region must be 1:8; got rl=$rl, r=$r"))

   end

end

function fill_FTable!(tbl::FTable, d::ellipseNh, r::Int)

   vmin, vmax, N = tbl.vmin, tbl.vmax, tbl.N
   h = (vmax - vmin) / (N - 1)

   P1, P2, P3 = tbl.P1, tbl.P2, tbl.P3

   @inbounds begin

      H = cld(r, 8)
      rl = r - 8 * (H - 1)

      Hj = H + 1

      T1  = d.T[1, H]
      T2  = d.T[2, H]
      T3  = d.T[3, H]
      T4  = d.T[4, H]
      T5  = d.T[5, H]
      T6  = d.T[6, H]
      T7  = d.T[7, H]
      T8  = d.T[8, H]
      T9  = d.T[9, H]
      T10 = d.T[10, H]
      T11 = d.T[11, H]
      T12 = d.T[12, H]
      T13 = d.T[13, H]

      # Outer ellipse constants.
      s1, c1 = sincos(d.tht[1])

      A1 = d.A[1]
      B1 = d.B[1]

      R11 = d.R1[1]
      R12 = d.R2[1]

      # Hole ellipse constants.
      sj, cj = sincos(d.tht[Hj])

      Aj = d.A[Hj]
      Bj = d.B[Hj]

      Rj1 = d.R1[Hj]
      Rj2 = d.R2[Hj]

      if rl == 1

         dθ = T2 - T1

         st5, ct5 = sincos(T5)
         st2, ct2 = sincos(T2)

         gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
         gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

         g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
         g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

         XRx = 0.5 * (gjx5 + g1x2)
         XRy = 0.5 * (gjy5 + g1y2)

         Xpx = XRx - T12
         Xpy = XRy - T13

         for i in 1:N

            vhat = vmin + (i - 1) * h

            Xxv = muladd(vhat, Xpx, T12)
            Xyv = muladd(vhat, Xpy, T13)

            θ = muladd(vhat, dθ, T1)
            st, ct = sincos(θ)

            Yxv = A1 + R11 * ct * c1 - R12 * st * s1
            Yyv = B1 + R11 * ct * s1 + R12 * st * c1

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv

         end

      elseif rl == 2

         dθo = T3 - T2
         dθh = T6 - T5

         for i in 1:N

            vhat = vmin + (i - 1) * h

            θo = muladd(vhat, dθo, T2)
            sto, cto = sincos(θo)

            gx = A1 + R11 * cto * c1 - R12 * sto * s1
            gy = B1 + R11 * cto * s1 + R12 * sto * c1

            θh = muladd(vhat, dθh, T5)
            sth, cth = sincos(θh)

            gjx = Aj + Rj1 * cth * cj - Rj2 * sth * sj
            gjy = Bj + Rj1 * cth * sj + Rj2 * sth * cj

            Xxv = 0.5 * (gx + gjx)
            Xyv = 0.5 * (gy + gjy)

            Yxv = gx
            Yyv = gy

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv

         end

      elseif rl == 3

         dθ = T4 - T3

         st6, ct6 = sincos(T6)
         st3, ct3 = sincos(T3)

         gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
         gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

         g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
         g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

         XLx = 0.5 * (gjx6 + g1x3)
         XLy = 0.5 * (gjy6 + g1y3)

         Xpx = T10 - XLx
         Xpy = T11 - XLy

         for i in 1:N

            vhat = vmin + (i - 1) * h

            Xxv = muladd(vhat, Xpx, XLx)
            Xyv = muladd(vhat, Xpy, XLy)

            θ = muladd(vhat, dθ, T3)
            st, ct = sincos(θ)

            Yxv = A1 + R11 * ct * c1 - R12 * st * s1
            Yyv = B1 + R11 * ct * s1 + R12 * st * c1

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv

         end

      elseif rl == 4

         dθo = T3 - T2
         dθh = T6 - T5

         for i in 1:N

            vhat = vmin + (i - 1) * h

            θo = muladd(vhat, dθo, T2)
            sto, cto = sincos(θo)

            gx = A1 + R11 * cto * c1 - R12 * sto * s1
            gy = B1 + R11 * cto * s1 + R12 * sto * c1

            θh = muladd(vhat, dθh, T5)
            sth, cth = sincos(θh)

            gjx = Aj + Rj1 * cth * cj - Rj2 * sth * sj
            gjy = Bj + Rj1 * cth * sj + Rj2 * sth * cj

            Xxv = 0.5 * (gx + gjx)
            Xyv = 0.5 * (gy + gjy)

            Yxv = gjx
            Yyv = gjy

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv

         end

      elseif rl == 5

         dθ = T7 - T6

         st6, ct6 = sincos(T6)
         st3, ct3 = sincos(T3)

         gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
         gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

         g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
         g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

         XLx = 0.5 * (gjx6 + g1x3)
         XLy = 0.5 * (gjy6 + g1y3)

         Xpx = T10 - XLx
         Xpy = T11 - XLy

         for i in 1:N

            vhat = vmin + (i - 1) * h

            Xxv = muladd(vhat, Xpx, XLx)
            Xyv = muladd(vhat, Xpy, XLy)

            θ = muladd(vhat, dθ, T6)
            st, ct = sincos(θ)

            Yxv = Aj + Rj1 * ct * cj - Rj2 * st * sj
            Yyv = Bj + Rj1 * ct * sj + Rj2 * st * cj

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv

         end

      elseif rl == 6

         dθ = T8 - T7

         Xpx = A1 - T10
         Xpy = B1 - T11

         for i in 1:N

            vhat = vmin + (i - 1) * h

            Xxv = muladd(vhat, Xpx, T10)
            Xyv = muladd(vhat, Xpy, T11)

            θ = muladd(vhat, dθ, T7)
            st, ct = sincos(θ)

            Yxv = Aj + Rj1 * ct * cj - Rj2 * st * sj
            Yyv = Bj + Rj1 * ct * sj + Rj2 * st * cj

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv

         end

      elseif rl == 7

         dθ = T9 - T8

         Xpx = T12 - A1
         Xpy = T13 - B1

         for i in 1:N

            vhat = vmin + (i - 1) * h

            Xxv = muladd(vhat, Xpx, A1)
            Xyv = muladd(vhat, Xpy, B1)

            θ = muladd(vhat, dθ, T8)
            st, ct = sincos(θ)

            Yxv = Aj + Rj1 * ct * cj - Rj2 * st * sj
            Yyv = Bj + Rj1 * ct * sj + Rj2 * st * cj

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv

         end

      elseif rl == 8

         dθ = T5 - T9 + 2.0π

         st5, ct5 = sincos(T5)
         st2, ct2 = sincos(T2)

         gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
         gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

         g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
         g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

         XRx = 0.5 * (gjx5 + g1x2)
         XRy = 0.5 * (gjy5 + g1y2)

         Xpx = XRx - T12
         Xpy = XRy - T13

         for i in 1:N

            vhat = vmin + (i - 1) * h

            Xxv = muladd(vhat, Xpx, T12)
            Xyv = muladd(vhat, Xpy, T13)

            θ = muladd(vhat, dθ, T9)
            st, ct = sincos(θ)

            Yxv = Aj + Rj1 * ct * cj - Rj2 * st * sj
            Yyv = Bj + Rj1 * ct * sj + Rj2 * st * cj

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv

         end

      else

         throw(ArgumentError("fill_FTable!(ellipseNh): local region must be 1:8; got rl=$rl, r=$r"))

      end

   end

   tbl.reg = r

   return tbl

end

@inline function f1I(t::Float64, s::Float64,
   d::ellipseNh, u::Float64, v::Float64, r::Int)

   t̂ = muladd(0.5, t, 0.5)

   Xxv = Xx(s, d, r)
   Yxv = Yx(s, d, r)

   return muladd(t̂, Yxv - Xxv, Xxv) - u

end

@inline function f2I(t::Float64, s::Float64,
   d::ellipseNh, u::Float64, v::Float64, r::Int)

   t̂ = muladd(0.5, t, 0.5)

   Xyv = Xy(s, d, r)
   Yyv = Yy(s, d, r)

   return muladd(t̂, Yyv - Xyv, Xyv) - v

end

@inline function JinvI(t::Float64, s::Float64,
   d::ellipseNh, u::Float64, v::Float64, r::Int)

   return jinvmap(d, t, s, r)

end

@inline function f_cont(v̂::Float64,
   d::ellipseNh, u::Float64, v::Float64, r::Int)

   @inbounds begin

      H = cld(r, 8)
      rl = r - 8 * (H - 1)

      Hj = H + 1

      T1  = d.T[1, H]
      T2  = d.T[2, H]
      T3  = d.T[3, H]
      T4  = d.T[4, H]
      T5  = d.T[5, H]
      T6  = d.T[6, H]
      T7  = d.T[7, H]
      T8  = d.T[8, H]
      T9  = d.T[9, H]
      T10 = d.T[10, H]
      T11 = d.T[11, H]
      T12 = d.T[12, H]
      T13 = d.T[13, H]

      s1, c1 = sincos(d.tht[1])
      sj, cj = sincos(d.tht[Hj])

      A1 = d.A[1]
      B1 = d.B[1]
      R11 = d.R1[1]
      R12 = d.R2[1]

      Aj = d.A[Hj]
      Bj = d.B[Hj]
      Rj1 = d.R1[Hj]
      Rj2 = d.R2[Hj]

      if rl == 1

         dθ = T2 - T1

         st5, ct5 = sincos(T5)
         st2, ct2 = sincos(T2)

         gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
         gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

         g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
         g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

         XRx = 0.5 * (gjx5 + g1x2)
         XRy = 0.5 * (gjy5 + g1y2)

         Xxv = muladd(v̂, XRx - T12, T12)
         Xyv = muladd(v̂, XRy - T13, T13)

         θ = muladd(v̂, dθ, T1)
         st, ct = sincos(θ)

         Yxv = A1 + R11 * ct * c1 - R12 * st * s1
         Yyv = B1 + R11 * ct * s1 + R12 * st * c1

      elseif rl == 2

         dθo = T3 - T2
         dθh = T6 - T5

         θo = muladd(v̂, dθo, T2)
         sto, cto = sincos(θo)

         gx = A1 + R11 * cto * c1 - R12 * sto * s1
         gy = B1 + R11 * cto * s1 + R12 * sto * c1

         θh = muladd(v̂, dθh, T5)
         sth, cth = sincos(θh)

         gjx = Aj + Rj1 * cth * cj - Rj2 * sth * sj
         gjy = Bj + Rj1 * cth * sj + Rj2 * sth * cj

         Xxv = 0.5 * (gx + gjx)
         Xyv = 0.5 * (gy + gjy)

         Yxv = gx
         Yyv = gy

      elseif rl == 3

         dθ = T4 - T3

         st6, ct6 = sincos(T6)
         st3, ct3 = sincos(T3)

         gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
         gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

         g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
         g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

         XLx = 0.5 * (gjx6 + g1x3)
         XLy = 0.5 * (gjy6 + g1y3)

         Xxv = muladd(v̂, T10 - XLx, XLx)
         Xyv = muladd(v̂, T11 - XLy, XLy)

         θ = muladd(v̂, dθ, T3)
         st, ct = sincos(θ)

         Yxv = A1 + R11 * ct * c1 - R12 * st * s1
         Yyv = B1 + R11 * ct * s1 + R12 * st * c1

      elseif rl == 4

         dθo = T3 - T2
         dθh = T6 - T5

         θo = muladd(v̂, dθo, T2)
         sto, cto = sincos(θo)

         gx = A1 + R11 * cto * c1 - R12 * sto * s1
         gy = B1 + R11 * cto * s1 + R12 * sto * c1

         θh = muladd(v̂, dθh, T5)
         sth, cth = sincos(θh)

         gjx = Aj + Rj1 * cth * cj - Rj2 * sth * sj
         gjy = Bj + Rj1 * cth * sj + Rj2 * sth * cj

         Xxv = 0.5 * (gx + gjx)
         Xyv = 0.5 * (gy + gjy)

         Yxv = gjx
         Yyv = gjy

      elseif rl == 5

         dθ = T7 - T6

         st6, ct6 = sincos(T6)
         st3, ct3 = sincos(T3)

         gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
         gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

         g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
         g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

         XLx = 0.5 * (gjx6 + g1x3)
         XLy = 0.5 * (gjy6 + g1y3)

         Xxv = muladd(v̂, T10 - XLx, XLx)
         Xyv = muladd(v̂, T11 - XLy, XLy)

         θ = muladd(v̂, dθ, T6)
         st, ct = sincos(θ)

         Yxv = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yyv = Bj + Rj1 * ct * sj + Rj2 * st * cj

      elseif rl == 6

         dθ = T8 - T7

         Xxv = muladd(v̂, A1 - T10, T10)
         Xyv = muladd(v̂, B1 - T11, T11)

         θ = muladd(v̂, dθ, T7)
         st, ct = sincos(θ)

         Yxv = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yyv = Bj + Rj1 * ct * sj + Rj2 * st * cj

      elseif rl == 7

         dθ = T9 - T8

         Xxv = muladd(v̂, T12 - A1, A1)
         Xyv = muladd(v̂, T13 - B1, B1)

         θ = muladd(v̂, dθ, T8)
         st, ct = sincos(θ)

         Yxv = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yyv = Bj + Rj1 * ct * sj + Rj2 * st * cj

      elseif rl == 8

         dθ = T5 - T9 + 2.0π

         st5, ct5 = sincos(T5)
         st2, ct2 = sincos(T2)

         gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
         gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

         g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
         g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

         XRx = 0.5 * (gjx5 + g1x2)
         XRy = 0.5 * (gjy5 + g1y2)

         Xxv = muladd(v̂, XRx - T12, T12)
         Xyv = muladd(v̂, XRy - T13, T13)

         θ = muladd(v̂, dθ, T9)
         st, ct = sincos(θ)

         Yxv = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yyv = Bj + Rj1 * ct * sj + Rj2 * st * cj

      else

         throw(ArgumentError("f_cont(ellipseNh): local region must be 1:8; got rl=$rl, r=$r"))

      end

   end

   return u * (Yyv - Xyv) - v * (Yxv - Xxv) +
          (Xyv * Yxv - Xxv * Yyv)

end

function mapinv(tbl::FTable, d::ellipseNh, u::Float64,
   v::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]

   r = p.reg

   # Stage 1: local 2D Newton.
   tN, sN = newtonR2D(f1I, f2I, JinvI,
      0.0, 0.0, 4, d, u, v, r; tol=1e-15)

   if tN !== :max

      x = xi_inv((1.0 + tN) / 2.0, p.ck0, p.ck1)
      y = xi_inv((1.0 + sN) / 2.0, p.tk0, p.tk1)

      return x, y

   end

   # Stage 2: BIS/table fallback.
   if tbl.reg != r
      display("Wrong table!")
      fill_FTable!(tbl, d, r)
   end

   rmi = tbl.rmi
   zxi = tbl.zxi

   n = find_roots!(rmi, tbl, d, u, v, r)

   @inbounds for i in 1:n

      # rmi[i] is in table coordinates [0,1].
      sref = 2.0 * rmi[i] - 1.0

      Xxv = Xx(sref, d, r)
      Xyv = Xy(sref, d, r)
      Yxv = Yx(sref, d, r)
      Yyv = Yy(sref, d, r)

      if abs(v - Xyv) < abs(u - Xxv)
         zxi[i] = (u - Xxv) / (Yxv - Xxv)
      else
         zxi[i] = (v - Xyv) / (Yyv - Xyv)
      end

   end

   best_cost = Inf
   best_idx = 0

   @inbounds for i in 1:n

      zy_norm = xi_inv(rmi[i], p.tk0, p.tk1)
      zx_norm = xi_inv(zxi[i], p.ck0, p.ck1)

      cost = max(abs(zy_norm), abs(zx_norm))

      if cost < best_cost
         best_cost = cost
         best_idx = i
      end

   end

   za = rmi[best_idx]
   zx = zxi[best_idx]

   t = xi_inv(zx, p.ck0, p.ck1)
   s = xi_inv(za, p.tk0, p.tk1)

   return t, s

end

"""
    ptconv(d::ellipseNh, t1::Float64, t2::Float64, idx::Int, ptdest::String)

Convert a parametric point between region ↔ patch coordinates.

Inputs:
- `t1`: first coordinate of the point.
- `t2`: second coordinate of the point.
- `idx`: region index if `ptdest == "to_pth"`, or patch index if `ptdest == "to_reg"`.
- `ptdest`: `"to_pth"` or `"to_reg"`.

Returns:
- `(t1out, t2out, out_idx)`, where `out_idx` is the patch index for `"to_pth"`
  and the region index for `"to_reg"`.
"""
function ptconv(d::ellipseNh, t1::Float64, t2::Float64,
   idx::Int, ptdest::String)

   if ptdest == "to_pth"

      r = idx

      @inbounds for k in 1:d.Npat

         p = d.pths[k]
         p.reg == r || continue

         t1k = xi_inv((t1 + 1.0) / 2.0, p.ck0, p.ck1)
         t2k = xi_inv((t2 + 1.0) / 2.0, p.tk0, p.tk1)

         if abs(t1k) <= 1.0 && abs(t2k) <= 1.0
            return t1k, t2k, k
         end

      end

      error("ptconv(to_pth): no patch in region $r contained the point.")

   elseif ptdest == "to_reg"

      k = idx
      p = d.pths[k]
      r = p.reg

      t1r = (p.ck1 - p.ck0) * t1 + (p.ck1 + p.ck0) - 1.0
      t2r = (p.tk1 - p.tk0) * t2 + (p.tk1 + p.tk0) - 1.0

      return t1r, t2r, r

   else

      error("ptconv: ptdest must be \"to_pth\" or \"to_reg\"")

   end

end

"""
    mapinv2(d::ellipseNh, t1::Float64, t2::Float64, k2::Int, k::Int)

Given a point `(u,v) = τₖ₂(t1,t2)` on patch `k2`, return its reference
coordinates on patch `k` in `[-1,1]^2`.

This differs from `mapinv` because we already know the point comes from
`(t1,t2)` on patch `k2`.
"""
function mapinv2(d::ellipseNh, t1::Float64, t2::Float64,
   k2::Int, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]

   tr1, tr2, _ = ptconv(d, t1, t2, k2, "to_reg")

   û = (tr1 + 1.0) / 2.0
   v̂ = (tr2 + 1.0) / 2.0

   t1k = xi_inv(û, p.ck0, p.ck1)
   t2k = xi_inv(v̂, p.tk0, p.tk1)

   return t1k, t2k

end

"""
    dfunc(d::ellipseNh, k::Int, t::Float64, s::Float64)

Boundary-distance factor raised to the method-dependent exponent.

Notes:
- `t` is expected to be `1 - t_actual`.
- The raw distance factor is

      1 - ck1 + (ck1 - ck0) * t / 2

- The exponent is `s - 1` if `s >= 0.5`, otherwise `s`.
"""
function dfunc(d::ellipseNh, k::Int, t::Float64, s::Float64)::Float64

   p = d.pths[k]

   hc = p.ck1 - p.ck0
   val = (1.0 - p.ck1) + hc * t / 2.0

   exp = s >= 0.5 ? (s - 1.0) : s

   return val^exp

end

function dfunc!(out::StridedArray{Float64}, d::ellipseNh, k::Int,
   t::StridedArray{Float64}, s::Float64)

   p = d.pths[k]

   exp = s >= 0.5 ? (s - 1.0) : s

   αc = 0.5 * (p.ck1 - p.ck0)
   βc = 1.0 - p.ck1

   @inbounds for i in eachindex(out, t)
      out[i] = muladd(αc, t[i], βc)^exp
   end

   return nothing

end

"""
    dfunc!(out, d::ellipseNh, k, t)

Raw boundary-distance factor, without raising to any power.
"""
function dfunc!(out::StridedArray{Float64}, d::ellipseNh, k::Int,
   t::StridedArray{Float64})

   p = d.pths[k]

   αc = 0.5 * (p.ck1 - p.ck0)
   βc = 1.0 - p.ck1

   @inbounds for i in eachindex(out, t)
      out[i] = muladd(αc, t[i], βc)
   end

   return nothing

end


"""
  diff_map!(out::Matrix{Float64}, Zx::Matrix{Float64},
   Zy::Matrix{Float64}, DJ::StridedArray{Float64},
   d::ellipseNh, u::Float64, v::Float64,
   u2::Matrix{Float64}, v2::Matrix{Float64},
   du::AbstractVector, dv::AbstractVector, k::Int;
   tol = 1e-3)

Fill `out` with ‖τ(u,v) - τ(u₂,v₂)‖ on patch `k`.

No allocations. Uses `mapxy_Dmap!` and a sixth-order Taylor correction
near `(u,v)`.
"""
function diff_map!(out::Matrix{Float64},
   Zx::Matrix{Float64}, Zy::Matrix{Float64}, DJ::StridedArray{Float64},
   d::ellipseNh, u::Float64, v::Float64,
   u2::Matrix{Float64}, v2::Matrix{Float64},
   du::AbstractVector, dv::AbstractVector, k::Int;
   tol::Float64=1e-3)

   nd_u = size(out, 1)
   nd_v = size(out, 2)

   p = d.pths[k]

   αc = 0.5 * (p.ck1 - p.ck0)
   βc = 0.5 * (p.ck0 + p.ck1)

   αt = 0.5 * (p.tk1 - p.tk0)
   βt = 0.5 * (p.tk0 + p.tk1)

   mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)

   xi1 = muladd(αc, u, βc)
   xi2 = muladd(αt, v, βt)

   reg = p.reg

   H = cld(reg, 8)
   r = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]
   T10 = d.T[10, H]
   T11 = d.T[11, H]
   T12 = d.T[12, H]
   T13 = d.T[13, H]

   s1, c1 = sincos(d.tht[1])
   sj, cj = sincos(d.tht[Hj])

   A1 = d.A[1]
   B1 = d.B[1]
   R11 = d.R1[1]
   R12 = d.R2[1]

   Aj = d.A[Hj]
   Bj = d.B[Hj]
   Rj1 = d.R1[Hj]
   Rj2 = d.R2[Hj]

   Xx = Xy = 0.0
   X1x = X1y = 0.0
   X2x = X2y = 0.0
   X3x = X3y = 0.0
   X4x = X4y = 0.0
   X5x = X5y = 0.0
   X6x = X6y = 0.0

   Yx = Yy = 0.0
   Y1x = Y1y = 0.0
   Y2x = Y2y = 0.0
   Y3x = Y3y = 0.0
   Y4x = Y4y = 0.0
   Y5x = Y5y = 0.0
   Y6x = Y6y = 0.0

   if r == 1

      dθ = T2 - T1
      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      XRx = 0.5 * (gjx5 + g1x2)
      XRy = 0.5 * (gjy5 + g1y2)

      X1x = XRx - T12
      X1y = XRy - T13

      Xx = muladd(xi2, X1x, T12)
      Xy = muladd(xi2, X1y, T13)

      θ = muladd(xi2, dθ, T1)
      st, ct = sincos(θ)

      Yx = A1 + R11 * ct * c1 - R12 * st * s1
      Yy = B1 + R11 * ct * s1 + R12 * st * c1

      g1x = -R11 * st * c1 - R12 * ct * s1
      g1y = -R11 * st * s1 + R12 * ct * c1

      g2x = A1 - Yx
      g2y = B1 - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   elseif r == 2 || r == 4

      dθo = T3 - T2
      dθh = T6 - T5

      dθo2 = dθo * dθo
      dθo3 = dθo2 * dθo
      dθo4 = dθo2 * dθo2
      dθo5 = dθo4 * dθo
      dθo6 = dθo3 * dθo3

      dθh2 = dθh * dθh
      dθh3 = dθh2 * dθh
      dθh4 = dθh2 * dθh2
      dθh5 = dθh4 * dθh
      dθh6 = dθh3 * dθh3

      θo = muladd(xi2, dθo, T2)
      sto, cto = sincos(θo)

      gox = A1 + R11 * cto * c1 - R12 * sto * s1
      goy = B1 + R11 * cto * s1 + R12 * sto * c1

      go1x = -R11 * sto * c1 - R12 * cto * s1
      go1y = -R11 * sto * s1 + R12 * cto * c1

      go2x = A1 - gox
      go2y = B1 - goy

      θh = muladd(xi2, dθh, T5)
      sth, cth = sincos(θh)

      ghx = Aj + Rj1 * cth * cj - Rj2 * sth * sj
      ghy = Bj + Rj1 * cth * sj + Rj2 * sth * cj

      gh1x = -Rj1 * sth * cj - Rj2 * cth * sj
      gh1y = -Rj1 * sth * sj + Rj2 * cth * cj

      gh2x = Aj - ghx
      gh2y = Bj - ghy

      Xx = 0.5 * (gox + ghx)
      Xy = 0.5 * (goy + ghy)

      X1x = 0.5 * (dθo * go1x + dθh * gh1x)
      X1y = 0.5 * (dθo * go1y + dθh * gh1y)

      X2x = 0.5 * (dθo2 * go2x + dθh2 * gh2x)
      X2y = 0.5 * (dθo2 * go2y + dθh2 * gh2y)

      X3x = -0.5 * (dθo3 * go1x + dθh3 * gh1x)
      X3y = -0.5 * (dθo3 * go1y + dθh3 * gh1y)

      X4x = -0.5 * (dθo4 * go2x + dθh4 * gh2x)
      X4y = -0.5 * (dθo4 * go2y + dθh4 * gh2y)

      X5x = 0.5 * (dθo5 * go1x + dθh5 * gh1x)
      X5y = 0.5 * (dθo5 * go1y + dθh5 * gh1y)

      X6x = 0.5 * (dθo6 * go2x + dθh6 * gh2x)
      X6y = 0.5 * (dθo6 * go2y + dθh6 * gh2y)

      if r == 2

         Yx = gox
         Yy = goy

         Y1x = dθo * go1x
         Y1y = dθo * go1y
         Y2x = dθo2 * go2x
         Y2y = dθo2 * go2y
         Y3x = -dθo3 * go1x
         Y3y = -dθo3 * go1y
         Y4x = -dθo4 * go2x
         Y4y = -dθo4 * go2y

         Y5x = dθo5 * go1x
         Y5y = dθo5 * go1y
         Y6x = dθo6 * go2x
         Y6y = dθo6 * go2y

      else

         Yx = ghx
         Yy = ghy

         Y1x = dθh * gh1x
         Y1y = dθh * gh1y
         Y2x = dθh2 * gh2x
         Y2y = dθh2 * gh2y
         Y3x = -dθh3 * gh1x
         Y3y = -dθh3 * gh1y
         Y4x = -dθh4 * gh2x
         Y4y = -dθh4 * gh2y

         Y5x = dθh5 * gh1x
         Y5y = dθh5 * gh1y
         Y6x = dθh6 * gh2x
         Y6y = dθh6 * gh2y

      end

   elseif r == 3

      dθ = T4 - T3
      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      XLx = 0.5 * (gjx6 + g1x3)
      XLy = 0.5 * (gjy6 + g1y3)

      X1x = T10 - XLx
      X1y = T11 - XLy

      Xx = muladd(xi2, X1x, XLx)
      Xy = muladd(xi2, X1y, XLy)

      θ = muladd(xi2, dθ, T3)
      st, ct = sincos(θ)

      Yx = A1 + R11 * ct * c1 - R12 * st * s1
      Yy = B1 + R11 * ct * s1 + R12 * st * c1

      g1x = -R11 * st * c1 - R12 * ct * s1
      g1y = -R11 * st * s1 + R12 * ct * c1

      g2x = A1 - Yx
      g2y = B1 - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   elseif r == 5

      dθ = T7 - T6
      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      XLx = 0.5 * (gjx6 + g1x3)
      XLy = 0.5 * (gjy6 + g1y3)

      X1x = T10 - XLx
      X1y = T11 - XLy

      Xx = muladd(xi2, X1x, XLx)
      Xy = muladd(xi2, X1y, XLy)

      θ = muladd(xi2, dθ, T6)
      st, ct = sincos(θ)

      Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
      Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

      g1x = -Rj1 * st * cj - Rj2 * ct * sj
      g1y = -Rj1 * st * sj + Rj2 * ct * cj

      g2x = Aj - Yx
      g2y = Bj - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   elseif r == 6

      dθ = T8 - T7
      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      X1x = A1 - T10
      X1y = B1 - T11

      Xx = muladd(xi2, X1x, T10)
      Xy = muladd(xi2, X1y, T11)

      θ = muladd(xi2, dθ, T7)
      st, ct = sincos(θ)

      Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
      Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

      g1x = -Rj1 * st * cj - Rj2 * ct * sj
      g1y = -Rj1 * st * sj + Rj2 * ct * cj

      g2x = Aj - Yx
      g2y = Bj - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   elseif r == 7

      dθ = T9 - T8
      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      X1x = T12 - A1
      X1y = T13 - B1

      Xx = muladd(xi2, X1x, A1)
      Xy = muladd(xi2, X1y, B1)

      θ = muladd(xi2, dθ, T8)
      st, ct = sincos(θ)

      Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
      Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

      g1x = -Rj1 * st * cj - Rj2 * ct * sj
      g1y = -Rj1 * st * sj + Rj2 * ct * cj

      g2x = Aj - Yx
      g2y = Bj - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   elseif r == 8

      dθ = T5 - T9 + 2.0π
      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      XRx = 0.5 * (gjx5 + g1x2)
      XRy = 0.5 * (gjy5 + g1y2)

      X1x = XRx - T12
      X1y = XRy - T13

      Xx = muladd(xi2, X1x, T12)
      Xy = muladd(xi2, X1y, T13)

      θ = muladd(xi2, dθ, T9)
      st, ct = sincos(θ)

      Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
      Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

      g1x = -Rj1 * st * cj - Rj2 * ct * sj
      g1y = -Rj1 * st * sj + Rj2 * ct * cj

      g2x = Aj - Yx
      g2y = Bj - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   else

      throw(ArgumentError("diff_map! local region must be 1:8; got r=$r, reg=$reg"))

   end

   αt2 = αt * αt
   αt3 = αt2 * αt
   αt4 = αt2 * αt2
   αt5 = αt4 * αt
   αt6 = αt3 * αt3

   dux = αc * (Yx - Xx)
   duy = αc * (Yy - Xy)

   duvx = αc * αt * (Y1x - X1x)
   duvy = αc * αt * (Y1y - X1y)

   duv2x = αc * αt2 * (Y2x - X2x)
   duv2y = αc * αt2 * (Y2y - X2y)

   duv3x = αc * αt3 * (Y3x - X3x)
   duv3y = αc * αt3 * (Y3y - X3y)

   duv4x = αc * αt4 * (Y4x - X4x)
   duv4y = αc * αt4 * (Y4y - X4y)

   duv5x = αc * αt5 * (Y5x - X5x)
   duv5y = αc * αt5 * (Y5y - X5y)

   dvx = αt * ((1.0 - xi1) * X1x + xi1 * Y1x)
   dvy = αt * ((1.0 - xi1) * X1y + xi1 * Y1y)

   dv2x = αt2 * ((1.0 - xi1) * X2x + xi1 * Y2x)
   dv2y = αt2 * ((1.0 - xi1) * X2y + xi1 * Y2y)

   dv3x = αt3 * ((1.0 - xi1) * X3x + xi1 * Y3x)
   dv3y = αt3 * ((1.0 - xi1) * X3y + xi1 * Y3y)

   dv4x = αt4 * ((1.0 - xi1) * X4x + xi1 * Y4x)
   dv4y = αt4 * ((1.0 - xi1) * X4y + xi1 * Y4y)

   dv5x = αt5 * ((1.0 - xi1) * X5x + xi1 * Y5x)
   dv5y = αt5 * ((1.0 - xi1) * X5y + xi1 * Y5y)

   dv6x = αt6 * ((1.0 - xi1) * X6x + xi1 * Y6x)
   dv6y = αt6 * ((1.0 - xi1) * X6y + xi1 * Y6y)

   tux, tvy = mapxy(d, u, v, k)

   @inbounds for j in 1:nd_v

      dvj = dv[j]
      dvj2 = dvj * dvj
      dvj3 = dvj2 * dvj
      dvj4 = dvj2 * dvj2
      dvj5 = dvj4 * dvj
      dvj6 = dvj3 * dvj3

      @inbounds for i in 1:nd_u

         uu = u2[i, j]
         vv = v2[i, j]

         if abs(u - uu) < tol && abs(v - vv) < tol

            dui = du[i]

            Dx = (dui * dux + dvj * dvx) -
                 (dui * dvj * duvx + dvj2 * dv2x / 2.0) +
                 (dui * dvj2 * duv2x / 2.0 + dvj3 * dv3x / 6.0) -
                 (dui * dvj3 * duv3x / 6.0 + dvj4 * dv4x / 24.0) +
                 (dui * dvj4 * duv4x / 24.0 + dvj5 * dv5x / 120.0) -
                 (dui * dvj5 * duv5x / 120.0 + dvj6 * dv6x / 720.0)

            Dy = (dui * duy + dvj * dvy) -
                 (dui * dvj * duvy + dvj2 * dv2y / 2.0) +
                 (dui * dvj2 * duv2y / 2.0 + dvj3 * dv3y / 6.0) -
                 (dui * dvj3 * duv3y / 6.0 + dvj4 * dv4y / 24.0) +
                 (dui * dvj4 * duv4y / 24.0 + dvj5 * dv5y / 120.0) -
                 (dui * dvj5 * duv5y / 120.0 + dvj6 * dv6y / 720.0)

            out[i, j] = hypot(Dx, Dy)

         else

            out[i, j] = hypot(tux - Zx[i, j], tvy - Zy[i, j])

         end

      end

   end

   return nothing

end

"""
    diff_map!(out::Vector{Float64},
      Zx::Vector{Float64}, Zy::Vector{Float64}, DJℓ::Vector{Float64},
      d::ellipseNh, u::Float64, v::Float64,
      u2::Float64, v2::Vector{Float64},
      du::Float64, dv::Vector{Float64}, k::Int;
      tol::Float64 = 1e-3)

Vector version of `diff_map!`.

Computes out[j] = ‖τ_k(u,v) - τ_k(u2,v2[j])‖

with Taylor fixup near `(u,v)`. No allocations.

Also fills Zx[j], Zy[j] = τ_k(u2,v2[j]) and DJℓ[j]= ‖∂τ/∂v(u2,v2[j])‖
"""
function diff_map!(out::Vector{Float64},
   Zx::Vector{Float64}, Zy::Vector{Float64}, DJℓ::Vector{Float64},
   d::ellipseNh, u::Float64, v::Float64,
   u2::Float64, v2::Vector{Float64},
   du::Float64, dv::Vector{Float64}, k::Int;
   tol::Float64=1e-3)

   nd_v = length(out)

   p = d.pths[k]

   αc = 0.5 * (p.ck1 - p.ck0)
   βc = 0.5 * (p.ck0 + p.ck1)

   αt = 0.5 * (p.tk1 - p.tk0)
   βt = 0.5 * (p.tk0 + p.tk1)

   xi1 = muladd(αc, u, βc)
   xi2 = muladd(αt, v, βt)

   xi1q = muladd(αc, u2, βc)

   reg = p.reg

   H = cld(reg, 8)
   rl = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]
   T10 = d.T[10, H]
   T11 = d.T[11, H]
   T12 = d.T[12, H]
   T13 = d.T[13, H]

   s1, c1 = sincos(d.tht[1])
   sj, cj = sincos(d.tht[Hj])

   A1 = d.A[1]
   B1 = d.B[1]
   R11 = d.R1[1]
   R12 = d.R2[1]

   Aj = d.A[Hj]
   Bj = d.B[Hj]
   Rj1 = d.R1[Hj]
   Rj2 = d.R2[Hj]

   # ------------------------------------------------------------
   # τ(u,v), the fixed evaluation point.
   # ------------------------------------------------------------
   tux, tvy = mapxy(d, u, v, k)

   # ------------------------------------------------------------
   # Fill Zx, Zy, DJℓ along the line (u2, v2[j]).
   # ------------------------------------------------------------
   if rl == 1

      dθ = T2 - T1

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      XRx = 0.5 * (gjx5 + g1x2)
      XRy = 0.5 * (gjy5 + g1y2)

      X1x = XRx - T12
      X1y = XRy - T13

      @inbounds for j in 1:nd_v

         xi2q = muladd(αt, v2[j], βt)

         Xxq = muladd(xi2q, X1x, T12)
         Xyq = muladd(xi2q, X1y, T13)

         θ = muladd(xi2q, dθ, T1)
         st, ct = sincos(θ)

         Yxq = A1 + R11 * ct * c1 - R12 * st * s1
         Yyq = B1 + R11 * ct * s1 + R12 * st * c1

         Y1xq = dθ * (-R11 * st * c1 - R12 * ct * s1)
         Y1yq = dθ * (-R11 * st * s1 + R12 * ct * c1)

         Zx[j] = muladd(xi1q, Yxq - Xxq, Xxq)
         Zy[j] = muladd(xi1q, Yyq - Xyq, Xyq)

         dvxq = αt * ((1.0 - xi1q) * X1x + xi1q * Y1xq)
         dvyq = αt * ((1.0 - xi1q) * X1y + xi1q * Y1yq)

         DJℓ[j] = hypot(dvxq, dvyq)

      end

   elseif rl == 2 || rl == 4

      dθo = T3 - T2
      dθh = T6 - T5

      @inbounds for j in 1:nd_v

         xi2q = muladd(αt, v2[j], βt)

         θo = muladd(xi2q, dθo, T2)
         sto, cto = sincos(θo)

         gox = A1 + R11 * cto * c1 - R12 * sto * s1
         goy = B1 + R11 * cto * s1 + R12 * sto * c1

         go1x = dθo * (-R11 * sto * c1 - R12 * cto * s1)
         go1y = dθo * (-R11 * sto * s1 + R12 * cto * c1)

         θh = muladd(xi2q, dθh, T5)
         sth, cth = sincos(θh)

         ghx = Aj + Rj1 * cth * cj - Rj2 * sth * sj
         ghy = Bj + Rj1 * cth * sj + Rj2 * sth * cj

         gh1x = dθh * (-Rj1 * sth * cj - Rj2 * cth * sj)
         gh1y = dθh * (-Rj1 * sth * sj + Rj2 * cth * cj)

         Xxq = 0.5 * (gox + ghx)
         Xyq = 0.5 * (goy + ghy)

         X1xq = 0.5 * (go1x + gh1x)
         X1yq = 0.5 * (go1y + gh1y)

         if rl == 2

            Yxq = gox
            Yyq = goy
            Y1xq = go1x
            Y1yq = go1y

         else

            Yxq = ghx
            Yyq = ghy
            Y1xq = gh1x
            Y1yq = gh1y

         end

         Zx[j] = muladd(xi1q, Yxq - Xxq, Xxq)
         Zy[j] = muladd(xi1q, Yyq - Xyq, Xyq)

         dvxq = αt * ((1.0 - xi1q) * X1xq + xi1q * Y1xq)
         dvyq = αt * ((1.0 - xi1q) * X1yq + xi1q * Y1yq)

         DJℓ[j] = hypot(dvxq, dvyq)

      end

   elseif rl == 3

      dθ = T4 - T3

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      XLx = 0.5 * (gjx6 + g1x3)
      XLy = 0.5 * (gjy6 + g1y3)

      X1x = T10 - XLx
      X1y = T11 - XLy

      @inbounds for j in 1:nd_v

         xi2q = muladd(αt, v2[j], βt)

         Xxq = muladd(xi2q, X1x, XLx)
         Xyq = muladd(xi2q, X1y, XLy)

         θ = muladd(xi2q, dθ, T3)
         st, ct = sincos(θ)

         Yxq = A1 + R11 * ct * c1 - R12 * st * s1
         Yyq = B1 + R11 * ct * s1 + R12 * st * c1

         Y1xq = dθ * (-R11 * st * c1 - R12 * ct * s1)
         Y1yq = dθ * (-R11 * st * s1 + R12 * ct * c1)

         Zx[j] = muladd(xi1q, Yxq - Xxq, Xxq)
         Zy[j] = muladd(xi1q, Yyq - Xyq, Xyq)

         dvxq = αt * ((1.0 - xi1q) * X1x + xi1q * Y1xq)
         dvyq = αt * ((1.0 - xi1q) * X1y + xi1q * Y1yq)

         DJℓ[j] = hypot(dvxq, dvyq)

      end

   elseif rl == 5

      dθ = T7 - T6

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      XLx = 0.5 * (gjx6 + g1x3)
      XLy = 0.5 * (gjy6 + g1y3)

      X1x = T10 - XLx
      X1y = T11 - XLy

      @inbounds for j in 1:nd_v

         xi2q = muladd(αt, v2[j], βt)

         Xxq = muladd(xi2q, X1x, XLx)
         Xyq = muladd(xi2q, X1y, XLy)

         θ = muladd(xi2q, dθ, T6)
         st, ct = sincos(θ)

         Yxq = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yyq = Bj + Rj1 * ct * sj + Rj2 * st * cj

         Y1xq = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         Y1yq = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         Zx[j] = muladd(xi1q, Yxq - Xxq, Xxq)
         Zy[j] = muladd(xi1q, Yyq - Xyq, Xyq)

         dvxq = αt * ((1.0 - xi1q) * X1x + xi1q * Y1xq)
         dvyq = αt * ((1.0 - xi1q) * X1y + xi1q * Y1yq)

         DJℓ[j] = hypot(dvxq, dvyq)

      end

   elseif rl == 6

      dθ = T8 - T7

      X1x = A1 - T10
      X1y = B1 - T11

      @inbounds for j in 1:nd_v

         xi2q = muladd(αt, v2[j], βt)

         Xxq = muladd(xi2q, X1x, T10)
         Xyq = muladd(xi2q, X1y, T11)

         θ = muladd(xi2q, dθ, T7)
         st, ct = sincos(θ)

         Yxq = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yyq = Bj + Rj1 * ct * sj + Rj2 * st * cj

         Y1xq = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         Y1yq = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         Zx[j] = muladd(xi1q, Yxq - Xxq, Xxq)
         Zy[j] = muladd(xi1q, Yyq - Xyq, Xyq)

         dvxq = αt * ((1.0 - xi1q) * X1x + xi1q * Y1xq)
         dvyq = αt * ((1.0 - xi1q) * X1y + xi1q * Y1yq)

         DJℓ[j] = hypot(dvxq, dvyq)

      end

   elseif rl == 7

      dθ = T9 - T8

      X1x = T12 - A1
      X1y = T13 - B1

      @inbounds for j in 1:nd_v

         xi2q = muladd(αt, v2[j], βt)

         Xxq = muladd(xi2q, X1x, A1)
         Xyq = muladd(xi2q, X1y, B1)

         θ = muladd(xi2q, dθ, T8)
         st, ct = sincos(θ)

         Yxq = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yyq = Bj + Rj1 * ct * sj + Rj2 * st * cj

         Y1xq = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         Y1yq = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         Zx[j] = muladd(xi1q, Yxq - Xxq, Xxq)
         Zy[j] = muladd(xi1q, Yyq - Xyq, Xyq)

         dvxq = αt * ((1.0 - xi1q) * X1x + xi1q * Y1xq)
         dvyq = αt * ((1.0 - xi1q) * X1y + xi1q * Y1yq)

         DJℓ[j] = hypot(dvxq, dvyq)

      end

   elseif rl == 8

      dθ = T5 - T9 + 2.0π

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      XRx = 0.5 * (gjx5 + g1x2)
      XRy = 0.5 * (gjy5 + g1y2)

      X1x = XRx - T12
      X1y = XRy - T13

      @inbounds for j in 1:nd_v

         xi2q = muladd(αt, v2[j], βt)

         Xxq = muladd(xi2q, X1x, T12)
         Xyq = muladd(xi2q, X1y, T13)

         θ = muladd(xi2q, dθ, T9)
         st, ct = sincos(θ)

         Yxq = Aj + Rj1 * ct * cj - Rj2 * st * sj
         Yyq = Bj + Rj1 * ct * sj + Rj2 * st * cj

         Y1xq = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         Y1yq = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         Zx[j] = muladd(xi1q, Yxq - Xxq, Xxq)
         Zy[j] = muladd(xi1q, Yyq - Xyq, Xyq)

         dvxq = αt * ((1.0 - xi1q) * X1x + xi1q * Y1xq)
         dvyq = αt * ((1.0 - xi1q) * X1y + xi1q * Y1yq)

         DJℓ[j] = hypot(dvxq, dvyq)

      end

   else

      throw(ArgumentError("diff_map!(Vector, ellipseNh): local region must be 1:8; got rl=$rl, reg=$reg"))

   end

   # ------------------------------------------------------------
   # Taylor data at the fixed point (u,v), through sixth order.
   # ------------------------------------------------------------
   Xx = Xy = 0.0
   X1x = X1y = 0.0
   X2x = X2y = 0.0
   X3x = X3y = 0.0
   X4x = X4y = 0.0
   X5x = X5y = 0.0
   X6x = X6y = 0.0

   Yx = Yy = 0.0
   Y1x = Y1y = 0.0
   Y2x = Y2y = 0.0
   Y3x = Y3y = 0.0
   Y4x = Y4y = 0.0
   Y5x = Y5y = 0.0
   Y6x = Y6y = 0.0

   if rl == 1

      dθ = T2 - T1
      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      XRx = 0.5 * (gjx5 + g1x2)
      XRy = 0.5 * (gjy5 + g1y2)

      X1x = XRx - T12
      X1y = XRy - T13

      Xx = muladd(xi2, X1x, T12)
      Xy = muladd(xi2, X1y, T13)

      θ = muladd(xi2, dθ, T1)
      st, ct = sincos(θ)

      Yx = A1 + R11 * ct * c1 - R12 * st * s1
      Yy = B1 + R11 * ct * s1 + R12 * st * c1

      g1x = -R11 * st * c1 - R12 * ct * s1
      g1y = -R11 * st * s1 + R12 * ct * c1

      g2x = A1 - Yx
      g2y = B1 - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   elseif rl == 2 || rl == 4

      dθo = T3 - T2
      dθh = T6 - T5

      dθo2 = dθo * dθo
      dθo3 = dθo2 * dθo
      dθo4 = dθo2 * dθo2
      dθo5 = dθo4 * dθo
      dθo6 = dθo3 * dθo3

      dθh2 = dθh * dθh
      dθh3 = dθh2 * dθh
      dθh4 = dθh2 * dθh2
      dθh5 = dθh4 * dθh
      dθh6 = dθh3 * dθh3

      θo = muladd(xi2, dθo, T2)
      sto, cto = sincos(θo)

      gox = A1 + R11 * cto * c1 - R12 * sto * s1
      goy = B1 + R11 * cto * s1 + R12 * sto * c1

      go1x = -R11 * sto * c1 - R12 * cto * s1
      go1y = -R11 * sto * s1 + R12 * cto * c1

      go2x = A1 - gox
      go2y = B1 - goy

      θh = muladd(xi2, dθh, T5)
      sth, cth = sincos(θh)

      ghx = Aj + Rj1 * cth * cj - Rj2 * sth * sj
      ghy = Bj + Rj1 * cth * sj + Rj2 * sth * cj

      gh1x = -Rj1 * sth * cj - Rj2 * cth * sj
      gh1y = -Rj1 * sth * sj + Rj2 * cth * cj

      gh2x = Aj - ghx
      gh2y = Bj - ghy

      Xx = 0.5 * (gox + ghx)
      Xy = 0.5 * (goy + ghy)

      X1x = 0.5 * (dθo * go1x + dθh * gh1x)
      X1y = 0.5 * (dθo * go1y + dθh * gh1y)

      X2x = 0.5 * (dθo2 * go2x + dθh2 * gh2x)
      X2y = 0.5 * (dθo2 * go2y + dθh2 * gh2y)

      X3x = -0.5 * (dθo3 * go1x + dθh3 * gh1x)
      X3y = -0.5 * (dθo3 * go1y + dθh3 * gh1y)

      X4x = -0.5 * (dθo4 * go2x + dθh4 * gh2x)
      X4y = -0.5 * (dθo4 * go2y + dθh4 * gh2y)

      X5x = 0.5 * (dθo5 * go1x + dθh5 * gh1x)
      X5y = 0.5 * (dθo5 * go1y + dθh5 * gh1y)

      X6x = 0.5 * (dθo6 * go2x + dθh6 * gh2x)
      X6y = 0.5 * (dθo6 * go2y + dθh6 * gh2y)

      if rl == 2

         Yx = gox
         Yy = goy

         Y1x = dθo * go1x
         Y1y = dθo * go1y
         Y2x = dθo2 * go2x
         Y2y = dθo2 * go2y
         Y3x = -dθo3 * go1x
         Y3y = -dθo3 * go1y
         Y4x = -dθo4 * go2x
         Y4y = -dθo4 * go2y

         Y5x = dθo5 * go1x
         Y5y = dθo5 * go1y
         Y6x = dθo6 * go2x
         Y6y = dθo6 * go2y

      else

         Yx = ghx
         Yy = ghy

         Y1x = dθh * gh1x
         Y1y = dθh * gh1y
         Y2x = dθh2 * gh2x
         Y2y = dθh2 * gh2y
         Y3x = -dθh3 * gh1x
         Y3y = -dθh3 * gh1y
         Y4x = -dθh4 * gh2x
         Y4y = -dθh4 * gh2y

         Y5x = dθh5 * gh1x
         Y5y = dθh5 * gh1y
         Y6x = dθh6 * gh2x
         Y6y = dθh6 * gh2y

      end

   elseif rl == 3

      dθ = T4 - T3
      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      XLx = 0.5 * (gjx6 + g1x3)
      XLy = 0.5 * (gjy6 + g1y3)

      X1x = T10 - XLx
      X1y = T11 - XLy

      Xx = muladd(xi2, X1x, XLx)
      Xy = muladd(xi2, X1y, XLy)

      θ = muladd(xi2, dθ, T3)
      st, ct = sincos(θ)

      Yx = A1 + R11 * ct * c1 - R12 * st * s1
      Yy = B1 + R11 * ct * s1 + R12 * st * c1

      g1x = -R11 * st * c1 - R12 * ct * s1
      g1y = -R11 * st * s1 + R12 * ct * c1

      g2x = A1 - Yx
      g2y = B1 - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   elseif rl == 5

      dθ = T7 - T6
      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      XLx = 0.5 * (gjx6 + g1x3)
      XLy = 0.5 * (gjy6 + g1y3)

      X1x = T10 - XLx
      X1y = T11 - XLy

      Xx = muladd(xi2, X1x, XLx)
      Xy = muladd(xi2, X1y, XLy)

      θ = muladd(xi2, dθ, T6)
      st, ct = sincos(θ)

      Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
      Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

      g1x = -Rj1 * st * cj - Rj2 * ct * sj
      g1y = -Rj1 * st * sj + Rj2 * ct * cj

      g2x = Aj - Yx
      g2y = Bj - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   elseif rl == 6

      dθ = T8 - T7
      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      X1x = A1 - T10
      X1y = B1 - T11

      Xx = muladd(xi2, X1x, T10)
      Xy = muladd(xi2, X1y, T11)

      θ = muladd(xi2, dθ, T7)
      st, ct = sincos(θ)

      Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
      Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

      g1x = -Rj1 * st * cj - Rj2 * ct * sj
      g1y = -Rj1 * st * sj + Rj2 * ct * cj

      g2x = Aj - Yx
      g2y = Bj - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   elseif rl == 7

      dθ = T9 - T8
      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      X1x = T12 - A1
      X1y = T13 - B1

      Xx = muladd(xi2, X1x, A1)
      Xy = muladd(xi2, X1y, B1)

      θ = muladd(xi2, dθ, T8)
      st, ct = sincos(θ)

      Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
      Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

      g1x = -Rj1 * st * cj - Rj2 * ct * sj
      g1y = -Rj1 * st * sj + Rj2 * ct * cj

      g2x = Aj - Yx
      g2y = Bj - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   elseif rl == 8

      dθ = T5 - T9 + 2.0π
      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      XRx = 0.5 * (gjx5 + g1x2)
      XRy = 0.5 * (gjy5 + g1y2)

      X1x = XRx - T12
      X1y = XRy - T13

      Xx = muladd(xi2, X1x, T12)
      Xy = muladd(xi2, X1y, T13)

      θ = muladd(xi2, dθ, T9)
      st, ct = sincos(θ)

      Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
      Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

      g1x = -Rj1 * st * cj - Rj2 * ct * sj
      g1y = -Rj1 * st * sj + Rj2 * ct * cj

      g2x = Aj - Yx
      g2y = Bj - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   else

      throw(ArgumentError("diff_map! local region must be 1:8; got rl=$rl, reg=$reg"))

   end

   αt2 = αt * αt
   αt3 = αt2 * αt
   αt4 = αt2 * αt2
   αt5 = αt4 * αt
   αt6 = αt3 * αt3

   dux = αc * (Yx - Xx)
   duy = αc * (Yy - Xy)

   duvx = αc * αt * (Y1x - X1x)
   duvy = αc * αt * (Y1y - X1y)

   duv2x = αc * αt2 * (Y2x - X2x)
   duv2y = αc * αt2 * (Y2y - X2y)

   duv3x = αc * αt3 * (Y3x - X3x)
   duv3y = αc * αt3 * (Y3y - X3y)

   duv4x = αc * αt4 * (Y4x - X4x)
   duv4y = αc * αt4 * (Y4y - X4y)

   duv5x = αc * αt5 * (Y5x - X5x)
   duv5y = αc * αt5 * (Y5y - X5y)

   dvx = αt * ((1.0 - xi1) * X1x + xi1 * Y1x)
   dvy = αt * ((1.0 - xi1) * X1y + xi1 * Y1y)

   dv2x = αt2 * ((1.0 - xi1) * X2x + xi1 * Y2x)
   dv2y = αt2 * ((1.0 - xi1) * X2y + xi1 * Y2y)

   dv3x = αt3 * ((1.0 - xi1) * X3x + xi1 * Y3x)
   dv3y = αt3 * ((1.0 - xi1) * X3y + xi1 * Y3y)

   dv4x = αt4 * ((1.0 - xi1) * X4x + xi1 * Y4x)
   dv4y = αt4 * ((1.0 - xi1) * X4y + xi1 * Y4y)

   dv5x = αt5 * ((1.0 - xi1) * X5x + xi1 * Y5x)
   dv5y = αt5 * ((1.0 - xi1) * X5y + xi1 * Y5y)

   dv6x = αt6 * ((1.0 - xi1) * X6x + xi1 * Y6x)
   dv6y = αt6 * ((1.0 - xi1) * X6y + xi1 * Y6y)

   @inbounds for j in 1:nd_v

      vv = v2[j]

      if abs(u - u2) < tol && abs(v - vv) < tol

         dvj = dv[j]
         dvj2 = dvj * dvj
         dvj3 = dvj2 * dvj
         dvj4 = dvj2 * dvj2
         dvj5 = dvj4 * dvj
         dvj6 = dvj3 * dvj3

         Dx = (du * dux + dvj * dvx) -
              (du * dvj * duvx + dvj2 * dv2x / 2.0) +
              (du * dvj2 * duv2x / 2.0 + dvj3 * dv3x / 6.0) -
              (du * dvj3 * duv3x / 6.0 + dvj4 * dv4x / 24.0) +
              (du * dvj4 * duv4x / 24.0 + dvj5 * dv5x / 120.0) -
              (du * dvj5 * duv5x / 120.0 + dvj6 * dv6x / 720.0)

         Dy = (du * duy + dvj * dvy) -
              (du * dvj * duvy + dvj2 * dv2y / 2.0) +
              (du * dvj2 * duv2y / 2.0 + dvj3 * dv3y / 6.0) -
              (du * dvj3 * duv3y / 6.0 + dvj4 * dv4y / 24.0) +
              (du * dvj4 * duv4y / 24.0 + dvj5 * dv5y / 120.0) -
              (du * dvj5 * duv5y / 120.0 + dvj6 * dv6y / 720.0)

         out[j] = hypot(Dx, Dy)

      else

         out[j] = hypot(tux - Zx[j], tvy - Zy[j])

      end

   end


   return nothing

end

"""
    diff_rmap!(out::Matrix{Float64}, Zx::Matrix{Float64}, Zy::Matrix{Float64},
               DJ::StridedArray{Float64}, d::ellipseNh,
               u::Float64, v::Float64,
               u2::Matrix{Float64}, v2::Matrix{Float64}, r::Matrix{Float64},
               du::AbstractVector, dv::AbstractVector, k::Int;
               tol::Float64 = 1e-3)

Compute ‖(τ(u,v) - τ(u₂,v₂)) / r‖ for the `k`-th patch, where

    u₂ = u - r .* du
    v₂ = v - r .* dv

No allocations. Uses `mapxy_Dmap!` and a sixth-order Taylor correction near `(u,v)`.
"""
function diff_rmap!(out::Matrix{Float64},
   Zx::Matrix{Float64}, Zy::Matrix{Float64}, DJ::StridedArray{Float64},
   d::ellipseNh, u::Float64, v::Float64,
   u2::Matrix{Float64}, v2::Matrix{Float64}, r::Matrix{Float64},
   du::AbstractVector, dv::AbstractVector, k::Int;
   tol::Float64=1e-3)

   nt = size(out, 1)
   nr = size(out, 2)

   p = d.pths[k]

   αc = 0.5 * (p.ck1 - p.ck0)
   βc = 0.5 * (p.ck0 + p.ck1)

   αt = 0.5 * (p.tk1 - p.tk0)
   βt = 0.5 * (p.tk0 + p.tk1)

   mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)

   xi1 = muladd(αc, u, βc)
   xi2 = muladd(αt, v, βt)

   reg = p.reg

   H = cld(reg, 8)
   rl = reg - 8 * (H - 1)

   Hj = H + 1

   T1 = d.T[1, H]
   T2 = d.T[2, H]
   T3 = d.T[3, H]
   T4 = d.T[4, H]
   T5 = d.T[5, H]
   T6 = d.T[6, H]
   T7 = d.T[7, H]
   T8 = d.T[8, H]
   T9 = d.T[9, H]
   T10 = d.T[10, H]
   T11 = d.T[11, H]
   T12 = d.T[12, H]
   T13 = d.T[13, H]

   s1, c1 = sincos(d.tht[1])
   sj, cj = sincos(d.tht[Hj])

   A1 = d.A[1]
   B1 = d.B[1]
   R11 = d.R1[1]
   R12 = d.R2[1]

   Aj = d.A[Hj]
   Bj = d.B[Hj]
   Rj1 = d.R1[Hj]
   Rj2 = d.R2[Hj]

   Xx = Xy = 0.0
   X1x = X1y = 0.0
   X2x = X2y = 0.0
   X3x = X3y = 0.0
   X4x = X4y = 0.0
   X5x = X5y = 0.0
   X6x = X6y = 0.0

   Yx = Yy = 0.0
   Y1x = Y1y = 0.0
   Y2x = Y2y = 0.0
   Y3x = Y3y = 0.0
   Y4x = Y4y = 0.0
   Y5x = Y5y = 0.0
   Y6x = Y6y = 0.0

   if rl == 1

      dθ = T2 - T1

      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      XRx = 0.5 * (gjx5 + g1x2)
      XRy = 0.5 * (gjy5 + g1y2)

      X1x = XRx - T12
      X1y = XRy - T13

      Xx = muladd(xi2, X1x, T12)
      Xy = muladd(xi2, X1y, T13)

      θ = muladd(xi2, dθ, T1)
      st, ct = sincos(θ)

      Yx = A1 + R11 * ct * c1 - R12 * st * s1
      Yy = B1 + R11 * ct * s1 + R12 * st * c1

      g1x = -R11 * st * c1 - R12 * ct * s1
      g1y = -R11 * st * s1 + R12 * ct * c1

      g2x = A1 - Yx
      g2y = B1 - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   elseif rl == 2 || rl == 4

      dθo = T3 - T2
      dθh = T6 - T5

      dθo2 = dθo * dθo
      dθo3 = dθo2 * dθo
      dθo4 = dθo2 * dθo2
      dθo5 = dθo4 * dθo
      dθo6 = dθo3 * dθo3

      dθh2 = dθh * dθh
      dθh3 = dθh2 * dθh
      dθh4 = dθh2 * dθh2
      dθh5 = dθh4 * dθh
      dθh6 = dθh3 * dθh3

      θo = muladd(xi2, dθo, T2)
      sto, cto = sincos(θo)

      gox = A1 + R11 * cto * c1 - R12 * sto * s1
      goy = B1 + R11 * cto * s1 + R12 * sto * c1

      go1x = -R11 * sto * c1 - R12 * cto * s1
      go1y = -R11 * sto * s1 + R12 * cto * c1

      go2x = A1 - gox
      go2y = B1 - goy

      θh = muladd(xi2, dθh, T5)
      sth, cth = sincos(θh)

      ghx = Aj + Rj1 * cth * cj - Rj2 * sth * sj
      ghy = Bj + Rj1 * cth * sj + Rj2 * sth * cj

      gh1x = -Rj1 * sth * cj - Rj2 * cth * sj
      gh1y = -Rj1 * sth * sj + Rj2 * cth * cj

      gh2x = Aj - ghx
      gh2y = Bj - ghy

      Xx = 0.5 * (gox + ghx)
      Xy = 0.5 * (goy + ghy)

      X1x = 0.5 * (dθo * go1x + dθh * gh1x)
      X1y = 0.5 * (dθo * go1y + dθh * gh1y)

      X2x = 0.5 * (dθo2 * go2x + dθh2 * gh2x)
      X2y = 0.5 * (dθo2 * go2y + dθh2 * gh2y)

      X3x = -0.5 * (dθo3 * go1x + dθh3 * gh1x)
      X3y = -0.5 * (dθo3 * go1y + dθh3 * gh1y)

      X4x = -0.5 * (dθo4 * go2x + dθh4 * gh2x)
      X4y = -0.5 * (dθo4 * go2y + dθh4 * gh2y)

      X5x = 0.5 * (dθo5 * go1x + dθh5 * gh1x)
      X5y = 0.5 * (dθo5 * go1y + dθh5 * gh1y)

      X6x = 0.5 * (dθo6 * go2x + dθh6 * gh2x)
      X6y = 0.5 * (dθo6 * go2y + dθh6 * gh2y)

      if rl == 2

         Yx = gox
         Yy = goy

         Y1x = dθo * go1x
         Y1y = dθo * go1y
         Y2x = dθo2 * go2x
         Y2y = dθo2 * go2y
         Y3x = -dθo3 * go1x
         Y3y = -dθo3 * go1y
         Y4x = -dθo4 * go2x
         Y4y = -dθo4 * go2y
         Y5x = dθo5 * go1x
         Y5y = dθo5 * go1y
         Y6x = dθo6 * go2x
         Y6y = dθo6 * go2y

      else

         Yx = ghx
         Yy = ghy

         Y1x = dθh * gh1x
         Y1y = dθh * gh1y
         Y2x = dθh2 * gh2x
         Y2y = dθh2 * gh2y
         Y3x = -dθh3 * gh1x
         Y3y = -dθh3 * gh1y
         Y4x = -dθh4 * gh2x
         Y4y = -dθh4 * gh2y
         Y5x = dθh5 * gh1x
         Y5y = dθh5 * gh1y
         Y6x = dθh6 * gh2x
         Y6y = dθh6 * gh2y

      end

   elseif rl == 3

      dθ = T4 - T3

      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      XLx = 0.5 * (gjx6 + g1x3)
      XLy = 0.5 * (gjy6 + g1y3)

      X1x = T10 - XLx
      X1y = T11 - XLy

      Xx = muladd(xi2, X1x, XLx)
      Xy = muladd(xi2, X1y, XLy)

      θ = muladd(xi2, dθ, T3)
      st, ct = sincos(θ)

      Yx = A1 + R11 * ct * c1 - R12 * st * s1
      Yy = B1 + R11 * ct * s1 + R12 * st * c1

      g1x = -R11 * st * c1 - R12 * ct * s1
      g1y = -R11 * st * s1 + R12 * ct * c1

      g2x = A1 - Yx
      g2y = B1 - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   elseif rl == 5

      dθ = T7 - T6

      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      st6, ct6 = sincos(T6)
      st3, ct3 = sincos(T3)

      gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
      gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

      g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
      g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

      XLx = 0.5 * (gjx6 + g1x3)
      XLy = 0.5 * (gjy6 + g1y3)

      X1x = T10 - XLx
      X1y = T11 - XLy

      Xx = muladd(xi2, X1x, XLx)
      Xy = muladd(xi2, X1y, XLy)

      θ = muladd(xi2, dθ, T6)
      st, ct = sincos(θ)

      Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
      Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

      g1x = -Rj1 * st * cj - Rj2 * ct * sj
      g1y = -Rj1 * st * sj + Rj2 * ct * cj

      g2x = Aj - Yx
      g2y = Bj - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   elseif rl == 6

      dθ = T8 - T7

      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      X1x = A1 - T10
      X1y = B1 - T11

      Xx = muladd(xi2, X1x, T10)
      Xy = muladd(xi2, X1y, T11)

      θ = muladd(xi2, dθ, T7)
      st, ct = sincos(θ)

      Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
      Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

      g1x = -Rj1 * st * cj - Rj2 * ct * sj
      g1y = -Rj1 * st * sj + Rj2 * ct * cj

      g2x = Aj - Yx
      g2y = Bj - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   elseif rl == 7

      dθ = T9 - T8

      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      X1x = T12 - A1
      X1y = T13 - B1

      Xx = muladd(xi2, X1x, A1)
      Xy = muladd(xi2, X1y, B1)

      θ = muladd(xi2, dθ, T8)
      st, ct = sincos(θ)

      Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
      Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

      g1x = -Rj1 * st * cj - Rj2 * ct * sj
      g1y = -Rj1 * st * sj + Rj2 * ct * cj

      g2x = Aj - Yx
      g2y = Bj - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   elseif rl == 8

      dθ = T5 - T9 + 2.0π

      dθ2 = dθ * dθ
      dθ3 = dθ2 * dθ
      dθ4 = dθ2 * dθ2
      dθ5 = dθ4 * dθ
      dθ6 = dθ3 * dθ3

      st5, ct5 = sincos(T5)
      st2, ct2 = sincos(T2)

      gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
      gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

      g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
      g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

      XRx = 0.5 * (gjx5 + g1x2)
      XRy = 0.5 * (gjy5 + g1y2)

      X1x = XRx - T12
      X1y = XRy - T13

      Xx = muladd(xi2, X1x, T12)
      Xy = muladd(xi2, X1y, T13)

      θ = muladd(xi2, dθ, T9)
      st, ct = sincos(θ)

      Yx = Aj + Rj1 * ct * cj - Rj2 * st * sj
      Yy = Bj + Rj1 * ct * sj + Rj2 * st * cj

      g1x = -Rj1 * st * cj - Rj2 * ct * sj
      g1y = -Rj1 * st * sj + Rj2 * ct * cj

      g2x = Aj - Yx
      g2y = Bj - Yy

      Y1x = dθ * g1x
      Y1y = dθ * g1y
      Y2x = dθ2 * g2x
      Y2y = dθ2 * g2y
      Y3x = -dθ3 * g1x
      Y3y = -dθ3 * g1y
      Y4x = -dθ4 * g2x
      Y4y = -dθ4 * g2y
      Y5x = dθ5 * g1x
      Y5y = dθ5 * g1y
      Y6x = dθ6 * g2x
      Y6y = dθ6 * g2y

   else

      throw(ArgumentError("diff_rmap!(ellipseNh): local region must be 1:8; got rl=$rl, reg=$reg"))

   end

   αt2 = αt * αt
   αt3 = αt2 * αt
   αt4 = αt2 * αt2
   αt5 = αt4 * αt
   αt6 = αt3 * αt3

   dux = αc * (Yx - Xx)
   duy = αc * (Yy - Xy)

   duvx = αc * αt * (Y1x - X1x)
   duvy = αc * αt * (Y1y - X1y)

   duv2x = αc * αt2 * (Y2x - X2x)
   duv2y = αc * αt2 * (Y2y - X2y)

   duv3x = αc * αt3 * (Y3x - X3x)
   duv3y = αc * αt3 * (Y3y - X3y)

   duv4x = αc * αt4 * (Y4x - X4x)
   duv4y = αc * αt4 * (Y4y - X4y)

   duv5x = αc * αt5 * (Y5x - X5x)
   duv5y = αc * αt5 * (Y5y - X5y)

   dvx = αt * ((1.0 - xi1) * X1x + xi1 * Y1x)
   dvy = αt * ((1.0 - xi1) * X1y + xi1 * Y1y)

   dv2x = αt2 * ((1.0 - xi1) * X2x + xi1 * Y2x)
   dv2y = αt2 * ((1.0 - xi1) * X2y + xi1 * Y2y)

   dv3x = αt3 * ((1.0 - xi1) * X3x + xi1 * Y3x)
   dv3y = αt3 * ((1.0 - xi1) * X3y + xi1 * Y3y)

   dv4x = αt4 * ((1.0 - xi1) * X4x + xi1 * Y4x)
   dv4y = αt4 * ((1.0 - xi1) * X4y + xi1 * Y4y)

   dv5x = αt5 * ((1.0 - xi1) * X5x + xi1 * Y5x)
   dv5y = αt5 * ((1.0 - xi1) * X5y + xi1 * Y5y)

   dv6x = αt6 * ((1.0 - xi1) * X6x + xi1 * Y6x)
   dv6y = αt6 * ((1.0 - xi1) * X6y + xi1 * Y6y)

   tux, tvy = mapxy(d, u, v, k)

   @inbounds for i in 1:nt

      dui = du[i]
      dvi = dv[i]

      @inbounds for j in 1:nr

         uu = u2[i, j]
         vv = v2[i, j]

         if abs(u - uu) < tol && abs(v - vv) < tol

            rij = r[i, j]

            r1 = dvi * rij
            r2 = r1 * r1
            r3 = r2 * r1
            r4 = r2 * r2
            r5 = r4 * r1

            Dx = (dui * dux + dvi * dvx) -
                 r1 * (dui * duvx + dvi * dv2x / 2.0) +
                 r2 * (dui * duv2x / 2.0 + dvi * dv3x / 6.0) -
                 r3 * (dui * duv3x / 6.0 + dvi * dv4x / 24.0) +
                 r4 * (dui * duv4x / 24.0 + dvi * dv5x / 120.0) -
                 r5 * (dui * duv5x / 120.0 + dvi * dv6x / 720.0)

            Dy = (dui * duy + dvi * dvy) -
                 r1 * (dui * duvy + dvi * dv2y / 2.0) +
                 r2 * (dui * duv2y / 2.0 + dvi * dv3y / 6.0) -
                 r3 * (dui * duv3y / 6.0 + dvi * dv4y / 24.0) +
                 r4 * (dui * duv4y / 24.0 + dvi * dv5y / 120.0) -
                 r5 * (dui * duv5y / 120.0 + dvi * dv6y / 720.0)

            out[i, j] = hypot(Dx, Dy)

         else

            out[i, j] = hypot(tux - Zx[i, j], tvy - Zy[i, j]) / r[i, j]

         end

      end

   end

   return nothing

end

"""
    Dwall(d::ellipseNh, u::Float64, v::Float64, k::Int) -> dvx, dvy

Derivative of the wall curve with respect to `v` at the side `u = ±1`
for the `k`-th patch.

Returns `(dvx, dvy)`.
"""
function Dwall(d::ellipseNh, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64}

   @inbounds begin

      p = d.pths[k]

      αu = 0.5 * (p.ck1 - p.ck0)
      βu = 0.5 * (p.ck0 + p.ck1)

      αv = 0.5 * (p.tk1 - p.tk0)
      βv = 0.5 * (p.tk0 + p.tk1)

      xi1 = muladd(αu, u, βu)
      xi2 = muladd(αv, v, βv)

      reg = p.reg

      H = cld(reg, 8)
      rl = reg - 8 * (H - 1)

      Hj = H + 1

      T1  = d.T[1, H]
      T2  = d.T[2, H]
      T3  = d.T[3, H]
      T4  = d.T[4, H]
      T5  = d.T[5, H]
      T6  = d.T[6, H]
      T7  = d.T[7, H]
      T8  = d.T[8, H]
      T9  = d.T[9, H]
      T10 = d.T[10, H]
      T11 = d.T[11, H]
      T12 = d.T[12, H]
      T13 = d.T[13, H]

      # Outer ellipse constants.
      s1, c1 = sincos(d.tht[1])

      A1 = d.A[1]
      B1 = d.B[1]

      R11 = d.R1[1]
      R12 = d.R2[1]

      # Hole ellipse constants.
      sj, cj = sincos(d.tht[Hj])

      Aj = d.A[Hj]
      Bj = d.B[Hj]

      Rj1 = d.R1[Hj]
      Rj2 = d.R2[Hj]

      if rl == 1

         dθ = T2 - T1

         st5, ct5 = sincos(T5)
         st2, ct2 = sincos(T2)

         gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
         gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

         g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
         g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

         XRx = 0.5 * (gjx5 + g1x2)
         XRy = 0.5 * (gjy5 + g1y2)

         Xpx = XRx - T12
         Xpy = XRy - T13

         θ = muladd(xi2, dθ, T1)
         st, ct = sincos(θ)

         dYx = dθ * (-R11 * st * c1 - R12 * ct * s1)
         dYy = dθ * (-R11 * st * s1 + R12 * ct * c1)

         dvx = (1.0 - xi1) * Xpx + xi1 * dYx
         dvy = (1.0 - xi1) * Xpy + xi1 * dYy

      elseif rl == 2

         dθo = T3 - T2
         dθh = T6 - T5

         θo = muladd(xi2, dθo, T2)
         sto, cto = sincos(θo)

         dgx = -R11 * sto * c1 - R12 * cto * s1
         dgy = -R11 * sto * s1 + R12 * cto * c1

         θh = muladd(xi2, dθh, T5)
         sth, cth = sincos(θh)

         dgjx = -Rj1 * sth * cj - Rj2 * cth * sj
         dgjy = -Rj1 * sth * sj + Rj2 * cth * cj

         dvx = 0.5 * ((1.0 - xi1) * dθh * dgjx + (1.0 + xi1) * dθo * dgx)
         dvy = 0.5 * ((1.0 - xi1) * dθh * dgjy + (1.0 + xi1) * dθo * dgy)

      elseif rl == 3

         dθ = T4 - T3

         st6, ct6 = sincos(T6)
         st3, ct3 = sincos(T3)

         gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
         gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

         g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
         g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

         XLx = 0.5 * (gjx6 + g1x3)
         XLy = 0.5 * (gjy6 + g1y3)

         Xpx = T10 - XLx
         Xpy = T11 - XLy

         θ = muladd(xi2, dθ, T3)
         st, ct = sincos(θ)

         dYx = dθ * (-R11 * st * c1 - R12 * ct * s1)
         dYy = dθ * (-R11 * st * s1 + R12 * ct * c1)

         dvx = (1.0 - xi1) * Xpx + xi1 * dYx
         dvy = (1.0 - xi1) * Xpy + xi1 * dYy

      elseif rl == 4

         dθo = T3 - T2
         dθh = T6 - T5

         θo = muladd(xi2, dθo, T2)
         sto, cto = sincos(θo)

         dgx = -R11 * sto * c1 - R12 * cto * s1
         dgy = -R11 * sto * s1 + R12 * cto * c1

         θh = muladd(xi2, dθh, T5)
         sth, cth = sincos(θh)

         dgjx = -Rj1 * sth * cj - Rj2 * cth * sj
         dgjy = -Rj1 * sth * sj + Rj2 * cth * cj

         dvx = 0.5 * ((1.0 - xi1) * dθo * dgx + (1.0 + xi1) * dθh * dgjx)
         dvy = 0.5 * ((1.0 - xi1) * dθo * dgy + (1.0 + xi1) * dθh * dgjy)

      elseif rl == 5

         dθ = T7 - T6

         st6, ct6 = sincos(T6)
         st3, ct3 = sincos(T3)

         gjx6 = Aj + Rj1 * ct6 * cj - Rj2 * st6 * sj
         gjy6 = Bj + Rj1 * ct6 * sj + Rj2 * st6 * cj

         g1x3 = A1 + R11 * ct3 * c1 - R12 * st3 * s1
         g1y3 = B1 + R11 * ct3 * s1 + R12 * st3 * c1

         XLx = 0.5 * (gjx6 + g1x3)
         XLy = 0.5 * (gjy6 + g1y3)

         Xpx = T10 - XLx
         Xpy = T11 - XLy

         θ = muladd(xi2, dθ, T6)
         st, ct = sincos(θ)

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         dvx = (1.0 - xi1) * Xpx + xi1 * dYx
         dvy = (1.0 - xi1) * Xpy + xi1 * dYy

      elseif rl == 6

         dθ = T8 - T7

         Xpx = A1 - T10
         Xpy = B1 - T11

         θ = muladd(xi2, dθ, T7)
         st, ct = sincos(θ)

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         dvx = (1.0 - xi1) * Xpx + xi1 * dYx
         dvy = (1.0 - xi1) * Xpy + xi1 * dYy

      elseif rl == 7

         dθ = T9 - T8

         Xpx = T12 - A1
         Xpy = T13 - B1

         θ = muladd(xi2, dθ, T8)
         st, ct = sincos(θ)

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         dvx = (1.0 - xi1) * Xpx + xi1 * dYx
         dvy = (1.0 - xi1) * Xpy + xi1 * dYy

      elseif rl == 8

         dθ = T5 - T9 + 2.0π

         st5, ct5 = sincos(T5)
         st2, ct2 = sincos(T2)

         gjx5 = Aj + Rj1 * ct5 * cj - Rj2 * st5 * sj
         gjy5 = Bj + Rj1 * ct5 * sj + Rj2 * st5 * cj

         g1x2 = A1 + R11 * ct2 * c1 - R12 * st2 * s1
         g1y2 = B1 + R11 * ct2 * s1 + R12 * st2 * c1

         XRx = 0.5 * (gjx5 + g1x2)
         XRy = 0.5 * (gjy5 + g1y2)

         Xpx = XRx - T12
         Xpy = XRy - T13

         θ = muladd(xi2, dθ, T9)
         st, ct = sincos(θ)

         dYx = dθ * (-Rj1 * st * cj - Rj2 * ct * sj)
         dYy = dθ * (-Rj1 * st * sj + Rj2 * ct * cj)

         dvx = (1.0 - xi1) * Xpx + xi1 * dYx
         dvy = (1.0 - xi1) * Xpy + xi1 * dYy

      else

         throw(ArgumentError("Dwall(ellipseNh): local region must be 1:8; got rl=$rl, reg=$reg"))

      end

   end

   return αv * dvx, αv * dvy

end

"""
    refine!(d::ellipseNh, Nc::Int, Nt::Int, K::Vector{Int})

Subdivide each patch in `K` into `Nc * Nt` subpatches, update `d.pths`,
`d.Npat`, recompute `d.kd`, and rebuild `d.Qpts` / `d.Qptsbd`.

Patches are re-sorted lexicographically by `(reg, ck0, ck1, tk0, tk1)`.
"""
function refine!(d::ellipseNh, Nc::Int, Nt::Int, K::Vector{Int})

   @assert Nc >= 1 && Nt >= 1 "Nc and Nt must be >= 1"
   @assert !isempty(K) "K, the set of patches to refine, is empty"

   children = Patch[]
   sizehint!(children, length(K) * Nc * Nt)

   @inbounds for k in K

      p = d.pths[k]

      cc = range(p.ck0, p.ck1; length=Nc + 1)
      tt = range(p.tk0, p.tk1; length=Nt + 1)

      for i in 1:Nc, j in 1:Nt
         push!(children, Patch(p.reg, cc[i], cc[i+1], tt[j], tt[j+1]))
      end

   end

   d.pths = vcat([d.pths[i] for i in 1:d.Npat if !(i in K)], children)

   sort!(d.pths, by = q -> (q.reg, q.ck0, q.ck1, q.tk0, q.tk1))

   d.Npat = length(d.pths)

   d.kd = [k for k in 1:d.Npat if d.pths[k].ck1 == 1.0]

   d.Qpts = Matrix{Float64}(undef, 8, d.Npat)

   @inbounds for k in 1:d.Npat
      @views V = d.Qpts[:, k]
      boundquad!(V, d, k)
   end

   d.Qptsbd = Matrix{Float64}(undef, 8, length(d.kd))

   @inbounds for (ℓ, k) in enumerate(d.kd)
      @views V = d.Qptsbd[:, ℓ]
      boundquadbd!(V, d, k)
   end

   return d

end

function Base.show(io::IO, d::ellipseNh)

   println(io, "ellipseNh with properties:")

   println(io, "  nh     = ", d.nh)

   println(io, "  A      = ", d.A)
   println(io, "  B      = ", d.B)

   println(io, "  R1     = ", d.R1)
   println(io, "  R2     = ", d.R2)

   println(io, "  tht    = ", d.tht)

   println(io, "  T      = ", size(d.T, 1), "×", size(d.T, 2), " Matrix")

   if length(d.kd) <= 8
      L = "Int[" * join(d.kd, ' ') * "]"
   else
      head = join(d.kd[1:4], ' ')
      tail = join(d.kd[end-4+1:end], ' ')
      L = "Int[$head … $tail]"
   end

   println(io, "  kd     = ", L)
   println(io, "  Npat   = ", d.Npat)
   println(io, "  pths   = Vector{Patch} (", length(d.pths), " patches)")
   println(io, "  Qpts   = ", size(d.Qpts, 1), "×", size(d.Qpts, 2), " Matrix")
   println(io, "  Qptsbd = ", size(d.Qptsbd, 1), "×", size(d.Qptsbd, 2), " Matrix")

end