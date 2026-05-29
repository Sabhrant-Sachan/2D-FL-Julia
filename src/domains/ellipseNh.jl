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

      # @inbounds for k in 1:Npat
      #    @views boundquad!(d.Qpts[:, k], d, k)
      # end

      # @inbounds for (ℓ, k) in enumerate(d.kd)
      #    @views boundquadbd!(d.Qptsbd[:, ℓ], d, k)
      # end

      return d
   end

end

#TODO