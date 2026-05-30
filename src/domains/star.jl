struct regionparams_star
   e₁x::Vector{Float64}
   e₁y::Vector{Float64}
   e₂x::Vector{Float64}
   e₂y::Vector{Float64}
   qx::Vector{Float64}
   qy::Vector{Float64}
   αx::Vector{Float64}
   αy::Vector{Float64}
   nme::Float64
end

mutable struct star <: abstractdomain

   A::Float64
   B::Float64
   R::Float64
   P::Int
   Ps::Float64
   Rc::Float64
   Lc::Float64
   L1::Float64
   L2::Float64
   tht1::Float64
   tht2::Float64
   nh::Int

   kd::Vector{Int}
   Npat::Int
   pths::Vector{Patch}
   Qpts::Matrix{Float64}
   Qptsbd::Matrix{Float64}

   RP::regionparams_star

   function star(; b,
      a=nothing, A=nothing, B=nothing,
      R=nothing, P=nothing, Ps=nothing,
      Rc=nothing, Lc=nothing,
      L1=nothing, L2=nothing,
      tht1=nothing, tht2=nothing,
      ck=nothing, tk=nothing)

      @assert isa(b, AbstractVector{Int}) && length(b) == 12
      @assert all(b .>= 1)

      nh = 0

      a = isnothing(a) ? ceil.(Int, 2 .* b ./ 3) : collect(Int, a)

      @assert length(a) == 12
      @assert all(a .>= 1)

      A = something(A, 0.0)
      B = something(B, 0.0)

      R = something(R, 1.0)
      P = something(P, 5)
      Ps = something(Ps, 1.0)

      @assert P >= 3
      @assert R > 0.0
      @assert Ps > 0.0

      tht1 = something(tht1, 0.22)
      tht2 = something(tht2, 0.19)

      L1 = something(L1, 0.35 * R)
      L2 = something(L2, 0.45 * R)

      Rc = something(Rc, 0.6 * R)
      Lc = something(Lc, 0.8 * Rc)

      @assert L1 > 0.0 && L2 > 0.0
      @assert Rc > 0.0 && Lc > 0.0

      # Total number of patches:
      #   regions 1:7 repeat P times,
      #   regions 8:12 are the central disc-like regions.
      Npat = P * dot(a[1:7], b[1:7]) + dot(a[8:12], b[8:12])

      # Template partitions, length 12.
      # For actual petal regions, we reuse template regions 1:7.
      # For center regions, we reuse template regions 8:12.
      ck = something(ck, [(0:a[k]) ./ a[k] for k in 1:12])
      tk = something(tk, [(0:b[k]) ./ b[k] for k in 1:12])

      @assert length(ck) == 12
      @assert length(tk) == 12

      # ----------------------------------------------------
      # Region parameters for all petals.
      #   tm = (2j - 1)π/P - θ₁
      #   tp = (2j - 3)π/P + θ₁
      #
      #   P2 = h(tm) * (cos(tm), sin(tm))
      #   P1 = h(tp) * (cos(tp), sin(tp))
      #
      #   e₁ = P2 - P1
      #   e₂ = (P2 + P1)/2
      #   q  = (e₁y, -e₁x)
      # ----------------------------------------------------
      e₁x = Vector{Float64}(undef, P)
      e₁y = Vector{Float64}(undef, P)
      e₂x = Vector{Float64}(undef, P)
      e₂y = Vector{Float64}(undef, P)
      qx = Vector{Float64}(undef, P)
      qy = Vector{Float64}(undef, P)

      # α[j] = Rc * (cos((2j-3)π/P), sin((2j-3)π/P)).
      # α[P+1] = α[1] is useful for region 7j-4.
      αx = Vector{Float64}(undef, P + 1)
      αy = Vector{Float64}(undef, P + 1)

      @inline hstar(θ::Float64) = 1.0 + Ps / 2.0 + Ps * cos(P * θ) / 2.0

      nme = 0.0

      @inbounds for j in 1:P

         tm = (2.0 * j - 1.0) * π / P - tht1
         tp = (2.0 * j - 3.0) * π / P + tht1

         hm = hstar(tm)
         hp = hstar(tp)

         P2x = hm * cos(tm)
         P2y = hm * sin(tm)

         P1x = hp * cos(tp)
         P1y = hp * sin(tp)

         e₁x[j] = P2x - P1x
         e₁y[j] = P2y - P1y

         e₂x[j] = 0.5 * (P2x + P1x)
         e₂y[j] = 0.5 * (P2y + P1y)

         qx[j] = e₁y[j]
         qy[j] = -e₁x[j]

         if j == 1
            nme = hypot(e₁x[j], e₁y[j])
         end

         ψ = (2.0 * j - 3.0) * π / P

         αx[j] = Rc * cos(ψ)
         αy[j] = Rc * sin(ψ)

      end

      αx[P+1] = αx[1]
      αy[P+1] = αy[1]

      RP = regionparams_star(e₁x, e₁y, e₂x, e₂y, qx, qy, αx, αy, nme)

      # ------------------------------------------------------------
      # Build patches.
      #
      # p.reg stores actual physical region:
      #   petal regions : 1, ..., 7P
      #   center regions: 7P+1, ..., 7P+5
      #
      # regloc is the local petal region 1:7.
      # rdisc  is the local central-disc region 1:5.
      # ------------------------------------------------------------

      pths = Vector{Patch}(undef, Npat)

      idx = 1

      @inbounds for jpet in 1:P

         for regloc in 1:7

            reg = regloc + 7 * (jpet - 1)

            ckℓ = ck[regloc]
            tkℓ = tk[regloc]

            for ic in 1:a[regloc]

               ck0 = ckℓ[ic]
               ck1 = ckℓ[ic+1]

               for it in 1:b[regloc]

                  tk0 = tkℓ[it]
                  tk1 = tkℓ[it+1]

                  pths[idx] = Patch(reg, ck0, ck1, tk0, tk1)

                  idx += 1

               end
            end
         end
      end

      # Central disc-like regions.
      # Template regions 8:12 become actual regions 7P+1 : 7P+5.

      @inbounds for rdisc in 1:5

         reg = 7 * P + rdisc

         # Central disc templates are stored in slots 8:12.
         regloc = rdisc + 7

         ckℓ = ck[regloc]
         tkℓ = tk[regloc]

         for ic in 1:a[regloc]

            ck0 = ckℓ[ic]
            ck1 = ckℓ[ic+1]

            for it in 1:b[regloc]

               tk0 = tkℓ[it]
               tk1 = tkℓ[it+1]

               pths[idx] = Patch(reg, ck0, ck1, tk0, tk1)

               idx += 1

            end
         end
      end

      # ------------------------------------------------------------
      # Boundary-touching patches.
      #
      # Outer boundary lives only on petal regions with ck1 == 1.
      #
      # Exclude:
      #   regloc == 4 : inner circular arc
      #   regloc == 7 : petal rectangle
      #
      # Central disc regions are not outer-boundary patches.
      # ------------------------------------------------------------

      kd = Vector{Int}(undef, Npat)
      nkd = 0

      @inbounds for k in 1:Npat

         reg = pths[k].reg

         if reg <= 7 * P

            regloc = mod1(reg, 7)

            if regloc != 4 && regloc != 7 && pths[k].ck1 == 1.0

               nkd += 1
               kd[nkd] = k

            end

         end
      end

      resize!(kd, nkd)

      Qpts = Matrix{Float64}(undef, 8, Npat)
      Qptsbd = Matrix{Float64}(undef, 8, nkd)

      d = new(A, B, R, P, Ps, Rc, Lc, L1, L2,
         tht1, tht2, nh, kd, Npat, pths, Qpts, Qptsbd, RP)

      @inbounds for k in 1:Npat
         @views V = d.Qpts[:, k]
         boundquad!(V, d, k)
      end

      @inbounds for (ℓ, k) in enumerate(d.kd)
         @views V = d.Qptsbd[:, ℓ]
         boundquadbd!(V, d, k)
      end

      return d

   end

end

function mapx(d::star, u::Float64, v::Float64, k::Int)::Float64

   x, _ = mapxy(d, u, v, k)

   return x

end

function mapy(d::star, u::Float64, v::Float64, k::Int)::Float64

   _, y = mapxy(d, u, v, k)

   return y

end

function mapxy(d::star, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]
   reg = p.reg

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   xi1 = muladd(0.5 * hc, u, 0.5 * (p.ck0 + p.ck1))
   xi2 = muladd(0.5 * ht, v, 0.5 * (p.tk0 + p.tk1))

   if reg <= 7 * d.P

      jpet = cld(reg, 7)
      regloc = mod1(reg, 7)

      e1x = d.RP.e₁x[jpet]
      e1y = d.RP.e₁y[jpet]
      e2x = d.RP.e₂x[jpet]
      e2y = d.RP.e₂y[jpet]
      qx = d.RP.qx[jpet]
      qy = d.RP.qy[jpet]
      nme = d.RP.nme

      if regloc == 7

         if d.L2 >= d.L1
            xi1 = muladd(0.5 * ht, u, 0.5 * (p.tk0 + p.tk1))
            xi2 = muladd(0.5 * hc, v, 0.5 * (p.ck0 + p.ck1))
         end

         Zx = d.A + d.R * e2x +
              (2.0 * xi2 - 1.0) * d.L1 * e1x / (2.0 * nme) +
              d.L2 * xi1 * qx / nme

         Zy = d.B + d.R * e2y +
              (2.0 * xi2 - 1.0) * d.L1 * e1y / (2.0 * nme) +
              d.L2 * xi1 * qy / nme

         return Zx, Zy

      end

      if regloc == 1

         Xx = d.A + d.R * e2x +
              (2.0 * xi2 - 1.0) * d.L1 * e1x / (2.0 * nme) +
              d.L2 * qx / nme

         Xy = d.B + d.R * e2y +
              (2.0 * xi2 - 1.0) * d.L1 * e1y / (2.0 * nme) +
              d.L2 * qy / nme

         th = xi2 * (2.0 * d.tht2) +
              2.0 * (jpet - 1) * π / d.P - d.tht2

      elseif regloc == 2

         Xx = d.A + d.R * e2x +
              d.L1 * e1x / (2.0 * nme) +
              (1.0 - xi2) * d.L2 * qx / nme

         Xy = d.B + d.R * e2y +
              d.L1 * e1y / (2.0 * nme) +
              (1.0 - xi2) * d.L2 * qy / nme

         th = xi2 * (π / d.P - d.tht2 - d.tht1) +
              2.0 * (jpet - 1) * π / d.P + d.tht2

      elseif regloc == 3

         alpx = d.RP.αx[jpet+1]
         alpy = d.RP.αy[jpet+1]

         Xx = d.A + xi2 * alpx +
              (1.0 - xi2) * (d.R * e2x + d.L1 * e1x / (2.0 * nme))

         Xy = d.B + xi2 * alpy +
              (1.0 - xi2) * (d.R * e2y + d.L1 * e1y / (2.0 * nme))

         th = xi2 * d.tht1 + (2.0 * jpet - 1.0) * π / d.P - d.tht1

      elseif regloc == 4

         Xx = d.A + d.R * e2x +
              (2.0 * xi2 - 1.0) * d.L1 * e1x / (2.0 * nme)

         Xy = d.B + d.R * e2y +
              (2.0 * xi2 - 1.0) * d.L1 * e1y / (2.0 * nme)

         th = xi2 * (2.0 * π / d.P) + (2.0 * jpet - 3.0) * π / d.P

         Yx = d.A + d.Rc * cos(th)
         Yy = d.B + d.Rc * sin(th)

         Zx = (1.0 - xi1) * Xx + xi1 * Yx
         Zy = (1.0 - xi1) * Xy + xi1 * Yy

         return Zx, Zy

      elseif regloc == 5

         alpx = d.RP.αx[jpet]
         alpy = d.RP.αy[jpet]

         Xx = d.A + (1.0 - xi2) * alpx +
              xi2 * (d.R * e2x - d.L1 * e1x / (2.0 * nme))

         Xy = d.B + (1.0 - xi2) * alpy +
              xi2 * (d.R * e2y - d.L1 * e1y / (2.0 * nme))

         th = xi2 * d.tht1 + (2.0 * jpet - 3.0) * π / d.P

      elseif regloc == 6

         Xx = d.A + d.R * e2x -
              d.L1 * e1x / (2.0 * nme) +
              xi2 * d.L2 * qx / nme

         Xy = d.B + d.R * e2y -
              d.L1 * e1y / (2.0 * nme) +
              xi2 * d.L2 * qy / nme

         th = xi2 * (π / d.P - d.tht2 - d.tht1) +
              (2.0 * jpet - 3.0) * π / d.P + d.tht1

      else

         throw(ArgumentError("mapxy(star): invalid local petal region regloc=$regloc"))

      end

      h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * th) / 2.0

      Yx = d.A + d.R * h * cos(th)
      Yy = d.B + d.R * h * sin(th)

      Zx = (1.0 - xi1) * Xx + xi1 * Yx
      Zy = (1.0 - xi1) * Xy + xi1 * Yy

      return Zx, Zy

   else

      rdisc = reg - 7 * d.P

      if rdisc == 5

         Zx = d.A - d.Lc / 2.0 + d.Lc * xi1
         Zy = d.B - d.Lc / 2.0 + d.Lc * xi2

         return Zx, Zy

      end

      if rdisc == 1

         Xx = d.A + d.Lc / 2.0
         Xy = d.B + d.Lc * (2.0 * xi2 - 1.0) / 2.0

      elseif rdisc == 2

         Xx = d.A + d.Lc * (1.0 - 2.0 * xi2) / 2.0
         Xy = d.B + d.Lc / 2.0

      elseif rdisc == 3

         Xx = d.A - d.Lc / 2.0
         Xy = d.B - d.Lc * (2.0 * xi2 - 1.0) / 2.0

      elseif rdisc == 4

         Xx = d.A - d.Lc * (1.0 - 2.0 * xi2) / 2.0
         Xy = d.B - d.Lc / 2.0

      else

         throw(ArgumentError("mapxy(star): invalid central region rdisc=$rdisc"))

      end

      thp = xi2 / 2.0 + (2.0 * rdisc - 3.0) / 4.0

      Yx = d.A + d.Rc * cospi(thp)
      Yy = d.B + d.Rc * sinpi(thp)

      Zx = (1.0 - xi1) * Xx + xi1 * Yx
      Zy = (1.0 - xi1) * Xy + xi1 * Yy

      return Zx, Zy

   end

end

function mapxy!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
   d::star, u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)

   p = d.pths[k]
   reg = p.reg

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   αc = 0.5 * hc
   αt = 0.5 * ht
   βc = p.ck0 + αc
   βt = p.tk0 + αt

   if reg <= 7 * d.P

      jpet = cld(reg, 7)
      regloc = mod1(reg, 7)

      e1x = d.RP.e₁x[jpet]
      e1y = d.RP.e₁y[jpet]
      e2x = d.RP.e₂x[jpet]
      e2y = d.RP.e₂y[jpet]
      qx = d.RP.qx[jpet]
      qy = d.RP.qy[jpet]
      nme = d.RP.nme

      if regloc == 7

         @inbounds for I in eachindex(Zx, Zy, u, v)

            if d.L2 >= d.L1
               xi1 = muladd(αt, u[I], βt)
               xi2 = muladd(αc, v[I], βc)
            else
               xi1 = muladd(αc, u[I], βc)
               xi2 = muladd(αt, v[I], βt)
            end

            Zx[I] = d.A + d.R * e2x +
                    (2.0 * xi2 - 1.0) * d.L1 * e1x / (2.0 * nme) +
                    d.L2 * xi1 * qx / nme

            Zy[I] = d.B + d.R * e2y +
                    (2.0 * xi2 - 1.0) * d.L1 * e1y / (2.0 * nme) +
                    d.L2 * xi1 * qy / nme

         end

         return nothing

      end

      if regloc == 1

         @inbounds for I in eachindex(Zx, Zy, u, v)

            xi1 = muladd(αc, u[I], βc)
            xi2 = muladd(αt, v[I], βt)

            Xx = d.A + d.R * e2x +
                 (2.0 * xi2 - 1.0) * d.L1 * e1x / (2.0 * nme) +
                 d.L2 * qx / nme

            Xy = d.B + d.R * e2y +
                 (2.0 * xi2 - 1.0) * d.L1 * e1y / (2.0 * nme) +
                 d.L2 * qy / nme

            th = xi2 * (2.0 * d.tht2) +
                 2.0 * (jpet - 1) * π / d.P - d.tht2

            h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * th) / 2.0

            Yx = d.A + d.R * h * cos(th)
            Yy = d.B + d.R * h * sin(th)

            Zx[I] = (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = (1.0 - xi1) * Xy + xi1 * Yy

         end

      elseif regloc == 2

         @inbounds for I in eachindex(Zx, Zy, u, v)

            xi1 = muladd(αc, u[I], βc)
            xi2 = muladd(αt, v[I], βt)

            Xx = d.A + d.R * e2x +
                 d.L1 * e1x / (2.0 * nme) +
                 (1.0 - xi2) * d.L2 * qx / nme

            Xy = d.B + d.R * e2y +
                 d.L1 * e1y / (2.0 * nme) +
                 (1.0 - xi2) * d.L2 * qy / nme

            th = xi2 * (π / d.P - d.tht2 - d.tht1) +
                 2.0 * (jpet - 1) * π / d.P + d.tht2

            h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * th) / 2.0

            Yx = d.A + d.R * h * cos(th)
            Yy = d.B + d.R * h * sin(th)

            Zx[I] = (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = (1.0 - xi1) * Xy + xi1 * Yy

         end

      elseif regloc == 3

         alpx = d.RP.αx[jpet+1]
         alpy = d.RP.αy[jpet+1]

         @inbounds for I in eachindex(Zx, Zy, u, v)

            xi1 = muladd(αc, u[I], βc)
            xi2 = muladd(αt, v[I], βt)

            Xx = d.A + xi2 * alpx +
                 (1.0 - xi2) * (d.R * e2x + d.L1 * e1x / (2.0 * nme))

            Xy = d.B + xi2 * alpy +
                 (1.0 - xi2) * (d.R * e2y + d.L1 * e1y / (2.0 * nme))

            th = xi2 * d.tht1 + (2.0 * jpet - 1.0) * π / d.P - d.tht1

            h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * th) / 2.0

            Yx = d.A + d.R * h * cos(th)
            Yy = d.B + d.R * h * sin(th)

            Zx[I] = (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = (1.0 - xi1) * Xy + xi1 * Yy

         end

      elseif regloc == 4

         @inbounds for I in eachindex(Zx, Zy, u, v)

            xi1 = muladd(αc, u[I], βc)
            xi2 = muladd(αt, v[I], βt)

            Xx = d.A + d.R * e2x +
                 (2.0 * xi2 - 1.0) * d.L1 * e1x / (2.0 * nme)

            Xy = d.B + d.R * e2y +
                 (2.0 * xi2 - 1.0) * d.L1 * e1y / (2.0 * nme)

            th = xi2 * (2.0 * π / d.P) + (2.0 * jpet - 3.0) * π / d.P

            Yx = d.A + d.Rc * cos(th)
            Yy = d.B + d.Rc * sin(th)

            Zx[I] = (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = (1.0 - xi1) * Xy + xi1 * Yy

         end

      elseif regloc == 5

         alpx = d.RP.αx[jpet]
         alpy = d.RP.αy[jpet]

         @inbounds for I in eachindex(Zx, Zy, u, v)

            xi1 = muladd(αc, u[I], βc)
            xi2 = muladd(αt, v[I], βt)

            Xx = d.A + (1.0 - xi2) * alpx +
                 xi2 * (d.R * e2x - d.L1 * e1x / (2.0 * nme))

            Xy = d.B + (1.0 - xi2) * alpy +
                 xi2 * (d.R * e2y - d.L1 * e1y / (2.0 * nme))

            th = xi2 * d.tht1 + (2.0 * jpet - 3.0) * π / d.P

            h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * th) / 2.0

            Yx = d.A + d.R * h * cos(th)
            Yy = d.B + d.R * h * sin(th)

            Zx[I] = (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = (1.0 - xi1) * Xy + xi1 * Yy

         end

      elseif regloc == 6

         @inbounds for I in eachindex(Zx, Zy, u, v)

            xi1 = muladd(αc, u[I], βc)
            xi2 = muladd(αt, v[I], βt)

            Xx = d.A + d.R * e2x -
                 d.L1 * e1x / (2.0 * nme) +
                 xi2 * d.L2 * qx / nme

            Xy = d.B + d.R * e2y -
                 d.L1 * e1y / (2.0 * nme) +
                 xi2 * d.L2 * qy / nme

            th = xi2 * (π / d.P - d.tht2 - d.tht1) +
                 (2.0 * jpet - 3.0) * π / d.P + d.tht1

            h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * th) / 2.0

            Yx = d.A + d.R * h * cos(th)
            Yy = d.B + d.R * h * sin(th)

            Zx[I] = (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = (1.0 - xi1) * Xy + xi1 * Yy

         end

      else

         throw(ArgumentError("mapxy!(star): invalid local petal region regloc=$regloc"))

      end

   else

      rdisc = reg - 7 * d.P

      if rdisc == 5

         @inbounds for I in eachindex(Zx, Zy, u, v)

            xi1 = muladd(αc, u[I], βc)
            xi2 = muladd(αt, v[I], βt)

            Zx[I] = d.A - d.Lc / 2.0 + d.Lc * xi1
            Zy[I] = d.B - d.Lc / 2.0 + d.Lc * xi2

         end

         return nothing

      end

      if rdisc == 1

         @inbounds for I in eachindex(Zx, Zy, u, v)

            xi1 = muladd(αc, u[I], βc)
            xi2 = muladd(αt, v[I], βt)

            Xx = d.A + d.Lc / 2.0
            Xy = d.B + d.Lc * (2.0 * xi2 - 1.0) / 2.0

            thp = xi2 / 2.0 + (2.0 * rdisc - 3.0) / 4.0

            Yx = d.A + d.Rc * cospi(thp)
            Yy = d.B + d.Rc * sinpi(thp)

            Zx[I] = (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = (1.0 - xi1) * Xy + xi1 * Yy

         end

      elseif rdisc == 2

         @inbounds for I in eachindex(Zx, Zy, u, v)

            xi1 = muladd(αc, u[I], βc)
            xi2 = muladd(αt, v[I], βt)

            Xx = d.A + d.Lc * (1.0 - 2.0 * xi2) / 2.0
            Xy = d.B + d.Lc / 2.0

            thp = xi2 / 2.0 + (2.0 * rdisc - 3.0) / 4.0

            Yx = d.A + d.Rc * cospi(thp)
            Yy = d.B + d.Rc * sinpi(thp)

            Zx[I] = (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = (1.0 - xi1) * Xy + xi1 * Yy

         end

      elseif rdisc == 3

         @inbounds for I in eachindex(Zx, Zy, u, v)

            xi1 = muladd(αc, u[I], βc)
            xi2 = muladd(αt, v[I], βt)

            Xx = d.A - d.Lc / 2.0
            Xy = d.B - d.Lc * (2.0 * xi2 - 1.0) / 2.0

            thp = xi2 / 2.0 + (2.0 * rdisc - 3.0) / 4.0

            Yx = d.A + d.Rc * cospi(thp)
            Yy = d.B + d.Rc * sinpi(thp)

            Zx[I] = (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = (1.0 - xi1) * Xy + xi1 * Yy

         end

      elseif rdisc == 4

         @inbounds for I in eachindex(Zx, Zy, u, v)

            xi1 = muladd(αc, u[I], βc)
            xi2 = muladd(αt, v[I], βt)

            Xx = d.A - d.Lc * (1.0 - 2.0 * xi2) / 2.0
            Xy = d.B - d.Lc / 2.0

            thp = xi2 / 2.0 + (2.0 * rdisc - 3.0) / 4.0

            Yx = d.A + d.Rc * cospi(thp)
            Yy = d.B + d.Rc * sinpi(thp)

            Zx[I] = (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = (1.0 - xi1) * Xy + xi1 * Yy

         end

      else

         throw(ArgumentError("mapxy!(star): invalid central region rdisc=$rdisc"))

      end

   end

   return nothing

end

function draw(d::star, flag=nothing; L::Int=33, show::Bool=true)

   clrs = (
      RGBf(0.8500, 0.3250, 0.0980),  # dark orange
      RGBf(0.4940, 0.1840, 0.5560),  # deep purple
      RGBf(0.4660, 0.6740, 0.1880),  # sea green
      RGBf(0.9290, 0.6270, 0.8900),  # light pink
      RGBf(0.3010, 0.7450, 0.9330),  # sky blue
      RGBf(1.0000, 0.7320, 0.4000),  # coral
      RGBf(0.6350, 0.0780, 0.1840),  # maroon
      RGBf(0.9880, 0.4350, 0.7880),  # light purple
      RGBf(0.0000, 0.5000, 0.0000),  # green
      RGBf(0.7700, 0.7060, 0.5220),  # light brown
      RGBf(0.0000, 0.0000, 1.0000),  # blue
      RGBf(0.0000, 0.0000, 1.0000),  # blue
   )

   colors = Vector{RGBf}(undef, 7 * d.P + 5)

   @inbounds for jpet in 1:d.P
      jj = 7 * (jpet - 1)

      colors[jj+1] = clrs[1]
      colors[jj+2] = clrs[2]
      colors[jj+3] = clrs[3]
      colors[jj+4] = clrs[4]
      colors[jj+5] = clrs[5]
      colors[jj+6] = clrs[6]
      colors[jj+7] = clrs[7]
   end

   jj = 7 * d.P + 1
   colors[jj] = clrs[8]
   colors[jj+1] = clrs[9]
   colors[jj+2] = clrs[10]
   colors[jj+3] = clrs[11]
   colors[jj+4] = clrs[12]

   return draw_geom(d, colors; flag=flag, L=L, show=show)
end

#-----------------------
"""
  gamx(d::star, t::Float64, k::Int)
  gamx!(out::StridedArray{Float64}, d::star, t::StridedArray{Float64}, k)
First (x) coordinate of the boundary parametrization γ(t).

- `gamx(d, t, k)` computes the **right boundary** of patch `k`, where `t ∈ [-1,1]`.

`t` may be a scalar or an array; the return has the same shape.
"""
function gamx(d::star, t::Float64, k::Int)::Float64

   p = d.pths[k]
   reg = p.reg

   if reg > 7 * d.P
      throw(ArgumentError("gamx(star): boundary gam is only for petal regions; got reg=$reg"))
   end

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ht = p.tk1 - p.tk0
   ξ = muladd(0.5 * ht, t, 0.5 * (p.tk0 + p.tk1))

   if regloc == 1

      θ = ξ * (2.0 * d.tht2) + 2.0 * (jpet - 1) * π / d.P - d.tht2

   elseif regloc == 2

      θ = ξ * (π / d.P - d.tht2 - d.tht1) + 2.0 * (jpet - 1) * π / d.P + d.tht2

   elseif regloc == 3

      θ = ξ * d.tht1 + (2.0 * jpet - 1.0) * π / d.P - d.tht1

   elseif regloc == 5

      θ = ξ * d.tht1 + (2.0 * jpet - 3.0) * π / d.P

   elseif regloc == 6

      θ = ξ * (π / d.P - d.tht2 - d.tht1) + (2.0 * jpet - 3.0) * π / d.P + d.tht1

   else

      throw(ArgumentError("gamx(star): expected boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"))

   end

   h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

   return d.A + d.R * h * cos(θ)

end

function gamx!(out::StridedArray{Float64}, d::star,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]
   reg = p.reg

   if reg > 7 * d.P
      throw(ArgumentError("gamx!(star): boundary gam is only for petal regions; got reg=$reg"))
   end

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = p.tk0 + αt

   @inbounds for I in eachindex(out, t)

      ξ = muladd(αt, t[I], βt)

      if regloc == 1

         θ = ξ * (2.0 * d.tht2) + 2.0 * (jpet - 1) * π / d.P - d.tht2

      elseif regloc == 2

         θ = ξ * (π / d.P - d.tht2 - d.tht1) + 2.0 * (jpet - 1) * π / d.P + d.tht2

      elseif regloc == 3

         θ = ξ * d.tht1 + (2.0 * jpet - 1.0) * π / d.P - d.tht1

      elseif regloc == 5

         θ = ξ * d.tht1 + (2.0 * jpet - 3.0) * π / d.P

      elseif regloc == 6

         θ = ξ * (π / d.P - d.tht2 - d.tht1) + (2.0 * jpet - 3.0) * π / d.P + d.tht1

      else

         throw(ArgumentError("gamx!(star): expected boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"))

      end

      h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

      out[I] = d.A + d.R * h * cos(θ)

   end

   return nothing

end

"""
   gamy(d::star, t::Float64, k::Int)
   gamy!(out::StridedArray{Float64}, d::star, t::StridedArray{Float64}, k)

Second (y) coordinate of the boundary parametrization γ(t).

- `gamy(d, t, k)` computes the right boundary of patch `k`, where `t ∈ [-1,1]`.

`t` can be a scalar `Float64` or an array; the return has the same shape.
"""
function gamy(d::star, t::Float64, k::Int)::Float64

   p = d.pths[k]
   reg = p.reg

   if reg > 7 * d.P
      throw(ArgumentError("gamy(star): boundary gam is only for petal regions; got reg=$reg"))
   end

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ht = p.tk1 - p.tk0
   ξ = muladd(0.5 * ht, t, 0.5 * (p.tk0 + p.tk1))

   if regloc == 1

      θ = ξ * (2.0 * d.tht2) + 2.0 * (jpet - 1) * π / d.P - d.tht2

   elseif regloc == 2

      θ = ξ * (π / d.P - d.tht2 - d.tht1) + 2.0 * (jpet - 1) * π / d.P + d.tht2

   elseif regloc == 3

      θ = ξ * d.tht1 + (2.0 * jpet - 1.0) * π / d.P - d.tht1

   elseif regloc == 5

      θ = ξ * d.tht1 + (2.0 * jpet - 3.0) * π / d.P

   elseif regloc == 6

      θ = ξ * (π / d.P - d.tht2 - d.tht1) + (2.0 * jpet - 3.0) * π / d.P + d.tht1

   else

      throw(ArgumentError("gamy(star): expected boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"))

   end

   h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

   return d.B + d.R * h * sin(θ)

end

function gamy!(out::StridedArray{Float64}, d::star,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]
   reg = p.reg

   if reg > 7 * d.P
      throw(ArgumentError("gamy!(star): boundary gam is only for petal regions; got reg=$reg"))
   end

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = p.tk0 + αt

   @inbounds for I in eachindex(out, t)

      ξ = muladd(αt, t[I], βt)

      if regloc == 1

         θ = ξ * (2.0 * d.tht2) + 2.0 * (jpet - 1) * π / d.P - d.tht2

      elseif regloc == 2

         θ = ξ * (π / d.P - d.tht2 - d.tht1) + 2.0 * (jpet - 1) * π / d.P + d.tht2

      elseif regloc == 3

         θ = ξ * d.tht1 + (2.0 * jpet - 1.0) * π / d.P - d.tht1

      elseif regloc == 5

         θ = ξ * d.tht1 + (2.0 * jpet - 3.0) * π / d.P

      elseif regloc == 6

         θ = ξ * (π / d.P - d.tht2 - d.tht1) + (2.0 * jpet - 3.0) * π / d.P + d.tht1

      else

         throw(ArgumentError("gamy!(star): expected boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"))

      end

      h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

      out[I] = d.B + d.R * h * sin(θ)

   end

   return nothing

end

function gam(d::star, t::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]
   reg = p.reg

   if reg > 7 * d.P
      throw(ArgumentError("gam(star): boundary gam is only for petal regions; got reg=$reg"))
   end

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ht = p.tk1 - p.tk0
   ξ = muladd(0.5 * ht, t, 0.5 * (p.tk0 + p.tk1))

   if regloc == 1

      θ = ξ * (2.0 * d.tht2) +
          2.0 * (jpet - 1) * π / d.P - d.tht2

   elseif regloc == 2

      θ = ξ * (π / d.P - d.tht2 - d.tht1) +
          2.0 * (jpet - 1) * π / d.P + d.tht2

   elseif regloc == 3

      θ = ξ * d.tht1 +
          (2.0 * jpet - 1.0) * π / d.P - d.tht1

   elseif regloc == 5

      θ = ξ * d.tht1 +
          (2.0 * jpet - 3.0) * π / d.P

   elseif regloc == 6

      θ = ξ * (π / d.P - d.tht2 - d.tht1) +
          (2.0 * jpet - 3.0) * π / d.P + d.tht1

   else

      throw(ArgumentError("gam(star): expected boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"))

   end

   h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

   return d.A + d.R * h * cos(θ), d.B + d.R * h * sin(θ)

end

function gam!(out::Vector{Float64}, d::star, t::Float64, k::Int)

   p = d.pths[k]
   reg = p.reg

   if reg > 7 * d.P
      throw(ArgumentError("gam!(star): boundary gam is only for petal regions; got reg=$reg"))
   end

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ht = p.tk1 - p.tk0
   ξ = muladd(0.5 * ht, t, 0.5 * (p.tk0 + p.tk1))

   if regloc == 1

      θ = ξ * (2.0 * d.tht2) + 2.0 * (jpet - 1) * π / d.P - d.tht2

   elseif regloc == 2

      θ = ξ * (π / d.P - d.tht2 - d.tht1) + 2.0 * (jpet - 1) * π / d.P + d.tht2

   elseif regloc == 3

      θ = ξ * d.tht1 + (2.0 * jpet - 1.0) * π / d.P - d.tht1

   elseif regloc == 5

      θ = ξ * d.tht1 + (2.0 * jpet - 3.0) * π / d.P

   elseif regloc == 6

      θ = ξ * (π / d.P - d.tht2 - d.tht1) + (2.0 * jpet - 3.0) * π / d.P + d.tht1

   else

      throw(ArgumentError("gam!(star): expected boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"))

   end

   h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

   out[1] = d.A + d.R * h * cos(θ)
   out[2] = d.B + d.R * h * sin(θ)

   return nothing

end

function gam!(out::Matrix{Float64}, d::star, t::Vector{Float64}, k::Int)

   p = d.pths[k]
   reg = p.reg

   if reg > 7 * d.P
      throw(ArgumentError("gam!(star): boundary gam is only for petal regions; got reg=$reg"))
   end

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = p.tk0 + αt

   @inbounds for j in eachindex(t)

      ξ = muladd(αt, t[j], βt)

      if regloc == 1

         θ = ξ * (2.0 * d.tht2) + 2.0 * (jpet - 1) * π / d.P - d.tht2

      elseif regloc == 2

         θ = ξ * (π / d.P - d.tht2 - d.tht1) + 2.0 * (jpet - 1) * π / d.P + d.tht2

      elseif regloc == 3

         θ = ξ * d.tht1 + (2.0 * jpet - 1.0) * π / d.P - d.tht1

      elseif regloc == 5

         θ = ξ * d.tht1 + (2.0 * jpet - 3.0) * π / d.P

      elseif regloc == 6

         θ = ξ * (π / d.P - d.tht2 - d.tht1) + (2.0 * jpet - 3.0) * π / d.P + d.tht1

      else

         throw(ArgumentError("gam!(star): expected boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"))

      end

      h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

      out[1, j] = d.A + d.R * h * cos(θ)
      out[2, j] = d.B + d.R * h * sin(θ)

   end

   return nothing

end

function drawbd(d::star, flag=true; L::Int=33, show::Bool=true)

   clrs = (
      RGBf(0.8500, 0.3250, 0.0980),  # dark orange
      RGBf(0.4940, 0.1840, 0.5560),  # deep purple
      RGBf(0.4660, 0.6740, 0.1880),  # sea green
      RGBf(0.9290, 0.6270, 0.8900),  # light pink
      RGBf(0.3010, 0.7450, 0.9330),  # sky blue
      RGBf(1.0000, 0.7320, 0.4000),  # coral
      RGBf(0.6350, 0.0780, 0.1840),  # maroon
      RGBf(0.9880, 0.4350, 0.7880),  # light purple
      RGBf(0.0000, 0.5000, 0.0000),  # green
      RGBf(0.7700, 0.7060, 0.5220),  # light brown
      RGBf(0.0000, 0.0000, 1.0000),  # blue
      RGBf(0.0000, 0.0000, 1.0000),  # blue
   )

   colors = Vector{RGBf}(undef, 7 * d.P + 5)

   @inbounds for jpet in 1:d.P
      jj = 7 * (jpet - 1)

      colors[jj+1] = clrs[1]
      colors[jj+2] = clrs[2]
      colors[jj+3] = clrs[3]
      colors[jj+4] = clrs[4]
      colors[jj+5] = clrs[5]
      colors[jj+6] = clrs[6]
      colors[jj+7] = clrs[7]
   end

   jj = 7 * d.P + 1
   colors[jj] = clrs[8]
   colors[jj+1] = clrs[9]
   colors[jj+2] = clrs[10]
   colors[jj+3] = clrs[11]
   colors[jj+4] = clrs[12]

   return drawbd_geom(d, colors; flag=flag, L=L, show=show)
end

"""
   dgamx(d::star, t::Float64, k::Int) -> Float64

Derivative of the first coordinate of the boundary parametrization.

- With `k`: derivative of the right boundary of patch `k`
  (`∂/∂t mapx(d, 1, t, k)` for `t ∈ [-1,1]`).
"""
function dgamx(d::star, t::Float64, k::Int)::Float64

   p = d.pths[k]
   reg = p.reg

   if reg > 7 * d.P
      throw(ArgumentError("dgamx(star): boundary derivative is only for petal regions; got reg=$reg"))
   end

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ht = p.tk1 - p.tk0
   ξ = muladd(0.5 * ht, t, 0.5 * (p.tk0 + p.tk1))

   if regloc == 1
      dθ = 2.0 * d.tht2
      θ = ξ * dθ + 2.0 * (jpet - 1) * π / d.P - d.tht2
   elseif regloc == 2
      dθ = π / d.P - d.tht2 - d.tht1
      θ = ξ * dθ + 2.0 * (jpet - 1) * π / d.P + d.tht2
   elseif regloc == 3
      dθ = d.tht1
      θ = ξ * dθ + (2.0 * jpet - 1.0) * π / d.P - d.tht1
   elseif regloc == 5
      dθ = d.tht1
      θ = ξ * dθ + (2.0 * jpet - 3.0) * π / d.P
   elseif regloc == 6
      dθ = π / d.P - d.tht2 - d.tht1
      θ = ξ * dθ + (2.0 * jpet - 3.0) * π / d.P + d.tht1
   else
      throw(ArgumentError("dgamx(star): expected boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"))
   end

   s, c = sincos(θ)

   h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0
   h1 = -0.5 * d.Ps * d.P * sin(d.P * θ)

   dgx_dθ = d.R * (h1 * c - h * s)

   return 0.5 * ht * dθ * dgx_dθ

end

function dgamx!(out::StridedArray{Float64}, d::star,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]
   reg = p.reg

   if reg > 7 * d.P
      throw(ArgumentError("dgamx!(star): boundary derivative is only for petal regions; got reg=$reg"))
   end

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = p.tk0 + αt

   if regloc == 1
      dθ = 2.0 * d.tht2
      bθ = 2.0 * (jpet - 1) * π / d.P - d.tht2
   elseif regloc == 2
      dθ = π / d.P - d.tht2 - d.tht1
      bθ = 2.0 * (jpet - 1) * π / d.P + d.tht2
   elseif regloc == 3
      dθ = d.tht1
      bθ = (2.0 * jpet - 1.0) * π / d.P - d.tht1
   elseif regloc == 5
      dθ = d.tht1
      bθ = (2.0 * jpet - 3.0) * π / d.P
   elseif regloc == 6
      dθ = π / d.P - d.tht2 - d.tht1
      bθ = (2.0 * jpet - 3.0) * π / d.P + d.tht1
   else
      throw(ArgumentError("dgamx!(star): expected boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"))
   end

   scale = αt * dθ

   @inbounds for I in eachindex(out, t)
      ξ = muladd(αt, t[I], βt)
      θ = muladd(dθ, ξ, bθ)

      s, c = sincos(θ)

      h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0
      h1 = -0.5 * d.Ps * d.P * sin(d.P * θ)

      out[I] = scale * d.R * (h1 * c - h * s)
   end

   return nothing

end

"""
   dgamy(d::star, t::Float64, k::Int) -> Float64

Derivative of the second coordinate of the boundary parametrization γ(t).

- With `k`: derivative of the right boundary of patch `k`
  (`∂/∂t mapy(d, 1, t, k)` for `t ∈ [-1,1]`).
"""
function dgamy(d::star, t::Float64, k::Int)::Float64

   p = d.pths[k]
   reg = p.reg

   if reg > 7 * d.P
      throw(ArgumentError("dgamy(star): boundary derivative is only for petal regions; got reg=$reg"))
   end

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ht = p.tk1 - p.tk0
   ξ = muladd(0.5 * ht, t, 0.5 * (p.tk0 + p.tk1))

   if regloc == 1
      dθ = 2.0 * d.tht2
      θ = ξ * dθ + 2.0 * (jpet - 1) * π / d.P - d.tht2
   elseif regloc == 2
      dθ = π / d.P - d.tht2 - d.tht1
      θ = ξ * dθ + 2.0 * (jpet - 1) * π / d.P + d.tht2
   elseif regloc == 3
      dθ = d.tht1
      θ = ξ * dθ + (2.0 * jpet - 1.0) * π / d.P - d.tht1
   elseif regloc == 5
      dθ = d.tht1
      θ = ξ * dθ + (2.0 * jpet - 3.0) * π / d.P
   elseif regloc == 6
      dθ = π / d.P - d.tht2 - d.tht1
      θ = ξ * dθ + (2.0 * jpet - 3.0) * π / d.P + d.tht1
   else
      throw(ArgumentError("dgamy(star): expected boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"))
   end

   s, c = sincos(θ)

   h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0
   h1 = -0.5 * d.Ps * d.P * sin(d.P * θ)

   dgy_dθ = d.R * (h1 * s + h * c)

   return 0.5 * ht * dθ * dgy_dθ

end

function dgamy!(out::StridedArray{Float64}, d::star,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]
   reg = p.reg

   if reg > 7 * d.P
      throw(ArgumentError("dgamy!(star): boundary derivative is only for petal regions; got reg=$reg"))
   end

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = p.tk0 + αt

   if regloc == 1
      dθ = 2.0 * d.tht2
      bθ = 2.0 * (jpet - 1) * π / d.P - d.tht2
   elseif regloc == 2
      dθ = π / d.P - d.tht2 - d.tht1
      bθ = 2.0 * (jpet - 1) * π / d.P + d.tht2
   elseif regloc == 3
      dθ = d.tht1
      bθ = (2.0 * jpet - 1.0) * π / d.P - d.tht1
   elseif regloc == 5
      dθ = d.tht1
      bθ = (2.0 * jpet - 3.0) * π / d.P
   elseif regloc == 6
      dθ = π / d.P - d.tht2 - d.tht1
      bθ = (2.0 * jpet - 3.0) * π / d.P + d.tht1
   else
      throw(ArgumentError("dgamy!(star): expected boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"))
   end

   scale = αt * dθ

   @inbounds for I in eachindex(out, t)
      ξ = muladd(αt, t[I], βt)
      θ = muladd(dθ, ξ, bθ)

      s, c = sincos(θ)

      h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0
      h1 = -0.5 * d.Ps * d.P * sin(d.P * θ)

      out[I] = scale * d.R * (h1 * s + h * c)
   end

   return nothing

end

function dgam(d::star, t::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]
   reg = p.reg

   if reg > 7 * d.P
      throw(ArgumentError("dgam(star): boundary derivative is only for petal regions; got reg=$reg"))
   end

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ht = p.tk1 - p.tk0
   ξ = muladd(0.5 * ht, t, 0.5 * (p.tk0 + p.tk1))

   if regloc == 1
      dθ = 2.0 * d.tht2
      θ = ξ * dθ + 2.0 * (jpet - 1) * π / d.P - d.tht2
   elseif regloc == 2
      dθ = π / d.P - d.tht2 - d.tht1
      θ = ξ * dθ + 2.0 * (jpet - 1) * π / d.P + d.tht2
   elseif regloc == 3
      dθ = d.tht1
      θ = ξ * dθ + (2.0 * jpet - 1.0) * π / d.P - d.tht1
   elseif regloc == 5
      dθ = d.tht1
      θ = ξ * dθ + (2.0 * jpet - 3.0) * π / d.P
   elseif regloc == 6
      dθ = π / d.P - d.tht2 - d.tht1
      θ = ξ * dθ + (2.0 * jpet - 3.0) * π / d.P + d.tht1
   else
      throw(ArgumentError("dgam(star): expected boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"))
   end

   s, c = sincos(θ)

   h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0
   h1 = -0.5 * d.Ps * d.P * sin(d.P * θ)

   scale = 0.5 * ht * dθ

   return scale * d.R * (h1 * c - h * s),
          scale * d.R * (h1 * s + h * c)

end

function dgam!(out::Vector{Float64}, d::star, t::Float64, k::Int)

   dx, dy = dgam(d, t, k)

   out[1] = dx
   out[2] = dy

   return nothing

end

function dgam!(out::Matrix{Float64}, d::star, t::Vector{Float64}, k::Int)

   p = d.pths[k]
   reg = p.reg

   if reg > 7 * d.P
      throw(ArgumentError("dgam!(star): boundary derivative is only for petal regions; got reg=$reg"))
   end

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = p.tk0 + αt

   if regloc == 1
      dθ = 2.0 * d.tht2
      bθ = 2.0 * (jpet - 1) * π / d.P - d.tht2
   elseif regloc == 2
      dθ = π / d.P - d.tht2 - d.tht1
      bθ = 2.0 * (jpet - 1) * π / d.P + d.tht2
   elseif regloc == 3
      dθ = d.tht1
      bθ = (2.0 * jpet - 1.0) * π / d.P - d.tht1
   elseif regloc == 5
      dθ = d.tht1
      bθ = (2.0 * jpet - 3.0) * π / d.P
   elseif regloc == 6
      dθ = π / d.P - d.tht2 - d.tht1
      bθ = (2.0 * jpet - 3.0) * π / d.P + d.tht1
   else
      throw(ArgumentError("dgam!(star): expected boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"))
   end

   scale = αt * dθ

   @inbounds for j in eachindex(t)
      ξ = muladd(αt, t[j], βt)
      θ = muladd(dθ, ξ, bθ)

      s, c = sincos(θ)

      h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0
      h1 = -0.5 * d.Ps * d.P * sin(d.P * θ)

      out[1, j] = scale * d.R * (h1 * c - h * s)
      out[2, j] = scale * d.R * (h1 * s + h * c)
   end

   return nothing

end

function gamp(d::star, t::Float64, k::Int)::Tuple{Float64,Float64}

   dx, dy = dgam(d, t, k)

   return dy, -dx

end

function gamp!(out::Matrix{Float64}, d::star, t::Vector{Float64}, k::Int)

   p = d.pths[k]
   reg = p.reg

   if reg > 7 * d.P
      throw(ArgumentError("gamp!(star): boundary derivative is only for petal regions; got reg=$reg"))
   end

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = p.tk0 + αt

   if regloc == 1
      dθ = 2.0 * d.tht2
      bθ = 2.0 * (jpet - 1) * π / d.P - d.tht2
   elseif regloc == 2
      dθ = π / d.P - d.tht2 - d.tht1
      bθ = 2.0 * (jpet - 1) * π / d.P + d.tht2
   elseif regloc == 3
      dθ = d.tht1
      bθ = (2.0 * jpet - 1.0) * π / d.P - d.tht1
   elseif regloc == 5
      dθ = d.tht1
      bθ = (2.0 * jpet - 3.0) * π / d.P
   elseif regloc == 6
      dθ = π / d.P - d.tht2 - d.tht1
      bθ = (2.0 * jpet - 3.0) * π / d.P + d.tht1
   else
      throw(ArgumentError("gamp!(star): expected boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"))
   end

   scale = αt * dθ

   @inbounds for j in eachindex(t)
      ξ = muladd(αt, t[j], βt)
      θ = muladd(dθ, ξ, bθ)

      s, c = sincos(θ)

      h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0
      h1 = -0.5 * d.Ps * d.P * sin(d.P * θ)

      dx = scale * d.R * (h1 * c - h * s)
      dy = scale * d.R * (h1 * s + h * c)

      out[1, j] = dy
      out[2, j] = -dx
   end

   return nothing

end

function nu!(out::Matrix{Float64}, d::star, t::Vector{Float64}, k::Int)

   p = d.pths[k]
   reg = p.reg

   if reg > 7 * d.P
      throw(ArgumentError("nu!(star): boundary derivative is only for petal regions; got reg=$reg"))
   end

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = p.tk0 + αt

   if regloc == 1
      dθ = 2.0 * d.tht2
      bθ = 2.0 * (jpet - 1) * π / d.P - d.tht2
   elseif regloc == 2
      dθ = π / d.P - d.tht2 - d.tht1
      bθ = 2.0 * (jpet - 1) * π / d.P + d.tht2
   elseif regloc == 3
      dθ = d.tht1
      bθ = (2.0 * jpet - 1.0) * π / d.P - d.tht1
   elseif regloc == 5
      dθ = d.tht1
      bθ = (2.0 * jpet - 3.0) * π / d.P
   elseif regloc == 6
      dθ = π / d.P - d.tht2 - d.tht1
      bθ = (2.0 * jpet - 3.0) * π / d.P + d.tht1
   else
      throw(ArgumentError("nu!(star): expected boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"))
   end

   scale = αt * dθ

   @inbounds for j in eachindex(t)
      ξ = muladd(αt, t[j], βt)
      θ = muladd(dθ, ξ, bθ)

      s, c = sincos(θ)

      h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0
      h1 = -0.5 * d.Ps * d.P * sin(d.P * θ)

      dx = scale * d.R * (h1 * c - h * s)
      dy = scale * d.R * (h1 * s + h * c)

      nx = dy
      ny = -dx
      nrm = hypot(nx, ny)

      out[1, j] = nx / nrm
      out[2, j] = ny / nrm
   end

   return nothing

end

# This function finds s such that γ_l(t) = γ_k(s),
# allowing s to lie outside [-1, 1].
function bdinv(d::star, t::Float64, l::Int, k::Int)::Float64

   @inbounds begin

      pl = d.pths[l]
      pk = d.pths[k]

      plr = pl.reg
      pkr = pk.reg

      jsrc = cld(plr, 7)
      jtar = cld(pkr, 7)

      rsrc = mod1(plr, 7)
      rtar = mod1(pkr, 7)

      # Map local t on patch l to its ξ coordinate.
      xil = muladd(0.5 * (pl.tk1 - pl.tk0), t, 0.5 * (pl.tk0 + pl.tk1))

      # Same physical angular region: same ξ, converted to target patch.
      if plr == pkr
         return xi_inv(xil, pk.tk0, pk.tk1)
      end

      # Source angle θ = as * ξ + bs.
      th = if rsrc == 1

         xil * (2.0 * d.tht2) + 2.0 * (jsrc - 1) * π / d.P - d.tht2

      elseif rsrc == 2

         xil * (π / d.P - d.tht2 - d.tht1) + 2.0 * (jsrc - 1) * π / d.P + d.tht2

      elseif rsrc == 3

         xil * d.tht1 + (2.0 * jsrc - 1.0) * π / d.P - d.tht1

      elseif rsrc == 5

         xil * d.tht1 + (2.0 * jsrc - 3.0) * π / d.P

      else # rsrc == 6

         xil * (π / d.P - d.tht2 - d.tht1) + (2.0 * jsrc - 3.0) * π / d.P + d.tht1

      end

      # Target angle map θ = ak * ξ + bk.
      ak, bk = if rtar == 1

         2.0 * d.tht2, 2.0 * (jtar - 1) * π / d.P - d.tht2

      elseif rtar == 2

         π / d.P - d.tht2 - d.tht1, 2.0 * (jtar - 1) * π / d.P + d.tht2

      elseif rtar == 3

         d.tht1, (2.0 * jtar - 1.0) * π / d.P - d.tht1

      elseif rtar == 5

         d.tht1, (2.0 * jtar - 3.0) * π / d.P

      else # rtar == 6

         π / d.P - d.tht2 - d.tht1, (2.0 * jtar - 3.0) * π / d.P + d.tht1

      end

      # Shift source θ by multiples of 2π closest to target angular interval.
      center = bk + 0.5 * ak
      m = round((center - th) / (2π))
      ths = th + 2π * m

      # Convert shifted angle to target ξ, then to local patch coordinate.
      xik = (ths - bk) / ak

      return xi_inv(xik, pk.tk0, pk.tk1)

   end

end

"""
   hderhigher(d::star, θ::Float64)

Compute h, h', ..., h⁽⁶⁾ for the star radial function

   h(θ) = 1 + Ps/2 + Ps*cos(Pθ)/2

Returns h, h1, h2, h3, h4, h5, h6
"""
function hderhigher(d::star, θ::Float64)

   P = Float64(d.P)
   a = 0.5 * d.Ps

   sP, cP = sincos(P * θ)

   P2 = P * P
   P3 = P2 * P
   P4 = P2 * P2
   P5 = P4 * P
   P6 = P3 * P3

   h  = 1.0 + a + a * cP
   h1 = -a * P  * sP
   h2 = -a * P2 * cP
   h3 =  a * P3 * sP
   h4 =  a * P4 * cP
   h5 = -a * P5 * sP
   h6 = -a * P6 * cP

   return h, h1, h2, h3, h4, h5, h6

end

"""
   DLP!(out::StridedArray{Float64}, d::star, t::Float64,
   tau::StridedArray{Float64}, k::Int,
   x::Vector{Float64}, G::Vector{Float64}, GP::Vector{Float64})

Same-panel double-layer kernel on the star boundary:

   K(t, τ) =((γ_k(τ) - γ_k(t)) ⋅ γᵖᵉʳᵖ_k(τ)) / ||γ_k(τ) - γ_k(t)||²

Uses the polar-form cancellation for γ(θ) = R h(θ)(cosθ, sinθ).
"""
function DLP!(out::StridedArray{Float64}, d::star, t::Float64,
   tau::StridedArray{Float64}, k::Int, x::Vector{Float64},
   G::Vector{Float64}, GP::Vector{Float64})

   p = d.pths[k]
   reg = p.reg

   jpet = cld(reg, 7)
   regloc = mod1(reg, 7)

   ak, bk = if regloc == 1

      2.0 * d.tht2,
      2.0 * (jpet - 1) * π / d.P - d.tht2

   elseif regloc == 2

      π / d.P - d.tht2 - d.tht1,
      2.0 * (jpet - 1) * π / d.P + d.tht2

   elseif regloc == 3

      d.tht1,
      (2.0 * jpet - 1.0) * π / d.P - d.tht1

   elseif regloc == 5

      d.tht1,
      (2.0 * jpet - 3.0) * π / d.P

   elseif regloc == 6

      π / d.P - d.tht2 - d.tht1,
      (2.0 * jpet - 3.0) * π / d.P + d.tht1

   else

      throw(ArgumentError(
         "DLP! star expects boundary local region 1,2,3,5,6; got reg=$reg, regloc=$regloc"
      ))

   end

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = p.tk0 + αt

   # t mapped from [-1,1] to [tk0,tk1], then to θ.
   tm = muladd(αt, t, βt)
   tht = muladd(ak, tm, bk)

   ht0, ht1, ht2, ht3, ht4, ht5, ht6 = hderhigher(d, tht)

   # dθ / dτ_reference
   λ = αt * ak

   @inbounds for i in eachindex(tau, out)

      τm = muladd(αt, tau[i], βt)
      thτ = muladd(ak, τm, bk)

      sPτ, cPτ = sincos(d.P * thτ)

      hτ = 1.0 + d.Ps / 2.0 + d.Ps * cPτ / 2.0
      hτ1 = -0.5 * d.Ps * d.P * sPτ

      Δ = 0.5 * (thτ - tht)

      sΔ, cΔ = sincos(Δ)

      if abs(Δ) < 2e-3

         Δ2 = Δ * Δ
         Δ3 = Δ2 * Δ
         Δ4 = Δ2 * Δ2

         q = ht1 + Δ * ht2 +
             Δ2 * ((2.0 / 3.0) * ht3 + ht1 / 6.0) +
             Δ3 * (ht4 / 3.0 + ht2 / 6.0) +
             Δ4 * ((2.0 / 15.0) * ht5 + ht3 / 9.0 + (7.0 / 360.0) * ht1)

         q2 = -ht2 +
              Δ * (-(4.0 / 3.0) * ht3 + (2.0 / 3.0) * ht1) +
              Δ2 * (-ht4 + ht2) +
              Δ3 * (-(8.0 / 15.0) * ht5 + (8.0 / 9.0) * ht3 + (4.0 / 45.0) * ht1) +
              Δ4 * (-(2.0 / 9.0) * ht6 + (5.0 / 9.0) * ht4 + ht2 / 9.0)

      else

         q = (hτ - ht0) / (2.0 * sΔ)
         q2 = (q - hτ1 * cΔ) / sΔ

      end

      num = ht0 * hτ + hτ * q2 + 2.0 * q * hτ1 * cΔ
      den = q * q + hτ * ht0

      out[i] = 0.5 * λ * num / den

   end

   return nothing

end

#-----------------------
function dergam(d::star, t::Float64)::Tuple{Float64,Float64,Float64,Float64}

   Pt = d.P * t

   h = 1.0 + d.Ps / 2.0 + d.Ps * cos(Pt) / 2.0
   dh = -d.Ps * d.P * sin(Pt) / 2.0

   st, ct = sincos(t)

   gx = d.R * h * ct
   gy = d.R * h * st

   dgx = d.R * (dh * ct - h * st)
   dgy = d.R * (dh * st + h * ct)

   return gx, gy, dgx, dgy

end

function Dmap!(out::StridedArray{Float64}, d::star,
   u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)

   p = d.pths[k]
   reg = p.reg

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   αu = 0.5 * hc
   αv = 0.5 * ht
   βu = p.ck0 + αu
   βv = p.tk0 + αv

   if reg <= 7 * d.P

      jpet = cld(reg, 7)
      regloc = mod1(reg, 7)

      e1x = d.RP.e₁x[jpet]
      e1y = d.RP.e₁y[jpet]
      e2x = d.RP.e₂x[jpet]
      e2y = d.RP.e₂y[jpet]
      qx = d.RP.qx[jpet]
      qy = d.RP.qy[jpet]
      nme = d.RP.nme

      if regloc == 7

         fill!(out, d.L1 * d.L2 * hc * ht / 4.0)

         return nothing

      end

      if regloc == 1

         @inbounds for I in eachindex(out, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = xi2 * (2.0 * d.tht2) +
                2.0 * (jpet - 1) * π / d.P - d.tht2

            gx, gy, dgx, dgy = dergam(d, θ)

            Xx = d.R * e2x +
                 (2.0 * xi2 - 1.0) * d.L1 * e1x / (2.0 * nme) +
                 d.L2 * qx / nme

            Xy = d.R * e2y +
                 (2.0 * xi2 - 1.0) * d.L1 * e1y / (2.0 * nme) +
                 d.L2 * qy / nme

            dux = αu * (gx - Xx)
            duy = αu * (gy - Xy)

            dvx = ht * d.L1 * (1.0 - xi1) * e1x / (2.0 * nme) +
                  ht * d.tht2 * xi1 * dgx

            dvy = ht * d.L1 * (1.0 - xi1) * e1y / (2.0 * nme) +
                  ht * d.tht2 * xi1 * dgy

            out[I] = abs(dux * dvy - dvx * duy)

         end

      elseif regloc == 2

         Δθ = π / d.P - d.tht2 - d.tht1

         @inbounds for I in eachindex(out, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = xi2 * Δθ + 2.0 * (jpet - 1) * π / d.P + d.tht2

            gx, gy, dgx, dgy = dergam(d, θ)

            Xx = d.R * e2x +
                 d.L1 * e1x / (2.0 * nme) +
                 (1.0 - xi2) * d.L2 * qx / nme

            Xy = d.R * e2y +
                 d.L1 * e1y / (2.0 * nme) +
                 (1.0 - xi2) * d.L2 * qy / nme

            dux = αu * (gx - Xx)
            duy = αu * (gy - Xy)

            dvx = -ht * d.L2 * (1.0 - xi1) * qx / (2.0 * nme) +
                  ht * Δθ * xi1 * dgx / 2.0

            dvy = -ht * d.L2 * (1.0 - xi1) * qy / (2.0 * nme) +
                  ht * Δθ * xi1 * dgy / 2.0

            out[I] = abs(dux * dvy - dvx * duy)

         end

      elseif regloc == 3

         alpx = d.RP.αx[jpet+1]
         alpy = d.RP.αy[jpet+1]

         Cx = d.R * e2x + d.L1 * e1x / (2.0 * nme)
         Cy = d.R * e2y + d.L1 * e1y / (2.0 * nme)

         @inbounds for I in eachindex(out, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = xi2 * d.tht1 +
                (2.0 * jpet - 1.0) * π / d.P - d.tht1

            gx, gy, dgx, dgy = dergam(d, θ)

            Xx = xi2 * alpx + (1.0 - xi2) * Cx
            Xy = xi2 * alpy + (1.0 - xi2) * Cy

            dux = αu * (gx - Xx)
            duy = αu * (gy - Xy)

            dvx = ht * (1.0 - xi1) * (alpx - Cx) / 2.0 +
                  ht * d.tht1 * xi1 * dgx / 2.0

            dvy = ht * (1.0 - xi1) * (alpy - Cy) / 2.0 +
                  ht * d.tht1 * xi1 * dgy / 2.0

            out[I] = abs(dux * dvy - dvx * duy)

         end

      elseif regloc == 4

         @inbounds for I in eachindex(out, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = xi2 * (2.0 * π / d.P) +
                (2.0 * jpet - 3.0) * π / d.P

            st, ct = sincos(θ)

            Rcct = d.Rc * ct
            Rcst = d.Rc * st

            Xx = d.R * e2x +
                 (2.0 * xi2 - 1.0) * d.L1 * e1x / (2.0 * nme)

            Xy = d.R * e2y +
                 (2.0 * xi2 - 1.0) * d.L1 * e1y / (2.0 * nme)

            dux = αu * (Rcct - Xx)
            duy = αu * (Rcst - Xy)

            dvx = ht * d.L1 * (1.0 - xi1) * e1x / (2.0 * nme) -
                  π * ht * d.Rc * xi1 * st / d.P

            dvy = ht * d.L1 * (1.0 - xi1) * e1y / (2.0 * nme) +
                  π * ht * d.Rc * xi1 * ct / d.P

            out[I] = abs(dux * dvy - dvx * duy)

         end

      elseif regloc == 5

         alpx = d.RP.αx[jpet]
         alpy = d.RP.αy[jpet]

         Cx = d.R * e2x - d.L1 * e1x / (2.0 * nme)
         Cy = d.R * e2y - d.L1 * e1y / (2.0 * nme)

         @inbounds for I in eachindex(out, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = xi2 * d.tht1 +
                (2.0 * jpet - 3.0) * π / d.P

            gx, gy, dgx, dgy = dergam(d, θ)

            Xx = (1.0 - xi2) * alpx + xi2 * Cx
            Xy = (1.0 - xi2) * alpy + xi2 * Cy

            dux = αu * (gx - Xx)
            duy = αu * (gy - Xy)

            dvx = ht * (1.0 - xi1) * (Cx - alpx) / 2.0 +
                  ht * d.tht1 * xi1 * dgx / 2.0

            dvy = ht * (1.0 - xi1) * (Cy - alpy) / 2.0 +
                  ht * d.tht1 * xi1 * dgy / 2.0

            out[I] = abs(dux * dvy - dvx * duy)

         end

      elseif regloc == 6

         Δθ = π / d.P - d.tht2 - d.tht1

         @inbounds for I in eachindex(out, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = xi2 * Δθ +
                (2.0 * jpet - 3.0) * π / d.P + d.tht1

            gx, gy, dgx, dgy = dergam(d, θ)

            Xx = d.R * e2x -
                 d.L1 * e1x / (2.0 * nme) +
                 xi2 * d.L2 * qx / nme

            Xy = d.R * e2y -
                 d.L1 * e1y / (2.0 * nme) +
                 xi2 * d.L2 * qy / nme

            dux = αu * (gx - Xx)
            duy = αu * (gy - Xy)

            dvx = ht * d.L2 * (1.0 - xi1) * qx / (2.0 * nme) +
                  ht * Δθ * xi1 * dgx / 2.0

            dvy = ht * d.L2 * (1.0 - xi1) * qy / (2.0 * nme) +
                  ht * Δθ * xi1 * dgy / 2.0

            out[I] = abs(dux * dvy - dvx * duy)

         end

      else

         throw(ArgumentError("Dmap!(star): invalid local petal region regloc=$regloc"))

      end

   else

      rdisc = reg - 7 * d.P

      if rdisc == 5

         fill!(out, d.Lc^2 * hc * ht / 4.0)

         return nothing

      end

      if rdisc == 1

         @inbounds for I in eachindex(out, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = π * xi2 / 2.0 - π / 4.0

            st, ct = sincos(θ)

            svt = d.Rc * st
            cvt = d.Rc * ct

            dux = αu * (cvt - d.Lc / 2.0)
            duy = αu * (svt - d.Lc * (2.0 * xi2 - 1.0) / 2.0)

            dvx = -π * ht * xi1 * svt / 4.0
            dvy = (1.0 - xi1) * d.Lc * ht / 2.0 +
                  π * ht * xi1 * cvt / 4.0

            out[I] = abs(dux * dvy - dvx * duy)

         end

      elseif rdisc == 2

         @inbounds for I in eachindex(out, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = π * xi2 / 2.0 + π / 4.0

            st, ct = sincos(θ)

            svt = d.Rc * st
            cvt = d.Rc * ct

            dux = αu * (cvt - d.Lc * (1.0 - 2.0 * xi2) / 2.0)
            duy = αu * (svt - d.Lc / 2.0)

            dvx = -(1.0 - xi1) * d.Lc * ht / 2.0 -
                  π * ht * xi1 * svt / 4.0
            dvy = π * ht * xi1 * cvt / 4.0

            out[I] = abs(dux * dvy - dvx * duy)

         end

      elseif rdisc == 3

         @inbounds for I in eachindex(out, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = π * xi2 / 2.0 + 3.0 * π / 4.0

            st, ct = sincos(θ)

            svt = d.Rc * st
            cvt = d.Rc * ct

            dux = αu * (cvt + d.Lc / 2.0)
            duy = αu * (svt + d.Lc * (2.0 * xi2 - 1.0) / 2.0)

            dvx = -π * ht * xi1 * svt / 4.0
            dvy = (xi1 - 1.0) * d.Lc * ht / 2.0 +
                  π * ht * xi1 * cvt / 4.0

            out[I] = abs(dux * dvy - dvx * duy)

         end

      elseif rdisc == 4

         @inbounds for I in eachindex(out, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = π * xi2 / 2.0 + 5.0 * π / 4.0

            st, ct = sincos(θ)

            svt = d.Rc * st
            cvt = d.Rc * ct

            dux = αu * (cvt + d.Lc * (1.0 - 2.0 * xi2) / 2.0)
            duy = αu * (svt + d.Lc / 2.0)

            dvx = (1.0 - xi1) * d.Lc * ht / 2.0 -
                  π * ht * xi1 * svt / 4.0
            dvy = π * ht * xi1 * cvt / 4.0

            out[I] = abs(dux * dvy - dvx * duy)

         end

      else

         throw(ArgumentError("Dmap!(star): invalid central region rdisc=$rdisc"))

      end

   end

   return nothing

end

"""
A combination of mapxy! and Dmap! function. No allocations.

The purpose of this function is to reduce repeated computations
related to cosine and sine.
"""
function mapxy_Dmap!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
   DJ::StridedArray{Float64}, d::star,
   u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)

   p = d.pths[k]
   reg = p.reg

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   αu = 0.5 * hc
   αv = 0.5 * ht
   βu = p.ck0 + αu
   βv = p.tk0 + αv

   if reg <= 7 * d.P

      jpet = cld(reg, 7)
      regloc = mod1(reg, 7)

      e1x = d.RP.e₁x[jpet]
      e1y = d.RP.e₁y[jpet]
      e2x = d.RP.e₂x[jpet]
      e2y = d.RP.e₂y[jpet]
      qx = d.RP.qx[jpet]
      qy = d.RP.qy[jpet]
      nme = d.RP.nme

      invn = 1.0 / nme

      # ------------------------------------------------------------
      # Petal rectangular region: affine map.
      # ------------------------------------------------------------
      if regloc == 7

         fill!(DJ, d.L1 * d.L2 * hc * ht / 4.0)

         @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

            if d.L2 >= d.L1
               xi1 = muladd(αv, u[I], βv)
               xi2 = muladd(αu, v[I], βu)
            else
               xi1 = muladd(αu, u[I], βu)
               xi2 = muladd(αv, v[I], βv)
            end

            Zx[I] = d.A + d.R * e2x +
                    (2.0 * xi2 - 1.0) * d.L1 * e1x * invn / 2.0 +
                    d.L2 * xi1 * qx * invn

            Zy[I] = d.B + d.R * e2y +
                    (2.0 * xi2 - 1.0) * d.L1 * e1y * invn / 2.0 +
                    d.L2 * xi1 * qy * invn

         end

         return nothing

      end

      # ------------------------------------------------------------
      # Petal curved/interpolated regions.
      # ------------------------------------------------------------
      if regloc == 1

         ak = 2.0 * d.tht2
         bk = 2.0 * (jpet - 1) * π / d.P - d.tht2

         @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = muladd(ak, xi2, bk)
            s, c = sincos(θ)
            sP, cP = sincos(d.P * θ)

            h = 1.0 + d.Ps / 2.0 + d.Ps * cP / 2.0
            h1 = -0.5 * d.Ps * d.P * sP

            Yx = d.R * h * c
            Yy = d.R * h * s

            dYx = d.R * (h1 * c - h * s)
            dYy = d.R * (h1 * s + h * c)

            Xx = d.R * e2x +
                 (2.0 * xi2 - 1.0) * d.L1 * e1x * invn / 2.0 +
                 d.L2 * qx * invn

            Xy = d.R * e2y +
                 (2.0 * xi2 - 1.0) * d.L1 * e1y * invn / 2.0 +
                 d.L2 * qy * invn

            dXx = d.L1 * e1x * invn
            dXy = d.L1 * e1y * invn

            Zx[I] = d.A + (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = d.B + (1.0 - xi1) * Xy + xi1 * Yy

            dux = αu * (Yx - Xx)
            duy = αu * (Yy - Xy)

            dvx = αv * ((1.0 - xi1) * dXx + xi1 * ak * dYx)
            dvy = αv * ((1.0 - xi1) * dXy + xi1 * ak * dYy)

            DJ[I] = abs(dux * dvy - dvx * duy)

         end

      elseif regloc == 2

         ak = π / d.P - d.tht2 - d.tht1
         bk = 2.0 * (jpet - 1) * π / d.P + d.tht2

         @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = muladd(ak, xi2, bk)
            s, c = sincos(θ)
            sP, cP = sincos(d.P * θ)

            h = 1.0 + d.Ps / 2.0 + d.Ps * cP / 2.0
            h1 = -0.5 * d.Ps * d.P * sP

            Yx = d.R * h * c
            Yy = d.R * h * s

            dYx = d.R * (h1 * c - h * s)
            dYy = d.R * (h1 * s + h * c)

            Xx = d.R * e2x +
                 d.L1 * e1x * invn / 2.0 +
                 (1.0 - xi2) * d.L2 * qx * invn

            Xy = d.R * e2y +
                 d.L1 * e1y * invn / 2.0 +
                 (1.0 - xi2) * d.L2 * qy * invn

            dXx = -d.L2 * qx * invn
            dXy = -d.L2 * qy * invn

            Zx[I] = d.A + (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = d.B + (1.0 - xi1) * Xy + xi1 * Yy

            dux = αu * (Yx - Xx)
            duy = αu * (Yy - Xy)

            dvx = αv * ((1.0 - xi1) * dXx + xi1 * ak * dYx)
            dvy = αv * ((1.0 - xi1) * dXy + xi1 * ak * dYy)

            DJ[I] = abs(dux * dvy - dvx * duy)

         end

      elseif regloc == 3

         ak = d.tht1
         bk = (2.0 * jpet - 1.0) * π / d.P - d.tht1

         alpx = d.RP.αx[jpet + 1]
         alpy = d.RP.αy[jpet + 1]

         Cx = d.R * e2x + d.L1 * e1x * invn / 2.0
         Cy = d.R * e2y + d.L1 * e1y * invn / 2.0

         dXx_const = alpx - Cx
         dXy_const = alpy - Cy

         @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = muladd(ak, xi2, bk)
            s, c = sincos(θ)
            sP, cP = sincos(d.P * θ)

            h = 1.0 + d.Ps / 2.0 + d.Ps * cP / 2.0
            h1 = -0.5 * d.Ps * d.P * sP

            Yx = d.R * h * c
            Yy = d.R * h * s

            dYx = d.R * (h1 * c - h * s)
            dYy = d.R * (h1 * s + h * c)

            Xx = xi2 * alpx + (1.0 - xi2) * Cx
            Xy = xi2 * alpy + (1.0 - xi2) * Cy

            Zx[I] = d.A + (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = d.B + (1.0 - xi1) * Xy + xi1 * Yy

            dux = αu * (Yx - Xx)
            duy = αu * (Yy - Xy)

            dvx = αv * ((1.0 - xi1) * dXx_const + xi1 * ak * dYx)
            dvy = αv * ((1.0 - xi1) * dXy_const + xi1 * ak * dYy)

            DJ[I] = abs(dux * dvy - dvx * duy)

         end

      elseif regloc == 4

         ak = 2.0 * π / d.P
         bk = (2.0 * jpet - 3.0) * π / d.P

         @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = muladd(ak, xi2, bk)
            s, c = sincos(θ)

            Yx = d.Rc * c
            Yy = d.Rc * s

            dYx = -d.Rc * s
            dYy =  d.Rc * c

            Xx = d.R * e2x +
                 (2.0 * xi2 - 1.0) * d.L1 * e1x * invn / 2.0

            Xy = d.R * e2y +
                 (2.0 * xi2 - 1.0) * d.L1 * e1y * invn / 2.0

            dXx = d.L1 * e1x * invn
            dXy = d.L1 * e1y * invn

            Zx[I] = d.A + (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = d.B + (1.0 - xi1) * Xy + xi1 * Yy

            dux = αu * (Yx - Xx)
            duy = αu * (Yy - Xy)

            dvx = αv * ((1.0 - xi1) * dXx + xi1 * ak * dYx)
            dvy = αv * ((1.0 - xi1) * dXy + xi1 * ak * dYy)

            DJ[I] = abs(dux * dvy - dvx * duy)

         end

      elseif regloc == 5

         ak = d.tht1
         bk = (2.0 * jpet - 3.0) * π / d.P

         alpx = d.RP.αx[jpet]
         alpy = d.RP.αy[jpet]

         Cx = d.R * e2x - d.L1 * e1x * invn / 2.0
         Cy = d.R * e2y - d.L1 * e1y * invn / 2.0

         dXx_const = Cx - alpx
         dXy_const = Cy - alpy

         @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = muladd(ak, xi2, bk)
            s, c = sincos(θ)
            sP, cP = sincos(d.P * θ)

            h = 1.0 + d.Ps / 2.0 + d.Ps * cP / 2.0
            h1 = -0.5 * d.Ps * d.P * sP

            Yx = d.R * h * c
            Yy = d.R * h * s

            dYx = d.R * (h1 * c - h * s)
            dYy = d.R * (h1 * s + h * c)

            Xx = (1.0 - xi2) * alpx + xi2 * Cx
            Xy = (1.0 - xi2) * alpy + xi2 * Cy

            Zx[I] = d.A + (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = d.B + (1.0 - xi1) * Xy + xi1 * Yy

            dux = αu * (Yx - Xx)
            duy = αu * (Yy - Xy)

            dvx = αv * ((1.0 - xi1) * dXx_const + xi1 * ak * dYx)
            dvy = αv * ((1.0 - xi1) * dXy_const + xi1 * ak * dYy)

            DJ[I] = abs(dux * dvy - dvx * duy)

         end

      elseif regloc == 6

         ak = π / d.P - d.tht2 - d.tht1
         bk = (2.0 * jpet - 3.0) * π / d.P + d.tht1

         @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = muladd(ak, xi2, bk)
            s, c = sincos(θ)
            sP, cP = sincos(d.P * θ)

            h = 1.0 + d.Ps / 2.0 + d.Ps * cP / 2.0
            h1 = -0.5 * d.Ps * d.P * sP

            Yx = d.R * h * c
            Yy = d.R * h * s

            dYx = d.R * (h1 * c - h * s)
            dYy = d.R * (h1 * s + h * c)

            Xx = d.R * e2x -
                 d.L1 * e1x * invn / 2.0 +
                 xi2 * d.L2 * qx * invn

            Xy = d.R * e2y -
                 d.L1 * e1y * invn / 2.0 +
                 xi2 * d.L2 * qy * invn

            dXx = d.L2 * qx * invn
            dXy = d.L2 * qy * invn

            Zx[I] = d.A + (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = d.B + (1.0 - xi1) * Xy + xi1 * Yy

            dux = αu * (Yx - Xx)
            duy = αu * (Yy - Xy)

            dvx = αv * ((1.0 - xi1) * dXx + xi1 * ak * dYx)
            dvy = αv * ((1.0 - xi1) * dXy + xi1 * ak * dYy)

            DJ[I] = abs(dux * dvy - dvx * duy)

         end

      else

         throw(ArgumentError("mapxy_Dmap!(star): invalid local petal region regloc=$regloc"))

      end

   else

      # ------------------------------------------------------------
      # Central disc-like regions.
      # ------------------------------------------------------------
      rdisc = reg - 7 * d.P

      if rdisc == 5

         fill!(DJ, d.Lc^2 * hc * ht / 4.0)

         @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            Zx[I] = d.A - d.Lc / 2.0 + d.Lc * xi1
            Zy[I] = d.B - d.Lc / 2.0 + d.Lc * xi2

         end

         return nothing

      elseif rdisc == 1

         @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = π * xi2 / 2.0 - π / 4.0
            s, c = sincos(θ)

            Yx = d.Rc * c
            Yy = d.Rc * s

            Xx = d.Lc / 2.0
            Xy = d.Lc * (2.0 * xi2 - 1.0) / 2.0

            dYx = -d.Rc * s
            dYy =  d.Rc * c

            dXx = 0.0
            dXy = d.Lc

            Zx[I] = d.A + (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = d.B + (1.0 - xi1) * Xy + xi1 * Yy

            dux = αu * (Yx - Xx)
            duy = αu * (Yy - Xy)

            dvx = αv * (xi1 * (π / 2.0) * dYx)
            dvy = αv * ((1.0 - xi1) * dXy + xi1 * (π / 2.0) * dYy)

            DJ[I] = abs(dux * dvy - dvx * duy)

         end

      elseif rdisc == 2

         @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = π * xi2 / 2.0 + π / 4.0
            s, c = sincos(θ)

            Yx = d.Rc * c
            Yy = d.Rc * s

            Xx = d.Lc * (1.0 - 2.0 * xi2) / 2.0
            Xy = d.Lc / 2.0

            dYx = -d.Rc * s
            dYy =  d.Rc * c

            dXx = -d.Lc
            dXy = 0.0

            Zx[I] = d.A + (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = d.B + (1.0 - xi1) * Xy + xi1 * Yy

            dux = αu * (Yx - Xx)
            duy = αu * (Yy - Xy)

            dvx = αv * ((1.0 - xi1) * dXx + xi1 * (π / 2.0) * dYx)
            dvy = αv * (xi1 * (π / 2.0) * dYy)

            DJ[I] = abs(dux * dvy - dvx * duy)

         end

      elseif rdisc == 3

         @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = π * xi2 / 2.0 + 3.0 * π / 4.0
            s, c = sincos(θ)

            Yx = d.Rc * c
            Yy = d.Rc * s

            Xx = -d.Lc / 2.0
            Xy = -d.Lc * (2.0 * xi2 - 1.0) / 2.0

            dYx = -d.Rc * s
            dYy =  d.Rc * c

            dXx = 0.0
            dXy = -d.Lc

            Zx[I] = d.A + (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = d.B + (1.0 - xi1) * Xy + xi1 * Yy

            dux = αu * (Yx - Xx)
            duy = αu * (Yy - Xy)

            dvx = αv * (xi1 * (π / 2.0) * dYx)
            dvy = αv * ((1.0 - xi1) * dXy + xi1 * (π / 2.0) * dYy)

            DJ[I] = abs(dux * dvy - dvx * duy)

         end

      elseif rdisc == 4

         @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            θ = π * xi2 / 2.0 + 5.0 * π / 4.0
            s, c = sincos(θ)

            Yx = d.Rc * c
            Yy = d.Rc * s

            Xx = -d.Lc * (1.0 - 2.0 * xi2) / 2.0
            Xy = -d.Lc / 2.0

            dYx = -d.Rc * s
            dYy =  d.Rc * c

            dXx = d.Lc
            dXy = 0.0

            Zx[I] = d.A + (1.0 - xi1) * Xx + xi1 * Yx
            Zy[I] = d.B + (1.0 - xi1) * Xy + xi1 * Yy

            dux = αu * (Yx - Xx)
            duy = αu * (Yy - Xy)

            dvx = αv * ((1.0 - xi1) * dXx + xi1 * (π / 2.0) * dYx)
            dvy = αv * (xi1 * (π / 2.0) * dYy)

            DJ[I] = abs(dux * dvy - dvx * duy)

         end

      else

         throw(ArgumentError("mapxy_Dmap!(star): invalid central region rdisc=$rdisc"))

      end

   end

   return nothing

end

function chk_map(d::star; n::Int=32, tol::Float64=5e-14)
   c1 = 3 * d.Ps^2 + 8 * d.Ps + 8
   c2 = 35 * d.Ps^4 + 160 * d.Ps^3 + 288 * d.Ps^2 + 256 * d.Ps + 128
   Iex = π * d.R^2 * ((d.A^2 + d.B^2) * c1 + d.R^2 * c2 / 32)/8
   return chkmap_geom(d, Iex; n=n, tol=tol)
end

"""
   jinvmap(d::star, u::Float64, v::Float64, r::Int) -> Tuple{Float64,Float64,Float64,Float64}

Inverse Jacobian of the region mapping at `(u,v) ∈ [-1,1]^2`.

Here `r` is the actual region number, not a patch index.

Returns `(J11, J12, J21, J22)` where

   J⁻¹ = [J11 J12;J21 J22]

For this region-level inverse, `hc = ht = 1`.
"""
function jinvmap(d::star, u::Float64, v::Float64, r::Int)

   # Region-level affine coordinates, since hc = ht = 1.
   xi1 = 0.5 * (u + 1.0)
   xi2 = 0.5 * (v + 1.0)

   αu = 0.5
   αv = 0.5

   if r <= 7 * d.P

      jpet = cld(r, 7)
      regloc = mod1(r, 7)

      e1x = d.RP.e₁x[jpet]
      e1y = d.RP.e₁y[jpet]
      e2x = d.RP.e₂x[jpet]
      e2y = d.RP.e₂y[jpet]
      qx = d.RP.qx[jpet]
      qy = d.RP.qy[jpet]
      nme = d.RP.nme
      invn = 1.0 / nme

      if regloc == 1

         ak = 2.0 * d.tht2
         bk = 2.0 * (jpet - 1) * π / d.P - d.tht2

         θ = muladd(ak, xi2, bk)
         s, c = sincos(θ)
         sP, cP = sincos(d.P * θ)

         h = 1.0 + d.Ps / 2.0 + d.Ps * cP / 2.0
         h1 = -0.5 * d.Ps * d.P * sP

         Yx = d.R * h * c
         Yy = d.R * h * s

         dYx = d.R * (h1 * c - h * s)
         dYy = d.R * (h1 * s + h * c)

         Xx = d.R * e2x +
              (2.0 * xi2 - 1.0) * d.L1 * e1x * invn / 2.0 +
              d.L2 * qx * invn

         Xy = d.R * e2y +
              (2.0 * xi2 - 1.0) * d.L1 * e1y * invn / 2.0 +
              d.L2 * qy * invn

         dXx = d.L1 * e1x * invn
         dXy = d.L1 * e1y * invn

         dux = αu * (Yx - Xx)
         duy = αu * (Yy - Xy)

         dvx = αv * ((1.0 - xi1) * dXx + xi1 * ak * dYx)
         dvy = αv * ((1.0 - xi1) * dXy + xi1 * ak * dYy)

      elseif regloc == 2

         ak = π / d.P - d.tht2 - d.tht1
         bk = 2.0 * (jpet - 1) * π / d.P + d.tht2

         θ = muladd(ak, xi2, bk)
         s, c = sincos(θ)
         sP, cP = sincos(d.P * θ)

         h = 1.0 + d.Ps / 2.0 + d.Ps * cP / 2.0
         h1 = -0.5 * d.Ps * d.P * sP

         Yx = d.R * h * c
         Yy = d.R * h * s

         dYx = d.R * (h1 * c - h * s)
         dYy = d.R * (h1 * s + h * c)

         Xx = d.R * e2x +
              d.L1 * e1x * invn / 2.0 +
              (1.0 - xi2) * d.L2 * qx * invn

         Xy = d.R * e2y +
              d.L1 * e1y * invn / 2.0 +
              (1.0 - xi2) * d.L2 * qy * invn

         dXx = -d.L2 * qx * invn
         dXy = -d.L2 * qy * invn

         dux = αu * (Yx - Xx)
         duy = αu * (Yy - Xy)

         dvx = αv * ((1.0 - xi1) * dXx + xi1 * ak * dYx)
         dvy = αv * ((1.0 - xi1) * dXy + xi1 * ak * dYy)

      elseif regloc == 3

         ak = d.tht1
         bk = (2.0 * jpet - 1.0) * π / d.P - d.tht1

         alpx = d.RP.αx[jpet+1]
         alpy = d.RP.αy[jpet+1]

         Cx = d.R * e2x + d.L1 * e1x * invn / 2.0
         Cy = d.R * e2y + d.L1 * e1y * invn / 2.0

         θ = muladd(ak, xi2, bk)
         s, c = sincos(θ)
         sP, cP = sincos(d.P * θ)

         h = 1.0 + d.Ps / 2.0 + d.Ps * cP / 2.0
         h1 = -0.5 * d.Ps * d.P * sP

         Yx = d.R * h * c
         Yy = d.R * h * s

         dYx = d.R * (h1 * c - h * s)
         dYy = d.R * (h1 * s + h * c)

         Xx = xi2 * alpx + (1.0 - xi2) * Cx
         Xy = xi2 * alpy + (1.0 - xi2) * Cy

         dXx = alpx - Cx
         dXy = alpy - Cy

         dux = αu * (Yx - Xx)
         duy = αu * (Yy - Xy)

         dvx = αv * ((1.0 - xi1) * dXx + xi1 * ak * dYx)
         dvy = αv * ((1.0 - xi1) * dXy + xi1 * ak * dYy)

      elseif regloc == 4

         ak = 2.0 * π / d.P
         bk = (2.0 * jpet - 3.0) * π / d.P

         θ = muladd(ak, xi2, bk)
         s, c = sincos(θ)

         Yx = d.Rc * c
         Yy = d.Rc * s

         dYx = -d.Rc * s
         dYy = d.Rc * c

         Xx = d.R * e2x +
              (2.0 * xi2 - 1.0) * d.L1 * e1x * invn / 2.0

         Xy = d.R * e2y +
              (2.0 * xi2 - 1.0) * d.L1 * e1y * invn / 2.0

         dXx = d.L1 * e1x * invn
         dXy = d.L1 * e1y * invn

         dux = αu * (Yx - Xx)
         duy = αu * (Yy - Xy)

         dvx = αv * ((1.0 - xi1) * dXx + xi1 * ak * dYx)
         dvy = αv * ((1.0 - xi1) * dXy + xi1 * ak * dYy)

      elseif regloc == 5

         ak = d.tht1
         bk = (2.0 * jpet - 3.0) * π / d.P

         alpx = d.RP.αx[jpet]
         alpy = d.RP.αy[jpet]

         Cx = d.R * e2x - d.L1 * e1x * invn / 2.0
         Cy = d.R * e2y - d.L1 * e1y * invn / 2.0

         θ = muladd(ak, xi2, bk)
         s, c = sincos(θ)
         sP, cP = sincos(d.P * θ)

         h = 1.0 + d.Ps / 2.0 + d.Ps * cP / 2.0
         h1 = -0.5 * d.Ps * d.P * sP

         Yx = d.R * h * c
         Yy = d.R * h * s

         dYx = d.R * (h1 * c - h * s)
         dYy = d.R * (h1 * s + h * c)

         Xx = (1.0 - xi2) * alpx + xi2 * Cx
         Xy = (1.0 - xi2) * alpy + xi2 * Cy

         dXx = Cx - alpx
         dXy = Cy - alpy

         dux = αu * (Yx - Xx)
         duy = αu * (Yy - Xy)

         dvx = αv * ((1.0 - xi1) * dXx + xi1 * ak * dYx)
         dvy = αv * ((1.0 - xi1) * dXy + xi1 * ak * dYy)

      elseif regloc == 6

         ak = π / d.P - d.tht2 - d.tht1
         bk = (2.0 * jpet - 3.0) * π / d.P + d.tht1

         θ = muladd(ak, xi2, bk)
         s, c = sincos(θ)
         sP, cP = sincos(d.P * θ)

         h = 1.0 + d.Ps / 2.0 + d.Ps * cP / 2.0
         h1 = -0.5 * d.Ps * d.P * sP

         Yx = d.R * h * c
         Yy = d.R * h * s

         dYx = d.R * (h1 * c - h * s)
         dYy = d.R * (h1 * s + h * c)

         Xx = d.R * e2x -
              d.L1 * e1x * invn / 2.0 +
              xi2 * d.L2 * qx * invn

         Xy = d.R * e2y -
              d.L1 * e1y * invn / 2.0 +
              xi2 * d.L2 * qy * invn

         dXx = d.L2 * qx * invn
         dXy = d.L2 * qy * invn

         dux = αu * (Yx - Xx)
         duy = αu * (Yy - Xy)

         dvx = αv * ((1.0 - xi1) * dXx + xi1 * ak * dYx)
         dvy = αv * ((1.0 - xi1) * dXy + xi1 * ak * dYy)

      else # regloc == 7, petal rectangle

         # τ_u = (L2/2nme) q, τ_v = (L1/2nme) e1.
         dux = d.L2 * qx * invn / 2.0
         duy = d.L2 * qy * invn / 2.0

         dvx = d.L1 * e1x * invn / 2.0
         dvy = d.L1 * e1y * invn / 2.0

      end

   else

      rdisc = r - 7 * d.P

      if rdisc == 1

         θ = π * xi2 / 2.0 - π / 4.0
         s, c = sincos(θ)

         Yx = d.Rc * c
         Yy = d.Rc * s

         Xx = d.Lc / 2.0
         Xy = d.Lc * (2.0 * xi2 - 1.0) / 2.0

         dYx = -d.Rc * s
         dYy = d.Rc * c

         dXx = 0.0
         dXy = d.Lc

         dux = αu * (Yx - Xx)
         duy = αu * (Yy - Xy)

         dvx = αv * (xi1 * (π / 2.0) * dYx)
         dvy = αv * ((1.0 - xi1) * dXy + xi1 * (π / 2.0) * dYy)

      elseif rdisc == 2

         θ = π * xi2 / 2.0 + π / 4.0
         s, c = sincos(θ)

         Yx = d.Rc * c
         Yy = d.Rc * s

         Xx = d.Lc * (1.0 - 2.0 * xi2) / 2.0
         Xy = d.Lc / 2.0

         dYx = -d.Rc * s
         dYy = d.Rc * c

         dXx = -d.Lc
         dXy = 0.0

         dux = αu * (Yx - Xx)
         duy = αu * (Yy - Xy)

         dvx = αv * ((1.0 - xi1) * dXx + xi1 * (π / 2.0) * dYx)
         dvy = αv * (xi1 * (π / 2.0) * dYy)

      elseif rdisc == 3

         θ = π * xi2 / 2.0 + 3.0 * π / 4.0
         s, c = sincos(θ)

         Yx = d.Rc * c
         Yy = d.Rc * s

         Xx = -d.Lc / 2.0
         Xy = -d.Lc * (2.0 * xi2 - 1.0) / 2.0

         dYx = -d.Rc * s
         dYy = d.Rc * c

         dXx = 0.0
         dXy = -d.Lc

         dux = αu * (Yx - Xx)
         duy = αu * (Yy - Xy)

         dvx = αv * (xi1 * (π / 2.0) * dYx)
         dvy = αv * ((1.0 - xi1) * dXy + xi1 * (π / 2.0) * dYy)

      elseif rdisc == 4

         θ = π * xi2 / 2.0 + 5.0 * π / 4.0
         s, c = sincos(θ)

         Yx = d.Rc * c
         Yy = d.Rc * s

         Xx = -d.Lc * (1.0 - 2.0 * xi2) / 2.0
         Xy = -d.Lc / 2.0

         dYx = -d.Rc * s
         dYy = d.Rc * c

         dXx = d.Lc
         dXy = 0.0

         dux = αu * (Yx - Xx)
         duy = αu * (Yy - Xy)

         dvx = αv * ((1.0 - xi1) * dXx + xi1 * (π / 2.0) * dYx)
         dvy = αv * (xi1 * (π / 2.0) * dYy)

      else # rdisc == 5, central square

         dux = d.Lc / 2.0
         duy = 0.0
         dvx = 0.0
         dvy = d.Lc / 2.0

      end

   end

   detJ = dux * dvy - dvx * duy
   invdetJ = 1.0 / detJ

   return dvy * invdetJ, -dvx * invdetJ, -duy * invdetJ, dux * invdetJ

end

"""
  mapinv(d::star, u::Float64, v::Float64, k::Int) -> Tuple{Float64,Float64}

Inverse mapping: given a physical point `(u,v)` on patch `k`, return
the reference coordinates `[t; s]` in `[-1,1]^2` such that `mapxy(d, t, s, k) = (u,v)`.
"""
@inline function Xx(s::Float64, d::star, r::Int)

   if r <= 7 * d.P

      jpet = cld(r, 7)
      regloc = mod1(r, 7)

      e1x = d.RP.e₁x[jpet]
      e2x = d.RP.e₂x[jpet]
      qx = d.RP.qx[jpet]
      nme = d.RP.nme
      invn = 1.0 / nme

      if regloc == 1

         return d.A + d.R * e2x +
                (2.0 * s - 1.0) * d.L1 * e1x * invn / 2.0 +
                d.L2 * qx * invn

      elseif regloc == 2

         return d.A + d.R * e2x +
                d.L1 * e1x * invn / 2.0 +
                (1.0 - s) * d.L2 * qx * invn

      elseif regloc == 3

         alpx = d.RP.αx[jpet+1]

         return d.A + s * alpx +
                (1.0 - s) * (d.R * e2x + d.L1 * e1x * invn / 2.0)

      elseif regloc == 4

         return d.A + d.R * e2x +
                (2.0 * s - 1.0) * d.L1 * e1x * invn / 2.0

      elseif regloc == 5

         alpx = d.RP.αx[jpet]

         return d.A + (1.0 - s) * alpx +
                s * (d.R * e2x - d.L1 * e1x * invn / 2.0)

      elseif regloc == 6

         return d.A + d.R * e2x -
                d.L1 * e1x * invn / 2.0 +
                s * d.L2 * qx * invn

      else

         throw(ArgumentError("Xx(star) expects non-rectangular petal region; got r=$r, regloc=$regloc"))

      end

   else

      rdisc = r - 7 * d.P

      if rdisc == 1

         return d.A + d.Lc / 2.0

      elseif rdisc == 2

         return d.A + d.Lc * (1.0 - 2.0 * s) / 2.0

      elseif rdisc == 3

         return d.A - d.Lc / 2.0

      elseif rdisc == 4

         return d.A - d.Lc * (1.0 - 2.0 * s) / 2.0

      else

         throw(ArgumentError("Xx(star) expects central region 1–4; got rdisc=$rdisc"))

      end

   end

end

@inline function Xy(s::Float64, d::star, r::Int)

   if r <= 7 * d.P

      jpet = cld(r, 7)
      regloc = mod1(r, 7)

      e1y = d.RP.e₁y[jpet]
      e2y = d.RP.e₂y[jpet]
      qy = d.RP.qy[jpet]
      nme = d.RP.nme
      invn = 1.0 / nme

      if regloc == 1

         return d.B + d.R * e2y +
                (2.0 * s - 1.0) * d.L1 * e1y * invn / 2.0 +
                d.L2 * qy * invn

      elseif regloc == 2

         return d.B + d.R * e2y +
                d.L1 * e1y * invn / 2.0 +
                (1.0 - s) * d.L2 * qy * invn

      elseif regloc == 3

         alpy = d.RP.αy[jpet+1]

         return d.B + s * alpy +
                (1.0 - s) * (d.R * e2y + d.L1 * e1y * invn / 2.0)

      elseif regloc == 4

         return d.B + d.R * e2y +
                (2.0 * s - 1.0) * d.L1 * e1y * invn / 2.0

      elseif regloc == 5

         alpy = d.RP.αy[jpet]

         return d.B + (1.0 - s) * alpy +
                s * (d.R * e2y - d.L1 * e1y * invn / 2.0)

      elseif regloc == 6

         return d.B + d.R * e2y -
                d.L1 * e1y * invn / 2.0 +
                s * d.L2 * qy * invn

      else

         throw(ArgumentError("Xy(star) expects non-rectangular petal region; got r=$r, regloc=$regloc"))

      end

   else

      rdisc = r - 7 * d.P

      if rdisc == 1

         return d.B + d.Lc * (2.0 * s - 1.0) / 2.0

      elseif rdisc == 2

         return d.B + d.Lc / 2.0

      elseif rdisc == 3

         return d.B - d.Lc * (2.0 * s - 1.0) / 2.0

      elseif rdisc == 4

         return d.B - d.Lc / 2.0

      else

         throw(ArgumentError("Xy(star) expects central region 1–4; got rdisc=$rdisc"))

      end

   end

end

@inline function Yx(s::Float64, d::star, r::Int)

   if r <= 7 * d.P

      jpet = cld(r, 7)
      regloc = mod1(r, 7)

      if regloc == 1

         θ = s * (2.0 * d.tht2) +
             2.0 * (jpet - 1) * π / d.P - d.tht2

         st, ct = sincos(θ)
         h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

         return d.A + d.R * h * ct

      elseif regloc == 2

         θ = s * (π / d.P - d.tht2 - d.tht1) +
             2.0 * (jpet - 1) * π / d.P + d.tht2

         st, ct = sincos(θ)
         h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

         return d.A + d.R * h * ct

      elseif regloc == 3

         θ = s * d.tht1 +
             (2.0 * jpet - 1.0) * π / d.P - d.tht1

         st, ct = sincos(θ)
         h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

         return d.A + d.R * h * ct

      elseif regloc == 4

         θ = s * (2.0 * π / d.P) +
             (2.0 * jpet - 3.0) * π / d.P

         return d.A + d.Rc * cos(θ)

      elseif regloc == 5

         θ = s * d.tht1 +
             (2.0 * jpet - 3.0) * π / d.P

         st, ct = sincos(θ)
         h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

         return d.A + d.R * h * ct

      elseif regloc == 6

         θ = s * (π / d.P - d.tht2 - d.tht1) +
             (2.0 * jpet - 3.0) * π / d.P + d.tht1

         st, ct = sincos(θ)
         h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

         return d.A + d.R * h * ct

      else

         throw(ArgumentError("Yx(star) expects non-rectangular petal region; got r=$r, regloc=$regloc"))

      end

   else

      rdisc = r - 7 * d.P

      if 1 <= rdisc <= 4

         return d.A + d.Rc * cospi(s / 2.0 + (2.0 * rdisc - 3.0) / 4.0)

      else

         throw(ArgumentError("Yx(star) expects central region 1–4; got rdisc=$rdisc"))

      end

   end

end

@inline function Yy(s::Float64, d::star, r::Int)

   if r <= 7 * d.P

      jpet = cld(r, 7)
      regloc = mod1(r, 7)

      if regloc == 1

         θ = s * (2.0 * d.tht2) +
             2.0 * (jpet - 1) * π / d.P - d.tht2

         st, ct = sincos(θ)
         h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

         return d.B + d.R * h * st

      elseif regloc == 2

         θ = s * (π / d.P - d.tht2 - d.tht1) +
             2.0 * (jpet - 1) * π / d.P + d.tht2

         st, ct = sincos(θ)
         h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

         return d.B + d.R * h * st

      elseif regloc == 3

         θ = s * d.tht1 +
             (2.0 * jpet - 1.0) * π / d.P - d.tht1

         st, ct = sincos(θ)
         h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

         return d.B + d.R * h * st

      elseif regloc == 4

         θ = s * (2.0 * π / d.P) +
             (2.0 * jpet - 3.0) * π / d.P

         return d.B + d.Rc * sin(θ)

      elseif regloc == 5

         θ = s * d.tht1 +
             (2.0 * jpet - 3.0) * π / d.P

         st, ct = sincos(θ)
         h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

         return d.B + d.R * h * st

      elseif regloc == 6

         θ = s * (π / d.P - d.tht2 - d.tht1) +
             (2.0 * jpet - 3.0) * π / d.P + d.tht1

         st, ct = sincos(θ)
         h = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

         return d.B + d.R * h * st

      else

         throw(ArgumentError("Yy(star) expects non-rectangular petal region; got r=$r, regloc=$regloc"))

      end

   else

      rdisc = r - 7 * d.P

      if 1 <= rdisc <= 4

         return d.B + d.Rc * sinpi(s / 2.0 + (2.0 * rdisc - 3.0) / 4.0)

      else

         throw(ArgumentError("Yy(star) expects central region 1–4; got rdisc=$rdisc"))

      end

   end

end

function fill_FTable!(tbl::FTable, d::star, r::Int)

   vmin, vmax, N = tbl.vmin, tbl.vmax, tbl.N
   htab = (vmax - vmin) / (N - 1)

   P1, P2, P3 = tbl.P1, tbl.P2, tbl.P3

   if r <= 7 * d.P

      jpet = cld(r, 7)
      regloc = mod1(r, 7)

      e1x = d.RP.e₁x[jpet]
      e1y = d.RP.e₁y[jpet]
      e2x = d.RP.e₂x[jpet]
      e2y = d.RP.e₂y[jpet]
      qx = d.RP.qx[jpet]
      qy = d.RP.qy[jpet]
      nme = d.RP.nme
      invn = 1.0 / nme

      if regloc == 1

         ak = 2.0 * d.tht2
         bk = 2.0 * (jpet - 1) * π / d.P - d.tht2

         @inbounds for i in 1:N
            s = vmin + (i - 1) * htab

            θ = muladd(ak, s, bk)
            st, ct = sincos(θ)
            hρ = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

            Xxv = d.A + d.R * e2x +
                  (2.0 * s - 1.0) * d.L1 * e1x * invn / 2.0 +
                  d.L2 * qx * invn

            Xyv = d.B + d.R * e2y +
                  (2.0 * s - 1.0) * d.L1 * e1y * invn / 2.0 +
                  d.L2 * qy * invn

            Yxv = d.A + d.R * hρ * ct
            Yyv = d.B + d.R * hρ * st

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv
         end

      elseif regloc == 2

         ak = π / d.P - d.tht2 - d.tht1
         bk = 2.0 * (jpet - 1) * π / d.P + d.tht2

         @inbounds for i in 1:N
            s = vmin + (i - 1) * htab

            θ = muladd(ak, s, bk)
            st, ct = sincos(θ)
            hρ = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

            Xxv = d.A + d.R * e2x +
                  d.L1 * e1x * invn / 2.0 +
                  (1.0 - s) * d.L2 * qx * invn

            Xyv = d.B + d.R * e2y +
                  d.L1 * e1y * invn / 2.0 +
                  (1.0 - s) * d.L2 * qy * invn

            Yxv = d.A + d.R * hρ * ct
            Yyv = d.B + d.R * hρ * st

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv
         end

      elseif regloc == 3

         ak = d.tht1
         bk = (2.0 * jpet - 1.0) * π / d.P - d.tht1

         alpx = d.RP.αx[jpet+1]
         alpy = d.RP.αy[jpet+1]

         Cx = d.A + d.R * e2x + d.L1 * e1x * invn / 2.0
         Cy = d.B + d.R * e2y + d.L1 * e1y * invn / 2.0

         @inbounds for i in 1:N
            s = vmin + (i - 1) * htab

            θ = muladd(ak, s, bk)
            st, ct = sincos(θ)
            hρ = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

            Xxv = d.A + s * alpx + (1.0 - s) * (Cx - d.A)
            Xyv = d.B + s * alpy + (1.0 - s) * (Cy - d.B)

            Yxv = d.A + d.R * hρ * ct
            Yyv = d.B + d.R * hρ * st

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv
         end

      elseif regloc == 4

         ak = 2.0 * π / d.P
         bk = (2.0 * jpet - 3.0) * π / d.P

         @inbounds for i in 1:N
            s = vmin + (i - 1) * htab

            θ = muladd(ak, s, bk)
            st, ct = sincos(θ)

            Xxv = d.A + d.R * e2x +
                  (2.0 * s - 1.0) * d.L1 * e1x * invn / 2.0

            Xyv = d.B + d.R * e2y +
                  (2.0 * s - 1.0) * d.L1 * e1y * invn / 2.0

            Yxv = d.A + d.Rc * ct
            Yyv = d.B + d.Rc * st

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv
         end

      elseif regloc == 5

         ak = d.tht1
         bk = (2.0 * jpet - 3.0) * π / d.P

         alpx = d.RP.αx[jpet]
         alpy = d.RP.αy[jpet]

         Cx = d.A + d.R * e2x - d.L1 * e1x * invn / 2.0
         Cy = d.B + d.R * e2y - d.L1 * e1y * invn / 2.0

         @inbounds for i in 1:N
            s = vmin + (i - 1) * htab

            θ = muladd(ak, s, bk)
            st, ct = sincos(θ)
            hρ = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

            Xxv = d.A + (1.0 - s) * alpx + s * (Cx - d.A)
            Xyv = d.B + (1.0 - s) * alpy + s * (Cy - d.B)

            Yxv = d.A + d.R * hρ * ct
            Yyv = d.B + d.R * hρ * st

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv
         end

      elseif regloc == 6

         ak = π / d.P - d.tht2 - d.tht1
         bk = (2.0 * jpet - 3.0) * π / d.P + d.tht1

         @inbounds for i in 1:N
            s = vmin + (i - 1) * htab

            θ = muladd(ak, s, bk)
            st, ct = sincos(θ)
            hρ = 1.0 + d.Ps / 2.0 + d.Ps * cos(d.P * θ) / 2.0

            Xxv = d.A + d.R * e2x -
                  d.L1 * e1x * invn / 2.0 +
                  s * d.L2 * qx * invn

            Xyv = d.B + d.R * e2y -
                  d.L1 * e1y * invn / 2.0 +
                  s * d.L2 * qy * invn

            Yxv = d.A + d.R * hρ * ct
            Yyv = d.B + d.R * hρ * st

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv
         end

      else

         throw(ArgumentError("fill_FTable!(star): expected non-rectangular petal region; got r=$r, regloc=$regloc"))

      end

   else

      rdisc = r - 7 * d.P

      if rdisc == 1

         @inbounds for i in 1:N
            s = vmin + (i - 1) * htab

            thp = s / 2.0 - 1.0 / 4.0

            Xxv = d.A + d.Lc / 2.0
            Xyv = d.B + d.Lc * (2.0 * s - 1.0) / 2.0

            Yxv = d.A + d.Rc * cospi(thp)
            Yyv = d.B + d.Rc * sinpi(thp)

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv
         end

      elseif rdisc == 2

         @inbounds for i in 1:N
            s = vmin + (i - 1) * htab

            thp = s / 2.0 + 1.0 / 4.0

            Xxv = d.A + d.Lc * (1.0 - 2.0 * s) / 2.0
            Xyv = d.B + d.Lc / 2.0

            Yxv = d.A + d.Rc * cospi(thp)
            Yyv = d.B + d.Rc * sinpi(thp)

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv
         end

      elseif rdisc == 3

         @inbounds for i in 1:N
            s = vmin + (i - 1) * htab

            thp = s / 2.0 + 3.0 / 4.0

            Xxv = d.A - d.Lc / 2.0
            Xyv = d.B - d.Lc * (2.0 * s - 1.0) / 2.0

            Yxv = d.A + d.Rc * cospi(thp)
            Yyv = d.B + d.Rc * sinpi(thp)

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv
         end

      elseif rdisc == 4

         @inbounds for i in 1:N
            s = vmin + (i - 1) * htab

            thp = s / 2.0 + 5.0 / 4.0

            Xxv = d.A - d.Lc * (1.0 - 2.0 * s) / 2.0
            Xyv = d.B - d.Lc / 2.0

            Yxv = d.A + d.Rc * cospi(thp)
            Yyv = d.B + d.Rc * sinpi(thp)

            P1[i] = Yyv - Xyv
            P2[i] = Yxv - Xxv
            P3[i] = Xyv * Yxv - Xxv * Yyv
         end

      else

         throw(ArgumentError("fill_FTable!(star): expected central region 1–4; got rdisc=$rdisc"))

      end

   end

   tbl.reg = r

   return tbl

end

@inline function f1I(t::Float64, s::Float64,
   d::star, u::Float64, v::Float64, r::Int)

   ŝ = (s + 1.0) / 2.0
   t̂ = (t + 1.0) / 2.0

   Xxv = Xx(ŝ, d, r)
   Yxv = Yx(ŝ, d, r)

   return (1.0 - t̂) * Xxv + t̂ * Yxv - u

end

@inline function f2I(t::Float64, s::Float64,
   d::star, u::Float64, v::Float64, r::Int)

   ŝ = (s + 1.0) / 2.0
   t̂ = (t + 1.0) / 2.0

   Xyv = Xy(ŝ, d, r)
   Yyv = Yy(ŝ, d, r)

   return (1.0 - t̂) * Xyv + t̂ * Yyv - v

end

@inline function JinvI(t::Float64, s::Float64,
   d::star, u::Float64, v::Float64, r::Int)

   return jinvmap(d, t, s, r)

end

@inline function f_cont(v̂::Float64, d::star,
   u::Float64, v::Float64, r::Int)

   Xxv = Xx(v̂, d, r)
   Xyv = Xy(v̂, d, r)

   Yxv = Yx(v̂, d, r)
   Yyv = Yy(v̂, d, r)

   dx = Yxv - Xxv
   dy = Yyv - Xyv

   term1 = muladd(u, dy, -v * dx)          # u*dy - v*dx
   term2 = muladd(Xyv, Yxv, -Xxv * Yyv)    # Xy*Yx - Xx*Yy

   return term1 + term2

end

function mapinv(tbl::FTable, d::star, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]
   r = p.reg

   # ----- Stage 1: 2D Newton -----
   tN, sN = newtonR2D(f1I, f2I, JinvI,
      0.0, 0.0, 4, d, u, v, r; tol=1e-15)

   if tN !== :max

      return xi_inv((1.0 + tN) / 2.0, p.ck0, p.ck1),
      xi_inv((1.0 + sN) / 2.0, p.tk0, p.tk1)

   end

   # ----- Stage 2: fallback 1D root search -----
   rmi = tbl.rmi
   zxi = tbl.zxi

   n = find_roots!(rmi, tbl, d, u, v, r)

   @inbounds for i in 1:n

      v̂ = rmi[i]

      Xxv = Xx(v̂, d, r)
      Xyv = Xy(v̂, d, r)

      Yxv = Yx(v̂, d, r)
      Yyv = Yy(v̂, d, r)

      if abs(v - Xyv) < abs(u - Xxv)
         zxi[i] = (u - Xxv) / (Yxv - Xxv)
      else
         zxi[i] = (v - Xyv) / (Yyv - Xyv)
      end

   end

   # Choose best candidate by normalized patch cost.
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

function mapinv(d::star, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]
   r = p.reg

   if r <= 7 * d.P

      jpet = cld(r, 7)
      regloc = mod1(r, 7)

      if regloc == 7

         e1x = d.RP.e₁x[jpet]
         e1y = d.RP.e₁y[jpet]

         e2x = d.RP.e₂x[jpet]
         e2y = d.RP.e₂y[jpet]

         qx = d.RP.qx[jpet]
         qy = d.RP.qy[jpet]

         nme = d.RP.nme
         invn = 1.0 / nme

         wx = u - d.A - d.R * e2x
         wy = v - d.B - d.R * e2y

         # Coordinates in the affine petal rectangle.
         ξq = (wx * qx + wy * qy) * (invn / d.L2)

         ξe = muladd(wx * e1x + wy * e1y, invn / d.L1, 0.5)

         if d.L2 >= d.L1

            # mapxy swaps axes for this rectangular patch.
            Z1 = xi_inv(ξq, p.tk0, p.tk1)
            Z2 = xi_inv(ξe, p.ck0, p.ck1)

         else

            Z1 = xi_inv(ξq, p.ck0, p.ck1)
            Z2 = xi_inv(ξe, p.tk0, p.tk1)

         end

         return Z1, Z2

      else

         error("mapinv(d::star, u, v, k): non-rectangular petal region requires mapinv(tbl, d, u, v, k). Got reg=$r, regloc=$regloc.")

      end

   else

      rdisc = r - 7 * d.P

      if rdisc == 5

         ξ1 = (u - d.A + d.Lc / 2.0) / d.Lc
         ξ2 = (v - d.B + d.Lc / 2.0) / d.Lc

         Z1 = xi_inv(ξ1, p.ck0, p.ck1)
         Z2 = xi_inv(ξ2, p.tk0, p.tk1)

         return Z1, Z2

      else

         error("mapinv(d::star, u, v, k): non-rectangular central region requires mapinv(tbl, d, u, v, k). Got reg=$r, rdisc=$rdisc.")

      end

   end

end

"""
   ptconv(d::star, t1::Float64, t2::Float64, idx::Int, ptdest::String)

Convert a parametric point between region ↔ patch coordinates.

Inputs:
- `t1` : first coordinate of point t
- `t2` : second coordinate of point t
- `idx`: region index if `"to_pth"`, or patch index if `"to_reg"`
- `ptdest`: `"to_pth"` or `"to_reg"`

Returns:
- `Tuple(Float64, Float64, out_idx::Int)`
where `out_idx` is patch index if `"to_pth"`, region index if `"to_reg"`.
"""
function ptconv(d::star, t1::Float64, t2::Float64, idx::Int, ptdest::String)

   if ptdest == "to_pth"

      r = idx

      for k in 1:d.Npat

         p = d.pths[k]
         p.reg == r || continue

         if r <= 7 * d.P && mod1(r, 7) == 7

            if d.L1 > d.L2
               t1k = xi_inv((t1 + 1.0) / 2.0, p.ck0, p.ck1)
               t2k = xi_inv((t2 + 1.0) / 2.0, p.tk0, p.tk1)
            else
               t1k = xi_inv((t1 + 1.0) / 2.0, p.tk0, p.tk1)
               t2k = xi_inv((t2 + 1.0) / 2.0, p.ck0, p.ck1)
            end

         else

            t1k = xi_inv((t1 + 1.0) / 2.0, p.ck0, p.ck1)
            t2k = xi_inv((t2 + 1.0) / 2.0, p.tk0, p.tk1)

         end

         if abs(t1k) ≤ 1.0 && abs(t2k) ≤ 1.0
            return t1k, t2k, k
         end

      end

      error("ptconv(to_pth): no patch in region $r contained the point.")

   elseif ptdest == "to_reg"

      k = idx
      p = d.pths[k]
      r = p.reg

      if r <= 7 * d.P && mod1(r, 7) == 7

         if d.L1 > d.L2
            t1r = (p.ck1 - p.ck0) * t1 + (p.ck1 + p.ck0) - 1.0
            t2r = (p.tk1 - p.tk0) * t2 + (p.tk1 + p.tk0) - 1.0
         else
            t1r = (p.tk1 - p.tk0) * t1 + (p.tk1 + p.tk0) - 1.0
            t2r = (p.ck1 - p.ck0) * t2 + (p.ck1 + p.ck0) - 1.0
         end

      else

         t1r = (p.ck1 - p.ck0) * t1 + (p.ck1 + p.ck0) - 1.0
         t2r = (p.tk1 - p.tk0) * t2 + (p.tk1 + p.tk0) - 1.0

      end

      return t1r, t2r, r

   else

      error("ptconv: ptdest must be \"to_pth\" or \"to_reg\"")

   end

end

"""
   mapinv2(d::star, t1::Float64, t2::Float64, k2::Int, k::Int) -> Tuple{Float64,Float64}

Given a point `(u,v) = τₖ₂(t1,t2)` on patch `k2`, return its reference
coordinates on patch `k` in `[-1,1]^2`.

This differs from `mapinv` because we already know `(u,v)` comes from
`(t1,t2)` on `k2`.

This should only be called when patches `k2` and `k` are in the same region.
"""
function mapinv2(d::star, t1::Float64, t2::Float64, k2::Int, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]

   tr1, tr2, _ = ptconv(d, t1, t2, k2, "to_reg")

   û = (tr1 + 1.0) / 2.0
   v̂ = (tr2 + 1.0) / 2.0

   if p.reg <= 7 * d.P && mod1(p.reg, 7) == 7 && d.L2 >= d.L1

      t1k = xi_inv(û, p.tk0, p.tk1)
      t2k = xi_inv(v̂, p.ck0, p.ck1)

   else

      t1k = xi_inv(û, p.ck0, p.ck1)
      t2k = xi_inv(v̂, p.tk0, p.tk1)

   end

   return t1k, t2k

end

"""
   dfunc(d::star, k::Int, t::Float64, s::Float64)::Float64

Factor that goes linearly to zero as we approach the outer star boundary.

- If patch `k` is rectangular or not an outer-boundary petal region, returns `1`.
- Otherwise uses

   val = 1 - ck1 + (ck1 - ck0) * t / 2

  and raises to `s-1` if `s ≥ 0.5`, else to `s`.

Notes:
- `t` is expected to be `1 - t_actual`.
"""
function dfunc(d::star, k::Int, t::Float64, s::Float64)::Float64

   p = d.pths[k]
   r = p.reg

   # Only petal regions touching the outer star boundary get the boundary factor.
   # Exclude:
   #   - central regions: r > 7P
   #   - petal rectangle: regloc == 7
   #   - inner circular petal region: regloc == 4
   if r > 7 * d.P || mod1(r, 7) == 7 || mod1(r, 7) == 4
      return 1.0
   end

   hc = p.ck1 - p.ck0
   val = (1.0 - p.ck1) + hc * t / 2.0
   exp = s ≥ 0.5 ? (s - 1.0) : s

   return val^exp

end

function dfunc!(out::StridedArray{Float64}, d::star, k::Int,
   t::StridedArray{Float64}, s::Float64)

   p = d.pths[k]
   r = p.reg

   # Only petal regions touching the outer star boundary get the boundary factor.
   if r > 7 * d.P || mod1(r, 7) == 7 || mod1(r, 7) == 4

      fill!(out, 1.0)

   else

      exp = s ≥ 0.5 ? (s - 1.0) : s
      αc = (p.ck1 - p.ck0) / 2.0
      βc = 1.0 - p.ck1

      @inbounds for i in eachindex(t)
         out[i] = muladd(αc, t[i], βc)^exp
      end

   end

   return nothing

end

function dfunc!(out::StridedArray{Float64}, d::star, k::Int,
   t::StridedArray{Float64})

   p = d.pths[k]
   r = p.reg

   # Only petal regions touching the outer star boundary get the boundary factor.
   if r > 7 * d.P || mod1(r, 7) == 7 || mod1(r, 7) == 4

      fill!(out, 1.0)

   else
      αc = (p.ck1 - p.ck0) / 2.0
      βc = 1.0 - p.ck1

      @inbounds for i in eachindex(t)
         out[i] = muladd(αc, t[i], βc)
      end

   end

   return nothing

end


"""
    gamderhigher(d::star, th::Float64)

Return centered star boundary values and derivatives through fourth order.

Returns gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy, d4gx, d4gy

where gx(th) = R*h(th)*cos(th), gy(th) = R*h(th)*sin(th),
and h(th) = 1 + Ps/2 + Ps*cos(P*th)/2.
"""
function gamderhigher(d::star, th::Float64)

   s, c = sincos(th)

   P = Float64(d.P)
   a = 0.5 * d.Ps

   sP, cP = sincos(P * th)

   P2 = P * P
   P3 = P2 * P
   P4 = P2 * P2

   h  = 1.0 + a + a * cP
   h1 = -a * P  * sP
   h2 = -a * P2 * cP
   h3 =  a * P3 * sP
   h4 =  a * P4 * cP

   gx = d.R * h * c
   gy = d.R * h * s

   dgx = d.R * (h1 * c - h * s)
   dgy = d.R * (h1 * s + h * c)

   d2gx = d.R * (h2 * c - 2.0 * h1 * s - h * c)
   d2gy = d.R * (h2 * s + 2.0 * h1 * c - h * s)

   d3gx = d.R * (h3 * c - 3.0 * h2 * s - 3.0 * h1 * c + h * s)
   d3gy = d.R * (h3 * s + 3.0 * h2 * c - 3.0 * h1 * s - h * c)

   d4gx = d.R * (h4 * c - 4.0 * h3 * s - 6.0 * h2 * c +
                  4.0 * h1 * s + h * c)

   d4gy = d.R * (h4 * s + 4.0 * h3 * c - 6.0 * h2 * s -
                  4.0 * h1 * c + h * s)

   return gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy, d4gx, d4gy

end

"""
    diff_map!(out::Matrix{Float64},
              Zx::Matrix{Float64}, Zy::Matrix{Float64},
              DJ::StridedArray{Float64},
              d::star, u::Float64, v::Float64,
              u2::Matrix{Float64}, v2::Matrix{Float64},
              du::AbstractVector, dv::AbstractVector, k::Int;
              tol::Float64 = 1e-4)

Fill `out` with ‖τ(u,v) - τ(u₂,v₂)‖ on patch `k`.

Here `(u,v)` are scalars, while `u2,v2` are matrix grids. The vectors
`du,dv` satisfy

    du[i] = u - u2[i,j]
    dv[j] = v - v2[i,j]

with the usual meshgrid convention.

Uses `mapxy_Dmap!` for the far evaluation and a fourth-order Taylor fixup near
`(u,v)` to avoid cancellation.
"""
function diff_map!(out::Matrix{Float64},
   Zx::Matrix{Float64}, Zy::Matrix{Float64}, DJ::StridedArray{Float64},
   d::star, u::Float64, v::Float64,
   u2::Matrix{Float64}, v2::Matrix{Float64},
   du::AbstractVector, dv::AbstractVector, k::Int;
   tol::Float64=1e-4)

   nd_u = size(out, 1)
   nd_v = size(out, 2)

   p = d.pths[k]
   reg = p.reg

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   αc = 0.5 * hc
   αt = 0.5 * ht
   βc = p.ck0 + αc
   βt = p.tk0 + αt

   mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)

   # ------------------------------------------------------------
   # Petal rectangular region: exact affine distance.
   # ------------------------------------------------------------
   if reg <= 7 * d.P && mod1(reg, 7) == 7

      jpet = cld(reg, 7)

      e1x = d.RP.e₁x[jpet]
      e1y = d.RP.e₁y[jpet]
      qx = d.RP.qx[jpet]
      qy = d.RP.qy[jpet]
      invn = 1.0 / d.RP.nme

      if d.L2 >= d.L1

         cU = d.L2 * αt * invn
         cV = d.L1 * αc * invn

      else

         cU = d.L2 * αc * invn
         cV = d.L1 * αt * invn

      end

      @inbounds for j in 1:nd_v
         dvj = dv[j]
         @inbounds for i in 1:nd_u
            dui = du[i]

            Dx = cU * dui * qx + cV * dvj * e1x
            Dy = cU * dui * qy + cV * dvj * e1y

            out[i, j] = hypot(Dx, Dy)
         end
      end

      return nothing

   end

   # ------------------------------------------------------------
   # Central square: exact affine distance.
   # ------------------------------------------------------------
   if reg == 7 * d.P + 5

      cU = d.Lc * αc
      cV = d.Lc * αt

      @inbounds for j in 1:nd_v
         dvj = dv[j]
         @inbounds for i in 1:nd_u
            dui = du[i]

            out[i, j] = hypot(cU * dui, cV * dvj)
         end
      end

      return nothing

   end

   # Scalar reference point mapped to region coordinates.
   xi1 = muladd(αc, u, βc)
   xi2 = muladd(αt, v, βt)

   if reg <= 7 * d.P

      jpet = cld(reg, 7)
      regloc = mod1(reg, 7)

      e1x = d.RP.e₁x[jpet]
      e1y = d.RP.e₁y[jpet]
      e2x = d.RP.e₂x[jpet]
      e2y = d.RP.e₂y[jpet]
      qx = d.RP.qx[jpet]
      qy = d.RP.qy[jpet]
      invn = 1.0 / d.RP.nme

      if regloc == 1

         ak = 2.0 * d.tht2
         bk = 2.0 * (jpet - 1) * π / d.P - d.tht2

         Xx0 = d.A + d.R * e2x +
               (2.0 * xi2 - 1.0) * d.L1 * e1x * invn / 2.0 +
               d.L2 * qx * invn

         Xy0 = d.B + d.R * e2y +
               (2.0 * xi2 - 1.0) * d.L1 * e1y * invn / 2.0 +
               d.L2 * qy * invn

         dXx = d.L1 * e1x * invn
         dXy = d.L1 * e1y * invn

         th = muladd(ak, xi2, bk)

         gx, gy, dgx, dgy,d2gx, d2gy,d3gx, 
         d3gy,d4gx, d4gy = gamderhigher(d, th)

         Yx0 = d.A + gx
         Yy0 = d.B + gy

      elseif regloc == 2

         ak = π / d.P - d.tht2 - d.tht1
         bk = 2.0 * (jpet - 1) * π / d.P + d.tht2

         Xx0 = d.A + d.R * e2x +
               d.L1 * e1x * invn / 2.0 +
               (1.0 - xi2) * d.L2 * qx * invn

         Xy0 = d.B + d.R * e2y +
               d.L1 * e1y * invn / 2.0 +
               (1.0 - xi2) * d.L2 * qy * invn

         dXx = -d.L2 * qx * invn
         dXy = -d.L2 * qy * invn

         th = muladd(ak, xi2, bk)

         gx, gy, dgx, dgy,d2gx, d2gy,d3gx, 
         d3gy,d4gx, d4gy = gamderhigher(d, th)

         Yx0 = d.A + gx
         Yy0 = d.B + gy

      elseif regloc == 3

         ak = d.tht1
         bk = (2.0 * jpet - 1.0) * π / d.P - d.tht1

         alpx = d.RP.αx[jpet+1]
         alpy = d.RP.αy[jpet+1]

         Cx = d.R * e2x + d.L1 * e1x * invn / 2.0
         Cy = d.R * e2y + d.L1 * e1y * invn / 2.0

         Xx0 = d.A + xi2 * alpx + (1.0 - xi2) * Cx
         Xy0 = d.B + xi2 * alpy + (1.0 - xi2) * Cy

         dXx = alpx - Cx
         dXy = alpy - Cy

         th = muladd(ak, xi2, bk)

         gx, gy, dgx, dgy,d2gx, d2gy,d3gx, 
         d3gy,d4gx, d4gy = gamderhigher(d, th)

         Yx0 = d.A + gx
         Yy0 = d.B + gy

      elseif regloc == 4

         ak = 2.0 * π / d.P
         bk = (2.0 * jpet - 3.0) * π / d.P

         Xx0 = d.A + d.R * e2x +
               (2.0 * xi2 - 1.0) * d.L1 * e1x * invn / 2.0

         Xy0 = d.B + d.R * e2y +
               (2.0 * xi2 - 1.0) * d.L1 * e1y * invn / 2.0

         dXx = d.L1 * e1x * invn
         dXy = d.L1 * e1y * invn

         th = muladd(ak, xi2, bk)
         sth, cth = sincos(th)

         Yx0 = d.A + d.Rc * cth
         Yy0 = d.B + d.Rc * sth

         dgx = -d.Rc * sth
         dgy = d.Rc * cth
         d2gx = -d.Rc * cth
         d2gy = -d.Rc * sth
         d3gx = d.Rc * sth
         d3gy = -d.Rc * cth
         d4gx = d.Rc * cth
         d4gy = d.Rc * sth

      elseif regloc == 5

         ak = d.tht1
         bk = (2.0 * jpet - 3.0) * π / d.P

         alpx = d.RP.αx[jpet]
         alpy = d.RP.αy[jpet]

         Cx = d.R * e2x - d.L1 * e1x * invn / 2.0
         Cy = d.R * e2y - d.L1 * e1y * invn / 2.0

         Xx0 = d.A + (1.0 - xi2) * alpx + xi2 * Cx
         Xy0 = d.B + (1.0 - xi2) * alpy + xi2 * Cy

         dXx = Cx - alpx
         dXy = Cy - alpy

         th = muladd(ak, xi2, bk)

         gx, gy, dgx, dgy,d2gx, d2gy,d3gx, 
         d3gy,d4gx, d4gy = gamderhigher(d, th)

         Yx0 = d.A + gx
         Yy0 = d.B + gy

      elseif regloc == 6

         ak = π / d.P - d.tht2 - d.tht1
         bk = (2.0 * jpet - 3.0) * π / d.P + d.tht1

         Xx0 = d.A + d.R * e2x -
               d.L1 * e1x * invn / 2.0 +
               xi2 * d.L2 * qx * invn

         Xy0 = d.B + d.R * e2y -
               d.L1 * e1y * invn / 2.0 +
               xi2 * d.L2 * qy * invn

         dXx = d.L2 * qx * invn
         dXy = d.L2 * qy * invn

         th = muladd(ak, xi2, bk)

         gx, gy, dgx, dgy,d2gx, d2gy,d3gx, 
         d3gy,d4gx, d4gy = gamderhigher(d, th)

         Yx0 = d.A + gx
         Yy0 = d.B + gy

      else

         throw(ArgumentError("diff_map! Taylor fixup expects non-rectangular star region; got reg=$reg, regloc=$regloc"))

      end

   else

      rdisc = reg - 7 * d.P

      ak = π / 2.0
      bk = (2.0 * rdisc - 3.0) * π / 4.0

      th = muladd(ak, xi2, bk)
      sth, cth = sincos(th)

      Yx0 = d.A + d.Rc * cth
      Yy0 = d.B + d.Rc * sth

      dgx = -d.Rc * sth
      dgy = d.Rc * cth
      d2gx = -d.Rc * cth
      d2gy = -d.Rc * sth
      d3gx = d.Rc * sth
      d3gy = -d.Rc * cth
      d4gx = d.Rc * cth
      d4gy = d.Rc * sth

      if rdisc == 1

         Xx0 = d.A + d.Lc / 2.0
         Xy0 = d.B + d.Lc * (2.0 * xi2 - 1.0) / 2.0

         dXx = 0.0
         dXy = d.Lc

      elseif rdisc == 2

         Xx0 = d.A + d.Lc * (1.0 - 2.0 * xi2) / 2.0
         Xy0 = d.B + d.Lc / 2.0

         dXx = -d.Lc
         dXy = 0.0

      elseif rdisc == 3

         Xx0 = d.A - d.Lc / 2.0
         Xy0 = d.B - d.Lc * (2.0 * xi2 - 1.0) / 2.0

         dXx = 0.0
         dXy = -d.Lc

      elseif rdisc == 4

         Xx0 = d.A - d.Lc * (1.0 - 2.0 * xi2) / 2.0
         Xy0 = d.B - d.Lc / 2.0

         dXx = d.Lc
         dXy = 0.0

      else

         throw(ArgumentError("diff_map! Taylor fixup expects central region 1–4 or square 5; got rdisc=$rdisc"))

      end

   end

   # ------------------------------------------------------------
   # Common derivative assembly.
   # X is linear in xi2, so higher X derivatives vanish.
   # ------------------------------------------------------------

   ak2 = ak * ak
   ak3 = ak2 * ak
   ak4 = ak2 * ak2

   αt2 = αt * αt
   αt3 = αt2 * αt
   αt4 = αt2 * αt2

   dux = αc * (Yx0 - Xx0)
   duy = αc * (Yy0 - Xy0)

   dvx = αt * ((1.0 - xi1) * dXx + xi1 * ak * dgx)
   dvy = αt * ((1.0 - xi1) * dXy + xi1 * ak * dgy)

   duvx = αc * αt * (ak * dgx - dXx)
   duvy = αc * αt * (ak * dgy - dXy)

   dv2x = αt2 * xi1 * ak2 * d2gx
   dv2y = αt2 * xi1 * ak2 * d2gy

   duv2x = αc * αt2 * ak2 * d2gx
   duv2y = αc * αt2 * ak2 * d2gy

   dv3x = αt3 * xi1 * ak3 * d3gx
   dv3y = αt3 * xi1 * ak3 * d3gy

   duv3x = αc * αt3 * ak3 * d3gx
   duv3y = αc * αt3 * ak3 * d3gy

   dv4x = αt4 * xi1 * ak4 * d4gx
   dv4y = αt4 * xi1 * ak4 * d4gy

   tux, tvy = mapxy(d, u, v, k)

   @inbounds for j in 1:nd_v

      dvj = dv[j]
      dv2 = dvj * dvj
      dv3 = dv2 * dvj
      dv4 = dv2 * dv2

      @inbounds for i in 1:nd_u

         uu = u2[i, j]
         vv = v2[i, j]

         if abs(u - uu) < tol && abs(v - vv) < tol

            dui = du[i]

            Dx = (dui * dux + dvj * dvx) -
                 (dui * dvj * duvx + dv2 * dv2x / 2.0) +
                 (dui * dv2 * duv2x / 2.0 + dv3 * dv3x / 6.0) -
                 (dui * dv3 * duv3x / 6.0 + dv4 * dv4x / 24.0)

            Dy = (dui * duy + dvj * dvy) -
                 (dui * dvj * duvy + dv2 * dv2y / 2.0) +
                 (dui * dv2 * duv2y / 2.0 + dv3 * dv3y / 6.0) -
                 (dui * dv3 * duv3y / 6.0 + dv4 * dv4y / 24.0)

            out[i, j] = hypot(Dx, Dy)

         else

            out[i, j] = hypot(tux - Zx[i, j], tvy - Zy[i, j])

         end

      end

   end

   return nothing

end

"""
    gamderhigher6(d::star, th::Float64)

Return centered star boundary values and derivatives through sixth order.

Returns gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy, d4gx, d4gy, d5gx, d5gy, d6gx, d6gy

where gx(th) = R*h(th)*cos(th), gy(th) = R*h(th)*sin(th).
"""
function gamderhigher6(d::star, th::Float64)

   s, c = sincos(th)

   h, h1, h2, h3, h4, h5, h6 = hderhigher(d, th)

   gx = d.R * h * c
   gy = d.R * h * s

   dgx = d.R * (h1 * c - h * s)
   dgy = d.R * (h1 * s + h * c)

   d2gx = d.R * (h2 * c - 2.0 * h1 * s - h * c)
   d2gy = d.R * (h2 * s + 2.0 * h1 * c - h * s)

   d3gx = d.R * (h3 * c - 3.0 * h2 * s - 3.0 * h1 * c + h * s)
   d3gy = d.R * (h3 * s + 3.0 * h2 * c - 3.0 * h1 * s - h * c)

   d4gx = d.R * (h4 * c - 4.0 * h3 * s - 6.0 * h2 * c +
                  4.0 * h1 * s + h * c)

   d4gy = d.R * (h4 * s + 4.0 * h3 * c - 6.0 * h2 * s -
                  4.0 * h1 * c + h * s)

   d5gx = d.R * (h5 * c - 5.0 * h4 * s - 10.0 * h3 * c +
                  10.0 * h2 * s + 5.0 * h1 * c - h * s)

   d5gy = d.R * (h5 * s + 5.0 * h4 * c - 10.0 * h3 * s -
                  10.0 * h2 * c + 5.0 * h1 * s + h * c)

   d6gx = d.R * (h6 * c - 6.0 * h5 * s - 15.0 * h4 * c +
                  20.0 * h3 * s + 15.0 * h2 * c -
                  6.0 * h1 * s - h * c)

   d6gy = d.R * (h6 * s + 6.0 * h5 * c - 15.0 * h4 * s -
                  20.0 * h3 * c + 15.0 * h2 * s +
                  6.0 * h1 * c - h * s)

   return gx, gy,
          dgx, dgy,
          d2gx, d2gy,
          d3gx, d3gy,
          d4gx, d4gy,
          d5gx, d5gy,
          d6gx, d6gy

end

"""
    diff_rmap!(out::Matrix{Float64}, Zx::Matrix{Float64}, Zy::Matrix{Float64},
               DJ::StridedArray{Float64}, d::star,
               u::Float64, v::Float64,
               u2::Matrix{Float64}, v2::Matrix{Float64}, r::Matrix{Float64},
               du::AbstractVector, dv::AbstractVector, k::Int;
               tol::Float64 = 1e-3)

Compute ‖(τ(u,v) - τ(u₂,v₂)) / r‖ for the `k`-th patch, where

    u₂ = u - r .* du
    v₂ = v - r .* dv

No allocations. Uses `mapxy_Dmap!` and a sixth-order Taylor correction near `(u,v)`.

For affine/rectangular regions, only `DJ` is updated; `Zx,Zy` are not updated.
"""
function diff_rmap!(out::Matrix{Float64},
   Zx::Matrix{Float64}, Zy::Matrix{Float64}, DJ::StridedArray{Float64},
   d::star, u::Float64, v::Float64,
   u2::Matrix{Float64}, v2::Matrix{Float64}, r::Matrix{Float64},
   du::AbstractVector, dv::AbstractVector, k::Int;
   tol::Float64=1e-3)

   nt = size(out, 1)
   nr = size(out, 2)

   # @assert length(du) == nt
   # @assert length(dv) == nt
   # @assert size(r) == (nt, nr)
   # @assert size(out) == (nt, nr)
   # @assert size(u2)  == size(out)
   # @assert size(v2)  == size(out)
   # @assert size(du) == size(dv)
   # @inbounds for i in 1:nt, j in 1:nr
   #   @assert isapprox(u2[i, j], u - r[i, j] * du[i]; rtol=0, atol=1e-14)
   #   @assert isapprox(v2[i, j], v - r[i, j] * dv[i]; rtol=0, atol=1e-14)
   # end

   p = d.pths[k]
   reg = p.reg

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   αc = 0.5 * hc
   αt = 0.5 * ht
   βc = p.ck0 + αc
   βt = p.tk0 + αt

   # ------------------------------------------------------------
   # Petal rectangular region: exact affine distance.
   # ------------------------------------------------------------
   if reg <= 7 * d.P && mod1(reg, 7) == 7

      jpet = cld(reg, 7)

      e1x = d.RP.e₁x[jpet]
      e1y = d.RP.e₁y[jpet]
      qx = d.RP.qx[jpet]
      qy = d.RP.qy[jpet]
      invn = 1.0 / d.RP.nme

      fill!(DJ, d.L1 * d.L2 * αc * αt)

      if d.L2 >= d.L1
         cU = d.L2 * αt * invn
         cV = d.L1 * αc * invn
      else
         cU = d.L2 * αc * invn
         cV = d.L1 * αt * invn
      end

      @inbounds for i in 1:nt
         dui = du[i]
         dvi = dv[i]

         hD = hypot(cU * dui * qx + cV * dvi * e1x,
            cU * dui * qy + cV * dvi * e1y)

         @inbounds for j in 1:nr
            out[i, j] = hD
         end
      end

      return nothing

   end

   # ------------------------------------------------------------
   # Central square: exact affine distance.
   # ------------------------------------------------------------
   if reg == 7 * d.P + 5

      cU = d.Lc * αc
      cV = d.Lc * αt

      fill!(DJ, cU * cV)

      @inbounds for i in 1:nt
         dui = du[i]
         dvi = dv[i]

         hD = hypot(cU * dui, cV * dvi)

         @inbounds for j in 1:nr
            out[i, j] = hD
         end
      end

      return nothing

   end

   # Non-affine regions: compute grid map/Jacobian.
   mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)

   # Scalar evaluation point in region coordinates.
   xi1 = muladd(αc, u, βc)
   xi2 = muladd(αt, v, βt)

   if reg <= 7 * d.P

      jpet = cld(reg, 7)
      regloc = mod1(reg, 7)

      e1x = d.RP.e₁x[jpet]
      e1y = d.RP.e₁y[jpet]
      e2x = d.RP.e₂x[jpet]
      e2y = d.RP.e₂y[jpet]
      qx = d.RP.qx[jpet]
      qy = d.RP.qy[jpet]
      invn = 1.0 / d.RP.nme

      if regloc == 1

         ak = 2.0 * d.tht2
         bk = 2.0 * (jpet - 1) * π / d.P - d.tht2

         Xx0 = d.A + d.R * e2x +
               (2.0 * xi2 - 1.0) * d.L1 * e1x * invn / 2.0 +
               d.L2 * qx * invn

         Xy0 = d.B + d.R * e2y +
               (2.0 * xi2 - 1.0) * d.L1 * e1y * invn / 2.0 +
               d.L2 * qy * invn

         dXx = d.L1 * e1x * invn
         dXy = d.L1 * e1y * invn

         th = muladd(ak, xi2, bk)

         gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy,
         d4gx, d4gy, d5gx, d5gy, d6gx, d6gy = gamderhigher6(d, th)

         Yx0 = d.A + gx
         Yy0 = d.B + gy

      elseif regloc == 2

         ak = π / d.P - d.tht2 - d.tht1
         bk = 2.0 * (jpet - 1) * π / d.P + d.tht2

         Xx0 = d.A + d.R * e2x +
               d.L1 * e1x * invn / 2.0 +
               (1.0 - xi2) * d.L2 * qx * invn

         Xy0 = d.B + d.R * e2y +
               d.L1 * e1y * invn / 2.0 +
               (1.0 - xi2) * d.L2 * qy * invn

         dXx = -d.L2 * qx * invn
         dXy = -d.L2 * qy * invn

         th = muladd(ak, xi2, bk)

         gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy,
         d4gx, d4gy, d5gx, d5gy, d6gx, d6gy = gamderhigher6(d, th)

         Yx0 = d.A + gx
         Yy0 = d.B + gy

      elseif regloc == 3

         ak = d.tht1
         bk = (2.0 * jpet - 1.0) * π / d.P - d.tht1

         alpx = d.RP.αx[jpet+1]
         alpy = d.RP.αy[jpet+1]

         Cx = d.R * e2x + d.L1 * e1x * invn / 2.0
         Cy = d.R * e2y + d.L1 * e1y * invn / 2.0

         Xx0 = d.A + xi2 * alpx + (1.0 - xi2) * Cx
         Xy0 = d.B + xi2 * alpy + (1.0 - xi2) * Cy

         dXx = alpx - Cx
         dXy = alpy - Cy

         th = muladd(ak, xi2, bk)

         gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy,
         d4gx, d4gy, d5gx, d5gy, d6gx, d6gy = gamderhigher6(d, th)

         Yx0 = d.A + gx
         Yy0 = d.B + gy

      elseif regloc == 4

         ak = 2.0 * π / d.P
         bk = (2.0 * jpet - 3.0) * π / d.P

         Xx0 = d.A + d.R * e2x +
               (2.0 * xi2 - 1.0) * d.L1 * e1x * invn / 2.0

         Xy0 = d.B + d.R * e2y +
               (2.0 * xi2 - 1.0) * d.L1 * e1y * invn / 2.0

         dXx = d.L1 * e1x * invn
         dXy = d.L1 * e1y * invn

         th = muladd(ak, xi2, bk)
         sth, cth = sincos(th)

         Yx0 = d.A + d.Rc * cth
         Yy0 = d.B + d.Rc * sth

         dgx = -d.Rc * sth
         dgy = d.Rc * cth
         d2gx = -d.Rc * cth
         d2gy = -d.Rc * sth
         d3gx = d.Rc * sth
         d3gy = -d.Rc * cth
         d4gx = d.Rc * cth
         d4gy = d.Rc * sth
         d5gx = -d.Rc * sth
         d5gy = d.Rc * cth
         d6gx = -d.Rc * cth
         d6gy = -d.Rc * sth

      elseif regloc == 5

         ak = d.tht1
         bk = (2.0 * jpet - 3.0) * π / d.P

         alpx = d.RP.αx[jpet]
         alpy = d.RP.αy[jpet]

         Cx = d.R * e2x - d.L1 * e1x * invn / 2.0
         Cy = d.R * e2y - d.L1 * e1y * invn / 2.0

         Xx0 = d.A + (1.0 - xi2) * alpx + xi2 * Cx
         Xy0 = d.B + (1.0 - xi2) * alpy + xi2 * Cy

         dXx = Cx - alpx
         dXy = Cy - alpy

         th = muladd(ak, xi2, bk)

         gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy,
         d4gx, d4gy, d5gx, d5gy, d6gx, d6gy = gamderhigher6(d, th)

         Yx0 = d.A + gx
         Yy0 = d.B + gy

      elseif regloc == 6

         ak = π / d.P - d.tht2 - d.tht1
         bk = (2.0 * jpet - 3.0) * π / d.P + d.tht1

         Xx0 = d.A + d.R * e2x -
               d.L1 * e1x * invn / 2.0 +
               xi2 * d.L2 * qx * invn

         Xy0 = d.B + d.R * e2y -
               d.L1 * e1y * invn / 2.0 +
               xi2 * d.L2 * qy * invn

         dXx = d.L2 * qx * invn
         dXy = d.L2 * qy * invn

         th = muladd(ak, xi2, bk)

         gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy,
         d4gx, d4gy, d5gx, d5gy, d6gx, d6gy = gamderhigher6(d, th)

         Yx0 = d.A + gx
         Yy0 = d.B + gy

      else

         throw(ArgumentError("diff_rmap! expects non-rectangular star petal region; got reg=$reg, regloc=$regloc"))

      end

   else

      rdisc = reg - 7 * d.P

      ak = π / 2.0
      bk = (2.0 * rdisc - 3.0) * π / 4.0

      th = muladd(ak, xi2, bk)
      sth, cth = sincos(th)

      Yx0 = d.A + d.Rc * cth
      Yy0 = d.B + d.Rc * sth

      dgx = -d.Rc * sth
      dgy = d.Rc * cth
      d2gx = -d.Rc * cth
      d2gy = -d.Rc * sth
      d3gx = d.Rc * sth
      d3gy = -d.Rc * cth
      d4gx = d.Rc * cth
      d4gy = d.Rc * sth
      d5gx = -d.Rc * sth
      d5gy = d.Rc * cth
      d6gx = -d.Rc * cth
      d6gy = -d.Rc * sth

      if rdisc == 1

         Xx0 = d.A + d.Lc / 2.0
         Xy0 = d.B + d.Lc * (2.0 * xi2 - 1.0) / 2.0

         dXx = 0.0
         dXy = d.Lc

      elseif rdisc == 2

         Xx0 = d.A + d.Lc * (1.0 - 2.0 * xi2) / 2.0
         Xy0 = d.B + d.Lc / 2.0

         dXx = -d.Lc
         dXy = 0.0

      elseif rdisc == 3

         Xx0 = d.A - d.Lc / 2.0
         Xy0 = d.B - d.Lc * (2.0 * xi2 - 1.0) / 2.0

         dXx = 0.0
         dXy = -d.Lc

      elseif rdisc == 4

         Xx0 = d.A - d.Lc * (1.0 - 2.0 * xi2) / 2.0
         Xy0 = d.B - d.Lc / 2.0

         dXx = d.Lc
         dXy = 0.0

      else

         throw(ArgumentError("diff_rmap! expects central region 1–4 or square 5; got rdisc=$rdisc"))

      end

   end

   # ------------------------------------------------------------
   # Common derivative assembly.
   # Derivatives are with respect to local reference u,v.
   # ------------------------------------------------------------

   q = αt * ak

   q2 = q * q
   q3 = q2 * q
   q4 = q2 * q2
   q5 = q4 * q
   q6 = q3 * q3

   dux = αc * (Yx0 - Xx0)
   duy = αc * (Yy0 - Xy0)

   dvx = αt * ((1.0 - xi1) * dXx + xi1 * ak * dgx)
   dvy = αt * ((1.0 - xi1) * dXy + xi1 * ak * dgy)

   duvx = αc * (q * dgx - αt * dXx)
   duvy = αc * (q * dgy - αt * dXy)

   dv2x = xi1 * q2 * d2gx
   dv2y = xi1 * q2 * d2gy

   duv2x = αc * q2 * d2gx
   duv2y = αc * q2 * d2gy

   dv3x = xi1 * q3 * d3gx
   dv3y = xi1 * q3 * d3gy

   duv3x = αc * q3 * d3gx
   duv3y = αc * q3 * d3gy

   dv4x = xi1 * q4 * d4gx
   dv4y = xi1 * q4 * d4gy

   duv4x = αc * q4 * d4gx
   duv4y = αc * q4 * d4gy

   dv5x = xi1 * q5 * d5gx
   dv5y = xi1 * q5 * d5gy

   duv5x = αc * q5 * d5gx
   duv5y = αc * q5 * d5gy

   dv6x = xi1 * q6 * d6gx
   dv6y = xi1 * q6 * d6gy

   tux = (1.0 - xi1) * Xx0 + xi1 * Yx0
   tvy = (1.0 - xi1) * Xy0 + xi1 * Yy0

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
    Dwall(d::star, u::Float64, v::Float64, k::Int) -> dvx, dvy

Derivative of the wall curve with respect to `v` at side `u = ±1`
for the `k`-th non-rectangular patch.

Returns `(dvx, dvy)`.
"""
function Dwall(d::star, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]
   reg = p.reg

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   αu = 0.5 * hc
   αv = 0.5 * ht
   βu = p.ck0 + αu
   βv = p.tk0 + αv

   xi1 = muladd(αu, u, βu)
   xi2 = muladd(αv, v, βv)

   if reg <= 7 * d.P

      jpet = cld(reg, 7)
      regloc = mod1(reg, 7)

      e1x = d.RP.e₁x[jpet]
      e1y = d.RP.e₁y[jpet]
      e2x = d.RP.e₂x[jpet]
      e2y = d.RP.e₂y[jpet]
      qx = d.RP.qx[jpet]
      qy = d.RP.qy[jpet]
      nme = d.RP.nme
      invn = inv(nme)

      if regloc == 1

         Δth = 2.0 * d.tht2
         th0 = 2.0 * (jpet - 1) * π / d.P - d.tht2
         th = muladd(Δth, xi2, th0)

         _, _, dgx, dgy = dergam(d, th)

         Xpx = d.L1 * e1x * invn
         Xpy = d.L1 * e1y * invn

         dvx = (1.0 - xi1) * Xpx + xi1 * Δth * dgx
         dvy = (1.0 - xi1) * Xpy + xi1 * Δth * dgy

      elseif regloc == 2

         Δth = π / d.P - d.tht2 - d.tht1
         th0 = 2.0 * (jpet - 1) * π / d.P + d.tht2
         th = muladd(Δth, xi2, th0)

         _, _, dgx, dgy = dergam(d, th)

         Xpx = -d.L2 * qx * invn
         Xpy = -d.L2 * qy * invn

         dvx = (1.0 - xi1) * Xpx + xi1 * Δth * dgx
         dvy = (1.0 - xi1) * Xpy + xi1 * Δth * dgy

      elseif regloc == 3

         Δth = d.tht1
         th0 = (2.0 * jpet - 1.0) * π / d.P - d.tht1
         th = muladd(Δth, xi2, th0)

         _, _, dgx, dgy = dergam(d, th)

         alpx = d.RP.αx[jpet + 1]
         alpy = d.RP.αy[jpet + 1]

         Cx = d.R * e2x + d.L1 * e1x * invn / 2.0
         Cy = d.R * e2y + d.L1 * e1y * invn / 2.0

         Xpx = alpx - Cx
         Xpy = alpy - Cy

         dvx = (1.0 - xi1) * Xpx + xi1 * Δth * dgx
         dvy = (1.0 - xi1) * Xpy + xi1 * Δth * dgy

      elseif regloc == 4

         Δth = 2.0 * π / d.P
         th0 = (2.0 * jpet - 3.0) * π / d.P
         th = muladd(Δth, xi2, th0)

         st, ct = sincos(th)

         Xpx = d.L1 * e1x * invn
         Xpy = d.L1 * e1y * invn

         Ypx = -d.Rc * st
         Ypy =  d.Rc * ct

         dvx = (1.0 - xi1) * Xpx + xi1 * Δth * Ypx
         dvy = (1.0 - xi1) * Xpy + xi1 * Δth * Ypy

      elseif regloc == 5

         Δth = d.tht1
         th0 = (2.0 * jpet - 3.0) * π / d.P
         th = muladd(Δth, xi2, th0)

         _, _, dgx, dgy = dergam(d, th)

         alpx = d.RP.αx[jpet]
         alpy = d.RP.αy[jpet]

         Cx = d.R * e2x - d.L1 * e1x * invn / 2.0
         Cy = d.R * e2y - d.L1 * e1y * invn / 2.0

         Xpx = Cx - alpx
         Xpy = Cy - alpy

         dvx = (1.0 - xi1) * Xpx + xi1 * Δth * dgx
         dvy = (1.0 - xi1) * Xpy + xi1 * Δth * dgy

      elseif regloc == 6

         Δth = π / d.P - d.tht2 - d.tht1
         th0 = (2.0 * jpet - 3.0) * π / d.P + d.tht1
         th = muladd(Δth, xi2, th0)

         _, _, dgx, dgy = dergam(d, th)

         Xpx = d.L2 * qx * invn
         Xpy = d.L2 * qy * invn

         dvx = (1.0 - xi1) * Xpx + xi1 * Δth * dgx
         dvy = (1.0 - xi1) * Xpy + xi1 * Δth * dgy

      else

         throw(ArgumentError("Dwall(star) is for non-rectangular petal regions; got reg=$reg, regloc=$regloc"))

      end

      return αv * dvx, αv * dvy

   else

      rdisc = reg - 7 * d.P

      if rdisc == 1

         th = π * xi2 / 2.0 - π / 4.0
         st, ct = sincos(th)

         Xpx = 0.0
         Xpy = d.Lc

         Ypx = -d.Rc * st
         Ypy =  d.Rc * ct

         dvx = (1.0 - xi1) * Xpx + xi1 * (π / 2.0) * Ypx
         dvy = (1.0 - xi1) * Xpy + xi1 * (π / 2.0) * Ypy

      elseif rdisc == 2

         th = π * xi2 / 2.0 + π / 4.0
         st, ct = sincos(th)

         Xpx = -d.Lc
         Xpy = 0.0

         Ypx = -d.Rc * st
         Ypy =  d.Rc * ct

         dvx = (1.0 - xi1) * Xpx + xi1 * (π / 2.0) * Ypx
         dvy = (1.0 - xi1) * Xpy + xi1 * (π / 2.0) * Ypy

      elseif rdisc == 3

         th = π * xi2 / 2.0 + 3.0 * π / 4.0
         st, ct = sincos(th)

         Xpx = 0.0
         Xpy = -d.Lc

         Ypx = -d.Rc * st
         Ypy =  d.Rc * ct

         dvx = (1.0 - xi1) * Xpx + xi1 * (π / 2.0) * Ypx
         dvy = (1.0 - xi1) * Xpy + xi1 * (π / 2.0) * Ypy

      elseif rdisc == 4

         th = π * xi2 / 2.0 + 5.0 * π / 4.0
         st, ct = sincos(th)

         Xpx = d.Lc
         Xpy = 0.0

         Ypx = -d.Rc * st
         Ypy =  d.Rc * ct

         dvx = (1.0 - xi1) * Xpx + xi1 * (π / 2.0) * Ypx
         dvy = (1.0 - xi1) * Xpy + xi1 * (π / 2.0) * Ypy

      else

         throw(ArgumentError("Dwall(star) is for non-rectangular central regions; got reg=$reg, rdisc=$rdisc"))

      end

      return αv * dvx, αv * dvy

   end

end

"""
    refine!(d::star, Nc::Int, Nt::Int, K::Vector{Int})

Subdivide each patch in `K` into `Nc * Nt` subpatches, update `d.pths`,
`d.Npat`, recompute `d.kd`, and rebuild `d.Qpts` / `d.Qptsbd`.
"""
function refine!(d::star, Nc::Int, Nt::Int, K::Vector{Int})

   @assert Nc ≥ 1 && Nt ≥ 1 "Nc and Nt must be ≥ 1"
   @assert !isempty(K) "K cannot be empty"

   Kset = Set(K)

   children = Patch[]
   sizehint!(children, length(K) * Nc * Nt)

   @inbounds for k in K

      p = d.pths[k]

      cc = range(p.ck0, p.ck1; length = Nc + 1)
      tt = range(p.tk0, p.tk1; length = Nt + 1)

      for i in 1:Nc, j in 1:Nt
         push!(children, Patch(p.reg, cc[i], cc[i+1], tt[j], tt[j+1]))
      end

   end

   d.pths = vcat([d.pths[i] for i in 1:d.Npat if !(i in Kset)], children)

   sort!(d.pths, by = q -> (q.reg, q.ck0, q.ck1, q.tk0, q.tk1))

   d.Npat = length(d.pths)

   d.kd = [
      k for k in 1:d.Npat
      if d.pths[k].reg <= 7 * d.P &&
         mod1(d.pths[k].reg, 7) != 4 &&
         mod1(d.pths[k].reg, 7) != 7 &&
         d.pths[k].ck1 == 1.0]

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

   return nothing

end

function Base.show(io::IO, d::star)

   println(io, "star with properties:")
   println(io, "  (A , B )  = (", d.A, ", ", d.B, ")")
   println(io, "  R         = ", d.R)
   println(io, "  P         = ", d.P)
   println(io, "  Ps        = ", d.Ps)
   println(io, "  Rc        = ", d.Rc)
   println(io, "  Lc        = ", d.Lc)
   println(io, "  (L₁, L₂)  = (", d.L1, ", ", d.L2, ")")
   println(io, "  (θ₁, θ₂)  = (", d.tht1, ", ", d.tht2, ")")
   println(io, "  No. holes = ", d.nh)

   if length(d.kd) <= 6
      L = "Int[" * join(d.kd, ' ') * "]"
   else
      head = join(d.kd[1:3], ' ')
      tail = join(d.kd[end-2:end], ' ')
      L = "Int[$head … $tail]"
   end

   println(io, "  kd     = ", L)
   println(io, "  Npat   = ", d.Npat)
   println(io, "  pths   = Vector{Patch} (", length(d.pths), " patches)")
   println(io, "  Qpts   = ", size(d.Qpts, 1), "×", size(d.Qpts, 2), " Matrix")
   println(io, "  Qptsbd = ", size(d.Qptsbd, 1), "×", size(d.Qptsbd, 2), " Matrix")

end
