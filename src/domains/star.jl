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