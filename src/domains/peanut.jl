"""
Mutable struct peanut

A     — x-coordinate of peanut center
B     — y-coordinate of peanut center
R     — scale/radius of peanut
P     — pinch factor, positive
L1    — length of y-axis side of each inscribed rectangle
L2    — length of x-axis side of each inscribed rectangle
nh    — number of holes
tht1  — construction angle, default 11π/100
tht2  — computed construction angle from vertical-line intersection
alpha — α = (3 + hP(π/2)) / 4

kd    — indices of patches touching boundary
Npat  — total number of patches
pths  — Vector{Patch}
Qpts  — 8 × Npat matrix
Qptsbd — 8 × length(kd) matrix
"""
mutable struct peanut <: abstractdomain

   A::Float64
   B::Float64
   R::Float64
   P::Float64
   L1::Float64
   L2::Float64
   nh::Int
   tht1::Float64
   tht2::Float64
   alpha::Float64
   kd::Vector{Int}
   Npat::Int
   pths::Vector{Patch}
   Qpts::Matrix{Float64}
   Qptsbd::Matrix{Float64}

   function peanut(; b,
      a = nothing,
      A = nothing, B = nothing,
      R = nothing, P = nothing,
      L1 = nothing, L2 = nothing,
      tht1 = nothing, ck = nothing, tk = nothing)

      @assert isa(b, AbstractVector{<:Int}) && length(b) == 12
      @assert all(b .> 0)

      nh = 0

      a = isnothing(a) ? ceil.(Int, 2 .* b ./ 3) : collect(Int, a)

      @assert length(a) == 12
      @assert all(a .> 0)

      A = Float64(something(A, 0.0))
      B = Float64(something(B, 0.0))
      R = Float64(something(R, 1.0))
      P = Float64(something(P, 0.5))
      tht1 = Float64(something(tht1, 11π / 100))

      @assert R > 0.0
      @assert P > 0.0
      @assert 0.0 < tht1 < π / 2

      L1 = Float64(something(L1, 0.5 * R))
      L2 = Float64(something(L2, L1))

      @assert L1 > 0.0
      @assert L2 > 0.0

      # hP(t) = sqrt(sqrt(P + cos(2t)^2) + cos(2t))
      hP(t) = sqrt(sqrt(P + cos(2t)^2) + cos(2t))

      # hP'(t)
      dhP(t) = begin
         c2 = cos(2t)
         -sin(2t) * (1 + c2 / sqrt(P + c2^2)) / hP(t)
      end

      alpha = (3 + hP(π / 2)) / 4

      # tht2 solves:
      # hP(t) * cos(t) - alpha + L2/(2R) = 0
      f(t) = hP(t) * cos(t) - alpha + L2 / (2R)
      df(t) = dhP(t) * cos(t) - hP(t) * sin(t)

      tht2 = newtonR1D(f, df, π / 4, 64) #Last entry means max iterations

      @assert tht1 < tht2 < π / 2 "Expected tht1 < tht2 < π/2; got tht1=$tht1, tht2=$tht2"

      Npat = dot(a, b)

      # ck[k] has length a[k] + 1
      # tk[k] has length b[k] + 1
      ck = something(ck, [(0:a[k]) ./ a[k] for k in 1:12])
      tk = something(tk, [(0:b[k]) ./ b[k] for k in 1:12])

      @assert length(ck) == 12
      @assert length(tk) == 12

      pths = Vector{Patch}(undef, Npat)

      kd = Vector{Int}(undef, Npat)
      nkd = 0

      idx = 1

      @inbounds for k in 1:12
         ak = a[k]
         bk = b[k]
         ckₖ = ck[k]
         tkₖ = tk[k]

         @assert length(ckₖ) == ak + 1
         @assert length(tkₖ) == bk + 1

         for i in 1:ak
            ck0 = Float64(ckₖ[i])
            ck1 = Float64(ckₖ[i + 1])

            for j in 1:bk
               tk0 = Float64(tkₖ[j])
               tk1 = Float64(tkₖ[j + 1])

               p = Patch(k, ck0, ck1, tk0, tk1)
               pths[idx] = p

               # Boundary patches are in regions 1:10 and have ck1 == 1.
               if p.reg <= 10 && p.ck1 == 1.0
                  nkd += 1
                  kd[nkd] = idx
               end

               idx += 1
            end
         end
      end

      resize!(kd, nkd)

      Qpts = Matrix{Float64}(undef, 8, Npat)
      Qptsbd = Matrix{Float64}(undef, 8, nkd)

      d = new(A, B, R, P, L1, L2, nh, tht1, tht2, alpha,
              kd, Npat, pths, Qpts, Qptsbd)

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

function mapx(d::peanut, u::Float64, v::Float64, k::Int)::Float64

   p = d.pths[k]

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   xi1 = hc * u / 2 + (p.ck1 + p.ck0) / 2
   xi2 = ht * v / 2 + (p.tk1 + p.tk0) / 2

   r = p.reg

   if r != 11 && r != 12

      X, theta = if r == 1

         d.A + d.R * d.alpha + d.L2 / 2,
         (2 * xi2 - 1) * d.tht1

      elseif r == 2

         d.A + d.R * d.alpha + d.L2 * (1 - 2 * xi2) / 2,
         xi2 * (d.tht2 - d.tht1) + d.tht1

      elseif r == 3

         d.A + (1 - xi2) * (d.R * d.alpha - d.L2 / 2),
         xi2 * (π / 2 - d.tht2) + d.tht2

      elseif r == 4

         d.A + xi2 * (d.L2 / 2 - d.R * d.alpha),
         xi2 * (π / 2 - d.tht2) + π / 2

      elseif r == 5

         d.A - d.R * d.alpha + d.L2 * (1 - 2 * xi2) / 2,
         xi2 * (d.tht2 - d.tht1) + π - d.tht2

      elseif r == 6

         d.A - d.R * d.alpha - d.L2 / 2,
         2 * xi2 * d.tht1 + π - d.tht1

      elseif r == 7

         d.A - d.R * d.alpha - d.L2 * (1 - 2 * xi2) / 2,
         xi2 * (d.tht2 - d.tht1) + π + d.tht1

      elseif r == 8

         d.A + (1 - xi2) * (d.L2 / 2 - d.R * d.alpha),
         xi2 * (π / 2 - d.tht2) + π + d.tht2

      elseif r == 9

         d.A + xi2 * (d.R * d.alpha - d.L2 / 2),
         xi2 * (π / 2 - d.tht2) - π / 2

      elseif r == 10

         d.A + d.R * d.alpha - d.L2 * (1 - 2 * xi2) / 2,
         xi2 * (d.tht2 - d.tht1) - d.tht2

      else

         throw(ArgumentError("mapx expects region 1–12; got reg=$r"))

      end

      st, ct = sincos(theta)
      f = ct^2 - st^2
      q = sqrt(d.P + f^2)
      hP = sqrt(f + q)

      Y = d.A + d.R * hP * ct

      return (1 - xi1) * X + xi1 * Y

   elseif r == 11

      if d.L2 >= d.L1
         xi1 = ht * u / 2 + (p.tk1 + p.tk0) / 2
      end

      return d.A + d.R * d.alpha - d.L2 / 2 + d.L2 * xi1

   else # r == 12

      if d.L2 >= d.L1
         xi1 = ht * u / 2 + (p.tk1 + p.tk0) / 2
      end

      return d.A - d.R * d.alpha - d.L2 / 2 + d.L2 * xi1

   end

end

function mapy(d::peanut, u::Float64, v::Float64, k::Int)::Float64

   p = d.pths[k]

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   xi1 = hc * u / 2 + (p.ck1 + p.ck0) / 2
   xi2 = ht * v / 2 + (p.tk1 + p.tk0) / 2

   r = p.reg

   if r != 11 && r != 12

      X, theta = if r == 1

         d.B + d.L1 * (2 * xi2 - 1) / 2,
         (2 * xi2 - 1) * d.tht1

      elseif r == 2

         d.B + d.L1 / 2,
         xi2 * (d.tht2 - d.tht1) + d.tht1

      elseif r == 3

         d.B,
         xi2 * (π / 2 - d.tht2) + d.tht2

      elseif r == 4

         d.B,
         xi2 * (π / 2 - d.tht2) + π / 2

      elseif r == 5

         d.B + d.L1 / 2,
         xi2 * (d.tht2 - d.tht1) + π - d.tht2

      elseif r == 6

         d.B - d.L1 * (2 * xi2 - 1) / 2,
         2 * xi2 * d.tht1 + π - d.tht1

      elseif r == 7

         d.B - d.L1 / 2,
         xi2 * (d.tht2 - d.tht1) + π + d.tht1

      elseif r == 8

         d.B,
         xi2 * (π / 2 - d.tht2) + π + d.tht2

      elseif r == 9

         d.B,
         xi2 * (π / 2 - d.tht2) - π / 2

      elseif r == 10

         d.B - d.L1 / 2,
         xi2 * (d.tht2 - d.tht1) - d.tht2

      else

         throw(ArgumentError("mapy expects region 1–12; got reg=$r"))

      end

      st, ct = sincos(theta)
      f = ct^2 - st^2
      q = sqrt(d.P + f^2)
      hP = sqrt(f + q)

      Y = d.B + d.R * hP * st

      return (1 - xi1) * X + xi1 * Y

   elseif r == 11 || r == 12

      if d.L2 >= d.L1
         xi2 = hc * v / 2 + (p.ck1 + p.ck0) / 2
      end

      return d.B - d.L1 / 2 + d.L1 * xi2

   else

      throw(ArgumentError("mapy expects region 1–12; got reg=$r"))

   end

end

function mapxy(d::peanut, u::Float64, v::Float64, k::Int)::Tuple{Float64, Float64}

   return mapx(d, u, v, k), mapy(d, u, v, k)

end

function mapxy!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
   d::peanut, u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0
   r = p.reg

   αu = hc / 2
   αv = ht / 2
   βu = p.ck0 + αu
   βv = p.tk0 + αv

   Rα = d.R * d.alpha
   L1h = d.L1 / 2
   L2h = d.L2 / 2

   if r == 1

      @inbounds for I in eachindex(u, v, Zx, Zy)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + Rα + L2h
         Xy = d.B + d.L1 * (2xi2 - 1) / 2

         theta = (2xi2 - 1) * d.tht1

         st, ct = sincos(theta)
         f = ct^2 - st^2
         hP = sqrt(f + sqrt(d.P + f^2))

         Yx = d.A + d.R * hP * ct
         Yy = d.B + d.R * hP * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)

      end

   elseif r == 2

      aθ = d.tht2 - d.tht1

      @inbounds for I in eachindex(u, v, Zx, Zy)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + Rα + d.L2 * (1 - 2xi2) / 2
         Xy = d.B + L1h

         theta = xi2 * aθ + d.tht1

         st, ct = sincos(theta)
         f = ct^2 - st^2
         hP = sqrt(f + sqrt(d.P + f^2))

         Yx = d.A + d.R * hP * ct
         Yy = d.B + d.R * hP * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)

      end

   elseif r == 3

      aθ = π / 2 - d.tht2
      Cx = Rα - L2h

      @inbounds for I in eachindex(u, v, Zx, Zy)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + (1 - xi2) * Cx
         Xy = d.B

         theta = xi2 * aθ + d.tht2

         st, ct = sincos(theta)
         f = ct^2 - st^2
         hP = sqrt(f + sqrt(d.P + f^2))

         Yx = d.A + d.R * hP * ct
         Yy = d.B + d.R * hP * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)

      end

   elseif r == 4

      aθ = π / 2 - d.tht2
      Cx = L2h - Rα

      @inbounds for I in eachindex(u, v, Zx, Zy)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + xi2 * Cx
         Xy = d.B

         theta = xi2 * aθ + π / 2

         st, ct = sincos(theta)
         f = ct^2 - st^2
         hP = sqrt(f + sqrt(d.P + f^2))

         Yx = d.A + d.R * hP * ct
         Yy = d.B + d.R * hP * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)

      end

   elseif r == 5

      aθ = d.tht2 - d.tht1

      @inbounds for I in eachindex(u, v, Zx, Zy)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A - Rα + d.L2 * (1 - 2xi2) / 2
         Xy = d.B + L1h

         theta = xi2 * aθ + π - d.tht2

         st, ct = sincos(theta)
         f = ct^2 - st^2
         hP = sqrt(f + sqrt(d.P + f^2))

         Yx = d.A + d.R * hP * ct
         Yy = d.B + d.R * hP * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)

      end

   elseif r == 6

      @inbounds for I in eachindex(u, v, Zx, Zy)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A - Rα - L2h
         Xy = d.B - d.L1 * (2xi2 - 1) / 2

         theta = 2xi2 * d.tht1 + π - d.tht1

         st, ct = sincos(theta)
         f = ct^2 - st^2
         hP = sqrt(f + sqrt(d.P + f^2))

         Yx = d.A + d.R * hP * ct
         Yy = d.B + d.R * hP * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)

      end

   elseif r == 7

      aθ = d.tht2 - d.tht1

      @inbounds for I in eachindex(u, v, Zx, Zy)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A - Rα - d.L2 * (1 - 2xi2) / 2
         Xy = d.B - L1h

         theta = xi2 * aθ + π + d.tht1

         st, ct = sincos(theta)
         f = ct^2 - st^2
         hP = sqrt(f + sqrt(d.P + f^2))

         Yx = d.A + d.R * hP * ct
         Yy = d.B + d.R * hP * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)

      end

   elseif r == 8

      aθ = π / 2 - d.tht2
      Cx = L2h - Rα

      @inbounds for I in eachindex(u, v, Zx, Zy)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + (1 - xi2) * Cx
         Xy = d.B

         theta = xi2 * aθ + π + d.tht2

         st, ct = sincos(theta)
         f = ct^2 - st^2
         hP = sqrt(f + sqrt(d.P + f^2))

         Yx = d.A + d.R * hP * ct
         Yy = d.B + d.R * hP * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)

      end

   elseif r == 9

      aθ = π / 2 - d.tht2
      Cx = Rα - L2h

      @inbounds for I in eachindex(u, v, Zx, Zy)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + xi2 * Cx
         Xy = d.B

         theta = xi2 * aθ - π / 2

         st, ct = sincos(theta)
         f = ct^2 - st^2
         hP = sqrt(f + sqrt(d.P + f^2))

         Yx = d.A + d.R * hP * ct
         Yy = d.B + d.R * hP * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)

      end

   elseif r == 10

      aθ = d.tht2 - d.tht1

      @inbounds for I in eachindex(u, v, Zx, Zy)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + Rα - d.L2 * (1 - 2xi2) / 2
         Xy = d.B - L1h

         theta = xi2 * aθ - d.tht2

         st, ct = sincos(theta)
         f = ct^2 - st^2
         hP = sqrt(f + sqrt(d.P + f^2))

         Yx = d.A + d.R * hP * ct
         Yy = d.B + d.R * hP * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)

      end

   elseif r == 11 || r == 12

      if d.L2 >= d.L1

         αu = ht / 2
         αv = hc / 2
         βu = p.tk0 + αu
         βv = p.ck0 + αv

      end

      x0 = if r == 11
         d.A + Rα - L2h
      else
         d.A - Rα - L2h
      end

      y0 = d.B - L1h

      @inbounds for I in eachindex(u, v, Zx, Zy)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Zx[I] = muladd(d.L2, xi1, x0)
         Zy[I] = muladd(d.L1, xi2, y0)

      end

   else

      throw(ArgumentError("mapxy! expects region 1–12; got reg=$r"))

   end

   return nothing

end

function draw(d::peanut, flag=nothing; L::Int=32, show::Bool=true)

   colors = (
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
      RGBf(0, 0, 1),                 # blue
      RGBf(0, 0, 1),                 # blue
   )

   return draw_geom(d, colors; flag=flag, L=L, show=show)

end

#-----------------------
"""
    gamx(d::peanut, t::Float64, k::Int)
    gamx!(out::StridedArray{Float64}, d::peanut, t::StridedArray{Float64}, k::Int)

First x-coordinate of the boundary parametrization gamma(t).

- `gamx(d, t, k)` computes the right boundary of patch `k`, where `t in [-1,1]`.
- Only defined for boundary regions `1:10`.

`t` may be a scalar or an array; the return/output has the same shape.
"""
function gamx(d::peanut, t::Float64, k::Int)::Float64

   p = d.pths[k]

   xi2 = (p.tk1 - p.tk0) * t / 2 + (p.tk0 + p.tk1) / 2

   theta = if p.reg == 1

      (2 * xi2 - 1) * d.tht1

   elseif p.reg == 2

      xi2 * (d.tht2 - d.tht1) + d.tht1

   elseif p.reg == 3

      xi2 * (π / 2 - d.tht2) + d.tht2

   elseif p.reg == 4

      xi2 * (π / 2 - d.tht2) + π / 2

   elseif p.reg == 5

      xi2 * (d.tht2 - d.tht1) + π - d.tht2

   elseif p.reg == 6

      2 * xi2 * d.tht1 + π - d.tht1

   elseif p.reg == 7

      xi2 * (d.tht2 - d.tht1) + π + d.tht1

   elseif p.reg == 8

      xi2 * (π / 2 - d.tht2) + π + d.tht2

   elseif p.reg == 9

      xi2 * (π / 2 - d.tht2) - π / 2

   elseif p.reg == 10

      xi2 * (d.tht2 - d.tht1) - d.tht2

   else

      throw(ArgumentError("gamx is only defined for boundary regions 1:10; got reg=$(p.reg)"))

   end

   st, ct = sincos(theta)

   f = ct^2 - st^2
   hP = sqrt(f + sqrt(d.P + f^2))

   return d.A + d.R * hP * ct

end

function gamx!(out::StridedArray{Float64}, d::peanut,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   if p.reg == 1

      ak = 2 * d.tht1
      bk = -d.tht1

   elseif p.reg == 2

      ak = d.tht2 - d.tht1
      bk = d.tht1

   elseif p.reg == 3

      ak = π / 2 - d.tht2
      bk = d.tht2

   elseif p.reg == 4

      ak = π / 2 - d.tht2
      bk = π / 2

   elseif p.reg == 5

      ak = d.tht2 - d.tht1
      bk = π - d.tht2

   elseif p.reg == 6

      ak = 2 * d.tht1
      bk = π - d.tht1

   elseif p.reg == 7

      ak = d.tht2 - d.tht1
      bk = π + d.tht1

   elseif p.reg == 8

      ak = π / 2 - d.tht2
      bk = π + d.tht2

   elseif p.reg == 9

      ak = π / 2 - d.tht2
      bk = -π / 2

   elseif p.reg == 10

      ak = d.tht2 - d.tht1
      bk = -d.tht2

   else

      throw(ArgumentError("gamx! is only defined for boundary regions 1:10; got reg=$(p.reg)"))

   end

   @inbounds for I in eachindex(t, out)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(ak, xi2, bk)

      st, ct = sincos(theta)

      f = ct^2 - st^2
      hP = sqrt(f + sqrt(d.P + f^2))

      out[I] = muladd(d.R * hP, ct, d.A)

   end

   return nothing

end

"""
    gamy(d::peanut, t::Float64, k::Int)
    gamy!(out::StridedArray{Float64}, d::peanut, t::StridedArray{Float64}, k::Int)

Second y-coordinate of the boundary parametrization gamma(t).

- `gamy(d, t, k)` computes the right boundary of patch `k`, where `t in [-1,1]`.
- Only defined for boundary regions `1:10`.

`t` may be a scalar or an array; the return/output has the same shape.
"""
function gamy(d::peanut, t::Float64, k::Int)::Float64

   p = d.pths[k]

   xi2 = (p.tk1 - p.tk0) * t / 2 + (p.tk0 + p.tk1) / 2

   theta = if p.reg == 1

      (2 * xi2 - 1) * d.tht1

   elseif p.reg == 2

      xi2 * (d.tht2 - d.tht1) + d.tht1

   elseif p.reg == 3

      xi2 * (π / 2 - d.tht2) + d.tht2

   elseif p.reg == 4

      xi2 * (π / 2 - d.tht2) + π / 2

   elseif p.reg == 5

      xi2 * (d.tht2 - d.tht1) + π - d.tht2

   elseif p.reg == 6

      2 * xi2 * d.tht1 + π - d.tht1

   elseif p.reg == 7

      xi2 * (d.tht2 - d.tht1) + π + d.tht1

   elseif p.reg == 8

      xi2 * (π / 2 - d.tht2) + π + d.tht2

   elseif p.reg == 9

      xi2 * (π / 2 - d.tht2) - π / 2

   elseif p.reg == 10

      xi2 * (d.tht2 - d.tht1) - d.tht2

   else

      throw(ArgumentError("gamy is only defined for boundary regions 1:10; got reg=$(p.reg)"))

   end

   st, ct = sincos(theta)

   f = ct^2 - st^2
   hP = sqrt(f + sqrt(d.P + f^2))

   return d.B + d.R * hP * st

end

function gamy!(out::StridedArray{Float64}, d::peanut,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   if p.reg == 1

      ak = 2 * d.tht1
      bk = -d.tht1

   elseif p.reg == 2

      ak = d.tht2 - d.tht1
      bk = d.tht1

   elseif p.reg == 3

      ak = π / 2 - d.tht2
      bk = d.tht2

   elseif p.reg == 4

      ak = π / 2 - d.tht2
      bk = π / 2

   elseif p.reg == 5

      ak = d.tht2 - d.tht1
      bk = π - d.tht2

   elseif p.reg == 6

      ak = 2 * d.tht1
      bk = π - d.tht1

   elseif p.reg == 7

      ak = d.tht2 - d.tht1
      bk = π + d.tht1

   elseif p.reg == 8

      ak = π / 2 - d.tht2
      bk = π + d.tht2

   elseif p.reg == 9

      ak = π / 2 - d.tht2
      bk = -π / 2

   elseif p.reg == 10

      ak = d.tht2 - d.tht1
      bk = -d.tht2

   else

      throw(ArgumentError("gamy! is only defined for boundary regions 1:10; got reg=$(p.reg)"))

   end

   @inbounds for I in eachindex(t, out)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(ak, xi2, bk)

      st, ct = sincos(theta)

      f = ct^2 - st^2
      hP = sqrt(f + sqrt(d.P + f^2))

      out[I] = muladd(d.R * hP, st, d.B)

   end

   return nothing

end

function gam(d::peanut, t::Float64, k::Int)::Tuple{Float64, Float64}

   p = d.pths[k]

   xi2 = (p.tk1 - p.tk0) * t / 2 + (p.tk0 + p.tk1) / 2

   theta = if p.reg == 1

      (2 * xi2 - 1) * d.tht1

   elseif p.reg == 2

      xi2 * (d.tht2 - d.tht1) + d.tht1

   elseif p.reg == 3

      xi2 * (π / 2 - d.tht2) + d.tht2

   elseif p.reg == 4

      xi2 * (π / 2 - d.tht2) + π / 2

   elseif p.reg == 5

      xi2 * (d.tht2 - d.tht1) + π - d.tht2

   elseif p.reg == 6

      2 * xi2 * d.tht1 + π - d.tht1

   elseif p.reg == 7

      xi2 * (d.tht2 - d.tht1) + π + d.tht1

   elseif p.reg == 8

      xi2 * (π / 2 - d.tht2) + π + d.tht2

   elseif p.reg == 9

      xi2 * (π / 2 - d.tht2) - π / 2

   elseif p.reg == 10

      xi2 * (d.tht2 - d.tht1) - d.tht2

   else

      throw(ArgumentError("gam is only defined for boundary regions 1:10; got reg=$(p.reg)"))

   end

   st, ct = sincos(theta)

   f = ct^2 - st^2
   hP = sqrt(f + sqrt(d.P + f^2))

   outx = muladd(d.R * hP, ct, d.A)
   outy = muladd(d.R * hP, st, d.B)

   return outx, outy

end

function gam!(out::Vector{Float64}, d::peanut, t::Float64, k::Int)

   x, y = gam(d, t, k)

   out[1] = x
   out[2] = y

   return nothing

end

function gam!(out::Matrix{Float64}, d::peanut, t::Vector{Float64}, k::Int)

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   if p.reg == 1

      ak = 2 * d.tht1
      bk = -d.tht1

   elseif p.reg == 2

      ak = d.tht2 - d.tht1
      bk = d.tht1

   elseif p.reg == 3

      ak = π / 2 - d.tht2
      bk = d.tht2

   elseif p.reg == 4

      ak = π / 2 - d.tht2
      bk = π / 2

   elseif p.reg == 5

      ak = d.tht2 - d.tht1
      bk = π - d.tht2

   elseif p.reg == 6

      ak = 2 * d.tht1
      bk = π - d.tht1

   elseif p.reg == 7

      ak = d.tht2 - d.tht1
      bk = π + d.tht1

   elseif p.reg == 8

      ak = π / 2 - d.tht2
      bk = π + d.tht2

   elseif p.reg == 9

      ak = π / 2 - d.tht2
      bk = -π / 2

   elseif p.reg == 10

      ak = d.tht2 - d.tht1
      bk = -d.tht2

   else

      throw(ArgumentError("gam! is only defined for boundary regions 1:10; got reg=$(p.reg)"))

   end

   @inbounds for I in eachindex(t)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(ak, xi2, bk)

      st, ct = sincos(theta)

      f = ct^2 - st^2
      hP = sqrt(f + sqrt(d.P + f^2))

      out[1, I] = muladd(d.R * hP, ct, d.A)
      out[2, I] = muladd(d.R * hP, st, d.B)

   end

   return nothing

end

function drawbd(d::peanut, flag=true; L::Int=33, show::Bool=true)

   colors = (
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
      RGBf(0, 0, 1),                 # blue
      RGBf(0, 0, 1),                 # blue
   )

   return drawbd_geom(d, colors; flag=flag, L=L, show=show)
end

"""
   dgamx(d::peanut, t::Float64, k::Int) -> Float64
   dgamx!(out::StridedArray{Float64}, d::peanut, t::StridedArray{Float64}, k::Int)

Derivative of the first coordinate of the boundary parametrization.

- With `k`: derivative of the right boundary of patch `k`,
  i.e. derivative of `gamx(d, t, k)` with respect to `t in [-1,1]`.
"""
function dgamx(d::peanut, t::Float64, k::Int)::Float64

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   xi2 = muladd(αt, t, βt)

   ak, bk = if p.reg == 1
      2 * d.tht1, -d.tht1
   elseif p.reg == 2
      d.tht2 - d.tht1, d.tht1
   elseif p.reg == 3
      π / 2 - d.tht2, d.tht2
   elseif p.reg == 4
      π / 2 - d.tht2, π / 2
   elseif p.reg == 5
      d.tht2 - d.tht1, π - d.tht2
   elseif p.reg == 6
      2 * d.tht1, π - d.tht1
   elseif p.reg == 7
      d.tht2 - d.tht1, π + d.tht1
   elseif p.reg == 8
      π / 2 - d.tht2, π + d.tht2
   elseif p.reg == 9
      π / 2 - d.tht2, -π / 2
   elseif p.reg == 10
      d.tht2 - d.tht1, -d.tht2
   else
      throw(ArgumentError("dgamx is only defined for boundary regions 1:10; got reg=$(p.reg)"))
   end

   theta = muladd(ak, xi2, bk)

   st, ct = sincos(theta)

   f = cos(2 * theta)
   q = sqrt(d.P + f^2)
   h = sqrt(f + q)

   dh = -2.0 * st * ct * h / q

   dgx = d.R * (dh * ct - h * st)

   return αt * ak * dgx

end

function dgamx!(out::StridedArray{Float64}, d::peanut,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   ak, bk = if p.reg == 1
      2 * d.tht1, -d.tht1
   elseif p.reg == 2
      d.tht2 - d.tht1, d.tht1
   elseif p.reg == 3
      π / 2 - d.tht2, d.tht2
   elseif p.reg == 4
      π / 2 - d.tht2, π / 2
   elseif p.reg == 5
      d.tht2 - d.tht1, π - d.tht2
   elseif p.reg == 6
      2 * d.tht1, π - d.tht1
   elseif p.reg == 7
      d.tht2 - d.tht1, π + d.tht1
   elseif p.reg == 8
      π / 2 - d.tht2, π + d.tht2
   elseif p.reg == 9
      π / 2 - d.tht2, -π / 2
   elseif p.reg == 10
      d.tht2 - d.tht1, -d.tht2
   else
      throw(ArgumentError("dgamx! is only defined for boundary regions 1:10; got reg=$(p.reg)"))
   end

   scale = αt * ak * d.R

   @inbounds for I in eachindex(t, out)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(ak, xi2, bk)

      st, ct = sincos(theta)

      f = cos(2 * theta)
      q = sqrt(d.P + f^2)
      h = sqrt(f + q)

      dh = -2.0 * st * ct * h /  q

      out[I] = scale * (dh * ct - h * st)

   end

   return nothing

end

"""
   dgamy(d::peanut, t::Float64, k::Int) -> Float64
   dgamy!(out::StridedArray{Float64}, d::peanut, t::StridedArray{Float64}, k::Int)

Derivative of the second coordinate of the boundary parametrization.
"""
function dgamy(d::peanut, t::Float64, k::Int)::Float64

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   xi2 = muladd(αt, t, βt)

   ak, bk = if p.reg == 1
      2 * d.tht1, -d.tht1
   elseif p.reg == 2
      d.tht2 - d.tht1, d.tht1
   elseif p.reg == 3
      π / 2 - d.tht2, d.tht2
   elseif p.reg == 4
      π / 2 - d.tht2, π / 2
   elseif p.reg == 5
      d.tht2 - d.tht1, π - d.tht2
   elseif p.reg == 6
      2 * d.tht1, π - d.tht1
   elseif p.reg == 7
      d.tht2 - d.tht1, π + d.tht1
   elseif p.reg == 8
      π / 2 - d.tht2, π + d.tht2
   elseif p.reg == 9
      π / 2 - d.tht2, -π / 2
   elseif p.reg == 10
      d.tht2 - d.tht1, -d.tht2
   else
      throw(ArgumentError("dgamy is only defined for boundary regions 1:10; got reg=$(p.reg)"))
   end

   theta = muladd(ak, xi2, bk)

   st, ct = sincos(theta)

   f = cos(2 * theta)
   q = sqrt(d.P + f^2)
   h = sqrt(f + q)

   dh = -2.0 * st * ct * h / q

   dgy = d.R * (dh * st + h * ct)

   return αt * ak * dgy

end

function dgamy!(out::StridedArray{Float64}, d::peanut,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   ak, bk = if p.reg == 1
      2 * d.tht1, -d.tht1
   elseif p.reg == 2
      d.tht2 - d.tht1, d.tht1
   elseif p.reg == 3
      π / 2 - d.tht2, d.tht2
   elseif p.reg == 4
      π / 2 - d.tht2, π / 2
   elseif p.reg == 5
      d.tht2 - d.tht1, π - d.tht2
   elseif p.reg == 6
      2 * d.tht1, π - d.tht1
   elseif p.reg == 7
      d.tht2 - d.tht1, π + d.tht1
   elseif p.reg == 8
      π / 2 - d.tht2, π + d.tht2
   elseif p.reg == 9
      π / 2 - d.tht2, -π / 2
   elseif p.reg == 10
      d.tht2 - d.tht1, -d.tht2
   else
      throw(ArgumentError("dgamy! is only defined for boundary regions 1:10; got reg=$(p.reg)"))
   end

   scale = αt * ak * d.R

   @inbounds for I in eachindex(t, out)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(ak, xi2, bk)

      st, ct = sincos(theta)

      f = cos(2 * theta)
      q = sqrt(d.P + f^2)
      h = sqrt(f + q)

      dh = -2.0 * st * ct * h / q

      out[I] = scale * (dh * st + h * ct)

   end

   return nothing

end

function dgam(d::peanut, t::Float64, k::Int)::Tuple{Float64, Float64}

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   xi2 = muladd(αt, t, βt)

   ak, bk = if p.reg == 1
      2 * d.tht1, -d.tht1
   elseif p.reg == 2
      d.tht2 - d.tht1, d.tht1
   elseif p.reg == 3
      π / 2 - d.tht2, d.tht2
   elseif p.reg == 4
      π / 2 - d.tht2, π / 2
   elseif p.reg == 5
      d.tht2 - d.tht1, π - d.tht2
   elseif p.reg == 6
      2 * d.tht1, π - d.tht1
   elseif p.reg == 7
      d.tht2 - d.tht1, π + d.tht1
   elseif p.reg == 8
      π / 2 - d.tht2, π + d.tht2
   elseif p.reg == 9
      π / 2 - d.tht2, -π / 2
   elseif p.reg == 10
      d.tht2 - d.tht1, -d.tht2
   else
      throw(ArgumentError("dgam is only defined for boundary regions 1:10; got reg=$(p.reg)"))
   end

   theta = muladd(ak, xi2, bk)

   st, ct = sincos(theta)

   f = cos(2 * theta)
   q = sqrt(d.P + f^2)
   h = sqrt(f + q)

   dh = -2.0 * st * ct * h / q

   scale = αt * ak * d.R

   outx = scale * (dh * ct - h * st)
   outy = scale * (dh * st + h * ct)

   return outx, outy

end

function dgam!(out::Vector{Float64}, d::peanut, t::Float64, k::Int)

   dx, dy = dgam(d, t, k)

   out[1] = dx
   out[2] = dy

   return nothing

end

function dgam!(out::Matrix{Float64}, d::peanut, t::Vector{Float64}, k::Int)

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   ak, bk = if p.reg == 1
      2 * d.tht1, -d.tht1
   elseif p.reg == 2
      d.tht2 - d.tht1, d.tht1
   elseif p.reg == 3
      π / 2 - d.tht2, d.tht2
   elseif p.reg == 4
      π / 2 - d.tht2, π / 2
   elseif p.reg == 5
      d.tht2 - d.tht1, π - d.tht2
   elseif p.reg == 6
      2 * d.tht1, π - d.tht1
   elseif p.reg == 7
      d.tht2 - d.tht1, π + d.tht1
   elseif p.reg == 8
      π / 2 - d.tht2, π + d.tht2
   elseif p.reg == 9
      π / 2 - d.tht2, -π / 2
   elseif p.reg == 10
      d.tht2 - d.tht1, -d.tht2
   else
      throw(ArgumentError("dgam! is only defined for boundary regions 1:10; got reg=$(p.reg)"))
   end

   scale = αt * ak * d.R

   @inbounds for I in eachindex(t)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(ak, xi2, bk)

      st, ct = sincos(theta)

      f = cos(2 * theta)
      q = sqrt(d.P + f^2)
      h = sqrt(f + q)

      dh = -2.0 * st * ct * h / q

      out[1, I] = scale * (dh * ct - h * st)
      out[2, I] = scale * (dh * st + h * ct)

   end

   return nothing

end

function gamp(d::peanut, t::Float64, k::Int)::Tuple{Float64, Float64}

   dx, dy = dgam(d, t, k)

   return dy, -dx

end

function gamp!(out::Matrix{Float64}, d::peanut, t::Vector{Float64}, k::Int)

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   ak, bk = if p.reg == 1
      2 * d.tht1, -d.tht1
   elseif p.reg == 2
      d.tht2 - d.tht1, d.tht1
   elseif p.reg == 3
      π / 2 - d.tht2, d.tht2
   elseif p.reg == 4
      π / 2 - d.tht2, π / 2
   elseif p.reg == 5
      d.tht2 - d.tht1, π - d.tht2
   elseif p.reg == 6
      2 * d.tht1, π - d.tht1
   elseif p.reg == 7
      d.tht2 - d.tht1, π + d.tht1
   elseif p.reg == 8
      π / 2 - d.tht2, π + d.tht2
   elseif p.reg == 9
      π / 2 - d.tht2, -π / 2
   elseif p.reg == 10
      d.tht2 - d.tht1, -d.tht2
   else
      throw(ArgumentError("gamp! is only defined for boundary regions 1:10; got reg=$(p.reg)"))
   end

   scale = αt * ak * d.R

   @inbounds for I in eachindex(t)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(ak, xi2, bk)

      st, ct = sincos(theta)

      f = cos(2 * theta)
      q = sqrt(d.P + f^2)
      h = sqrt(f + q)

      dh = -2.0 * st * ct * h / q

      dx = scale * (dh * ct - h * st)
      dy = scale * (dh * st + h * ct)

      out[1, I] = dy
      out[2, I] = -dx

   end

   return nothing

end

function nu!(out::Matrix{Float64}, d::peanut, t::Vector{Float64}, k::Int)

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   ak, bk = if p.reg == 1
      2 * d.tht1, -d.tht1
   elseif p.reg == 2
      d.tht2 - d.tht1, d.tht1
   elseif p.reg == 3
      π / 2 - d.tht2, d.tht2
   elseif p.reg == 4
      π / 2 - d.tht2, π / 2
   elseif p.reg == 5
      d.tht2 - d.tht1, π - d.tht2
   elseif p.reg == 6
      2 * d.tht1, π - d.tht1
   elseif p.reg == 7
      d.tht2 - d.tht1, π + d.tht1
   elseif p.reg == 8
      π / 2 - d.tht2, π + d.tht2
   elseif p.reg == 9
      π / 2 - d.tht2, -π / 2
   elseif p.reg == 10
      d.tht2 - d.tht1, -d.tht2
   else
      throw(ArgumentError("nu! is only defined for boundary regions 1:10; got reg=$(p.reg)"))
   end

   scale = αt * ak * d.R

   @inbounds for I in eachindex(t)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(ak, xi2, bk)

      st, ct = sincos(theta)

      f = cos(2 * theta)
      q = sqrt(d.P + f^2)
      h = sqrt(f + q)

      dh = -2.0 * st * ct * h / q

      dx = scale * (dh * ct - h * st)
      dy = scale * (dh * st + h * ct)

      S = hypot(dy, dx)

      out[1, I] = dy / S
      out[2, I] = -dx / S

   end

   return nothing

end

# This function finds s such that γ_l(t) = γ_k(s),
# allowing s to lie outside [-1, 1].
function bdinv(d::peanut, t::Float64, l::Int, k::Int)::Float64

   @inbounds begin

      pl = d.pths[l]
      pk = d.pths[k]

      plr = pl.reg
      pkr = pk.reg

      # Map local t on patch l to its xi coordinate.
      xil = muladd(0.5 * (pl.tk1 - pl.tk0), t, 0.5 * (pl.tk0 + pl.tk1))

      # If both patches are in the same angular region, the same xi works.
      if plr == pkr
         return xi_inv(xil, pk.tk0, pk.tk1)
      end

      # Compute physical angle theta corresponding to γ_l(t).
      th = if plr == 1
         xil * 2 * d.tht1 - d.tht1
      elseif plr == 2
         xil * (d.tht2 - d.tht1) + d.tht1
      elseif plr == 3
         xil * (π / 2 - d.tht2) + d.tht2
      elseif plr == 4
         xil * (π / 2 - d.tht2) + π / 2
      elseif plr == 5
         xil * (d.tht2 - d.tht1) + π - d.tht2
      elseif plr == 6
         xil * 2 * d.tht1 + π - d.tht1
      elseif plr == 7
         xil * (d.tht2 - d.tht1) + π + d.tht1
      elseif plr == 8
         xil * (π / 2 - d.tht2) + π + d.tht2
      elseif plr == 9
         xil * (π / 2 - d.tht2) - π / 2
      elseif plr == 10
         xil * (d.tht2 - d.tht1) - d.tht2
      else
         throw(ArgumentError("bdinv source region must be 1–10; got reg=$plr"))
      end

      # Get angular map for target patch k:
      #     theta = ak * xi + bk
      ak, bk = if pkr == 1
         2 * d.tht1, -d.tht1
      elseif pkr == 2
         d.tht2 - d.tht1, d.tht1
      elseif pkr == 3
         π / 2 - d.tht2, d.tht2
      elseif pkr == 4
         π / 2 - d.tht2, π / 2
      elseif pkr == 5
         d.tht2 - d.tht1, π - d.tht2
      elseif pkr == 6
         2 * d.tht1, π - d.tht1
      elseif pkr == 7
         d.tht2 - d.tht1, π + d.tht1
      elseif pkr == 8
         π / 2 - d.tht2, π + d.tht2
      elseif pkr == 9
         π / 2 - d.tht2, -π / 2
      elseif pkr == 10
         d.tht2 - d.tht1, -d.tht2
      else
         throw(ArgumentError("bdinv target region must be 1–10; got reg=$pkr"))
      end

      # Shift th by multiples of 2π so it is closest to the target region.
      # The target angular interval is [bk, bk + ak].
      center = bk + 0.5 * ak
      m = round((center - th) / (2π))
      ths = th + 2π * m

      # Convert angle to target xi, then xi to local patch coordinate s.
      xik = (ths - bk) / ak

      return xi_inv(xik, pk.tk0, pk.tk1)

   end

end

"""
   hderhigher(d::peanut, θ::Float64)

Compute h, h', h'', h⁽³⁾, h⁽⁴⁾, h⁽⁵⁾, h⁽⁶⁾ for the peanut radial function

   h(θ) = sqrt(cos(2θ) + sqrt(P + cos(2θ)^2))

where P = d.P.

Returns

   h, h1, h2, h3, h4, h5, h6

with derivatives taken with respect to θ.
"""
function hderhigher(d::peanut, θ::Float64)

   P = d.P

   # f = cos(2θ)
   f = cos(2θ)
   f1 = -2.0 * sin(2θ)
   f2 = -4.0 * f
   f3 = -4.0 * f1
   f4 = 16.0 * f
   f5 = 16.0 * f1
   f6 = -64.0 * f

   # l = P + f^2
   l = P + f^2

   l1 = 2.0 * f * f1

   l2 = 2.0 * f1^2 + 2.0 * f * f2

   l3 = 6.0 * f1 * f2 + 2.0 * f * f3

   l4 = 6.0 * f2^2 + 8.0 * f1 * f3 + 2.0 * f * f4

   l5 = 20.0 * f2 * f3 + 10.0 * f1 * f4 + 2.0 * f * f5

   l6 = 20.0 * f3^2 +
        30.0 * f2 * f4 +
        12.0 * f1 * f5 +
        2.0 * f * f6

   # q = sqrt(l)
   q = sqrt(l)

   q1 = l1 / (2.0 * q)

   q2 = (l2 - 2.0 * q1^2) / (2.0 * q)

   q3 = (l3 - 6.0 * q1 * q2) / (2.0 * q)

   q4 = (l4 - 6.0 * q2^2 - 8.0 * q1 * q3) / (2.0 * q)

   q5 = (l5 - 10.0 * q1 * q4 - 20.0 * q2 * q3) / (2.0 * q)

   q6 = (l6 -
         12.0 * q1 * q5 -
         30.0 * q2 * q4 -
         20.0 * q3^2) / (2.0 * q)

   # h = sqrt(f + q)
   m = f + q
   m1 = f1 + q1
   m2 = f2 + q2
   m3 = f3 + q3
   m4 = f4 + q4
   m5 = f5 + q5
   m6 = f6 + q6

   h = sqrt(m)

   h1 = m1 / (2.0 * h)

   h2 = (m2 - 2.0 * h1^2) / (2.0 * h)

   h3 = (m3 - 6.0 * h1 * h2) / (2.0 * h)

   h4 = (m4 - 6.0 * h2^2 - 8.0 * h1 * h3) / (2.0 * h)

   h5 = (m5 - 10.0 * h1 * h4 - 20.0 * h2 * h3) / (2.0 * h)

   h6 = (m6 -
         12.0 * h1 * h5 -
         30.0 * h2 * h4 -
         20.0 * h3^2) / (2.0 * h)

   return h, h1, h2, h3, h4, h5, h6

end

"""
    DLP!(out::StridedArray{Float64}, d::peanut, t::Float64,
         tau::StridedArray{Float64}, k::Int,
         x::Vector{Float64}, G::Vector{Float64}, GP::Vector{Float64})

Same-panel double-layer kernel on the peanut boundary:

    K(t, τ) =
    ((γ_k(τ) - γ_k(t)) ⋅ γᵖᵉʳᵖ_k(τ)) /
    ||γ_k(τ) - γ_k(t)||²

Uses the analytically simplified polar form where the removable sin²(Δ)
singularity is canceled over the whole panel.

The Taylor expansions for q and q₂ are used when Δ is small.
"""
function DLP!(out::StridedArray{Float64}, d::peanut, t::Float64,
   tau::StridedArray{Float64}, k::Int, x::Vector{Float64},
   G::Vector{Float64}, GP::Vector{Float64})

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = p.tk0 + αt

   reg = p.reg

   ak, bk = if reg == 1
      2 * d.tht1, -d.tht1
   elseif reg == 2
      d.tht2 - d.tht1, d.tht1
   elseif reg == 3
      π / 2 - d.tht2, d.tht2
   elseif reg == 4
      π / 2 - d.tht2, π / 2
   elseif reg == 5
      d.tht2 - d.tht1, π - d.tht2
   elseif reg == 6
      2 * d.tht1, π - d.tht1
   elseif reg == 7
      d.tht2 - d.tht1, π + d.tht1
   elseif reg == 8
      π / 2 - d.tht2, π + d.tht2
   elseif reg == 9
      π / 2 - d.tht2, -π / 2
   elseif reg == 10
      d.tht2 - d.tht1, -d.tht2
   else
      throw(ArgumentError("DLP! self term defined only for boundary regions 1–10; got reg=$reg"))
   end

   # t mapped from [-1,1] to [tk0,tk1], then to theta.
   tm = muladd(αt, t, βt)
   tht = muladd(ak, tm, bk)

   # h and derivatives at theta_t.
   ht0, ht1, ht2, ht3, ht4, ht5, ht6 = hderhigher(d, tht)

   # dtheta / dtau
   λ = αt * ak

   @inbounds for i in eachindex(tau, out)

      τm = muladd(αt, tau[i], βt)
      thτ = muladd(ak, τm, bk)

      # hτ and hτ'
      fτ = cos(2 * thτ)
      qτ = sqrt(d.P + fτ^2)
      hτ = sqrt(fτ + qτ)

      fτ1 = -2.0 * sin(2 * thτ)
      hτ1 = fτ1 * hτ / (2.0 * qτ)

      Δ = 0.5 * (thτ - tht)

      sΔ, cΔ = sincos(Δ)
      # (2e-3)^5 = 3.2e-14
      if abs(Δ) < 2e-3

         Δ2 = Δ * Δ
         Δ3 = Δ2 * Δ
         Δ4 = Δ2 * Δ2

         q = ht1 + Δ * ht2 +
             Δ2 * ((2 / 3) * ht3 + ht1 / 6) +
             Δ3 * (ht4 / 3 + ht2 / 6) +
             Δ4 * ((2 / 15) * ht5 + ht3 / 9 + (7 / 360) * ht1)

         q2 = -ht2 +
              Δ * (-(4 / 3) * ht3 + (2 / 3) * ht1) +
              Δ2 * (-ht4 + ht2) +
              Δ3 * (-(8 / 15) * ht5 + (8 / 9) * ht3 + (4 / 45) * ht1) +
              Δ4 * (-(2 / 9) * ht6 + (5 / 9) * ht4 + ht2 / 9)

      else

         q = (hτ - ht0) / (2.0 * sΔ)
         q2 = (q - hτ1 * cΔ) / sΔ

      end

      num = ht0 * hτ + hτ * q2 + 2.0 * q * hτ1 * cΔ
      den = q^2 + hτ * ht0

      out[i] = 0.5 * λ * num / den

   end

   return nothing

end
#-----------------------

"""
   gamder(d::peanut, θ::Float64)

Return centered peanut boundary values and first derivatives with respect to `θ`.

Returns gx, gy, dgx, dgy

where gx(θ) = R*h(θ)*cos(θ), gy(θ) = R*h(θ)*sin(θ)

and h(θ) = sqrt(cos(2θ) + sqrt(P + cos(2θ)^2)).
"""
function gamder(d::peanut, θ::Float64)::Tuple{Float64,Float64,Float64,Float64}

   st, ct = sincos(θ)

   f = cos(2 * θ)
   q = sqrt(d.P + f^2)
   h = sqrt(f + q)

   dh = -2.0 * st * ct * h / q

   gx = d.R * h * ct
   gy = d.R * h * st

   dgx = d.R * (dh * ct - h * st)
   dgy = d.R * (dh * st + h * ct)

   return gx, gy, dgx, dgy

end

"""
    Dmap!(DJ::StridedArray{Float64}, d::peanut,
          u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)

Fill `DJ` with the absolute determinant of the Jacobian of the patch map.

No allocations. `u`, `v`, and `DJ` must have compatible indexing.
"""
function Dmap!(DJ::StridedArray{Float64}, d::peanut,
   u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   αc = hc / 2
   αt = ht / 2
   βc = p.ck0 + αc
   βt = p.tk0 + αt

   r = p.reg

   if r == 11 || r == 12

      fill!(DJ, d.L1 * d.L2 * hc * ht / 4)

      return nothing

   end

   Rα = d.R * d.alpha
   L1h = d.L1 / 2
   L2h = d.L2 / 2

   if r == 1

      ak = 2 * d.tht1
      bk = -d.tht1

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = Rα + L2h
         Xy = d.L1 * (2vhat - 1) / 2

         dXx = 0.0
         dXy = d.L1

         th = muladd(ak, vhat, bk)

         st, ct = sincos(th)
         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)

         f1 = -2.0 * sin(2 * th)
         dh = f1 * h / (2.0 * q)

         gx = d.R * h * ct
         gy = d.R * h * st
         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * dXx + uhat * αt * ak * dgx
         dvy = (1 - uhat) * αt * dXy + uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 2

      ak = d.tht2 - d.tht1
      bk = d.tht1

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = Rα + d.L2 * (1 - 2vhat) / 2
         Xy = L1h

         dXx = -d.L2
         dXy = 0.0

         th = muladd(ak, vhat, bk)

         st, ct = sincos(th)
         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)

         f1 = -2.0 * sin(2 * th)
         dh = f1 * h / (2.0 * q)

         gx = d.R * h * ct
         gy = d.R * h * st
         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * dXx + uhat * αt * ak * dgx
         dvy = (1 - uhat) * αt * dXy + uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 3

      ak = π / 2 - d.tht2
      bk = d.tht2
      Cx = Rα - L2h

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = (1 - vhat) * Cx
         Xy = 0.0

         dXx = -Cx
         dXy = 0.0

         th = muladd(ak, vhat, bk)

         st, ct = sincos(th)
         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)

         f1 = -2.0 * sin(2 * th)
         dh = f1 * h / (2.0 * q)

         gx = d.R * h * ct
         gy = d.R * h * st
         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * dXx + uhat * αt * ak * dgx
         dvy = (1 - uhat) * αt * dXy + uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 4

      ak = π / 2 - d.tht2
      bk = π / 2
      Cx = L2h - Rα

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = vhat * Cx
         Xy = 0.0

         dXx = Cx
         dXy = 0.0

         th = muladd(ak, vhat, bk)

         st, ct = sincos(th)
         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)

         f1 = -2.0 * sin(2 * th)
         dh = f1 * h / (2.0 * q)

         gx = d.R * h * ct
         gy = d.R * h * st
         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * dXx + uhat * αt * ak * dgx
         dvy = (1 - uhat) * αt * dXy + uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 5

      ak = d.tht2 - d.tht1
      bk = π - d.tht2

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = -Rα + d.L2 * (1 - 2vhat) / 2
         Xy = L1h

         dXx = -d.L2
         dXy = 0.0

         th = muladd(ak, vhat, bk)

         st, ct = sincos(th)
         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)

         f1 = -2.0 * sin(2 * th)
         dh = f1 * h / (2.0 * q)

         gx = d.R * h * ct
         gy = d.R * h * st
         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * dXx + uhat * αt * ak * dgx
         dvy = (1 - uhat) * αt * dXy + uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 6

      ak = 2 * d.tht1
      bk = π - d.tht1

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = -Rα - L2h
         Xy = -d.L1 * (2vhat - 1) / 2

         dXx = 0.0
         dXy = -d.L1

         th = muladd(ak, vhat, bk)

         st, ct = sincos(th)
         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)

         f1 = -2.0 * sin(2 * th)
         dh = f1 * h / (2.0 * q)

         gx = d.R * h * ct
         gy = d.R * h * st
         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * dXx + uhat * αt * ak * dgx
         dvy = (1 - uhat) * αt * dXy + uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 7

      ak = d.tht2 - d.tht1
      bk = π + d.tht1

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = -Rα - d.L2 * (1 - 2vhat) / 2
         Xy = -L1h

         dXx = d.L2
         dXy = 0.0

         th = muladd(ak, vhat, bk)

         st, ct = sincos(th)
         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)

         f1 = -2.0 * sin(2 * th)
         dh = f1 * h / (2.0 * q)

         gx = d.R * h * ct
         gy = d.R * h * st
         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * dXx + uhat * αt * ak * dgx
         dvy = (1 - uhat) * αt * dXy + uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 8

      ak = π / 2 - d.tht2
      bk = π + d.tht2
      Cx = L2h - Rα

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = (1 - vhat) * Cx
         Xy = 0.0

         dXx = -Cx
         dXy = 0.0

         th = muladd(ak, vhat, bk)

         st, ct = sincos(th)
         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)

         f1 = -2.0 * sin(2 * th)
         dh = f1 * h / (2.0 * q)

         gx = d.R * h * ct
         gy = d.R * h * st
         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * dXx + uhat * αt * ak * dgx
         dvy = (1 - uhat) * αt * dXy + uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 9

      ak = π / 2 - d.tht2
      bk = -π / 2
      Cx = Rα - L2h

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = vhat * Cx
         Xy = 0.0

         dXx = Cx
         dXy = 0.0

         th = muladd(ak, vhat, bk)

         st, ct = sincos(th)
         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)

         f1 = -2.0 * sin(2 * th)
         dh = f1 * h / (2.0 * q)

         gx = d.R * h * ct
         gy = d.R * h * st
         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * dXx + uhat * αt * ak * dgx
         dvy = (1 - uhat) * αt * dXy + uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 10

      ak = d.tht2 - d.tht1
      bk = -d.tht2

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = Rα - d.L2 * (1 - 2vhat) / 2
         Xy = -L1h

         dXx = d.L2
         dXy = 0.0

         th = muladd(ak, vhat, bk)

         st, ct = sincos(th)
         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)

         f1 = -2.0 * sin(2 * th)
         dh = f1 * h / (2.0 * q)

         gx = d.R * h * ct
         gy = d.R * h * st
         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * dXx + uhat * αt * ak * dgx
         dvy = (1 - uhat) * αt * dXy + uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   else

      throw(ArgumentError("Dmap! for peanut expects region 1–12; got reg=$r"))

   end

   return nothing

end

"""
A combination of mapxy! and Dmap! function (No allocations!)

The purpose of this function is to reduce repeated computations related
to cosine, sine, hP, and hP'.
"""
function mapxy_Dmap!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
   DJ::StridedArray{Float64}, d::peanut, u::StridedArray{Float64},
   v::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   αc = hc / 2
   αt = ht / 2
   βc = p.ck0 + αc
   βt = p.tk0 + αt

   r = p.reg

   Rα = d.R * d.alpha
   L1h = d.L1 / 2
   L2h = d.L2 / 2

   if r == 11 || r == 12

      if d.L2 >= d.L1
         αu = αt
         βu = βt
         αv = αc
         βv = βc
      else
         αu = αc
         βu = βc
         αv = αt
         βv = βt
      end

      x0 = if r == 11
         d.A + Rα - L2h
      else
         d.A - Rα - L2h
      end

      y0 = d.B - L1h
      detJ = d.L1 * d.L2 * hc * ht / 4

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)
         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Zx[I] = muladd(d.L2, xi1, x0)
         Zy[I] = muladd(d.L1, xi2, y0)
         DJ[I] = detJ
      end

      return nothing
   end

   if r == 1

      ak = 2 * d.tht1
      bk = -d.tht1
      Xxc = Rα + L2h

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)
         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xyc = d.L1 * (2vhat - 1) / 2

         th = muladd(ak, vhat, bk)
         st, ct = sincos(th)

         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)
         dh = -2.0 * st * ct * h / q

         Yxc = d.R * h * ct
         Yyc = d.R * h * st

         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         Zx[I] = d.A + muladd(uhat, Yxc - Xxc, Xxc)
         Zy[I] = d.B + muladd(uhat, Yyc - Xyc, Xyc)

         dux = αc * (Yxc - Xxc)
         duy = αc * (Yyc - Xyc)

         dvx = uhat * αt * ak * dgx
         dvy = (1 - uhat) * αt * d.L1 + uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)
      end

   elseif r == 2

      ak = d.tht2 - d.tht1
      bk = d.tht1
      Xyc = L1h

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)
         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xxc = Rα + d.L2 * (1 - 2vhat) / 2

         th = muladd(ak, vhat, bk)
         st, ct = sincos(th)

         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)
         dh = -2.0 * st * ct * h / q

         Yxc = d.R * h * ct
         Yyc = d.R * h * st

         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         Zx[I] = d.A + muladd(uhat, Yxc - Xxc, Xxc)
         Zy[I] = d.B + muladd(uhat, Yyc - Xyc, Xyc)

         dux = αc * (Yxc - Xxc)
         duy = αc * (Yyc - Xyc)

         dvx = -(1 - uhat) * αt * d.L2 + uhat * αt * ak * dgx
         dvy = uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)
      end

   elseif r == 3

      ak = π / 2 - d.tht2
      bk = d.tht2
      Cx = Rα - L2h

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)
         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xxc = (1 - vhat) * Cx

         th = muladd(ak, vhat, bk)
         st, ct = sincos(th)

         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)
         dh = -2.0 * st * ct * h / q

         Yxc = d.R * h * ct
         Yyc = d.R * h * st

         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         Zx[I] = d.A + muladd(uhat, Yxc - Xxc, Xxc)
         Zy[I] = d.B + uhat * Yyc

         dux = αc * (Yxc - Xxc)
         duy = αc * Yyc

         dvx = -(1 - uhat) * αt * Cx + uhat * αt * ak * dgx
         dvy = uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)
      end

   elseif r == 4

      ak = π / 2 - d.tht2
      bk = π / 2
      Cx = L2h - Rα

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)
         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xxc = vhat * Cx

         th = muladd(ak, vhat, bk)
         st, ct = sincos(th)

         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)
         dh = -2.0 * st * ct * h / q

         Yxc = d.R * h * ct
         Yyc = d.R * h * st

         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         Zx[I] = d.A + muladd(uhat, Yxc - Xxc, Xxc)
         Zy[I] = d.B + uhat * Yyc

         dux = αc * (Yxc - Xxc)
         duy = αc * Yyc

         dvx = (1 - uhat) * αt * Cx + uhat * αt * ak * dgx
         dvy = uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)
      end

   elseif r == 5

      ak = d.tht2 - d.tht1
      bk = π - d.tht2
      Xyc = L1h

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)
         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xxc = -Rα + d.L2 * (1 - 2vhat) / 2

         th = muladd(ak, vhat, bk)
         st, ct = sincos(th)

         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)
         dh = -2.0 * st * ct * h / q

         Yxc = d.R * h * ct
         Yyc = d.R * h * st

         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         Zx[I] = d.A + muladd(uhat, Yxc - Xxc, Xxc)
         Zy[I] = d.B + muladd(uhat, Yyc - Xyc, Xyc)

         dux = αc * (Yxc - Xxc)
         duy = αc * (Yyc - Xyc)

         dvx = -(1 - uhat) * αt * d.L2 + uhat * αt * ak * dgx
         dvy = uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)
      end

   elseif r == 6

      ak = 2 * d.tht1
      bk = π - d.tht1
      Xxc = -Rα - L2h

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)
         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xyc = -d.L1 * (2vhat - 1) / 2

         th = muladd(ak, vhat, bk)
         st, ct = sincos(th)

         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)
         dh = -2.0 * st * ct * h / q

         Yxc = d.R * h * ct
         Yyc = d.R * h * st

         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         Zx[I] = d.A + muladd(uhat, Yxc - Xxc, Xxc)
         Zy[I] = d.B + muladd(uhat, Yyc - Xyc, Xyc)

         dux = αc * (Yxc - Xxc)
         duy = αc * (Yyc - Xyc)

         dvx = uhat * αt * ak * dgx
         dvy = -(1 - uhat) * αt * d.L1 + uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)
      end

   elseif r == 7

      ak = d.tht2 - d.tht1
      bk = π + d.tht1
      Xyc = -L1h

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)
         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xxc = -Rα - d.L2 * (1 - 2vhat) / 2

         th = muladd(ak, vhat, bk)
         st, ct = sincos(th)

         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)
         dh = -2.0 * st * ct * h / q

         Yxc = d.R * h * ct
         Yyc = d.R * h * st

         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         Zx[I] = d.A + muladd(uhat, Yxc - Xxc, Xxc)
         Zy[I] = d.B + muladd(uhat, Yyc - Xyc, Xyc)

         dux = αc * (Yxc - Xxc)
         duy = αc * (Yyc - Xyc)

         dvx = (1 - uhat) * αt * d.L2 + uhat * αt * ak * dgx
         dvy = uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)
      end

   elseif r == 8

      ak = π / 2 - d.tht2
      bk = π + d.tht2
      Cx = L2h - Rα

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)
         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xxc = (1 - vhat) * Cx

         th = muladd(ak, vhat, bk)
         st, ct = sincos(th)

         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)
         dh = -2.0 * st * ct * h / q

         Yxc = d.R * h * ct
         Yyc = d.R * h * st

         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         Zx[I] = d.A + muladd(uhat, Yxc - Xxc, Xxc)
         Zy[I] = d.B + uhat * Yyc

         dux = αc * (Yxc - Xxc)
         duy = αc * Yyc

         dvx = -(1 - uhat) * αt * Cx + uhat * αt * ak * dgx
         dvy = uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)
      end

   elseif r == 9

      ak = π / 2 - d.tht2
      bk = -π / 2
      Cx = Rα - L2h

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)
         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xxc = vhat * Cx

         th = muladd(ak, vhat, bk)
         st, ct = sincos(th)

         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)
         dh = -2.0 * st * ct * h / q

         Yxc = d.R * h * ct
         Yyc = d.R * h * st

         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         Zx[I] = d.A + muladd(uhat, Yxc - Xxc, Xxc)
         Zy[I] = d.B + uhat * Yyc

         dux = αc * (Yxc - Xxc)
         duy = αc * Yyc

         dvx = (1 - uhat) * αt * Cx + uhat * αt * ak * dgx
         dvy = uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)
      end

   elseif r == 10

      ak = d.tht2 - d.tht1
      bk = -d.tht2
      Xyc = -L1h

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)
         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xxc = Rα - d.L2 * (1 - 2vhat) / 2

         th = muladd(ak, vhat, bk)
         st, ct = sincos(th)

         f = cos(2 * th)
         q = sqrt(d.P + f^2)
         h = sqrt(f + q)
         dh = -2.0 * st * ct * h / q

         Yxc = d.R * h * ct
         Yyc = d.R * h * st

         dgx = d.R * (dh * ct - h * st)
         dgy = d.R * (dh * st + h * ct)

         Zx[I] = d.A + muladd(uhat, Yxc - Xxc, Xxc)
         Zy[I] = d.B + muladd(uhat, Yyc - Xyc, Xyc)

         dux = αc * (Yxc - Xxc)
         duy = αc * (Yyc - Xyc)

         dvx = (1 - uhat) * αt * d.L2 + uhat * αt * ak * dgx
         dvy = uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)
      end

   else

      throw(ArgumentError("mapxy_Dmap! for peanut expects region 1–12; got reg=$r"))

   end

   return nothing

end

function chk_map(d::peanut; n::Int=32, tol::Float64=5e-14)
   Iex = d.R^2 * (π * d.R^2 * (1 + d.P) / 2 + 2 * (d.A^2 + d.B^2) *
         sqrt(1 + d.P) * ellipe(1 / (1 + d.P)))
   return chkmap_geom(d, Iex; n=n, tol=tol)
end

"""
    jinvmap(d::peanut, u::Float64, v::Float64, r::Int)

Inverse Jacobian of the peanut region mapping at a single reference point
`(u,v) ∈ [-1,1]^2`.

Returns `(J11, J12, J21, J22)`, the entries of `J⁻¹`.

Only the region index `r` is needed because this is used on the region-level
map, so `hc = ht = 1`.
"""
function jinvmap(d::peanut, u::Float64, v::Float64, r::Int)

   uhat = (u + 1) / 2
   vhat = (v + 1) / 2

   Rα = d.R * d.alpha
   L1h = d.L1 / 2
   L2h = d.L2 / 2

   if r == 1

      ak = 2 * d.tht1
      th = muladd(ak, vhat, -d.tht1)

      st, ct = sincos(th)
      f = cos(2 * th)
      q = sqrt(d.P + f^2)
      h = sqrt(f + q)
      dh = -2.0 * st * ct * h / q

      gx = d.R * h * ct
      gy = d.R * h * st
      dgx = d.R * (dh * ct - h * st)
      dgy = d.R * (dh * st + h * ct)

      dux = (gx - Rα - L2h) / 2
      duy = (gy - d.L1 * (2vhat - 1) / 2) / 2

      dvx = uhat * ak * dgx / 2
      dvy = (1 - uhat) * d.L1 / 2 + uhat * ak * dgy / 2

   elseif r == 2

      ak = d.tht2 - d.tht1
      th = muladd(ak, vhat, d.tht1)

      st, ct = sincos(th)
      f = cos(2 * th)
      q = sqrt(d.P + f^2)
      h = sqrt(f + q)
      dh = -2.0 * st * ct * h / q

      gx = d.R * h * ct
      gy = d.R * h * st
      dgx = d.R * (dh * ct - h * st)
      dgy = d.R * (dh * st + h * ct)

      dux = (gx - Rα - d.L2 * (1 - 2vhat) / 2) / 2
      duy = (gy - L1h) / 2

      dvx = -(1 - uhat) * d.L2 / 2 + uhat * ak * dgx / 2
      dvy = uhat * ak * dgy / 2

   elseif r == 3

      ak = π / 2 - d.tht2
      Cx = Rα - L2h
      th = muladd(ak, vhat, d.tht2)

      st, ct = sincos(th)
      f = cos(2 * th)
      q = sqrt(d.P + f^2)
      h = sqrt(f + q)
      dh = -2.0 * st * ct * h / q

      gx = d.R * h * ct
      gy = d.R * h * st
      dgx = d.R * (dh * ct - h * st)
      dgy = d.R * (dh * st + h * ct)

      dux = (gx - (1 - vhat) * Cx) / 2
      duy = gy / 2

      dvx = -(1 - uhat) * Cx / 2 + uhat * ak * dgx / 2
      dvy = uhat * ak * dgy / 2

   elseif r == 4

      ak = π / 2 - d.tht2
      Cx = L2h - Rα
      th = muladd(ak, vhat, π / 2)

      st, ct = sincos(th)
      f = cos(2 * th)
      q = sqrt(d.P + f^2)
      h = sqrt(f + q)
      dh = -2.0 * st * ct * h / q

      gx = d.R * h * ct
      gy = d.R * h * st
      dgx = d.R * (dh * ct - h * st)
      dgy = d.R * (dh * st + h * ct)

      dux = (gx - vhat * Cx) / 2
      duy = gy / 2

      dvx = (1 - uhat) * Cx / 2 + uhat * ak * dgx / 2
      dvy = uhat * ak * dgy / 2

   elseif r == 5

      ak = d.tht2 - d.tht1
      th = muladd(ak, vhat, π - d.tht2)

      st, ct = sincos(th)
      f = cos(2 * th)
      q = sqrt(d.P + f^2)
      h = sqrt(f + q)
      dh = -2.0 * st * ct * h / q

      gx = d.R * h * ct
      gy = d.R * h * st
      dgx = d.R * (dh * ct - h * st)
      dgy = d.R * (dh * st + h * ct)

      dux = (gx + Rα - d.L2 * (1 - 2vhat) / 2) / 2
      duy = (gy - L1h) / 2

      dvx = -(1 - uhat) * d.L2 / 2 + uhat * ak * dgx / 2
      dvy = uhat * ak * dgy / 2

   elseif r == 6

      ak = 2 * d.tht1
      th = muladd(ak, vhat, π - d.tht1)

      st, ct = sincos(th)
      f = cos(2 * th)
      q = sqrt(d.P + f^2)
      h = sqrt(f + q)
      dh = -2.0 * st * ct * h / q

      gx = d.R * h * ct
      gy = d.R * h * st
      dgx = d.R * (dh * ct - h * st)
      dgy = d.R * (dh * st + h * ct)

      dux = (gx + Rα + L2h) / 2
      duy = (gy + d.L1 * (2vhat - 1) / 2) / 2

      dvx = uhat * ak * dgx / 2
      dvy = -(1 - uhat) * d.L1 / 2 + uhat * ak * dgy / 2

   elseif r == 7

      ak = d.tht2 - d.tht1
      th = muladd(ak, vhat, π + d.tht1)

      st, ct = sincos(th)
      f = cos(2 * th)
      q = sqrt(d.P + f^2)
      h = sqrt(f + q)
      dh = -2.0 * st * ct * h / q

      gx = d.R * h * ct
      gy = d.R * h * st
      dgx = d.R * (dh * ct - h * st)
      dgy = d.R * (dh * st + h * ct)

      dux = (gx + Rα + d.L2 * (1 - 2vhat) / 2) / 2
      duy = (gy + L1h) / 2

      dvx = (1 - uhat) * d.L2 / 2 + uhat * ak * dgx / 2
      dvy = uhat * ak * dgy / 2

   elseif r == 8

      ak = π / 2 - d.tht2
      Cx = L2h - Rα
      th = muladd(ak, vhat, π + d.tht2)

      st, ct = sincos(th)
      f = cos(2 * th)
      q = sqrt(d.P + f^2)
      h = sqrt(f + q)
      dh = -2.0 * st * ct * h / q

      gx = d.R * h * ct
      gy = d.R * h * st
      dgx = d.R * (dh * ct - h * st)
      dgy = d.R * (dh * st + h * ct)

      dux = (gx - (1 - vhat) * Cx) / 2
      duy = gy / 2

      dvx = -(1 - uhat) * Cx / 2 + uhat * ak * dgx / 2
      dvy = uhat * ak * dgy / 2

   elseif r == 9

      ak = π / 2 - d.tht2
      Cx = Rα - L2h
      th = muladd(ak, vhat, -π / 2)

      st, ct = sincos(th)
      f = cos(2 * th)
      q = sqrt(d.P + f^2)
      h = sqrt(f + q)
      dh = -2.0 * st * ct * h / q

      gx = d.R * h * ct
      gy = d.R * h * st
      dgx = d.R * (dh * ct - h * st)
      dgy = d.R * (dh * st + h * ct)

      dux = (gx - vhat * Cx) / 2
      duy = gy / 2

      dvx = (1 - uhat) * Cx / 2 + uhat * ak * dgx / 2
      dvy = uhat * ak * dgy / 2

   elseif r == 10

      ak = d.tht2 - d.tht1
      th = muladd(ak, vhat, -d.tht2)

      st, ct = sincos(th)
      f = cos(2 * th)
      q = sqrt(d.P + f^2)
      h = sqrt(f + q)
      dh = -2.0 * st * ct * h / q

      gx = d.R * h * ct
      gy = d.R * h * st
      dgx = d.R * (dh * ct - h * st)
      dgy = d.R * (dh * st + h * ct)

      dux = (gx - Rα + d.L2 * (1 - 2vhat) / 2) / 2
      duy = (gy + L1h) / 2

      dvx = (1 - uhat) * d.L2 / 2 + uhat * ak * dgx / 2
      dvy = uhat * ak * dgy / 2

   else

      throw(ArgumentError("jinvmap for peanut expects region 1–10; got reg=$r"))

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
  mapinv(d::peanut, u::Float64, v::Float64, k::Int) -> Tuple{Float64,Float64}

Inverse mapping: given a physical point `(u,v)` on patch `k`, return
the reference coordinates `[t; s]` in `[-1,1]^2` such that `mapxy(d, t, s, k) = (u,v)`.
"""
@inline function Xx(s::Float64, d::peanut, r::Int)

   Rα = d.R * d.alpha

   if r == 1
      return d.A + Rα + d.L2 / 2

   elseif r == 2
      return d.A + Rα + d.L2 * (1 - 2s) / 2

   elseif r == 3
      return d.A + (1 - s) * (Rα - d.L2 / 2)

   elseif r == 4
      return d.A + s * (d.L2 / 2 - Rα)

   elseif r == 5
      return d.A - Rα + d.L2 * (1 - 2s) / 2

   elseif r == 6
      return d.A - Rα - d.L2 / 2

   elseif r == 7
      return d.A - Rα - d.L2 * (1 - 2s) / 2

   elseif r == 8
      return d.A + (1 - s) * (d.L2 / 2 - Rα)

   elseif r == 9
      return d.A + s * (Rα - d.L2 / 2)

   elseif r == 10
      return d.A + Rα - d.L2 * (1 - 2s) / 2

   else
      throw(ArgumentError("Xx for peanut expects region 1–10; got r=$r"))
   end

end

@inline function Xy(s::Float64, d::peanut, r::Int)

   if r == 1
      return d.B + d.L1 * (2s - 1) / 2

   elseif r == 2
      return d.B + d.L1 / 2

   elseif r == 3
      return d.B

   elseif r == 4
      return d.B

   elseif r == 5
      return d.B + d.L1 / 2

   elseif r == 6
      return d.B - d.L1 * (2s - 1) / 2

   elseif r == 7
      return d.B - d.L1 / 2

   elseif r == 8
      return d.B

   elseif r == 9
      return d.B

   elseif r == 10
      return d.B - d.L1 / 2

   else
      throw(ArgumentError("Xy for peanut expects region 1–10; got r=$r"))
   end

end

@inline function Yx(s::Float64, d::peanut, r::Int)

   th = if r == 1
      (2s - 1) * d.tht1

   elseif r == 2
      s * (d.tht2 - d.tht1) + d.tht1

   elseif r == 3
      s * (π / 2 - d.tht2) + d.tht2

   elseif r == 4
      s * (π / 2 - d.tht2) + π / 2

   elseif r == 5
      s * (d.tht2 - d.tht1) + π - d.tht2

   elseif r == 6
      2s * d.tht1 + π - d.tht1

   elseif r == 7
      s * (d.tht2 - d.tht1) + π + d.tht1

   elseif r == 8
      s * (π / 2 - d.tht2) + π + d.tht2

   elseif r == 9
      s * (π / 2 - d.tht2) - π / 2

   elseif r == 10
      s * (d.tht2 - d.tht1) - d.tht2

   else
      throw(ArgumentError("Yx for peanut expects region 1–10; got r=$r"))
   end

   st, ct = sincos(th)

   f = cos(2 * th)
   hP = sqrt(f + sqrt(d.P + f^2))

   return d.A + d.R * hP * ct

end

@inline function Yy(s::Float64, d::peanut, r::Int)

   th = if r == 1
      (2s - 1) * d.tht1

   elseif r == 2
      s * (d.tht2 - d.tht1) + d.tht1

   elseif r == 3
      s * (π / 2 - d.tht2) + d.tht2

   elseif r == 4
      s * (π / 2 - d.tht2) + π / 2

   elseif r == 5
      s * (d.tht2 - d.tht1) + π - d.tht2

   elseif r == 6
      2s * d.tht1 + π - d.tht1

   elseif r == 7
      s * (d.tht2 - d.tht1) + π + d.tht1

   elseif r == 8
      s * (π / 2 - d.tht2) + π + d.tht2

   elseif r == 9
      s * (π / 2 - d.tht2) - π / 2

   elseif r == 10
      s * (d.tht2 - d.tht1) - d.tht2

   else
      throw(ArgumentError("Yy for peanut expects region 1–10; got r=$r"))
   end

   st, ct = sincos(th)

   f = cos(2 * th)
   hP = sqrt(f + sqrt(d.P + f^2))

   return d.B + d.R * hP * st

end

function fill_FTable!(tbl::FTable, d::peanut, r::Int)

   vmin, vmax, N = tbl.vmin, tbl.vmax, tbl.N
   h = (vmax - vmin) / (N - 1)

   P1, P2, P3 = tbl.P1, tbl.P2, tbl.P3

   Rα = d.R * d.alpha
   L1h = d.L1 / 2
   L2h = d.L2 / 2

   if r == 1

      ak = 2 * d.tht1
      bk = -d.tht1
      Xxv = d.A + Rα + L2h

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * h

         Xyv = d.B + d.L1 * (2v̂ - 1) / 2

         th = muladd(ak, v̂, bk)
         st, ct = sincos(th)
         f = cos(2 * th)
         hP = sqrt(f + sqrt(d.P + f^2))

         Yxv = d.A + d.R * hP * ct
         Yyv = d.B + d.R * hP * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 2

      ak = d.tht2 - d.tht1
      bk = d.tht1
      Xyv = d.B + L1h

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * h

         Xxv = d.A + Rα + d.L2 * (1 - 2v̂) / 2

         th = muladd(ak, v̂, bk)
         st, ct = sincos(th)
         f = cos(2 * th)
         hP = sqrt(f + sqrt(d.P + f^2))

         Yxv = d.A + d.R * hP * ct
         Yyv = d.B + d.R * hP * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 3

      ak = π / 2 - d.tht2
      bk = d.tht2
      Cx = Rα - L2h
      Xyv = d.B

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * h

         Xxv = d.A + (1 - v̂) * Cx

         th = muladd(ak, v̂, bk)
         st, ct = sincos(th)
         f = cos(2 * th)
         hP = sqrt(f + sqrt(d.P + f^2))

         Yxv = d.A + d.R * hP * ct
         Yyv = d.B + d.R * hP * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 4

      ak = π / 2 - d.tht2
      bk = π / 2
      Cx = L2h - Rα
      Xyv = d.B

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * h

         Xxv = d.A + v̂ * Cx

         th = muladd(ak, v̂, bk)
         st, ct = sincos(th)
         f = cos(2 * th)
         hP = sqrt(f + sqrt(d.P + f^2))

         Yxv = d.A + d.R * hP * ct
         Yyv = d.B + d.R * hP * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 5

      ak = d.tht2 - d.tht1
      bk = π - d.tht2
      Xyv = d.B + L1h

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * h

         Xxv = d.A - Rα + d.L2 * (1 - 2v̂) / 2

         th = muladd(ak, v̂, bk)
         st, ct = sincos(th)
         f = cos(2 * th)
         hP = sqrt(f + sqrt(d.P + f^2))

         Yxv = d.A + d.R * hP * ct
         Yyv = d.B + d.R * hP * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 6

      ak = 2 * d.tht1
      bk = π - d.tht1
      Xxv = d.A - Rα - L2h

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * h

         Xyv = d.B - d.L1 * (2v̂ - 1) / 2

         th = muladd(ak, v̂, bk)
         st, ct = sincos(th)
         f = cos(2 * th)
         hP = sqrt(f + sqrt(d.P + f^2))

         Yxv = d.A + d.R * hP * ct
         Yyv = d.B + d.R * hP * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 7

      ak = d.tht2 - d.tht1
      bk = π + d.tht1
      Xyv = d.B - L1h

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * h

         Xxv = d.A - Rα - d.L2 * (1 - 2v̂) / 2

         th = muladd(ak, v̂, bk)
         st, ct = sincos(th)
         f = cos(2 * th)
         hP = sqrt(f + sqrt(d.P + f^2))

         Yxv = d.A + d.R * hP * ct
         Yyv = d.B + d.R * hP * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 8

      ak = π / 2 - d.tht2
      bk = π + d.tht2
      Cx = L2h - Rα
      Xyv = d.B

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * h

         Xxv = d.A + (1 - v̂) * Cx

         th = muladd(ak, v̂, bk)
         st, ct = sincos(th)
         f = cos(2 * th)
         hP = sqrt(f + sqrt(d.P + f^2))

         Yxv = d.A + d.R * hP * ct
         Yyv = d.B + d.R * hP * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 9

      ak = π / 2 - d.tht2
      bk = -π / 2
      Cx = Rα - L2h
      Xyv = d.B

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * h

         Xxv = d.A + v̂ * Cx

         th = muladd(ak, v̂, bk)
         st, ct = sincos(th)
         f = cos(2 * th)
         hP = sqrt(f + sqrt(d.P + f^2))

         Yxv = d.A + d.R * hP * ct
         Yyv = d.B + d.R * hP * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 10

      ak = d.tht2 - d.tht1
      bk = -d.tht2
      Xyv = d.B - L1h

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * h

         Xxv = d.A + Rα - d.L2 * (1 - 2v̂) / 2

         th = muladd(ak, v̂, bk)
         st, ct = sincos(th)
         f = cos(2 * th)
         hP = sqrt(f + sqrt(d.P + f^2))

         Yxv = d.A + d.R * hP * ct
         Yyv = d.B + d.R * hP * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   else

      throw(ArgumentError("fill_FTable! for peanut expects region 1–10; got r=$r"))

   end

   tbl.reg = r
   return tbl

end

@inline function f1I(t::Float64, s::Float64,
   d::peanut, u::Float64, v::Float64, r::Int)

   t̂ = (t + 1) / 2
   ŝ = (s + 1) / 2

   Xxv = Xx(ŝ, d, r)
   Yxv = Yx(ŝ, d, r)

   return (1 - t̂) * Xxv + t̂ * Yxv - u

end

@inline function f2I(t::Float64, s::Float64,
   d::peanut, u::Float64, v::Float64, r::Int)

   t̂ = (t + 1) / 2
   ŝ = (s + 1) / 2

   Xyv = Xy(ŝ, d, r)
   Yyv = Yy(ŝ, d, r)

   return (1 - t̂) * Xyv + t̂ * Yyv - v

end

@inline function JinvI(t::Float64, s::Float64,
   d::peanut, u::Float64, v::Float64, r::Int)

   return jinvmap(d, t, s, r)

end

@inline function f_cont(v̂::Float64, d::peanut,
   u::Float64, v::Float64, r::Int)

   Rα = d.R * d.alpha
   L1h = d.L1 / 2
   L2h = d.L2 / 2

   if r == 1

      Xxv = d.A + Rα + L2h
      Xyv = d.B + d.L1 * (2v̂ - 1) / 2
      th = (2v̂ - 1) * d.tht1

   elseif r == 2

      Xxv = d.A + Rα + d.L2 * (1 - 2v̂) / 2
      Xyv = d.B + L1h
      th = v̂ * (d.tht2 - d.tht1) + d.tht1

   elseif r == 3

      Xxv = d.A + (1 - v̂) * (Rα - L2h)
      Xyv = d.B
      th = v̂ * (π / 2 - d.tht2) + d.tht2

   elseif r == 4

      Xxv = d.A + v̂ * (L2h - Rα)
      Xyv = d.B
      th = v̂ * (π / 2 - d.tht2) + π / 2

   elseif r == 5

      Xxv = d.A - Rα + d.L2 * (1 - 2v̂) / 2
      Xyv = d.B + L1h
      th = v̂ * (d.tht2 - d.tht1) + π - d.tht2

   elseif r == 6

      Xxv = d.A - Rα - L2h
      Xyv = d.B - d.L1 * (2v̂ - 1) / 2
      th = 2v̂ * d.tht1 + π - d.tht1

   elseif r == 7

      Xxv = d.A - Rα - d.L2 * (1 - 2v̂) / 2
      Xyv = d.B - L1h
      th = v̂ * (d.tht2 - d.tht1) + π + d.tht1

   elseif r == 8

      Xxv = d.A + (1 - v̂) * (L2h - Rα)
      Xyv = d.B
      th = v̂ * (π / 2 - d.tht2) + π + d.tht2

   elseif r == 9

      Xxv = d.A + v̂ * (Rα - L2h)
      Xyv = d.B
      th = v̂ * (π / 2 - d.tht2) - π / 2

   elseif r == 10

      Xxv = d.A + Rα - d.L2 * (1 - 2v̂) / 2
      Xyv = d.B - L1h
      th = v̂ * (d.tht2 - d.tht1) - d.tht2

   else

      throw(ArgumentError("f_cont for peanut expects region 1–10; got r=$r"))

   end

   st, ct = sincos(th)
   f = cos(2 * th)
   hP = sqrt(f + sqrt(d.P + f^2))

   Yxv = d.A + d.R * hP * ct
   Yyv = d.B + d.R * hP * st

   dx = Yxv - Xxv
   dy = Yyv - Xyv

   term1 = muladd(u, dy, -v * dx)
   term2 = muladd(Xyv, Yxv, -Xxv * Yyv)

   return term1 + term2

end

function mapinv(tbl::FTable, d::peanut, u::Float64,
   v::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]
   r = p.reg

   # ----- Stage 1: 2D Newton -----
   tN, sN = newtonR2D(f1I, f2I, JinvI,
      0.0, 0.0, 4, d, u, v, r; tol = 1e-15)

   if tN !== :max

      x = xi_inv((1 + tN) / 2, p.ck0, p.ck1)
      y = xi_inv((1 + sN) / 2, p.tk0, p.tk1)

      return x, y

   end

   # ----- Stage 2: BIS / FTable method -----
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

# Separate inverse for the two rectangular regions.
function mapinv(d::peanut, u::Float64,
   v::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]
   r = p.reg

   if r == 11

      xξ = (u - d.A - d.R * d.alpha + d.L2 / 2) / d.L2
      yξ = (v - d.B + d.L1 / 2) / d.L1

      if d.L2 >= d.L1
         Z1 = xi_inv(xξ, p.tk0, p.tk1)
         Z2 = xi_inv(yξ, p.ck0, p.ck1)
      else
         Z1 = xi_inv(xξ, p.ck0, p.ck1)
         Z2 = xi_inv(yξ, p.tk0, p.tk1)
      end

      return Z1, Z2

   elseif r == 12

      xξ = (u - d.A + d.R * d.alpha + d.L2 / 2) / d.L2
      yξ = (v - d.B + d.L1 / 2) / d.L1

      if d.L2 >= d.L1
         Z1 = xi_inv(xξ, p.tk0, p.tk1)
         Z2 = xi_inv(yξ, p.ck0, p.ck1)
      else
         Z1 = xi_inv(xξ, p.ck0, p.ck1)
         Z2 = xi_inv(yξ, p.tk0, p.tk1)
      end

      return Z1, Z2

   end

   throw(ArgumentError("mapinv(d::peanut) is only for rectangular regions 11 and 12; got reg=$r"))

end

"""
    ptconv(d::peanut, t1::Float64, t2::Float64, idx::Int, ptdest::String)

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
function ptconv(d::peanut, t1::Float64, t2::Float64, idx::Int, ptdest::String)

   if ptdest == "to_pth"

      r = idx

      for k in 1:d.Npat

         p = d.pths[k]
         p.reg == r || continue

         if r != 11 && r != 12

            t1k = xi_inv((t1 + 1) / 2, p.ck0, p.ck1)
            t2k = xi_inv((t2 + 1) / 2, p.tk0, p.tk1)

         else

            if d.L1 > d.L2

               t1k = xi_inv((t1 + 1) / 2, p.ck0, p.ck1)
               t2k = xi_inv((t2 + 1) / 2, p.tk0, p.tk1)

            else

               t1k = xi_inv((t1 + 1) / 2, p.tk0, p.tk1)
               t2k = xi_inv((t2 + 1) / 2, p.ck0, p.ck1)

            end

         end

         if abs(t1k) ≤ 1 && abs(t2k) ≤ 1
            return t1k, t2k, k
         end

      end

      error("ptconv(to_pth): no patch in region $r contained the point.")

   elseif ptdest == "to_reg"

      k = idx
      p = d.pths[k]
      r = p.reg

      if r != 11 && r != 12

         t1r = (p.ck1 - p.ck0) * t1 + (p.ck1 + p.ck0) - 1
         t2r = (p.tk1 - p.tk0) * t2 + (p.tk1 + p.tk0) - 1

      else

         if d.L1 > d.L2

            t1r = (p.ck1 - p.ck0) * t1 + (p.ck1 + p.ck0) - 1
            t2r = (p.tk1 - p.tk0) * t2 + (p.tk1 + p.tk0) - 1

         else

            t1r = (p.tk1 - p.tk0) * t1 + (p.tk1 + p.tk0) - 1
            t2r = (p.ck1 - p.ck0) * t2 + (p.ck1 + p.ck0) - 1

         end

      end

      return t1r, t2r, r

   else

      error("ptconv: ptdest must be \"to_pth\" or \"to_reg\"")

   end

end

"""
    mapinv2(d::peanut, t1::Float64, t2::Float64, k2::Int, k::Int) -> Tuple{Float64,Float64}

Given a point `(u,v) = τₖ₂(t1,t2)` on patch `k2`, return its reference
coordinates on patch `k` in `[-1,1]^2`.

This differs from `mapinv` because we already know `(u,v)` comes from
`(t1,t2)` on `k2`.

This should only be called when patches `k2` and `k` are in the same region.
"""
function mapinv2(d::peanut, t1::Float64, t2::Float64,
   k2::Int, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]

   # Convert the point from patch coordinates on k2 to region coordinates.
   tr1, tr2, _ = ptconv(d, t1, t2, k2, "to_reg")

   # Map from [-1,1] to [0,1].
   û = (tr1 + 1) / 2
   v̂ = (tr2 + 1) / 2

   # Rectangular regions 11 and 12 use the same axis-swap convention
   # as mapxy/mapinv when L2 >= L1.
   if d.L2 >= d.L1 && (p.reg == 11 || p.reg == 12)

      t1k = xi_inv(û, p.tk0, p.tk1)
      t2k = xi_inv(v̂, p.ck0, p.ck1)

   else

      t1k = xi_inv(û, p.ck0, p.ck1)
      t2k = xi_inv(v̂, p.tk0, p.tk1)

   end

   return t1k, t2k

end

"""
Factor that goes linearly to zero as we approach the boundary.

- If patch `k` is in rectangular region 11 or 12, returns `1`.
- Otherwise uses

    d = 1 - ck1 + (ck1 - ck0) * t / 2

  then raises elementwise to the power `s-1` if `s ≥ 0.5`, else to `s`.

Notes:
- `t` is expected to be `1 - t_actual`.
- `t` may be a scalar `Float64` or any `StridedArray{Float64}`.
"""
function dfunc(d::peanut, k::Int, t::Float64, s::Float64)::Float64

   p = d.pths[k]

   if p.reg == 11 || p.reg == 12
      return 1.0
   end

   hc = p.ck1 - p.ck0
   val = (1.0 - p.ck1) + hc * t / 2
   exp = s ≥ 0.5 ? (s - 1) : s

   return val^exp

end

function dfunc!(out::StridedArray{Float64}, d::peanut, k::Int,
   t::StridedArray{Float64}, s::Float64)

   p = d.pths[k]

   if p.reg == 11 || p.reg == 12

      fill!(out, 1.0)

   else

      exp = s ≥ 0.5 ? (s - 1) : s
      αc = (p.ck1 - p.ck0) / 2
      βc = 1.0 - p.ck1

      @inbounds for i in eachindex(t)
         out[i] = muladd(αc, t[i], βc)^exp
      end

   end

   return nothing

end

"""
    gamderhigher(d::peanut, th::Float64)

Return centered peanut boundary coordinates and derivatives up to 4th order.

Returns
    gx, gy, dgx, dgy, d2gx, d2gy,
    d3gx, d3gy, d4gx, d4gy

where gx(th) = R*h(th)*cos(th), gy(th) = R*h(th)*sin(th)
and h(th) = sqrt(cos(2th) + sqrt(P + cos(2th)^2)).
"""
function gamderhigher(d::peanut, th::Float64)

   st, ct = sincos(th)

   # f = cos(2th)
   f = cos(2 * th)
   f1 = -2.0 * sin(2 * th)
   f2 = -4.0 * f
   f3 = -4.0 * f1
   f4 = 16.0 * f

   # l = P + f^2
   l = d.P + f^2
   l1 = 2.0 * f * f1
   l2 = 2.0 * f1^2 + 2.0 * f * f2
   l3 = 6.0 * f1 * f2 + 2.0 * f * f3
   l4 = 6.0 * f2^2 + 8.0 * f1 * f3 + 2.0 * f * f4

   # q = sqrt(l)
   q = sqrt(l)

   q1 = l1 / (2.0 * q)

   q2 = (l2 - 2.0 * q1^2) / (2.0 * q)

   q3 = (l3 - 6.0 * q1 * q2) / (2.0 * q)

   q4 = (l4 - 6.0 * q2^2 - 8.0 * q1 * q3) / (2.0 * q)

   # h = sqrt(f + q)
   m = f + q
   m1 = f1 + q1
   m2 = f2 + q2
   m3 = f3 + q3
   m4 = f4 + q4

   h = sqrt(m)

   h1 = m1 / (2.0 * h)

   h2 = (m2 - 2.0 * h1^2) / (2.0 * h)

   h3 = (m3 - 6.0 * h1 * h2) / (2.0 * h)

   h4 = (m4 - 6.0 * h2^2 - 8.0 * h1 * h3) / (2.0 * h)

   gx = d.R * h * ct
   gy = d.R * h * st

   dgx = d.R * (h1 * ct - h * st)
   dgy = d.R * (h1 * st + h * ct)

   d2gx = d.R * (h2 * ct - 2.0 * h1 * st - h * ct)
   d2gy = d.R * (h2 * st + 2.0 * h1 * ct - h * st)

   d3gx = d.R * (h3 * ct - 3.0 * h2 * st - 3.0 * h1 * ct + h * st)
   d3gy = d.R * (h3 * st + 3.0 * h2 * ct - 3.0 * h1 * st - h * ct)

   d4gx = d.R * (h4 * ct - 4.0 * h3 * st - 6.0 * h2 * ct + 4.0 * h1 * st + h * ct)
   d4gy = d.R * (h4 * st + 4.0 * h3 * ct - 6.0 * h2 * st - 4.0 * h1 * ct + h * st)

   return gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy, d4gx, d4gy

end

"""
    diff_map!(out, Zx, Zy, DJ, d::peanut, u, v, u2, v2, du, dv, k; tol=1e-4)

Fill `out` with `||τ(u,v) - τ(u₂,v₂)||` on patch `k`.

No allocations. Uses `mapxy_Dmap!` for the far values and Taylor expansion
near `(u,v)` to avoid cancellation.
"""
function diff_map!(out::Matrix{Float64},
   Zx::Matrix{Float64}, Zy::Matrix{Float64}, DJ::StridedArray{Float64},
   d::peanut, u::Float64, v::Float64,
   u2::Matrix{Float64}, v2::Matrix{Float64},
   du::AbstractVector, dv::AbstractVector, k::Int;
   tol::Float64=1e-4)

   nd_u = size(out, 1)
   nd_v = size(out, 2)

   p = d.pths[k]

   αc = (p.ck1 - p.ck0) / 2
   αt = (p.tk1 - p.tk0) / 2

   # Important: fill Zx, Zy, and DJ before any early return.
   mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)

   if p.reg == 11 || p.reg == 12

      if d.L2 >= d.L1
         cDx = d.L2 * αt
         cDy = d.L1 * αc
      else
         cDx = d.L2 * αc
         cDy = d.L1 * αt
      end

      @inbounds for j in 1:nd_v
         dvj = dv[j]
         @inbounds for i in 1:nd_u
            out[i, j] = hypot(cDx * du[i], cDy * dvj)
         end
      end

      return nothing

   end

   uhat = muladd(αc, u, p.ck0 + αc)
   vhat = muladd(αt, v, p.tk0 + αt)

   Rα = d.R * d.alpha
   L1h = d.L1 / 2
   L2h = d.L2 / 2

   if p.reg == 1

      ak = 2 * d.tht1
      bk = -d.tht1

      Xx = Rα + L2h
      Xy = d.L1 * (2vhat - 1) / 2

      dXx = 0.0
      dXy = d.L1

   elseif p.reg == 2

      ak = d.tht2 - d.tht1
      bk = d.tht1

      Xx = Rα + d.L2 * (1 - 2vhat) / 2
      Xy = L1h

      dXx = -d.L2
      dXy = 0.0

   elseif p.reg == 3

      ak = π / 2 - d.tht2
      bk = d.tht2

      Xx = (1 - vhat) * (Rα - L2h)
      Xy = 0.0

      dXx = -(Rα - L2h)
      dXy = 0.0

   elseif p.reg == 4

      ak = π / 2 - d.tht2
      bk = π / 2

      Xx = vhat * (L2h - Rα)
      Xy = 0.0

      dXx = L2h - Rα
      dXy = 0.0

   elseif p.reg == 5

      ak = d.tht2 - d.tht1
      bk = π - d.tht2

      Xx = -Rα + d.L2 * (1 - 2vhat) / 2
      Xy = L1h

      dXx = -d.L2
      dXy = 0.0

   elseif p.reg == 6

      ak = 2 * d.tht1
      bk = π - d.tht1

      Xx = -Rα - L2h
      Xy = -d.L1 * (2vhat - 1) / 2

      dXx = 0.0
      dXy = -d.L1

   elseif p.reg == 7

      ak = d.tht2 - d.tht1
      bk = π + d.tht1

      Xx = -Rα - d.L2 * (1 - 2vhat) / 2
      Xy = -L1h

      dXx = d.L2
      dXy = 0.0

   elseif p.reg == 8

      ak = π / 2 - d.tht2
      bk = π + d.tht2

      Xx = (1 - vhat) * (L2h - Rα)
      Xy = 0.0

      dXx = -(L2h - Rα)
      dXy = 0.0

   elseif p.reg == 9

      ak = π / 2 - d.tht2
      bk = -π / 2

      Xx = vhat * (Rα - L2h)
      Xy = 0.0

      dXx = Rα - L2h
      dXy = 0.0

   elseif p.reg == 10

      ak = d.tht2 - d.tht1
      bk = -d.tht2

      Xx = Rα - d.L2 * (1 - 2vhat) / 2
      Xy = -L1h

      dXx = d.L2
      dXy = 0.0

   else

      throw(ArgumentError("diff_map! for peanut expects region 1–12; got reg=$(p.reg)"))

   end

   th = muladd(ak, vhat, bk)

   gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy, d4gx, d4gy =
      gamderhigher(d, th)

   q = αt * ak
   q2 = q * q
   q3 = q2 * q
   q4 = q2 * q2

   dux = αc * (gx - Xx)
   duy = αc * (gy - Xy)

   dvx = (1 - uhat) * αt * dXx + uhat * q * dgx
   dvy = (1 - uhat) * αt * dXy + uhat * q * dgy

   duvx = αc * (q * dgx - αt * dXx)
   duvy = αc * (q * dgy - αt * dXy)

   dv2x = uhat * q2 * d2gx
   dv2y = uhat * q2 * d2gy

   duv2x = αc * q2 * d2gx
   duv2y = αc * q2 * d2gy

   dv3x = uhat * q3 * d3gx
   dv3y = uhat * q3 * d3gy

   duv3x = αc * q3 * d3gx
   duv3y = αc * q3 * d3gy

   dv4x = uhat * q4 * d4gx
   dv4y = uhat * q4 * d4gy

   tux, tvy = mapxy(d, u, v, k)

   @inbounds for j in 1:nd_v

      dvj = dv[j]
      dvj2 = dvj * dvj
      dvj3 = dvj2 * dvj
      dvj4 = dvj2 * dvj2

      @inbounds for i in 1:nd_u

         uu = u2[i, j]
         vv = v2[i, j]

         if abs(u - uu) < tol && abs(v - vv) < tol

            dui = du[i]

            Dx = (dui * dux + dvj * dvx) -
                 (dui * dvj * duvx + dvj2 * dv2x / 2) +
                 (dui * dvj2 * duv2x / 2 + dvj3 * dv3x / 6) -
                 (dui * dvj3 * duv3x / 6 + dvj4 * dv4x / 24)

            Dy = (dui * duy + dvj * dvy) -
                 (dui * dvj * duvy + dvj2 * dv2y / 2) +
                 (dui * dvj2 * duv2y / 2 + dvj3 * dv3y / 6) -
                 (dui * dvj3 * duv3y / 6 + dvj4 * dv4y / 24)

            out[i, j] = hypot(Dx, Dy)

         else

            out[i, j] = hypot(tux - Zx[i, j], tvy - Zy[i, j])

         end

      end

   end

   return nothing

end

"""
    gamderhigher6(d::peanut, th::Float64)

Return centered peanut boundary coordinates and derivatives up to 6th order.

Returns gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy, d4gx, d4gy, d5gx, d5gy, d6gx, d6gy.
"""
function gamderhigher6(d::peanut, th::Float64)

   st, ct = sincos(th)

   h, h1, h2, h3, h4, h5, h6 = hderhigher(d, th)

   gx = d.R * h * ct
   gy = d.R * h * st

   dgx = d.R * (h1 * ct - h * st)
   dgy = d.R * (h1 * st + h * ct)

   d2gx = d.R * (h2 * ct - 2.0 * h1 * st - h * ct)
   d2gy = d.R * (h2 * st + 2.0 * h1 * ct - h * st)

   d3gx = d.R * (h3 * ct - 3.0 * h2 * st - 3.0 * h1 * ct + h * st)
   d3gy = d.R * (h3 * st + 3.0 * h2 * ct - 3.0 * h1 * st - h * ct)

   d4gx = d.R * (h4 * ct - 4.0 * h3 * st - 6.0 * h2 * ct +
                 4.0 * h1 * st + h * ct)

   d4gy = d.R * (h4 * st + 4.0 * h3 * ct - 6.0 * h2 * st -
                 4.0 * h1 * ct + h * st)

   d5gx = d.R * (h5 * ct - 5.0 * h4 * st - 10.0 * h3 * ct +
                 10.0 * h2 * st + 5.0 * h1 * ct - h * st)

   d5gy = d.R * (h5 * st + 5.0 * h4 * ct - 10.0 * h3 * st -
                 10.0 * h2 * ct + 5.0 * h1 * st + h * ct)

   d6gx = d.R * (h6 * ct - 6.0 * h5 * st - 15.0 * h4 * ct +
                 20.0 * h3 * st + 15.0 * h2 * ct -
                 6.0 * h1 * st - h * ct)

   d6gy = d.R * (h6 * st + 6.0 * h5 * ct - 15.0 * h4 * st -
                 20.0 * h3 * ct + 15.0 * h2 * st +
                 6.0 * h1 * ct - h * st)

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
             DJ::StridedArray{Float64}, d::peanut,
             u::Float64, v::Float64,
             u2::Matrix{Float64}, v2::Matrix{Float64}, r::Matrix{Float64},
             du::AbstractVector, dv::AbstractVector, k::Int;
             tol = 1e-3)

Compute

    ||τ(u,v) - τ(u₂,v₂)|| / r

for the `k`-th patch, where

    u₂ = u - r .* du
    v₂ = v - r .* dv

No allocations. Uses `mapxy_Dmap!` and Taylor correction near `(u,v)`.
"""
function diff_rmap!(out::Matrix{Float64},
   Zx::Matrix{Float64}, Zy::Matrix{Float64}, DJ::StridedArray{Float64},
   d::peanut, u::Float64, v::Float64,
   u2::Matrix{Float64}, v2::Matrix{Float64}, r::Matrix{Float64},
   du::AbstractVector, dv::AbstractVector, k::Int;
   tol::Float64 = 1e-3)

   nt = size(out, 1)
   nr = size(out, 2)

   p = d.pths[k]

   αc = (p.ck1 - p.ck0) / 2
   αt = (p.tk1 - p.tk0) / 2

   if p.reg == 11 || p.reg == 12

      if d.L2 >= d.L1
         cDx = d.L2 * αt
         cDy = d.L1 * αc
      else
         cDx = d.L2 * αc
         cDy = d.L1 * αt
      end

      fill!(DJ, cDx * cDy)

      @inbounds for i in 1:nt
         dui = du[i]
         dvi = dv[i]
         hD = hypot(cDx * dui, cDy * dvi)

         @inbounds for j in 1:nr
            out[i, j] = hD
         end
      end

      return nothing

   end

   mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)

   uhat = muladd(αc, u, p.ck0 + αc)
   vhat = muladd(αt, v, p.tk0 + αt)

   Rα = d.R * d.alpha
   L1h = d.L1 / 2
   L2h = d.L2 / 2

   if p.reg == 1

      ak = 2 * d.tht1
      bk = -d.tht1

      Xx = Rα + L2h
      Xy = d.L1 * (2vhat - 1) / 2

      dXx = 0.0
      dXy = d.L1

   elseif p.reg == 2

      ak = d.tht2 - d.tht1
      bk = d.tht1

      Xx = Rα + d.L2 * (1 - 2vhat) / 2
      Xy = L1h

      dXx = -d.L2
      dXy = 0.0

   elseif p.reg == 3

      ak = π / 2 - d.tht2
      bk = d.tht2

      Xx = (1 - vhat) * (Rα - L2h)
      Xy = 0.0

      dXx = -(Rα - L2h)
      dXy = 0.0

   elseif p.reg == 4

      ak = π / 2 - d.tht2
      bk = π / 2

      Xx = vhat * (L2h - Rα)
      Xy = 0.0

      dXx = L2h - Rα
      dXy = 0.0

   elseif p.reg == 5

      ak = d.tht2 - d.tht1
      bk = π - d.tht2

      Xx = -Rα + d.L2 * (1 - 2vhat) / 2
      Xy = L1h

      dXx = -d.L2
      dXy = 0.0

   elseif p.reg == 6

      ak = 2 * d.tht1
      bk = π - d.tht1

      Xx = -Rα - L2h
      Xy = -d.L1 * (2vhat - 1) / 2

      dXx = 0.0
      dXy = -d.L1

   elseif p.reg == 7

      ak = d.tht2 - d.tht1
      bk = π + d.tht1

      Xx = -Rα - d.L2 * (1 - 2vhat) / 2
      Xy = -L1h

      dXx = d.L2
      dXy = 0.0

   elseif p.reg == 8

      ak = π / 2 - d.tht2
      bk = π + d.tht2

      Xx = (1 - vhat) * (L2h - Rα)
      Xy = 0.0

      dXx = -(L2h - Rα)
      dXy = 0.0

   elseif p.reg == 9

      ak = π / 2 - d.tht2
      bk = -π / 2

      Xx = vhat * (Rα - L2h)
      Xy = 0.0

      dXx = Rα - L2h
      dXy = 0.0

   elseif p.reg == 10

      ak = d.tht2 - d.tht1
      bk = -d.tht2

      Xx = Rα - d.L2 * (1 - 2vhat) / 2
      Xy = -L1h

      dXx = d.L2
      dXy = 0.0

   else

      throw(ArgumentError("diff_rmap! for peanut expects region 1–12; got reg=$(p.reg)"))

   end

   th = muladd(ak, vhat, bk)

   gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy, 
   d4gx, d4gy, d5gx, d5gy, d6gx, d6gy = gamderhigher6(d, th)

   q = αt * ak

   q2 = q * q
   q3 = q2 * q
   q4 = q2 * q2
   q5 = q4 * q
   q6 = q3 * q3

   dux = αc * (gx - Xx)
   duy = αc * (gy - Xy)

   dvx = (1.0 - uhat) * αt * dXx + uhat * q * dgx
   dvy = (1.0 - uhat) * αt * dXy + uhat * q * dgy

   duvx = αc * (q * dgx - αt * dXx)
   duvy = αc * (q * dgy - αt * dXy)

   dv2x = uhat * q2 * d2gx
   dv2y = uhat * q2 * d2gy

   duv2x = αc * q2 * d2gx
   duv2y = αc * q2 * d2gy

   dv3x = uhat * q3 * d3gx
   dv3y = uhat * q3 * d3gy

   duv3x = αc * q3 * d3gx
   duv3y = αc * q3 * d3gy

   dv4x = uhat * q4 * d4gx
   dv4y = uhat * q4 * d4gy

   duv4x = αc * q4 * d4gx
   duv4y = αc * q4 * d4gy

   dv5x = uhat * q5 * d5gx
   dv5y = uhat * q5 * d5gy

   duv5x = αc * q5 * d5gx
   duv5y = αc * q5 * d5gy

   dv6x = uhat * q6 * d6gx
   dv6y = uhat * q6 * d6gy

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
    Dwall(d::peanut, u::Float64, v::Float64, k::Int) -> dvx, dvy

Derivative of the wall curve w.r.t. `v` at the singular side `u = ±1`
for the `k`-th non-rectangular patch.

Returns a Tuple `(dvx, dvy)`.
"""
function Dwall(d::peanut, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]
   reg = p.reg

   if reg == 11 || reg == 12
      throw(ArgumentError("Dwall is for non-rectangular patches; got reg=$reg"))
   end

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   αu = 0.5 * hc
   αv = 0.5 * ht

   βu = p.ck0 + αu
   βv = p.tk0 + αv

   xi1 = muladd(αu, u, βu)
   xi2 = muladd(αv, v, βv)

   Rα = d.R * d.alpha
   L2h = d.L2 / 2

   if reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1
      th = muladd(Δth, xi2, th0)

      _, _, dgx, dgy = gamder(d, th)

      dvx = xi1 * Δth * dgx
      dvy = (1.0 - xi1) * d.L1 + xi1 * Δth * dgy

   elseif reg == 2

      Δth = d.tht2 - d.tht1
      th0 = d.tht1
      th = muladd(Δth, xi2, th0)

      _, _, dgx, dgy = gamder(d, th)

      dvx = -(1.0 - xi1) * d.L2 + xi1 * Δth * dgx
      dvy = xi1 * Δth * dgy

   elseif reg == 3

      Δth = π / 2 - d.tht2
      th0 = d.tht2
      th = muladd(Δth, xi2, th0)

      _, _, dgx, dgy = gamder(d, th)

      dvx = -(1.0 - xi1) * (Rα - L2h) + xi1 * Δth * dgx
      dvy = xi1 * Δth * dgy

   elseif reg == 4

      Δth = π / 2 - d.tht2
      th0 = π / 2
      th = muladd(Δth, xi2, th0)

      _, _, dgx, dgy = gamder(d, th)

      dvx = (1.0 - xi1) * (L2h - Rα) + xi1 * Δth * dgx
      dvy = xi1 * Δth * dgy

   elseif reg == 5

      Δth = d.tht2 - d.tht1
      th0 = π - d.tht2
      th = muladd(Δth, xi2, th0)

      _, _, dgx, dgy = gamder(d, th)

      dvx = -(1.0 - xi1) * d.L2 + xi1 * Δth * dgx
      dvy = xi1 * Δth * dgy

   elseif reg == 6

      Δth = 2 * d.tht1
      th0 = π - d.tht1
      th = muladd(Δth, xi2, th0)

      _, _, dgx, dgy = gamder(d, th)

      dvx = xi1 * Δth * dgx
      dvy = -(1.0 - xi1) * d.L1 + xi1 * Δth * dgy

   elseif reg == 7

      Δth = d.tht2 - d.tht1
      th0 = π + d.tht1
      th = muladd(Δth, xi2, th0)

      _, _, dgx, dgy = gamder(d, th)

      dvx = (1.0 - xi1) * d.L2 + xi1 * Δth * dgx
      dvy = xi1 * Δth * dgy

   elseif reg == 8

      Δth = π / 2 - d.tht2
      th0 = π + d.tht2
      th = muladd(Δth, xi2, th0)

      _, _, dgx, dgy = gamder(d, th)

      dvx = -(1.0 - xi1) * (L2h - Rα) + xi1 * Δth * dgx
      dvy = xi1 * Δth * dgy

   elseif reg == 9

      Δth = π / 2 - d.tht2
      th0 = -π / 2
      th = muladd(Δth, xi2, th0)

      _, _, dgx, dgy = gamder(d, th)

      dvx = (1.0 - xi1) * (Rα - L2h) + xi1 * Δth * dgx
      dvy = xi1 * Δth * dgy

   elseif reg == 10

      Δth = d.tht2 - d.tht1
      th0 = -d.tht2
      th = muladd(Δth, xi2, th0)

      _, _, dgx, dgy = gamder(d, th)

      dvx = (1.0 - xi1) * d.L2 + xi1 * Δth * dgx
      dvy = xi1 * Δth * dgy

   else

      throw(ArgumentError("Dwall is defined only for regions 1–10; got reg=$reg"))

   end

   return αv * dvx, αv * dvy

end

"""
   refine!(d::peanut, Nc::Int, Nt::Int, K::AbstractVector{<:Int})

Subdivide each patch in `K` into `Nc * Nt` subpatches, split in `c` and `t`.
Updates `d.pths`, `d.Npat`, recomputes `d.kd`, and rebuilds `d.Qpts` / `d.Qptsbd`.

Notes
- Patches are re-sorted lexicographically by `(reg, ck0, ck1, tk0, tk1)`.
- `Qpts` uses `boundquad!(V, d, k)`.
- `Qptsbd` uses `boundquadbd!(V, d, k)` for patches in `d.kd`.
"""
function refine!(d::peanut, Nc::Int, Nt::Int, K::Vector{Int})

   @assert Nc ≥ 1 && Nt ≥ 1 "Nc and Nt must be ≥ 1"
   @assert !isempty(K) "K (set of patches to refine) is empty"

   # Create the subdivided children.
   children = Patch[]
   sizehint!(children, length(K) * Nc * Nt)

   for k in K

      p = d.pths[k]

      cc = range(p.ck0, p.ck1; length = Nc + 1)
      tt = range(p.tk0, p.tk1; length = Nt + 1)

      for i in 1:Nc, j in 1:Nt
         push!(children, Patch(p.reg, cc[i], cc[i + 1], tt[j], tt[j + 1]))
      end

   end

   # Merge and sort.
   d.pths = vcat([d.pths[i] for i in 1:d.Npat if !(i in K)], children)

   sort!(d.pths, by = q -> (q.reg, q.ck0, q.ck1, q.tk0, q.tk1))

   # Update counts.
   d.Npat = length(d.pths)

   # Boundary-touching patches are regions 1:10 with ck1 == 1.
   d.kd = [k for k in 1:d.Npat if d.pths[k].reg <= 10 && d.pths[k].ck1 == 1.0]

   # Qpts: 8 × Npat.
   d.Qpts = Matrix{Float64}(undef, 8, d.Npat)

   @inbounds for k in 1:d.Npat
      @views V = d.Qpts[:, k]
      boundquad!(V, d, k)
   end

   # Qptsbd: 8 × length(kd).
   d.Qptsbd = Matrix{Float64}(undef, 8, length(d.kd))

   @inbounds for (l, k) in enumerate(d.kd)
      @views V = d.Qptsbd[:, l]
      boundquadbd!(V, d, k)
   end

   return d

end

function Base.show(io::IO, d::peanut)

   println(io, "peanut with properties:")
   println(io, "  (A, B)   = (", d.A, ", ", d.B, ")")
   println(io, "  R        = ", d.R)
   println(io, "  P        = ", d.P)
   println(io, "  tht1     = ", d.tht1)
   println(io, "  tht2     = ", d.tht2)
   println(io, "  alpha    = ", d.alpha)
   println(io, "  (L1, L2) = (", d.L1, ", ", d.L2, ")")
   println(io, "  No. of holes = ", d.nh)

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