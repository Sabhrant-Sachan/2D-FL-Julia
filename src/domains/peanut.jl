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

      # No padding, unlike Matlab.
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

      if abs(Δ) < 2e-3

         Δ2 = Δ * Δ
         Δ3 = Δ2 * Δ
         Δ4 = Δ2 * Δ2

         q = ht1 +
             Δ * ht2 +
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