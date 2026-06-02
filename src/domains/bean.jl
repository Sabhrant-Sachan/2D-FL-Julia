"""
Mutable struct bean

A  — x-coordinate of bean center (Float64)

B  — y-coordinate of bean center (Float64)

R — Scale of the bean (Float64)

L₁ — length of rectangle's y-axis side (Float64)

L₂ — length of rectangle's x-axis side (Float64)

nₕ — number of holes (nonnegative Int) 

θ₁ — A parameter used to constuct the patches in bean. (Float64)

θ₂ — A parameter used to constuct the patches in bean. (Float64)

θ₃ — A parameter used to constuct the patches in bean. (Float64)

kd — indices of patches touching boundary (Vector{Int})  

Nₚₐₜ — total number of patches (≥ 1)  (Int)

pths — Npat stuctures of type Patch. The first
       row `reg` represents region in which the patch is present.
       2,3 row represents partition ck values.
       4,5 row represents partition tk values.
       With this five data values, the patch and the mapping
       is uniquely found and can be used in further 
       calculations. Definitions of ck and tk given in
       constructor of class. It is a matrix of size 5*Npat. 

Qpts — 8*Npat matrix (filled by bd_quad)  

Qptsbd — 8*Npat matrix (filled by bd_quadbd)

"""

struct regionparams_bean

   e₁x::Float64
   e₁y::Float64
   e₂x::Float64
   e₂y::Float64

   q₁x::Float64
   q₁y::Float64
   q₂x::Float64
   q₂y::Float64

   e₃x::Float64
   e₃y::Float64
   e₄x::Float64
   e₄y::Float64
   # common norm of e₁,e₂,q₁,q₂
   nme::Float64   

end

mutable struct bean <: abstractdomain

   A::Float64
   B::Float64
   R::Float64
   L1::Float64
   L2::Float64
   tht1::Float64
   tht2::Float64
   tht3::Float64
   nh::Int

   kd::Vector{Int}
   Npat::Int
   pths::Vector{Patch}
   Qpts::Matrix{Float64}
   Qptsbd::Matrix{Float64}

   RP::regionparams_bean

   function bean(; b, a = nothing,
      A = nothing, B = nothing,
      R = nothing, L1 = nothing, L2 = nothing,
      tht1 = nothing, tht2 = nothing, tht3 = nothing,
      ck = nothing, tk = nothing)

      @assert isa(b, AbstractVector{Int}) && length(b) == 12
      @assert all(b .> 0)

      nh = 0

      a = isnothing(a) ? ceil.(Int, 2 .* b ./ 3) : collect(Int, a)

      @assert length(a) == 12
      @assert all(a .> 0)

      A = Float64(something(A, 0.0))
      B = Float64(something(B, 0.0))
      R = Float64(something(R, 1.0))

      L1 = Float64(something(L1, 0.3 * R))
      L2 = Float64(something(L2, 0.25 * R))

      tht1 = Float64(something(tht1, 0.23))
      tht2 = Float64(something(tht2, 0.225))
      tht3 = Float64(something(tht3, 0.35))

      @assert R > 0.0
      @assert L1 > 0.0
      @assert L2 > 0.0

      Npat = dot(a, b)

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

               if p.reg != 11 && p.reg != 12 && p.ck1 == 1.0
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

      # ------------------------------------------------------------
      # Construction parameters for bean regions.
      # ------------------------------------------------------------
      θ = π / 4 + tht1

      s, c = sincos(θ)

      s2, c2 = s * s, c * c

      s4, c4 = s2 * s2, c2 * c2

      C = cos(2 * tht1)

      e₁x, e₁y = C * s2, C * c2
      e₂x, e₂y = C * c2, C * s2

      q₁x, q₁y = -C * c2, C * s2
      q₂x, q₂y =  C * s2, -C * c2

      e₃x, e₃y = c4, s4
      e₄x, e₄y = s4, c4

      nme = hypot(e₁x, e₁y)

      RP = regionparams_bean(
         e₁x, e₁y, e₂x, e₂y,
         q₁x, q₁y, q₂x, q₂y,
         e₃x, e₃y, e₄x, e₄y,
         nme)

      d = new(A, B, R, L1, L2, tht1, tht2, tht3, nh,
              kd, Npat, pths, Qpts, Qptsbd, RP)

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

function mapx(d::bean, u::Float64, v::Float64, k::Int)::Float64
   x, _ = mapxy(d, u, v, k)
   return x
end

function mapy(d::bean, u::Float64, v::Float64, k::Int)::Float64
   _, y = mapxy(d, u, v, k)
   return y
end

function mapxy(d::bean, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]
   rp = d.RP

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   xi1 = muladd(hc / 2, u, (p.ck0 + p.ck1) / 2)
   xi2 = muladd(ht / 2, v, (p.tk0 + p.tk1) / 2)

   r = p.reg

   if r == 11 || r == 12

      if d.L2 >= d.L1
         xi1 = muladd(ht / 2, u, (p.tk0 + p.tk1) / 2)
         xi2 = muladd(hc / 2, v, (p.ck0 + p.ck1) / 2)
      end

      if r == 11
         Xx = d.A + d.R * rp.e₄x +
              (2xi2 - 1) * d.L1 * rp.e₂x / (2 * rp.nme) +
              d.L2 * xi1 * rp.q₂x / rp.nme

         Xy = d.B + d.R * rp.e₄y +
              (2xi2 - 1) * d.L1 * rp.e₂y / (2 * rp.nme) +
              d.L2 * xi1 * rp.q₂y / rp.nme
      else
         Xx = d.A + d.R * rp.e₃x +
              (2xi2 - 1) * d.L1 * rp.e₁x / (2 * rp.nme) +
              d.L2 * xi1 * rp.q₁x / rp.nme

         Xy = d.B + d.R * rp.e₃y +
              (2xi2 - 1) * d.L1 * rp.e₁y / (2 * rp.nme) +
              d.L2 * xi1 * rp.q₁y / rp.nme
      end

      return Xx, Xy
   end

   if r == 1

      Xx = d.A + d.R * rp.e₄x +
           (2xi2 - 1) * d.L1 * rp.e₂x / (2 * rp.nme) +
           d.L2 * rp.q₂x / rp.nme

      Xy = d.B + d.R * rp.e₄y +
           (2xi2 - 1) * d.L1 * rp.e₂y / (2 * rp.nme) +
           d.L2 * rp.q₂y / rp.nme

      th = xi2 * (d.tht2 + d.tht3) - d.tht3

   elseif r == 2

      Xx = d.A + d.R * rp.e₄x +
           d.L1 * rp.e₂x / (2 * rp.nme) +
           (1 - xi2) * d.L2 * rp.q₂x / rp.nme

      Xy = d.B + d.R * rp.e₄y +
           d.L1 * rp.e₂y / (2 * rp.nme) +
           (1 - xi2) * d.L2 * rp.q₂y / rp.nme

      th = xi2 * (π / 4 - d.tht1 - d.tht2) + d.tht2

   elseif r == 3

      Xx = d.A + d.R * xi2 / 4 + d.R * (1 - xi2) * rp.e₄x
      Xy = d.B + d.R * xi2 / 4 + d.R * (1 - xi2) * rp.e₄y

      th = xi2 * d.tht1 + π / 4 - d.tht1

   elseif r == 4

      Xx = d.A + d.R * (1 - xi2) / 4 + d.R * xi2 * rp.e₃x
      Xy = d.B + d.R * (1 - xi2) / 4 + d.R * xi2 * rp.e₃y

      th = xi2 * d.tht1 + π / 4

   elseif r == 5

      Xx = d.A + d.R * rp.e₃x +
           d.L1 * rp.e₁x / (2 * rp.nme) +
           xi2 * d.L2 * rp.q₁x / rp.nme

      Xy = d.B + d.R * rp.e₃y +
           d.L1 * rp.e₁y / (2 * rp.nme) +
           xi2 * d.L2 * rp.q₁y / rp.nme

      th = xi2 * (π / 4 - d.tht1 - d.tht2) + π / 4 + d.tht1

   elseif r == 6

      Xx = d.A + d.R * rp.e₃x +
           (1 - 2xi2) * d.L1 * rp.e₁x / (2 * rp.nme) +
           d.L2 * rp.q₁x / rp.nme

      Xy = d.B + d.R * rp.e₃y +
           (1 - 2xi2) * d.L1 * rp.e₁y / (2 * rp.nme) +
           d.L2 * rp.q₁y / rp.nme

      th = xi2 * (d.tht2 + d.tht3) + π / 2 - d.tht2

   elseif r == 7

      Xx = d.A + d.R * rp.e₃x -
           d.L1 * rp.e₁x / (2 * rp.nme) +
           (1 - xi2) * d.L2 * rp.q₁x / rp.nme

      Xy = d.B + d.R * rp.e₃y -
           d.L1 * rp.e₁y / (2 * rp.nme) +
           (1 - xi2) * d.L2 * rp.q₁y / rp.nme

      th = xi2 * (π / 4 - d.tht1 - d.tht3) + π / 2 + d.tht3

   elseif r == 8

      Xx = d.A + d.R * xi2 / 4 + d.R * (1 - xi2) * rp.e₃x
      Xy = d.B + d.R * xi2 / 4 + d.R * (1 - xi2) * rp.e₃y

      th = xi2 * d.tht1 + 3π / 4 - d.tht1

   elseif r == 9

      Xx = d.A + d.R * (1 - xi2) / 4 + d.R * xi2 * rp.e₄x
      Xy = d.B + d.R * (1 - xi2) / 4 + d.R * xi2 * rp.e₄y

      th = xi2 * d.tht1 + 3π / 4

   elseif r == 10

      Xx = d.A + d.R * rp.e₄x -
           d.L1 * rp.e₂x / (2 * rp.nme) +
           xi2 * d.L2 * rp.q₂x / rp.nme

      Xy = d.B + d.R * rp.e₄y -
           d.L1 * rp.e₂y / (2 * rp.nme) +
           xi2 * d.L2 * rp.q₂y / rp.nme

      th = xi2 * (π / 4 - d.tht1 - d.tht3) + 3π / 4 + d.tht1

   else

      throw(ArgumentError("mapxy expects region 1–12; got reg=$r"))

   end

   st, ct = sincos(th)
   h = st^3 + ct^3

   Yx = d.A + d.R * h * ct
   Yy = d.B + d.R * h * st

   Zx = muladd(xi1, Yx - Xx, Xx)
   Zy = muladd(xi1, Yy - Xy, Xy)

   return Zx, Zy

end

function mapxy!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
   d::bean, u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)

   p = d.pths[k]
   rp = d.RP

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   αu = hc / 2
   αv = ht / 2
   βu = p.ck0 + αu
   βv = p.tk0 + αv

   r = p.reg

   if r == 1

      ak = d.tht2 + d.tht3
      bk = -d.tht3

      @inbounds for I in eachindex(u, v, Zx, Zy)
         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + d.R * rp.e₄x +
              (2xi2 - 1) * d.L1 * rp.e₂x / (2 * rp.nme) +
              d.L2 * rp.q₂x / rp.nme

         Xy = d.B + d.R * rp.e₄y +
              (2xi2 - 1) * d.L1 * rp.e₂y / (2 * rp.nme) +
              d.L2 * rp.q₂y / rp.nme

         th = muladd(ak, xi2, bk)
         st, ct = sincos(th)
         h = st^3 + ct^3

         Yx = d.A + d.R * h * ct
         Yy = d.B + d.R * h * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)
      end

   elseif r == 2

      ak = π / 4 - d.tht1 - d.tht2
      bk = d.tht2

      @inbounds for I in eachindex(u, v, Zx, Zy)
         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + d.R * rp.e₄x +
              d.L1 * rp.e₂x / (2 * rp.nme) +
              (1 - xi2) * d.L2 * rp.q₂x / rp.nme

         Xy = d.B + d.R * rp.e₄y +
              d.L1 * rp.e₂y / (2 * rp.nme) +
              (1 - xi2) * d.L2 * rp.q₂y / rp.nme

         th = muladd(ak, xi2, bk)
         st, ct = sincos(th)
         h = st^3 + ct^3

         Yx = d.A + d.R * h * ct
         Yy = d.B + d.R * h * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)
      end

   elseif r == 3

      ak = d.tht1
      bk = π / 4 - d.tht1

      @inbounds for I in eachindex(u, v, Zx, Zy)
         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + d.R * xi2 / 4 + d.R * (1 - xi2) * rp.e₄x
         Xy = d.B + d.R * xi2 / 4 + d.R * (1 - xi2) * rp.e₄y

         th = muladd(ak, xi2, bk)
         st, ct = sincos(th)
         h = st^3 + ct^3

         Yx = d.A + d.R * h * ct
         Yy = d.B + d.R * h * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)
      end

   elseif r == 4

      ak = d.tht1
      bk = π / 4

      @inbounds for I in eachindex(u, v, Zx, Zy)
         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + d.R * (1 - xi2) / 4 + d.R * xi2 * rp.e₃x
         Xy = d.B + d.R * (1 - xi2) / 4 + d.R * xi2 * rp.e₃y

         th = muladd(ak, xi2, bk)
         st, ct = sincos(th)
         h = st^3 + ct^3

         Yx = d.A + d.R * h * ct
         Yy = d.B + d.R * h * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)
      end

   elseif r == 5

      ak = π / 4 - d.tht1 - d.tht2
      bk = π / 4 + d.tht1

      @inbounds for I in eachindex(u, v, Zx, Zy)
         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + d.R * rp.e₃x +
              d.L1 * rp.e₁x / (2 * rp.nme) +
              xi2 * d.L2 * rp.q₁x / rp.nme

         Xy = d.B + d.R * rp.e₃y +
              d.L1 * rp.e₁y / (2 * rp.nme) +
              xi2 * d.L2 * rp.q₁y / rp.nme

         th = muladd(ak, xi2, bk)
         st, ct = sincos(th)
         h = st^3 + ct^3

         Yx = d.A + d.R * h * ct
         Yy = d.B + d.R * h * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)
      end

   elseif r == 6

      ak = d.tht2 + d.tht3
      bk = π / 2 - d.tht2

      @inbounds for I in eachindex(u, v, Zx, Zy)
         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + d.R * rp.e₃x +
              (1 - 2xi2) * d.L1 * rp.e₁x / (2 * rp.nme) +
              d.L2 * rp.q₁x / rp.nme

         Xy = d.B + d.R * rp.e₃y +
              (1 - 2xi2) * d.L1 * rp.e₁y / (2 * rp.nme) +
              d.L2 * rp.q₁y / rp.nme

         th = muladd(ak, xi2, bk)
         st, ct = sincos(th)
         h = st^3 + ct^3

         Yx = d.A + d.R * h * ct
         Yy = d.B + d.R * h * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)
      end

   elseif r == 7

      ak = π / 4 - d.tht1 - d.tht3
      bk = π / 2 + d.tht3

      @inbounds for I in eachindex(u, v, Zx, Zy)
         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + d.R * rp.e₃x -
              d.L1 * rp.e₁x / (2 * rp.nme) +
              (1 - xi2) * d.L2 * rp.q₁x / rp.nme

         Xy = d.B + d.R * rp.e₃y -
              d.L1 * rp.e₁y / (2 * rp.nme) +
              (1 - xi2) * d.L2 * rp.q₁y / rp.nme

         th = muladd(ak, xi2, bk)
         st, ct = sincos(th)
         h = st^3 + ct^3

         Yx = d.A + d.R * h * ct
         Yy = d.B + d.R * h * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)
      end

   elseif r == 8

      ak = d.tht1
      bk = 3π / 4 - d.tht1

      @inbounds for I in eachindex(u, v, Zx, Zy)
         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + d.R * xi2 / 4 + d.R * (1 - xi2) * rp.e₃x
         Xy = d.B + d.R * xi2 / 4 + d.R * (1 - xi2) * rp.e₃y

         th = muladd(ak, xi2, bk)
         st, ct = sincos(th)
         h = st^3 + ct^3

         Yx = d.A + d.R * h * ct
         Yy = d.B + d.R * h * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)
      end

   elseif r == 9

      ak = d.tht1
      bk = 3π / 4

      @inbounds for I in eachindex(u, v, Zx, Zy)
         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + d.R * (1 - xi2) / 4 + d.R * xi2 * rp.e₄x
         Xy = d.B + d.R * (1 - xi2) / 4 + d.R * xi2 * rp.e₄y

         th = muladd(ak, xi2, bk)
         st, ct = sincos(th)
         h = st^3 + ct^3

         Yx = d.A + d.R * h * ct
         Yy = d.B + d.R * h * st

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)
      end

   elseif r == 10

      ak = π / 4 - d.tht1 - d.tht3
      bk = 3π / 4 + d.tht1

      @inbounds for I in eachindex(u, v, Zx, Zy)
         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + d.R * rp.e₄x -
              d.L1 * rp.e₂x / (2 * rp.nme) +
              xi2 * d.L2 * rp.q₂x / rp.nme

         Xy = d.B + d.R * rp.e₄y -
              d.L1 * rp.e₂y / (2 * rp.nme) +
              xi2 * d.L2 * rp.q₂y / rp.nme

         th = muladd(ak, xi2, bk)
         st, ct = sincos(th)
         h = st^3 + ct^3

         Yx = d.A + d.R * h * ct
         Yy = d.B + d.R * h * st

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

      if r == 11

         @inbounds for I in eachindex(u, v, Zx, Zy)
            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            Zx[I] = d.A + d.R * rp.e₄x +
                    (2xi2 - 1) * d.L1 * rp.e₂x / (2 * rp.nme) +
                    d.L2 * xi1 * rp.q₂x / rp.nme

            Zy[I] = d.B + d.R * rp.e₄y +
                    (2xi2 - 1) * d.L1 * rp.e₂y / (2 * rp.nme) +
                    d.L2 * xi1 * rp.q₂y / rp.nme
         end

      else # r == 12

         @inbounds for I in eachindex(u, v, Zx, Zy)
            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            Zx[I] = d.A + d.R * rp.e₃x +
                    (2xi2 - 1) * d.L1 * rp.e₁x / (2 * rp.nme) +
                    d.L2 * xi1 * rp.q₁x / rp.nme

            Zy[I] = d.B + d.R * rp.e₃y +
                    (2xi2 - 1) * d.L1 * rp.e₁y / (2 * rp.nme) +
                    d.L2 * xi1 * rp.q₁y / rp.nme
         end

      end

   else

      throw(ArgumentError("mapxy! expects region 1–12; got reg=$r"))

   end

   return nothing

end

function draw(d::bean, flag = nothing; L::Int = 33, show::Bool = true)

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
  gamx(d::bean, t::Float64, k::Int)
  gamx!(out::StridedArray{Float64}, d::bean, t::StridedArray{Float64}, k)
First (x) coordinate of the boundary parametrization γ(t).

- `gamx(d, t, k)` computes the **right boundary** of patch `k`, where `t ∈ [-1,1]`.

`t` may be a scalar or an array; the return has the same shape.
"""
function gamx(d::bean, t::Float64, k::Int)::Float64

   p = d.pths[k]
   r = p.reg

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt
   ξ = muladd(αt, t, βt)

   θ = if r == 1
      ξ * (d.tht2 + d.tht3) - d.tht3

   elseif r == 2
      ξ * (π / 4 - d.tht1 - d.tht2) + d.tht2

   elseif r == 3
      ξ * d.tht1 + π / 4 - d.tht1

   elseif r == 4
      ξ * d.tht1 + π / 4

   elseif r == 5
      ξ * (π / 4 - d.tht1 - d.tht2) + π / 4 + d.tht1

   elseif r == 6
      ξ * (d.tht2 + d.tht3) + π / 2 - d.tht2

   elseif r == 7
      ξ * (π / 4 - d.tht1 - d.tht3) + π / 2 + d.tht3

   elseif r == 8
      ξ * d.tht1 + 3π / 4 - d.tht1

   elseif r == 9
      ξ * d.tht1 + 3π / 4

   elseif r == 10
      ξ * (π / 4 - d.tht1 - d.tht3) + 3π / 4 + d.tht1

   else
      throw(ArgumentError("gamx for bean is defined only for boundary regions 1–10; got reg=$r"))
   end

   st, ct = sincos(θ)
   h = st^3 + ct^3

   return d.A + d.R * h * ct

end

function gamx!(out::StridedArray{Float64}, d::bean,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]
   r = p.reg

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   ak, bk = if r == 1
      d.tht2 + d.tht3, -d.tht3
   elseif r == 2
      π / 4 - d.tht1 - d.tht2, d.tht2
   elseif r == 3
      d.tht1, π / 4 - d.tht1
   elseif r == 4
      d.tht1, π / 4
   elseif r == 5
      π / 4 - d.tht1 - d.tht2, π / 4 + d.tht1
   elseif r == 6
      d.tht2 + d.tht3, π / 2 - d.tht2
   elseif r == 7
      π / 4 - d.tht1 - d.tht3, π / 2 + d.tht3
   elseif r == 8
      d.tht1, 3π / 4 - d.tht1
   elseif r == 9
      d.tht1, 3π / 4
   elseif r == 10
      π / 4 - d.tht1 - d.tht3, 3π / 4 + d.tht1
   else
      throw(ArgumentError("gamx! for bean is defined only for boundary regions 1–10; got reg=$r"))
   end

   @inbounds for I in eachindex(t, out)
      ξ = muladd(αt, t[I], βt)
      θ = muladd(ak, ξ, bk)

      st, ct = sincos(θ)
      h = st^3 + ct^3

      out[I] = d.A + d.R * h * ct
   end

   return nothing

end
"""
    gamy(d::bean, t::Float64, k::Int)
    gamy!(out::StridedArray{Float64}, d::bean, t::StridedArray{Float64}, k)

Second (y) coordinate of the boundary parametrization γ(t).

- `gamy(d, t, k)` computes the right boundary of patch `k`, where `t ∈ [-1,1]`.

`t` can be a scalar `Float64` or an array; the return has the same shape.
"""
function gamy(d::bean, t::Float64, k::Int)::Float64

   p = d.pths[k]
   r = p.reg

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt
   ξ = muladd(αt, t, βt)

   θ = if r == 1
      ξ * (d.tht2 + d.tht3) - d.tht3

   elseif r == 2
      ξ * (π / 4 - d.tht1 - d.tht2) + d.tht2

   elseif r == 3
      ξ * d.tht1 + π / 4 - d.tht1

   elseif r == 4
      ξ * d.tht1 + π / 4

   elseif r == 5
      ξ * (π / 4 - d.tht1 - d.tht2) + π / 4 + d.tht1

   elseif r == 6
      ξ * (d.tht2 + d.tht3) + π / 2 - d.tht2

   elseif r == 7
      ξ * (π / 4 - d.tht1 - d.tht3) + π / 2 + d.tht3

   elseif r == 8
      ξ * d.tht1 + 3π / 4 - d.tht1

   elseif r == 9
      ξ * d.tht1 + 3π / 4

   elseif r == 10
      ξ * (π / 4 - d.tht1 - d.tht3) + 3π / 4 + d.tht1

   else
      throw(ArgumentError("gamy for bean is defined only for boundary regions 1–10; got reg=$r"))
   end

   st, ct = sincos(θ)
   h = st^3 + ct^3

   return d.B + d.R * h * st

end

function gamy!(out::StridedArray{Float64}, d::bean,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]
   r = p.reg

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   ak, bk = if r == 1
      d.tht2 + d.tht3, -d.tht3
   elseif r == 2
      π / 4 - d.tht1 - d.tht2, d.tht2
   elseif r == 3
      d.tht1, π / 4 - d.tht1
   elseif r == 4
      d.tht1, π / 4
   elseif r == 5
      π / 4 - d.tht1 - d.tht2, π / 4 + d.tht1
   elseif r == 6
      d.tht2 + d.tht3, π / 2 - d.tht2
   elseif r == 7
      π / 4 - d.tht1 - d.tht3, π / 2 + d.tht3
   elseif r == 8
      d.tht1, 3π / 4 - d.tht1
   elseif r == 9
      d.tht1, 3π / 4
   elseif r == 10
      π / 4 - d.tht1 - d.tht3, 3π / 4 + d.tht1
   else
      throw(ArgumentError("gamy! for bean is defined only for boundary regions 1–10; got reg=$r"))
   end

   @inbounds for I in eachindex(t, out)
      ξ = muladd(αt, t[I], βt)
      θ = muladd(ak, ξ, bk)

      st, ct = sincos(θ)
      h = st^3 + ct^3

      out[I] = d.B + d.R * h * st
   end

   return nothing

end

function gam(d::bean, t::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]
   r = p.reg

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt
   ξ = muladd(αt, t, βt)

   θ = if r == 1
      ξ * (d.tht2 + d.tht3) - d.tht3

   elseif r == 2
      ξ * (π / 4 - d.tht1 - d.tht2) + d.tht2

   elseif r == 3
      ξ * d.tht1 + π / 4 - d.tht1

   elseif r == 4
      ξ * d.tht1 + π / 4

   elseif r == 5
      ξ * (π / 4 - d.tht1 - d.tht2) + π / 4 + d.tht1

   elseif r == 6
      ξ * (d.tht2 + d.tht3) + π / 2 - d.tht2

   elseif r == 7
      ξ * (π / 4 - d.tht1 - d.tht3) + π / 2 + d.tht3

   elseif r == 8
      ξ * d.tht1 + 3π / 4 - d.tht1

   elseif r == 9
      ξ * d.tht1 + 3π / 4

   elseif r == 10
      ξ * (π / 4 - d.tht1 - d.tht3) + 3π / 4 + d.tht1

   else
      throw(ArgumentError("gam for bean is defined only for boundary regions 1–10; got reg=$r"))
   end

   st, ct = sincos(θ)
   h = st^3 + ct^3

   gx = d.A + d.R * h * ct
   gy = d.B + d.R * h * st

   return gx, gy

end

function gam!(out::Vector{Float64}, d::bean, t::Float64, k::Int)

   p = d.pths[k]
   r = p.reg

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt
   ξ = muladd(αt, t, βt)

   ak, bk = if r == 1
      d.tht2 + d.tht3, -d.tht3
   elseif r == 2
      π / 4 - d.tht1 - d.tht2, d.tht2
   elseif r == 3
      d.tht1, π / 4 - d.tht1
   elseif r == 4
      d.tht1, π / 4
   elseif r == 5
      π / 4 - d.tht1 - d.tht2, π / 4 + d.tht1
   elseif r == 6
      d.tht2 + d.tht3, π / 2 - d.tht2
   elseif r == 7
      π / 4 - d.tht1 - d.tht3, π / 2 + d.tht3
   elseif r == 8
      d.tht1, 3π / 4 - d.tht1
   elseif r == 9
      d.tht1, 3π / 4
   elseif r == 10
      π / 4 - d.tht1 - d.tht3, 3π / 4 + d.tht1
   else
      throw(ArgumentError("gam! for bean is defined only for boundary regions 1–10; got reg=$r"))
   end

   θ = muladd(ak, ξ, bk)
   st, ct = sincos(θ)
   h = st^3 + ct^3

   out[1] = d.A + d.R * h * ct
   out[2] = d.B + d.R * h * st

   return nothing

end

function gam!(out::Matrix{Float64}, d::bean, t::Vector{Float64}, k::Int)

   p = d.pths[k]
   r = p.reg

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   ak, bk = if r == 1
      d.tht2 + d.tht3, -d.tht3

   elseif r == 2
      π / 4 - d.tht1 - d.tht2, d.tht2

   elseif r == 3
      d.tht1, π / 4 - d.tht1

   elseif r == 4
      d.tht1, π / 4

   elseif r == 5
      π / 4 - d.tht1 - d.tht2, π / 4 + d.tht1

   elseif r == 6
      d.tht2 + d.tht3, π / 2 - d.tht2

   elseif r == 7
      π / 4 - d.tht1 - d.tht3, π / 2 + d.tht3

   elseif r == 8
      d.tht1, 3π / 4 - d.tht1

   elseif r == 9
      d.tht1, 3π / 4

   elseif r == 10
      π / 4 - d.tht1 - d.tht3, 3π / 4 + d.tht1

   else
      throw(ArgumentError("gam! for bean is defined only for boundary regions 1–10; got reg=$r"))
   end

   @inbounds for I in eachindex(t)
      ξ = muladd(αt, t[I], βt)
      θ = muladd(ak, ξ, bk)

      st, ct = sincos(θ)
      h = st^3 + ct^3

      out[1, I] = d.A + d.R * h * ct
      out[2, I] = d.B + d.R * h * st
   end

   return nothing

end

function drawbd(d::bean, flag = true; L::Int = 33, show::Bool = true)

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

    return drawbd_geom(d, colors; flag = flag, L = L, show = show)
end

"""
    dgamx(d::bean, t::Float64, k::Int) -> Float64

Derivative of the first coordinate of the boundary parametrization.

With `k`: derivative of the right boundary of patch `k`,
i.e. `∂/∂t mapx(d, 1, t, k)` for `t ∈ [-1,1]`.
"""
function dgamx(d::bean, t::Float64, k::Int)::Float64

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = ht / 2
   βt = p.tk0 + αt

   if p.reg == 1
      Δth = d.tht2 + d.tht3
      th0 = -d.tht3

   elseif p.reg == 2
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = d.tht2

   elseif p.reg == 3
      Δth = d.tht1
      th0 = π / 4 - d.tht1

   elseif p.reg == 4
      Δth = d.tht1
      th0 = π / 4

   elseif p.reg == 5
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = π / 4 + d.tht1

   elseif p.reg == 6
      Δth = d.tht2 + d.tht3
      th0 = π / 2 - d.tht2

   elseif p.reg == 7
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = π / 2 + d.tht3

   elseif p.reg == 8
      Δth = d.tht1
      th0 = 3π / 4 - d.tht1

   elseif p.reg == 9
      Δth = d.tht1
      th0 = 3π / 4

   elseif p.reg == 10
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = 3π / 4 + d.tht1

   else
      throw(ArgumentError("dgamx is defined only for regions 1–10; got reg=$(p.reg)"))
   end

   ξ = muladd(αt, t, βt)
   θ = muladd(Δth, ξ, th0)

   s, c = sincos(θ)

   h = s^3 + c^3
   h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

   scale = 0.5 * ht * Δth

   return scale * d.R * (h1 * c - h * s)

end

function dgamx!(out::StridedArray{Float64}, d::bean,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = ht / 2
   βt = p.tk0 + αt

   if p.reg == 1
      Δth = d.tht2 + d.tht3
      th0 = -d.tht3

   elseif p.reg == 2
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = d.tht2

   elseif p.reg == 3
      Δth = d.tht1
      th0 = π / 4 - d.tht1

   elseif p.reg == 4
      Δth = d.tht1
      th0 = π / 4

   elseif p.reg == 5
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = π / 4 + d.tht1

   elseif p.reg == 6
      Δth = d.tht2 + d.tht3
      th0 = π / 2 - d.tht2

   elseif p.reg == 7
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = π / 2 + d.tht3

   elseif p.reg == 8
      Δth = d.tht1
      th0 = 3π / 4 - d.tht1

   elseif p.reg == 9
      Δth = d.tht1
      th0 = 3π / 4

   elseif p.reg == 10
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = 3π / 4 + d.tht1

   else
      throw(ArgumentError("dgamx! is defined only for regions 1–10; got reg=$(p.reg)"))
   end

   scale = 0.5 * ht * Δth * d.R

   @inbounds for I in eachindex(out, t)

      ξ = muladd(αt, t[I], βt)
      θ = muladd(Δth, ξ, th0)

      s, c = sincos(θ)

      h = s^3 + c^3
      h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

      out[I] = scale * (h1 * c - h * s)

   end

   return nothing

end

"""
    dgamy(d::bean, t::Float64, k::Int) -> Float64

Derivative of the second coordinate of the boundary parametrization.

With `k`: derivative of the right boundary of patch `k`,
i.e. `∂/∂t mapy(d, 1, t, k)` for `t ∈ [-1,1]`.
"""
function dgamy(d::bean, t::Float64, k::Int)::Float64

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = ht / 2
   βt = p.tk0 + αt

   if p.reg == 1
      Δth = d.tht2 + d.tht3
      th0 = -d.tht3

   elseif p.reg == 2
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = d.tht2

   elseif p.reg == 3
      Δth = d.tht1
      th0 = π / 4 - d.tht1

   elseif p.reg == 4
      Δth = d.tht1
      th0 = π / 4

   elseif p.reg == 5
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = π / 4 + d.tht1

   elseif p.reg == 6
      Δth = d.tht2 + d.tht3
      th0 = π / 2 - d.tht2

   elseif p.reg == 7
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = π / 2 + d.tht3

   elseif p.reg == 8
      Δth = d.tht1
      th0 = 3π / 4 - d.tht1

   elseif p.reg == 9
      Δth = d.tht1
      th0 = 3π / 4

   elseif p.reg == 10
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = 3π / 4 + d.tht1

   else
      throw(ArgumentError("dgamy is defined only for regions 1–10; got reg=$(p.reg)"))
   end

   ξ = muladd(αt, t, βt)
   θ = muladd(Δth, ξ, th0)

   s, c = sincos(θ)

   h = s^3 + c^3
   h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

   scale = 0.5 * ht * Δth

   return scale * d.R * (h1 * s + h * c)

end

function dgamy!(out::StridedArray{Float64}, d::bean,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = ht / 2
   βt = p.tk0 + αt

   if p.reg == 1
      Δth = d.tht2 + d.tht3
      th0 = -d.tht3

   elseif p.reg == 2
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = d.tht2

   elseif p.reg == 3
      Δth = d.tht1
      th0 = π / 4 - d.tht1

   elseif p.reg == 4
      Δth = d.tht1
      th0 = π / 4

   elseif p.reg == 5
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = π / 4 + d.tht1

   elseif p.reg == 6
      Δth = d.tht2 + d.tht3
      th0 = π / 2 - d.tht2

   elseif p.reg == 7
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = π / 2 + d.tht3

   elseif p.reg == 8
      Δth = d.tht1
      th0 = 3π / 4 - d.tht1

   elseif p.reg == 9
      Δth = d.tht1
      th0 = 3π / 4

   elseif p.reg == 10
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = 3π / 4 + d.tht1

   else
      throw(ArgumentError("dgamy! is defined only for regions 1–10; got reg=$(p.reg)"))
   end

   scale = 0.5 * ht * Δth * d.R

   @inbounds for I in eachindex(out, t)

      ξ = muladd(αt, t[I], βt)
      θ = muladd(Δth, ξ, th0)

      s, c = sincos(θ)

      h = s^3 + c^3
      h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

      out[I] = scale * (h1 * s + h * c)

   end

   return nothing

end

function dgam(d::bean, t::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = ht / 2
   βt = p.tk0 + αt

   if p.reg == 1
      Δth = d.tht2 + d.tht3
      th0 = -d.tht3

   elseif p.reg == 2
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = d.tht2

   elseif p.reg == 3
      Δth = d.tht1
      th0 = π / 4 - d.tht1

   elseif p.reg == 4
      Δth = d.tht1
      th0 = π / 4

   elseif p.reg == 5
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = π / 4 + d.tht1

   elseif p.reg == 6
      Δth = d.tht2 + d.tht3
      th0 = π / 2 - d.tht2

   elseif p.reg == 7
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = π / 2 + d.tht3

   elseif p.reg == 8
      Δth = d.tht1
      th0 = 3π / 4 - d.tht1

   elseif p.reg == 9
      Δth = d.tht1
      th0 = 3π / 4

   elseif p.reg == 10
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = 3π / 4 + d.tht1

   else
      throw(ArgumentError("dgam is defined only for regions 1–10; got reg=$(p.reg)"))
   end

   ξ = muladd(αt, t, βt)
   θ = muladd(Δth, ξ, th0)

   s, c = sincos(θ)

   h = s^3 + c^3
   h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

   scale = 0.5 * ht * Δth * d.R

   dx = scale * (h1 * c - h * s)
   dy = scale * (h1 * s + h * c)

   return dx, dy

end

function dgam!(out::Vector{Float64}, d::bean, t::Float64, k::Int)

   dx, dy = dgam(d, t, k)

   out[1] = dx
   out[2] = dy

   return nothing

end

function dgam!(out::Matrix{Float64}, d::bean,
   t::Vector{Float64}, k::Int)

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = ht / 2
   βt = p.tk0 + αt

   if p.reg == 1
      Δth = d.tht2 + d.tht3
      th0 = -d.tht3

   elseif p.reg == 2
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = d.tht2

   elseif p.reg == 3
      Δth = d.tht1
      th0 = π / 4 - d.tht1

   elseif p.reg == 4
      Δth = d.tht1
      th0 = π / 4

   elseif p.reg == 5
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = π / 4 + d.tht1

   elseif p.reg == 6
      Δth = d.tht2 + d.tht3
      th0 = π / 2 - d.tht2

   elseif p.reg == 7
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = π / 2 + d.tht3

   elseif p.reg == 8
      Δth = d.tht1
      th0 = 3π / 4 - d.tht1

   elseif p.reg == 9
      Δth = d.tht1
      th0 = 3π / 4

   elseif p.reg == 10
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = 3π / 4 + d.tht1

   else
      throw(ArgumentError("dgam! is defined only for regions 1–10; got reg=$(p.reg)"))
   end

   scale = 0.5 * ht * Δth * d.R

   @inbounds for I in eachindex(t)

      ξ = muladd(αt, t[I], βt)
      θ = muladd(Δth, ξ, th0)

      s, c = sincos(θ)

      h = s^3 + c^3
      h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

      out[1, I] = scale * (h1 * c - h * s)
      out[2, I] = scale * (h1 * s + h * c)

   end

   return nothing

end

function gamp(d::bean, t::Float64, k::Int)::Tuple{Float64,Float64}

   dx, dy = dgam(d, t, k)

   return dy, -dx

end

function gamp!(out::Matrix{Float64}, d::bean,
   t::Vector{Float64}, k::Int)

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = ht / 2
   βt = p.tk0 + αt

   if p.reg == 1
      Δth = d.tht2 + d.tht3
      th0 = -d.tht3

   elseif p.reg == 2
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = d.tht2

   elseif p.reg == 3
      Δth = d.tht1
      th0 = π / 4 - d.tht1

   elseif p.reg == 4
      Δth = d.tht1
      th0 = π / 4

   elseif p.reg == 5
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = π / 4 + d.tht1

   elseif p.reg == 6
      Δth = d.tht2 + d.tht3
      th0 = π / 2 - d.tht2

   elseif p.reg == 7
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = π / 2 + d.tht3

   elseif p.reg == 8
      Δth = d.tht1
      th0 = 3π / 4 - d.tht1

   elseif p.reg == 9
      Δth = d.tht1
      th0 = 3π / 4

   elseif p.reg == 10
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = 3π / 4 + d.tht1

   else
      throw(ArgumentError("gamp! is defined only for regions 1–10; got reg=$(p.reg)"))
   end

   scale = 0.5 * ht * Δth * d.R

   @inbounds for I in eachindex(t)

      ξ = muladd(αt, t[I], βt)
      θ = muladd(Δth, ξ, th0)

      s, c = sincos(θ)

      h = s^3 + c^3
      h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

      dx = scale * (h1 * c - h * s)
      dy = scale * (h1 * s + h * c)

      out[1, I] = dy
      out[2, I] = -dx

   end

   return nothing

end

function nu!(out::Matrix{Float64}, d::bean,
   t::Vector{Float64}, k::Int)

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = ht / 2
   βt = p.tk0 + αt

   if p.reg == 1
      Δth = d.tht2 + d.tht3
      th0 = -d.tht3

   elseif p.reg == 2
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = d.tht2

   elseif p.reg == 3
      Δth = d.tht1
      th0 = π / 4 - d.tht1

   elseif p.reg == 4
      Δth = d.tht1
      th0 = π / 4

   elseif p.reg == 5
      Δth = π / 4 - d.tht1 - d.tht2
      th0 = π / 4 + d.tht1

   elseif p.reg == 6
      Δth = d.tht2 + d.tht3
      th0 = π / 2 - d.tht2

   elseif p.reg == 7
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = π / 2 + d.tht3

   elseif p.reg == 8
      Δth = d.tht1
      th0 = 3π / 4 - d.tht1

   elseif p.reg == 9
      Δth = d.tht1
      th0 = 3π / 4

   elseif p.reg == 10
      Δth = π / 4 - d.tht1 - d.tht3
      th0 = 3π / 4 + d.tht1

   else
      throw(ArgumentError("nu! is defined only for regions 1–10; got reg=$(p.reg)"))
   end

   scale = 0.5 * ht * Δth * d.R

   @inbounds for I in eachindex(t)

      ξ = muladd(αt, t[I], βt)
      θ = muladd(Δth, ξ, th0)

      s, c = sincos(θ)

      h = s^3 + c^3
      h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

      dx = scale * (h1 * c - h * s)
      dy = scale * (h1 * s + h * c)

      S = hypot(dx, dy)

      out[1, I] = dy / S
      out[2, I] = -dx / S

   end

   return nothing

end

# This function finds s such that γ_l(t) = γ_k(s),
# allowing s to lie outside [-1, 1].
function bdinv(d::bean, t::Float64, l::Int, k::Int)::Float64

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
         xil * (d.tht2 + d.tht3) - d.tht3

      elseif plr == 2
         xil * (π / 4 - d.tht1 - d.tht2) + d.tht2

      elseif plr == 3
         xil * d.tht1 + π / 4 - d.tht1

      elseif plr == 4
         xil * d.tht1 + π / 4

      elseif plr == 5
         xil * (π / 4 - d.tht1 - d.tht2) + π / 4 + d.tht1

      elseif plr == 6
         xil * (d.tht2 + d.tht3) + π / 2 - d.tht2

      elseif plr == 7
         xil * (π / 4 - d.tht1 - d.tht3) + π / 2 + d.tht3

      elseif plr == 8
         xil * d.tht1 + 3π / 4 - d.tht1

      elseif plr == 9
         xil * d.tht1 + 3π / 4

      elseif plr == 10
         xil * (π / 4 - d.tht1 - d.tht3) + 3π / 4 + d.tht1

      else
         throw(ArgumentError("bdinv source region must be 1–10; got reg=$plr"))
      end

      # Get angular map for target patch k:
      #     theta = ak * xi + bk
      ak, bk = if pkr == 1
         d.tht2 + d.tht3, -d.tht3

      elseif pkr == 2
         π / 4 - d.tht1 - d.tht2, d.tht2

      elseif pkr == 3
         d.tht1, π / 4 - d.tht1

      elseif pkr == 4
         d.tht1, π / 4

      elseif pkr == 5
         π / 4 - d.tht1 - d.tht2, π / 4 + d.tht1

      elseif pkr == 6
         d.tht2 + d.tht3, π / 2 - d.tht2

      elseif pkr == 7
         π / 4 - d.tht1 - d.tht3, π / 2 + d.tht3

      elseif pkr == 8
         d.tht1, 3π / 4 - d.tht1

      elseif pkr == 9
         d.tht1, 3π / 4

      elseif pkr == 10
         π / 4 - d.tht1 - d.tht3, 3π / 4 + d.tht1

      else
         throw(ArgumentError("bdinv target region must be 1–10; got reg=$pkr"))
      end

      # Shift th by multiples of 2π so it is closest to the target region.
      center = bk + 0.5 * ak
      m = round((center - th) / (2π))
      ths = th + 2π * m

      # Convert angle to target xi, then xi to local patch coordinate s.
      xik = (ths - bk) / ak

      return xi_inv(xik, pk.tk0, pk.tk1)

   end

end

"""
   hderhigher(d::bean, θ::Float64)

Compute h, h', h'', h⁽³⁾, h⁽⁴⁾, h⁽⁵⁾, h⁽⁶⁾ for the bean radial function

   h(θ) = sin(θ)^3 + cos(θ)^3

Returns

   h, h1, h2, h3, h4, h5, h6

with derivatives taken with respect to θ.
"""
function hderhigher(d::bean, θ::Float64)

   s, c = sincos(θ)

   spc = s + c
   smc = s - c
   sc = s * c

   h  = spc * (1.0 - sc)

   h1 = 3.0 * sc * smc

   h2 = -3.0 * spc * (1.0 - 3.0 * sc)

   h3 = -3.0 * smc * (2.0 + 9.0 * sc)

   h4 = 3.0 * spc * (7.0 - 27.0 * sc)

   h5 = 3.0 * smc * (20.0 + 81.0 * sc)

   h6 = -3.0 * spc * (61.0 - 243.0 * sc)

   return h, h1, h2, h3, h4, h5, h6

end

"""
    DLP!(out::StridedArray{Float64}, d::bean, t::Float64,
         tau::StridedArray{Float64}, k::Int,
         x::Vector{Float64}, G::Vector{Float64}, GP::Vector{Float64})

Same-panel double-layer kernel on the bean boundary:

    K(t, τ) =
    ((γ_k(τ) - γ_k(t)) ⋅ γᵖᵉʳᵖ_k(τ)) /
    ||γ_k(τ) - γ_k(t)||²

Uses the analytically simplified polar form where the removable sin²(Δ)
singularity is canceled over the whole panel.

The Taylor expansions for q and q₂ are used when Δ is small.
"""
function DLP!(out::StridedArray{Float64}, d::bean, t::Float64,
   tau::StridedArray{Float64}, k::Int, x::Vector{Float64},
   G::Vector{Float64}, GP::Vector{Float64})

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = p.tk0 + αt

   reg = p.reg

   ak, bk = if reg == 1
      d.tht2 + d.tht3, -d.tht3

   elseif reg == 2
      π / 4 - d.tht1 - d.tht2, d.tht2

   elseif reg == 3
      d.tht1, π / 4 - d.tht1

   elseif reg == 4
      d.tht1, π / 4

   elseif reg == 5
      π / 4 - d.tht1 - d.tht2, π / 4 + d.tht1

   elseif reg == 6
      d.tht2 + d.tht3, π / 2 - d.tht2

   elseif reg == 7
      π / 4 - d.tht1 - d.tht3, π / 2 + d.tht3

   elseif reg == 8
      d.tht1, 3π / 4 - d.tht1

   elseif reg == 9
      d.tht1, 3π / 4

   elseif reg == 10
      π / 4 - d.tht1 - d.tht3, 3π / 4 + d.tht1

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

      sτ, cτ = sincos(thτ)

      spcτ = sτ + cτ
      smcτ = sτ - cτ
      scτ = sτ * cτ

      hτ = spcτ * (1.0 - scτ)
      hτ1 = 3.0 * scτ * smcτ

      Δ = 0.5 * (thτ - tht)

      sΔ, cΔ = sincos(Δ)

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
    gamder(d::bean, θ::Float64)

Return centered bean boundary values and first derivatives with respect to `θ`.

Returns gx, gy, dgx, dgy where

    gx(θ) = R*h(θ)*cos(θ), gy(θ) = R*h(θ)*sin(θ)

and h(θ) = sin(θ)^3 + cos(θ)^3.
"""
function gamder(d::bean, θ::Float64)::Tuple{Float64,Float64,Float64,Float64}

   s, c = sincos(θ)

   h = s^3 + c^3
   h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

   gx = d.R * h * c
   gy = d.R * h * s

   dgx = d.R * (h1 * c - h * s)
   dgy = d.R * (h1 * s + h * c)

   return gx, gy, dgx, dgy

end

"""
    Dmap!(DJ::StridedArray{Float64}, d::bean,
          u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)

Fill `DJ` with the absolute determinant of the Jacobian of the bean patch map.

No allocations. `u`, `v`, and `DJ` must have compatible indexing.
"""
function Dmap!(DJ::StridedArray{Float64}, d::bean,
   u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)

   p = d.pths[k]
   rp = d.RP

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

   if r == 1

      ak = d.tht2 + d.tht3
      bk = -d.tht3

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * rp.e₄x +
              (2vhat - 1) * d.L1 * rp.e₂x / (2 * rp.nme) +
              d.L2 * rp.q₂x / rp.nme

         Xy = d.R * rp.e₄y +
              (2vhat - 1) * d.L1 * rp.e₂y / (2 * rp.nme) +
              d.L2 * rp.q₂y / rp.nme

         th = muladd(ak, vhat, bk)
         gx, gy, dgx, dgy = gamder(d, th)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * d.L1 * rp.e₂x / rp.nme +
               uhat * αt * ak * dgx

         dvy = (1 - uhat) * αt * d.L1 * rp.e₂y / rp.nme +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 2

      ak = π / 4 - d.tht1 - d.tht2
      bk = d.tht2

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * rp.e₄x +
              d.L1 * rp.e₂x / (2 * rp.nme) +
              (1 - vhat) * d.L2 * rp.q₂x / rp.nme

         Xy = d.R * rp.e₄y +
              d.L1 * rp.e₂y / (2 * rp.nme) +
              (1 - vhat) * d.L2 * rp.q₂y / rp.nme

         th = muladd(ak, vhat, bk)
         gx, gy, dgx, dgy = gamder(d, th)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = -(1 - uhat) * αt * d.L2 * rp.q₂x / rp.nme +
               uhat * αt * ak * dgx

         dvy = -(1 - uhat) * αt * d.L2 * rp.q₂y / rp.nme +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 3

      ak = d.tht1
      bk = π / 4 - d.tht1

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * vhat / 4 + d.R * (1 - vhat) * rp.e₄x
         Xy = d.R * vhat / 4 + d.R * (1 - vhat) * rp.e₄y

         th = muladd(ak, vhat, bk)
         gx, gy, dgx, dgy = gamder(d, th)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * d.R * (1 / 4 - rp.e₄x) +
               uhat * αt * ak * dgx

         dvy = (1 - uhat) * αt * d.R * (1 / 4 - rp.e₄y) +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 4

      ak = d.tht1
      bk = π / 4

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * (1 - vhat) / 4 + d.R * vhat * rp.e₃x
         Xy = d.R * (1 - vhat) / 4 + d.R * vhat * rp.e₃y

         th = muladd(ak, vhat, bk)
         gx, gy, dgx, dgy = gamder(d, th)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * d.R * (rp.e₃x - 1 / 4) +
               uhat * αt * ak * dgx

         dvy = (1 - uhat) * αt * d.R * (rp.e₃y - 1 / 4) +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 5

      ak = π / 4 - d.tht1 - d.tht2
      bk = π / 4 + d.tht1

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * rp.e₃x +
              d.L1 * rp.e₁x / (2 * rp.nme) +
              vhat * d.L2 * rp.q₁x / rp.nme

         Xy = d.R * rp.e₃y +
              d.L1 * rp.e₁y / (2 * rp.nme) +
              vhat * d.L2 * rp.q₁y / rp.nme

         th = muladd(ak, vhat, bk)
         gx, gy, dgx, dgy = gamder(d, th)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * d.L2 * rp.q₁x / rp.nme +
               uhat * αt * ak * dgx

         dvy = (1 - uhat) * αt * d.L2 * rp.q₁y / rp.nme +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 6

      ak = d.tht2 + d.tht3
      bk = π / 2 - d.tht2

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * rp.e₃x +
              (1 - 2vhat) * d.L1 * rp.e₁x / (2 * rp.nme) +
              d.L2 * rp.q₁x / rp.nme

         Xy = d.R * rp.e₃y +
              (1 - 2vhat) * d.L1 * rp.e₁y / (2 * rp.nme) +
              d.L2 * rp.q₁y / rp.nme

         th = muladd(ak, vhat, bk)
         gx, gy, dgx, dgy = gamder(d, th)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = -(1 - uhat) * αt * d.L1 * rp.e₁x / rp.nme +
               uhat * αt * ak * dgx

         dvy = -(1 - uhat) * αt * d.L1 * rp.e₁y / rp.nme +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 7

      ak = π / 4 - d.tht1 - d.tht3
      bk = π / 2 + d.tht3

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * rp.e₃x -
              d.L1 * rp.e₁x / (2 * rp.nme) +
              (1 - vhat) * d.L2 * rp.q₁x / rp.nme

         Xy = d.R * rp.e₃y -
              d.L1 * rp.e₁y / (2 * rp.nme) +
              (1 - vhat) * d.L2 * rp.q₁y / rp.nme

         th = muladd(ak, vhat, bk)
         gx, gy, dgx, dgy = gamder(d, th)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = -(1 - uhat) * αt * d.L2 * rp.q₁x / rp.nme +
               uhat * αt * ak * dgx

         dvy = -(1 - uhat) * αt * d.L2 * rp.q₁y / rp.nme +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 8

      ak = d.tht1
      bk = 3π / 4 - d.tht1

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * vhat / 4 + d.R * (1 - vhat) * rp.e₃x
         Xy = d.R * vhat / 4 + d.R * (1 - vhat) * rp.e₃y

         th = muladd(ak, vhat, bk)
         gx, gy, dgx, dgy = gamder(d, th)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * d.R * (1 / 4 - rp.e₃x) +
               uhat * αt * ak * dgx

         dvy = (1 - uhat) * αt * d.R * (1 / 4 - rp.e₃y) +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 9

      ak = d.tht1
      bk = 3π / 4

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * (1 - vhat) / 4 + d.R * vhat * rp.e₄x
         Xy = d.R * (1 - vhat) / 4 + d.R * vhat * rp.e₄y

         th = muladd(ak, vhat, bk)
         gx, gy, dgx, dgy = gamder(d, th)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * d.R * (rp.e₄x - 1 / 4) +
               uhat * αt * ak * dgx

         dvy = (1 - uhat) * αt * d.R * (rp.e₄y - 1 / 4) +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 10

      ak = π / 4 - d.tht1 - d.tht3
      bk = 3π / 4 + d.tht1

      @inbounds for I in eachindex(u, v, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * rp.e₄x -
              d.L1 * rp.e₂x / (2 * rp.nme) +
              vhat * d.L2 * rp.q₂x / rp.nme

         Xy = d.R * rp.e₄y -
              d.L1 * rp.e₂y / (2 * rp.nme) +
              vhat * d.L2 * rp.q₂y / rp.nme

         th = muladd(ak, vhat, bk)
         gx, gy, dgx, dgy = gamder(d, th)

         dux = αc * (gx - Xx)
         duy = αc * (gy - Xy)

         dvx = (1 - uhat) * αt * d.L2 * rp.q₂x / rp.nme +
               uhat * αt * ak * dgx

         dvy = (1 - uhat) * αt * d.L2 * rp.q₂y / rp.nme +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   else

      throw(ArgumentError("Dmap! for bean expects region 1–12; got reg=$r"))

   end

   return nothing

end

"""
A combination of mapxy! and Dmap! function (No allocations!)

The purpose of this function is to reduce repeated computations related
to cosine, sine, h, and h'.
"""
function mapxy_Dmap!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
   DJ::StridedArray{Float64}, d::bean, u::StridedArray{Float64},
   v::StridedArray{Float64}, k::Int)

   p = d.pths[k]
   rp = d.RP

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   αc = hc / 2
   αt = ht / 2
   βc = p.ck0 + αc
   βt = p.tk0 + αt

   r = p.reg
   invn = 1.0 / rp.nme

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

      detJ = d.L1 * d.L2 * hc * ht / 4

      if r == 11

         @inbounds for I in eachindex(u, v, Zx, Zy, DJ)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            Zx[I] = d.A + d.R * rp.e₄x +
                    (2xi2 - 1) * d.L1 * rp.e₂x * invn / 2 +
                    d.L2 * xi1 * rp.q₂x * invn

            Zy[I] = d.B + d.R * rp.e₄y +
                    (2xi2 - 1) * d.L1 * rp.e₂y * invn / 2 +
                    d.L2 * xi1 * rp.q₂y * invn

            DJ[I] = detJ

         end

      else # r == 12

         @inbounds for I in eachindex(u, v, Zx, Zy, DJ)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            Zx[I] = d.A + d.R * rp.e₃x +
                    (2xi2 - 1) * d.L1 * rp.e₁x * invn / 2 +
                    d.L2 * xi1 * rp.q₁x * invn

            Zy[I] = d.B + d.R * rp.e₃y +
                    (2xi2 - 1) * d.L1 * rp.e₁y * invn / 2 +
                    d.L2 * xi1 * rp.q₁y * invn

            DJ[I] = detJ

         end

      end

      return nothing

   end

   if r == 1

      ak = d.tht2 + d.tht3
      bk = -d.tht3

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * rp.e₄x +
              (2vhat - 1) * d.L1 * rp.e₂x * invn / 2 +
              d.L2 * rp.q₂x * invn

         Xy = d.R * rp.e₄y +
              (2vhat - 1) * d.L1 * rp.e₂y * invn / 2 +
              d.L2 * rp.q₂y * invn

         th = muladd(ak, vhat, bk)
         s, c = sincos(th)

         h = s^3 + c^3
         h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

         Yx = d.R * h * c
         Yy = d.R * h * s

         dgx = d.R * (h1 * c - h * s)
         dgy = d.R * (h1 * s + h * c)

         Zx[I] = d.A + muladd(uhat, Yx - Xx, Xx)
         Zy[I] = d.B + muladd(uhat, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = (1 - uhat) * αt * d.L1 * rp.e₂x * invn +
               uhat * αt * ak * dgx

         dvy = (1 - uhat) * αt * d.L1 * rp.e₂y * invn +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 2

      ak = π / 4 - d.tht1 - d.tht2
      bk = d.tht2

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * rp.e₄x +
              d.L1 * rp.e₂x * invn / 2 +
              (1 - vhat) * d.L2 * rp.q₂x * invn

         Xy = d.R * rp.e₄y +
              d.L1 * rp.e₂y * invn / 2 +
              (1 - vhat) * d.L2 * rp.q₂y * invn

         th = muladd(ak, vhat, bk)
         s, c = sincos(th)

         h = s^3 + c^3
         h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

         Yx = d.R * h * c
         Yy = d.R * h * s

         dgx = d.R * (h1 * c - h * s)
         dgy = d.R * (h1 * s + h * c)

         Zx[I] = d.A + muladd(uhat, Yx - Xx, Xx)
         Zy[I] = d.B + muladd(uhat, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = -(1 - uhat) * αt * d.L2 * rp.q₂x * invn +
               uhat * αt * ak * dgx

         dvy = -(1 - uhat) * αt * d.L2 * rp.q₂y * invn +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 3

      ak = d.tht1
      bk = π / 4 - d.tht1

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * vhat / 4 + d.R * (1 - vhat) * rp.e₄x
         Xy = d.R * vhat / 4 + d.R * (1 - vhat) * rp.e₄y

         th = muladd(ak, vhat, bk)
         s, c = sincos(th)

         h = s^3 + c^3
         h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

         Yx = d.R * h * c
         Yy = d.R * h * s

         dgx = d.R * (h1 * c - h * s)
         dgy = d.R * (h1 * s + h * c)

         Zx[I] = d.A + muladd(uhat, Yx - Xx, Xx)
         Zy[I] = d.B + muladd(uhat, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = (1 - uhat) * αt * d.R * (1 / 4 - rp.e₄x) +
               uhat * αt * ak * dgx

         dvy = (1 - uhat) * αt * d.R * (1 / 4 - rp.e₄y) +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 4

      ak = d.tht1
      bk = π / 4

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * (1 - vhat) / 4 + d.R * vhat * rp.e₃x
         Xy = d.R * (1 - vhat) / 4 + d.R * vhat * rp.e₃y

         th = muladd(ak, vhat, bk)
         s, c = sincos(th)

         h = s^3 + c^3
         h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

         Yx = d.R * h * c
         Yy = d.R * h * s

         dgx = d.R * (h1 * c - h * s)
         dgy = d.R * (h1 * s + h * c)

         Zx[I] = d.A + muladd(uhat, Yx - Xx, Xx)
         Zy[I] = d.B + muladd(uhat, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = (1 - uhat) * αt * d.R * (rp.e₃x - 1 / 4) +
               uhat * αt * ak * dgx

         dvy = (1 - uhat) * αt * d.R * (rp.e₃y - 1 / 4) +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 5

      ak = π / 4 - d.tht1 - d.tht2
      bk = π / 4 + d.tht1

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * rp.e₃x +
              d.L1 * rp.e₁x * invn / 2 +
              vhat * d.L2 * rp.q₁x * invn

         Xy = d.R * rp.e₃y +
              d.L1 * rp.e₁y * invn / 2 +
              vhat * d.L2 * rp.q₁y * invn

         th = muladd(ak, vhat, bk)
         s, c = sincos(th)

         h = s^3 + c^3
         h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

         Yx = d.R * h * c
         Yy = d.R * h * s

         dgx = d.R * (h1 * c - h * s)
         dgy = d.R * (h1 * s + h * c)

         Zx[I] = d.A + muladd(uhat, Yx - Xx, Xx)
         Zy[I] = d.B + muladd(uhat, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = (1 - uhat) * αt * d.L2 * rp.q₁x * invn +
               uhat * αt * ak * dgx

         dvy = (1 - uhat) * αt * d.L2 * rp.q₁y * invn +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 6

      ak = d.tht2 + d.tht3
      bk = π / 2 - d.tht2

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * rp.e₃x +
              (1 - 2vhat) * d.L1 * rp.e₁x * invn / 2 +
              d.L2 * rp.q₁x * invn

         Xy = d.R * rp.e₃y +
              (1 - 2vhat) * d.L1 * rp.e₁y * invn / 2 +
              d.L2 * rp.q₁y * invn

         th = muladd(ak, vhat, bk)
         s, c = sincos(th)

         h = s^3 + c^3
         h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

         Yx = d.R * h * c
         Yy = d.R * h * s

         dgx = d.R * (h1 * c - h * s)
         dgy = d.R * (h1 * s + h * c)

         Zx[I] = d.A + muladd(uhat, Yx - Xx, Xx)
         Zy[I] = d.B + muladd(uhat, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = -(1 - uhat) * αt * d.L1 * rp.e₁x * invn +
               uhat * αt * ak * dgx

         dvy = -(1 - uhat) * αt * d.L1 * rp.e₁y * invn +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 7

      ak = π / 4 - d.tht1 - d.tht3
      bk = π / 2 + d.tht3

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * rp.e₃x -
              d.L1 * rp.e₁x * invn / 2 +
              (1 - vhat) * d.L2 * rp.q₁x * invn

         Xy = d.R * rp.e₃y -
              d.L1 * rp.e₁y * invn / 2 +
              (1 - vhat) * d.L2 * rp.q₁y * invn

         th = muladd(ak, vhat, bk)
         s, c = sincos(th)

         h = s^3 + c^3
         h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

         Yx = d.R * h * c
         Yy = d.R * h * s

         dgx = d.R * (h1 * c - h * s)
         dgy = d.R * (h1 * s + h * c)

         Zx[I] = d.A + muladd(uhat, Yx - Xx, Xx)
         Zy[I] = d.B + muladd(uhat, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = -(1 - uhat) * αt * d.L2 * rp.q₁x * invn +
               uhat * αt * ak * dgx

         dvy = -(1 - uhat) * αt * d.L2 * rp.q₁y * invn +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 8

      ak = d.tht1
      bk = 3π / 4 - d.tht1

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * vhat / 4 + d.R * (1 - vhat) * rp.e₃x
         Xy = d.R * vhat / 4 + d.R * (1 - vhat) * rp.e₃y

         th = muladd(ak, vhat, bk)
         s, c = sincos(th)

         h = s^3 + c^3
         h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

         Yx = d.R * h * c
         Yy = d.R * h * s

         dgx = d.R * (h1 * c - h * s)
         dgy = d.R * (h1 * s + h * c)

         Zx[I] = d.A + muladd(uhat, Yx - Xx, Xx)
         Zy[I] = d.B + muladd(uhat, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = (1 - uhat) * αt * d.R * (1 / 4 - rp.e₃x) +
               uhat * αt * ak * dgx

         dvy = (1 - uhat) * αt * d.R * (1 / 4 - rp.e₃y) +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 9

      ak = d.tht1
      bk = 3π / 4

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * (1 - vhat) / 4 + d.R * vhat * rp.e₄x
         Xy = d.R * (1 - vhat) / 4 + d.R * vhat * rp.e₄y

         th = muladd(ak, vhat, bk)
         s, c = sincos(th)

         h = s^3 + c^3
         h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

         Yx = d.R * h * c
         Yy = d.R * h * s

         dgx = d.R * (h1 * c - h * s)
         dgy = d.R * (h1 * s + h * c)

         Zx[I] = d.A + muladd(uhat, Yx - Xx, Xx)
         Zy[I] = d.B + muladd(uhat, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = (1 - uhat) * αt * d.R * (rp.e₄x - 1 / 4) +
               uhat * αt * ak * dgx

         dvy = (1 - uhat) * αt * d.R * (rp.e₄y - 1 / 4) +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   elseif r == 10

      ak = π / 4 - d.tht1 - d.tht3
      bk = 3π / 4 + d.tht1

      @inbounds for I in eachindex(u, v, Zx, Zy, DJ)

         uhat = muladd(αc, u[I], βc)
         vhat = muladd(αt, v[I], βt)

         Xx = d.R * rp.e₄x -
              d.L1 * rp.e₂x * invn / 2 +
              vhat * d.L2 * rp.q₂x * invn

         Xy = d.R * rp.e₄y -
              d.L1 * rp.e₂y * invn / 2 +
              vhat * d.L2 * rp.q₂y * invn

         th = muladd(ak, vhat, bk)
         s, c = sincos(th)

         h = s^3 + c^3
         h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

         Yx = d.R * h * c
         Yy = d.R * h * s

         dgx = d.R * (h1 * c - h * s)
         dgy = d.R * (h1 * s + h * c)

         Zx[I] = d.A + muladd(uhat, Yx - Xx, Xx)
         Zy[I] = d.B + muladd(uhat, Yy - Xy, Xy)

         dux = αc * (Yx - Xx)
         duy = αc * (Yy - Xy)

         dvx = (1 - uhat) * αt * d.L2 * rp.q₂x * invn +
               uhat * αt * ak * dgx

         dvy = (1 - uhat) * αt * d.L2 * rp.q₂y * invn +
               uhat * αt * ak * dgy

         DJ[I] = abs(dux * dvy - dvx * duy)

      end

   else

      throw(ArgumentError("mapxy_Dmap! for bean expects region 1–12; got reg=$r"))

   end

   return nothing

end

function chk_map(d::bean; n::Int=32, tol::Float64=5e-14)
   Iex = d.R^2 * π * (d.R^2 * 123 / 1024 + 5 * (d.A^2 + d.B^2) / 16
         + 3 * (d.A + d.B) * d.R / 16)
   return chkmap_geom(d, Iex; n=n, tol=tol)
end

"""
    jinvmap(d::bean, u::Float64, v::Float64, r::Int)

Inverse Jacobian of the bean region mapping at a single reference point
`(u,v) ∈ [-1,1]^2`.

Returns `(J11, J12, J21, J22)`, the entries of `J⁻¹`.

Only the region index `r` is needed because this is used on the region-level
map, so `hc = ht = 1`.
"""
function jinvmap(d::bean, u::Float64, v::Float64, r::Int)

   rp = d.RP
   invn = 1.0 / rp.nme

   uhat = (u + 1) / 2
   vhat = (v + 1) / 2

   if r == 1

      ak = d.tht2 + d.tht3
      th = muladd(ak, vhat, -d.tht3)

      s, c = sincos(th)
      h = s^3 + c^3
      h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

      gx = d.R * h * c
      gy = d.R * h * s
      dgx = d.R * (h1 * c - h * s)
      dgy = d.R * (h1 * s + h * c)

      Xx = d.R * rp.e₄x +
           (2vhat - 1) * d.L1 * rp.e₂x * invn / 2 +
           d.L2 * rp.q₂x * invn

      Xy = d.R * rp.e₄y +
           (2vhat - 1) * d.L1 * rp.e₂y * invn / 2 +
           d.L2 * rp.q₂y * invn

      dux = (gx - Xx) / 2
      duy = (gy - Xy) / 2

      dvx = (1 - uhat) * d.L1 * rp.e₂x * invn / 2 +
            uhat * ak * dgx / 2

      dvy = (1 - uhat) * d.L1 * rp.e₂y * invn / 2 +
            uhat * ak * dgy / 2

   elseif r == 2

      ak = π / 4 - d.tht1 - d.tht2
      th = muladd(ak, vhat, d.tht2)

      s, c = sincos(th)
      h = s^3 + c^3
      h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

      gx = d.R * h * c
      gy = d.R * h * s
      dgx = d.R * (h1 * c - h * s)
      dgy = d.R * (h1 * s + h * c)

      Xx = d.R * rp.e₄x +
           d.L1 * rp.e₂x * invn / 2 +
           (1 - vhat) * d.L2 * rp.q₂x * invn

      Xy = d.R * rp.e₄y +
           d.L1 * rp.e₂y * invn / 2 +
           (1 - vhat) * d.L2 * rp.q₂y * invn

      dux = (gx - Xx) / 2
      duy = (gy - Xy) / 2

      dvx = -(1 - uhat) * d.L2 * rp.q₂x * invn / 2 +
            uhat * ak * dgx / 2

      dvy = -(1 - uhat) * d.L2 * rp.q₂y * invn / 2 +
            uhat * ak * dgy / 2

   elseif r == 3

      ak = d.tht1
      th = muladd(ak, vhat, π / 4 - d.tht1)

      s, c = sincos(th)
      h = s^3 + c^3
      h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

      gx = d.R * h * c
      gy = d.R * h * s
      dgx = d.R * (h1 * c - h * s)
      dgy = d.R * (h1 * s + h * c)

      Xx = d.R * vhat / 4 + d.R * (1 - vhat) * rp.e₄x
      Xy = d.R * vhat / 4 + d.R * (1 - vhat) * rp.e₄y

      dux = (gx - Xx) / 2
      duy = (gy - Xy) / 2

      dvx = (1 - uhat) * d.R * (1 / 4 - rp.e₄x) / 2 +
            uhat * ak * dgx / 2

      dvy = (1 - uhat) * d.R * (1 / 4 - rp.e₄y) / 2 +
            uhat * ak * dgy / 2

   elseif r == 4

      ak = d.tht1
      th = muladd(ak, vhat, π / 4)

      s, c = sincos(th)
      h = s^3 + c^3
      h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

      gx = d.R * h * c
      gy = d.R * h * s
      dgx = d.R * (h1 * c - h * s)
      dgy = d.R * (h1 * s + h * c)

      Xx = d.R * (1 - vhat) / 4 + d.R * vhat * rp.e₃x
      Xy = d.R * (1 - vhat) / 4 + d.R * vhat * rp.e₃y

      dux = (gx - Xx) / 2
      duy = (gy - Xy) / 2

      dvx = (1 - uhat) * d.R * (rp.e₃x - 1 / 4) / 2 +
            uhat * ak * dgx / 2

      dvy = (1 - uhat) * d.R * (rp.e₃y - 1 / 4) / 2 +
            uhat * ak * dgy / 2

   elseif r == 5

      ak = π / 4 - d.tht1 - d.tht2
      th = muladd(ak, vhat, π / 4 + d.tht1)

      s, c = sincos(th)
      h = s^3 + c^3
      h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

      gx = d.R * h * c
      gy = d.R * h * s
      dgx = d.R * (h1 * c - h * s)
      dgy = d.R * (h1 * s + h * c)

      Xx = d.R * rp.e₃x +
           d.L1 * rp.e₁x * invn / 2 +
           vhat * d.L2 * rp.q₁x * invn

      Xy = d.R * rp.e₃y +
           d.L1 * rp.e₁y * invn / 2 +
           vhat * d.L2 * rp.q₁y * invn

      dux = (gx - Xx) / 2
      duy = (gy - Xy) / 2

      dvx = (1 - uhat) * d.L2 * rp.q₁x * invn / 2 +
            uhat * ak * dgx / 2

      dvy = (1 - uhat) * d.L2 * rp.q₁y * invn / 2 +
            uhat * ak * dgy / 2

   elseif r == 6

      ak = d.tht2 + d.tht3
      th = muladd(ak, vhat, π / 2 - d.tht2)

      s, c = sincos(th)
      h = s^3 + c^3
      h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

      gx = d.R * h * c
      gy = d.R * h * s
      dgx = d.R * (h1 * c - h * s)
      dgy = d.R * (h1 * s + h * c)

      Xx = d.R * rp.e₃x +
           (1 - 2vhat) * d.L1 * rp.e₁x * invn / 2 +
           d.L2 * rp.q₁x * invn

      Xy = d.R * rp.e₃y +
           (1 - 2vhat) * d.L1 * rp.e₁y * invn / 2 +
           d.L2 * rp.q₁y * invn

      dux = (gx - Xx) / 2
      duy = (gy - Xy) / 2

      dvx = -(1 - uhat) * d.L1 * rp.e₁x * invn / 2 +
            uhat * ak * dgx / 2

      dvy = -(1 - uhat) * d.L1 * rp.e₁y * invn / 2 +
            uhat * ak * dgy / 2

   elseif r == 7

      ak = π / 4 - d.tht1 - d.tht3
      th = muladd(ak, vhat, π / 2 + d.tht3)

      s, c = sincos(th)
      h = s^3 + c^3
      h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

      gx = d.R * h * c
      gy = d.R * h * s
      dgx = d.R * (h1 * c - h * s)
      dgy = d.R * (h1 * s + h * c)

      Xx = d.R * rp.e₃x -
           d.L1 * rp.e₁x * invn / 2 +
           (1 - vhat) * d.L2 * rp.q₁x * invn

      Xy = d.R * rp.e₃y -
           d.L1 * rp.e₁y * invn / 2 +
           (1 - vhat) * d.L2 * rp.q₁y * invn

      dux = (gx - Xx) / 2
      duy = (gy - Xy) / 2

      dvx = -(1 - uhat) * d.L2 * rp.q₁x * invn / 2 +
            uhat * ak * dgx / 2

      dvy = -(1 - uhat) * d.L2 * rp.q₁y * invn / 2 +
            uhat * ak * dgy / 2

   elseif r == 8

      ak = d.tht1
      th = muladd(ak, vhat, 3π / 4 - d.tht1)

      s, c = sincos(th)
      h = s^3 + c^3
      h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

      gx = d.R * h * c
      gy = d.R * h * s
      dgx = d.R * (h1 * c - h * s)
      dgy = d.R * (h1 * s + h * c)

      Xx = d.R * vhat / 4 + d.R * (1 - vhat) * rp.e₃x
      Xy = d.R * vhat / 4 + d.R * (1 - vhat) * rp.e₃y

      dux = (gx - Xx) / 2
      duy = (gy - Xy) / 2

      dvx = (1 - uhat) * d.R * (1 / 4 - rp.e₃x) / 2 +
            uhat * ak * dgx / 2

      dvy = (1 - uhat) * d.R * (1 / 4 - rp.e₃y) / 2 +
            uhat * ak * dgy / 2

   elseif r == 9

      ak = d.tht1
      th = muladd(ak, vhat, 3π / 4)

      s, c = sincos(th)
      h = s^3 + c^3
      h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

      gx = d.R * h * c
      gy = d.R * h * s
      dgx = d.R * (h1 * c - h * s)
      dgy = d.R * (h1 * s + h * c)

      Xx = d.R * (1 - vhat) / 4 + d.R * vhat * rp.e₄x
      Xy = d.R * (1 - vhat) / 4 + d.R * vhat * rp.e₄y

      dux = (gx - Xx) / 2
      duy = (gy - Xy) / 2

      dvx = (1 - uhat) * d.R * (rp.e₄x - 1 / 4) / 2 +
            uhat * ak * dgx / 2

      dvy = (1 - uhat) * d.R * (rp.e₄y - 1 / 4) / 2 +
            uhat * ak * dgy / 2

   elseif r == 10

      ak = π / 4 - d.tht1 - d.tht3
      th = muladd(ak, vhat, 3π / 4 + d.tht1)

      s, c = sincos(th)
      h = s^3 + c^3
      h1 = 3.0 * s^2 * c - 3.0 * c^2 * s

      gx = d.R * h * c
      gy = d.R * h * s
      dgx = d.R * (h1 * c - h * s)
      dgy = d.R * (h1 * s + h * c)

      Xx = d.R * rp.e₄x -
           d.L1 * rp.e₂x * invn / 2 +
           vhat * d.L2 * rp.q₂x * invn

      Xy = d.R * rp.e₄y -
           d.L1 * rp.e₂y * invn / 2 +
           vhat * d.L2 * rp.q₂y * invn

      dux = (gx - Xx) / 2
      duy = (gy - Xy) / 2

      dvx = (1 - uhat) * d.L2 * rp.q₂x * invn / 2 +
            uhat * ak * dgx / 2

      dvy = (1 - uhat) * d.L2 * rp.q₂y * invn / 2 +
            uhat * ak * dgy / 2

   else

      throw(ArgumentError("jinvmap for bean expects region 1–10; got reg=$r"))

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
  mapinv(d::bean, u::Float64, v::Float64, k::Int) -> Tuple{Float64,Float64}

Inverse mapping: given a physical point `(u,v)` on patch `k`, return
the reference coordinates `[t; s]` in `[-1,1]^2` such that `mapxy(d, t, s, k) = (u,v)`.
"""
@inline function Xx(s::Float64, d::bean, r::Int)

   rp = d.RP
   invn = 1.0 / rp.nme

   if r == 1
      return d.A + d.R * rp.e₄x +
             (2s - 1) * d.L1 * rp.e₂x * invn / 2 +
             d.L2 * rp.q₂x * invn

   elseif r == 2
      return d.A + d.R * rp.e₄x +
             d.L1 * rp.e₂x * invn / 2 +
             (1 - s) * d.L2 * rp.q₂x * invn

   elseif r == 3
      return d.A + d.R * s / 4 + d.R * (1 - s) * rp.e₄x

   elseif r == 4
      return d.A + d.R * (1 - s) / 4 + d.R * s * rp.e₃x

   elseif r == 5
      return d.A + d.R * rp.e₃x +
             d.L1 * rp.e₁x * invn / 2 +
             s * d.L2 * rp.q₁x * invn

   elseif r == 6
      return d.A + d.R * rp.e₃x +
             (1 - 2s) * d.L1 * rp.e₁x * invn / 2 +
             d.L2 * rp.q₁x * invn

   elseif r == 7
      return d.A + d.R * rp.e₃x -
             d.L1 * rp.e₁x * invn / 2 +
             (1 - s) * d.L2 * rp.q₁x * invn

   elseif r == 8
      return d.A + d.R * s / 4 + d.R * (1 - s) * rp.e₃x

   elseif r == 9
      return d.A + d.R * (1 - s) / 4 + d.R * s * rp.e₄x

   elseif r == 10
      return d.A + d.R * rp.e₄x -
             d.L1 * rp.e₂x * invn / 2 +
             s * d.L2 * rp.q₂x * invn

   else
      throw(ArgumentError("Xx for bean expects region 1–10; got r=$r"))
   end

end

@inline function Xy(s::Float64, d::bean, r::Int)

   rp = d.RP
   invn = 1.0 / rp.nme

   if r == 1
      return d.B + d.R * rp.e₄y +
             (2s - 1) * d.L1 * rp.e₂y * invn / 2 +
             d.L2 * rp.q₂y * invn

   elseif r == 2
      return d.B + d.R * rp.e₄y +
             d.L1 * rp.e₂y * invn / 2 +
             (1 - s) * d.L2 * rp.q₂y * invn

   elseif r == 3
      return d.B + d.R * s / 4 + d.R * (1 - s) * rp.e₄y

   elseif r == 4
      return d.B + d.R * (1 - s) / 4 + d.R * s * rp.e₃y

   elseif r == 5
      return d.B + d.R * rp.e₃y +
             d.L1 * rp.e₁y * invn / 2 +
             s * d.L2 * rp.q₁y * invn

   elseif r == 6
      return d.B + d.R * rp.e₃y +
             (1 - 2s) * d.L1 * rp.e₁y * invn / 2 +
             d.L2 * rp.q₁y * invn

   elseif r == 7
      return d.B + d.R * rp.e₃y -
             d.L1 * rp.e₁y * invn / 2 +
             (1 - s) * d.L2 * rp.q₁y * invn

   elseif r == 8
      return d.B + d.R * s / 4 + d.R * (1 - s) * rp.e₃y

   elseif r == 9
      return d.B + d.R * (1 - s) / 4 + d.R * s * rp.e₄y

   elseif r == 10
      return d.B + d.R * rp.e₄y -
             d.L1 * rp.e₂y * invn / 2 +
             s * d.L2 * rp.q₂y * invn

   else
      throw(ArgumentError("Xy for bean expects region 1–10; got r=$r"))
   end

end

@inline function Yx(s::Float64, d::bean, r::Int)

   θ = if r == 1
      s * (d.tht2 + d.tht3) - d.tht3
   elseif r == 2
      s * (π / 4 - d.tht1 - d.tht2) + d.tht2
   elseif r == 3
      s * d.tht1 + π / 4 - d.tht1
   elseif r == 4
      s * d.tht1 + π / 4
   elseif r == 5
      s * (π / 4 - d.tht1 - d.tht2) + π / 4 + d.tht1
   elseif r == 6
      s * (d.tht2 + d.tht3) + π / 2 - d.tht2
   elseif r == 7
      s * (π / 4 - d.tht1 - d.tht3) + π / 2 + d.tht3
   elseif r == 8
      s * d.tht1 + 3π / 4 - d.tht1
   elseif r == 9
      s * d.tht1 + 3π / 4
   elseif r == 10
      s * (π / 4 - d.tht1 - d.tht3) + 3π / 4 + d.tht1
   else
      throw(ArgumentError("Yx for bean expects region 1–10; got r=$r"))
   end

   st, ct = sincos(θ)
   h = st^3 + ct^3

   return d.A + d.R * h * ct

end

@inline function Yy(s::Float64, d::bean, r::Int)

   θ = if r == 1
      s * (d.tht2 + d.tht3) - d.tht3
   elseif r == 2
      s * (π / 4 - d.tht1 - d.tht2) + d.tht2
   elseif r == 3
      s * d.tht1 + π / 4 - d.tht1
   elseif r == 4
      s * d.tht1 + π / 4
   elseif r == 5
      s * (π / 4 - d.tht1 - d.tht2) + π / 4 + d.tht1
   elseif r == 6
      s * (d.tht2 + d.tht3) + π / 2 - d.tht2
   elseif r == 7
      s * (π / 4 - d.tht1 - d.tht3) + π / 2 + d.tht3
   elseif r == 8
      s * d.tht1 + 3π / 4 - d.tht1
   elseif r == 9
      s * d.tht1 + 3π / 4
   elseif r == 10
      s * (π / 4 - d.tht1 - d.tht3) + 3π / 4 + d.tht1
   else
      throw(ArgumentError("Yy for bean expects region 1–10; got r=$r"))
   end

   st, ct = sincos(θ)
   h = st^3 + ct^3

   return d.B + d.R * h * st

end

function fill_FTable!(tbl::FTable, d::bean, r::Int)

   vmin, vmax, N = tbl.vmin, tbl.vmax, tbl.N
   hstep = (vmax - vmin) / (N - 1)

   P1, P2, P3 = tbl.P1, tbl.P2, tbl.P3

   rp = d.RP
   invn = 1.0 / rp.nme

   if r == 1

      ak = d.tht2 + d.tht3
      bk = -d.tht3

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * hstep

         Xxv = d.A + d.R * rp.e₄x +
               (2v̂ - 1) * d.L1 * rp.e₂x * invn / 2 +
               d.L2 * rp.q₂x * invn

         Xyv = d.B + d.R * rp.e₄y +
               (2v̂ - 1) * d.L1 * rp.e₂y * invn / 2 +
               d.L2 * rp.q₂y * invn

         θ = muladd(ak, v̂, bk)
         st, ct = sincos(θ)
         h = st^3 + ct^3

         Yxv = d.A + d.R * h * ct
         Yyv = d.B + d.R * h * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 2

      ak = π / 4 - d.tht1 - d.tht2
      bk = d.tht2

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * hstep

         Xxv = d.A + d.R * rp.e₄x +
               d.L1 * rp.e₂x * invn / 2 +
               (1 - v̂) * d.L2 * rp.q₂x * invn

         Xyv = d.B + d.R * rp.e₄y +
               d.L1 * rp.e₂y * invn / 2 +
               (1 - v̂) * d.L2 * rp.q₂y * invn

         θ = muladd(ak, v̂, bk)
         st, ct = sincos(θ)
         h = st^3 + ct^3

         Yxv = d.A + d.R * h * ct
         Yyv = d.B + d.R * h * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 3

      ak = d.tht1
      bk = π / 4 - d.tht1

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * hstep

         Xxv = d.A + d.R * v̂ / 4 + d.R * (1 - v̂) * rp.e₄x
         Xyv = d.B + d.R * v̂ / 4 + d.R * (1 - v̂) * rp.e₄y

         θ = muladd(ak, v̂, bk)
         st, ct = sincos(θ)
         h = st^3 + ct^3

         Yxv = d.A + d.R * h * ct
         Yyv = d.B + d.R * h * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 4

      ak = d.tht1
      bk = π / 4

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * hstep

         Xxv = d.A + d.R * (1 - v̂) / 4 + d.R * v̂ * rp.e₃x
         Xyv = d.B + d.R * (1 - v̂) / 4 + d.R * v̂ * rp.e₃y

         θ = muladd(ak, v̂, bk)
         st, ct = sincos(θ)
         h = st^3 + ct^3

         Yxv = d.A + d.R * h * ct
         Yyv = d.B + d.R * h * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 5

      ak = π / 4 - d.tht1 - d.tht2
      bk = π / 4 + d.tht1

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * hstep

         Xxv = d.A + d.R * rp.e₃x +
               d.L1 * rp.e₁x * invn / 2 +
               v̂ * d.L2 * rp.q₁x * invn

         Xyv = d.B + d.R * rp.e₃y +
               d.L1 * rp.e₁y * invn / 2 +
               v̂ * d.L2 * rp.q₁y * invn

         θ = muladd(ak, v̂, bk)
         st, ct = sincos(θ)
         h = st^3 + ct^3

         Yxv = d.A + d.R * h * ct
         Yyv = d.B + d.R * h * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 6

      ak = d.tht2 + d.tht3
      bk = π / 2 - d.tht2

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * hstep

         Xxv = d.A + d.R * rp.e₃x +
               (1 - 2v̂) * d.L1 * rp.e₁x * invn / 2 +
               d.L2 * rp.q₁x * invn

         Xyv = d.B + d.R * rp.e₃y +
               (1 - 2v̂) * d.L1 * rp.e₁y * invn / 2 +
               d.L2 * rp.q₁y * invn

         θ = muladd(ak, v̂, bk)
         st, ct = sincos(θ)
         h = st^3 + ct^3

         Yxv = d.A + d.R * h * ct
         Yyv = d.B + d.R * h * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 7

      ak = π / 4 - d.tht1 - d.tht3
      bk = π / 2 + d.tht3

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * hstep

         Xxv = d.A + d.R * rp.e₃x -
               d.L1 * rp.e₁x * invn / 2 +
               (1 - v̂) * d.L2 * rp.q₁x * invn

         Xyv = d.B + d.R * rp.e₃y -
               d.L1 * rp.e₁y * invn / 2 +
               (1 - v̂) * d.L2 * rp.q₁y * invn

         θ = muladd(ak, v̂, bk)
         st, ct = sincos(θ)
         h = st^3 + ct^3

         Yxv = d.A + d.R * h * ct
         Yyv = d.B + d.R * h * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 8

      ak = d.tht1
      bk = 3π / 4 - d.tht1

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * hstep

         Xxv = d.A + d.R * v̂ / 4 + d.R * (1 - v̂) * rp.e₃x
         Xyv = d.B + d.R * v̂ / 4 + d.R * (1 - v̂) * rp.e₃y

         θ = muladd(ak, v̂, bk)
         st, ct = sincos(θ)
         h = st^3 + ct^3

         Yxv = d.A + d.R * h * ct
         Yyv = d.B + d.R * h * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 9

      ak = d.tht1
      bk = 3π / 4

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * hstep

         Xxv = d.A + d.R * (1 - v̂) / 4 + d.R * v̂ * rp.e₄x
         Xyv = d.B + d.R * (1 - v̂) / 4 + d.R * v̂ * rp.e₄y

         θ = muladd(ak, v̂, bk)
         st, ct = sincos(θ)
         h = st^3 + ct^3

         Yxv = d.A + d.R * h * ct
         Yyv = d.B + d.R * h * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   elseif r == 10

      ak = π / 4 - d.tht1 - d.tht3
      bk = 3π / 4 + d.tht1

      @inbounds for i in 1:N
         v̂ = vmin + (i - 1) * hstep

         Xxv = d.A + d.R * rp.e₄x -
               d.L1 * rp.e₂x * invn / 2 +
               v̂ * d.L2 * rp.q₂x * invn

         Xyv = d.B + d.R * rp.e₄y -
               d.L1 * rp.e₂y * invn / 2 +
               v̂ * d.L2 * rp.q₂y * invn

         θ = muladd(ak, v̂, bk)
         st, ct = sincos(θ)
         h = st^3 + ct^3

         Yxv = d.A + d.R * h * ct
         Yyv = d.B + d.R * h * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv
      end

   else
      throw(ArgumentError("fill_FTable! for bean expects region 1–10; got r=$r"))
   end

   tbl.reg = r
   return tbl

end

@inline function f1I(t::Float64, s::Float64,
  d::bean, u::Float64, v::Float64, r::Int)
  ŝ = (s + 1) / 2
  t̂ = (t + 1) / 2

  Xxv, Yxv = Xx(ŝ, d, r),  Yx(ŝ, d, r)

  return (1 - t̂) * Xxv + t̂ * Yxv - u
end

@inline function f2I(t::Float64, s::Float64,
  d::bean, u::Float64, v::Float64, r::Int)
    ŝ = (s + 1) / 2
    t̂ = (t + 1) / 2

    Xyv, Yyv = Xy(ŝ, d, r), Yy(ŝ, d, r)

    return (1 - t̂) * Xyv + t̂ * Yyv - v
end

@inline function f_cont(v̂::Float64, d::bean,
   u::Float64, v::Float64, r::Int)

   rp = d.RP
   invn = 1.0 / rp.nme

   if r == 1

      Xxv = d.A + d.R * rp.e₄x +
            (2v̂ - 1) * d.L1 * rp.e₂x * invn / 2 +
            d.L2 * rp.q₂x * invn

      Xyv = d.B + d.R * rp.e₄y +
            (2v̂ - 1) * d.L1 * rp.e₂y * invn / 2 +
            d.L2 * rp.q₂y * invn

      θ = v̂ * (d.tht2 + d.tht3) - d.tht3

   elseif r == 2

      Xxv = d.A + d.R * rp.e₄x +
            d.L1 * rp.e₂x * invn / 2 +
            (1 - v̂) * d.L2 * rp.q₂x * invn

      Xyv = d.B + d.R * rp.e₄y +
            d.L1 * rp.e₂y * invn / 2 +
            (1 - v̂) * d.L2 * rp.q₂y * invn

      θ = v̂ * (π / 4 - d.tht1 - d.tht2) + d.tht2

   elseif r == 3

      Xxv = d.A + d.R * v̂ / 4 + d.R * (1 - v̂) * rp.e₄x
      Xyv = d.B + d.R * v̂ / 4 + d.R * (1 - v̂) * rp.e₄y

      θ = v̂ * d.tht1 + π / 4 - d.tht1

   elseif r == 4

      Xxv = d.A + d.R * (1 - v̂) / 4 + d.R * v̂ * rp.e₃x
      Xyv = d.B + d.R * (1 - v̂) / 4 + d.R * v̂ * rp.e₃y

      θ = v̂ * d.tht1 + π / 4

   elseif r == 5

      Xxv = d.A + d.R * rp.e₃x +
            d.L1 * rp.e₁x * invn / 2 +
            v̂ * d.L2 * rp.q₁x * invn

      Xyv = d.B + d.R * rp.e₃y +
            d.L1 * rp.e₁y * invn / 2 +
            v̂ * d.L2 * rp.q₁y * invn

      θ = v̂ * (π / 4 - d.tht1 - d.tht2) + π / 4 + d.tht1

   elseif r == 6

      Xxv = d.A + d.R * rp.e₃x +
            (1 - 2v̂) * d.L1 * rp.e₁x * invn / 2 +
            d.L2 * rp.q₁x * invn

      Xyv = d.B + d.R * rp.e₃y +
            (1 - 2v̂) * d.L1 * rp.e₁y * invn / 2 +
            d.L2 * rp.q₁y * invn

      θ = v̂ * (d.tht2 + d.tht3) + π / 2 - d.tht2

   elseif r == 7

      Xxv = d.A + d.R * rp.e₃x -
            d.L1 * rp.e₁x * invn / 2 +
            (1 - v̂) * d.L2 * rp.q₁x * invn

      Xyv = d.B + d.R * rp.e₃y -
            d.L1 * rp.e₁y * invn / 2 +
            (1 - v̂) * d.L2 * rp.q₁y * invn

      θ = v̂ * (π / 4 - d.tht1 - d.tht3) + π / 2 + d.tht3

   elseif r == 8

      Xxv = d.A + d.R * v̂ / 4 + d.R * (1 - v̂) * rp.e₃x
      Xyv = d.B + d.R * v̂ / 4 + d.R * (1 - v̂) * rp.e₃y

      θ = v̂ * d.tht1 + 3π / 4 - d.tht1

   elseif r == 9

      Xxv = d.A + d.R * (1 - v̂) / 4 + d.R * v̂ * rp.e₄x
      Xyv = d.B + d.R * (1 - v̂) / 4 + d.R * v̂ * rp.e₄y

      θ = v̂ * d.tht1 + 3π / 4

   elseif r == 10

      Xxv = d.A + d.R * rp.e₄x -
            d.L1 * rp.e₂x * invn / 2 +
            v̂ * d.L2 * rp.q₂x * invn

      Xyv = d.B + d.R * rp.e₄y -
            d.L1 * rp.e₂y * invn / 2 +
            v̂ * d.L2 * rp.q₂y * invn

      θ = v̂ * (π / 4 - d.tht1 - d.tht3) + 3π / 4 + d.tht1

   else
      throw(ArgumentError("f_cont for bean expects region 1–10; got r=$r"))
   end

   st, ct = sincos(θ)
   h = st^3 + ct^3

   Yxv = d.A + d.R * h * ct
   Yyv = d.B + d.R * h * st

   dx = Yxv - Xxv
   dy = Yyv - Xyv

   term1 = muladd(u, dy, -v * dx)
   term2 = muladd(Xyv, Yxv, -Xxv * Yyv)

   return term1 + term2

end

@inline function JinvI(t::Float64, s::Float64,
  d::bean, u::Float64, v::Float64, r::Int)
  return jinvmap(d, t, s, r)  # returns J11,J12,J21,J22
end


function mapinv(tbl::FTable, d::bean, u::Float64,
   v::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]
   r = p.reg

   if !(1 <= r <= 10)
      throw(ArgumentError("mapinv(tbl, d::bean, ...) expects region 1–10; got reg=$r"))
   end

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

# Separate inverse for rectangular regions 11 and 12.
function mapinv(d::bean, u::Float64,
   v::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]
   r = p.reg
   rp = d.RP

   if r == 11

      wx = u - d.A - d.R * rp.e₄x
      wy = v - d.B - d.R * rp.e₄y

      xi1 = (wx * rp.q₂x + wy * rp.q₂y) / (d.L2 * rp.nme)
      xi2 = 0.5 + (wx * rp.e₂x + wy * rp.e₂y) / (d.L1 * rp.nme)

      if d.L2 >= d.L1
         Z1 = xi_inv(xi1, p.tk0, p.tk1)
         Z2 = xi_inv(xi2, p.ck0, p.ck1)
      else
         Z1 = xi_inv(xi1, p.ck0, p.ck1)
         Z2 = xi_inv(xi2, p.tk0, p.tk1)
      end

      return Z1, Z2

   elseif r == 12

      wx = u - d.A - d.R * rp.e₃x
      wy = v - d.B - d.R * rp.e₃y

      xi1 = (wx * rp.q₁x + wy * rp.q₁y) / (d.L2 * rp.nme)
      xi2 = 0.5 + (wx * rp.e₁x + wy * rp.e₁y) / (d.L1 * rp.nme)

      if d.L2 >= d.L1
         Z1 = xi_inv(xi1, p.tk0, p.tk1)
         Z2 = xi_inv(xi2, p.ck0, p.ck1)
      else
         Z1 = xi_inv(xi1, p.ck0, p.ck1)
         Z2 = xi_inv(xi2, p.tk0, p.tk1)
      end

      return Z1, Z2

   end

   throw(ArgumentError("mapinv(d::bean, u, v, k) is only for rectangular regions 11 and 12; got reg=$r"))

end

"""
    ptconv(d::bean, t1::Float64, t2::Float64, idx::Int, ptdest::String)

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
function ptconv(d::bean, t1::Float64, t2::Float64, idx::Int, ptdest::String)

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
    mapinv2(d::bean, t1::Float64, t2::Float64, k2::Int, k::Int) -> Tuple{Float64,Float64}

Given a point `(u,v) = τₖ₂(t1,t2)` on patch `k2`, return its reference
coordinates on patch `k` in `[-1,1]^2`.

This differs from `mapinv` because we already know `(u,v)` comes from
`(t1,t2)` on `k2`.

This should only be called when patches `k2` and `k` are in the same region.
"""
function mapinv2(d::bean, t1::Float64, t2::Float64,
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
function dfunc(d::bean, k::Int, t::Float64, s::Float64)::Float64

   p = d.pths[k]

   if p.reg == 11 || p.reg == 12
      return 1.0
   end

   hc = p.ck1 - p.ck0
   val = (1.0 - p.ck1) + hc * t / 2
   exp = s ≥ 0.5 ? (s - 1) : s

   return val^exp

end

function dfunc!(out::StridedArray{Float64}, d::bean, k::Int,
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

function dfunc!(out::StridedArray{Float64}, d::bean, k::Int, t::StridedArray{Float64})

   p = d.pths[k]

   if p.reg == 11 || p.reg == 12

      fill!(out, 1.0)

   else
      αc = (p.ck1 - p.ck0) / 2
      βc = 1.0 - p.ck1

      @inbounds for i in eachindex(t)
         out[i] = muladd(αc, t[i], βc)
      end

   end

   return nothing

end

"""
   gamderhigher(d::bean, th::Float64)

Return centered bean boundary values and derivatives through sixth order.

Returns gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy, d4gx, d4gy, d5gx, d5gy, d6gx, d6gy

where gx(th) = R*h(th)*cos(th), gy(th) = R*h(th)*sin(th),
and h(th) = sin(th)^3 + cos(th)^3.
"""
function gamderhigher(d::bean, th::Float64)

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
  diff_map!(out::Matrix{Float64}, Zx::Matrix{Float64}, Zy::Matrix{Float64},
            DJ::StridedArray{Float64}, d::bean,
            u::Float64, v::Float64,
            u2::Matrix{Float64}, v2::Matrix{Float64},
            du::AbstractVector, dv::AbstractVector, k::Int;
            tol = 1e-3)

Fill `out` with ‖τ(u,v) - τ(u₂,v₂)‖ on patch `k`.

No allocations. Uses `mapxy_Dmap!` and Taylor correction near `(u,v)`.
"""
function diff_map!(out::Matrix{Float64},
   Zx::Matrix{Float64}, Zy::Matrix{Float64}, DJ::StridedArray{Float64},
   d::bean, u::Float64, v::Float64,
   u2::Matrix{Float64}, v2::Matrix{Float64},
   du::AbstractVector, dv::AbstractVector, k::Int;
   tol::Float64=1e-3)

   nd_u = size(out, 1)
   nd_v = size(out, 2)

   p = d.pths[k]
   rp = d.RP

   αc = (p.ck1 - p.ck0) / 2
   αt = (p.tk1 - p.tk0) / 2

   invn = 1.0 / rp.nme

   mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)

   if p.reg == 11 || p.reg == 12

      if d.L2 >= d.L1
         αu = αt
         αv = αc
      else
         αu = αc
         αv = αt
      end

      if p.reg == 11

         qx = rp.q₂x
         qy = rp.q₂y
         ex = rp.e₂x
         ey = rp.e₂y

      else # p.reg == 12

         qx = rp.q₁x
         qy = rp.q₁y
         ex = rp.e₁x
         ey = rp.e₁y

      end

      @inbounds for j in 1:nd_v
         dvj = dv[j]

         @inbounds for i in 1:nd_u
            dui = du[i]

            Dx = d.L2 * αu * dui * qx * invn +
                 d.L1 * αv * dvj * ex * invn

            Dy = d.L2 * αu * dui * qy * invn +
                 d.L1 * αv * dvj * ey * invn

            out[i, j] = hypot(Dx, Dy)
         end
      end

      return nothing

   end

   uhat = muladd(αc, u, p.ck0 + αc)
   vhat = muladd(αt, v, p.tk0 + αt)

   if p.reg == 1

      ak = d.tht2 + d.tht3
      bk = -d.tht3

      Xx = d.R * rp.e₄x +
           (2.0 * vhat - 1.0) * d.L1 * rp.e₂x * invn / 2.0 +
           d.L2 * rp.q₂x * invn

      Xy = d.R * rp.e₄y +
           (2.0 * vhat - 1.0) * d.L1 * rp.e₂y * invn / 2.0 +
           d.L2 * rp.q₂y * invn

      dXx = d.L1 * rp.e₂x * invn
      dXy = d.L1 * rp.e₂y * invn

   elseif p.reg == 2

      ak = π / 4 - d.tht1 - d.tht2
      bk = d.tht2

      Xx = d.R * rp.e₄x +
           d.L1 * rp.e₂x * invn / 2.0 +
           (1.0 - vhat) * d.L2 * rp.q₂x * invn

      Xy = d.R * rp.e₄y +
           d.L1 * rp.e₂y * invn / 2.0 +
           (1.0 - vhat) * d.L2 * rp.q₂y * invn

      dXx = -d.L2 * rp.q₂x * invn
      dXy = -d.L2 * rp.q₂y * invn

   elseif p.reg == 3

      ak = d.tht1
      bk = π / 4 - d.tht1

      Xx = d.R * vhat / 4.0 + d.R * (1.0 - vhat) * rp.e₄x
      Xy = d.R * vhat / 4.0 + d.R * (1.0 - vhat) * rp.e₄y

      dXx = d.R * (1.0 / 4.0 - rp.e₄x)
      dXy = d.R * (1.0 / 4.0 - rp.e₄y)

   elseif p.reg == 4

      ak = d.tht1
      bk = π / 4

      Xx = d.R * (1.0 - vhat) / 4.0 + d.R * vhat * rp.e₃x
      Xy = d.R * (1.0 - vhat) / 4.0 + d.R * vhat * rp.e₃y

      dXx = d.R * (rp.e₃x - 1.0 / 4.0)
      dXy = d.R * (rp.e₃y - 1.0 / 4.0)

   elseif p.reg == 5

      ak = π / 4 - d.tht1 - d.tht2
      bk = π / 4 + d.tht1

      Xx = d.R * rp.e₃x +
           d.L1 * rp.e₁x * invn / 2.0 +
           vhat * d.L2 * rp.q₁x * invn

      Xy = d.R * rp.e₃y +
           d.L1 * rp.e₁y * invn / 2.0 +
           vhat * d.L2 * rp.q₁y * invn

      dXx = d.L2 * rp.q₁x * invn
      dXy = d.L2 * rp.q₁y * invn

   elseif p.reg == 6

      ak = d.tht2 + d.tht3
      bk = π / 2 - d.tht2

      Xx = d.R * rp.e₃x +
           (1.0 - 2.0 * vhat) * d.L1 * rp.e₁x * invn / 2.0 +
           d.L2 * rp.q₁x * invn

      Xy = d.R * rp.e₃y +
           (1.0 - 2.0 * vhat) * d.L1 * rp.e₁y * invn / 2.0 +
           d.L2 * rp.q₁y * invn

      dXx = -d.L1 * rp.e₁x * invn
      dXy = -d.L1 * rp.e₁y * invn

   elseif p.reg == 7

      ak = π / 4 - d.tht1 - d.tht3
      bk = π / 2 + d.tht3

      Xx = d.R * rp.e₃x -
           d.L1 * rp.e₁x * invn / 2.0 +
           (1.0 - vhat) * d.L2 * rp.q₁x * invn

      Xy = d.R * rp.e₃y -
           d.L1 * rp.e₁y * invn / 2.0 +
           (1.0 - vhat) * d.L2 * rp.q₁y * invn

      dXx = -d.L2 * rp.q₁x * invn
      dXy = -d.L2 * rp.q₁y * invn

   elseif p.reg == 8

      ak = d.tht1
      bk = 3π / 4 - d.tht1

      Xx = d.R * vhat / 4.0 + d.R * (1.0 - vhat) * rp.e₃x
      Xy = d.R * vhat / 4.0 + d.R * (1.0 - vhat) * rp.e₃y

      dXx = d.R * (1.0 / 4.0 - rp.e₃x)
      dXy = d.R * (1.0 / 4.0 - rp.e₃y)

   elseif p.reg == 9

      ak = d.tht1
      bk = 3π / 4

      Xx = d.R * (1.0 - vhat) / 4.0 + d.R * vhat * rp.e₄x
      Xy = d.R * (1.0 - vhat) / 4.0 + d.R * vhat * rp.e₄y

      dXx = d.R * (rp.e₄x - 1.0 / 4.0)
      dXy = d.R * (rp.e₄y - 1.0 / 4.0)

   elseif p.reg == 10

      ak = π / 4 - d.tht1 - d.tht3
      bk = 3π / 4 + d.tht1

      Xx = d.R * rp.e₄x -
           d.L1 * rp.e₂x * invn / 2.0 +
           vhat * d.L2 * rp.q₂x * invn

      Xy = d.R * rp.e₄y -
           d.L1 * rp.e₂y * invn / 2.0 +
           vhat * d.L2 * rp.q₂y * invn

      dXx = d.L2 * rp.q₂x * invn
      dXy = d.L2 * rp.q₂y * invn

   else

      throw(ArgumentError("diff_map! for bean expects region 1–12; got reg=$(p.reg)"))

   end

   th = muladd(ak, vhat, bk)

   gx, gy,
   dgx, dgy,
   d2gx, d2gy,
   d3gx, d3gy,
   d4gx, d4gy,
   d5gx, d5gy,
   d6gx, d6gy = gamderhigher(d, th)

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
  diff_rmap!(out::Matrix{Float64}, Zx::Matrix{Float64}, Zy::Matrix{Float64},
             DJ::StridedArray{Float64}, d::bean,
             u::Float64, v::Float64,
             u2::Matrix{Float64}, v2::Matrix{Float64},
             r::Matrix{Float64}, du::AbstractVector,
             dv::AbstractVector, k::Int;
             tol = 1e-3)

Fill `out` with ‖τ(u,v) - τ(u₂,v₂)‖ / r on patch `k`.

No allocations. Uses `mapxy_Dmap!` and Taylor correction near `(u,v)`.
"""
function diff_rmap!(out::Matrix{Float64},
   Zx::Matrix{Float64}, Zy::Matrix{Float64}, DJ::StridedArray{Float64},
   d::bean, u::Float64, v::Float64,
   u2::Matrix{Float64}, v2::Matrix{Float64},
   r::Matrix{Float64}, du::AbstractVector,
   dv::AbstractVector, k::Int;
   tol::Float64 = 1e-3)

   nd_u = size(out, 1)
   nd_v = size(out, 2)

   p = d.pths[k]
   rp = d.RP

   αc = (p.ck1 - p.ck0) / 2
   αt = (p.tk1 - p.tk0) / 2

   invn = 1.0 / rp.nme

   # Rectangular regions are affine, so the scaled difference is exact.
   if p.reg == 11 || p.reg == 12

      fill!(DJ, d.L1 * d.L2 * αc * αt)

      if d.L2 >= d.L1
         cV = d.L1 * αc * invn
         cU = d.L2 * αt * invn
      else
         cV = d.L1 * αt * invn
         cU = d.L2 * αc * invn
      end

      if p.reg == 11
         qx = rp.q₂x
         qy = rp.q₂y
         ex = rp.e₂x
         ey = rp.e₂y
      else
         qx = rp.q₁x
         qy = rp.q₁y
         ex = rp.e₁x
         ey = rp.e₁y
      end

      @inbounds for i in 1:nd_u
         dui = du[i]
         dvi = dv[i]

         @inbounds for j in 1:nd_v

            Dx = (cV * dvi) * ex + (cU * dui) * qx
            Dy = (cV * dvi) * ey + (cU * dui) * qy

            out[i, j] = hypot(Dx, Dy)
         end
      end

      return nothing

   end

   mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)

   uhat = muladd(αc, u, p.ck0 + αc)
   vhat = muladd(αt, v, p.tk0 + αt)

   if p.reg == 1

      ak = d.tht2 + d.tht3
      bk = -d.tht3

      Xx = d.R * rp.e₄x +
           (2.0 * vhat - 1.0) * d.L1 * rp.e₂x * invn / 2.0 +
           d.L2 * rp.q₂x * invn

      Xy = d.R * rp.e₄y +
           (2.0 * vhat - 1.0) * d.L1 * rp.e₂y * invn / 2.0 +
           d.L2 * rp.q₂y * invn

      dXx = d.L1 * rp.e₂x * invn
      dXy = d.L1 * rp.e₂y * invn

   elseif p.reg == 2

      ak = π / 4 - d.tht1 - d.tht2
      bk = d.tht2

      Xx = d.R * rp.e₄x +
           d.L1 * rp.e₂x * invn / 2.0 +
           (1.0 - vhat) * d.L2 * rp.q₂x * invn

      Xy = d.R * rp.e₄y +
           d.L1 * rp.e₂y * invn / 2.0 +
           (1.0 - vhat) * d.L2 * rp.q₂y * invn

      dXx = -d.L2 * rp.q₂x * invn
      dXy = -d.L2 * rp.q₂y * invn

   elseif p.reg == 3

      ak = d.tht1
      bk = π / 4 - d.tht1

      Xx = d.R * vhat / 4.0 + d.R * (1.0 - vhat) * rp.e₄x
      Xy = d.R * vhat / 4.0 + d.R * (1.0 - vhat) * rp.e₄y

      dXx = d.R * (1.0 / 4.0 - rp.e₄x)
      dXy = d.R * (1.0 / 4.0 - rp.e₄y)

   elseif p.reg == 4

      ak = d.tht1
      bk = π / 4

      Xx = d.R * (1.0 - vhat) / 4.0 + d.R * vhat * rp.e₃x
      Xy = d.R * (1.0 - vhat) / 4.0 + d.R * vhat * rp.e₃y

      dXx = d.R * (rp.e₃x - 1.0 / 4.0)
      dXy = d.R * (rp.e₃y - 1.0 / 4.0)

   elseif p.reg == 5

      ak = π / 4 - d.tht1 - d.tht2
      bk = π / 4 + d.tht1

      Xx = d.R * rp.e₃x +
           d.L1 * rp.e₁x * invn / 2.0 +
           vhat * d.L2 * rp.q₁x * invn

      Xy = d.R * rp.e₃y +
           d.L1 * rp.e₁y * invn / 2.0 +
           vhat * d.L2 * rp.q₁y * invn

      dXx = d.L2 * rp.q₁x * invn
      dXy = d.L2 * rp.q₁y * invn

   elseif p.reg == 6

      ak = d.tht2 + d.tht3
      bk = π / 2 - d.tht2

      Xx = d.R * rp.e₃x +
           (1.0 - 2.0 * vhat) * d.L1 * rp.e₁x * invn / 2.0 +
           d.L2 * rp.q₁x * invn

      Xy = d.R * rp.e₃y +
           (1.0 - 2.0 * vhat) * d.L1 * rp.e₁y * invn / 2.0 +
           d.L2 * rp.q₁y * invn

      dXx = -d.L1 * rp.e₁x * invn
      dXy = -d.L1 * rp.e₁y * invn

   elseif p.reg == 7

      ak = π / 4 - d.tht1 - d.tht3
      bk = π / 2 + d.tht3

      Xx = d.R * rp.e₃x -
           d.L1 * rp.e₁x * invn / 2.0 +
           (1.0 - vhat) * d.L2 * rp.q₁x * invn

      Xy = d.R * rp.e₃y -
           d.L1 * rp.e₁y * invn / 2.0 +
           (1.0 - vhat) * d.L2 * rp.q₁y * invn

      dXx = -d.L2 * rp.q₁x * invn
      dXy = -d.L2 * rp.q₁y * invn

   elseif p.reg == 8

      ak = d.tht1
      bk = 3π / 4 - d.tht1

      Xx = d.R * vhat / 4.0 + d.R * (1.0 - vhat) * rp.e₃x
      Xy = d.R * vhat / 4.0 + d.R * (1.0 - vhat) * rp.e₃y

      dXx = d.R * (1.0 / 4.0 - rp.e₃x)
      dXy = d.R * (1.0 / 4.0 - rp.e₃y)

   elseif p.reg == 9

      ak = d.tht1
      bk = 3π / 4

      Xx = d.R * (1.0 - vhat) / 4.0 + d.R * vhat * rp.e₄x
      Xy = d.R * (1.0 - vhat) / 4.0 + d.R * vhat * rp.e₄y

      dXx = d.R * (rp.e₄x - 1.0 / 4.0)
      dXy = d.R * (rp.e₄y - 1.0 / 4.0)

   elseif p.reg == 10

      ak = π / 4 - d.tht1 - d.tht3
      bk = 3π / 4 + d.tht1

      Xx = d.R * rp.e₄x -
           d.L1 * rp.e₂x * invn / 2.0 +
           vhat * d.L2 * rp.q₂x * invn

      Xy = d.R * rp.e₄y -
           d.L1 * rp.e₂y * invn / 2.0 +
           vhat * d.L2 * rp.q₂y * invn

      dXx = d.L2 * rp.q₂x * invn
      dXy = d.L2 * rp.q₂y * invn

   else

      throw(ArgumentError("diff_rmap! for bean expects region 1–12; got reg=$(p.reg)"))

   end

   th = muladd(ak, vhat, bk)

   gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy,
   d4gx, d4gy, d5gx, d5gy, d6gx, d6gy = gamderhigher(d, th)

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

   @inbounds for i in 1:nd_u

      dui = du[i]
      dvi = dv[i]

      @inbounds for j in 1:nd_v

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
   Dwall(d::bean, u::Float64, v::Float64, k::Int) -> dvx, dvy

Derivative of the wall curve w.r.t. `v` at the singular side `u = ±1`
for the `k`-th non-rectangular patch.

Returns a Tuple `(dvx, dvy)`.
"""
function Dwall(d::bean, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64}

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

   rp = d.RP
   invn = 1.0 / rp.nme

   if reg == 1

      ak = d.tht2 + d.tht3
      th = muladd(ak, xi2, -d.tht3)

      _, _, dgx, dgy = gamder(d, th)

      dvx = (1.0 - xi1) * d.L1 * rp.e₂x * invn + xi1 * ak * dgx
      dvy = (1.0 - xi1) * d.L1 * rp.e₂y * invn + xi1 * ak * dgy

   elseif reg == 2

      ak = π / 4 - d.tht1 - d.tht2
      th = muladd(ak, xi2, d.tht2)

      _, _, dgx, dgy = gamder(d, th)

      dvx = -(1.0 - xi1) * d.L2 * rp.q₂x * invn + xi1 * ak * dgx
      dvy = -(1.0 - xi1) * d.L2 * rp.q₂y * invn + xi1 * ak * dgy

   elseif reg == 3

      ak = d.tht1
      th = muladd(ak, xi2, π / 4 - d.tht1)

      _, _, dgx, dgy = gamder(d, th)

      dvx = (1.0 - xi1) * d.R * (1.0 / 4 - rp.e₄x) + xi1 * ak * dgx
      dvy = (1.0 - xi1) * d.R * (1.0 / 4 - rp.e₄y) + xi1 * ak * dgy

   elseif reg == 4

      ak = d.tht1
      th = muladd(ak, xi2, π / 4)

      _, _, dgx, dgy = gamder(d, th)

      dvx = (1.0 - xi1) * d.R * (rp.e₃x - 1.0 / 4) + xi1 * ak * dgx
      dvy = (1.0 - xi1) * d.R * (rp.e₃y - 1.0 / 4) + xi1 * ak * dgy

   elseif reg == 5

      ak = π / 4 - d.tht1 - d.tht2
      th = muladd(ak, xi2, π / 4 + d.tht1)

      _, _, dgx, dgy = gamder(d, th)

      dvx = (1.0 - xi1) * d.L2 * rp.q₁x * invn + xi1 * ak * dgx
      dvy = (1.0 - xi1) * d.L2 * rp.q₁y * invn + xi1 * ak * dgy

   elseif reg == 6

      ak = d.tht2 + d.tht3
      th = muladd(ak, xi2, π / 2 - d.tht2)

      _, _, dgx, dgy = gamder(d, th)

      dvx = -(1.0 - xi1) * d.L1 * rp.e₁x * invn + xi1 * ak * dgx
      dvy = -(1.0 - xi1) * d.L1 * rp.e₁y * invn + xi1 * ak * dgy

   elseif reg == 7

      ak = π / 4 - d.tht1 - d.tht3
      th = muladd(ak, xi2, π / 2 + d.tht3)

      _, _, dgx, dgy = gamder(d, th)

      dvx = -(1.0 - xi1) * d.L2 * rp.q₁x * invn + xi1 * ak * dgx
      dvy = -(1.0 - xi1) * d.L2 * rp.q₁y * invn + xi1 * ak * dgy

   elseif reg == 8

      ak = d.tht1
      th = muladd(ak, xi2, 3π / 4 - d.tht1)

      _, _, dgx, dgy = gamder(d, th)

      dvx = (1.0 - xi1) * d.R * (1.0 / 4 - rp.e₃x) + xi1 * ak * dgx
      dvy = (1.0 - xi1) * d.R * (1.0 / 4 - rp.e₃y) + xi1 * ak * dgy

   elseif reg == 9

      ak = d.tht1
      th = muladd(ak, xi2, 3π / 4)

      _, _, dgx, dgy = gamder(d, th)

      dvx = (1.0 - xi1) * d.R * (rp.e₄x - 1.0 / 4) + xi1 * ak * dgx
      dvy = (1.0 - xi1) * d.R * (rp.e₄y - 1.0 / 4) + xi1 * ak * dgy

   elseif reg == 10

      ak = π / 4 - d.tht1 - d.tht3
      th = muladd(ak, xi2, 3π / 4 + d.tht1)

      _, _, dgx, dgy = gamder(d, th)

      dvx = (1.0 - xi1) * d.L2 * rp.q₂x * invn + xi1 * ak * dgx
      dvy = (1.0 - xi1) * d.L2 * rp.q₂y * invn + xi1 * ak * dgy

   else

      throw(ArgumentError("Dwall is defined only for regions 1–10; got reg=$reg"))

   end

   return αv * dvx, αv * dvy

end

"""
   refine!(d::bean, Nc::Int, Nt::Int, K::Vector{Int})

Subdivide each patch in `K` into `Nc * Nt` subpatches, split in `c` and `t`.
Updates `d.pths`, `d.Npat`, recomputes `d.kd`, and rebuilds `d.Qpts` / `d.Qptsbd`.

Notes
- Patches are re-sorted lexicographically by `(reg, ck0, ck1, tk0, tk1)`.
- `Qpts` uses `boundquad!(V, d, k)`.
- `Qptsbd` uses `boundquadbd!(V, d, k)` for patches in `d.kd`.
"""
function refine!(d::bean, Nc::Int, Nt::Int, K::Vector{Int})

   @assert Nc ≥ 1 && Nt ≥ 1 "Nc and Nt must be ≥ 1"
   @assert !isempty(K) "K (set of patches to refine) is empty"

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

   d.pths = vcat([d.pths[i] for i in 1:d.Npat if !(i in K)], children)
   sort!(d.pths, by = q -> (q.reg, q.ck0, q.ck1, q.tk0, q.tk1))

   d.Npat = length(d.pths)

   d.kd = [
      k for k in 1:d.Npat
      if d.pths[k].reg != 11 && d.pths[k].reg != 12 && d.pths[k].ck1 == 1.0
   ]

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

function Base.show(io::IO, d::bean)

   println(io, "bean with properties:")
   println(io, "  (A, B)   = (", d.A, ", ", d.B, ")")
   println(io, "  R        = ", d.R)
   println(io, "  (L1, L2) = (", d.L1, ", ", d.L2, ")")
   println(io, "  tht1     = ", d.tht1)
   println(io, "  tht2     = ", d.tht2)
   println(io, "  tht3     = ", d.tht3)
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