"""
Mutable struct squircle

A  — x-coordinate of squircle center (Float64)

B  — y-coordinate of squircle center (Float64)

P  — A parameter bigger than 2. By default 4

θ₁ — A parameter used to constuct the patches in squircle. (Float64)

R₁ — Scale of the squircle in x axis (Float64)

R₂ — Scale of the squircle in y axis (Float64)

L₁ — Length of the y axis side of the rectangle (Float64)

L₂ — Length of the x axis side of the rectangle (Float64)

Together, L1×L2 are the dimensions of the inscribed rectangle. 

nₕ — number of holes (nonnegative Int) 

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
mutable struct squircle <: abstractdomain

   A::Float64
   B::Float64
   P::Float64
   R1::Float64
   R2::Float64
   L1::Float64
   L2::Float64
   nh::Int
   kd::Vector{Int}
   tht1::Float64
   Npat::Int
   pths::Vector{Patch}
   Qpts::Matrix{Float64}
   Qptsbd::Matrix{Float64}

   function squircle(; b,
      a=nothing, A=nothing, B=nothing,
      P=nothing, R1=nothing, R2=nothing,
      L1=nothing, L2=nothing,
      tht1=nothing,
      ck=nothing, tk=nothing)

      """
      Required:
        - b :: Vector{Int} of length 5, positive entries

      Optional:
        - a     :: Vector{Int} length 5, defaults to ceil.(2*b/3)
        - A, B  :: center of the squircle, defaults to 0.0
        - P     :: squircle exponent, defaults to 4.0, should be > 2
        - R1    :: x-axis scale, defaults to 1.0
        - R2    :: y-axis scale, defaults to 1.0
        - L1    :: inscribed rectangle y-side length, defaults to R2
        - L2    :: inscribed rectangle x-side length, defaults to R1
        - tht1  :: patch construction angle, defaults to π/4
        - ck    :: vector-of-vectors, ck[k] length a[k] + 1
        - tk    :: vector-of-vectors, tk[k] length b[k] + 1
      """

      @assert isa(b, AbstractVector{<:Int}) && length(b) == 5
      @assert all(b .> 0)

      nh = 0

      a = isnothing(a) ? ceil.(Int, 2 .* b ./ 3) : collect(Int, a)

      @assert length(a) == 5
      @assert all(a .> 0)

      A = Float64(something(A, 0.0))
      B = Float64(something(B, 0.0))
      P = Float64(something(P, 4.0))
      R1 = Float64(something(R1, 1.0))
      R2 = Float64(something(R2, 1.0))
      tht1 = Float64(something(tht1, π / 4))

      @assert P >= 2.0
      @assert R1 > 0.0
      @assert R2 > 0.0

      L1 = Float64(something(L1, R2))
      L2 = Float64(something(L2, R1))

      @assert L1 > 0.0
      @assert L2 > 0.0

      Npat = dot(a, b)

      # ck[k] has length a[k] + 1
      # tk[k] has length b[k] + 1
      ck = something(ck, [(0:a[k]) ./ a[k] for k in 1:5])
      tk = something(tk, [(0:b[k]) ./ b[k] for k in 1:5])

      @assert length(ck) == 5
      @assert length(tk) == 5

      pths = Vector{Patch}(undef, Npat)

      # Upper bound; shrink later.
      kd = Vector{Int}(undef, Npat)
      nkd = 0

      idx = 1

      @inbounds for k in 1:5
         ak = a[k]
         bk = b[k]
         ckₖ = ck[k]
         tkₖ = tk[k]

         @assert length(ckₖ) == ak + 1
         @assert length(tkₖ) == bk + 1

         for i in 1:ak
            ck0 = Float64(ckₖ[i])
            ck1 = Float64(ckₖ[i+1])

            for j in 1:bk
               tk0 = Float64(tkₖ[j])
               tk1 = Float64(tkₖ[j+1])

               p = Patch(k, ck0, ck1, tk0, tk1)
               pths[idx] = p

               if p.reg != 5 && p.ck1 == 1.0
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

      d = new(A, B, P, R1, R2, L1, L2, nh, kd, tht1, Npat, pths, Qpts, Qptsbd)

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

function mapx(d::squircle, u::Float64, v::Float64, k::Int)::Float64

   # p denotes information of kth patch
   p = d.pths[k]

   # Difference in c values
   hc = p.ck1 - p.ck0

   xi1 = hc * u / 2 + (p.ck1 + p.ck0) / 2

   # Difference in t values
   ht = p.tk1 - p.tk0

   xi2 = ht * v / 2 + (p.tk1 + p.tk0) / 2

   if p.reg != 5

      X, thet = if p.reg == 1

         d.A + d.L2 / 2, xi2 * 2 * d.tht1 - d.tht1

      elseif p.reg == 2

         d.A + d.L2 * (1 - 2 * xi2) / 2,
         xi2 * (π - 2 * d.tht1) + d.tht1

      elseif p.reg == 3

         d.A - d.L2 / 2,
         xi2 * 2 * d.tht1 + π - d.tht1

      else # p.reg == 4

         d.A - d.L2 * (1 - 2 * xi2) / 2,
         xi2 * (π - 2 * d.tht1) + π + d.tht1

      end

      st, ct = sincos(thet)

      rho = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

      Y = d.A + rho * ct

      return (1 - xi1) * X + xi1 * Y

   else

      xi1 = d.L2 >= d.L1 ? ht * u / 2 + (p.tk1 + p.tk0) / 2 : xi1

      return d.A - d.L2 / 2 + d.L2 * xi1

   end

end

function mapy(d::squircle, u::Float64, v::Float64, k::Int)::Float64

   # p denotes information of kth patch
   p = d.pths[k]

   # Difference in c values
   hc = p.ck1 - p.ck0

   xi1 = hc * u / 2 + (p.ck1 + p.ck0) / 2

   # Difference in t values
   ht = p.tk1 - p.tk0

   xi2 = ht * v / 2 + (p.tk1 + p.tk0) / 2

   if p.reg != 5

      X, thet = if p.reg == 1

         d.B + d.L1 * (2 * xi2 - 1) / 2,
         xi2 * 2 * d.tht1 - d.tht1

      elseif p.reg == 2

         d.B + d.L1 / 2,
         xi2 * (π - 2 * d.tht1) + d.tht1

      elseif p.reg == 3

         d.B - d.L1 * (2 * xi2 - 1) / 2,
         xi2 * 2 * d.tht1 + π - d.tht1

      else # p.reg == 4

         d.B - d.L1 / 2,
         xi2 * (π - 2 * d.tht1) + π + d.tht1

      end

      st, ct = sincos(thet)

      rho = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

      Y = d.B + rho * st

      return (1 - xi1) * X + xi1 * Y

   else

      xi2 = d.L2 >= d.L1 ? hc * v / 2 + (p.ck1 + p.ck0) / 2 : xi2

      return d.B - d.L1 / 2 + d.L1 * xi2

   end

end

function mapxy(d::squircle, u::Float64, v::Float64, k::Int)::Tuple{Float64, Float64}

   return mapx(d, u, v, k), mapy(d, u, v, k)

end

function mapxy!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
   d::squircle, u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)

   # @assert size(Zx) == size(Zy) == size(u) == size(v)

   p = d.pths[k]
   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   if p.reg == 1

      αu = hc / 2
      αv = ht / 2
      βu = p.ck0 + αu
      βv = p.tk0 + αv

      @inbounds for I in eachindex(u, v, Zx, Zy)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + d.L2 / 2
         Xy = d.B + d.L1 * (2xi2 - 1) / 2

         thet = xi2 * 2 * d.tht1 - d.tht1

         st, ct = sincos(thet)

         rho = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

         Yx = muladd(rho, ct, d.A)
         Yy = muladd(rho, st, d.B)

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)

      end

   elseif p.reg == 2

      αu = hc / 2
      αv = ht / 2
      βu = p.ck0 + αu
      βv = p.tk0 + αv

      @inbounds for I in eachindex(u, v, Zx, Zy)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A + d.L2 * (1 - 2xi2) / 2
         Xy = d.B + d.L1 / 2

         thet = xi2 * (π - 2 * d.tht1) + d.tht1

         st, ct = sincos(thet)

         rho = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

         Yx = muladd(rho, ct, d.A)
         Yy = muladd(rho, st, d.B)

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)

      end

   elseif p.reg == 3

      αu = hc / 2
      αv = ht / 2
      βu = p.ck0 + αu
      βv = p.tk0 + αv

      @inbounds for I in eachindex(u, v, Zx, Zy)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A - d.L2 / 2
         Xy = d.B - d.L1 * (2xi2 - 1) / 2

         thet = xi2 * 2 * d.tht1 + π - d.tht1

         st, ct = sincos(thet)

         rho = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

         Yx = muladd(rho, ct, d.A)
         Yy = muladd(rho, st, d.B)

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)

      end

   elseif p.reg == 4

      αu = hc / 2
      αv = ht / 2
      βu = p.ck0 + αu
      βv = p.tk0 + αv

      @inbounds for I in eachindex(u, v, Zx, Zy)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Xx = d.A - d.L2 * (1 - 2xi2) / 2
         Xy = d.B - d.L1 / 2

         thet = xi2 * (π - 2 * d.tht1) + π + d.tht1

         st, ct = sincos(thet)

         rho = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

         Yx = muladd(rho, ct, d.A)
         Yy = muladd(rho, st, d.B)

         Zx[I] = muladd(xi1, Yx - Xx, Xx)
         Zy[I] = muladd(xi1, Yy - Xy, Xy)

      end

   else # p.reg == 5

      if d.L2 >= d.L1

         αu = ht / 2
         αv = hc / 2
         βu = p.tk0 + αu
         βv = p.ck0 + αv

         @inbounds for I in eachindex(u, v, Zx, Zy)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            Zx[I] = d.A - d.L2 / 2 + d.L2 * xi1
            Zy[I] = d.B - d.L1 / 2 + d.L1 * xi2

         end

      else

         αu = hc / 2
         αv = ht / 2
         βu = p.ck0 + αu
         βv = p.tk0 + αv

         @inbounds for I in eachindex(u, v, Zx, Zy)

            xi1 = muladd(αu, u[I], βu)
            xi2 = muladd(αv, v[I], βv)

            Zx[I] = d.A - d.L2 / 2 + d.L2 * xi1
            Zy[I] = d.B - d.L1 / 2 + d.L1 * xi2

         end

      end

   end

   return nothing

end

function draw(d::squircle, flag=nothing; L::Int=32, show::Bool=true)

  colors = (
    RGBf(1, 0, 0),         # red
    RGBf(1, 0.5, 0.25),    # orange
    RGBf(0.58, 0, 0.82),   # purple
    RGBf(0, 1, 0),         # green
    RGBf(0, 0, 1),         # blue
  )

  return draw_geom(d, colors; flag=flag, L=L, show=show)
end

#-----------------------
"""
   gamx(d::squircle, t::Float64, k::Int)
   gamx!(out::StridedArray{Float64}, d::squircle, t::StridedArray{Float64}, k::Int)

First x-coordinate of the boundary parametrization gamma(t).

- `gamx(d, t, k)` computes the right boundary of patch `k`, where `t in [-1,1]`.

`t` may be a scalar or an array; the return/output has the same shape.
"""
function gamx(d::squircle, t::Float64, k::Int)::Float64

   p = d.pths[k]

   xi2 = (p.tk1 - p.tk0) * t / 2 + (p.tk0 + p.tk1) / 2

   theta = if p.reg == 1

      xi2 * 2 * d.tht1 - d.tht1

   elseif p.reg == 2

      xi2 * (pi - 2 * d.tht1) + d.tht1

   elseif p.reg == 3

      xi2 * 2 * d.tht1 + pi - d.tht1

   else # p.reg == 4

      xi2 * (pi - 2 * d.tht1) + pi + d.tht1

   end

   st, ct = sincos(theta)

   hP = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

   return d.A + hP * ct

end

function gamx!(out::StridedArray{Float64}, d::squircle,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   if p.reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

   elseif p.reg == 2

      Δth = pi - 2 * d.tht1
      th0 = d.tht1

   elseif p.reg == 3

      Δth = 2 * d.tht1
      th0 = pi - d.tht1

   else # p.reg == 4

      Δth = pi - 2 * d.tht1
      th0 = pi + d.tht1

   end

   @inbounds for I in eachindex(t, out)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(Δth, xi2, th0)

      st, ct = sincos(theta)

      hP = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

      out[I] = muladd(hP, ct, d.A)

   end

   return nothing

end

"""
   gamy(d::squircle, t::Float64, k::Int)
   gamy!(out::StridedArray{Float64}, d::squircle, t::StridedArray{Float64}, k::Int)

Second y-coordinate of the boundary parametrization gamma(t).

- `gamy(d, t, k)` computes the right boundary of patch `k`, where `t in [-1,1]`.

`t` may be a scalar or an array; the return/output has the same shape.
"""
function gamy(d::squircle, t::Float64, k::Int)::Float64

   p = d.pths[k]

   xi2 = (p.tk1 - p.tk0) * t / 2 + (p.tk0 + p.tk1) / 2

   theta = if p.reg == 1

      xi2 * 2 * d.tht1 - d.tht1

   elseif p.reg == 2

      xi2 * (pi - 2 * d.tht1) + d.tht1

   elseif p.reg == 3

      xi2 * 2 * d.tht1 + pi - d.tht1

   else # p.reg == 4

      xi2 * (pi - 2 * d.tht1) + pi + d.tht1

   end

   st, ct = sincos(theta)

   hP = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

   return d.B + hP * st

end

function gamy!(out::StridedArray{Float64}, d::squircle,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   if p.reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

   elseif p.reg == 2

      Δth = pi - 2 * d.tht1
      th0 = d.tht1

   elseif p.reg == 3

      Δth = 2 * d.tht1
      th0 = pi - d.tht1

   else # p.reg == 4

      Δth = pi - 2 * d.tht1
      th0 = pi + d.tht1

   end

   @inbounds for I in eachindex(t, out)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(Δth, xi2, th0)

      st, ct = sincos(theta)

      hP = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

      out[I] = muladd(hP, st, d.B)

   end

   return nothing

end

function gam(d::squircle, t::Float64, k::Int)::Tuple{Float64, Float64}

   p = d.pths[k]

   xi2 = (p.tk1 - p.tk0) * t / 2 + (p.tk0 + p.tk1) / 2

   theta = if p.reg == 1

      xi2 * 2 * d.tht1 - d.tht1

   elseif p.reg == 2

      xi2 * (pi - 2 * d.tht1) + d.tht1

   elseif p.reg == 3

      xi2 * 2 * d.tht1 + pi - d.tht1

   else # p.reg == 4

      xi2 * (pi - 2 * d.tht1) + pi + d.tht1

   end

   st, ct = sincos(theta)

   hP = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

   outx = muladd(hP, ct, d.A)
   outy = muladd(hP, st, d.B)

   return outx, outy

end

function gam!(out::Vector{Float64}, d::squircle, t::Float64, k::Int)

   p = d.pths[k]

   xi2 = (p.tk1 - p.tk0) * t / 2 + (p.tk0 + p.tk1) / 2

   theta = if p.reg == 1

      xi2 * 2 * d.tht1 - d.tht1

   elseif p.reg == 2

      xi2 * (pi - 2 * d.tht1) + d.tht1

   elseif p.reg == 3

      xi2 * 2 * d.tht1 + pi - d.tht1

   else # p.reg == 4

      xi2 * (pi - 2 * d.tht1) + pi + d.tht1

   end

   st, ct = sincos(theta)

   hP = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

   out[1] = muladd(hP, ct, d.A)
   out[2] = muladd(hP, st, d.B)

   return nothing

end

function gam!(out::Matrix{Float64}, d::squircle, t::Vector{Float64}, k::Int)

   # out is a Matrix with 2 rows and length(t) columns

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   if p.reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

   elseif p.reg == 2

      Δth = pi - 2 * d.tht1
      th0 = d.tht1

   elseif p.reg == 3

      Δth = 2 * d.tht1
      th0 = pi - d.tht1

   else # p.reg == 4

      Δth = pi - 2 * d.tht1
      th0 = pi + d.tht1

   end

   @inbounds for I in eachindex(t)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(Δth, xi2, th0)

      st, ct = sincos(theta)

      hP = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

      out[1, I] = muladd(hP, ct, d.A)
      out[2, I] = muladd(hP, st, d.B)

   end

   return nothing

end

function drawbd(d::squircle, flag = true; L::Int = 33, show::Bool = true)

    colors = (
        RGBf(1, 0, 0),         # red
        RGBf(1, 0.5, 0.25),    # orange
        RGBf(0.58, 0, 0.82),   # purple
        RGBf(0, 1, 0),         # green
        RGBf(0, 0, 1),         # blue
    )

    return drawbd_geom(d, colors; flag = flag, L = L, show = show)
end

"""
   dgamx(d::squircle, t::Float64, k::Int) -> Float64
   dgamx!(out::StridedArray{Float64}, d::squircle, t::StridedArray{Float64}, k::Int)

Derivative of the first coordinate of the boundary parametrization.

- With `k`: derivative of the right boundary of patch `k`,
  i.e. derivative of `gamx(d, t, k)` with respect to `t in [-1,1]`.
"""
function dgamx(d::squircle, t::Float64, k::Int)::Float64

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   if p.reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

   elseif p.reg == 2

      Δth = pi - 2 * d.tht1
      th0 = d.tht1

   elseif p.reg == 3

      Δth = 2 * d.tht1
      th0 = pi - d.tht1

   else

      Δth = pi - 2 * d.tht1
      th0 = pi + d.tht1

   end

   R1P = d.R1^d.P
   R2P = d.R2^d.P

   xi2 = muladd(αt, t, βt)
   theta = muladd(Δth, xi2, th0)

   st, ct = sincos(theta)

   ct1 = ct * abs(ct)^(d.P - 2) / R1P
   st1 = st * abs(st)^(d.P - 2) / R2P

   gg = ct1 * ct + st1 * st
   dgg = st * ct1 - ct * st1

   h = gg^(-1 / d.P)
   dh = h * dgg / gg

   dgx = dh * ct - h * st

   return αt * Δth * dgx

end

function dgamx!(out::StridedArray{Float64}, d::squircle,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   if p.reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

   elseif p.reg == 2

      Δth = pi - 2 * d.tht1
      th0 = d.tht1

   elseif p.reg == 3

      Δth = 2 * d.tht1
      th0 = pi - d.tht1

   else

      Δth = pi - 2 * d.tht1
      th0 = pi + d.tht1

   end

   scale = αt * Δth
   R1P = d.R1^d.P
   R2P = d.R2^d.P
   invP = -1 / d.P

   @inbounds for I in eachindex(t, out)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(Δth, xi2, th0)

      st, ct = sincos(theta)

      ct1 = ct * abs(ct)^(d.P - 2) / R1P
      st1 = st * abs(st)^(d.P - 2) / R2P

      gg = ct1 * ct + st1 * st
      dgg = st * ct1 - ct * st1

      h = gg^invP
      dh = h * dgg / gg

      dgx = dh * ct - h * st

      out[I] = scale * dgx

   end

   return nothing

end

"""
   dgamy(d::squircle, t::Float64, k::Int) -> Float64
   dgamy!(out::StridedArray{Float64}, d::squircle, t::StridedArray{Float64}, k::Int)

Derivative of the second coordinate of the boundary parametrization.
"""
function dgamy(d::squircle, t::Float64, k::Int)::Float64

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   if p.reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

   elseif p.reg == 2

      Δth = pi - 2 * d.tht1
      th0 = d.tht1

   elseif p.reg == 3

      Δth = 2 * d.tht1
      th0 = pi - d.tht1

   else

      Δth = pi - 2 * d.tht1
      th0 = pi + d.tht1

   end

   R1P = d.R1^d.P
   R2P = d.R2^d.P

   xi2 = muladd(αt, t, βt)
   theta = muladd(Δth, xi2, th0)

   st, ct = sincos(theta)

   ct1 = ct * abs(ct)^(d.P - 2) / R1P
   st1 = st * abs(st)^(d.P - 2) / R2P

   gg = ct1 * ct + st1 * st
   dgg = st * ct1 - ct * st1

   h = gg^(-1 / d.P)
   dh = h * dgg / gg

   dgy = dh * st + h * ct

   return αt * Δth * dgy

end

function dgamy!(out::StridedArray{Float64}, d::squircle,
   t::StridedArray{Float64}, k::Int)

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   if p.reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

   elseif p.reg == 2

      Δth = pi - 2 * d.tht1
      th0 = d.tht1

   elseif p.reg == 3

      Δth = 2 * d.tht1
      th0 = pi - d.tht1

   else

      Δth = pi - 2 * d.tht1
      th0 = pi + d.tht1

   end

   scale = αt * Δth
   R1P = d.R1^d.P
   R2P = d.R2^d.P
   invP = -1 / d.P

   @inbounds for I in eachindex(t, out)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(Δth, xi2, th0)

      st, ct = sincos(theta)

      ct1 = ct * abs(ct)^(d.P - 2) / R1P
      st1 = st * abs(st)^(d.P - 2) / R2P

      gg = ct1 * ct + st1 * st
      dgg = st * ct1 - ct * st1

      h = gg^invP
      dh = h * dgg / gg

      dgy = dh * st + h * ct

      out[I] = scale * dgy

   end

   return nothing

end

function dgam(d::squircle, t::Float64, k::Int)::Tuple{Float64, Float64}

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   if p.reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

   elseif p.reg == 2

      Δth = pi - 2 * d.tht1
      th0 = d.tht1

   elseif p.reg == 3

      Δth = 2 * d.tht1
      th0 = pi - d.tht1

   else

      Δth = pi - 2 * d.tht1
      th0 = pi + d.tht1

   end

   R1P = d.R1^d.P
   R2P = d.R2^d.P

   xi2 = muladd(αt, t, βt)
   theta = muladd(Δth, xi2, th0)

   st, ct = sincos(theta)

   ct1 = ct * abs(ct)^(d.P - 2) / R1P
   st1 = st * abs(st)^(d.P - 2) / R2P

   gg = ct1 * ct + st1 * st
   dgg = st * ct1 - ct * st1

   h = gg^(-1 / d.P)
   dh = h * dgg / gg

   scale = αt * Δth

   outx = scale * (dh * ct - h * st)
   outy = scale * (dh * st + h * ct)

   return outx, outy

end

function dgam!(out::Vector{Float64}, d::squircle, t::Float64, k::Int)

   dx, dy = dgam(d, t, k)

   out[1] = dx
   out[2] = dy

   return nothing

end

function dgam!(out::Matrix{Float64}, d::squircle, t::Vector{Float64}, k::Int)

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   if p.reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

   elseif p.reg == 2

      Δth = pi - 2 * d.tht1
      th0 = d.tht1

   elseif p.reg == 3

      Δth = 2 * d.tht1
      th0 = pi - d.tht1

   else

      Δth = pi - 2 * d.tht1
      th0 = pi + d.tht1

   end

   scale = αt * Δth
   R1P = d.R1^d.P
   R2P = d.R2^d.P
   invP = -1 / d.P

   @inbounds for I in eachindex(t)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(Δth, xi2, th0)

      st, ct = sincos(theta)

      ct1 = ct * abs(ct)^(d.P - 2) / R1P
      st1 = st * abs(st)^(d.P - 2) / R2P

      gg = ct1 * ct + st1 * st
      dgg = st * ct1 - ct * st1

      h = gg^invP
      dh = h * dgg / gg

      out[1, I] = scale * (dh * ct - h * st)
      out[2, I] = scale * (dh * st + h * ct)

   end

   return nothing

end

function gamp(d::squircle, t::Float64, k::Int)::Tuple{Float64,Float64}

   dx, dy = dgam(d, t, k)

   return dy, -dx

end

function gamp!(out::Matrix{Float64}, d::squircle, t::Vector{Float64}, k::Int)

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   if p.reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

   elseif p.reg == 2

      Δth = pi - 2 * d.tht1
      th0 = d.tht1

   elseif p.reg == 3

      Δth = 2 * d.tht1
      th0 = pi - d.tht1

   else

      Δth = pi - 2 * d.tht1
      th0 = pi + d.tht1

   end

   scale = αt * Δth
   R1P = d.R1^d.P
   R2P = d.R2^d.P
   invP = -1 / d.P

   @inbounds for I in eachindex(t)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(Δth, xi2, th0)

      st, ct = sincos(theta)

      ct1 = ct * abs(ct)^(d.P - 2) / R1P
      st1 = st * abs(st)^(d.P - 2) / R2P

      gg = ct1 * ct + st1 * st
      dgg = st * ct1 - ct * st1

      h = gg ^ invP
      dh = h * dgg / gg

      dx = scale * (dh * ct - h * st)
      dy = scale * (dh * st + h * ct)

      out[1, I] = dy
      out[2, I] = -dx

   end

   return nothing

end

function nu!(out::Matrix{Float64}, d::squircle, t::Vector{Float64}, k::Int)

   p = d.pths[k]

   αt = (p.tk1 - p.tk0) / 2
   βt = p.tk0 + αt

   if p.reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

   elseif p.reg == 2

      Δth = pi - 2 * d.tht1
      th0 = d.tht1

   elseif p.reg == 3

      Δth = 2 * d.tht1
      th0 = pi - d.tht1

   else

      @assert p.reg == 4 "nu! is only defined for boundary regions 1:4"
      Δth = pi - 2 * d.tht1
      th0 = pi + d.tht1

   end

   scale = αt * Δth
   R1P = d.R1^d.P
   R2P = d.R2^d.P
   invP = -1 / d.P

   @inbounds for I in eachindex(t)

      xi2 = muladd(αt, t[I], βt)
      theta = muladd(Δth, xi2, th0)

      st, ct = sincos(theta)

      ct1 = ct * abs(ct)^(d.P - 2) / R1P
      st1 = st * abs(st)^(d.P - 2) / R2P

      gg = ct1 * ct + st1 * st
      dgg = st * ct1 - ct * st1

      h = gg^invP
      dh = h * dgg / gg

      dx = scale * (dh * ct - h * st)
      dy = scale * (dh * st + h * ct)

      # gamp = (dy, -dx)
      S = hypot(dy, dx)

      out[1, I] = dy / S
      out[2, I] = -dx / S

   end

   return nothing

end

# This function finds s such that γ_l(t) = γ_k(s),
# allowing s to lie outside [-1, 1].
function bdinv(d::squircle, t::Float64, l::Int, k::Int)::Float64
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
         xil * (π - 2 * d.tht1) + d.tht1
      elseif plr == 3
         xil * 2 * d.tht1 + π - d.tht1
      elseif plr == 4
         xil * (π - 2 * d.tht1) + π + d.tht1
      else
         throw(ArgumentError("bdinv source region must be 1–4; got reg=$plr"))
      end

      # Get angular map for target patch k:
      #     theta = ak * xi + bk
      ak, bk = if pkr == 1
         2 * d.tht1, -d.tht1
      elseif pkr == 2
         π - 2 * d.tht1, d.tht1
      elseif pkr == 3
         2 * d.tht1, π - d.tht1
      elseif pkr == 4
         π - 2 * d.tht1, π + d.tht1
      else
         throw(ArgumentError("bdinv target region must be 1–4; got reg=$pkr"))
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
   hderhigher(d::squircle, θ::Float64)

Compute h, h', h'', h⁽³⁾, h⁽⁴⁾, h⁽⁵⁾, h⁽⁶⁾ for the squircle radial function

   h(θ) = ( |cos(θ)|^P / R1^P + |sin(θ)|^P / R2^P )^(-1/P)
"""
function hderhigher(d::squircle, θ::Float64)

   P = d.P

   st, ct = sincos(θ)

   R1P = d.R1^P
   R2P = d.R2^P

   if P == 2.0

      A = 1.0 / R1P
      B = 1.0 / R2P
      Δ = B - A

      g = A * ct^2 + B * st^2
      g1 = 2.0 * st * ct * Δ
      g2 = 2.0 * (ct^2 - st^2) * Δ
      g3 = -8.0 * st * ct * Δ
      g4 = -8.0 * (ct^2 - st^2) * Δ
      g5 = 32.0 * st * ct * Δ
      g6 = 32.0 * (ct^2 - st^2) * Δ

   elseif P == 4.0

      A = 1.0 / R1P
      B = 1.0 / R2P

      st2 = st * st
      ct2 = ct * ct
      st4 = st2 * st2
      ct4 = ct2 * ct2

      g = A * ct4 + B * st4

      g1 = -4.0 * st * ct * (A * ct2 - B * st2)

      g2 = -4.0 * (
         A * ct4 +
         B * st4 -
         3.0 * (A + B) * st2 * ct2
      )

      g3 = 8.0 * st * ct * (
         A * (5.0 * ct2 - 3.0 * st2) +
         B * (3.0 * ct2 - 5.0 * st2)
      )

      g4 = 8.0 * (
         A * (3.0 * st4 - 24.0 * st2 * ct2 + 5.0 * ct4) +
         B * (5.0 * st4 - 24.0 * st2 * ct2 + 3.0 * ct4)
      )

      g5 = -32.0 * st * ct * (
         A * (17.0 * ct2 - 15.0 * st2) +
         B * (15.0 * ct2 - 17.0 * st2)
      )

      g6 = -32.0 * (
         A * (128.0 * st4 - 130.0 * st2 + 17.0) +
         B * (128.0 * st4 - 126.0 * st2 + 15.0)
      )

   else

      ct0 = abs(ct)^(P - 6) / R1P
      ct1 = abs(ct)^(P - 4) / R1P
      ct2 = ct * ct * ct1
      ct3 = ct * ct * ct2

      st0 = abs(st)^(P - 6) / R2P
      st1 = abs(st)^(P - 4) / R2P
      st2 = st * st * st1
      st3 = st * st * st2

      ct_sq = ct * ct
      st_sq = st * st

      g = ct3 + st3

      g1 = -P * st * ct * (ct2 - st2)

      g2 = P * (P - 1) * (ct2 + st2) - P^2 * g

      g3 = -P * (P - 1) * (P - 2) * st * ct * (ct1 - st1) -
           P^2 * g1

      g4 = P * (P - 1) * (P - 2) * (P - 3) * (ct1 + st1) -
           P * (P - 1) * (P - 2)^2 * (ct2 + st2) -
           P^2 * g2

      g5 = P * (P - 1) * (P - 2) * (P - 3) * (P - 4) *
           st * ct * (st0 - ct0) -
           P * (P - 1) * (P - 2)^3 *
           st * ct * (st1 - ct1) -
           P^2 * g3

      C1 = P * (P - 1) * (P - 2) * (P - 3) * (P - 4)
      C2 = P * (P - 1) * (P - 2)^3

      F0 = (ct_sq - st_sq) * (st0 - ct0) +
           (P - 6) * (ct_sq * st0 + st_sq * ct0)

      F1 = (ct_sq - st_sq) * (st1 - ct1) +
           (P - 4) * (ct_sq * st1 + st_sq * ct1)

      g6 = C1 * F0 - C2 * F1 - P^2 * g4

   end

   h = g^(-1 / P)

   h1 = -(h * g1) / (P * g)

   h2 = -((P + 1) * h1 * g1 + h * g2) / (P * g)

   h3 = -((2P + 1) * h2 * g1 +
          (P + 2) * h1 * g2 +
          h * g3) / (P * g)

   h4 = -((3P + 1) * h3 * g1 +
          (3P + 3) * h2 * g2 +
          (P + 3) * h1 * g3 +
          h * g4) / (P * g)

   h5 = -((4P + 1) * h4 * g1 +
          (6P + 4) * h3 * g2 +
          (4P + 6) * h2 * g3 +
          (P + 4) * h1 * g4 +
          h * g5) / (P * g)

   h6 = -((5P + 1) * h5 * g1 +
          (10P + 5) * h4 * g2 +
          (10P + 10) * h3 * g3 +
          (5P + 10) * h2 * g4 +
          (P + 5) * h1 * g5 +
          h * g6) / (P * g)

   return h, h1, h2, h3, h4, h5, h6

end

"""
    DLP(d::squircle, t::Float64, tau::StridedArray{Float64}, k::Int)

Double-layer kernel on the boundary:
K(t, τ) = ((γ_k(τ) - γ_k(t))⋅  γᵖᵉʳᵖ_k(τ)) / ‖γ_k(τ) - γ_k(t)‖² 
and for k == l the limiting value is taken for patch k.
The array method returns an array with the ***same shape*** as `tau`.
"""
function DLP!(out::StridedArray{Float64}, d::squircle, t::Float64,
   tau::StridedArray{Float64}, k::Int, x::Vector{Float64},
   G::Vector{Float64}, GP::Vector{Float64})

   p = d.pths[k]

   ht = p.tk1 - p.tk0
   αt = 0.5 * ht
   βt = p.tk0 + αt

   reg = p.reg

   if reg == 1

      ak = 2 * d.tht1
      bk = -d.tht1

   elseif reg == 2

      ak = π - 2 * d.tht1
      bk = d.tht1

   elseif reg == 3

      ak = 2 * d.tht1
      bk = π - d.tht1

   elseif reg == 4

      ak = π - 2 * d.tht1
      bk = π + d.tht1

   else

      throw(ArgumentError("DLP! self term defined only for boundary regions 1–4; got reg=$reg"))

   end

   P = d.P
   R1P = d.R1^P
   R2P = d.R2^P
   invP = -1 / P

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

      stτ, ctτ = sincos(thτ)

      ct1 = ctτ * abs(ctτ)^(P - 2) / R1P
      st1 = stτ * abs(stτ)^(P - 2) / R2P

      ggτ = ct1 * ctτ + st1 * stτ
      dggτ = stτ * ct1 - ctτ * st1

      hτ = ggτ^invP
      hτ1 = hτ * dggτ / ggτ

      Δ = 0.5 * (thτ - tht)

      sΔ, cΔ = sincos(Δ)

      #This tolerance needs to be higher if P is large
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

         q = (hτ - ht0) / (2 * sΔ)
         q2 = (q - hτ1 * cΔ) / sΔ

      end

      num = ht0 * hτ + hτ * q2 + 2 * q * hτ1 * cΔ
      den = q^2 + hτ * ht0

      out[i] = 0.5 * λ * num / den

   end

   return nothing

end

#-----------------------

"""
   gamder(d::squircle, theta::Float64) -> Tuple{Float64, Float64, Float64, Float64}
   gamder!(out::Vector{Float64}, d::squircle, theta::Float64)

Evaluate the squircle boundary parametrization and its theta-derivative.

Returns `(gx, gy, dgx, dgy)`, where

   gx(theta) = h(theta) * cos(theta)
   gy(theta) = h(theta) * sin(theta)

and `dgx`, `dgy` are derivatives with respect to `theta`.

The center `(d.A, d.B)` is not included here. Add it outside when needed:
   x = d.A + gx
   y = d.B + gy

gamder! writes `(gx, gy, dgx, dgy)` into `out`. (Requires `length(out) >= 4`).
"""
function gamder(d::squircle, theta::Float64)::Tuple{Float64, Float64, Float64, Float64}

   st, ct = sincos(theta)

   R1P = d.R1^d.P
   R2P = d.R2^d.P

   ct1 = ct * abs(ct)^(d.P - 2) / R1P
   st1 = st * abs(st)^(d.P - 2) / R2P

   gg = ct1 * ct + st1 * st
   dgg = st * ct1 - ct * st1

   h = gg^(-1 / d.P)
   dh = h * dgg / gg

   gx = h * ct
   gy = h * st

   dgx = dh * ct - h * st
   dgy = dh * st + h * ct

   return gx, gy, dgx, dgy

end

function gamder!(out::Vector{Float64}, d::squircle, theta::Float64)

   st, ct = sincos(theta)

   R1P = d.R1^d.P
   R2P = d.R2^d.P

   ct1 = ct * abs(ct)^(d.P - 2) / R1P
   st1 = st * abs(st)^(d.P - 2) / R2P

   gg = ct1 * ct + st1 * st
   dgg = st * ct1 - ct * st1

   h = gg^(-1 / d.P)
   dh = h * dgg / gg

   out[1] = h * ct
   out[2] = h * st
   out[3] = dh * ct - h * st
   out[4] = dh * st + h * ct

   return nothing

end

function Dmap!(out::StridedArray{Float64}, d::squircle,
   u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)

   p = d.pths[k]
   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0
   reg = p.reg

   αu = hc / 2
   αv = ht / 2
   βu = p.ck0 + αu      # û = αu * u + βu
   βv = p.tk0 + αv      # v̂ = αv * v + βv
   αuv = αu * αv

   if reg == 5
      # constant over the rectangular patch
      fill!(out, d.L1 * d.L2 * αuv)
      return nothing
   end

   if reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

      L2h = d.L2 / 2
      L1h = d.L1 / 2

      @inbounds for I in eachindex(out, u, v)

         û = muladd(αu, u[I], βu)
         v̂ = muladd(αv, v[I], βv)

         th = muladd(Δth, v̂, th0)

         gx, gy, dgx, dgy = gamder(d, th)

         # X = (L2/2, L1*(2v̂ - 1)/2)
         dux = gx - L2h
         duy = gy - d.L1 * (2v̂ - 1) / 2

         dvx = û * Δth * dgx
         dvy = (1 - û) * d.L1 + û * Δth * dgy

         out[I] = αuv * abs(dux * dvy - dvx * duy)

      end

   elseif reg == 2

      Δth = pi - 2 * d.tht1
      th0 = d.tht1

      L2h = d.L2 / 2
      L1h = d.L1 / 2

      @inbounds for I in eachindex(out, u, v)

         û = muladd(αu, u[I], βu)
         v̂ = muladd(αv, v[I], βv)

         th = muladd(Δth, v̂, th0)

         gx, gy, dgx, dgy = gamder(d, th)

         # X = (L2*(1 - 2v̂)/2, L1/2)
         dux = gx - d.L2 * (1 - 2v̂) / 2
         duy = gy - L1h

         dvx = (û - 1) * d.L2 + û * Δth * dgx
         dvy = û * Δth * dgy

         out[I] = αuv * abs(dux * dvy - dvx * duy)

      end

   elseif reg == 3

      Δth = 2 * d.tht1
      th0 = pi - d.tht1

      L2h = d.L2 / 2

      @inbounds for I in eachindex(out, u, v)

         û = muladd(αu, u[I], βu)
         v̂ = muladd(αv, v[I], βv)

         th = muladd(Δth, v̂, th0)

         gx, gy, dgx, dgy = gamder(d, th)

         # X = (-L2/2, -L1*(2v̂ - 1)/2)
         dux = gx + L2h
         duy = gy + d.L1 * (2v̂ - 1) / 2

         dvx = û * Δth * dgx
         dvy = (û - 1) * d.L1 + û * Δth * dgy

         out[I] = αuv * abs(dux * dvy - dvx * duy)

      end

   elseif reg == 4

      Δth = pi - 2 * d.tht1
      th0 = pi + d.tht1

      L1h = d.L1 / 2

      @inbounds for I in eachindex(out, u, v)

         û = muladd(αu, u[I], βu)
         v̂ = muladd(αv, v[I], βv)

         th = muladd(Δth, v̂, th0)

         gx, gy, dgx, dgy = gamder(d, th)

         # X = (-L2*(1 - 2v̂)/2, -L1/2)
         dux = gx + d.L2 * (1 - 2v̂) / 2
         duy = gy + L1h

         dvx = (1 - û) * d.L2 + û * Δth * dgx
         dvy = û * Δth * dgy

         out[I] = αuv * abs(dux * dvy - dvx * duy)

      end

   else

      throw(ArgumentError("Dmap! region must be 1–5; got reg=$reg"))

   end

   return nothing

end

"""
A combination of mapxy! and Dmap! function (No allocations!)
The purpose of this function is to reduce computations
related to cosine and sine's.
"""
function mapxy_Dmap!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
   DJ::StridedArray{Float64}, d::squircle, u::StridedArray{Float64},
   v::StridedArray{Float64}, k::Int)

   p = d.pths[k]
   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0
   reg = p.reg

   αu = hc / 2
   αv = ht / 2
   βu = p.ck0 + αu
   βv = p.tk0 + αv
   αuv = αu * αv

   L1h = d.L1 / 2
   L2h = d.L2 / 2

   if reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

      Xx = d.A + L2h
      Xy0 = d.B - L1h

      @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

         û = muladd(αu, u[I], βu)
         v̂ = muladd(αv, v[I], βv)

         th = muladd(Δth, v̂, th0)

         gx, gy, dgx, dgy = gamder(d, th)

         Yx = d.A + gx
         Yy = d.B + gy

         Xy = muladd(d.L1, v̂, Xy0)

         Zx[I] = muladd(û, Yx - Xx, Xx)
         Zy[I] = muladd(û, Yy - Xy, Xy)

         dux = gx - L2h
         duy = gy - d.L1 * (2v̂ - 1) / 2

         dvx = û * Δth * dgx
         dvy = (1 - û) * d.L1 + û * Δth * dgy

         DJ[I] = αuv * abs(dux * dvy - dvx * duy)

      end

   elseif reg == 2

      Δth = pi - 2 * d.tht1
      th0 = d.tht1

      Xx0 = d.A + L2h
      Xy = d.B + L1h

      @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

         û = muladd(αu, u[I], βu)
         v̂ = muladd(αv, v[I], βv)

         th = muladd(Δth, v̂, th0)

         gx, gy, dgx, dgy = gamder(d, th)

         Yx = d.A + gx
         Yy = d.B + gy

         Xx = muladd(-d.L2, v̂, Xx0)

         Zx[I] = muladd(û, Yx - Xx, Xx)
         Zy[I] = muladd(û, Yy - Xy, Xy)

         dux = gx - d.L2 * (1 - 2v̂) / 2
         duy = gy - L1h

         dvx = (û - 1) * d.L2 + û * Δth * dgx
         dvy = û * Δth * dgy

         DJ[I] = αuv * abs(dux * dvy - dvx * duy)

      end

   elseif reg == 3

      Δth = 2 * d.tht1
      th0 = pi - d.tht1

      Xx = d.A - L2h
      Xy0 = d.B + L1h

      @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

         û = muladd(αu, u[I], βu)
         v̂ = muladd(αv, v[I], βv)

         th = muladd(Δth, v̂, th0)

         gx, gy, dgx, dgy = gamder(d, th)

         Yx = d.A + gx
         Yy = d.B + gy

         Xy = muladd(-d.L1, v̂, Xy0)

         Zx[I] = muladd(û, Yx - Xx, Xx)
         Zy[I] = muladd(û, Yy - Xy, Xy)

         dux = gx + L2h
         duy = gy + d.L1 * (2v̂ - 1) / 2

         dvx = û * Δth * dgx
         dvy = (û - 1) * d.L1 + û * Δth * dgy

         DJ[I] = αuv * abs(dux * dvy - dvx * duy)

      end

   elseif reg == 4

      Δth = pi - 2 * d.tht1
      th0 = pi + d.tht1

      Xx0 = d.A - L2h
      Xy = d.B - L1h

      @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

         û = muladd(αu, u[I], βu)
         v̂ = muladd(αv, v[I], βv)

         th = muladd(Δth, v̂, th0)

         gx, gy, dgx, dgy = gamder(d, th)

         Yx = d.A + gx
         Yy = d.B + gy

         Xx = muladd(d.L2, v̂, Xx0)

         Zx[I] = muladd(û, Yx - Xx, Xx)
         Zy[I] = muladd(û, Yy - Xy, Xy)

         dux = gx + d.L2 * (1 - 2v̂) / 2
         duy = gy + L1h

         dvx = (1 - û) * d.L2 + û * Δth * dgx
         dvy = û * Δth * dgy

         DJ[I] = αuv * abs(dux * dvy - dvx * duy)

      end

   elseif reg == 5

      fill!(DJ, d.L1 * d.L2 * αuv)

      if d.L2 >= d.L1

         αu = ht / 2
         αv = hc / 2
         βu = p.tk0 + αu
         βv = p.ck0 + αv

      end

      x0 = d.A - L2h
      y0 = d.B - L1h

      @inbounds for I in eachindex(Zx, Zy, DJ, u, v)

         xi1 = muladd(αu, u[I], βu)
         xi2 = muladd(αv, v[I], βv)

         Zx[I] = muladd(d.L2, xi1, x0)
         Zy[I] = muladd(d.L1, xi2, y0)

      end

   else

      throw(ArgumentError("mapxy_Dmap! region must be 1–5; got reg=$reg"))

   end

   return nothing

end

function chk_map(d::squircle; n::Int = 32, tol::Float64 = 5e-14)
    Iex = d.R1 * d.R2 * gamma(1 + 1/d.P)*(
        (d.R1^2 + d.R2^2) * gamma(3.0/d.P) / gamma(4.0/d.P) 
    + 2.0*(d.A^2 + d.B^2) * gamma(1.0/d.P) / gamma(2.0/d.P)
    )
    return chkmap_geom(d, Iex; n = n, tol = tol)
end

"""
    jinvmap(d::squircle, u::Float64, v::Float64, r::Int) -> tuple

Inverse Jacobian of the mapping at a single reference point `(u,v) ∈ [-1,1]^2`
for the given region `r`.

Returns a 4-tuple `(J11, J12, J21, J22)` containing the components of `J⁻¹`.

Only the region number is needed because this is used on the reference region
with `ht = hc = 1`.
"""
function jinvmap(d::squircle, u::Float64, v::Float64, r::Int)

   û = (u + 1) / 2
   v̂ = (v + 1) / 2

   if r == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

      th = muladd(Δth, v̂, th0)

      gx, gy, dgx, dgy = gamder(d, th)

      dux = (gx - d.L2 / 2) / 2
      dvx = û * Δth * dgx / 2

      duy = (gy - d.L1 * (2v̂ - 1) / 2) / 2
      dvy = (1 - û) * d.L1 / 2 + û * Δth * dgy / 2

   elseif r == 2

      Δth = pi - 2 * d.tht1
      th0 = d.tht1

      th = muladd(Δth, v̂, th0)

      gx, gy, dgx, dgy = gamder(d, th)

      dux = (gx - d.L2 * (1 - 2v̂) / 2) / 2
      dvx = (û - 1) * d.L2 / 2 + û * Δth * dgx / 2

      duy = (gy - d.L1 / 2) / 2
      dvy = û * Δth * dgy / 2

   elseif r == 3

      Δth = 2 * d.tht1
      th0 = pi - d.tht1

      th = muladd(Δth, v̂, th0)

      gx, gy, dgx, dgy = gamder(d, th)

      dux = (gx + d.L2 / 2) / 2
      dvx = û * Δth * dgx / 2

      duy = (gy + d.L1 * (2v̂ - 1) / 2) / 2
      dvy = (û - 1) * d.L1 / 2 + û * Δth * dgy / 2

   elseif r == 4

      Δth = pi - 2 * d.tht1
      th0 = pi + d.tht1

      th = muladd(Δth, v̂, th0)

      gx, gy, dgx, dgy = gamder(d, th)

      dux = (gx + d.L2 * (1 - 2v̂) / 2) / 2
      dvx = (1 - û) * d.L2 / 2 + û * Δth * dgx / 2

      duy = (gy + d.L1 / 2) / 2
      dvy = û * Δth * dgy / 2

   else

      throw(ArgumentError("jinvmap for squircle expects region 1–4; got r=$r"))

   end

   detJ = dux * dvy - dvx * duy
   invdet = 1.0 / detJ

   # inverse of [dux dvx; duy dvy] is (1/detJ)[dvy -dvx; -duy dux]
   J11 = invdet * dvy
   J12 = -invdet * dvx
   J21 = -invdet * duy
   J22 = invdet * dux

   return J11, J12, J21, J22

end

"""
  mapinv(tbl::FTable, d::squircle, u::Float64, v::Float64, k::Int)
      -> Tuple{Float64,Float64}

Inverse mapping: given a physical point `(u,v)` and patch `k`, return
the reference coordinates `(t, s)` such that `mapxy(d, t, s, k) = (u,v)`.

Strategy:
1. If `k` is the center rectangle, use the affine inverse.
2. Otherwise, try 2D Newton on the reference region.
3. If Newton fails, use the 1D root fallback based on `FTable`.
"""
@inline function Xx(s::Float64, d::squircle, r::Int)
   if r == 1
      return d.A + d.L2 / 2
   elseif r == 2
      return d.A + d.L2 * (1 - 2s) / 2
   elseif r == 3
      return d.A - d.L2 / 2
   elseif r == 4
      return d.A - d.L2 * (1 - 2s) / 2
   else
      throw(ArgumentError("Xx expects region 1–4; got r=$r"))
   end
end

@inline function Xy(s::Float64, d::squircle, r::Int)
   if r == 1
      return d.B + d.L1 * (2s - 1) / 2
   elseif r == 2
      return d.B + d.L1 / 2
   elseif r == 3
      return d.B - d.L1 * (2s - 1) / 2
   elseif r == 4
      return d.B - d.L1 / 2
   else
      throw(ArgumentError("Xy expects region 1–4; got r=$r"))
   end
end

@inline function Yx(s::Float64, d::squircle, r::Int)

   th = if r == 1

      s * 2 * d.tht1 - d.tht1

   elseif r == 2

      s * (π - 2 * d.tht1) + d.tht1

   elseif r == 3

      s * 2 * d.tht1 + π - d.tht1

   elseif r == 4

      s * (π - 2 * d.tht1) + π + d.tht1

   else

      throw(ArgumentError("Yx expects region 1–4; got r=$r"))

   end

   st, ct = sincos(th)

   hP = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

   return d.A + hP * ct

end

@inline function Yy(s::Float64, d::squircle, r::Int)

   th = if r == 1

      s * 2 * d.tht1 - d.tht1

   elseif r == 2

      s * (π - 2 * d.tht1) + d.tht1

   elseif r == 3

      s * 2 * d.tht1 + π - d.tht1

   elseif r == 4

      s * (π - 2 * d.tht1) + π + d.tht1

   else

      throw(ArgumentError("Yy expects region 1–4; got r=$r"))

   end

   st, ct = sincos(th)

   hP = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

   return d.B + hP * st

end

function fill_FTable!(tbl::FTable, d::squircle, r::Int)

   vmin, vmax, N = tbl.vmin, tbl.vmax, tbl.N
   h = (vmax - vmin) / (N - 1)

   P1, P2, P3 = tbl.P1, tbl.P2, tbl.P3

   if r == 1

      Xx_const = d.A + d.L2 / 2

      @inbounds for i in 1:N

         v̂ = vmin + (i - 1) * h

         Xxv = Xx_const
         Xyv = d.B + d.L1 * (2v̂ - 1) / 2

         th = v̂ * 2 * d.tht1 - d.tht1
         st, ct = sincos(th)
         hP = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

         Yxv = d.A + hP * ct
         Yyv = d.B + hP * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv

      end

   elseif r == 2

      Xy_const = d.B + d.L1 / 2

      @inbounds for i in 1:N

         v̂ = vmin + (i - 1) * h

         Xxv = d.A + d.L2 * (1 - 2v̂) / 2
         Xyv = Xy_const

         th = v̂ * (π - 2 * d.tht1) + d.tht1
         st, ct = sincos(th)
         hP = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

         Yxv = d.A + hP * ct
         Yyv = d.B + hP * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv

      end

   elseif r == 3

      Xx_const = d.A - d.L2 / 2

      @inbounds for i in 1:N

         v̂ = vmin + (i - 1) * h

         Xxv = Xx_const
         Xyv = d.B - d.L1 * (2v̂ - 1) / 2

         th = v̂ * 2 * d.tht1 + π - d.tht1
         st, ct = sincos(th)
         hP = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

         Yxv = d.A + hP * ct
         Yyv = d.B + hP * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv

      end

   elseif r == 4

      Xy_const = d.B - d.L1 / 2

      @inbounds for i in 1:N

         v̂ = vmin + (i - 1) * h

         Xxv = d.A - d.L2 * (1 - 2v̂) / 2
         Xyv = Xy_const

         th = v̂ * (π - 2 * d.tht1) + π + d.tht1
         st, ct = sincos(th)
         hP = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

         Yxv = d.A + hP * ct
         Yyv = d.B + hP * st

         P1[i] = Yyv - Xyv
         P2[i] = Yxv - Xxv
         P3[i] = Xyv * Yxv - Xxv * Yyv

      end

   else

      throw(ArgumentError("fill_FTable! expects region 1–4; got r=$r"))

   end

   tbl.reg = r

   return tbl

end

@inline function f1I(t::Float64, s::Float64,
   d::squircle, u::Float64, v::Float64, r::Int)

   ŝ = (s + 1) / 2
   t̂ = (t + 1) / 2

   Xxv = Xx(ŝ, d, r)
   Yxv = Yx(ŝ, d, r)

   return (1 - t̂) * Xxv + t̂ * Yxv - u

end

@inline function f2I(t::Float64, s::Float64,
   d::squircle, u::Float64, v::Float64, r::Int)

   ŝ = (s + 1) / 2
   t̂ = (t + 1) / 2

   Xyv = Xy(ŝ, d, r)
   Yyv = Yy(ŝ, d, r)

   return (1 - t̂) * Xyv + t̂ * Yyv - v

end

@inline function JinvI(t::Float64, s::Float64,
   d::squircle, u::Float64, v::Float64, r::Int)

   return jinvmap(d, t, s, r)

end

@inline function f_cont(v̂::Float64,
   d::squircle, u::Float64, v::Float64, r::Int)

   if r == 1

      Xxv = d.A + d.L2 / 2
      Xyv = d.B + d.L1 * (2v̂ - 1) / 2
      th = v̂ * 2 * d.tht1 - d.tht1

   elseif r == 2

      Xxv = d.A + d.L2 * (1 - 2v̂) / 2
      Xyv = d.B + d.L1 / 2
      th = v̂ * (π - 2 * d.tht1) + d.tht1

   elseif r == 3

      Xxv = d.A - d.L2 / 2
      Xyv = d.B - d.L1 * (2v̂ - 1) / 2
      th = v̂ * 2 * d.tht1 + π - d.tht1

   elseif r == 4

      Xxv = d.A - d.L2 * (1 - 2v̂) / 2
      Xyv = d.B - d.L1 / 2
      th = v̂ * (π - 2 * d.tht1) + π + d.tht1

   else

      throw(ArgumentError("f_cont expects region 1–4; got r=$r"))

   end

   st, ct = sincos(th)

   hP = ((abs(ct) / d.R1)^d.P + (abs(st) / d.R2)^d.P)^(-1 / d.P)

   Yxv = d.A + hP * ct
   Yyv = d.B + hP * st

   dx = Yxv - Xxv
   dy = Yyv - Xyv

   term1 = muladd(u, dy, -v * dx)
   term2 = muladd(Xyv, Yxv, -Xxv * Yyv)

   return term1 + term2

end

function mapinv(tbl::FTable, d::squircle, u::Float64,
   v::Float64, k::Int)::Tuple{Float64, Float64}

   p = d.pths[k]
   r = p.reg

   # ----- curved patches reg = 1..4 -----
   # Stage 1: 2D Newton via subroutines.
   tN, sN = newtonR2D(f1I, f2I, JinvI,
      0.0, 0.0, 4, d, u, v, r; tol = 1e-15)

   if tN !== :max

      return xi_inv((1 + tN) / 2, p.ck0, p.ck1),
             xi_inv((1 + sN) / 2, p.tk0, p.tk1)

   end

   # Stage 2: BIS method inversion.

   if tbl.reg != r
      fill_FTable!(tbl, d, r)
   end

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

function mapinv(d::squircle, u::Float64,
   v::Float64, k::Int)::Tuple{Float64,Float64}

   p = d.pths[k]

   if d.L2 >= d.L1

      Z1 = xi_inv((u - d.A + d.L2 / 2) / d.L2, p.tk0, p.tk1)
      Z2 = xi_inv((v - d.B + d.L1 / 2) / d.L1, p.ck0, p.ck1)

   else

      Z1 = xi_inv((u - d.A + d.L2 / 2) / d.L2, p.ck0, p.ck1)
      Z2 = xi_inv((v - d.B + d.L1 / 2) / d.L1, p.tk0, p.tk1)

   end

   return Z1, Z2
end
"""
    ptconv(d::squircle,  t1::Float64, t2::Float64, idx::Int, ptdest::String)

Convert a parametric point between region ↔ patch coordinates.

Inputs:
- `t1` : 1st cooridnate of point t
- `t2` : 2nd cooridnate of point t
- `idx`: region index (if `"to_pth"`) or patch index (if `"to_reg"`)
- `ptdest`: `"to_pth"` or `"to_reg"`

Returns:
- `Tuple(Float64, Float64, out_idx::Int)`
  where `out_idx` is patch index if `"to_pth"`, region index if `"to_reg"`.
"""
function ptconv(d::squircle, t1::Float64, t2::Float64, idx::Int, ptdest::String)
  #@assert length(t) == 2 "t must be length-2 vector [t1,t2]"

  if ptdest == "to_pth"
    r = idx

    for k in 1:d.Npat
      p = d.pths[k]
      p.reg == r || continue

      if r != 5
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
    #@assert 1 ≤ k ≤ d.Npat "patch index out of range"
    p = d.pths[k]
    r = p.reg

    if r != 5
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
    mapinv2(d::squircle, t1::Float64, t2::Float64, k2::Int, k::Int) -> Tuple{Float64,Float64}

Given a point `(u,v) = τₖ₂(t1,t2)` on patch `k2`, return its reference
coordinates on patch `k` in `[-1,1]^2`.

This differs from `mapinv` because we already know `(u,v)` comes from `(t1,t2)` on `k2`.
"""
function mapinv2(d::squircle, t1::Float64, t2::Float64, k2::Int, k::Int)::Tuple{Float64,Float64}
    p = d.pths[k]  # Fields reg, ck0, ck1, tk0, tk1

    # Convert the given (t1,t2,k2) to the "regional" normalized coords of the *point*
    # Expecting something that returns two values in [-1,1]:
    #   tr1 ≈ u_ref, tr2 ≈ v_ref  (for the *global/regional* parameterization)
    # Adjust this call if your ptconv signature differs.
    tr1, tr2, _ = ptconv(d, t1, t2, k2, "to_reg")   

    # Map from [-1,1] -> [0,1]
    û = (tr1 + 1) / 2
    v̂ = (tr2 + 1) / 2

    # Axis swap for the center rectangle if L2 ≥ L1
    if d.L2 >= d.L1 && p.reg == 5
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
- If patch `k` is in region 5 → returns ones with same shape as `t`.
- Otherwise uses: `d = 1 - ck1 + (ck1 - ck0) * t / 2`, then
raises elementwise to the power `s-1` if `s ≥ 0.5`, else to `s`.

Notes:
- `t` is expected to be **1 - t_actual**
- `t` may be a scalar `Float64` or any `StridedArray`; the return has the same shape.
"""

function dfunc(d::squircle, k::Int, t::Float64, s::Float64)::Float64

  p = d.pths[k]

  if p.reg == 5
    return 1.0
  end

  hc = p.ck1 - p.ck0
  val = (1.0 - p.ck1) + hc * t / 2
  exp = s ≥ 0.5 ? (s - 1) : s
  return val^exp

end

function dfunc!(out::StridedArray{Float64}, d::squircle, k::Int, t::StridedArray{Float64}, s::Float64)

  p = d.pths[k]

  if p.reg == 5
    fill!(out, 1.0)
  else
    exp = s ≥ 0.5 ? (s - 1) : s
    αc = (p.ck1 - p.ck0) / 2
    βc = 1.0 - p.ck1
    @inbounds for i in eachindex(t)
      out[i] = (muladd(αc,t[i],βc))^exp
    end
  end

  return nothing

end

function dfunc!(out::StridedArray{Float64}, d::squircle, k::Int, t::StridedArray{Float64})

  p = d.pths[k]

  if p.reg == 5
    fill!(out, 1.0)
  else
    αc = (p.ck1 - p.ck0) / 2
    βc = 1.0 - p.ck1
    @inbounds for i in eachindex(t)
      out[i] = muladd(αc,t[i],βc)
    end
  end

  return nothing

end

"""
    gamderhigher(d::squircle, th::Float64)

Return centered squircle boundary derivatives through order 4.

Returns gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy, d4gx, d4gy

where

    gx(th) = hP(th) * cos(th)
    gy(th) = hP(th) * sin(th)

and derivatives are with respect to `th`.

Supports P == 2, P == 4 and P > 4. 
"""
function gamderhigher(d::squircle, th::Float64)

   P = d.P

   st, ct = sincos(th)

   if P == 2.0
      A1 = 1.0 / d.R1^P
      B1 = 1.0 / d.R2^P
      Δ = B1 - A1

      gg = A1 * ct^2 + B1 * st^2
      dgg = 2.0 * st * ct * Δ
      d2gg = 2.0 * (ct^2 - st^2) * Δ
      d3gg = -8.0 * st * ct * Δ
      d4gg = -8.0 * (ct^2 - st^2) * Δ

   elseif P == 4.0
      A1 = 1.0 / d.R1^P
      B1 = 1.0 / d.R2^P

      st2 = st * st
      ct2 = ct * ct
      st4 = st2 * st2
      ct4 = ct2 * ct2

      gg = A1 * ct4 + B1 * st4

      dgg = -4.0 * st * ct * (A1 * ct2 - B1 * st2)

      d2gg = -4.0 * (A1 * ct4 + B1 * st4 -
         3.0 * (A1 + B1) * st2 * ct2)

      d3gg = 8.0 * st * ct * (
                A1 * (5.0 * ct2 - 3.0 * st2) +
                B1 * (3.0 * ct2 - 5.0 * st2))

      d4gg = 8.0 * (
         A1 * (3.0 * st4 - 24.0 * st2 * ct2 + 5.0 * ct4) +
         B1 * (5.0 * st4 - 24.0 * st2 * ct2 + 3.0 * ct4))

   else
      ct1 = abs(ct)^(P - 4) / d.R1^P
      ct2 = ct * ct * ct1
      ct3 = ct * ct * ct2

      st1 = abs(st)^(P - 4) / d.R2^P
      st2 = st * st * st1
      st3 = st * st * st2

      gg = ct3 + st3

      dgg = -P * st * ct * (ct2 - st2)

      d2gg = P * (P - 1) * (ct2 + st2) - P^2 * gg

      d3gg = -P * (P - 1) * (P - 2) * st * ct * (ct1 - st1) -
             P^2 * dgg

      d4gg = P * (P - 1) * (P - 2) * (P - 3) * (ct1 + st1) -
             P * (P - 1) * (P - 2)^2 * (ct2 + st2) -
             P^2 * d2gg

   end

   h = gg^(-1 / P)

   dh = -(h * dgg) / (P * gg)

   d2h = -((P + 1) * dh * dgg + h * d2gg) / (P * gg)

   d3h = -((2P + 1) * d2h * dgg +
           (P + 2) * dh * d2gg +
           h * d3gg) / (P * gg)

   d4h = -((3P + 1) * d3h * dgg +
           (3P + 3) * d2h * d2gg +
           (P + 3) * dh * d3gg +
           h * d4gg) / (P * gg)

   gx = h * ct
   gy = h * st

   dgx = dh * ct - h * st
   dgy = dh * st + h * ct

   d2gx = d2h * ct - 2 * dh * st - h * ct
   d2gy = d2h * st + 2 * dh * ct - h * st

   d3gx = d3h * ct - 3 * d2h * st - 3 * dh * ct + h * st
   d3gy = d3h * st + 3 * d2h * ct - 3 * dh * st - h * ct

   d4gx = d4h * ct - 4 * d3h * st - 6 * d2h * ct + 4 * dh * st + h * ct
   d4gy = d4h * st + 4 * d3h * ct - 6 * d2h * st - 4 * dh * ct + h * st

   return gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy, d4gx, d4gy

end

"""
  diff_map!(out::Matrix{Float64}, Zx::Matrix{Float64}, Zy::Matrix{Float64},
            DJ::StridedArray{Float64}, d::squircle,
            u::Float64, v::Float64,
            u2::Matrix{Float64}, v2::Matrix{Float64},
            du::AbstractVector, dv::AbstractVector, k::Int;
            tol = 1e-4)

Fill `out` with ‖τ(u,v) - τ(u₂,v₂)‖ on patch `k`.

No allocations. Uses `mapxy_Dmap!` and Taylor correction near `(u,v)`.
"""
function diff_map!(out::Matrix{Float64},
   Zx::Matrix{Float64}, Zy::Matrix{Float64}, DJ::StridedArray{Float64},
   d::squircle, u::Float64, v::Float64,
   u2::Matrix{Float64}, v2::Matrix{Float64},
   du::AbstractVector, dv::AbstractVector, k::Int;
   tol::Float64 = 1e-4)

   nd_u = size(out, 1)
   nd_v = size(out, 2)

   p = d.pths[k]

   αc = (p.ck1 - p.ck0) / 2
   αt = (p.tk1 - p.tk0) / 2

   mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)

   if p.reg == 5

      if d.L2 >= d.L1
         cDx, cDy = d.L2 * αt, d.L1 * αc
      else
         cDx, cDy = d.L2 * αc, d.L1 * αt
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

   if p.reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

      Xx = d.L2 / 2
      Xy = d.L1 * (2vhat - 1) / 2
      dXx = 0.0
      dXy = d.L1

   elseif p.reg == 2

      Δth = π - 2 * d.tht1
      th0 = d.tht1

      Xx = d.L2 * (1 - 2vhat) / 2
      Xy = d.L1 / 2
      dXx = -d.L2
      dXy = 0.0

   elseif p.reg == 3

      Δth = 2 * d.tht1
      th0 = π - d.tht1

      Xx = -d.L2 / 2
      Xy = -d.L1 * (2vhat - 1) / 2
      dXx = 0.0
      dXy = -d.L1

   elseif p.reg == 4

      Δth = π - 2 * d.tht1
      th0 = π + d.tht1

      Xx = -d.L2 * (1 - 2vhat) / 2
      Xy = -d.L1 / 2
      dXx = d.L2
      dXy = 0.0

   end

   th = muladd(Δth, vhat, th0)

   gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy, d4gx, d4gy =
      gamderhigher(d, th)

   q = αt * Δth
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
    gamderhigher6(d::squircle, th::Float64)

Return centered squircle boundary derivatives through order 6.

Returns gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy, d4gx, d4gy, 
d5gx, d5gy, d6gx, d6gy

where gx(th) = h(th) * cos(th), gy(th) = h(th) * sin(th).

and derivatives are with respect to `th`.
"""
function gamderhigher6(d::squircle, th::Float64)

   st, ct = sincos(th)

   h, h1, h2, h3, h4, h5, h6 = hderhigher(d, th)

   gx = h * ct
   gy = h * st

   dgx = h1 * ct - h * st
   dgy = h1 * st + h * ct

   d2gx = h2 * ct - 2.0 * h1 * st - h * ct
   d2gy = h2 * st + 2.0 * h1 * ct - h * st

   d3gx = h3 * ct - 3.0 * h2 * st - 3.0 * h1 * ct + h * st
   d3gy = h3 * st + 3.0 * h2 * ct - 3.0 * h1 * st - h * ct

   d4gx = h4 * ct - 4.0 * h3 * st - 6.0 * h2 * ct +
          4.0 * h1 * st + h * ct

   d4gy = h4 * st + 4.0 * h3 * ct - 6.0 * h2 * st -
          4.0 * h1 * ct + h * st

   d5gx = h5 * ct - 5.0 * h4 * st - 10.0 * h3 * ct +
          10.0 * h2 * st + 5.0 * h1 * ct - h * st

   d5gy = h5 * st + 5.0 * h4 * ct - 10.0 * h3 * st -
          10.0 * h2 * ct + 5.0 * h1 * st + h * ct

   d6gx = h6 * ct - 6.0 * h5 * st - 15.0 * h4 * ct +
          20.0 * h3 * st + 15.0 * h2 * ct -
          6.0 * h1 * st - h * ct

   d6gy = h6 * st + 6.0 * h5 * ct - 15.0 * h4 * st -
          20.0 * h3 * ct + 15.0 * h2 * st +
          6.0 * h1 * ct - h * st

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
             DJ::StridedArray{Float64}, d::squircle,
             u::Float64, v::Float64,
             u2::Matrix{Float64}, v2::Matrix{Float64}, r::Matrix{Float64},
             du::AbstractVector, dv::AbstractVector, k::Int;
             tol = 1e-3)

Compute ‖(τ(u,v) - τ(u₂,v₂)) / r‖ for the `k`-th patch, where

    u₂ = u - r .* du
    v₂ = v - r .* dv

No allocations. Uses `mapxy_Dmap!` and Taylor correction near `(u,v)`.
"""
function diff_rmap!(out::Matrix{Float64},
   Zx::Matrix{Float64}, Zy::Matrix{Float64}, DJ::StridedArray{Float64},
   d::squircle, u::Float64, v::Float64,
   u2::Matrix{Float64}, v2::Matrix{Float64}, r::Matrix{Float64},
   du::AbstractVector, dv::AbstractVector, k::Int;
   tol::Float64 = 1e-3)

   nt = size(out, 1)
   nr = size(out, 2)

   p = d.pths[k]

   αc = (p.ck1 - p.ck0) / 2
   αt = (p.tk1 - p.tk0) / 2

   if p.reg == 5

      if d.L2 >= d.L1
         cDx, cDy = d.L2 * αt, d.L1 * αc
      else
         cDx, cDy = d.L2 * αc, d.L1 * αt
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

   if p.reg == 1

      Δth = 2 * d.tht1
      th0 = -d.tht1

      Xx = d.L2 / 2
      Xy = d.L1 * (2vhat - 1) / 2
      dXx = 0.0
      dXy = d.L1

   elseif p.reg == 2

      Δth = π - 2 * d.tht1
      th0 = d.tht1

      Xx = d.L2 * (1 - 2vhat) / 2
      Xy = d.L1 / 2
      dXx = -d.L2
      dXy = 0.0

   elseif p.reg == 3

      Δth = 2 * d.tht1
      th0 = π - d.tht1

      Xx = -d.L2 / 2
      Xy = -d.L1 * (2vhat - 1) / 2
      dXx = 0.0
      dXy = -d.L1

   elseif p.reg == 4

      Δth = π - 2 * d.tht1
      th0 = π + d.tht1

      Xx = -d.L2 * (1 - 2vhat) / 2
      Xy = -d.L1 / 2
      dXx = d.L2
      dXy = 0.0

   else

      throw(ArgumentError("diff_rmap! expects region 1–5; got reg=$(p.reg)"))

   end

   th = muladd(Δth, vhat, th0)

   gx, gy, dgx, dgy, d2gx, d2gy, d3gx, d3gy, 
   d4gx, d4gy, d5gx, d5gy, d6gx, d6gy = gamderhigher6(d, th)

   q = αt * Δth

   q2 = q * q
   q3 = q2 * q
   q4 = q2 * q2
   q5 = q4 * q
   q6 = q3 * q3

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
    Dwall(d::squircle, u::Float64, v::Float64, k::Int) -> dvx, dvy

Derivative of the wall curve with respect to `v` at the side `u = ±1`
for the `k`-th patch.

Returns a Tuple `(dvx, dvy)`.
"""
function Dwall(d::squircle, u::Float64, v::Float64, k::Int)::Tuple{Float64, Float64}

   p = d.pths[k]

   hc = p.ck1 - p.ck0
   ht = p.tk1 - p.tk0

   # affine maps: [-1,1] → [ck0,ck1] and [tk0,tk1]
   v̂ = ht * (v + 1) / 2 + p.tk0
   û = hc * (u + 1) / 2 + p.ck0

   if p.reg == 1

      th = v̂ * 2 * d.tht1 - d.tht1
      Δth = 2 * d.tht1

      _, _, dgx, dgy = gamder(d, th)

      dvx = ht * û * Δth * dgx / 2
      dvy = (1 - û) * d.L1 * ht / 2 + ht * û * Δth * dgy / 2

   elseif p.reg == 2

      th = v̂ * (π - 2 * d.tht1) + d.tht1
      Δth = π - 2 * d.tht1

      _, _, dgx, dgy = gamder(d, th)

      dvx = -(1 - û) * d.L2 * ht / 2 + ht * û * Δth * dgx / 2
      dvy = ht * û * Δth * dgy / 2

   elseif p.reg == 3

      th = v̂ * 2 * d.tht1 + π - d.tht1
      Δth = 2 * d.tht1

      _, _, dgx, dgy = gamder(d, th)

      dvx = ht * û * Δth * dgx / 2
      dvy = (û - 1) * d.L1 * ht / 2 + ht * û * Δth * dgy / 2

   elseif p.reg == 4

      th = v̂ * (π - 2 * d.tht1) + π + d.tht1
      Δth = π - 2 * d.tht1

      _, _, dgx, dgy = gamder(d, th)

      dvx = (1 - û) * d.L2 * ht / 2 + ht * û * Δth * dgx / 2
      dvy = ht * û * Δth * dgy / 2

   else

      throw(ArgumentError("Dwall expects region 1–4; got reg=$(p.reg)"))

   end

   return dvx, dvy

end

"""
   refine!(d::squircle, Nc::Int, Nt::Int, K::AbstractVector{<:Int})

Subdivide each patch in `K` into `Nc * Nt` subpatches, split in `c` and `t`.
Updates `d.pths`, `d.Npat`, recomputes `d.kd`, and rebuilds `d.Qpts` / `d.Qptsbd`.

Notes
- Patches are re-sorted lexicographically by `(reg, ck0, ck1, tk0, tk1)`.
- `Qpts` uses `boundquad!(V, d, k)`.
- `Qptsbd` uses `boundquadbd!(V, d, k)` for patches in `d.kd`.
"""
function refine!(d::squircle, Nc::Int, Nt::Int, K::Vector{Int})
   @assert Nc ≥ 1 && Nt ≥ 1 "Nc and Nt must be ≥ 1"
   @assert !isempty(K) "K (set of patches to refine) is empty"

   # create the subdivided children
   children = Patch[]
   sizehint!(children, length(K) * Nc * Nt)

   for k in K
      p = d.pths[k]
      cc = range(p.ck0, p.ck1; length=Nc + 1)
      tt = range(p.tk0, p.tk1; length=Nt + 1)

      for i in 1:Nc, j in 1:Nt
         push!(children, Patch(p.reg, cc[i], cc[i+1], tt[j], tt[j+1]))
      end
   end

   # merge and sort
   d.pths = vcat([d.pths[i] for i in 1:d.Npat if !(i in K)], children)
   sort!(d.pths, by=q -> (q.reg, q.ck0, q.ck1, q.tk0, q.tk1))

   # update counts
   d.Npat = length(d.pths)

   # recompute boundary-touching patch indices
   d.kd = [k for k in 1:d.Npat if d.pths[k].reg != 5 && d.pths[k].ck1 == 1.0]

   # Qpts: 8 × Npat
   d.Qpts = Matrix{Float64}(undef, 8, d.Npat)

   @inbounds for k in 1:d.Npat
      @views V = d.Qpts[:, k]
      boundquad!(V, d, k)
   end

   # Qptsbd: 8 × length(kd)
   d.Qptsbd = Matrix{Float64}(undef, 8, length(d.kd))

   @inbounds for (l, k) in enumerate(d.kd)
      @views V = d.Qptsbd[:, l]
      boundquadbd!(V, d, k)
   end

   return d
end

function Base.show(io::IO, d::squircle)

   println(io, "squircle with properties:")
   println(io, "  (A, B)   = (", d.A, ", ", d.B, ")")
   println(io, "  P        = ", d.P)
   println(io, "  tht1     = ", d.tht1)
   println(io, "  (R1, R2) = (", d.R1, ", ", d.R2, ")")
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

