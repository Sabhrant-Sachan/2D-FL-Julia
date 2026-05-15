"""
Mutable struct disc

A  — x-coordinate of disc center (Float64)

B  — y-coordinate of disc center (Float64)

R  — radius (> 0)  (Float64)

L₁ — length of rectangle's y-axis side (Float64)

L₂ — length of rectangle's x-axis side (Float64)

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
mutable struct disc <: abstractdomain

  A::Float64
  B::Float64
  R::Float64
  L1::Float64
  L2::Float64
  nh::Int
  kd::Vector{Int}
  Npat::Int
  pths::Vector{Patch}   # 5 * Npat
  Qpts::Matrix{Float64}   # 8 * Npat
  Qptsbd::Matrix{Float64} # 8 * Npat

  function disc(; b,
    a=nothing, A=nothing, B=nothing,
    R=nothing, L1=nothing, L2=nothing,
    ck=nothing, tk=nothing)
    #Construct an intsance for this structure

    """
    Required:
      - b :: Vector{Int} of length 5 (Col), positive entries
    
    Optional (all keywords):
      - a  :: Vector{Int} length 5 (defaults to ceil.(2*b/3))
      - A,B:: Center of the disc (defaults 0.0)
      - R  :: Float64 (default 1.0)
      - L1 :: Float64 (default 0.8*R)
      - L2 :: Float64 (default 0.8*R)
      - ck :: vector{vector{Float64}} (default equispaced per a)
      - tk :: vector{vector{Float64}} (default equispaced per b)
    """
    @assert isa(b, AbstractVector{<:Int}) && length(b) == 5

    #The disc has obviously no holes
    nh = 0

    a = isnothing(a) ? ceil.(Int, 2 .* b ./ 3) : collect(Int, a)

    #By default, a unit disc at origin
    A = something(A, 0)
    B = something(B, 0)
    R = something(R, 1.0)

    #By default, a sqaure of size length 0.8 units is inscribed
    L1 = something(L1, 0.8 * R)
    L2 = something(L2, 0.8 * R)

    #Sum up all the subpatches

    Npat = dot(a,b)

    # ck, tk as vector-of-vectors with *no padding*:
    #   ck[k] has length a[k] + 1
    #   tk[k] has length b[k] + 1
    ck = something(ck, [(0:a[k]) ./ a[k] for k in 1:5])
    tk = something(tk, [(0:b[k]) ./ b[k] for k in 1:5])

    # Preallocate pths and kd
    pths = Vector{Patch}(undef, Npat)

    #Starting with the assumption that all interioir 
    #patches touch the boundary
    kd = Vector{Int}(undef, Npat)   # upper bound length
    nkd = 0                         # actual count

    # Fill pths and kd in a single pass
    idx = 1
    @inbounds for k in 1:5
      akℓ = a[k]
      bkℓ = b[k]
      ckℓ = ck[k]
      tkℓ = tk[k]

      @inbounds for i in 1:akℓ
        ck1 = ckℓ[i]
        ck2 = ckℓ[i+1]
        @inbounds for j in 1:bkℓ
          tk1 = tkℓ[j]
          tk2 = tkℓ[j+1]

          p = Patch(k, ck1, ck2, tk1, tk2)
          pths[idx] = p

          if p.reg != 5 && p.ck1 == 1.0
            nkd += 1
            kd[nkd] = idx
          end

          idx += 1
        end
      end
    end

    # Shrink kd to actual size
    resize!(kd, nkd)

    # Allocate quadrature arrays (no need for zeros; we'll overwrite everything)
    Qpts = Matrix{Float64}(undef, 8, Npat)
    Qptsbd = Matrix{Float64}(undef, 8, nkd)

    d = new(A, B, R, L1, L2, nh, kd, Npat, pths, Qpts, Qptsbd)

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

function mapx(d::disc, u::Float64, v::Float64, k::Int)::Float64

  # p denotes information of kth patch
  p = d.pths[k]

  # Difference in c values
  hc = p.ck1 - p.ck0

  xi1 = hc * u / 2 + (p.ck1 + p.ck0) / 2

  # Difference in t values
  ht = p.tk1 - p.tk0

  xi2 = ht * v / 2 + (p.tk1 + p.tk0) / 2

  if p.reg != 5

    X = if p.reg == 1

      d.A + d.L2 / 2

    elseif p.reg == 3

      d.A - d.L2 / 2

    elseif p.reg == 2

      d.A + d.L2 * (1 - 2 * xi2) / 2

    else # p.reg == 4

      d.A - d.L2 * (1 - 2 * xi2) / 2

    end

    Y = d.A + d.R * cospi(xi2 / 2 + (2 * p.reg - 3) / 4)

    return (1 - xi1) * X + xi1 * Y

  else

    xi1 = d.L2 >= d.L1 ? ht * u / 2 + (p.tk1 + p.tk0) / 2 : xi1

    return d.A - d.L2 / 2 + d.L2 * xi1

  end

end

function mapy(d::disc, u::Float64, v::Float64, k::Int)::Float64

  # p denotes information of kth patch
  p = d.pths[k]

  # Difference in c values
  hc = p.ck1 - p.ck0

  xi1 = hc * u / 2 + (p.ck1 + p.ck0) / 2

  # Difference in t values
  ht = p.tk1 - p.tk0

  xi2 = ht * v / 2 + (p.tk1 + p.tk0) / 2

  if p.reg != 5

    X = if p.reg == 2

      d.B + d.L1 / 2

    elseif p.reg == 4

      d.B - d.L1 / 2

    elseif p.reg == 1

      d.B + d.L1 * (2 * xi2 - 1) / 2

    else # p.reg == 3

      d.B - d.L1 * (2 * xi2 - 1) / 2

    end

    Y = d.B + d.R * sinpi(xi2 / 2 + (2 * p.reg - 3) / 4)

    return (1 - xi1) * X + xi1 * Y

  else

    xi2 = d.L2 >= d.L1 ? hc * v / 2 + (p.ck1 + p.ck0) / 2 : xi2

    return d.B - d.L1 / 2 + d.L1 * xi2

  end

end

function mapxy(d::disc, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64}

  return mapx(d, u, v, k), mapy(d, u, v, k)

end

function mapxy!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
  d::disc, u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)
  #@assert size(Zx) == size(Zy) == size(u) == size(v)
  p = d.pths[k]
  hc = p.ck1 - p.ck0
  ht = p.tk1 - p.tk0

  if p.reg != 5
    # xi1 = map u from [-1,1] → [ck0,ck1]; xi2 = map v → [tk0,tk1]
    αu = hc / 2
    αv = ht / 2
    βu = p.ck0 + αu
    βv = p.tk0 + αv

    # Corner line (Xx,Xy) for this region
    # Iterate once over a single index range I that is 
    # valid for all four arrays.
    @inbounds for I in eachindex(u, v, Zx, Zy)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      # fixed edge point for given xi2
      Xx = if p.reg == 1
        d.A + d.L2 / 2
      elseif p.reg == 2
        d.A + d.L2 * (1 - 2xi2) / 2
      elseif p.reg == 3
        d.A - d.L2 / 2
      else # p.reg == 4
        d.A - d.L2 * (1 - 2xi2) / 2
      end

      Xy = if p.reg == 1
        d.B + d.L1 * (2xi2 - 1) / 2
      elseif p.reg == 2
        d.B + d.L1 / 2
      elseif p.reg == 3
        d.B - d.L1 * (2xi2 - 1) / 2
      else # p.reg == 4
        d.B - d.L1 / 2
      end

      # circular arc point for given xi2
      t = π *(xi2 / 2 + (2 * p.reg - 3) / 4)

      s, c = sincos(t)

      Yx = muladd(d.R, c, d.A) #d.A + d.R * c
      Yy = muladd(d.R, s, d.B) #d.B + d.R * s

      # blend: (1 - xi1) * X + xi1 * Y
      Zx[I] = muladd(xi1, (Yx - Xx), Xx)
      Zy[I] = muladd(xi1, (Yy - Xy), Xy)
    end

  else
    # rectangular patch (reg==5)
    # if L2 >= L1, swap affine maps per your code path
    if d.L2 >= d.L1
      αu = ht / 2
      αv = hc / 2
      βu = p.tk0 + αu
      βv = p.ck0 + αv
      @inbounds  for I in eachindex(u, v, Zx, Zy)
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
      @inbounds  for I in eachindex(u, v, Zx, Zy)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)
        Zx[I] = d.A - d.L2 / 2 + d.L2 * xi1
        Zy[I] = d.B - d.L1 / 2 + d.L1 * xi2
      end

    end
  end

  return nothing

end

function draw(d::disc, flag=nothing; L::Int=32, show::Bool=true)

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
    gamx(d::disc, t::Float64, k::Int)
    gamx(d::disc, t::StridedArray{Float64}, k::Int)

First (x) coordinate of the boundary parametrization γ(t).

- `gamx(d, t, k)` computes the **right boundary** of patch `k`, where `t ∈ [-1,1]`.

`t` may be a scalar or an array; the return has the same shape.
"""
function gamx(d::disc, t::Float64, k::Int)::Float64
    p = d.pths[k]
    # Map t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
    xi2 = (p.tk1 - p.tk0) * t / 2 + (p.tk0 + p.tk1) / 2
    phase = xi2 / 2 + (2*p.reg - 3) / 4
    return d.A + d.R * cospi(phase)
end

function gamx!(out::StridedArray{Float64}, d::disc, t::StridedArray{Float64}, k::Int)
  
  p  = d.pths[k]
  αt = (p.tk1 - p.tk0) / 2
  βt = p.tk0 + αt
  c  = pi * (2 * p.reg - 3) / 4
  @inbounds  for I in eachindex(t)
    # Map t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
    τ = pi * muladd(αt, t[I], βt) / 2 + c
    out[I] = muladd(d.R, cos(τ), d.A) 
  end

  return nothing
end

"""
    gamy(d::disc, t::Float64, k::Int)
    gamy(d::disc, t::StridedArray{Float64}, k::Int)

Second (y) coordinate of the boundary parametrization γ(t).

- `gamy(d, t, k)` computes the right boundary of patch `k`, where `t ∈ [-1,1]`.

`t` can be a scalar `Float64` or any `StridedArray{Float64}`; the output matches the input shape.
"""
function gamy(d::disc, t::Float64, k::Int)::Float64
    p = d.pths[k]
    # Map t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
    xi2   = (p.tk1 - p.tk0) * t / 2 + (p.tk0 + p.tk1) / 2
    t = xi2 / 2 + (2*p.reg - 3) / 4
    return d.B + d.R * sinpi(t)
end

function gamy!(out::StridedArray{Float64}, d::disc, t::StridedArray{Float64}, k::Int)

  # Map t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
  p = d.pths[k]
  αt = (p.tk1 - p.tk0) / 2
  βt = p.tk0 + αt
  c  = pi * (2 * p.reg - 3) / 4
  @inbounds  for I in eachindex(t)
    τ = pi * muladd(αt, t[I], βt) / 2 + c
    out[I] = muladd(d.R, sin(τ), d.B)
  end

  return nothing
end

function gam(d::disc, t::Float64, k::Int)::Tuple{Float64,Float64}
 
  p  = d.pths[k]
  αt = (p.tk1 - p.tk0) / 2
  βt = p.tk0 + αt
  c  = pi * (2 * p.reg - 3) / 4
  τ = pi * muladd(αt, t, βt) / 2 + c
  outx = muladd(d.R, cos(τ), d.A)
  outy = muladd(d.R, sin(τ), d.B)

  return outx, outy
end

function gam!(out::Vector{Float64}, d::disc, t::Float64, k::Int)
 
  p  = d.pths[k]
  αt = (p.tk1 - p.tk0) / 2
  βt = p.tk0 + αt
  c  = pi * (2 * p.reg - 3) / 4
  τ = pi * muladd(αt, t, βt) / 2 + c
  out[1] = muladd(d.R, cos(τ), d.A)
  out[2] = muladd(d.R, sin(τ), d.B)

  return nothing
end

function gam!(out::Matrix{Float64}, d::disc, t::Vector{Float64}, k::Int)
  #Out is a Matrix of same columns as length of t and 2 rows
  p  = d.pths[k]
  αt = (p.tk1 - p.tk0) / 2
  βt = p.tk0 + αt
  c  = π * (2 * p.reg - 3) / 4

  @inbounds  for I in eachindex(t)
    τ = π * muladd(αt, t[I], βt) / 2 + c
    out[1, I] = muladd(d.R, cos(τ), d.A) #x coord
    out[2, I] = muladd(d.R, sin(τ), d.B) #y coord
  end

  return nothing
end

function drawbd(d::disc, flag = true; L::Int = 33, show::Bool = true)

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
    dgamx(d::disc, t::Float64, k::Int) -> Float64
    dgamx(d::disc, t::StridedArray{Float64}, k::Int) -> StridedArray{Float64}

Derivative of the first coordinate of the boundary parametrization.

- With `k`: derivative of the **right boundary** of patch `k`
  (i.e., `∂/∂t mapx(d, 1, t, k)` for `t ∈ [-1,1]`).
  Uses  `t = ξ₂/2 + (2*reg - 3)/4`, where
  `ξ₂ = (tk1 - tk0)*t/2 + (tk0 + tk1)/2`. Result is
  `-(π*ht*R/4) * sinpi(t)` with `ht = tk1 - tk0`.
"""
function dgamx(d::disc, t::Float64, k::Int)::Float64
    p  = d.pths[k]
    ht = p.tk1 - p.tk0
    t = ht * t / 2 + (p.tk0 + p.tk1) / 2
    t = t / 2 + (2*p.reg - 3) / 4
    return -(π * ht * d.R / 4) * sinpi(t)
end

function dgamx!(out::StridedArray{Float64}, d::disc, t::StridedArray{Float64}, k::Int)
  
  p  = d.pths[k]
  αt = (p.tk1 - p.tk0) / 2
  βt = p.tk0 + αt
  c₁ =  π * (2 * p.reg - 3) / 4
  c₂ = -π * αt * d.R / 2
  @inbounds  for I in eachindex(t)
    # Map t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
    τ = pi * muladd(αt, t[I], βt) / 2 + c₁
    out[I] = c₂ * sin(τ)
  end

  return nothing
end

"""
    dgamy(d::disc, t::Float64, k::Int) -> Float64
    dgamy(d::disc, t::StridedArray{Float64}, k::Int) -> StridedArray{Float64}

Derivative of the second coordinate of the boundary parametrization γ(t).

- With `k`: derivative of the **right boundary** of patch `k`
  (`∂/∂t mapy(d, 1, t, k)` for `t ∈ [-1,1]`).
  If `ht = tk1 - tk0` and `t` is re-used for the phase
  `t = ( (ht*t/2 + (tk0+tk1)/2) / 2 ) + (2*reg - 3)/4`,
  then the result is `(π*ht*R/4) * cospi(t)`.
"""
function dgamy(d::disc, t::Float64, k::Int)::Float64
    p  = d.pths[k]
    ht = p.tk1 - p.tk0
    t = ht*(t+1)/2 + p.tk0
    t = t / 2 + (2*p.reg - 3) / 4
    return (π * ht * d.R / 4) * cospi(t)
end

function dgamy!(out::StridedArray{Float64}, d::disc, t::StridedArray{Float64}, k::Int)
  
  p  = d.pths[k]
  αt = (p.tk1 - p.tk0) / 2
  βt = p.tk0 + αt
  c₁ =  π * (2 * p.reg - 3) / 4
  c₂ =  π * αt * d.R / 2
  @inbounds  for I in eachindex(t)
    τ = pi * muladd(αt, t[I], βt) / 2 + c₁
    out[I] = c₂ * cos(τ)
  end

  return nothing
end

function dgam(d::disc, t::Float64, k::Int)::Tuple{Float64,Float64}

  p = d.pths[k]
  αt = (p.tk1 - p.tk0) / 2
  βt = p.tk0 + αt
  c₁ = π * (2 * p.reg - 3) / 4
  c₂ = π * αt * d.R / 2
  τ = π * muladd(αt, t, βt) / 2 + c₁
  s, c = sincos(τ)
  outx = -c₂ * s
  outy = c₂ * c

  return outx, outy
end

function dgam!(out::Vector{Float64}, d::disc, t::Float64, k::Int)

  p = d.pths[k]
  αt = (p.tk1 - p.tk0) / 2
  βt = p.tk0 + αt
  c₁ = π * (2 * p.reg - 3) / 4
  c₂ = π * αt * d.R / 2
  τ = π * muladd(αt, t, βt) / 2 + c₁
  s, c = sincos(τ)
  out[1] = -c₂ * s
  out[2] = c₂ * c

  return nothing
end

function dgam!(out::Matrix{Float64}, d::disc, t::Vector{Float64}, k::Int)

  p = d.pths[k]
  αt = (p.tk1 - p.tk0) / 2
  βt = p.tk0 + αt
  c₁ = π * (2 * p.reg - 3) / 4
  c₂ = π * αt * d.R / 2
  @inbounds  for I in eachindex(t)
    τ = pi * muladd(αt, t[I], βt) / 2 + c₁
    s, c = sincos(τ)
    out[1, I] = -c₂ * s
    out[2, I] = c₂ * c
  end

  return nothing
end

function gamp(d::disc, t::Float64, k::Int)::Tuple{Float64,Float64}
    return dgamy(d, t, k), -dgamx(d, t, k)
end

function gamp!(out::Matrix{Float64}, d::disc, t::Vector{Float64}, k::Int)

  p = d.pths[k]
  αt = (p.tk1 - p.tk0) / 2
  βt = p.tk0 + αt
  c₁ =  π * (2 * p.reg - 3) / 4
  c₂ =  π * αt * d.R / 2
  
  @inbounds  for I in eachindex(t)
    τ = pi * muladd(αt, t[I], βt) / 2 + c₁
    s, c = sincos(τ)
    out[2I-1] = c₂ * c
    out[2I]   = c₂ * s
  end

  return nothing
end

function nu!(out::Matrix{Float64}, d::disc, t::Vector{Float64}, k::Int) 
  p = d.pths[k]
  αt = (p.tk1 - p.tk0) / 2
  βt = p.tk0 + αt
  c₁ =  π * (2 * p.reg - 3) / 4

  @inbounds  for I in eachindex(t)
    τ = pi * muladd(αt, t[I], βt) / 2 + c₁
    s, c = sincos(τ)
    #Divided by c2, so that norm is 1.
    #But for other domains, divide by norm
    out[1, I] = c
    out[2, I] = s 
  end
  return nothing
end

"""
    DLP(d::disc, t::Float64, l::Int, tau::Float64, k::Int) -> Float64
    DLP(d::disc, t::Float64, l::Int, tau::StridedArray{Float64}, k::Int)

Double-layer kernel on the boundary:
K(t, τ) = ((γ_k(τ) - γ_l(t)) ⋅ γ'_k(τ)) / ‖γ_k(τ) - γ_l(t)‖² for k ≠ l,
and for k == l the limiting value is (π * ht / 8), where ht = tk1 - tk0
for patch k. The array method returns an array with the ***same shape*** as `tau`.
"""
function DLP(d::disc, t::Float64, l::Int, tau::Float64, k::Int)
    if k != l
        x  = gam(d, t, l)            # 2-vector
        g  = gam(d, tau, k)          # 2-vector
        gp = gamp(d, tau, k)         # 2-vector
        δ  = g .- x
        num = dot(δ, gp)
        den = δ[1]^2 + δ[2]^2
        return num / den
    else
        p  = d.pths[k]
        ht = p.tk1 - p.tk0
        return (π * ht) / 8
    end
end

function DLP!(out::StridedArray{Float64}, d::disc, t::Float64, l::Int, 
  tau::StridedArray{Float64}, k::Int, x::Vector{Float64}, G::Vector{Float64},
  GP::Vector{Float64})
  if k != l
    gam!(x, d, t, l)

    @inbounds  for i in eachindex(tau)
      ti = tau[i]
      gam!(G, d, ti, k)
      gamp!(GP, d, ti, k)

      dx = G[1] - x[1]
      dy = G[2] - x[2]
      # num = dot(Δ, GP); den = dot(Δ, Δ)
      num = muladd(dx, GP[1], dy * GP[2])
      den = muladd(dx, dx, dy * dy)
      out[i] = num / den
    end

  else
    p = d.pths[k]
    ht = p.tk1 - p.tk0
    fill!(out, (π * ht) / 8)
  end

  return nothing
end

#-----------------------
function Dmap!(out::StridedArray{Float64}, d::disc, 
  u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)

  p = d.pths[k]
  hc = p.ck1 - p.ck0
  ht = p.tk1 - p.tk0
  reg = p.reg

  αu = hc / 2
  αv = ht / 2
  βu = p.ck0 + αu      # û = αu * u + βu
  βv = p.tk0 + αv      # v̂ = αv * v + βv

  if p.reg == 5
    # constant over the patch
    fill!(out, d.L1 * d.L2 * αu * αv)
    return nothing
  end

  @inbounds for I in eachindex(out)
    û = muladd(αu, u[I], βu)
    v̂ = muladd(αv, v[I], βv)

    vt = π * v̂ / 2 + (2 * reg - 3) * π / 4
    svt, cvt = sincos(vt)
    svt *= d.R
    cvt *= d.R

    if reg == 1
      dux = cvt - d.L2 / 2
      dvx = -(π * û * svt) / 2
      duy = svt - d.L1 * (2v̂ - 1) / 2
      dvy = (1 - û) * d.L1 + (π * û * cvt) / 2
    elseif reg == 2
      dux = (cvt - d.L2 * (1 - 2v̂) / 2)
      dvx = (û - 1) * d.L2 - (π * û * svt) / 2
      duy = svt - d.L1 / 2
      dvy = (π * û * cvt) / 2
    elseif reg == 3
      dux = cvt + d.L2 / 2
      dvx = -(π * û * svt) / 2
      duy = svt + d.L1 * (2v̂ - 1) / 2
      dvy = (û - 1) * d.L1 + (π * û * cvt) / 2
    else # reg == 4
      dux = cvt + d.L2 * (1 - 2v̂) / 2
      dvx = (1 - û) * d.L2 - (π * û * svt) / 2
      duy = svt + d.L1 / 2
      dvy = (π * û * cvt) / 2
    end

    out[I] = αu * αv * abs(dux * dvy - dvx * duy)

  end

  return nothing
end

"""
A combination of mapxy! and Dmap! function (No allocations!)
The purpose of this function is to reduce computations  
related to cosine and sine's.  
"""
function mapxy_Dmap!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
  DJ::StridedArray{Float64}, d::disc, u::StridedArray{Float64},
  v::StridedArray{Float64}, k::Int)

  p   = d.pths[k]
  hc  = p.ck1 - p.ck0
  ht  = p.tk1 - p.tk0
  reg = p.reg

  αu = 0.5 * hc
  αv = 0.5 * ht
  βu = p.ck0 + αu
  βv = p.tk0 + αv
  αuv = αu * αv
  # --- curved patches reg 1..4: select formula ONCE, then loop with no if ---

  # vt = π*v̂/2 + (2*reg-3)*π/4 = (π/2)*v̂ + θ₀
  θ₀ = (2.0 * reg - 3.0) * (π / 4.0)
  c1 = (π / 2.0)

  # Common helpers
  inv2 = 0.5
  L1h = d.L1 * inv2
  L2h = d.L2 * inv2

  # We’ll implement the four region-specific affine pieces as scalars.
  # For each reg, Xx/Xy are affine in v̂, and dux/dvx/duy/dvy follow the
  # exact same structure with different signs and which term uses (2v̂-1) vs (1-2v̂).

  if reg == 1
    Xx₀ = d.A + L2h
    Xy₀ = d.B - L1h
    # Xy = B + L1*(2v̂-1)/2 = (B - L1/2) + L1*v̂ = Xy₀ + L1*v̂
    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      û = muladd(αu, u[I], βu)
      v̂ = muladd(αv, v[I], βv)

      vt = muladd(c1, v̂, θ₀)
      svt, cvt = sincos(vt)
      svtR = d.R * svt
      cvtR = d.R * cvt

      Xx = Xx₀
      Xy = muladd(d.L1, v̂, Xy₀)

      # Jacobian pieces (match your original)
      dux = cvtR - L2h
      dvx = -(c1 * û * svtR)
      duy = svtR - (muladd(d.L1, v̂, -L1h))     # svtR - L1*(2v̂-1)/2
      dvy = muladd(1.0 - û, d.L1, c1 * û * cvtR)

      Yx = d.A + cvtR
      Yy = d.B + svtR

      Zx[I] = muladd(û, Yx - Xx, Xx)
      Zy[I] = muladd(û, Yy - Xy, Xy)

      DJ[I] = αuv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 2
    # Xx = A + L2*(1-2v̂)/2 = (A + L2/2) - L2*v̂
    Xx₀ = d.A + L2h
    Xy₀ = d.B + L1h

    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      û = muladd(αu, u[I], βu)
      v̂ = muladd(αv, v[I], βv)

      vt = muladd(c1, v̂, θ₀)
      svt, cvt = sincos(vt)
      svtR = d.R * svt
      cvtR = d.R * cvt

      Xx = muladd(-d.L2, v̂, Xx₀)
      Xy = Xy₀

      dux = cvtR - (muladd(-d.L2, v̂, Xx₀) - d.A) # same as cvt - Xx shift; keep as original below
      dux = cvtR - (L2h - d.L2 * v̂)            # cvtR - L2*(1-2v̂)/2
      dvx = muladd(û - 1.0, d.L2, -(c1 * û * svtR))
      duy = svtR - L1h
      dvy = c1 * û * cvtR

      Yx = d.A + cvtR
      Yy = d.B + svtR

      Zx[I] = muladd(û, Yx - Xx, Xx)
      Zy[I] = muladd(û, Yy - Xy, Xy)

      DJ[I] = αuv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 3
    Xx₀ = d.A - L2h
    Xy₀ = d.B + L1h                 # because Xy = B - L1*(2v̂-1)/2 = (B + L1/2) - L1*v̂

    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      û = muladd(αu, u[I], βu)
      v̂ = muladd(αv, v[I], βv)

      vt = muladd(c1, v̂, θ₀)
      svt, cvt = sincos(vt)
      svtR = d.R * svt
      cvtR = d.R * cvt

      Xx = Xx₀
      Xy = muladd(-d.L1, v̂, Xy₀)

      dux = cvtR + L2h
      dvx = -(c1 * û * svtR)
      duy = svtR + (muladd(d.L1, v̂, -L1h))     # svtR + L1*(2v̂-1)/2
      dvy = muladd(û - 1.0, d.L1, c1 * û * cvtR)

      Yx = d.A + cvtR
      Yy = d.B + svtR

      Zx[I] = muladd(û, Yx - Xx, Xx)
      Zy[I] = muladd(û, Yy - Xy, Xy)

      DJ[I] = αuv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 4
    # Xx = A - L2*(1-2v̂)/2 = (A - L2/2) + L2*v̂
    Xx₀ = d.A - L2h
    Xy₀ = d.B - L1h

    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      û = muladd(αu, u[I], βu)
      v̂ = muladd(αv, v[I], βv)

      vt = muladd(c1, v̂, θ₀)
      svt, cvt = sincos(vt)
      svtR = d.R * svt
      cvtR = d.R * cvt

      Xx = muladd(d.L2, v̂, Xx₀)
      Xy = Xy₀

      dux = cvtR + (L2h - d.L2 * v̂)            # cvtR + L2*(1-2v̂)/2
      dvx = muladd(1.0 - û, d.L2, -(c1 * û * svtR))
      duy = svtR + L1h
      dvy = c1 * û * cvtR

      Yx = d.A + cvtR
      Yy = d.B + svtR

      Zx[I] = muladd(û, Yx - Xx, Xx)
      Zy[I] = muladd(û, Yy - Xy, Xy)

      DJ[I] = αuv * abs(dux * dvy - dvx * duy)
    end
  else
    # --- rectangular patch: no per-point branches ---
    fill!(DJ, d.L1 * d.L2 * αuv)

    if d.L2 >= d.L1
      # swap affine maps: xi1 from t/u, xi2 from c/v
      αu = 0.5 * ht
      αv = 0.5 * hc
      βu = p.tk0 + αu
      βv = p.ck0 + αv
    end

    x0 = d.A - 0.5 * d.L2
    y0 = d.B - 0.5 * d.L1

    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)
      Zx[I] = muladd(d.L2, xi1, x0)
      Zy[I] = muladd(d.L1, xi2, y0)
    end
  end

  return nothing
end

function chk_map(d::disc; n::Int = 32, tol::Float64 = 5e-14)
    Iex = π * d.R^4 / 2 + d.R^2 * (d.A^2 + d.B^2) * π
    return chkmap_geom(d, Iex; n = n, tol = tol)
end

"""
    jinvmap!(d::disc, u::Float64, v::Float64, r::Int) -> tuple

Inverse Jacobian of the mapping at a single reference point `(u,v) ∈ [-1,1]^2`
for the given patch. Returns a 4 Tuple which are components of matrix `J⁻¹`.
Only the region of the patch is needed as inverse are only computed over them.
ht = hc = 1
Notes:
- For the rectangular center region (`reg == 5`), the Jacobian is constant
  (affine map), so the inverse is a simple diagonal inverse (with the
  axis-swap handled when `L2 ≥ L1`, as in `mapxy`/`Dmap`).
"""
# Return (J11, J12, J21, J22) of the inverse Jacobian
function jinvmap(d::disc, u::Float64, v::Float64, r::Int)

  û = (u + 1) / 2
  v̂ = (v + 1) / 2

  vt = π * v̂ / 2 + (2 * r - 3) * π / 4
  svt = d.R * sin(vt)
  cvt = d.R * cos(vt)

  if r == 1
    dux = (cvt - d.L2 / 2) / 2
    dvx = -π * û * svt / 4
    duy = (svt - d.L1 * (2 * v̂ - 1) / 2) / 2
    dvy = (1 - û) * d.L1 / 2 + π * û * cvt / 4

  elseif r == 2
    dux = (cvt - d.L2 * (1 - 2 * v̂) / 2) / 2
    dvx = -(1 - û) * d.L2 / 2 - π * û * svt / 4
    duy = (svt - d.L1 / 2) / 2
    dvy = π * û * cvt / 4

  elseif r == 3
    dux = (cvt + d.L2 / 2) / 2
    dvx = -π * û * svt / 4
    duy = (svt + d.L1 * (2 * v̂ - 1) / 2) / 2
    dvy = (û - 1) * d.L1 / 2 + π * û * cvt / 4

  else # r == 4
    dux = (cvt + d.L2 * (1 - 2 * v̂) / 2) / 2
    dvx = (1 - û) * d.L2 / 2 - π * û * svt / 4
    duy = (svt + d.L1 / 2) / 2
    dvy = π * û * cvt / 4
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
  mapinv(d::disc, u::Float64, v::Float64, k::Int) -> Tuple{Float64,Float64}

Inverse mapping: given a physical point `(u,v)` on patch `k`, return
the reference coordinates `[t; s]` in `[-1,1]^2` such that `mapxy(d, t, s, k) = (u,v)`.

Strategy:
1. 2D Newton using `jinvmap`. If it fails to converge, fall back to a 1D root
   in `t ∈ [-30,30]` (`ApproxFun`) and reconstruct `s`.
"""

@inline function Xx(s::Float64, d::disc, r::Int)
    # ŝ ∈ [0,1]
    if r == 1
        return d.A + d.L2 / 2
    elseif r == 2
        return d.A + d.L2 * (1 - 2*s) / 2
    elseif r == 3
        return d.A - d.L2 / 2
    else # r == 4
        return d.A - d.L2 * (1 - 2*s) / 2
    end
end

@inline function Xy(s::Float64, d::disc, r::Int)
    if r == 1
        return d.B + d.L1 * (2*s - 1) / 2
    elseif r == 2
        return d.B + d.L1 / 2
    elseif r == 3
        return d.B - d.L1 * (2*s - 1) / 2
    else # r == 4
        return d.B - d.L1 / 2
    end
end

@inline function Yx(s::Float64, d::disc, r::Int)
    return d.A + d.R * cospi(s / 2 + (2*r - 3) / 4)
end

@inline function Yy(s::Float64, d::disc, r::Int)
    return d.B + d.R * sinpi(s / 2 + (2*r - 3) / 4)
end

function fill_FTable!(tbl::FTable, d::disc, r::Int)
  vmin, vmax, N = tbl.vmin, tbl.vmax, tbl.N
  h = (vmax - vmin) / (N - 1)

  P1, P2, P3 = tbl.P1, tbl.P2, tbl.P3

  if r == 1
    Xx_const = d.A + d.L2 / 2
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xy = d.B + d.L1 * (2 * v̂ - 1) / 2
      s, c = sincos(π * (v̂ / 2 - 1 / 4))
      Yx = d.A + d.R * c
      Yy = d.B + d.R * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx_const
      P3[i] = Xy * Yx - Xx_const * Yy
    end

  elseif r == 2
    Xy_const = d.B + d.L1 / 2
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A + d.L2 * (1 - 2 * v̂) / 2
      s, c = sincos(π * (v̂ / 2 + 1 / 4))
      Yx = d.A + d.R * c
      Yy = d.B + d.R * s
      P1[i] = Yy - Xy_const
      P2[i] = Yx - Xx
      P3[i] = Xy_const * Yx - Xx * Yy
    end

  elseif r == 3
    Xx_const = d.A - d.L2 / 2
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xy = d.B - d.L1 * (2 * v̂ - 1) / 2
      s, c = sincos(π * (v̂ / 2 + 3 / 4))
      Yx = d.A + d.R * c
      Yy = d.B + d.R * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx_const
      P3[i] = Xy * Yx - Xx_const * Yy
    end

  else  # r == 4
    Xy_const = d.B - d.L1 / 2
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A - d.L2 * (1 - 2 * v̂) / 2
      s, c = sincos(π * (v̂ / 2 + 5 / 4))
      Yx = d.A + d.R * c
      Yy = d.B + d.R * s
      P1[i] = Yy - Xy_const
      P2[i] = Yx - Xx
      P3[i] = Xy_const * Yx - Xx * Yy
    end
  end

  tbl.reg = r
  return tbl
end

@inline function f1I(t::Float64, s::Float64,
  d::disc, u::Float64, v::Float64, r::Int)
  ŝ = (s + 1) / 2
  t̂ = (t + 1) / 2

  Xxv, Yxv = Xx(ŝ, d, r),  Yx(ŝ, d, r)

  return (1 - t̂) * Xxv + t̂ * Yxv - u
end

@inline function f2I(t::Float64, s::Float64,
  d::disc, u::Float64, v::Float64, r::Int)
    ŝ = (s + 1) / 2
    t̂ = (t + 1) / 2

    Xyv, Yyv = Xy(ŝ, d, r), Yy(ŝ, d, r)

    return (1 - t̂) * Xyv + t̂ * Yyv - v
end

@inline function JinvI(t::Float64, s::Float64,
  d::disc, u::Float64, v::Float64, r::Int)
  return jinvmap(d, t, s, r)  # returns J11,J12,J21,J22
end

#Used in find_roots! function
@inline function f_cont(v̂::Float64,
  d::disc, u::Float64, v::Float64, r::Int)

  Xxv, Xyv = if r == 1
    d.A + d.L2 / 2, d.B + d.L1 * (2 * v̂ - 1) / 2
  elseif r == 2
    d.A + d.L2 * (1 - 2 * v̂) / 2, d.B + d.L1 / 2
  elseif r == 3
    d.A - d.L2 / 2, d.B - d.L1 * (2 * v̂ - 1) / 2
  else
    d.A - d.L2 * (1 - 2 * v̂) / 2, d.B - d.L1 / 2
  end

  phase = π * (v̂ / 2 + (2 * r - 3) / 4)

  Yxv = d.A + d.R * cos(phase)
  Yyv = d.B + d.R * sin(phase)

  return u * (Yyv - Xyv) - v * (Yxv - Xxv) + (Xyv*Yxv - Xxv*Yyv)
end

function mapinv(tbl::FTable, d::disc, u::Float64, 
  v::Float64, k::Int)::Tuple{Float64,Float64}

  p = d.pths[k]

  r = p.reg

  # central rectangle (affine inverse; watch axis swap per your map)
  if r == 5
    if d.L2 >= d.L1
      Z1 = xi_inv(((u - d.A + d.L2 / 2) / d.L2), p.tk0, p.tk1)
      Z2 = xi_inv(((v - d.B + d.L1 / 2) / d.L1), p.ck0, p.ck1)
    else
      Z1 = xi_inv(((u - d.A + d.L2 / 2) / d.L2), p.ck0, p.ck1)
      Z2 = xi_inv(((v - d.B + d.L1 / 2) / d.L1), p.tk0, p.tk1)
    end
    return Z1, Z2
  end

  # ----- curved patches reg = 1..4 -----
  # --- Stage 1: 2D Newton via Subroutines ---
  tN, sN = newtonR2D(f1I, f2I, JinvI,
    0.0, 0.0, 4, d, u, v, r; tol=1e-15)

  if tN !== :max

    x = xi_inv((1 + tN) / 2, p.ck0, p.ck1)
    y = xi_inv((1 + sN) / 2, p.tk0, p.tk1)

    return x, y

  end

  #--- Stage 2: BIS method inversion ---

  # ensure the table is filled for this region
  if tbl.reg != r
    display("hi")
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

  # choose best candidate by cost
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

  # final reference coords in [-1,1]
  t = xi_inv(zx, p.ck0, p.ck1)   # from zx
  s = xi_inv(za, p.tk0, p.tk1)   # from v̂ (root)

  return t, s

end

"""
    ptconv(d::disc,  t1::Float64, t2::Float64, idx::Int, ptdest::String)

Convert a parametric point between region ↔ patch coordinates.

Inputs:
- `t1` : 1st cooridnate of point t
- `t1` : 2nd cooridnate of point t
- `idx`: region index (if `"to_pth"`) or patch index (if `"to_reg"`)
- `ptdest`: `"to_pth"` or `"to_reg"`

Returns:
- `Tuple(Float64, Float64, out_idx::Int)`
  where `out_idx` is patch index if `"to_pth"`, region index if `"to_reg"`.
"""
function ptconv(d::disc, t1::Float64, t2::Float64, idx::Int, ptdest::String)
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
    mapinv2(d::disc, t1::Float64, t2::Float64, k2::Int, k::Int) -> Tuple{Float64,Float64}

Given a point `(u,v) = τₖ₂(t1,t2)` on patch `k2`, return its reference
coordinates on patch `k` in `[-1,1]^2`.

This differs from `mapinv` because we already know `(u,v)` comes from `(t1,t2)` on `k2`.

Test the mapinv and mapinv2 function by running the following code:
for k in 2:d.Npat
    z1, z2 = mapxy(d, 0, 0, k)
    Zx, Zy = if d.pths[k].reg == d.pths[k-1].reg
        mapinv2(d, 0, 0, k, k - 1)
    else
        mapinv(d, z1, z2, k - 1)
    end
    E1 = abs(z1 - mapx(d, Zx, Zy, k - 1))
    E2 = abs(z2 - mapy(d, Zx, Zy, k - 1))
    println(E1, " \t ", E2)
end

"""
function mapinv2(d::disc, t1::Float64, t2::Float64, k2::Int, k::Int)::Tuple{Float64,Float64}
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

function dfunc(d::disc, k::Int, t::Float64, s::Float64)::Float64

  p = d.pths[k]

  if p.reg == 5
    return 1.0
  end

  hc = p.ck1 - p.ck0
  val = (1.0 - p.ck1) + hc * t / 2
  exp = s ≥ 0.5 ? (s - 1) : s
  return val^exp

end

function dfunc!(out::StridedArray{Float64}, d::disc, k::Int, t::StridedArray{Float64}, s::Float64)

  p = d.pths[k]

  if p.reg == 5
    fill!(out, 1.0)
  else
    exp = s ≥ 0.5 ? (s - 1) : s
    αc = (p.ck1 - p.ck0) / 2
    βc = 1.0 - p.ck1
    @inbounds for i in eachindex(t)
      out[i] = (muladd(αc,t[i],1.0 - p.ck1))^exp
    end
  end

  return nothing

end

"""
  diff_map!(out::Matrix{Float64}, Zx::Matrix{Float64}, Zy::Matrix{Float64}, DJ::StridedArray{Float64},
            d::disc, u::Float64, v::Float64,
            u2::Matrix{Float64}, v2::Matrix{Float64},
            du::AbstractVector, dv::AbstractVector, k::Int;
            tol = 1e-4)

Fill `out` with ‖τ(u,v) - τ(u₂,v₂)‖ on patch `k`, where `(u,v)` are scalars and
`u2,v2` are matrices (meshgrid-like) with size `(length(du), length(dv))`.

No allocations! This map computes within itself
mapxy!(Zx, Zy, d, u2, v2, k)
  DLP!(DJ,     d, u2, v2, k)
it uses the combined function mapxy_Dmap!
"""
function diff_map!(out::Matrix{Float64},
  Zx::Matrix{Float64}, Zy::Matrix{Float64}, DJ::StridedArray{Float64},
  d::disc, u::Float64, v::Float64,
  u2::Matrix{Float64}, v2::Matrix{Float64},
  du::AbstractVector, dv::AbstractVector, k::Int;
  tol::Float64=1e-4)
  #This only works when nd_u = nd_v
  #nd_u, nd_v = length(du), length(dv)
  nd_u = size(out, 1)
  nd_v = size(out, 2)
  #@assert size(out) == (nd_u, nd_v)
  #@assert size(u2)  == size(out)
  #@assert size(v2)  == size(out)

  p = d.pths[k]
  αc = (p.ck1 - p.ck0) / 2     # Δc/2
  αt = (p.tk1 - p.tk0) / 2     # Δt/2

  mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)
  
  # --- rectangular patch (reg == 5): closed form, no mapping needed
  if p.reg == 5
    if d.L2 >= d.L1
      cDx, cDy = d.L2 * αt, d.L1 * αc
    else
      cDx, cDy = d.L2 * αc, d.L1 * αt
    end

    @inbounds for j in 1:nd_v
      dj = dv[j]
      @inbounds for i in 1:nd_u
        out[i, j] = hypot(cDx * du[i], cDy * dj)
      end
    end

    return nothing
  end

  # affine maps: [-1,1] → [ck0,ck1] × [tk0,tk1]
  û = muladd(αc, u, p.ck0 + αc)
  v̂ = muladd(αt, v, p.tk0 + αt)

  # angle/trig
  vt = muladd(π / 2, v̂, (2p.reg - 3) * (π / 4))
  svt, cvt = sincos(vt)
  svt *= d.R
  cvt *= d.R

  # ===== coefficients simplified with αc, αt =====
  # common scalars
  L1 = d.L1
  L2 = d.L2

  πa = π * αt
  πa_u = πa * û        # π * αt * û
  πa_αc = πa * αc       # π * αt * αc

  πa2 = πa * πa        # π^2 * αt^2
  πa3 = πa2 * πa       # π^3 * αt^3
  πa4 = πa2 * πa2      # π^4 * αt^4

  # x-direction higher derivatives (common across regions)
  C2x = -(πa2 * cvt) / 4        # factor in dv2x, duv2x
  C3x = (πa3 * svt) / 8        # factor in dv3x, duv3x
  C4x = (πa4 * cvt) / 16       # dv4x

  dv2x = C2x * û
  duv2x = C2x * αc
  dv3x = C3x * û
  duv3x = C3x * αc
  dv4x = C4x * û

  # y-direction higher derivatives (common across regions)
  C2y = -(πa2 * svt) / 4
  C3y = -(πa3 * cvt) / 8
  C4y = (πa4 * svt) / 16

  dv2y = C2y * û
  duv2y = C2y * αc
  dv3y = C3y * û
  duv3y = C3y * αc
  dv4y = C4y * û

  # some other shared helpers
  half_L1 = L1 / 2
  half_L2 = L2 / 2
  two_v̂_minus1 = 2v̂ - 1
  one_minus_û = 1 - û

  # now split on region
  if p.reg == 1
    # x
    dux = αc * (cvt - half_L2)
    dvx = -(πa_u * svt) / 2
    duvx = -(πa_αc * svt) / 2

    # y
    duy = αc * (svt - L1 * two_v̂_minus1 / 2)
    dvy = one_minus_û * L1 * αt + (πa_u * cvt) / 2
    duvy = (πa_αc * cvt) / 2 - αt * αc * L1

  elseif p.reg == 2
    # x
    dux = αc * (cvt - L2 * (1 - 2v̂) / 2)
    dvx = -one_minus_û * L2 * αt - (πa_u * svt) / 2
    duvx = αt * αc * L2 - (πa_αc * svt) / 2

    # y
    duy = αc * (svt - half_L1)
    dvy = (πa_u * cvt) / 2
    duvy = (πa_αc * cvt) / 2

  elseif p.reg == 3
    # x
    dux = αc * (cvt + half_L2)
    dvx = -(πa_u * svt) / 2
    duvx = -(πa_αc * svt) / 2

    # y
    duy = αc * (svt + L1 * two_v̂_minus1 / 2)
    dvy = -one_minus_û * L1 * αt + (πa_u * cvt) / 2
    duvy = (πa_αc * cvt) / 2 + αt * αc * L1

  else # p.reg == 4
    # x
    dux = αc * (cvt + L2 * (1 - 2v̂) / 2)
    dvx = one_minus_û * L2 * αt - (πa_u * svt) / 2
    duvx = -(αt * αc * L2) - (πa_αc * svt) / 2

    # y
    duy = αc * (svt + half_L1)
    dvy = (πa_u * cvt) / 2
    duvy = (πa_αc * cvt) / 2
  end
  # --- general case (reg ∈ {1,2,3,4})
  # Map the reference scalar point and the whole grid
  tux, tvy = mapxy(d, u, v, k)

  @inbounds for j in 1:nd_v
    dvj = dv[j]
    dvj2 = dvj * dvj
    dvj3 = dvj2 * dvj
    dvj4 = dvj2 * dvj2
    @inbounds for i in 1:nd_u
      uu = u2[i, j]
      vv = v2[i, j]
      if (abs(u - uu) < tol) && (abs(v - vv) < tol)
        dui = du[i]
        # near the evaluation point: Taylor fixup
        #Below computes the x and y coordinate of:
        #τ(u,v) - τ(u-du,v-dv) = dτᵤ(u,v) * du + dτᵥ(u,v) * dv - (...)
        Dx = (dui * dux + dvj * dvx) -
             (dui * dvj * duvx + dvj2 * (dv2x/2)) +
             (dui * dvj2 * (duv2x/2) + dvj3 * (dv3x/6)) -
             (dui * dvj3 * (duv3x/6) + dvj4 * (dv4x/24))

        Dy = (dui * duy + dvj * dvy) -
             (dui * dvj * duvy + dvj2 * (dv2y/2)) +
             (dui * dvj2 * (duv2y/2) + dvj3 * (dv3y/6)) -
             (dui * dvj3 * (duv3y/6) + dvj4 * (dv4y/24))

        out[i, j] = hypot(Dx, Dy)
      else
        # far: direct geometric difference
        out[i, j] = hypot(tux - Zx[i, j], tvy - Zy[i, j])
      end
    end
  end

  return nothing
end

"""
  diff_rmap!(out::Matrix{Float64}, x::Matrix{Float64}, Zy::Matrix{Float64}, DJ::StridedArray{Float64}
            d::disc, u::Float64, v::Float64,
            u2::Matrix{Float64}, v2::Matrix{Float64}, r::Matrix{Float64},
            du::AbstractVector, dv::AbstractVector, k::Int;
            tol = 1e-4)

Compute ‖(τ(u,v) - τ(u₂,v₂)) / r‖ for the `k`-th patch, with
`u₂ = u - r .* du`, `v₂ = v - r .* dv`.

- `(u,v)` are scalars.
- `u2, v2, r` are matrices of the **same size**.
- `du, dv` are vectors.
No allocations! This map computes within itself
mapxy!(Zx, Zy, d, u2, v2, k)
  DLP!(DJ,     d, u2, v2, k)
it uses the combined function mapxy_Dmap!
  Unlike diff_map! function, Zx, Zy in regions
  with affine mappings are not updated!
  But the Jacobian DJ is always updated.
"""
function diff_rmap!(out::Matrix{Float64},
  Zx::Matrix{Float64}, Zy::Matrix{Float64}, DJ::StridedArray{Float64},
  d::disc, u::Float64, v::Float64,
  u2::Matrix{Float64}, v2::Matrix{Float64}, r::Matrix{Float64},
  du::AbstractVector, dv::AbstractVector, k::Int;
  tol::Float64=1e-5)

  nt = size(out, 1)   # angular
  nr = size(out, 2)   # radial
  #passed all these assert tests!
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
  αc = (p.ck1 - p.ck0) / 2     # Δc/2
  αt = (p.tk1 - p.tk0) / 2     # Δt/2

  # --- rectangular patch (reg == 5): closed form, no mapping needed
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
        for j in 1:nr
            out[i, j] = hD
        end
    end

    return nothing
  end

  mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)

  # affine maps: [-1,1] → [ck0,ck1] × [tk0,tk1]
  û = muladd(αc, u, p.ck0 + αc)
  v̂ = muladd(αt, v, p.tk0 + αt)

  # angle/trig
  vt = muladd(π / 2, v̂, (2p.reg - 3) * (π / 4))
  svt, cvt = sincos(vt)
  svt *= d.R
  cvt *= d.R

  # ===== coefficients simplified with αc, αt =====
  # common scalars
  L1 = d.L1
  L2 = d.L2

  πa = π * αt
  πa_u = πa * û        # π * αt * û
  πa_αc = πa * αc       # π * αt * αc

  πa2 = πa * πa        # π^2 * αt^2
  πa3 = πa2 * πa       # π^3 * αt^3
  πa4 = πa2 * πa2      # π^4 * αt^4

  # x-direction higher derivatives (common across regions)
  C2x = -(πa2 * cvt) / 4        # factor in dv2x, duv2x
  C3x = (πa3 * svt) / 8        # factor in dv3x, duv3x
  C4x = (πa4 * cvt) / 16       # dv4x

  dv2x = C2x * û
  duv2x = C2x * αc
  dv3x = C3x * û
  duv3x = C3x * αc
  dv4x = C4x * û

  # y-direction higher derivatives (common across regions)
  C2y = -(πa2 * svt) / 4
  C3y = -(πa3 * cvt) / 8
  C4y = (πa4 * svt) / 16

  dv2y = C2y * û
  duv2y = C2y * αc
  dv3y = C3y * û
  duv3y = C3y * αc
  dv4y = C4y * û

  # some other shared helpers
  half_L1 = L1 / 2
  half_L2 = L2 / 2
  two_v̂_minus1 = 2v̂ - 1
  one_minus_û = 1 - û

  # now split on region
  if p.reg == 1
    # x
    dux = αc * (cvt - half_L2)
    dvx = -(πa_u * svt) / 2
    duvx = -(πa_αc * svt) / 2

    # y
    duy = αc * (svt - L1 * two_v̂_minus1 / 2)
    dvy = one_minus_û * L1 * αt + (πa_u * cvt) / 2
    duvy = (πa_αc * cvt) / 2 - αt * αc * L1

  elseif p.reg == 2
    # x
    dux = αc * (cvt - L2 * (1 - 2v̂) / 2)
    dvx = -one_minus_û * L2 * αt - (πa_u * svt) / 2
    duvx = αt * αc * L2 - (πa_αc * svt) / 2

    # y
    duy = αc * (svt - half_L1)
    dvy = (πa_u * cvt) / 2
    duvy = (πa_αc * cvt) / 2

  elseif p.reg == 3
    # x
    dux = αc * (cvt + half_L2)
    dvx = -(πa_u * svt) / 2
    duvx = -(πa_αc * svt) / 2

    # y
    duy = αc * (svt + L1 * two_v̂_minus1 / 2)
    dvy = -one_minus_û * L1 * αt + (πa_u * cvt) / 2
    duvy = (πa_αc * cvt) / 2 + αt * αc * L1

  else # p.reg == 4
    # x
    dux = αc * (cvt + L2 * (1 - 2v̂) / 2)
    dvx = one_minus_û * L2 * αt - (πa_u * svt) / 2
    duvx = -(αt * αc * L2) - (πa_αc * svt) / 2

    # y
    duy = αc * (svt + half_L1)
    dvy = (πa_u * cvt) / 2
    duvy = (πa_αc * cvt) / 2
  end
  # --- general case (reg ∈ {1,2,3,4})
  # Map the reference scalar point and the whole grid
  tux, tvy = mapxy(d, u, v, k)

  @inbounds for i in 1:nt
    dvi = dv[i]
    dui = du[i]
    @inbounds for j in 1:nr
      uu = u2[i, j]
      vv = v2[i, j]
      if (abs(u - uu) < tol) && (abs(v - vv) < tol)
        r1 = dvi *  r[i, j]
        r2 = r1 * r1
        r3 = r2 * r1
        # near the evaluation point: Taylor fixup
        #Below computes the x and y coordinate of:
        #(τ(u,v) - τ(u-r*du,v-r*dv))/r = dτᵤ(u,v) * du +  dτᵥ(u,v) * dv - r*(...)
        Dx = (dui * dux + dvi * dvx) -
             r1 * (dui * duvx + dvi * (dv2x / 2)) +
             r2 * (dui * (duv2x / 2) + dvi * (dv3x / 6)) -
             r3 * (dui * (duv3x / 6) + dvi * (dv4x / 24))

        Dy = (dui * duy + dvi * dvy) -
             r1 * (dui * duvy + dvi * (dv2y / 2)) +
             r2 * (dui * (duv2y / 2) + dvi * (dv3y / 6)) -
             r3 * (dui * (duv3y / 6) + dvi * (dv4y / 24))

        out[i, j] = hypot(Dx, Dy)
      else
        # far: direct geometric difference
        out[i, j] = hypot(tux - Zx[i, j], tvy - Zy[i, j]) / r[i, j]
      end
    end
  end

  return nothing
end


"""
    Dwall(d::disc, u::Float64, v::Float64, k::Int) -> dvx, dvy

Derivative of the wall curve w.r.t. `v` at the singular side `u = ±1`
for the `k`-th (non-quadrilateral) patch.

Returns a Tuple `dvx, dvy`.
"""
function Dwall(d::disc, u::Float64, v::Float64, k::Int)::Tuple{Float64, Float64}

    p  = d.pths[k]                    

    # Δ’s
    hc = p.ck1 - p.ck0
    ht = p.tk1 - p.tk0

    # affine maps: [-1,1] → [ck0,ck1] and [tk0,tk1]
    v̂ = ht*(v + 1)/2 + p.tk0
    û = hc*(u + 1)/2 + p.ck0

    # angular map and radius-scaled trig
    vt  = π*v̂/2 + (2*p.reg - 3)*π/4
    svt = d.R * sin(vt)
    cvt = d.R * cos(vt)

    if p.reg == 1
        dvx = -π * ht * û * svt / 4
        dvy = (1 - û) * d.L1 * ht / 2 + π * ht * û * cvt / 4
    elseif p.reg == 2
        dvx = -(1 - û) * d.L2 * ht / 2 - π * ht * û * svt / 4
        dvy =  π * ht * û * cvt / 4
    elseif p.reg == 3
        dvx = -π * ht * û * svt / 4
        dvy = (û - 1) * d.L1 * ht / 2 + π * ht * û * cvt / 4
    else # p.reg == 4
        dvx =  (1 - û) * d.L2 * ht / 2 - π * ht * û * svt / 4
        dvy =  π * ht * û * cvt / 4
    end

    return dvx, dvy
end


"""
    boundquad!(P::SubArray{Float64}, d::disc, k::Int) -> nothing

Find a **bounding quadrilateral** for the k-th patch.

Mutates P into in `[x1,y1,x2,y2,x3,y3,x4,y4]` in **counter-clockwise** order

If the patch is already a quadrilateral, return its fourcorners directly; 
otherwise:
1) build the two opposite sides `u = -1` and `u = +1`,
2) along each side, intersect the **tangent line** at `τ(u,t)` with the
   two straight sides to find the extremal intersections in parameter space,
3) compute the Float64-space intersection points and assemble P.
"""
#|| always requires the left operand to be a Bool
#|| can return non-boolean values.
#false || "ok"    # returns "ok"
#true || "nope"   # returns true
function boundquad!(P::SubArray{Float64}, d::disc, k::Int)

  p = d.pths[k]

  # Four corners (tuples immediately unpacked into scalars)
  L1p1x, L1p1y = mapm1(d, -1.0, k)   # (u=-1, v=-1)
  L1p2x, L1p2y = mapp1(d, -1.0, k)   # (u=+1, v=-1)
  L2p1x, L2p1y = mapm1(d, 1.0, k)   # (u=-1, v=+1)
  L2p2x, L2p2y = mapp1(d, 1.0, k)   # (u=+1, v=+1)

  if p.reg == 5
    # already a quadrilateral: corners in CCW order
    P[1] = L1p1x
    P[2] = L1p1y
    P[3] = L1p2x
    P[4] = L1p2y
    P[5] = L2p2x
    P[6] = L2p2y
    P[7] = L2p1x
    P[8] = L2p1y
    return nothing
  end
  # Common t grid
  N = 101
  tpts = range(-1.0, 1.0; length=N)

  # ---------- side away from boundary: u = -1 ----------
  best_score1 = Inf
  best_idx1 = 0

  for (i, t) in enumerate(tpts)
    # point on u = -1 wall
    Cx, Cy = mapm1(d, t, k)
    dvx, dvy = Dwall(d, -1.0, t, k)

    # z1, z2 = parameters returned by tanginterp (u, v or similar)
    z1, z2 = tanginterp(
      L1p1x, L1p1y, L1p2x, L1p2y,
      Cx, Cy, dvx, dvy,
      L2p1x, L2p1y, L2p2x, L2p2y)

    ok, s = SVum1(z1, z2)
    if ok && s < best_score1
      best_score1 = s
      best_idx1 = i
    end
  end

  if best_idx1 != 0
    to = tpts[best_idx1]
    # recompute geometry at best t
    Cx, Cy = mapm1(d, to, k)
    dvx, dvy = Dwall(d, -1.0, to, k)

    P1x, P1y, P2x, P2y = tanginterx(
      L1p1x, L1p1y, L1p2x, L1p2y,
      Cx, Cy, dvx, dvy,
      L2p1x, L2p1y, L2p2x, L2p2y)
  else
    # fallback
    P1x, P1y = L1p1x, L1p1y
    P2x, P2y = L2p1x, L2p1y
  end

  # ---------- side closer to boundary: u = +1 ----------
  best_score2 = Inf
  best_idx2 = 0

  for (i, t) in enumerate(tpts)
    Cx, Cy = mapp1(d, t, k)
    dvx, dvy = Dwall(d, 1.0, t, k)

    z1, z2 = tanginterp(
      L1p1x, L1p1y, L1p2x, L1p2y,
      Cx, Cy, dvx, dvy,
      L2p1x, L2p1y, L2p2x, L2p2y)

    ok, s = SVup1(z1, z2)
    if ok && s < best_score2
      best_score2 = s
      best_idx2 = i
    end
  end

  if best_idx2 != 0
    to = tpts[best_idx2]
    Cx, Cy = mapp1(d, to, k)
    dvx, dvy = Dwall(d, 1.0, to, k)

    P3x, P3y, P4x, P4y = tanginterx(
      L1p1x, L1p1y, L1p2x, L1p2y,
      Cx, Cy, dvx, dvy,
      L2p1x, L2p1y, L2p2x, L2p2y)
  else
    P3x, P3y = L1p2x, L1p2y
    P4x, P4y = L2p2x, L2p2y
  end

  # assemble in CCW order
  P[1] = P1x
  P[2] = P1y
  P[3] = P3x
  P[4] = P3y
  P[5] = P4x
  P[6] = P4y
  P[7] = P2x
  P[8] = P2y

  return nothing
end


"""
    refine!(d::disc, Nc::Int, Nt::Int, K::AbstractVector{<:Int})

Subdivide each patch in `K` into `Nc * Nt` subpatches (split in `c` and `t`),
update `d.pths`, `d.Npat`, recompute `d.kd`, and rebuild `d.Qpts` / `d.Qptsbd`.

Notes
- Patches are re-sorted lexicographically by `(reg, ck0, ck1, tk0, tk1)`
- `Qpts` uses `boundquad(d,k)`; `Qptsbd` uses `boundquadbd(d,k)` for patches in `d.kd`.
"""
function refine!(d::disc, Nc::Int, Nt::Int, K::Vector{Int})
    @assert Nc ≥ 1 && Nt ≥ 1 "Nc and Nt must be ≥ 1"
    @assert !isempty(K) "K (set of patches to refine) is empty"

    # create the subdivided children
    children = Patch[]
    sizehint!(children, length(K) * Nc * Nt)

    for k in K
        p = d.pths[k]
        cc = range(p.ck0, p.ck1; length = Nc + 1)
        tt = range(p.tk0, p.tk1; length = Nt + 1)
        for i in 1:Nc, j in 1:Nt
            push!(children, Patch(p.reg, cc[i], cc[i+1], tt[j], tt[j+1]))
        end
    end

    # merge and sort
    d.pths = vcat([d.pths[i] for i in 1:d.Npat if !(i in K)], children)
    sort!(d.pths, by = q -> (q.reg, q.ck0, q.ck1, q.tk0, q.tk1))

    # update counts
    d.Npat = length(d.pths)

    # --- recompute kd (boundary-touching patches) and quadrilaterals
    d.kd = [k for k in 1:d.Npat if d.pths[k].reg != 5 && d.pths[k].ck1 == 1.0]

    # Qpts: 8 × Npat, each column is boundquad of that patch
    d.Qpts = zeros(Float64, 8, d.Npat)
    @inbounds for k in 1:d.Npat
      @views V = d.Qpts[:, k]
      boundquad!(V, d, k)
    end
 
    # Qptsbd: 8 × |kd|, boundary quads for boundary patches (in kd order)
    d.Qptsbd = zeros(Float64, 8, length(d.kd))
    @inbounds for (ℓ, k) in enumerate(d.kd)
      @views V = d.Qptsbd[:, ℓ]
      boundquadbd!(V, d, k)
    end

    return d
end

function Base.show(io::IO, d::disc)

  println(io, "disc with properties:")
  println(io, "  (A , B )  = (", d.A,",", d.B,")")
  println(io, "  (L₁, L₂)  = (", d.L1,",", d.L2,")")
  println(io, "  Radius of disc = ", d.R)
  println(io, "  No. of holes   = ", d.nh)

  if length(d.kd) <= 6
    L = "Int[" * join(d.kd, ' ') * "]"
  else
    head = join(d.kd[1:3], ' ')
    tail = join(d.kd[end-3+1:end], ' ')
    L = "Int[$head … $tail]"
  end
  
  println(io, "  kd     = ", L)
  println(io, "  Npat   = ", d.Npat)
  println(io, "  pths   = Vector{Patch} (", length(d.pths), " patches)")
  println(io, "  Qpts   = ", size(d.Qpts, 1), "×", size(d.Qpts, 2), " Matrix")
  println(io, "  Qptsbd = ", size(d.Qptsbd, 1), "×", size(d.Qptsbd, 2), " Matrix")

end
