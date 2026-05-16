"""
Mutable struct kite

A  — x-coordinate of kite center (Float64)

B  — y-coordinate of kite center (Float64)

P  — Distortion of kite from an ellipse (Float64)

R₁ — Scale of the kite in x axis (Float64)

R₂ — Scale of the kite in y axis (Float64)

L₁ — length of rectangle's y-axis side (Float64)

L₂ — length of rectangle's x-axis side (Float64)

L₃ — Decides placement of rectangles on lines (Float64)

nₕ — number of holes (nonnegative Int) 

θ₁ — A parameter used to constuct the patches in kite. (Float64)

θ₂ — A parameter used to constuct the patches in kite. (Float64)

θ₃ — A parameter used to constuct the patches in kite. (Float64)

θ₄ — A parameter used to constuct the patches in kite. (Float64)

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
struct regionparams_kite
  e₁x::Float64
  e₁y::Float64
  e₂x::Float64
  e₂y::Float64
  q₁x::Float64
  q₁y::Float64
  q₂x::Float64
  q₂y::Float64
  α₁x::Float64
  α₁y::Float64
  α₂x::Float64
  α₂y::Float64
  nme::Float64 #Norm of e

end

mutable struct kite <: abstractdomain

  A::Float64
  B::Float64
  P::Float64
  R1::Float64
  R2::Float64
  L1::Float64
  L2::Float64
  L3::Float64
  tht1::Float64
  tht2::Float64
  tht3::Float64
  tht4::Float64
  nh::Int
  kd::Vector{Int}
  Npat::Int
  pths::Vector{Patch}   # 5 * Npat
  Qpts::Matrix{Float64}   # 8 * Npat
  Qptsbd::Matrix{Float64} # 8 * Npat
  RP::regionparams_kite

  function kite(; b,
    a=nothing, A=nothing, B=nothing,
    P=nothing, R1=nothing, R2=nothing, 
    L1=nothing, L2=nothing, L3=nothing, 
    tht1=nothing, tht2=nothing, 
    tht3=nothing, tht4=nothing, 
    ck=nothing, tk=nothing)
    #Construct an intsance for this structure

    """
    Required:
      - b :: Vector{Int} of length 12 (Col), positive entries
    
    Optional (all keywords):
      - a  :: Vector{Int} length 12 (defaults to ceil.(2*b/3))
      - A,B:: Center of the kite (defaults (0,0))
      - R1 :: Float64 (default 1.0)
      - R2 :: Float64 (default 1.5)
      - P  :: Float64 (default 1.3)
      - L1 :: Float64 (default 0.35*R2)
      - L2 :: Float64 (default 0.35*R1)
      - L3 :: Float64 (default 0.46)
      - θ₁ :: Float64 (default 0.90)
      - θ₂ :: Float64 (default 0.25)
      - θ₃ :: Float64 (default 1.25)
      - θ₄ :: Float64 (default 2.40)
      - ck :: vector{vector{Float64}} (default equispaced per a)
      - tk :: vector{vector{Float64}} (default equispaced per b)
    """
    @assert isa(b, AbstractVector{Int}) && length(b) == 12

    #The kite has obviously no holes
    nh = 0

    a = isnothing(a) ? ceil.(Int, 2 .* b ./ 3) : collect(Int, a)

    #By default, kite at origin
    A = something(A, 0.0)
    B = something(B, 0.0)

    P = something(P, 1.3)

    R1 = something(R1, 1//1)
    R2 = something(R2, 3//2)

    L1 = something(L1, 7//20 * R1)
    L2 = something(L2, 7//20 * R2)
    L3 = something(L3, 0.46)

    tht1= something(tht1, 0.90)
    tht2= something(tht2, 0.25)
    tht3= something(tht3, 1.25)
    tht4= something(tht4, 2.40)

    #Sum up all the subpatches

    Npat = dot(a,b)

    # ck, tk as vector-of-vectors with *no padding*:
    #   ck[k] has length a[k] + 1
    #   tk[k] has length b[k] + 1
    ck = something(ck, [(0:a[k]) ./ a[k] for k in 1:12])
    tk = something(tk, [(0:b[k]) ./ b[k] for k in 1:12])

    # Preallocate pths and kd
    pths = Vector{Patch}(undef, Npat)

    #Starting with the assumption that all interioir 
    #patches touch the boundary
    kd = Vector{Int}(undef, Npat)   # upper bound length
    nkd = 0                         # actual count

    # Fill pths and kd in a single pass
    idx = 1
    @inbounds for k in 1:12
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

          if (p.reg != 11 && p.reg != 12) && p.ck1 == 1.0
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

    tp = (tht1 + tht2) / 2
    tm = (tht1 - tht2) / 2
    # e₁
    e₁x = 2 * R1 * cos(tp) * cos(tm) - P * sin(2 * tm) * sin(2 * tp)
    e₁y = 2 * R2 * cos(tp) * sin(tm)

    e₂x, e₂y = e₁x, -e₁y    # e₂
    nme = sqrt(e₁x^2 + e₁y^2) # ||e₁||=||e₂||
    q₁x, q₁y = e₂y, e₂x     # q₁
    q₂x, q₂y = q₁x, -q₁y    # q₂

    #α₁ and α₂
    α₁x = L3 * e₂x - R1 * cos(tht2) - P * (sin(tht2)^2)
    α₁y = L3 * e₂y - R2 * sin(tht2)
    α₂x, α₂y = α₁x, L3 * e₁y + R2 * sin(tht2)

    RP = regionparams_kite(e₁x, e₁y, e₂x, e₂y, q₁x, q₁y, q₂x, q₂y, α₁x, α₁y, α₂x, α₂y, nme)

    d = new(A, B, P, R1, R2, L1, L2, L3, tht1, tht2, tht3, tht4, nh, kd, Npat, pths, Qpts, Qptsbd, RP)
    
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

function mapx(d::kite, u::Float64, v::Float64, k::Int)::Float64

  # p denotes information of kth patch
  p = d.pths[k]

  # Difference in c values
  hc = p.ck1 - p.ck0

  xi1 = hc * u / 2 + (p.ck1 + p.ck0) / 2

  # Difference in t values
  ht = p.tk1 - p.tk0

  xi2 = ht * v / 2 + (p.tk1 + p.tk0) / 2

  ex = d.RP.e₁x
  ey = d.RP.e₁y
  nme= d.RP.nme
  qx = d.RP.q₁x
  αx = d.RP.α₁x

  if p.reg != 11 && p.reg != 12

    if p.reg == 1

      X = d.A + αx + (2 * xi2 - 1) * d.L1 * ex / (2*nme) + d.L2 * qx / nme

      th = xi2 * (d.tht4 - d.tht3) - d.tht4

    elseif p.reg == 2

      X = d.A + αx + d.L1 * ex / (2*nme) + (1 - xi2) * d.L2 * qx / nme

      th = xi2 * (d.tht3 - d.tht1) - d.tht3

    elseif p.reg == 3

      X = d.A + (1 - xi2) * αx

      th = xi2 * d.tht1 - d.tht1

    elseif p.reg == 4

      X = d.A + xi2 * αx 

      th = xi2 * d.tht1

    elseif p.reg == 5

      X = d.A + αx + d.L1 * ex/(2*nme) + xi2 * d.L2 * qx / nme

      th = xi2 * (d.tht3 - d.tht1) + d.tht1;

    elseif p.reg == 6

      X = d.A + αx - (2 * xi2 - 1) * d.L1 * ex / (2*nme) + d.L2 * qx / nme

      th = xi2 * (d.tht4 - d.tht3) + d.tht3

    elseif p.reg == 7

      X = d.A + αx - d.L1 * ex / (2 * nme) + (1 - xi2) * d.L2 * qx / nme

      th = xi2 * (π - d.tht2 - d.tht4) + d.tht4

    elseif p.reg == 8

      X = d.A + (1 - xi2) * αx 

      th = xi2 * d.tht2 + π - d.tht2

    elseif p.reg == 9

      X = d.A + xi2 * αx

      th = xi2 * d.tht2 + π

    elseif p.reg == 10

      X = d.A + αx - d.L1 * ex / (2*nme) + d.L2 * xi2 * qx / nme

      th = xi2 * (π - d.tht2 - d.tht4) + π + d.tht2

    end

    s, c = sincos(th)

    Y = d.A + d.R1 * c - d.P * s * s

    return (1 - xi1) * X + xi1 * Y

  else

    if d.L2 >= d.L1

      xi1 = ht * u / 2 + (p.tk1 + p.tk0) / 2

      xi2 = hc * v / 2 + (p.ck1 + p.ck0) / 2

    end

    return d.A + αx + d.L2 * xi1 * qx / nme +(2 * xi2 - 1) * d.L1 * ex/(2*nme)

  end

end

function mapy(d::kite, u::Float64, v::Float64, k::Int)::Float64

  # p denotes information of kth patch
  p = d.pths[k]

  # Difference in c values
  hc = p.ck1 - p.ck0
  xi1 = hc * u / 2 + (p.ck1 + p.ck0) / 2

  # Difference in t values
  ht = p.tk1 - p.tk0
  xi2 = ht * v / 2 + (p.tk1 + p.tk0) / 2

  # Pull precomputed region params
  ey  = d.RP.e₁y
  nme = d.RP.nme
  qy  = d.RP.q₁y          
  α1y = d.RP.α₁y
  α2y = d.RP.α₂y

  if p.reg != 11 && p.reg != 12

    if p.reg == 1

      X = d.B + α1y + (2 * xi2 - 1) * d.L1 * (-ey) / (2 * nme) + d.L2 * (-qy) / nme

      th = xi2 * (d.tht4 - d.tht3) - d.tht4

    elseif p.reg == 2

      X = d.B + α1y + d.L1 * (-ey) / (2 * nme) + (1 - xi2) * d.L2 * (-qy) / nme

      th = xi2 * (d.tht3 - d.tht1) - d.tht3

    elseif p.reg == 3

      X = d.B + (1 - xi2) * α1y

      th = xi2 * d.tht1 - d.tht1

    elseif p.reg == 4

      X = d.B + xi2 * α2y

      th = xi2 * d.tht1

    elseif p.reg == 5

      X = d.B + α2y + d.L1 * ey / (2 * nme) + xi2 * d.L2 * qy / nme

      th = xi2 * (d.tht3 - d.tht1) + d.tht1

    elseif p.reg == 6

      X = d.B + α2y - (2 * xi2 - 1) * d.L1 * ey / (2 * nme) + d.L2 * qy / nme

      th = xi2 * (d.tht4 - d.tht3) + d.tht3

    elseif p.reg == 7

      X = d.B + α2y - d.L1 * ey / (2 * nme) + (1 - xi2) * d.L2 * qy / nme

      th = xi2 * (π - d.tht2 - d.tht4) + d.tht4

    elseif p.reg == 8

      X = d.B + (1 - xi2) * α2y

      th = xi2 * d.tht2 + π - d.tht2

    elseif p.reg == 9

      X = d.B + xi2 * α1y

      th = xi2 * d.tht2 + π

    elseif p.reg == 10

      X = d.B + α1y - d.L1 * (-ey) / (2 * nme) + d.L2 * xi2 * (-qy) / nme

      th = xi2 * (π - d.tht2 - d.tht4) + π + d.tht2

    end

    s = sin(th)
    Y = d.B + d.R2 * s

    return (1 - xi1) * X + xi1 * Y

  elseif p.reg == 11

    if d.L2 >= d.L1
      xi1 = ht * u / 2 + (p.tk1 + p.tk0) / 2
      xi2 = hc * v / 2 + (p.ck1 + p.ck0) / 2
    end

    return d.B + α1y + d.L2 * xi1 * (-qy) / nme + (2 * xi2 - 1) * d.L1 * (-ey) / (2 * nme)

  else # p.reg == 12

    if d.L2 >= d.L1
      xi1 = ht * u / 2 + (p.tk1 + p.tk0) / 2
      xi2 = hc * v / 2 + (p.ck1 + p.ck0) / 2
    end

    return d.B + α2y + d.L2 * xi1 * qy / nme + (2 * xi2 - 1) * d.L1 * ey / (2 * nme)

  end

end

function mapxy(d::kite, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64}

  p = d.pths[k]

  # base affine maps into (xi1, xi2)
  hc  = p.ck1 - p.ck0
  ht  = p.tk1 - p.tk0
  αu = hc / 2
  αv = ht / 2
  βu = p.ck0 + αu
  βv = p.tk0 + αv
  xi1 = muladd(αu, u, βu)
  xi2 = muladd(αv, v, βv)

  # precomputed geometry params
  ex  = d.RP.e₁x
  ey  = d.RP.e₁y
  nme = d.RP.nme
  qx  = d.RP.q₁x
  qy  = d.RP.q₁y
  αx  = d.RP.α₁x
  α1y = d.RP.α₁y
  α2y = d.RP.α₂y

  reg = p.reg

  if reg != 11 && reg != 12

    if reg == 1
      Xx = d.A + αx + (2*xi2 - 1) * d.L1 * ex / (2*nme) + d.L2 * qx / nme
      Xy = d.B + α1y + (2*xi2 - 1) * d.L1 * (-ey) / (2*nme) + d.L2 * (-qy) / nme
      th = xi2 * (d.tht4 - d.tht3) - d.tht4

    elseif reg == 2
      Xx = d.A + αx + d.L1 * ex / (2*nme) + (1 - xi2) * d.L2 * qx / nme
      Xy = d.B + α1y + d.L1 * (-ey) / (2*nme) + (1 - xi2) * d.L2 * (-qy) / nme
      th = xi2 * (d.tht3 - d.tht1) - d.tht3

    elseif reg == 3
      Xx = d.A + (1 - xi2) * αx
      Xy = d.B + (1 - xi2) * α1y
      th = xi2 * d.tht1 - d.tht1

    elseif reg == 4
      Xx = d.A + xi2 * αx
      Xy = d.B + xi2 * α2y
      th = xi2 * d.tht1

    elseif reg == 5
      Xx = d.A + αx + d.L1 * ex / (2*nme) + xi2 * d.L2 * qx / nme
      Xy = d.B + α2y + d.L1 * ey / (2*nme) + xi2 * d.L2 * qy / nme
      th = xi2 * (d.tht3 - d.tht1) + d.tht1

    elseif reg == 6
      Xx = d.A + αx - (2*xi2 - 1) * d.L1 * ex / (2*nme) + d.L2 * qx / nme
      Xy = d.B + α2y - (2*xi2 - 1) * d.L1 * ey / (2*nme) + d.L2 * qy / nme
      th = xi2 * (d.tht4 - d.tht3) + d.tht3

    elseif reg == 7
      Xx = d.A + αx - d.L1 * ex / (2*nme) + (1 - xi2) * d.L2 * qx / nme
      Xy = d.B + α2y - d.L1 * ey / (2*nme) + (1 - xi2) * d.L2 * qy / nme
      th = xi2 * (π - d.tht2 - d.tht4) + d.tht4

    elseif reg == 8
      Xx = d.A + (1 - xi2) * αx
      Xy = d.B + (1 - xi2) * α2y
      th = xi2 * d.tht2 + π - d.tht2

    elseif reg == 9
      Xx = d.A + xi2 * αx
      Xy = d.B + xi2 * α1y
      th = xi2 * d.tht2 + π

    elseif reg == 10
      Xx = d.A + αx - d.L1 * ex / (2*nme) + d.L2 * xi2 * qx / nme
      Xy = d.B + α1y - d.L1 * (-ey) / (2*nme) + d.L2 * xi2 * (-qy) / nme
      th = xi2 * (π - d.tht2 - d.tht4) + π + d.tht2
    end

    s, c = sincos(th)

    Yx = d.A + d.R1 * c - d.P * s * s
    Yy = d.B + d.R2 * s

    # (1-xi1)*X + xi1*Y  ==  X + xi1*(Y-X)
    Zx = muladd(xi1, Yx - Xx, Xx)
    Zy = muladd(xi1, Yy - Xy, Xy)

    return Zx, Zy

  else
    # regions 11/12 (lines), possibly swap parameterization if L2 >= L1
    if d.L2 >= d.L1
      xi1 = muladd(ht * 0.5, u, (p.tk1 + p.tk0) * 0.5)
      xi2 = muladd(hc * 0.5, v, (p.ck1 + p.ck0) * 0.5)
    end

    Zx = d.A + αx  + d.L2 * xi1 * qx / nme + (2*xi2 - 1) * d.L1 * ex / (2*nme)

    if reg == 11
      Zy = d.B + α1y + d.L2 * xi1 * (-qy) / nme + (2*xi2 - 1) * d.L1 * (-ey) / (2*nme)
    else # reg == 12
      Zy = d.B + α2y + d.L2 * xi1 * ( qy) / nme + (2*xi2 - 1) * d.L1 * ( ey) / (2*nme)
    end

    return Zx, Zy
  end
end

function mapxy!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
  d::kite, u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)
  p = d.pths[k]
  reg = p.reg

  # Affine maps to xi1, xi2 (default: xi1 from c/u, xi2 from t/v)
  hc = p.ck1 - p.ck0
  ht = p.tk1 - p.tk0
  αu = 0.5 * hc
  αv = 0.5 * ht
  βu = p.ck0 + αu
  βv = p.tk0 + αv

  # Geometry
  ex = d.RP.e₁x
  ey = d.RP.e₁y
  nme = d.RP.nme
  qx = d.RP.q₁x
  qy = d.RP.q₁y
  αx = d.RP.α₁x
  α1y = d.RP.α₁y
  α2y = d.RP.α₂y

  # Common scalings
  invn = inv(nme)
  hinvn = 0.5 * invn
  sL1x = d.L1 * ex * hinvn          # L1*ex/(2*nme)
  sL1y = d.L1 * ey * hinvn          # L1*ey/(2*nme)
  sL2x = d.L2 * qx * invn               # L2*qx/nme
  sL2y = d.L2 * qy * invn               # L2*qy/nme

  if reg == 1
    Xx₀ = d.A + αx + sL2x
    Xy₀ = d.B + α1y - sL2y
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

    @inbounds for I in eachindex(Zx, Zy, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      t = muladd(2.0, xi2, -1.0)
      Xx = muladd(t, sL1x, Xx₀)
      Xy = muladd(t, -sL1y, Xy₀)
      th = muladd(xi2, Δth, th0)

      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s

      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)
    end

  elseif reg == 2
    Xx₀ = d.A + αx + sL1x
    Xy₀ = d.B + α1y - sL1y
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

    @inbounds for I in eachindex(Zx, Zy, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      Xx = muladd(1.0 - xi2, sL2x, Xx₀)
      Xy = muladd(1.0 - xi2, -sL2y, Xy₀)
      th = muladd(xi2, Δth, th0)

      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s

      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)
    end

  elseif reg == 3
    Δth = d.tht1
    th0 = -d.tht1

    @inbounds for I in eachindex(Zx, Zy, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      Xx = muladd(1.0 - xi2, αx, d.A)
      Xy = muladd(1.0 - xi2, α1y, d.B)
      th = muladd(xi2, Δth, th0)

      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s

      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)
    end

  elseif reg == 4
    Δth = d.tht1
    th0 = 0.0

    @inbounds for I in eachindex(Zx, Zy, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      Xx = muladd(xi2, αx, d.A)
      Xy = muladd(xi2, α2y, d.B)
      th = muladd(xi2, Δth, th0)

      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s

      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)
    end

  elseif reg == 5
    Xx₀ = d.A + αx + sL1x
    Xy₀ = d.B + α2y + sL1y
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

    @inbounds for I in eachindex(Zx, Zy, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      Xx = muladd(xi2, sL2x, Xx₀)
      Xy = muladd(xi2, sL2y, Xy₀)
      th = muladd(xi2, Δth, th0)

      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s

      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)
    end

  elseif reg == 6
    Xx₀ = d.A + αx + sL2x
    Xy₀ = d.B + α2y + sL2y
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

    @inbounds for I in eachindex(Zx, Zy, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      t = muladd(2.0, xi2, -1.0)
      Xx = muladd(-t, sL1x, Xx₀)
      Xy = muladd(-t, sL1y, Xy₀)
      th = muladd(xi2, Δth, th0)

      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s

      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)
    end

  elseif reg == 7
    Xx₀ = d.A + αx - sL1x
    Xy₀ = d.B + α2y - sL1y
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

    @inbounds for I in eachindex(Zx, Zy, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      Xx = muladd(1.0 - xi2, sL2x, Xx₀)
      Xy = muladd(1.0 - xi2, sL2y, Xy₀)
      th = muladd(xi2, Δth, th0)

      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s

      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)
    end

  elseif reg == 8
    Δth = d.tht2
    th0 = π - d.tht2

    @inbounds for I in eachindex(Zx, Zy, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      Xx = muladd(1.0 - xi2, αx, d.A)
      Xy = muladd(1.0 - xi2, α2y, d.B)
      th = muladd(xi2, Δth, th0)

      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s

      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)
    end

  elseif reg == 9
    Δth = d.tht2
    th0 = π

    @inbounds for I in eachindex(Zx, Zy, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      Xx = muladd(xi2, αx, d.A)
      Xy = muladd(xi2, α1y, d.B)
      th = muladd(xi2, Δth, th0)

      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s

      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)
    end

  elseif reg == 10
    Xx₀ = d.A + αx - sL1x
    Xy₀ = d.B + α1y + sL1y
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

    @inbounds for I in eachindex(Zx, Zy, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      Xx = muladd(xi2, sL2x, Xx₀)
      Xy = muladd(xi2, -sL2y, Xy₀)
      th = muladd(xi2, Δth, th0)

      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s

      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)
    end

  elseif reg == 11
    # swap parameterization if L2 >= L1:
    # xi1 from t/u, xi2 from c/v
    if d.L2 >= d.L1
      αu = 0.5 * ht
      αv = 0.5 * hc
      βu = p.tk0 + αu
      βv = p.ck0 + αv
    end

    @inbounds for I in eachindex(Zx, Zy, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)
      t = muladd(2.0, xi2, -1.0)

      Zx[I] = d.A + αx  + (d.L2 * xi1) * qx * invn + t * sL1x

      Zy[I] = d.B + α1y + (d.L2 * xi1) * (-qy) * invn + t * (-sL1y)
    end

  elseif reg == 12
    # swap parameterization if L2 >= L1:
    # xi1 from t/u, xi2 from c/v
    if d.L2 >= d.L1
      αu = 0.5 * ht
      αv = 0.5 * hc
      βu = p.tk0 + αu
      βv = p.ck0 + αv
    end

    @inbounds for I in eachindex(Zx, Zy, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      t = muladd(2.0, xi2, -1.0)

      Zx[I] = d.A + αx  + (d.L2 * xi1) * qx * invn + t * sL1x

      Zy[I] = d.B + α2y + (d.L2 * xi1) * (qy) * invn + t * (sL1y)
    end

  end

  return nothing
end

function draw(d::kite, flag = nothing; L::Int = 33, show::Bool = true)

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
  gamx(d::kite, t::Float64, k::Int)
  gamx!(out::StridedArray{Float64}, d::kite, t::StridedArray{Float64}, k)
First (x) coordinate of the boundary parametrization γ(t).

- `gamx(d, t, k)` computes the **right boundary** of patch `k`, where `t ∈ [-1,1]`.

`t` may be a scalar or an array; the return has the same shape.
"""
function gamx(d::kite, t::Float64, k::Int)::Float64
  p = d.pths[k]
  # Map t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
  xi2 = muladd(0.5 * (p.tk1 - p.tk0), t, 0.5 * (p.tk0 + p.tk1))

  th = if p.reg == 1
    xi2 * (d.tht4 - d.tht3) - d.tht4

  elseif p.reg == 2
    xi2 * (d.tht3 - d.tht1) - d.tht3

  elseif p.reg == 3
    xi2 * d.tht1 - d.tht1

  elseif p.reg == 4
    xi2 * d.tht1

  elseif p.reg == 5
    xi2 * (d.tht3 - d.tht1) + d.tht1

  elseif p.reg == 6
    xi2 * (d.tht4 - d.tht3) + d.tht3

  elseif p.reg == 7
    xi2 * (π - d.tht2 - d.tht4) + d.tht4

  elseif p.reg == 8
    xi2 * d.tht2 + π - d.tht2

  elseif p.reg == 9
    xi2 * d.tht2 + π

  elseif p.reg == 10
    xi2 * (π - d.tht2 - d.tht4) + π + d.tht2

  else
    throw(ArgumentError("gamx is defined only for regions 1–10; got reg=$(p.reg)"))

  end

  s, c = sincos(th)

  return d.A + d.R1 * c - d.P * s * s
end

function gamx!(out::StridedArray{Float64}, d::kite, t::StridedArray{Float64}, k::Int)
  p = d.pths[k]

  ht = p.tk1 - p.tk0
  αt = ht / 2
  βt = p.tk0 + αt

  if p.reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

  elseif p.reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

  elseif p.reg == 3
    Δth = d.tht1
    th0 = -d.tht1

  elseif p.reg == 4
    Δth = d.tht1
    th0 = 0.0

  elseif p.reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

  elseif p.reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

  elseif p.reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

  elseif p.reg == 8
    Δth = d.tht2
    th0 = π - d.tht2

  elseif p.reg == 9
    Δth = d.tht2
    th0 = π

  elseif p.reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

  end

  @inbounds for I in eachindex(out, t)
    xi2 = muladd(αt, t[I], βt)

    s, c = sincos(muladd(xi2, Δth, th0))

    out[I] = d.A + d.R1 * c - d.P * s * s

  end

  return nothing
end
"""
    gamy(d::kite, t::Float64, k::Int)
    gamy!(out::StridedArray{Float64}, d::kite, t::StridedArray{Float64}, k)

Second (y) coordinate of the boundary parametrization γ(t).

- `gamy(d, t, k)` computes the right boundary of patch `k`, where `t ∈ [-1,1]`.

`t` can be a scalar `Float64` or an array; the return has the same shape.
"""
function gamy(d::kite, t::Float64, k::Int)::Float64
  p = d.pths[k]

  # Map t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
  xi2 = (p.tk1 - p.tk0) * t / 2 + (p.tk0 + p.tk1) / 2

  th = if p.reg == 1
    xi2 * (d.tht4 - d.tht3) - d.tht4
  elseif p.reg == 2
    xi2 * (d.tht3 - d.tht1) - d.tht3
  elseif p.reg == 3
    xi2 * d.tht1 - d.tht1
  elseif p.reg == 4
    xi2 * d.tht1
  elseif p.reg == 5
    xi2 * (d.tht3 - d.tht1) + d.tht1
  elseif p.reg == 6
    xi2 * (d.tht4 - d.tht3) + d.tht3
  elseif p.reg == 7
    xi2 * (π - d.tht2 - d.tht4) + d.tht4
  elseif p.reg == 8
    xi2 * d.tht2 + π - d.tht2
  elseif p.reg == 9
    xi2 * d.tht2 + π
  elseif p.reg == 10
    xi2 * (π - d.tht2 - d.tht4) + π + d.tht2
  else
    throw(ArgumentError("gamy is defined only for regions 1–10; got reg=$(p.reg)"))
  end

  return d.B + d.R2 * sin(th)
end

function gamy!(out::StridedArray{Float64}, d::kite, t::StridedArray{Float64}, k::Int)
  p = d.pths[k]

  ht = p.tk1 - p.tk0
  αt = ht / 2
  βt = p.tk0 + αt

  if p.reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

  elseif p.reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

  elseif p.reg == 3
    Δth = d.tht1
    th0 = -d.tht1

  elseif p.reg == 4
    Δth = d.tht1
    th0 = 0.0

  elseif p.reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

  elseif p.reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

  elseif p.reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

  elseif p.reg == 8
    Δth = d.tht2
    th0 = π - d.tht2

  elseif p.reg == 9
    Δth = d.tht2
    th0 = π

  elseif p.reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

  end

  @inbounds for I in eachindex(out, t)
    xi2 = muladd(αt, t[I], βt)
    out[I] = d.B + d.R2 * sin(muladd(xi2, Δth, th0))

  end

  return nothing
end

function gam(d::kite, t::Float64, k::Int)::Tuple{Float64,Float64}

  p = d.pths[k]
  # Map t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
  xi2 = muladd(0.5 * (p.tk1 - p.tk0), t, 0.5 * (p.tk0 + p.tk1))

  th = if p.reg == 1
    xi2 * (d.tht4 - d.tht3) - d.tht4

  elseif p.reg == 2
    xi2 * (d.tht3 - d.tht1) - d.tht3

  elseif p.reg == 3
    xi2 * d.tht1 - d.tht1

  elseif p.reg == 4
    xi2 * d.tht1

  elseif p.reg == 5
    xi2 * (d.tht3 - d.tht1) + d.tht1

  elseif p.reg == 6
    xi2 * (d.tht4 - d.tht3) + d.tht3

  elseif p.reg == 7
    xi2 * (π - d.tht2 - d.tht4) + d.tht4

  elseif p.reg == 8
    xi2 * d.tht2 + π - d.tht2

  elseif p.reg == 9
    xi2 * d.tht2 + π

  elseif p.reg == 10
    xi2 * (π - d.tht2 - d.tht4) + π + d.tht2

  else
    throw(ArgumentError("gam is defined only for regions 1–10; got reg=$(p.reg)"))

  end

  s, c = sincos(th)

  Yx = d.A + d.R1 * c - d.P * s * s

  Yy = d.B + d.R2 * s

  return Yx, Yy

end

function gam!(out::Vector{Float64}, d::kite, t::Float64, k::Int)

  p = d.pths[k]
  # Map t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
  xi2 = muladd(0.5 * (p.tk1 - p.tk0), t, 0.5 * (p.tk0 + p.tk1))

  th = if p.reg == 1
    xi2 * (d.tht4 - d.tht3) - d.tht4

  elseif p.reg == 2
    xi2 * (d.tht3 - d.tht1) - d.tht3

  elseif p.reg == 3
    xi2 * d.tht1 - d.tht1

  elseif p.reg == 4
    xi2 * d.tht1

  elseif p.reg == 5
    xi2 * (d.tht3 - d.tht1) + d.tht1

  elseif p.reg == 6
    xi2 * (d.tht4 - d.tht3) + d.tht3

  elseif p.reg == 7
    xi2 * (π - d.tht2 - d.tht4) + d.tht4

  elseif p.reg == 8
    xi2 * d.tht2 + π - d.tht2

  elseif p.reg == 9
    xi2 * d.tht2 + π

  elseif p.reg == 10
    xi2 * (π - d.tht2 - d.tht4) + π + d.tht2

  else
    throw(ArgumentError("gam is defined only for regions 1–10; got reg=$(p.reg)"))

  end

  s, c = sincos(th)

  out[1] = d.A + d.R1 * c - d.P * s * s

  out[2] = d.B + d.R2 * s

  return nothing

end

function gam!(out::Matrix{Float64}, d::kite, t::Vector{Float64}, k::Int)

  p = d.pths[k]

  ht = p.tk1 - p.tk0
  αt = ht / 2
  βt = p.tk0 + αt

  if p.reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

  elseif p.reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

  elseif p.reg == 3
    Δth = d.tht1
    th0 = -d.tht1

  elseif p.reg == 4
    Δth = d.tht1
    th0 = 0.0

  elseif p.reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

  elseif p.reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

  elseif p.reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

  elseif p.reg == 8
    Δth = d.tht2
    th0 = π - d.tht2

  elseif p.reg == 9
    Δth = d.tht2
    th0 = π

  elseif p.reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

  end

  @inbounds for I in eachindex(t)
    xi2 = muladd(αt, t[I], βt)
    s, c = sincos(muladd(xi2, Δth, th0))
    out[1, I] = d.A + d.R1 * c - d.P * s * s
    out[2, I] = d.B + d.R2 * s
  end

  return nothing

end

function drawbd(d::kite, flag = true; L::Int = 33, show::Bool = true)

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
    dgamx(d::kite, t::Float64, k::Int) -> Float64

Derivative of the first coordinate of the boundary parametrization.

- With `k`: derivative of the **right boundary** of patch `k`
  (i.e., `∂/∂t mapx(d, 1, t, k)` for `t ∈ [-1,1]`).
"""
function dgamx(d::kite, t::Float64, k::Int)::Float64
  p = d.pths[k]

  ht = p.tk1 - p.tk0
  αt = ht / 2
  βt = p.tk0 + αt

  if p.reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

  elseif p.reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

  elseif p.reg == 3
    Δth = d.tht1
    th0 = -d.tht1

  elseif p.reg == 4
    Δth = d.tht1
    th0 = 0.0

  elseif p.reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

  elseif p.reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

  elseif p.reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

  elseif p.reg == 8
    Δth = d.tht2
    th0 = π - d.tht2

  elseif p.reg == 9
    Δth = d.tht2
    th0 = π

  elseif p.reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

  else
    throw(ArgumentError("dgamx is defined only for regions 1–10; got reg=$(p.reg)"))

  end

  xi2 = muladd(αt, t, βt)

  s, c = sincos(muladd(xi2, Δth, th0))

  return ht * Δth * (- d.R1 * s - 2 * d.P * s * c) / 2
end

function dgamx!(out::StridedArray{Float64}, d::kite, t::StridedArray{Float64}, k::Int)

  p = d.pths[k]

  ht = p.tk1 - p.tk0
  αt = ht / 2
  βt = p.tk0 + αt

  if p.reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

  elseif p.reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

  elseif p.reg == 3
    Δth = d.tht1
    th0 = -d.tht1

  elseif p.reg == 4
    Δth = d.tht1
    th0 = 0.0

  elseif p.reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

  elseif p.reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

  elseif p.reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

  elseif p.reg == 8
    Δth = d.tht2
    th0 = π - d.tht2

  elseif p.reg == 9
    Δth = d.tht2
    th0 = π

  elseif p.reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

  else
    throw(ArgumentError("dgamx! is defined only for regions 1–10; got reg=$(p.reg)"))

  end

  @inbounds for I in eachindex(out, t)

    xi2 = muladd(αt, t[I], βt)

    s, c = sincos(muladd(xi2, Δth, th0))

    out[I] = ht * Δth * (-d.R1 * s - 2 * d.P * s * c) / 2

  end

  return nothing
end

"""
    dgamy(d::kite, t::Float64, k::Int) -> Float64

Derivative of the second coordinate of the boundary parametrization γ(t).

- With `k`: derivative of the **right boundary** of patch `k`
  (`∂/∂t mapy(d, 1, t, k)` for `t ∈ [-1,1]`).
"""
function dgamy(d::kite, t::Float64, k::Int)::Float64
  p = d.pths[k]

  ht = p.tk1 - p.tk0
  αt = 0.5 * ht
  βt = p.tk0 + αt

  if p.reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

  elseif p.reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

  elseif p.reg == 3
    Δth = d.tht1
    th0 = -d.tht1

  elseif p.reg == 4
    Δth = d.tht1
    th0 = 0.0

  elseif p.reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

  elseif p.reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

  elseif p.reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

  elseif p.reg == 8
    Δth = d.tht2
    th0 = π - d.tht2

  elseif p.reg == 9
    Δth = d.tht2
    th0 = π

  elseif p.reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

  else
    throw(ArgumentError("dgamy is defined only for regions 1–10; got reg=$(p.reg)"))

  end

  xi2 = muladd(αt, t, βt)

  # y(θ) = B + R2*sin(θ)
  # dy/dθ = R2*cos(θ)
  # dθ/dt = (ht/2)*Δth
  return (0.5 * ht * Δth) * (d.R2 * cos(muladd(xi2, Δth, th0)))
end

function dgamy!(out::StridedArray{Float64}, d::kite, t::StridedArray{Float64}, k::Int)
  p = d.pths[k]

  ht = p.tk1 - p.tk0
  αt = 0.5 * ht
  βt = p.tk0 + αt

  if p.reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

  elseif p.reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

  elseif p.reg == 3
    Δth = d.tht1
    th0 = -d.tht1

  elseif p.reg == 4
    Δth = d.tht1
    th0 = 0.0

  elseif p.reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

  elseif p.reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

  elseif p.reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

  elseif p.reg == 8
    Δth = d.tht2
    th0 = π - d.tht2

  elseif p.reg == 9
    Δth = d.tht2
    th0 = π

  elseif p.reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

  else
    throw(ArgumentError("dgamy! is defined only for regions 1–10; got reg=$(p.reg)"))

  end

  sc = 0.5 * ht * Δth * d.R2

  @inbounds for I in eachindex(out, t)
    xi2 = muladd(αt, t[I], βt)
    θ = muladd(xi2, Δth, th0)
    out[I] = sc * cos(θ)
  end

  return nothing
end

function dgam(d::kite, t::Float64, k::Int)::Tuple{Float64,Float64}
  p = d.pths[k]

  ht = p.tk1 - p.tk0
  αt = 0.5 * ht
  βt = p.tk0 + αt

  if p.reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

  elseif p.reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

  elseif p.reg == 3
    Δth = d.tht1
    th0 = -d.tht1

  elseif p.reg == 4
    Δth = d.tht1
    th0 = 0.0

  elseif p.reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

  elseif p.reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

  elseif p.reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

  elseif p.reg == 8
    Δth = d.tht2
    th0 = π - d.tht2

  elseif p.reg == 9
    Δth = d.tht2
    th0 = π

  elseif p.reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

  else
    throw(ArgumentError("dgam is defined only for regions 1–10; got reg=$(p.reg)"))
    
  end

  xi2 = muladd(αt, t, βt)
  θ = muladd(xi2, Δth, th0)

  s, c = sincos(θ)

  scale = 0.5 * ht * Δth

  # dx/dθ = -R1*sinθ - P*sin(2θ) = -R1*s - 2P*s*c
  dYx = scale * (-d.R1 * s - 2.0 * d.P * s * c)

  # dy/dθ = R2*cosθ
  dYy = scale * (d.R2 * c)

  return dYx, dYy
end

function dgam!(out::Vector{Float64}, d::kite, t::Float64, k::Int)
  p = d.pths[k]

  ht = p.tk1 - p.tk0
  αt = 0.5 * ht
  βt = p.tk0 + αt

  if p.reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

  elseif p.reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

  elseif p.reg == 3
    Δth = d.tht1
    th0 = -d.tht1

  elseif p.reg == 4
    Δth = d.tht1
    th0 = 0.0

  elseif p.reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

  elseif p.reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

  elseif p.reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

  elseif p.reg == 8
    Δth = d.tht2
    th0 = π - d.tht2

  elseif p.reg == 9
    Δth = d.tht2
    th0 = π

  elseif p.reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

  else
    throw(ArgumentError("dgam! is defined only for regions 1–10; got reg=$(p.reg)"))
    
  end

  xi2 = muladd(αt, t, βt)
  θ = muladd(xi2, Δth, th0)

  s, c = sincos(θ)

  scale = 0.5 * ht * Δth

  # dx/dθ = -R1*sinθ - P*sin(2θ) = -R1*s - 2P*s*c
  out[1] = scale * (-d.R1 * s - 2.0 * d.P * s * c)

  # dy/dθ = R2*cosθ
  out[2] = scale * (d.R2 * c)

  return nothing
end

function dgam!(out::Matrix{Float64}, d::kite, t::Vector{Float64}, k::Int)

  p = d.pths[k]

  ht = p.tk1 - p.tk0
  αt = 0.5 * ht
  βt = p.tk0 + αt

  if p.reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

  elseif p.reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

  elseif p.reg == 3
    Δth = d.tht1
    th0 = -d.tht1

  elseif p.reg == 4
    Δth = d.tht1
    th0 = 0.0

  elseif p.reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

  elseif p.reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

  elseif p.reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

  elseif p.reg == 8
    Δth = d.tht2
    th0 = π - d.tht2

  elseif p.reg == 9
    Δth = d.tht2
    th0 = π

  elseif p.reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

  else
    throw(ArgumentError("dgam! is defined only for regions 1–10; got reg=$(p.reg)"))

  end

  scale = 0.5 * ht * Δth

  @inbounds for I in eachindex(t)
    xi2 = muladd(αt, t[I], βt)
    θ = muladd(xi2, Δth, th0)

    s, c = sincos(θ)

    out[1, I] = scale * (-d.R1 * s - 2.0 * d.P * s * c)
    out[2, I] = scale * (d.R2 * c)
  end

  return nothing
end

function gamp(d::kite, t::Float64, k::Int)::Tuple{Float64,Float64}
    return dgamy(d, t, k), -dgamx(d, t, k)
end

function gamp!(out::Matrix{Float64}, d::kite, t::Vector{Float64}, k::Int)
  #Outputs dy, -dx, the sign is changed in dgamx coeffs!
  p = d.pths[k]

  ht = p.tk1 - p.tk0
  αt = 0.5 * ht
  βt = p.tk0 + αt

  if p.reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

  elseif p.reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

  elseif p.reg == 3
    Δth = d.tht1
    th0 = -d.tht1

  elseif p.reg == 4
    Δth = d.tht1
    th0 = 0.0

  elseif p.reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

  elseif p.reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

  elseif p.reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

  elseif p.reg == 8
    Δth = d.tht2
    th0 = π - d.tht2

  elseif p.reg == 9
    Δth = d.tht2
    th0 = π

  elseif p.reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

  else
    throw(ArgumentError("gamp! is defined only for regions 1–10; got reg=$(p.reg)"))

  end

  scale = 0.5 * ht * Δth

  @inbounds for I in eachindex(t)
    xi2 = muladd(αt, t[I], βt)
    θ = muladd(xi2, Δth, th0)

    s, c = sincos(θ)

    out[1, I] = scale * (d.R2 * c)                    # dy/dt
    out[2, I] = scale * (d.R1 * s + 2.0 * d.P * s * c)# -dx/dt
  end

  return nothing
end

function nu!(out::Matrix{Float64}, d::kite, t::Vector{Float64}, k::Int)
  p = d.pths[k]

  ht = p.tk1 - p.tk0
  αt = 0.5 * ht
  βt = p.tk0 + αt

  if p.reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

  elseif p.reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

  elseif p.reg == 3
    Δth = d.tht1
    th0 = -d.tht1

  elseif p.reg == 4
    Δth = d.tht1
    th0 = 0.0

  elseif p.reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

  elseif p.reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

  elseif p.reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

  elseif p.reg == 8
    Δth = d.tht2
    th0 = π - d.tht2

  elseif p.reg == 9
    Δth = d.tht2
    th0 = π

  elseif p.reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

  else
    throw(ArgumentError("nu! is defined only for regions 1–10; got reg=$(p.reg)"))

  end

  scale = 0.5 * ht * Δth

  @inbounds for I in eachindex(t)
    xi2 = muladd(αt, t[I], βt)
    θ = muladd(xi2, Δth, th0)

    s, c = sincos(θ)

    # dx/dθ = -R1*sinθ - P*sin(2θ) = -R1*s - 2P*s*c
    dx = scale * (-d.R1 * s - 2.0 * d.P * s * c)
    dy = scale * (d.R2 * c)

    S = sqrt(dx^2 + dy^2)

    out[1, I] = dy / S
    out[2, I] = -dx / S
  end

  return nothing
end

# This function finds s such that γ_l(t) = γ_k(s),
# allowing s to lie outside [-1, 1].
function bdinv(d::kite, t::Float64, l::Int, k::Int)::Float64
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
      xil * (d.tht4 - d.tht3) - d.tht4
    elseif plr == 2
      xil * (d.tht3 - d.tht1) - d.tht3
    elseif plr == 3
      xil * d.tht1 - d.tht1
    elseif plr == 4
      xil * d.tht1
    elseif plr == 5
      xil * (d.tht3 - d.tht1) + d.tht1
    elseif plr == 6
      xil * (d.tht4 - d.tht3) + d.tht3
    elseif plr == 7
      xil * (π - d.tht2 - d.tht4) + d.tht4
    elseif plr == 8
      xil * d.tht2 + π - d.tht2
    elseif plr == 9
      xil * d.tht2 + π
    elseif plr == 10
      xil * (π - d.tht2 - d.tht4) + π + d.tht2
    else
      throw(ArgumentError("bdinv source region must be 1–10; got reg=$plr"))
    end

    # Get angular map for target patch k:
    #     theta = ak * xi + bk
    ak, bk = if pkr == 1
      d.tht4 - d.tht3, -d.tht4
    elseif pkr == 2
      d.tht3 - d.tht1, -d.tht3
    elseif pkr == 3
      d.tht1, -d.tht1
    elseif pkr == 4
      d.tht1, 0.0
    elseif pkr == 5
      d.tht3 - d.tht1, d.tht1
    elseif pkr == 6
      d.tht4 - d.tht3, d.tht3
    elseif pkr == 7
      π - d.tht2 - d.tht4, d.tht4
    elseif pkr == 8
      d.tht2, π - d.tht2
    elseif pkr == 9
      d.tht2, π
    elseif pkr == 10
      π - d.tht2 - d.tht4, π + d.tht2
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
    DLP(d::kite, t::Float64, tau::StridedArray{Float64}, k::Int)

Double-layer kernel on the boundary:
K(t, τ) = ((γ_k(τ) - γ_k(t))⋅  γᵖᵉʳᵖ_k(τ)) / ‖γ_k(τ) - γ_k(t)‖² 
and for k == l the limiting value is taken for patch k.
The array method returns an array with the ***same shape*** as `tau`.
"""
function DLP!(out::StridedArray{Float64}, d::kite, t::Float64,
  tau::StridedArray{Float64}, k::Int, x::Vector{Float64},
  G::Vector{Float64}, GP::Vector{Float64})

  p = d.pths[k]

  ht = p.tk1 - p.tk0
  αt = 0.5 * ht
  βt = p.tk0 + αt

  reg = p.reg

  if reg == 1
    ak = d.tht4 - d.tht3
    bk = -d.tht4

  elseif reg == 2
    ak = d.tht3 - d.tht1
    bk = -d.tht3

  elseif reg == 3
    ak = d.tht1
    bk = -d.tht1

  elseif reg == 4
    ak = d.tht1
    bk = 0.0

  elseif reg == 5
    ak = d.tht3 - d.tht1
    bk = d.tht1

  elseif reg == 6
    ak = d.tht4 - d.tht3
    bk = d.tht3

  elseif reg == 7
    ak = π - d.tht2 - d.tht4
    bk = d.tht4

  elseif reg == 8
    ak = d.tht2
    bk = π - d.tht2

  elseif reg == 9
    ak = d.tht2
    bk = π

  elseif reg == 10
    ak = π - d.tht2 - d.tht4
    bk = π + d.tht2

  else
    throw(ArgumentError("DLP! self term defined only for regions 1–10; got reg=$reg"))

  end

  R1 = d.R1
  R2 = d.R2
  P = d.P

  # t mapped to [tk0, tk1]
  tm = muladd(αt, t, βt)
  pref = 0.25 * ht * ak * R2
  halfak = 0.5 * ak

  @inbounds for i in eachindex(tau, out)
    τm = muladd(αt, tau[i], βt)

    tt = muladd(halfak, tm + τm, bk)
    Δ = halfak * (τm - tm)

    st, ct = sincos(tt)
    sΔ, cΔ = sincos(Δ)

    ct2 = ct * ct
    s2t = 2.0 * st * ct

    c_thet_tau = muladd(ct, cΔ, -st * sΔ)

    vv1 = muladd(P * s2t, cΔ, R1 * st)
    vv2 = R2 * ct

    den = muladd(vv1, vv1, vv2 * vv2)
    num = muladd(2.0 * P * c_thet_tau, ct2, R1)

    out[i] = pref * num / den
  end

  return nothing
end

#-----------------------
function dergam(d::kite, t::Float64)::Tuple{Float64,Float64,Float64,Float64} 

  s,c = sincos(t)

  gx, gy   = d.R1 * c - d.P * s * s, d.R2 * s

  dgx, dgy = -d.R1 * s - 2.0 * d.P * s * c, d.R2 * c 

  return gx, gy, dgx, dgy

end

function Dmap!(out::StridedArray{Float64}, d::kite,
  u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)

  p = d.pths[k]
  hc = p.ck1 - p.ck0
  ht = p.tk1 - p.tk0
  reg = p.reg

  αu = 0.5 * hc
  αv = 0.5 * ht
  βu = p.ck0 + αu      # û = αu*u + (ck0+ck1)/2
  βv = p.tk0 + αv      # v̂ = αv*v + (tk0+tk1)/2

  if reg == 11 || reg == 12
    fill!(out, d.L1 * d.L2 * αu * αv)
    return nothing
  end

  # Precomputed geometry params 
  e1x = d.RP.e₁x
  e1y = d.RP.e₁y
  e2x = d.RP.e₂x
  e2y = d.RP.e₂y
  nme = d.RP.nme
  q1x = d.RP.q₁x
  q1y = d.RP.q₁y
  q2x = d.RP.q₂x
  q2y = d.RP.q₂y
  α1x = d.RP.α₁x
  α1y = d.RP.α₁y
  α2x = d.RP.α₂x
  α2y = d.RP.α₂y

  if reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

    @inbounds for I in eachindex(out, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)
      th = muladd(xi2, Δth, th0)

      gx, gy, dgx, dgy = dergam(d, th)

      dux = gx - α1x - (2 * xi2 - 1) * d.L1 * e2x / (2 * nme) - d.L2 * q2x / nme
      duy = gy - α1y - (2 * xi2 - 1) * d.L1 * e2y / (2 * nme) - d.L2 * q2y / nme

      dvx = d.L1 * (1 - xi1) * e2x / nme + Δth * xi1 * dgx
      dvy = d.L1 * (1 - xi1) * e2y / nme + Δth * xi1 * dgy

      out[I] = αu * αv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

    @inbounds for I in eachindex(out, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)
      th = muladd(xi2, Δth, th0)

      gx, gy, dgx, dgy = dergam(d, th)

      dux = gx - α1x - d.L1 * e2x / (2 * nme) - (1 - xi2) * d.L2 * q2x / nme
      duy = gy - α1y - d.L1 * e2y / (2 * nme) - (1 - xi2) * d.L2 * q2y / nme

      dvx = -d.L2 * (1 - xi1) * q2x / nme + Δth * xi1 * dgx
      dvy = -d.L2 * (1 - xi1) * q2y / nme + Δth * xi1 * dgy

      out[I] = αu * αv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 3
    Δth = d.tht1
    th0 = -d.tht1

    @inbounds for I in eachindex(out, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)
      th = muladd(xi2, Δth, th0)

      gx, gy, dgx, dgy = dergam(d, th)

      dux = gx - (1 - xi2) * α1x
      duy = gy - (1 - xi2) * α1y

      dvx = -(1 - xi1) * α1x + Δth * xi1 * dgx
      dvy = -(1 - xi1) * α1y + Δth * xi1 * dgy

      out[I] = αu * αv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 4
    Δth = d.tht1
    th0 = 0.0

    @inbounds for I in eachindex(out, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)
      th = muladd(xi2, Δth, th0)

      gx, gy, dgx, dgy = dergam(d, th)

      dux = gx - xi2 * α2x
      duy = gy - xi2 * α2y

      dvx = (1 - xi1) * α2x + Δth * xi1 * dgx
      dvy = (1 - xi1) * α2y + Δth * xi1 * dgy

      out[I] = αu * αv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

    @inbounds for I in eachindex(out, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)
      th = muladd(xi2, Δth, th0)

      gx, gy, dgx, dgy = dergam(d, th)

      dux = gx - α2x - d.L1 * e1x / (2 * nme) - xi2 * d.L2 * q1x / nme
      duy = gy - α2y - d.L1 * e1y / (2 * nme) - xi2 * d.L2 * q1y / nme

      dvx = d.L2 * (1 - xi1) * q1x / nme + Δth * xi1 * dgx
      dvy = d.L2 * (1 - xi1) * q1y / nme + Δth * xi1 * dgy

      out[I] = αu * αv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

    @inbounds for I in eachindex(out, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)
      th = muladd(xi2, Δth, th0)

      gx, gy, dgx, dgy = dergam(d, th)

      dux = gx - α2x + (2 * xi2 - 1) * d.L1 * e1x / (2 * nme) - d.L2 * q1x / nme
      duy = gy - α2y + (2 * xi2 - 1) * d.L1 * e1y / (2 * nme) - d.L2 * q1y / nme

      dvx = -d.L1 * (1 - xi1) * e1x / nme + Δth * xi1 * dgx
      dvy = -d.L1 * (1 - xi1) * e1y / nme + Δth * xi1 * dgy

      out[I] = αu * αv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

    @inbounds for I in eachindex(out, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)
      th = muladd(xi2, Δth, th0)

      gx, gy, dgx, dgy = dergam(d, th)

      dux = gx - α2x + d.L1 * e1x / (2 * nme) - (1 - xi2) * d.L2 * q1x / nme
      duy = gy - α2y + d.L1 * e1y / (2 * nme) - (1 - xi2) * d.L2 * q1y / nme

      dvx = -d.L2 * (1 - xi1) * q1x / nme + Δth * xi1 * dgx
      dvy = -d.L2 * (1 - xi1) * q1y / nme + Δth * xi1 * dgy

      out[I] = αu * αv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 8
    Δth = d.tht2
    th0 = π - d.tht2

    @inbounds for I in eachindex(out, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)
      th = muladd(xi2, Δth, th0)

      gx, gy, dgx, dgy = dergam(d, th)

      dux = gx - (1 - xi2) * α2x
      duy = gy - (1 - xi2) * α2y

      dvx = -(1 - xi1) * α2x + Δth * xi1 * dgx
      dvy = -(1 - xi1) * α2y + Δth * xi1 * dgy

      out[I] = αu * αv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 9
    Δth = d.tht2
    th0 = π

    @inbounds for I in eachindex(out, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)
      th = muladd(xi2, Δth, th0)

      gx, gy, dgx, dgy = dergam(d, th)

      dux = gx - xi2 * α1x
      duy = gy - xi2 * α1y

      dvx = (1 - xi1) * α1x + Δth * xi1 * dgx
      dvy = (1 - xi1) * α1y + Δth * xi1 * dgy

      out[I] = αu * αv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

    @inbounds for I in eachindex(out, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)
      th = muladd(xi2, Δth, th0)

      gx, gy, dgx, dgy = dergam(d, th)

      dux = gx - α1x + d.L1 * e2x / (2 * nme) - xi2 * d.L2 * q2x / nme
      duy = gy - α1y + d.L1 * e2y / (2 * nme) - xi2 * d.L2 * q2y / nme

      dvx = d.L2 * (1 - xi1) * q2x / nme + Δth * xi1 * dgx
      dvy = d.L2 * (1 - xi1) * q2y / nme + Δth * xi1 * dgy

      out[I] = αu * αv * abs(dux * dvy - dvx * duy)
    end

  else
    throw(ArgumentError("Dmap! is defined only for regions 1–12; got reg=$reg"))
  end

  return nothing
end

"""
A combination of mapxy! and Dmap! function (No allocations!)
The purpose of this function is to reduce computations  
related to cosine and sine's.  
"""
function mapxy_Dmap!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
  DJ::StridedArray{Float64}, d::kite, u::StridedArray{Float64}, 
  v::StridedArray{Float64}, k::Int)
  p = d.pths[k]
  reg = p.reg

  # Affine maps to xi1, xi2 (default: xi1 from c/u, xi2 from t/v)
  hc = p.ck1 - p.ck0
  ht = p.tk1 - p.tk0
  αu = 0.5 * hc
  αv = 0.5 * ht
  βu = p.ck0 + αu
  βv = p.tk0 + αv

  # Precomputed geometry params 
  e1x = d.RP.e₁x
  e1y = d.RP.e₁y
  e2x = d.RP.e₂x
  e2y = d.RP.e₂y
  nme = d.RP.nme
  q1x = d.RP.q₁x
  q1y = d.RP.q₁y
  q2x = d.RP.q₂x
  q2y = d.RP.q₂y
  α1x = d.RP.α₁x
  α1y = d.RP.α₁y
  α2x = d.RP.α₂x
  α2y = d.RP.α₂y

  # Common scalings
  invn  = inv(nme)
  hinvn = 0.5 * invn
  αuv   = αu * αv  

  if reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

    # mapxy X precomputes
    sL1x = d.L1 * e1x * hinvn
    sL1y = d.L1 * e1y * hinvn
    sL2x = d.L2 * q1x * invn
    sL2y = d.L2 * q1y * invn

    Xx₀ = d.A + α1x + sL2x
    Xy₀ = d.B + α1y - sL2y

    # Dmap constants for this region (uses e2, q2, α1)
    cL1x = d.L1 * e2x * hinvn         
    cL1y = d.L1 * e2y * hinvn
    cL2x = d.L2 * q2x * invn
    cL2y = d.L2 * q2y * invn

    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      t2m1 = muladd(2.0, xi2, -1.0)
      Xx = muladd(t2m1,  sL1x, Xx₀)
      Xy = muladd(t2m1, -sL1y, Xy₀)

      th = muladd(xi2, Δth, th0)
      s, c = sincos(th)

      # Y and Z (mapxy)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)

      # dgx,dgy from same sincos(th)
      dgx = -d.R1 * s - 2.0 * d.P * s * c
      dgy =  d.R2 * c

      # Dmap unscaled partials (matching MATLAB after mapping u,v -> xi1,xi2)
      dux = Yx - (d.A + α1x) - t2m1 * cL1x - cL2x
      duy = Yy - (d.B + α1y) - t2m1 * cL1y - cL2y

      dvx = (d.L1 * (1.0 - xi1) * e2x) * invn + Δth * xi1 * dgx
      dvy = (d.L1 * (1.0 - xi1) * e2y) * invn + Δth * xi1 * dgy

      DJ[I] = αuv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

    sL1x = d.L1 * e1x * hinvn
    sL1y = d.L1 * e1y * hinvn
    sL2x = d.L2 * q1x * invn
    sL2y = d.L2 * q1y * invn

    Xx₀ = d.A + α1x + sL1x
    Xy₀ = d.B + α1y - sL1y

    cL1x = d.L1 * e2x * hinvn
    cL1y = d.L1 * e2y * hinvn
    cL2x = d.L2 * q2x * invn
    cL2y = d.L2 * q2y * invn

    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      one_m_xi2 = 1.0 - xi2

      Xx = muladd(one_m_xi2, sL2x, Xx₀)
      Xy = muladd(one_m_xi2, -sL2y, Xy₀)

      th = muladd(xi2, Δth, th0)
      s, c = sincos(th)

      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)

      dgx = -d.R1 * s - 2.0 * d.P * s * c
      dgy = d.R2 * c

      dux = Yx - (d.A + α1x) - cL1x - one_m_xi2 * cL2x
      duy = Yy - (d.B + α1y) - cL1y - one_m_xi2 * cL2y

      dvx = -(d.L2 * (1.0 - xi1) * q2x) * invn + Δth * xi1 * dgx
      dvy = -(d.L2 * (1.0 - xi1) * q2y) * invn + Δth * xi1 * dgy

      DJ[I] = αuv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 3
    Δth = d.tht1
    th0 = -d.tht1

    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      one_m_xi2 = 1.0 - xi2

      Xx = muladd(one_m_xi2, α1x, d.A)
      Xy = muladd(one_m_xi2, α1y, d.B)

      th = muladd(xi2, Δth, th0)
      s, c = sincos(th)

      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)

      dgx = -d.R1 * s - 2.0 * d.P * s * c
      dgy =  d.R2 * c

      dux = (Yx - d.A) - one_m_xi2 * α1x
      duy = (Yy - d.B) - one_m_xi2 * α1y

      dvx = -(1.0 - xi1) * α1x + Δth * xi1 * dgx
      dvy = -(1.0 - xi1) * α1y + Δth * xi1 * dgy

      DJ[I] = αuv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 4
    Δth = d.tht1
    th0 = 0.0

    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      Xx = muladd(xi2, α1x, d.A)
      Xy = muladd(xi2, α2y, d.B)

      th = muladd(xi2, Δth, th0)
      s, c = sincos(th)

      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)

      dgx = -d.R1 * s - 2.0 * d.P * s * c
      dgy =  d.R2 * c

      dux = (Yx - d.A) - xi2 * α2x
      duy = (Yy - d.B) - xi2 * α2y

      dvx = (1.0 - xi1) * α2x + Δth * xi1 * dgx
      dvy = (1.0 - xi1) * α2y + Δth * xi1 * dgy

      DJ[I] = αuv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

    sL1x = d.L1 * e1x * hinvn
    sL1y = d.L1 * e1y * hinvn
    sL2x = d.L2 * q1x * invn
    sL2y = d.L2 * q1y * invn

    Xx₀ = d.A + α1x + sL1x
    Xy₀ = d.B + α2y + sL1y

    cL1x = d.L1 * e1x * hinvn
    cL1y = d.L1 * e1y * hinvn
    cL2x = d.L2 * q1x * invn
    cL2y = d.L2 * q1y * invn

    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      Xx = muladd(xi2, sL2x, Xx₀)
      Xy = muladd(xi2, sL2y, Xy₀)

      th = muladd(xi2, Δth, th0)
      s, c = sincos(th)

      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)

      dgx = -d.R1 * s - 2.0 * d.P * s * c
      dgy =  d.R2 * c

      dux = Yx - (d.A + α2x) - cL1x - xi2 * cL2x
      duy = Yy - (d.B + α2y) - cL1y - xi2 * cL2y

      dvx = (d.L2 * (1.0 - xi1) * q1x) * invn + Δth * xi1 * dgx
      dvy = (d.L2 * (1.0 - xi1) * q1y) * invn + Δth * xi1 * dgy

      DJ[I] = αuv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

    sL1x = d.L1 * e1x * hinvn
    sL1y = d.L1 * e1y * hinvn
    sL2x = d.L2 * q1x * invn
    sL2y = d.L2 * q1y * invn

    Xx₀ = d.A + α1x + sL2x
    Xy₀ = d.B + α2y + sL2y

    cL1x = d.L1 * e1x * hinvn
    cL1y = d.L1 * e1y * hinvn
    cL2x = d.L2 * q1x * invn
    cL2y = d.L2 * q1y * invn

    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      t2m1 = muladd(2.0, xi2, -1.0)
      Xx = muladd(-t2m1, sL1x, Xx₀)
      Xy = muladd(-t2m1, sL1y, Xy₀)

      th = muladd(xi2, Δth, th0)
      s, c = sincos(th)

      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)

      dgx = -d.R1 * s - 2.0 * d.P * s * c
      dgy = d.R2 * c

      dux = Yx - (d.A + α2x) + t2m1 * cL1x - cL2x
      duy = Yy - (d.B + α2y) + t2m1 * cL1y - cL2y

      dvx = -(d.L1 * (1.0 - xi1) * e1x) * invn + Δth * xi1 * dgx
      dvy = -(d.L1 * (1.0 - xi1) * e1y) * invn + Δth * xi1 * dgy

      DJ[I] = αuv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

    sL1x = d.L1 * e1x * hinvn
    sL1y = d.L1 * e1y * hinvn
    sL2x = d.L2 * q1x * invn
    sL2y = d.L2 * q1y * invn

    Xx₀ = d.A + α1x - sL1x
    Xy₀ = d.B + α2y - sL1y

    cL1x = d.L1 * e1x * hinvn
    cL1y = d.L1 * e1y * hinvn
    cL2x = d.L2 * q1x * invn
    cL2y = d.L2 * q1y * invn

    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      one_m_xi2 = 1.0 - xi2
      Xx = muladd(one_m_xi2, sL2x, Xx₀)
      Xy = muladd(one_m_xi2, sL2y, Xy₀)

      th = muladd(xi2, Δth, th0)
      s, c = sincos(th)

      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)

      dgx = -d.R1 * s - 2.0 * d.P * s * c
      dgy =  d.R2 * c

      dux = Yx - (d.A + α2x) + cL1x - one_m_xi2 * cL2x
      duy = Yy - (d.B + α2y) + cL1y - one_m_xi2 * cL2y

      dvx = -(d.L2 * (1.0 - xi1) * q1x) * invn + Δth * xi1 * dgx
      dvy = -(d.L2 * (1.0 - xi1) * q1y) * invn + Δth * xi1 * dgy

      DJ[I] = αuv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 8
    Δth = d.tht2
    th0 = π - d.tht2

    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      Xx = muladd(1.0 - xi2, α1x, d.A)
      Xy = muladd(1.0 - xi2, α2y, d.B)

      th = muladd(xi2, Δth, th0)
      s, c = sincos(th)

      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)

      dgx = -d.R1 * s - 2.0 * d.P * s * c
      dgy =  d.R2 * c

      dux = (Yx - d.A) - (1.0 - xi2) * α2x
      duy = (Yy - d.B) - (1.0 - xi2) * α2y

      dvx = -(1.0 - xi1) * α2x + Δth * xi1 * dgx
      dvy = -(1.0 - xi1) * α2y + Δth * xi1 * dgy

      DJ[I] = αuv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 9
    Δth = d.tht2
    th0 = π

    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      Xx = muladd(xi2, α1x, d.A)
      Xy = muladd(xi2, α1y, d.B)

      th = muladd(xi2, Δth, th0)
      s, c = sincos(th)

      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)

      dgx = -d.R1 * s - 2.0 * d.P * s * c
      dgy =  d.R2 * c

      dux = (Yx - d.A) - xi2 * α1x
      duy = (Yy - d.B) - xi2 * α1y

      dvx = (1.0 - xi1) * α1x + Δth * xi1 * dgx
      dvy = (1.0 - xi1) * α1y + Δth * xi1 * dgy

      DJ[I] = αuv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

    sL1x = d.L1 * e1x * hinvn
    sL1y = d.L1 * e1y * hinvn
    sL2x = d.L2 * q1x * invn
    sL2y = d.L2 * q1y * invn

    Xx₀ = d.A + α1x - sL1x
    Xy₀ = d.B + α1y + sL1y

    cL1x = d.L1 * e2x * hinvn
    cL1y = d.L1 * e2y * hinvn
    cL2x = d.L2 * q2x * invn
    cL2y = d.L2 * q2y * invn

    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      Xx = muladd(xi2, sL2x, Xx₀)
      Xy = muladd(xi2, -sL2y, Xy₀)

      th = muladd(xi2, Δth, th0)
      s, c = sincos(th)

      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      Zx[I] = muladd(xi1, Yx - Xx, Xx)
      Zy[I] = muladd(xi1, Yy - Xy, Xy)

      dgx = -d.R1 * s - 2.0 * d.P * s * c
      dgy =  d.R2 * c

      dux = Yx - (d.A + α1x) + cL1x - xi2 * cL2x
      duy = Yy - (d.B + α1y) + cL1y - xi2 * cL2y

      dvx = (d.L2 * (1.0 - xi1) * q2x) * invn + Δth * xi1 * dgx
      dvy = (d.L2 * (1.0 - xi1) * q2y) * invn + Δth * xi1 * dgy

      DJ[I] = αuv * abs(dux * dvy - dvx * duy)
    end

  elseif reg == 11
    # swap parameterization if L2 >= L1:
    # xi1 from t/u, xi2 from c/v
    DJconst = d.L1 * d.L2 * αuv

    if d.L2 >= d.L1
      αu = 0.5 * ht
      αv = 0.5 * hc
      βu = p.tk0 + αu
      βv = p.ck0 + αv
    end

    sL1x = d.L1 * e2x * hinvn
    sL1y = d.L1 * e2y * hinvn
    sL2x = d.L2 * q2x * invn
    sL2y = d.L2 * q2y * invn

    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)
      t = muladd(2.0, xi2, -1.0)

      Zx[I] = d.A + α1x + xi1 * sL2x + t * sL1x
      Zy[I] = d.B + α1y + xi1 * sL2y + t * sL1y
      DJ[I] = DJconst
    end

  elseif reg == 12
    # swap parameterization if L2 >= L1:
    # xi1 from t/u, xi2 from c/v
    DJconst = d.L1 * d.L2 * αuv

    if d.L2 >= d.L1
      αu = 0.5 * ht
      αv = 0.5 * hc
      βu = p.tk0 + αu
      βv = p.ck0 + αv
    end

    sL1x = d.L1 * e1x * hinvn
    sL1y = d.L1 * e1y * hinvn
    sL2x = d.L2 * q1x * invn
    sL2y = d.L2 * q1y * invn


    @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      t = muladd(2.0, xi2, -1.0)

      Zx[I] = d.A + α2x + xi1 * sL2x + t * sL1x
      Zy[I] = d.B + α2y + xi1 * sL2y + t * sL1y
      DJ[I] = DJconst
    end

  else
    throw(ArgumentError("mapxy_Dmap! expected reg in 1–12; got reg=$reg"))

  end

  return nothing
end

function chk_map(d::kite; n::Int=32, tol::Float64=5e-14)
  #f!(F, x, y) = fill!(F, 1.0)
  #Iex =  π * d.R1 * d.R2 
  Iex = π * (
    (d.A^2 + d.B^2) * d.R1 * d.R2 -
    d.A * d.P * d.R1 * d.R2 / 2 +
    d.R1 * d.R2 * (d.R1^2 + d.R2^2 + d.P^2 / 2) / 4
  )

  return chkmap_geom(d, Iex; n=n, tol=tol)
end

"""
    jinvmap(d::kite, u::Float64, v::Float64, r::Int) -> tuple

Inverse Jacobian of the mapping at a single reference point `(u,v) ∈ [-1,1]^2`
for the given patch. Returns a 4 Tuple which are components of matrix `J⁻¹`.
Only the region of the patch is needed as inverse are only computed over them.
ht = hc = 1
Notes:
- For the rectangular region, the Jacobian is constant
  (affine map), so the inverse is a simple diagonal inverse (with the
  axis-swap handled when `L2 ≥ L1`, as in `mapxy`/`Dmap`).
"""
# Return (J11, J12, J21, J22) of the inverse Jacobian
function jinvmap(d::kite, u::Float64, v::Float64, r::Int)

  û = (u + 1) / 2
  v̂ = (v + 1) / 2

  # Precomputed geometry params 
  e1x = d.RP.e₁x
  e1y = d.RP.e₁y
  e2x = d.RP.e₂x
  e2y = d.RP.e₂y
  nme = d.RP.nme
  q1x = d.RP.q₁x
  q1y = d.RP.q₁y
  q2x = d.RP.q₂x
  q2y = d.RP.q₂y
  α1x = d.RP.α₁x
  α1y = d.RP.α₁y
  α2x = d.RP.α₂x
  α2y = d.RP.α₂y

  if r == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4

    th = muladd(v̂, Δth, th0)

    gx, gy, dgx, dgy = dergam(d, th)

    dux = gx - α1x - (2 * v̂ - 1) * d.L1 * e2x / (2 * nme) - d.L2 * q2x / nme
    duy = gy - α1y - (2 * v̂ - 1) * d.L1 * e2y / (2 * nme) - d.L2 * q2y / nme

    dvx = d.L1 * (1 - û) * e2x / nme + Δth * û * dgx
    dvy = d.L1 * (1 - û) * e2y / nme + Δth * û * dgy

  elseif r == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3

    th = muladd(v̂, Δth, th0)

    gx, gy, dgx, dgy = dergam(d, th)

    dux = gx - α1x - d.L1 * e2x / (2 * nme) - (1 - v̂) * d.L2 * q2x / nme
    duy = gy - α1y - d.L1 * e2y / (2 * nme) - (1 - v̂) * d.L2 * q2y / nme

    dvx = -d.L2 * (1 - û) * q2x / nme + Δth * û * dgx
    dvy = -d.L2 * (1 - û) * q2y / nme + Δth * û * dgy

  elseif r == 3
    Δth = d.tht1
    th0 = -d.tht1

    th = muladd(v̂, Δth, th0)

    gx, gy, dgx, dgy = dergam(d, th)

    dux = gx - (1 - v̂) * α1x
    duy = gy - (1 - v̂) * α1y

    dvx = -(1 - û) * α1x + Δth * û * dgx
    dvy = -(1 - û) * α1y + Δth * û * dgy

  elseif r == 4
    Δth = d.tht1
    th0 = 0.0

    th = muladd(v̂, Δth, th0)

    gx, gy, dgx, dgy = dergam(d, th)

    dux = gx - v̂ * α2x
    duy = gy - v̂ * α2y

    dvx = (1 - û) * α2x + Δth * û * dgx
    dvy = (1 - û) * α2y + Δth * û * dgy

  elseif r == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1

    th = muladd(v̂, Δth, th0)

    gx, gy, dgx, dgy = dergam(d, th)

    dux = gx - α2x - d.L1 * e1x / (2 * nme) - v̂ * d.L2 * q1x / nme
    duy = gy - α2y - d.L1 * e1y / (2 * nme) - v̂ * d.L2 * q1y / nme

    dvx = d.L2 * (1 - û) * q1x / nme + Δth * û * dgx
    dvy = d.L2 * (1 - û) * q1y / nme + Δth * û * dgy

  elseif r == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3

    th = muladd(v̂, Δth, th0)

    gx, gy, dgx, dgy = dergam(d, th)

    dux = gx - α2x + (2 * v̂ - 1) * d.L1 * e1x / (2 * nme) - d.L2 * q1x / nme
    duy = gy - α2y + (2 * v̂ - 1) * d.L1 * e1y / (2 * nme) - d.L2 * q1y / nme

    dvx = -d.L1 * (1 - û) * e1x / nme + Δth * û * dgx
    dvy = -d.L1 * (1 - û) * e1y / nme + Δth * û * dgy

  elseif r == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4

    th = muladd(v̂, Δth, th0)

    gx, gy, dgx, dgy = dergam(d, th)

    dux = gx - α2x + d.L1 * e1x / (2 * nme) - (1 - v̂) * d.L2 * q1x / nme
    duy = gy - α2y + d.L1 * e1y / (2 * nme) - (1 - v̂) * d.L2 * q1y / nme

    dvx = -d.L2 * (1 - û) * q1x / nme + Δth * û * dgx
    dvy = -d.L2 * (1 - û) * q1y / nme + Δth * û * dgy


  elseif r == 8
    Δth = d.tht2
    th0 = π - d.tht2

    th = muladd(v̂, Δth, th0)

    gx, gy, dgx, dgy = dergam(d, th)

    dux = gx - (1 - v̂) * α2x
    duy = gy - (1 - v̂) * α2y

    dvx = -(1 - û) * α2x + Δth * û * dgx
    dvy = -(1 - û) * α2y + Δth * û * dgy


  elseif r == 9
    Δth = d.tht2
    th0 = π

    th = muladd(v̂, Δth, th0)

    gx, gy, dgx, dgy = dergam(d, th)

    dux = gx - v̂ * α1x
    duy = gy - v̂ * α1y

    dvx = (1 - û) * α1x + Δth * û * dgx
    dvy = (1 - û) * α1y + Δth * û * dgy

  elseif r == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2

    th = muladd(v̂, Δth, th0)

    gx, gy, dgx, dgy = dergam(d, th)

    dux = gx - α1x + d.L1 * e2x / (2 * nme) - v̂ * d.L2 * q2x / nme
    duy = gy - α1y + d.L1 * e2y / (2 * nme) - v̂ * d.L2 * q2y / nme

    dvx = d.L2 * (1 - û) * q2x / nme + Δth * û * dgx
    dvy = d.L2 * (1 - û) * q2y / nme + Δth * û * dgy
  else
    throw(ArgumentError("jinvmap is defined only for regions 1–10; got reg=$r"))

  end

  detJ = dux * dvy - dvx * duy
  invdet = 2.0 / detJ

  # inverse of [dux dvx; duy dvy] is (1/detJ)[dvy -dvx; -duy dux]
  J11 = invdet * dvy
  J12 = -invdet * dvx
  J21 = -invdet * duy
  J22 = invdet * dux

  return J11, J12, J21, J22

end

"""
  mapinv(d::kite, u::Float64, v::Float64, k::Int) -> Tuple{Float64,Float64}

Inverse mapping: given a physical point `(u,v)` on patch `k`, return
the reference coordinates `[t; s]` in `[-1,1]^2` such that `mapxy(d, t, s, k) = (u,v)`.
"""
@inline function Xx(s::Float64, d::kite, r::Int)

  ex = d.RP.e₁x
  nme = d.RP.nme
  qx = d.RP.q₁x
  αx = d.RP.α₁x

  if r == 1
    return d.A + αx + (2 * s - 1) * d.L1 * ex / (2 * nme) + d.L2 * qx / nme

  elseif r == 2
    return d.A + αx + d.L1 * ex / (2 * nme) + (1 - s) * d.L2 * qx / nme

  elseif r == 3
    return d.A + (1 - s) * αx

  elseif r == 4
    return d.A + s * αx

  elseif r == 5
    return d.A + αx + d.L1 * ex / (2 * nme) + s * d.L2 * qx / nme

  elseif r == 6
    return d.A + αx - (2 * s - 1) * d.L1 * ex / (2 * nme) + d.L2 * qx / nme

  elseif r == 7
    return d.A + αx - d.L1 * ex / (2 * nme) + (1 - s) * d.L2 * qx / nme

  elseif r == 8
    return d.A + (1 - s) * αx

  elseif r == 9
    return d.A + s * αx

  elseif r == 10
    return d.A + αx - d.L1 * ex / (2 * nme) + d.L2 * s * qx / nme

  end

end

@inline function Xy(s::Float64, d::kite, r::Int)

  ey  = d.RP.e₁y
  nme = d.RP.nme
  qy  = d.RP.q₁y          
  α1y = d.RP.α₁y
  α2y = d.RP.α₂y

  if r == 1

    return d.B + α1y + (2 * s - 1) * d.L1 * (-ey) / (2 * nme) + d.L2 * (-qy) / nme

  elseif r == 2

    return d.B + α1y + d.L1 * (-ey) / (2 * nme) + (1 - s) * d.L2 * (-qy) / nme

  elseif r == 3

    return d.B + (1 - s) * α1y

  elseif r == 4

    return d.B + s * α2y

  elseif r == 5

    return d.B + α2y + d.L1 * ey / (2 * nme) + s * d.L2 * qy / nme

  elseif r == 6

    return d.B + α2y - (2 * s - 1) * d.L1 * ey / (2 * nme) + d.L2 * qy / nme

  elseif r == 7

    return d.B + α2y - d.L1 * ey / (2 * nme) + (1 - s) * d.L2 * qy / nme

  elseif r == 8

    return d.B + (1 - s) * α2y

  elseif r == 9

    return d.B + s * α1y

  elseif r == 10

    return d.B + α1y - d.L1 * (-ey) / (2 * nme) + d.L2 * s * (-qy) / nme

  end
end

@inline function Yx(s::Float64, d::kite, r::Int)

  th = if r == 1
    s * (d.tht4 - d.tht3) - d.tht4

  elseif r == 2
    s * (d.tht3 - d.tht1) - d.tht3

  elseif r == 3
    s * d.tht1 - d.tht1

  elseif r == 4
    s * d.tht1

  elseif r == 5
    s * (d.tht3 - d.tht1) + d.tht1

  elseif r == 6
    s * (d.tht4 - d.tht3) + d.tht3

  elseif r == 7
    s * (π - d.tht2 - d.tht4) + d.tht4

  elseif r == 8
    s * d.tht2 + π - d.tht2

  elseif r == 9
    s * d.tht2 + π

  elseif r == 10
    s * (π - d.tht2 - d.tht4) + π + d.tht2

  end

  si, co = sincos(th)

  return d.A + d.R1 * co - d.P * si * si

end

@inline function Yy(s::Float64, d::kite, r::Int)
  th = if r == 1
    s * (d.tht4 - d.tht3) - d.tht4

  elseif r == 2
    s * (d.tht3 - d.tht1) - d.tht3

  elseif r == 3
    s * d.tht1 - d.tht1

  elseif r == 4
    s * d.tht1

  elseif r == 5
    s * (d.tht3 - d.tht1) + d.tht1

  elseif r == 6
    s * (d.tht4 - d.tht3) + d.tht3

  elseif r == 7
    s * (π - d.tht2 - d.tht4) + d.tht4

  elseif r == 8
    s * d.tht2 + π - d.tht2

  elseif r == 9
    s * d.tht2 + π

  elseif r == 10
    s * (π - d.tht2 - d.tht4) + π + d.tht2

  end

  return d.B + d.R2 * sin(th)
end

function fill_FTable!(tbl::FTable, d::kite, r::Int)
  vmin, vmax, N = tbl.vmin, tbl.vmax, tbl.N
  h = (vmax - vmin) / (N - 1)

  P1, P2, P3 = tbl.P1, tbl.P2, tbl.P3

  ex = d.RP.e₁x
  nme = d.RP.nme
  qx = d.RP.q₁x
  αx = d.RP.α₁x
  ey  = d.RP.e₁y
  qy  = d.RP.q₁y          
  α1y = d.RP.α₁y
  α2y = d.RP.α₂y

  if r == 1
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A + αx + (2 * v̂ - 1) * d.L1 * ex / (2 * nme) + d.L2 * qx / nme
      Xy = d.B + α1y + (2 * v̂ - 1) * d.L1 * (-ey) / (2 * nme) + d.L2 * (-qy) / nme
      th = v̂ * (d.tht4 - d.tht3) - d.tht4
      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx
      P3[i] = Xy * Yx - Xx * Yy
    end
  elseif r == 2
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A + αx + d.L1 * ex / (2 * nme) + (1 - v̂) * d.L2 * qx / nme
      Xy = d.B + α1y + d.L1 * (-ey) / (2 * nme) + (1 - v̂) * d.L2 * (-qy) / nme
      th = v̂ * (d.tht3 - d.tht1) - d.tht3
      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx
      P3[i] = Xy * Yx - Xx * Yy
    end

  elseif r == 3
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A + (1 - v̂) * αx
      Xy = d.B + (1 - v̂) * α1y
      th = v̂ * d.tht1 - d.tht1
      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx
      P3[i] = Xy * Yx - Xx * Yy
    end

  elseif r == 4
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A + v̂ * αx
      Xy = d.B + v̂ * α2y
      th = v̂ * d.tht1
      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx
      P3[i] = Xy * Yx - Xx * Yy
    end
  elseif r == 5
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A + αx + d.L1 * ex / (2 * nme) + v̂ * d.L2 * qx / nme
      Xy = d.B + α2y + d.L1 * ey / (2 * nme) + v̂ * d.L2 * qy / nme
      th = v̂ * (d.tht3 - d.tht1) + d.tht1
      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx
      P3[i] = Xy * Yx - Xx * Yy
    end
  elseif r == 6
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A + αx - (2 * v̂ - 1) * d.L1 * ex / (2 * nme) + d.L2 * qx / nme
      Xy = d.B + α2y - (2 * v̂ - 1) * d.L1 * ey / (2 * nme) + d.L2 * qy / nme
      th = v̂ * (d.tht4 - d.tht3) + d.tht3
      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx
      P3[i] = Xy * Yx - Xx * Yy
    end
  elseif r == 7
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A + αx - d.L1 * ex / (2 * nme) + (1 - v̂) * d.L2 * qx / nme
      Xy = d.B + α2y - d.L1 * ey / (2 * nme) + (1 - v̂) * d.L2 * qy / nme
      th = v̂ * (π - d.tht2 - d.tht4) + d.tht4
      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx
      P3[i] = Xy * Yx - Xx * Yy
    end
  elseif r == 8
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A + (1 - v̂) * αx
      Xy = d.B + (1 - v̂) * α2y
      th = v̂ * d.tht2 + π - d.tht2
      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx
      P3[i] = Xy * Yx - Xx * Yy
    end
  elseif r == 9
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A + v̂ * αx
      Xy = d.B + v̂ * α1y
      th = v̂ * d.tht2 + π
      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx
      P3[i] = Xy * Yx - Xx * Yy
    end
  elseif r == 10
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A + αx - d.L1 * ex / (2 * nme) + d.L2 * v̂ * qx / nme
      Xy = d.B + α1y - d.L1 * (-ey) / (2 * nme) + d.L2 * v̂ * (-qy) / nme
      th = v̂ * (π - d.tht2 - d.tht4) + π + d.tht2
      s, c = sincos(th)
      Yx = d.A + d.R1 * c - d.P * s * s
      Yy = d.B + d.R2 * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx
      P3[i] = Xy * Yx - Xx * Yy
    end
  end

  tbl.reg = r
  return tbl
end

@inline function f1I(t::Float64, s::Float64,
  d::kite, u::Float64, v::Float64, r::Int)
  ŝ = (s + 1) / 2
  t̂ = (t + 1) / 2

  Xxv, Yxv = Xx(ŝ, d, r),  Yx(ŝ, d, r)

  return (1 - t̂) * Xxv + t̂ * Yxv - u
end

@inline function f2I(t::Float64, s::Float64,
  d::kite, u::Float64, v::Float64, r::Int)
    ŝ = (s + 1) / 2
    t̂ = (t + 1) / 2

    Xyv, Yyv = Xy(ŝ, d, r), Yy(ŝ, d, r)

    return (1 - t̂) * Xyv + t̂ * Yyv - v
end

@inline function JinvI(t::Float64, s::Float64,
  d::kite, u::Float64, v::Float64, r::Int)
  return jinvmap(d, t, s, r)  # returns J11,J12,J21,J22
end

@inline function f_cont(v̂::Float64, d::kite, u::Float64, v::Float64, r::Int)

  # Geometry
  ex = d.RP.e₁x
  ey = d.RP.e₁y
  nme = d.RP.nme
  qx = d.RP.q₁x
  qy = d.RP.q₁y
  αx = d.RP.α₁x
  α1y = d.RP.α₁y
  α2y = d.RP.α₂y

  if r == 1
    Xxv = d.A + αx + (2 * v̂ - 1) * d.L1 * ex / (2 * nme) + d.L2 * qx / nme
    Xyv = d.B + α1y + (2 * v̂ - 1) * d.L1 * (-ey) / (2 * nme) + d.L2 * (-qy) / nme
    th = v̂ * (d.tht4 - d.tht3) - d.tht4
  elseif r == 2
    Xxv = d.A + αx + d.L1 * ex / (2 * nme) + (1 - v̂) * d.L2 * qx / nme
    Xyv = d.B + α1y + d.L1 * (-ey) / (2 * nme) + (1 - v̂) * d.L2 * (-qy) / nme
    th = v̂ * (d.tht3 - d.tht1) - d.tht3
  elseif r == 3
    Xxv = d.A + (1 - v̂) * αx
    Xyv = d.B + (1 - v̂) * α1y
    th = v̂ * d.tht1 - d.tht1
  elseif r == 4
    Xxv = d.A + v̂ * αx
    Xyv = d.B + v̂ * α2y
    th = v̂ * d.tht1
  elseif r == 5
    Xxv = d.A + αx + d.L1 * ex / (2 * nme) + v̂ * d.L2 * qx / nme
    Xyv = d.B + α2y + d.L1 * ey / (2 * nme) + v̂ * d.L2 * qy / nme
    th = v̂ * (d.tht3 - d.tht1) + d.tht1
  elseif r == 6
    Xxv = d.A + αx - (2 * v̂ - 1) * d.L1 * ex / (2 * nme) + d.L2 * qx / nme
    Xyv = d.B + α2y - (2 * v̂ - 1) * d.L1 * ey / (2 * nme) + d.L2 * qy / nme
    th = v̂ * (d.tht4 - d.tht3) + d.tht3
  elseif r == 7
    Xxv = d.A + αx - d.L1 * ex / (2 * nme) + (1 - v̂) * d.L2 * qx / nme
    Xyv = d.B + α2y - d.L1 * ey / (2 * nme) + (1 - v̂) * d.L2 * qy / nme
    th = v̂ * (π - d.tht2 - d.tht4) + d.tht4
  elseif r == 8
    Xxv = d.A + (1 - v̂) * αx
    Xyv = d.B + (1 - v̂) * α2y
    th = v̂ * d.tht2 + π - d.tht2
  elseif r == 9
    Xxv = d.A + v̂ * αx
    Xyv = d.B + v̂ * α1y
    th = v̂ * d.tht2 + π
  elseif r == 10
    Xxv = d.A + αx - d.L1 * ex / (2 * nme) + d.L2 * v̂ * qx / nme
    Xyv = d.B + α1y - d.L1 * (-ey) / (2 * nme) + d.L2 * v̂ * (-qy) / nme
    th = v̂ * (π - d.tht2 - d.tht4) + π + d.tht2
  end

  s, c = sincos(th)

  Yxv = d.A + d.R1 * c - d.P * s * s
  Yyv = d.B + d.R2 * s

  dx = Yxv - Xxv
  dy = Yyv - Xyv

  term1 = muladd(u, dy, -v * dx)              # u*dy - v*dx
  term2 = muladd(Xyv, Yxv, -Xxv * Yyv)        # det
  return term1 + term2
end

function mapinv(tbl::FTable, d::kite, u::Float64,v::Float64, k::Int)::Tuple{Float64,Float64}

  p = d.pths[k]
  r = p.reg

  # rectangles
  if r == 11 
    # Precomputed geometry params
    e2x = d.RP.e₂x
    e2y = d.RP.e₂y
    nme = d.RP.nme        
    q2x = d.RP.q₂x
    q2y = d.RP.q₂y
    α1x = d.RP.α₁x
    α1y = d.RP.α₁y

    wx = u - d.A - α1x
    wy = v - d.B - α1y

    invn = inv(nme)

    t_q2 = (wx*q2x + wy*q2y) * (invn / d.L2)

    t_e2 = muladd(wx*e2x + wy*e2y, invn / d.L1, 0.5)

    if d.L2 >= d.L1
        Z1 = xi_inv(t_q2, p.tk0, p.tk1)
        Z2 = xi_inv(t_e2, p.ck0, p.ck1)
    else
        Z1 = xi_inv(t_q2, p.ck0, p.ck1)
        Z2 = xi_inv(t_e2, p.tk0, p.tk1)
    end

    return Z1, Z2

  elseif r == 12

    e1x = d.RP.e₁x
    e1y = d.RP.e₁y
    q1x = d.RP.q₁x
    q1y = d.RP.q₁y
    nme = d.RP.nme       
    α2x = d.RP.α₂x      
    α2y = d.RP.α₂y     

    wx = u - d.A - α2x
    wy = v - d.B - α2y

    invn = inv(nme)

    t_q1 = (wx*q1x + wy*q1y) * (invn / d.L2)

    t_e1 = muladd((wx*e1x + wy*e1y), invn / d.L1, 0.5)

    if d.L2 >= d.L1
        Z1 = xi_inv(t_q1, p.tk0, p.tk1)
        Z2 = xi_inv(t_e1, p.ck0, p.ck1)
    else
        Z1 = xi_inv(t_q1, p.ck0, p.ck1)
        Z2 = xi_inv(t_e1, p.tk0, p.tk1)
    end

    return Z1, Z2

  end

  # ----- curved patches reg = 1..10 -----
  # --- Stage 1: 2D Newton via Subroutines ---
  tN, sN = newtonR2D(f1I, f2I, JinvI,
    0.0, 0.0, 4, d, u, v, r; tol=1e-15)

  if tN !== :max

    return xi_inv((1 + tN) / 2, p.ck0, p.ck1), xi_inv((1 + sN) / 2, p.tk0, p.tk1)

  end

  #--- Stage 2: BIS method inversion ---

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
    ptconv(d::kite,  t1::Float64, t2::Float64, idx::Int, ptdest::String)

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
function ptconv(d::kite, t1::Float64, t2::Float64, idx::Int, ptdest::String)
  #@assert length(t) == 2 "t must be length-2 vector [t1,t2]"

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
    #@assert 1 ≤ k ≤ d.Npat "patch index out of range"
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
    mapinv2(d::kite, t1::Float64, t2::Float64, k2::Int, k::Int) -> Tuple{Float64,Float64}

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
function mapinv2(d::kite, t1::Float64, t2::Float64, k2::Int, k::Int)::Tuple{Float64,Float64}
    p = d.pths[k]  # Fields reg, ck0, ck1, tk0, tk1

    # Convert the given (t1,t2,k2) to the "regional" normalized coords of the *point*
    # Expecting something that returns two values in [-1,1]:
    #   tr1 ≈ u_ref, tr2 ≈ v_ref  (for the *global/regional* parameterization)
    tr1, tr2, _ = ptconv(d, t1, t2, k2, "to_reg")   

    # Map from [-1,1] -> [0,1]
    û = (tr1 + 1) / 2
    v̂ = (tr2 + 1) / 2

    # Axis swap for the center rectangle if L2 ≥ L1
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
- If patch `k` is in region 5 → returns ones with same shape as `t`.
- Otherwise uses: `d = 1 - ck1 + (ck1 - ck0) * t / 2`, then
raises elementwise to the power `s-1` if `s ≥ 0.5`, else to `s`.

Notes:
- `t` is expected to be **1 - t_actual**
- `t` may be a scalar `Float64` or any `StridedArray{Float64}`; the return has the same shape.
"""
function dfunc(d::kite, k::Int, t::Float64, s::Float64)::Float64

  p = d.pths[k]

  if p.reg == 11 || p.reg == 12
    return 1.0
  end

  hc = p.ck1 - p.ck0
  val = (1.0 - p.ck1) + hc * t / 2
  exp = s ≥ 0.5 ? (s - 1) : s
  return val^exp

end

function dfunc!(out::StridedArray{Float64}, d::kite, k::Int, t::StridedArray{Float64}, s::Float64)

  p = d.pths[k]

  if p.reg == 11 || p.reg == 12
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

"""
  diff_map!(out::Matrix{Float64}, Zx::Matrix{Float64}, Zy::Matrix{Float64}, 
  DJ::Matrix{Float64}, d::kite, u::Float64, v::Float64,
  u2::Matrix{Float64}, v2::Matrix{Float64}, du::Vector{Float64}, 
  dv::Vector{Float64}, k::Int; tol = 1e-4)

Fill `out` with ‖τ(u,v) - τ(u₂,v₂)‖ on patch `k`, where `(u,v)` are scalars and
`u2,v2` are matrices (meshgrid-like) with size `(length(du), length(dv))`.

No allocations! This map computes within itself
mapxy!(Zx, Zy, d, u2, v2, k)
  DLP!(DJ,     d, u2, v2, k)
it uses the combined function mapxy_Dmap!
"""
function diff_map!(out::Matrix{Float64},
  Zx::Matrix{Float64}, Zy::Matrix{Float64}, DJ::Matrix{Float64},
  d::kite, u::Float64, v::Float64,
  u2::Matrix{Float64}, v2::Matrix{Float64},
  du::Vector{Float64}, dv::Vector{Float64}, k::Int;
  tol::Float64=1e-5)

  nd_u = size(out, 1)
  nd_v = size(out, 2)

  p = d.pths[k]
  reg = p.reg

  # affine scalings for reference -> (ck0,ck1)/(tk0,tk1)
  hc = p.ck1 - p.ck0
  ht = p.tk1 - p.tk0
  αu = 0.5 * hc
  αv = 0.5 * ht
  βu = p.ck0 + αu
  βv = p.tk0 + αv

  # --- Precomputed geometry params
  e1x = d.RP.e₁x
  e1y = d.RP.e₁y
  e2x = d.RP.e₂x
  e2y = d.RP.e₂y
  nme = d.RP.nme
  q1x = d.RP.q₁x
  q1y = d.RP.q₁y
  q2x = d.RP.q₂x
  q2y = d.RP.q₂y
  α1x = d.RP.α₁x
  α1y = d.RP.α₁y
  α2x = d.RP.α₂x
  α2y = d.RP.α₂y

  invn = inv(nme)

  mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)
  # ----------------------------------------------
  # Special rectangular patches reg==11 or reg==12
  # ----------------------------------------------
  if reg == 11

    if d.L2 >= d.L1
      cV = d.L1 * αu * invn   # multiplies dv
      cU = d.L2 * αv * invn   # multiplies du
    else
      cV = d.L1 * αv * invn
      cU = d.L2 * αu * invn
    end

    @inbounds for j in 1:nd_v
      dvj = dv[j]
      @inbounds for i in 1:nd_u
        dui = du[i]
        Dx =  (cV * dvj) * e2x - (cU * dui) * q2x
        Dy =  (cV * dvj) * e2y - (cU * dui) * q2y
        out[i, j] = hypot(Dx, Dy)
      end
    end

    return nothing

  elseif reg == 12
    # MATLAB uses e1,q1 with same algebra
    if d.L2 >= d.L1
      cV = d.L1 * αu * invn
      cU = d.L2 * αv * invn
    else
      cV = d.L1 * αv * invn
      cU = d.L2 * αu * invn
    end

    @inbounds for j in 1:nd_v
      dvj = dv[j]
      @inbounds for i in 1:nd_u
        dui = du[i]
        Dx =  (cV * dvj) * e1x - (cU * dui) * q1x
        Dy =  (cV * dvj) * e1y - (cU * dui) * q1y
        out[i, j] = hypot(Dx, Dy)
      end
    end

    return nothing

  end

  # ------------------------
  # General case reg ∈ 1:10
  # ------------------------

  tux, tvy = mapxy(d, u, v, k)

  # mapped scalar (u,v) to (xi1,xi2)
  xi1 = muladd(αu, u, βu)
  xi2 = muladd(αv, v, βv)

  if reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4
  elseif reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3
  elseif reg == 3
    Δth = d.tht1
    th0 = -d.tht1
  elseif reg == 4
    Δth = d.tht1
    th0 = 0.0
  elseif reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1
  elseif reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3
  elseif reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4
  elseif reg == 8
    Δth = d.tht2
    th0 = π - d.tht2
  elseif reg == 9
    Δth = d.tht2
    th0 = π
  else # reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2
  end

  th = muladd(xi2, Δth, th0)

  st, ct = sincos(th)
  st2 = 2.0 * st * ct
  ct2 = cos(2.0 * th)

  # x-derivs
  gx  = d.R1 * ct - d.P * (st * st)
  dgx = -d.R1 * st - d.P * st2
  d2gx = -d.R1 * ct - 2.0 * d.P * ct2
  d3gx =  d.R1 * st + 4.0 * d.P * st2
  d4gx =  d.R1 * ct + 8.0 * d.P * ct2

  # y-derivs
  gy  = d.R2 * st
  dgy = d.R2 * ct
  d2gy = -gy
  d3gy = -dgy
  d4gy =  gy

  # powers of αv
  αv2 = αv * αv
  αv3 = αv2 * αv
  αv4 = αv2 * αv2

  if reg == 1
    dux = αu * (gx - α1x - (2.0 * xi2 - 1.0) * d.L1 * e2x * (0.5 * invn) - d.L2 * q2x * invn)
    duy = αu * (gy - α1y - (2.0 * xi2 - 1.0) * d.L1 * e2y * (0.5 * invn) - d.L2 * q2y * invn)

    dvx  = αv * ( d.L1 * (1.0 - xi1) * e2x * invn + Δth * xi1 * dgx )
    dvy  = αv * ( d.L1 * (1.0 - xi1) * e2y * invn + Δth * xi1 * dgy )

    duvx = αu * αv * ( -d.L1 * e2x * invn + Δth * dgx )
    duvy = αu * αv * ( -d.L1 * e2y * invn + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 2
    dux = αu * (gx - α1x - d.L1 * e2x * (0.5 * invn) - (1.0 - xi2) * d.L2 * q2x * invn)
    duy = αu * (gy - α1y - d.L1 * e2y * (0.5 * invn) - (1.0 - xi2) * d.L2 * q2y * invn)

    dvx  = αv * ( -d.L2 * (1.0 - xi1) * q2x * invn + Δth * xi1 * dgx )
    dvy  = αv * ( -d.L2 * (1.0 - xi1) * q2y * invn + Δth * xi1 * dgy )

    duvx = αu * αv * ( d.L2 * q2x * invn + Δth * dgx )
    duvy = αu * αv * ( d.L2 * q2y * invn + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 3
    dux = αu * (gx - (1.0 - xi2) * α1x)
    duy = αu * (gy - (1.0 - xi2) * α1y)

    dvx  = αv * ( -(1.0 - xi1) * α1x + Δth * xi1 * dgx )
    dvy  = αv * ( -(1.0 - xi1) * α1y + Δth * xi1 * dgy )

    duvx = αu * αv * ( α1x + Δth * dgx )
    duvy = αu * αv * ( α1y + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 4
    dux = αu * (gx - xi2 * α2x)
    duy = αu * (gy - xi2 * α2y)

    dvx  = αv * ( (1.0 - xi1) * α2x + Δth * xi1 * dgx )
    dvy  = αv * ( (1.0 - xi1) * α2y + Δth * xi1 * dgy )

    duvx = αu * αv * ( -α2x + Δth * dgx )
    duvy = αu * αv * ( -α2y + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 5
    dux = αu * (gx - α2x - d.L1 * e1x * (0.5 * invn) - xi2 * d.L2 * q1x * invn)
    duy = αu * (gy - α2y - d.L1 * e1y * (0.5 * invn) - xi2 * d.L2 * q1y * invn)

    dvx  = αv * ( d.L2 * (1.0 - xi1) * q1x * invn + Δth * xi1 * dgx )
    dvy  = αv * ( d.L2 * (1.0 - xi1) * q1y * invn + Δth * xi1 * dgy )

    duvx = αu * αv * ( -d.L2 * q1x * invn + Δth * dgx )
    duvy = αu * αv * ( -d.L2 * q1y * invn + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 6
    dux = αu * (gx - α2x + (2.0 * xi2 - 1.0) * d.L1 * e1x * (0.5 * invn) - d.L2 * q1x * invn)
    duy = αu * (gy - α2y + (2.0 * xi2 - 1.0) * d.L1 * e1y * (0.5 * invn) - d.L2 * q1y * invn)

    dvx  = αv * ( -d.L1 * (1.0 - xi1) * e1x * invn + Δth * xi1 * dgx )
    dvy  = αv * ( -d.L1 * (1.0 - xi1) * e1y * invn + Δth * xi1 * dgy )

    duvx = αu * αv * ( d.L1 * e1x * invn + Δth * dgx )
    duvy = αu * αv * ( d.L1 * e1y * invn + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 7
    dux = αu * (gx - α2x + d.L1 * e1x * (0.5 * invn) - (1.0 - xi2) * d.L2 * q1x * invn)
    duy = αu * (gy - α2y + d.L1 * e1y * (0.5 * invn) - (1.0 - xi2) * d.L2 * q1y * invn)

    dvx  = αv * ( -d.L2 * (1.0 - xi1) * q1x * invn + Δth * xi1 * dgx )
    dvy  = αv * ( -d.L2 * (1.0 - xi1) * q1y * invn + Δth * xi1 * dgy )

    duvx = αu * αv * ( d.L2 * q1x * invn + Δth * dgx )
    duvy = αu * αv * ( d.L2 * q1y * invn + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 8
    dux = αu * (gx - (1.0 - xi2) * α2x)
    duy = αu * (gy - (1.0 - xi2) * α2y)

    dvx  = αv * ( -(1.0 - xi1) * α2x + Δth * xi1 * dgx )
    dvy  = αv * ( -(1.0 - xi1) * α2y + Δth * xi1 * dgy )

    duvx = αu * αv * ( α2x + Δth * dgx )
    duvy = αu * αv * ( α2y + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 9
    dux = αu * (gx - xi2 * α1x)
    duy = αu * (gy - xi2 * α1y)

    dvx  = αv * ( (1.0 - xi1) * α1x + Δth * xi1 * dgx )
    dvy  = αv * ( (1.0 - xi1) * α1y + Δth * xi1 * dgy )

    duvx = αu * αv * ( -α1x + Δth * dgx )
    duvy = αu * αv * ( -α1y + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  else # reg == 10
    dux = αu * (gx - α1x + d.L1 * e2x * (0.5 * invn) - xi2 * d.L2 * q2x * invn)
    duy = αu * (gy - α1y + d.L1 * e2y * (0.5 * invn) - xi2 * d.L2 * q2y * invn)

    dvx  = αv * ( d.L2 * (1.0 - xi1) * q2x * invn + Δth * xi1 * dgx )
    dvy  = αv * ( d.L2 * (1.0 - xi1) * q2y * invn + Δth * xi1 * dgy )

    duvx = αu * αv * ( -d.L2 * q2x * invn + Δth * dgx )
    duvy = αu * αv * ( -d.L2 * q2y * invn + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)
  end

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

end

"""
  diff_rmap!(out::Matrix{Float64}, x::Matrix{Float64}, Zy::Matrix{Float64}, 
  DJ::StridedArray{Float64}, d::kite, u::Float64, v::Float64,
  u2::Matrix{Float64}, v2::Matrix{Float64}, r::Matrix{Float64},
  du::Vector{Float64}, dv::Vector{Float64}, k::Int; tol = 1e-4)

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
  d::kite, u::Float64, v::Float64,
  u2::Matrix{Float64}, v2::Matrix{Float64}, r::Matrix{Float64},
  du::Vector{Float64}, dv::Vector{Float64}, k::Int;
  tol::Float64=1e-4)

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
  reg = p.reg
  # affine scalings for reference -> (ck0,ck1)/(tk0,tk1)
  hc = p.ck1 - p.ck0
  ht = p.tk1 - p.tk0
  αu = 0.5 * hc
  αv = 0.5 * ht
  βu = p.ck0 + αu
  βv = p.tk0 + αv

  # --- Precomputed geometry params
  e1x = d.RP.e₁x
  e1y = d.RP.e₁y
  e2x = d.RP.e₂x
  e2y = d.RP.e₂y
  nme = d.RP.nme
  q1x = d.RP.q₁x
  q1y = d.RP.q₁y
  q2x = d.RP.q₂x
  q2y = d.RP.q₂y
  α1x = d.RP.α₁x
  α1y = d.RP.α₁y
  α2x = d.RP.α₂x
  α2y = d.RP.α₂y

  invn = inv(nme)

  if reg == 11

    fill!(DJ, d.L1 * d.L2 * αu * αv)

    if d.L2 >= d.L1
      cV = d.L1 * αu * invn   # multiplies dv
      cU = d.L2 * αv * invn   # multiplies du
    else
      cV = d.L1 * αv * invn
      cU = d.L2 * αu * invn
    end

    @inbounds for i in 1:nt
      dui = du[i]
      dvi = dv[i]
      Dx = (cV * dvi) * e2x - (cU * dui) * q2x
      Dy = (cV * dvi) * e2y - (cU * dui) * q2y
      hD = hypot(Dx, Dy)
      for j in 1:nr
        out[i, j] = hD
      end
    end

    return nothing

  elseif reg == 12

    fill!(DJ, d.L1 * d.L2 * αu * αv)

    # MATLAB uses e1,q1 with same algebra
    if d.L2 >= d.L1
      cV = d.L1 * αu * invn
      cU = d.L2 * αv * invn
    else
      cV = d.L1 * αv * invn
      cU = d.L2 * αu * invn
    end

    @inbounds for i in 1:nt
      dui = du[i]
      dvi = dv[i]
      Dx = (cV * dvi) * e1x - (cU * dui) * q1x
      Dy = (cV * dvi) * e1y - (cU * dui) * q1y
      hD = hypot(Dx, Dy)
      for j in 1:nr
        out[i, j] = hD
      end
    end

    return nothing

  end

  # ------------------------
  # General case reg ∈ 1:10
  # ------------------------

  mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)

  tux, tvy = mapxy(d, u, v, k)

  # mapped scalar (u,v) to (xi1,xi2)
  xi1 = muladd(αu, u, βu)
  xi2 = muladd(αv, v, βv)

  if reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4
  elseif reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3
  elseif reg == 3
    Δth = d.tht1
    th0 = -d.tht1
  elseif reg == 4
    Δth = d.tht1
    th0 = 0.0
  elseif reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1
  elseif reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3
  elseif reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4
  elseif reg == 8
    Δth = d.tht2
    th0 = π - d.tht2
  elseif reg == 9
    Δth = d.tht2
    th0 = π
  else # reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2
  end

  th = muladd(xi2, Δth, th0)

  st, ct = sincos(th)
  st2 = 2.0 * st * ct
  ct2 = cos(2.0 * th)

  # x-derivs
  gx  = d.R1 * ct - d.P * (st * st)
  dgx = -d.R1 * st - d.P * st2
  d2gx = -d.R1 * ct - 2.0 * d.P * ct2
  d3gx =  d.R1 * st + 4.0 * d.P * st2
  d4gx =  d.R1 * ct + 8.0 * d.P * ct2

  # y-derivs
  gy  = d.R2 * st
  dgy = d.R2 * ct
  d2gy = -gy
  d3gy = -dgy
  d4gy =  gy

  # powers of αv
  αv2 = αv * αv
  αv3 = αv2 * αv
  αv4 = αv2 * αv2

  if reg == 1
    dux = αu * (gx - α1x - (2.0 * xi2 - 1.0) * d.L1 * e2x * (0.5 * invn) - d.L2 * q2x * invn)
    duy = αu * (gy - α1y - (2.0 * xi2 - 1.0) * d.L1 * e2y * (0.5 * invn) - d.L2 * q2y * invn)

    dvx  = αv * ( d.L1 * (1.0 - xi1) * e2x * invn + Δth * xi1 * dgx )
    dvy  = αv * ( d.L1 * (1.0 - xi1) * e2y * invn + Δth * xi1 * dgy )

    duvx = αu * αv * ( -d.L1 * e2x * invn + Δth * dgx )
    duvy = αu * αv * ( -d.L1 * e2y * invn + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 2
    dux = αu * (gx - α1x - d.L1 * e2x * (0.5 * invn) - (1.0 - xi2) * d.L2 * q2x * invn)
    duy = αu * (gy - α1y - d.L1 * e2y * (0.5 * invn) - (1.0 - xi2) * d.L2 * q2y * invn)

    dvx  = αv * ( -d.L2 * (1.0 - xi1) * q2x * invn + Δth * xi1 * dgx )
    dvy  = αv * ( -d.L2 * (1.0 - xi1) * q2y * invn + Δth * xi1 * dgy )

    duvx = αu * αv * ( d.L2 * q2x * invn + Δth * dgx )
    duvy = αu * αv * ( d.L2 * q2y * invn + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 3
    dux = αu * (gx - (1.0 - xi2) * α1x)
    duy = αu * (gy - (1.0 - xi2) * α1y)

    dvx  = αv * ( -(1.0 - xi1) * α1x + Δth * xi1 * dgx )
    dvy  = αv * ( -(1.0 - xi1) * α1y + Δth * xi1 * dgy )

    duvx = αu * αv * ( α1x + Δth * dgx )
    duvy = αu * αv * ( α1y + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 4
    dux = αu * (gx - xi2 * α2x)
    duy = αu * (gy - xi2 * α2y)

    dvx  = αv * ( (1.0 - xi1) * α2x + Δth * xi1 * dgx )
    dvy  = αv * ( (1.0 - xi1) * α2y + Δth * xi1 * dgy )

    duvx = αu * αv * ( -α2x + Δth * dgx )
    duvy = αu * αv * ( -α2y + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 5
    dux = αu * (gx - α2x - d.L1 * e1x * (0.5 * invn) - xi2 * d.L2 * q1x * invn)
    duy = αu * (gy - α2y - d.L1 * e1y * (0.5 * invn) - xi2 * d.L2 * q1y * invn)

    dvx  = αv * ( d.L2 * (1.0 - xi1) * q1x * invn + Δth * xi1 * dgx )
    dvy  = αv * ( d.L2 * (1.0 - xi1) * q1y * invn + Δth * xi1 * dgy )

    duvx = αu * αv * ( -d.L2 * q1x * invn + Δth * dgx )
    duvy = αu * αv * ( -d.L2 * q1y * invn + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 6
    dux = αu * (gx - α2x + (2.0 * xi2 - 1.0) * d.L1 * e1x * (0.5 * invn) - d.L2 * q1x * invn)
    duy = αu * (gy - α2y + (2.0 * xi2 - 1.0) * d.L1 * e1y * (0.5 * invn) - d.L2 * q1y * invn)

    dvx  = αv * ( -d.L1 * (1.0 - xi1) * e1x * invn + Δth * xi1 * dgx )
    dvy  = αv * ( -d.L1 * (1.0 - xi1) * e1y * invn + Δth * xi1 * dgy )

    duvx = αu * αv * ( d.L1 * e1x * invn + Δth * dgx )
    duvy = αu * αv * ( d.L1 * e1y * invn + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 7
    dux = αu * (gx - α2x + d.L1 * e1x * (0.5 * invn) - (1.0 - xi2) * d.L2 * q1x * invn)
    duy = αu * (gy - α2y + d.L1 * e1y * (0.5 * invn) - (1.0 - xi2) * d.L2 * q1y * invn)

    dvx  = αv * ( -d.L2 * (1.0 - xi1) * q1x * invn + Δth * xi1 * dgx )
    dvy  = αv * ( -d.L2 * (1.0 - xi1) * q1y * invn + Δth * xi1 * dgy )

    duvx = αu * αv * ( d.L2 * q1x * invn + Δth * dgx )
    duvy = αu * αv * ( d.L2 * q1y * invn + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 8
    dux = αu * (gx - (1.0 - xi2) * α2x)
    duy = αu * (gy - (1.0 - xi2) * α2y)

    dvx  = αv * ( -(1.0 - xi1) * α2x + Δth * xi1 * dgx )
    dvy  = αv * ( -(1.0 - xi1) * α2y + Δth * xi1 * dgy )

    duvx = αu * αv * ( α2x + Δth * dgx )
    duvy = αu * αv * ( α2y + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  elseif reg == 9
    dux = αu * (gx - xi2 * α1x)
    duy = αu * (gy - xi2 * α1y)

    dvx  = αv * ( (1.0 - xi1) * α1x + Δth * xi1 * dgx )
    dvy  = αv * ( (1.0 - xi1) * α1y + Δth * xi1 * dgy )

    duvx = αu * αv * ( -α1x + Δth * dgx )
    duvy = αu * αv * ( -α1y + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)

  else # reg == 10
    dux = αu * (gx - α1x + d.L1 * e2x * (0.5 * invn) - xi2 * d.L2 * q2x * invn)
    duy = αu * (gy - α1y + d.L1 * e2y * (0.5 * invn) - xi2 * d.L2 * q2y * invn)

    dvx  = αv * ( d.L2 * (1.0 - xi1) * q2x * invn + Δth * xi1 * dgx )
    dvy  = αv * ( d.L2 * (1.0 - xi1) * q2y * invn + Δth * xi1 * dgy )

    duvx = αu * αv * ( -d.L2 * q2x * invn + Δth * dgx )
    duvy = αu * αv * ( -d.L2 * q2y * invn + Δth * dgy )

    dv2x  = αv2 * (Δth^2) * (xi1 * d2gx)
    dv2y  = αv2 * (Δth^2) * (xi1 * d2gy)
    duv2x = αu * αv2 * (Δth^2) * d2gx
    duv2y = αu * αv2 * (Δth^2) * d2gy

    dv3x  = αv3 * (Δth^3) * (xi1 * d3gx)
    dv3y  = αv3 * (Δth^3) * (xi1 * d3gy)
    duv3x = αu * αv3 * (Δth^3) * d3gx
    duv3y = αu * αv3 * (Δth^3) * d3gy

    dv4x  = αv4 * (Δth^4) * (xi1 * d4gx)
    dv4y  = αv4 * (Δth^4) * (xi1 * d4gy)
  end

  @inbounds for i in 1:nt
    dvi = dv[i]
    dui = du[i]
    @inbounds for j in 1:nr
      uu = u2[i, j]
      vv = v2[i, j]
      if (abs(u - uu) < tol) && (abs(v - vv) < tol)
        r1 = dvi * r[i, j]
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

end


"""
    Dwall(d::kite, u::Float64, v::Float64, k::Int) -> dvx, dvy

Derivative of the wall curve w.r.t. `v` at the singular side `u = ±1`
for the `k`-th (non-quadrilateral) patch.
Returns a Tuple `dvx, dvy`.
"""
function Dwall(d::kite, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64}

  p   = d.pths[k]
  reg = p.reg

  if reg == 11 || reg == 12
    throw(ArgumentError("Dwall is for non-quadrilateral patches; got reg=$reg"))
  end

  # Affine maps: xi1 = αu*u + βu, xi2 = αv*v + βv
  hc = p.ck1 - p.ck0
  ht = p.tk1 - p.tk0
  αu = 0.5 * hc
  αv = 0.5 * ht
  βu = p.ck0 + αu
  βv = p.tk0 + αv

  xi1 = muladd(αu, u, βu)
  xi2 = muladd(αv, v, βv)

  # Precomputed geometry params
  e1x = d.RP.e₁x
  e1y = d.RP.e₁y
  e2x = d.RP.e₂x
  e2y = d.RP.e₂y
  nme = d.RP.nme
  q1x = d.RP.q₁x
  q1y = d.RP.q₁y
  q2x = d.RP.q₂x
  q2y = d.RP.q₂y
  α1x = d.RP.α₁x
  α1y = d.RP.α₁y
  α2x = d.RP.α₂x
  α2y = d.RP.α₂y

  invn = inv(nme)

  if reg == 1
    Δth = d.tht4 - d.tht3
    th0 = -d.tht4
    th  = muladd(xi2, Δth, th0)

    _, _, dgx, dgy = dergam(d, th)

    dvx = d.L1 * (1.0 - xi1) * e2x * invn + (Δth * xi1) * dgx
    dvy = d.L1 * (1.0 - xi1) * e2y * invn + (Δth * xi1) * dgy

  elseif reg == 2
    Δth = d.tht3 - d.tht1
    th0 = -d.tht3
    th  = muladd(xi2, Δth, th0)

    _, _, dgx, dgy = dergam(d, th)

    dvx = -d.L2 * (1.0 - xi1) * q2x * invn + (Δth * xi1) * dgx
    dvy = -d.L2 * (1.0 - xi1) * q2y * invn + (Δth * xi1) * dgy

  elseif reg == 3
    Δth = d.tht1
    th0 = -d.tht1
    th  = muladd(xi2, Δth, th0)

    _, _, dgx, dgy = dergam(d, th)

    dvx = -(1.0 - xi1) * α1x + (Δth * xi1) * dgx
    dvy = -(1.0 - xi1) * α1y + (Δth * xi1) * dgy

  elseif reg == 4
    Δth = d.tht1
    th0 = 0.0
    th  = muladd(xi2, Δth, th0)

    _, _, dgx, dgy = dergam(d, th)

    dvx = (1.0 - xi1) * α2x + (Δth * xi1) * dgx
    dvy = (1.0 - xi1) * α2y + (Δth * xi1) * dgy

  elseif reg == 5
    Δth = d.tht3 - d.tht1
    th0 = d.tht1
    th  = muladd(xi2, Δth, th0)

    _, _, dgx, dgy = dergam(d, th)

    dvx = d.L2 * (1.0 - xi1) * q1x * invn + (Δth * xi1) * dgx
    dvy = d.L2 * (1.0 - xi1) * q1y * invn + (Δth * xi1) * dgy

  elseif reg == 6
    Δth = d.tht4 - d.tht3
    th0 = d.tht3
    th  = muladd(xi2, Δth, th0)

    _, _, dgx, dgy = dergam(d, th)

    dvx = -d.L1 * (1.0 - xi1) * e1x * invn + (Δth * xi1) * dgx
    dvy = -d.L1 * (1.0 - xi1) * e1y * invn + (Δth * xi1) * dgy

  elseif reg == 7
    Δth = π - d.tht2 - d.tht4
    th0 = d.tht4
    th  = muladd(xi2, Δth, th0)

    _, _, dgx, dgy = dergam(d, th)

    dvx = -d.L2 * (1.0 - xi1) * q1x * invn + (Δth * xi1) * dgx
    dvy = -d.L2 * (1.0 - xi1) * q1y * invn + (Δth * xi1) * dgy

  elseif reg == 8
    Δth = d.tht2
    th0 = π - d.tht2
    th  = muladd(xi2, Δth, th0)

    _, _, dgx, dgy = dergam(d, th)

    dvx = -(1.0 - xi1) * α2x + (Δth * xi1) * dgx
    dvy = -(1.0 - xi1) * α2y + (Δth * xi1) * dgy

  elseif reg == 9
    Δth = d.tht2
    th0 = π
    th  = muladd(xi2, Δth, th0)

    _, _, dgx, dgy = dergam(d, th)

    dvx = (1.0 - xi1) * α1x + (Δth * xi1) * dgx
    dvy = (1.0 - xi1) * α1y + (Δth * xi1) * dgy

  elseif reg == 10
    Δth = π - d.tht2 - d.tht4
    th0 = π + d.tht2
    th  = muladd(xi2, Δth, th0)

    _, _, dgx, dgy = dergam(d, th)

    dvx = d.L2 * (1.0 - xi1) * q2x * invn + (Δth * xi1) * dgx
    dvy = d.L2 * (1.0 - xi1) * q2y * invn + (Δth * xi1) * dgy

  else
    throw(ArgumentError("Dwall is defined only for regions 1–10; got reg=$reg"))
  end

  return αv * dvx, αv * dvy
end

"""
    boundquad!(P::SubArray{Float64}, d::kite, k::Int) -> nothing

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
function boundquad!(P::SubArray{Float64}, d::kite, k::Int)

  p = d.pths[k]

  # Four corners (tuples immediately unpacked into scalars)
  L1p1x, L1p1y = mapm1(d, -1.0, k)   # (u=-1, v=-1)
  L1p2x, L1p2y = mapp1(d, -1.0, k)   # (u=+1, v=-1)
  L2p1x, L2p1y = mapm1(d, 1.0, k)   # (u=-1, v=+1)
  L2p2x, L2p2y = mapp1(d, 1.0, k)   # (u=+1, v=+1)

  if p.reg == 11 || p.reg == 12
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
    refine!(d::kite, Nc::Int, Nt::Int, K::AbstractVector{<:Int})

Subdivide each patch in `K` into `Nc * Nt` subpatches (split in `c` and `t`),
update `d.pths`, `d.Npat`, recompute `d.kd`, and rebuild `d.Qpts` / `d.Qptsbd`.

Notes
- Patches are re-sorted lexicographically by `(reg, ck0, ck1, tk0, tk1)`
- `Qpts` uses `boundquad(d,k)`; `Qptsbd` uses `boundquadbd(d,k)` for patches in `d.kd`.
"""
function refine!(d::kite, Nc::Int, Nt::Int, K::Vector{Int})
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
    d.kd = [k for k in 1:d.Npat if (d.pths[k].reg != 11 && d.pths[k].reg != 12) && d.pths[k].ck1 == 1.0]

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

    return nothing
end


function Base.show(io::IO, d::kite)

  println(io, "kite with properties:")
  println(io, "  (A , B )  = (", d.A,",", d.B,")")
  println(io, "  (R₁, R₂)  = (", d.R1,",", d.R2,")")
  println(io, "  (L₁, L₂)  = (", d.L1,",", d.L2,")")
  println(io, "  (θ₁, θ₂)  = (", d.tht1,",", d.tht2,")")

  println(io, "  (θ₃, θ₄, L₃)  = (", d.tht3,",", d.tht4,",", d.L3,")")
  println(io, "  Distortion    = ", d.P)
  println(io, "  No. of holes  = ", d.nh)

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
