"""
Mutable struct ellipse

A  — x-coordinate of ellipse center (Float64)

B  — y-coordinate of ellipse center (Float64)

R₁ — Radius of the ellipse in x axis (Float64)

R₂ — Radius of the ellipse in y axis (Float64)

Cθ — Cos of angle of rotation of the ellipse. 0 ≤ θ < π. (Float64)

Sθ — Sin of angle of rotation of the ellipse. 0 ≤ θ < π. (Float64)

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
mutable struct ellipse <: abstractdomain

  A::Float64
  B::Float64
  R1::Float64
  R2::Float64
  L1::Float64
  L2::Float64
  nh::Int
  kd::Vector{Int}
  Cθ::Float64
  Sθ::Float64
  Npat::Int
  pths::Vector{Patch}   # 5 * Npat
  Qpts::Matrix{Float64}   # 8 * Npat
  Qptsbd::Matrix{Float64} # 8 * Npat

  function ellipse(; b,
    a=nothing, A=nothing, B=nothing,
    R1=nothing, R2=nothing, L1=nothing, 
    L2=nothing, tht=nothing, ck=nothing, 
    tk=nothing)
    #Construct an intsance for this structure

    """
    Required:
      - b :: Vector{Int} of length 5 (Col), positive entries
    
    Optional (all keywords):
      - a  :: Vector{Int} length 5 (defaults to ceil.(2*b/3))
      - A,B:: Center of the ellipse (defaults 0.0)
      - R1 :: Float64 (default 1.0)
      - R2 :: Float64 (default 1.5)
      - L1 :: Float64 (default 0.8*R2)
      - L2 :: Float64 (default 0.8*R1)
      - θ  :: Float64 (default 0.0)
      - ck :: vector{vector{Float64}} (default equispaced per a)
      - tk :: vector{vector{Float64}} (default equispaced per b)
    """
    @assert isa(b, AbstractVector{Int}) && length(b) == 5

    #The ellipse has obviously no holes
    nh = 0

    a = isnothing(a) ? ceil.(Int, 2 .* b ./ 3) : collect(Int, a)

    #By default, ellipse at origin
    A = something(A, 0.0)
    B = something(B, 0.0)

    R1 = something(R1, 1//1)
    R2 = something(R2, 3//2)

    L1 = something(L1, 4//5 * R2)
    L2 = something(L2, 4//5 * R1)

    tht= something(tht, 0.0)

    Sθ, Cθ = sincos(tht)

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

    d = new(A, B, R1, R2, L1, L2, nh, kd, Cθ, Sθ, Npat, pths, Qpts, Qptsbd)

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

function mapx(d::ellipse, u::Float64, v::Float64, k::Int)::Float64

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

      d.A + d.L2 * d.Cθ / 2 - d.L1 * d.Sθ * (2*xi2 - 1) / 2

    elseif p.reg == 3

      d.A - d.L2 * d.Cθ / 2 + d.L1 * d.Sθ * (2*xi2 - 1) / 2

    elseif p.reg == 2

      d.A + d.L2 * d.Cθ * (1 - 2*xi2) / 2 - d.L1 * d.Sθ / 2

    else # p.reg == 4

      d.A - d.L2 * d.Cθ * (1 - 2*xi2) / 2 + d.L1 * d.Sθ / 2

    end

    th = π * (xi2 / 2 + (2 * p.reg - 3) / 4)

    s, c = sincos(th)

    Y = d.A + d.R1 * d.Cθ * c - d.R2 * d.Sθ * s

    return (1 - xi1) * X + xi1 * Y

  else

    if d.L2 >= d.L1 
      
      xi1 = ht * u / 2 + (p.tk1 + p.tk0) / 2

      xi2 = hc * v / 2 + (p.ck1 + p.ck0) / 2

    end

    return d.A + d.L2*(xi1 - 0.5)*d.Cθ - d.L1*(xi2 - 0.5)*d.Sθ

  end

end

function mapy(d::ellipse, u::Float64, v::Float64, k::Int)::Float64

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

      d.B + d.L2 * d.Sθ / 2 + d.L1 * d.Cθ * (2*xi2 - 1) / 2

    elseif p.reg == 3

      d.B - d.L2 * d.Sθ / 2 - d.L1 * d.Cθ * (2*xi2 - 1) / 2

    elseif p.reg == 2

      d.B + d.L2 * d.Sθ * (1 - 2*xi2) / 2 + d.L1 * d.Cθ / 2

    else # p.reg == 4

      d.B - d.L2 * d.Sθ * (1 - 2*xi2) / 2 - d.L1 * d.Cθ / 2

    end

    th = π * (xi2 / 2 + (2 * p.reg - 3) / 4)

    s, c = sincos(th)

    Y = d.B + d.R1 * d.Sθ * c + d.R2 * d.Cθ * s

    return (1 - xi1) * X + xi1 * Y

  else

    if d.L2 >= d.L1 
      
      xi1 = ht * u / 2 + (p.tk1 + p.tk0) / 2

      xi2 = hc * v / 2 + (p.ck1 + p.ck0) / 2

    end

    return d.B + d.L2*(xi1 - 0.5)*d.Sθ + d.L1*(xi2 - 0.5)*d.Cθ

  end

end

function mapxy(d::ellipse, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64}

  p = d.pths[k]

  # Difference in c values
  hc = p.ck1 - p.ck0

  xi1 = hc * u / 2 + (p.ck1 + p.ck0) / 2

  # Difference in t values
  ht = p.tk1 - p.tk0

  xi2 = ht * v / 2 + (p.tk1 + p.tk0) / 2

  if p.reg != 5

    if p.reg == 1

      Xx = d.A + d.L2 * d.Cθ / 2 - d.L1 * d.Sθ * (2 * xi2 - 1) / 2
      Xy = d.B + d.L2 * d.Sθ / 2 + d.L1 * d.Cθ * (2 * xi2 - 1) / 2

    elseif p.reg == 2

      Xx = d.A + d.L2 * d.Cθ * (1 - 2 * xi2) / 2 - d.L1 * d.Sθ / 2
      Xy = d.B + d.L2 * d.Sθ * (1 - 2 * xi2) / 2 + d.L1 * d.Cθ / 2

    elseif p.reg == 3

      Xx = d.A - d.L2 * d.Cθ / 2 + d.L1 * d.Sθ * (2 * xi2 - 1) / 2
      Xy = d.B - d.L2 * d.Sθ / 2 - d.L1 * d.Cθ * (2 * xi2 - 1) / 2

    else # p.reg == 4

      Xx = d.A - d.L2 * d.Cθ * (1 - 2 * xi2) / 2 + d.L1 * d.Sθ / 2
      Xy = d.B - d.L2 * d.Sθ * (1 - 2 * xi2) / 2 - d.L1 * d.Cθ / 2

    end

    th = π * (xi2 / 2 + (2 * p.reg - 3) / 4)

    s, c = sincos(th)

    Yx = d.A + d.R1 * d.Cθ * c - d.R2 * d.Sθ * s
    Yy = d.B + d.R1 * d.Sθ * c + d.R2 * d.Cθ * s

    Zx = (1 - xi1) * Xx + xi1 * Yx
    Zy = (1 - xi1) * Xy + xi1 * Yy

  else

    if d.L2 >= d.L1

      xi1 = ht * u / 2 + (p.tk1 + p.tk0) / 2

      xi2 = hc * v / 2 + (p.ck1 + p.ck0) / 2

    end

    Zx = d.A + d.L2 * (xi1 - 0.5) * d.Cθ - d.L1 * (xi2 - 0.5) * d.Sθ
    Zy = d.B + d.L2 * (xi1 - 0.5) * d.Sθ + d.L1 * (xi2 - 0.5) * d.Cθ

  end

  return Zx, Zy
end

function mapxy!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
  d::ellipse, u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)
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
    c₁₁ = d.R1 * d.Cθ
    c₁₂ = d.R1 * d.Sθ
    c₂₁ = d.R2 * d.Cθ
    c₂₂ = -d.R2 * d.Sθ
    # Corner line (Xx,Xy) for this region
    # Iterate once over a single index range I that is 
    # valid for all four arrays.
    @inbounds for I in eachindex(u, v, Zx, Zy)
      xi1 = muladd(αu, u[I], βu)
      xi2 = muladd(αv, v[I], βv)

      # fixed edge point for given xi2
      Xx = if p.reg == 1
        d.A + d.L2 * d.Cθ / 2 - d.L1 * d.Sθ * (2*xi2 - 1) / 2
      elseif p.reg == 2
        d.A + d.L2 * d.Cθ * (1 - 2*xi2) / 2 - d.L1 * d.Sθ / 2
      elseif p.reg == 3
        d.A - d.L2 * d.Cθ / 2 + d.L1 * d.Sθ * (2*xi2 - 1) / 2
      else # p.reg == 4
        d.A - d.L2 * d.Cθ * (1 - 2*xi2) / 2 + d.L1 * d.Sθ / 2
      end
    
      Xy = if p.reg == 1
        d.B + d.L2 * d.Sθ / 2 + d.L1 * d.Cθ * (2*xi2 - 1) / 2
      elseif p.reg == 2
        d.B + d.L2 * d.Sθ * (1 - 2*xi2) / 2 + d.L1 * d.Cθ / 2
      elseif p.reg == 3
        d.B - d.L2 * d.Sθ / 2 - d.L1 * d.Cθ * (2*xi2 - 1) / 2
      else # p.reg == 4
        d.B - d.L2 * d.Sθ * (1 - 2*xi2) / 2 - d.L1 * d.Cθ / 2
      end

      # circular arc point for given xi2
      t = π *(xi2 / 2 + (2 * p.reg - 3) / 4)

      s, c = sincos(t)

      #d.A + d.R1 * d.Cθ * c - d.R2 * d.Sθ * s
      Yx = muladd(c₂₂, s, muladd(c₁₁, c, d.A)) 
      #d.B + d.R1 * d.Sθ * c + d.R2 * d.Cθ * s
      Yy = muladd(c₂₁, s, muladd(c₁₂, c, d.B))  

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
        Zx[I] = d.A + d.L2*(xi1 - 0.5)*d.Cθ - d.L1*(xi2 - 0.5)*d.Sθ
        Zy[I] = d.B + d.L2*(xi1 - 0.5)*d.Sθ + d.L1*(xi2 - 0.5)*d.Cθ
      end

    else
      αu = hc / 2
      αv = ht / 2
      βu = p.ck0 + αu
      βv = p.tk0 + αv
      @inbounds  for I in eachindex(u, v, Zx, Zy)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)
        Zx[I] = d.A + d.L2*(xi1 - 0.5)*d.Cθ - d.L1*(xi2 - 0.5)*d.Sθ
        Zy[I] = d.B + d.L2*(xi1 - 0.5)*d.Sθ + d.L1*(xi2 - 0.5)*d.Cθ
      end

    end
  end

  return nothing

end

function draw(d::ellipse, flag=nothing)
  # number of samples for the patch boundary
  L = 33

  #This is treated as a column vector
  t = collect(range(-1, 1, L))

  # vertical edges (left and right)
  #Fill column vector [-1, 1] 1*L times
  T1x = repeat([-1.0, 1.0], 1, L)   # 2 × L
  T1y = repeat(t', 2, 1)        # 2 × L

  # horizontal edges (top and bottom)
  T2x = repeat(t, 1, 2)         # L × 2
  T2y = repeat([-1.0, 1.0]', L, 1)  # L × 2

  # define 5 colors
  colors = (
    RGBf(1, 0, 0),         # red
    RGBf(1, 0.5, 0.25),    # orange
    RGBf(0.58, 0, 0.82),   # purple
    RGBf(0, 1, 0),         # green
    RGBf(0, 0, 1),         # blue
  )

  fig = Figure(size=(650, 650))
  ax = Axis(fig[1, 1];
    aspect=DataAspect(),
    xlabel=latexstring("\$x\$"),
    ylabel=latexstring("\$y\$"),
    xlabelsize=22,      # label font sizes
    ylabelsize=22,
    xticklabelsize=18,  # tick label sizes
    yticklabelsize=18,
  )

  flag = isnothing(flag) ? 0 : 1

  Zx1 = similar(T1x)
  Zy1 = similar(T1y)
  Zx2 = similar(T2x)
  Zy2 = similar(T2y)

  for p in 1:d.Npat

    k = d.pths[p].reg   # region index

    # left/right boundaries
    # Zx1, Zy1 = mapxy(d, T1x, T1y, p)
    mapxy!(Zx1, Zy1, d, T1x, T1y, p)
    lines!(ax, vec(Zx1[1, :]), vec(Zy1[1, :]), color=colors[k], linewidth=2)
    lines!(ax, vec(Zx1[2, :]), vec(Zy1[2, :]), color=colors[k], linewidth=2)

    # top/bottom boundaries
    # Zx2, Zy2 = mapxy(d, T2x, T2y, p)
    mapxy!(Zx2, Zy2, d, T2x, T2y, p)
    lines!(ax, vec(Zx2[:, 1]), vec(Zy2[:, 1]), color=colors[k], linewidth=2)
    lines!(ax, vec(Zx2[:, 2]), vec(Zy2[:, 2]), color=colors[k], linewidth=2)

    if flag == 1
      # write patch number at its center
      cx, cy = mapxy(d, 0.0, 0.0, p)
      lab = latexstring("\$" * string(p) * "\$")
      text!(ax, cx, cy;
        text=lab, align=(:center, :center),
        fontsize=16, color=:black)
    end

  end

  #xlims!(ax, d.A - d.R1 - 0.1, d.A + d.R1 + 0.1)
  #ylims!(ax, d.B - d.R2 - 0.1, d.B + d.R2 + 0.1)

  display(GLMakie.Screen(),fig)
  return (fig, ax)
  #savefig(plt, "myplot.svg")
end

#-----------------------

"""
    gamx(d::ellipse, t::Float64, k::Int)
    gamx(d::ellipse, t::StridedArray{Float64}, k::Int)

First (x) coordinate of the boundary parametrization γ(t).

- `gamx(d, t, k)` computes the **right boundary** of patch `k`, where `t ∈ [-1,1]`.

`t` may be a scalar or an array; the return has the same shape.
"""
function gamx(d::ellipse, t::Float64, k::Int)::Float64
    p = d.pths[k]
    # Map t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
    xi2 = (p.tk1 - p.tk0) * t / 2 + (p.tk0 + p.tk1) / 2
    t = π*(xi2 / 2 + (2*p.reg - 3) / 4)
    return d.A + d.R1 * d.Cθ *cos(t) - d.R2 * d.Sθ *sin(t)
end

function gamx!(out::StridedArray{Float64}, d::ellipse, t::StridedArray{Float64}, k::Int)
    p  = d.pths[k]
    αt = (p.tk1 - p.tk0) / 2
    βt = p.tk0 + αt
    C  = π * (2 * p.reg - 3) / 4
    c₁₁ = d.R1 * d.Cθ
    c₂₂ =-d.R2 * d.Sθ

    @inbounds for I in eachindex(out, t)
        xi2 = muladd(αt, t[I], βt)      # ξ₂
        τ   = π * (xi2 / 2) + C        
        s, c = sincos(τ)
        out[I] = muladd(c₂₂, s, muladd(c₁₁, c, d.A))
    end
    return nothing
end

"""
    gamy(d::ellipse, t::Float64, k::Int)
    gamy(d::ellipse, t::StridedArray{Float64}, k::Int)

Second (y) coordinate of the boundary parametrization γ(t).

- `gamy(d, t, k)` computes the right boundary of patch `k`, where `t ∈ [-1,1]`.

`t` can be a scalar `Float64` or any `StridedArray{Float64}`; the output matches the input shape.
"""
function gamy(d::ellipse, t::Float64, k::Int)::Float64
    p = d.pths[k]
    # Map t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
    xi2   = (p.tk1 - p.tk0) * t / 2 + (p.tk0 + p.tk1) / 2
    t = π *( xi2 / 2 + (2*p.reg - 3) / 4 )
    return d.B + d.R1 * d.Sθ *cos(t) + d.R2 * d.Cθ *sin(t)
end

function gamy!(out::StridedArray{Float64}, d::ellipse, t::StridedArray{Float64}, k::Int)
    p  = d.pths[k]
    αt = (p.tk1 - p.tk0) / 2
    βt = p.tk0 + αt
    C  = π * (2 * p.reg - 3) / 4
    c₁₂ = d.R1 * d.Sθ
    c₂₁ = d.R2 * d.Cθ

    @inbounds for I in eachindex(out, t)
        xi2 = muladd(αt, t[I], βt)      # ξ₂
        τ   = π * (xi2 / 2) + C        
        s, c = sincos(τ)
        out[I] = muladd(c₂₁, s, muladd(c₁₂, c, d.B))  
    end
    return nothing
end

function gam(d::ellipse, t::Float64, k::Int)::Tuple{Float64,Float64}
 
  p  = d.pths[k]
  αt = (p.tk1 - p.tk0) / 2
  βt = p.tk0 + αt
  C  = π * (2 * p.reg - 3) / 4
  τ = π * muladd(αt, t, βt) / 2 + C

  c₁₁ = d.R1 * d.Cθ
  c₁₂ = d.R1 * d.Sθ
  c₂₁ = d.R2 * d.Cθ
  c₂₂ = -d.R2 * d.Sθ

  s, c = sincos(τ)

  #d.A + d.R1 * d.Cθ * c - d.R2 * d.Sθ * s
  Yx = muladd(c₂₂, s, muladd(c₁₁, c, d.A))
  #d.B + d.R1 * d.Sθ * c + d.R2 * d.Cθ * s
  Yy = muladd(c₂₁, s, muladd(c₁₂, c, d.B))

  return Yx, Yy
end

function gam!(out::Vector{Float64}, d::ellipse, t::Float64, k::Int)
 
  p  = d.pths[k]
  αt = (p.tk1 - p.tk0) / 2
  βt = p.tk0 + αt
  C  = π * (2 * p.reg - 3) / 4
  τ = π * muladd(αt, t, βt) / 2 + C

  c₁₁ = d.R1 * d.Cθ
  c₁₂ = d.R1 * d.Sθ
  c₂₁ = d.R2 * d.Cθ
  c₂₂ = -d.R2 * d.Sθ

  s, c = sincos(τ)

  out[1] = muladd(c₂₂, s, muladd(c₁₁, c, d.A))
  out[2] = muladd(c₂₁, s, muladd(c₁₂, c, d.B))

  return nothing
end

function gam!(out::Matrix{Float64}, d::ellipse, t::Vector{Float64}, k::Int)
  #Out is a Matrix of same columns as length of t and 2 rows
  p  = d.pths[k]
  αt = (p.tk1 - p.tk0) / 2
  βt = p.tk0 + αt
  C  = π * (2 * p.reg - 3) / 4

  c₁₁ = d.R1 * d.Cθ
  c₁₂ = d.R1 * d.Sθ
  c₂₁ = d.R2 * d.Cθ
  c₂₂ = -d.R2 * d.Sθ

  @inbounds  for I in eachindex(t)
    τ = π * muladd(αt, t[I], βt) / 2 + C

    s, c = sincos(τ)

    out[1, I] = muladd(c₂₂, s, muladd(c₁₁, c, d.A)) #x coord
    out[2, I] = muladd(c₂₁, s, muladd(c₁₂, c, d.B)) #y coord
  end

  return nothing
end

function drawbd(d::ellipse)
  # Number of sampling points
  L = 33
  t = collect(range(-1, 1, L))

  # define 5 colors
  clr = (
    RGBf(1, 0, 0),         # red
    RGBf(1, 0.5, 0.25),    # orange
    RGBf(0.58, 0, 0.82),   # purple
    RGBf(0, 1, 0),         # green
    RGBf(0, 0, 1),         # blue
  )

  fig = Figure(size=(650, 650))

  ax = Axis(fig[1, 1];
    aspect=DataAspect(),
    xlabel=latexstring("\$x\$"),
    ylabel=latexstring("\$y\$"),
    xlabelsize=22,      # label font sizes
    ylabelsize=22,
    xticklabelsize=18,  # tick label sizes
    yticklabelsize=18,
  )

  ll = 1

  Tx = similar(t)
  Ty = similar(t)

  for p in d.kd
    k = d.pths[p].reg  # region index

    # right boundaries of patch
    #Tx = gamx(d, t, p)
    #Ty = gamy(d, t, p)
    gamx!(Tx, d, t, p)
    gamy!(Ty, d, t, p)
    # Plot boundary curve in colour of region
    lines!(ax, Tx, Ty; color=clr[k], linewidth=2)

    cx, cy = mapxy(d, 0.0, 0.0, p)
    lab = latexstring("\$" * string(p) * "\$")
    text!(ax, cx, cy; text=lab, align=(:center, :center),
      fontsize=16, color=:black)

    ll += 1
  end

  display(GLMakie.Screen(),fig)
  return (fig, ax)
  #savefig(plt, "myplot.svg")
end

"""
    dgamx(d::ellipse, t::Float64, k::Int) -> Float64
    dgamx(d::ellipse, t::StridedArray{Float64}, k::Int) -> StridedArray{Float64}

Derivative of the first coordinate of the boundary parametrization.

- With `k`: derivative of the **right boundary** of patch `k`
  (i.e., `∂/∂t mapx(d, 1, t, k)` for `t ∈ [-1,1]`).
"""
function dgamx(d::ellipse, t::Float64, k::Int)::Float64
    p  = d.pths[k]
    ht = p.tk1 - p.tk0
    t = ht * (t + 1) / 2 + p.tk0
    τ = π * (t / 2 + (2*p.reg - 3) / 4)
    s, c = sincos(τ)
    return -π * ht * ( d.R1 * d.Cθ * s + d.R2 * d.Sθ * c ) /4
end

function dgamx!(out::StridedArray{Float64}, d::ellipse, t::StridedArray{Float64}, k::Int)
    p  = d.pths[k]
    ht = p.tk1 - p.tk0
    
    αt = ht / 2
    βt = p.tk0 + αt

    C  = π * (2 * p.reg - 3) / 4

    c₁ = -(π * ht / 4) * d.R1 * d.Cθ
    c₂ = -(π * ht / 4) * d.R2 * d.Sθ

    @inbounds for I in eachindex(out, t)
        τ = π * muladd(αt, t[I], βt) / 2 + C

        s, c = sincos(τ)

        out[I] = muladd(c₂, c, c₁ * s)  
    end
    return nothing
end
"""
    dgamy(d::ellipse, t::Float64, k::Int) -> Float64
    dgamy(d::ellipse, t::StridedArray{Float64}, k::Int) -> StridedArray{Float64}

Derivative of the second coordinate of the boundary parametrization γ(t).

- With `k`: derivative of the **right boundary** of patch `k`
  (`∂/∂t mapy(d, 1, t, k)` for `t ∈ [-1,1]`).
"""
function dgamy(d::ellipse, t::Float64, k::Int)::Float64
    p  = d.pths[k]
    ht = p.tk1 - p.tk0
    t = ht * (t + 1) / 2 + p.tk0
    τ = π * (t / 2 + (2*p.reg - 3) / 4)
    s, c = sincos(τ)
    return (π * ht / 4) * (-d.R1 * d.Sθ * s + d.R2 * d.Cθ * c)
end

function dgamy!(out::StridedArray{Float64}, d::ellipse, t::StridedArray{Float64}, k::Int)
    p  = d.pths[k]
    ht = p.tk1 - p.tk0
    
    αt = ht / 2
    βt = p.tk0 + αt

    C  = π * (2 * p.reg - 3) / 4

    c₁ =-(π * ht / 4) * d.R1 * d.Sθ
    c₂ = (π * ht / 4) * d.R2 * d.Cθ

    @inbounds for I in eachindex(out, t)
        τ = π * muladd(αt, t[I], βt) / 2 + C

        s, c = sincos(τ)

        out[I] = muladd(c₂, c, c₁ * s)  
    end
    return nothing
end

function dgam(d::ellipse, t::Float64, k::Int)::Tuple{Float64,Float64}

  p = d.pths[k]
  ht = p.tk1 - p.tk0

  αt = ht / 2
  βt = p.tk0 + αt

  C = π * (2 * p.reg - 3) / 4

  τ = π * muladd(αt, t, βt) / 2 + C

  s, c = sincos(τ)

  # dgamx coefficients
  c₁ = -(π * ht / 4) * d.R1 * d.Cθ
  c₂ = -(π * ht / 4) * d.R2 * d.Sθ
  # dgamy coefficients
  c₃ = -(π * ht / 4) * d.R1 * d.Sθ
  c₄ = (π * ht / 4) * d.R2 * d.Cθ

  outx = muladd(c₂, c, c₁ * s)  
  outy = muladd(c₄, c, c₃ * s)  
  
  return outx, outy
end

function dgam!(out::Vector{Float64}, d::ellipse, t::Float64, k::Int)

  p = d.pths[k]
  ht = p.tk1 - p.tk0

  αt = ht / 2
  βt = p.tk0 + αt

  C = π * (2 * p.reg - 3) / 4

  τ = π * muladd(αt, t, βt) / 2 + C

  s, c = sincos(τ)

  # dgamx coefficients
  c₁ = -(π * ht / 4) * d.R1 * d.Cθ
  c₂ = -(π * ht / 4) * d.R2 * d.Sθ
  # dgamy coefficients
  c₃ = -(π * ht / 4) * d.R1 * d.Sθ
  c₄ = (π * ht / 4) * d.R2 * d.Cθ

  out[1] = muladd(c₂, c, c₁ * s)  
  out[2] = muladd(c₄, c, c₃ * s)  

  return nothing
end

function dgam!(out::Matrix{Float64}, d::ellipse, t::Vector{Float64}, k::Int)

  p = d.pths[k]
  ht = p.tk1 - p.tk0

  αt = ht / 2
  βt = p.tk0 + αt

  C = π * (2 * p.reg - 3) / 4

  # dgamx coefficients
  c₁ = -(π * ht / 4) * d.R1 * d.Cθ
  c₂ = -(π * ht / 4) * d.R2 * d.Sθ
  # dgamy coefficients
  c₃ = -(π * ht / 4) * d.R1 * d.Sθ
  c₄ = (π * ht / 4) * d.R2 * d.Cθ
  @inbounds  for I in eachindex(t)
    τ = π * muladd(αt, t[I], βt) / 2 + C

    s, c = sincos(τ)

    out[1, I] = muladd(c₂, c, c₁ * s)  
    out[2, I] = muladd(c₄, c, c₃ * s)  
  end

  return nothing
end

function gamp!(out::Matrix{Float64}, d::ellipse, t::Vector{Float64}, k::Int)
  #Outputs dy, -dx, the sign is changed in dgamx coeffs!
  p = d.pths[k]
  ht = p.tk1 - p.tk0

  αt = ht / 2
  βt = p.tk0 + αt

  C = π * (2 * p.reg - 3) / 4

  # dgamx coefficients
  c₁ = (π * ht / 4) * d.R1 * d.Cθ
  c₂ = (π * ht / 4) * d.R2 * d.Sθ
  # dgamy coefficients
  c₃ = -(π * ht / 4) * d.R1 * d.Sθ
  c₄ = (π * ht / 4) * d.R2 * d.Cθ
  @inbounds  for I in eachindex(t)
    τ = π * muladd(αt, t[I], βt) / 2 + C

    s, c = sincos(τ)

    out[1, I] = muladd(c₄, c, c₃ * s)
    out[2, I] = muladd(c₂, c, c₁ * s)
  end

  return nothing
end

function nu!(out::Vector{Float64}, d::ellipse, t::Float64, k::Int) 
  p = d.pths[k]
  ht = p.tk1 - p.tk0

  αt = ht / 2
  βt = p.tk0 + αt

  C = π * (2 * p.reg - 3) / 4

  τ = π * muladd(αt, t, βt) / 2 + C

  s, c = sincos(τ)

  # minus dgamx coefficients
  c₁ = (π * ht / 4) * d.R1 * d.Cθ
  c₂ = (π * ht / 4) * d.R2 * d.Sθ
  # dgamy coefficients
  c₃ = -(π * ht / 4) * d.R1 * d.Sθ
  c₄ = (π * ht / 4) * d.R2 * d.Cθ

  dx = muladd(c₂, c, c₁ * s)  
  dy = muladd(c₄, c, c₃ * s)  
  S  = sqrt(dx^2+dy^2)

  out[1] = dy/S
  out[2] = dx/S

  return nothing
end

function nu!(out::Matrix{Float64}, d::ellipse, t::Vector{Float64}, k::Int) 
  #Outputs dy, -dx, the sign is changed in dgamx coeffs!
  p = d.pths[k]
  ht = p.tk1 - p.tk0

  αt = ht / 2
  βt = p.tk0 + αt

  C = π * (2 * p.reg - 3) / 4

  # minus dgamx coefficients
  c₁ = (π * ht / 4) * d.R1 * d.Cθ
  c₂ = (π * ht / 4) * d.R2 * d.Sθ
  # dgamy coefficients
  c₃ = -(π * ht / 4) * d.R1 * d.Sθ
  c₄ = (π * ht / 4) * d.R2 * d.Cθ
  @inbounds  for I in eachindex(t)
    τ = π * muladd(αt, t[I], βt) / 2 + C

    s, c = sincos(τ)

    dx = muladd(c₂, c, c₁ * s)  
    dy = muladd(c₄, c, c₃ * s)  
    S  = sqrt(dx^2+dy^2)

    out[1, I] = dy/S
    out[2, I] = dx/S
  end

  return nothing
end

"""
    DLP(d::ellipse, t::Float64, l::Int, tau::StridedArray{Float64}, k::Int) -> StridedArray{Float64}

Double-layer kernel on the boundary:
K(t, τ) = ((γ_k(τ) - γ_l(t)) ⋅ γ'_k(τ)) / ‖γ_k(τ) - γ_l(t)‖² for k ≠ l,
and for k == l the limiting value is taken for patch k.
The array method returns an array with the ***same shape*** as `tau`.
"""
function DLP!(out::StridedArray{Float64}, d::ellipse, t::Float64, l::Int, 
  tau::StridedArray{Float64},k::Int, x::Vector{Float64}, 
  G::Vector{Float64}, GP::Vector{Float64})
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

    aτ, bτ = ht / 2, p.tk0 + ht / 2

    t̂ = muladd(aτ, t, bτ)

    C = π * (2 * p.reg - 3) / 4

    pref  = ht * π * d.R1 * d.R2 / 8

    @inbounds for I in eachindex(tau, out)

      τ = muladd(aτ, tau[I], bτ)

      tt = (π/4) * (t̂ + τ) + C

      s, c = sincos(tt)
      # r1^2 sin^2 + r2^2 cos^2
      denom = muladd(d.R1^2, s*s, d.R2^2 * (c*c)) 

      out[I] = pref / denom
    end

  end

  return nothing
end

#-----------------------
function Dmap!(out::StridedArray{Float64}, d::ellipse, 
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

    vt = π * (v̂ / 2 + (2 * reg - 3) / 4)
    svt, cvt = sincos(vt)

    if reg == 1
      dux = d.L1 * (2 * v̂ - 1) * d.Sθ / 2 - d.L2 * d.Cθ / 2 + d.R1 * d.Cθ * cvt - d.R2 * d.Sθ * svt
      dvx = (û - 1) * d.L1 * d.Sθ - π * û * (d.R1 * d.Cθ * svt + d.R2 * d.Sθ * cvt) / 2
      duy = d.L1 * (1 - 2 * v̂) * d.Cθ / 2 - d.L2 * d.Sθ / 2 + d.R1 * d.Sθ * cvt + d.R2 * d.Cθ * svt
      dvy = (1 - û) * d.L1 * d.Cθ + π * û * (-d.R1 * d.Sθ * svt + d.R2 * d.Cθ * cvt) / 2

    elseif reg == 2
      dux = d.L2 * (2 * v̂ - 1) * d.Cθ / 2 + d.L1 * d.Sθ / 2 + d.R1 * d.Cθ * cvt - d.R2 * d.Sθ * svt
      dvx = (û - 1) * d.L2 * d.Cθ - π * û * (d.R1 * d.Cθ * svt + d.R2 * d.Sθ * cvt) / 2
      duy = d.L2 * (2 * v̂ - 1) * d.Sθ / 2 - d.L1 * d.Cθ / 2 + d.R1 * d.Sθ * cvt + d.R2 * d.Cθ * svt
      dvy = (û - 1) * d.L2 * d.Sθ + π * û * (-d.R1 * d.Sθ * svt + d.R2 * d.Cθ * cvt) / 2

    elseif reg == 3
      dux = -d.L1 * (2 * v̂ - 1) * d.Sθ / 2 + d.L2 * d.Cθ / 2 + d.R1 * d.Cθ * cvt - d.R2 * d.Sθ * svt
      dvx = (1 - û) * d.L1 * d.Sθ - π * û * (d.R1 * d.Cθ * svt + d.R2 * d.Sθ * cvt) / 2
      duy = -d.L1 * (1 - 2 * v̂) * d.Cθ / 2 + d.L2 * d.Sθ / 2 + d.R1 * d.Sθ * cvt + d.R2 * d.Cθ * svt
      dvy = (û - 1) * d.L1 * d.Cθ + π * û * (-d.R1 * d.Sθ * svt + d.R2 * d.Cθ * cvt) / 2

    else # reg == 4
      dux = -d.L2 * (2 * v̂ - 1) * d.Cθ / 2 - d.L1 * d.Sθ / 2 + d.R1 * d.Cθ * cvt - d.R2 * d.Sθ * svt
      dvx = (1 - û) * d.L2 * d.Cθ - π * û * (d.R1 * d.Cθ * svt + d.R2 * d.Sθ * cvt) / 2
      duy = -d.L2 * (2 * v̂ - 1) * d.Sθ / 2 + d.L1 * d.Cθ / 2 + d.R1 * d.Sθ * cvt + d.R2 * d.Cθ * svt
      dvy = (1 - û) * d.L2 * d.Sθ + π * û * (-d.R1 * d.Sθ * svt + d.R2 * d.Cθ * cvt) / 2
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
  DJ::StridedArray{Float64}, d::ellipse, u::StridedArray{Float64}, 
  v::StridedArray{Float64}, k::Int)
  p = d.pths[k]
  hc = p.ck1 - p.ck0
  ht = p.tk1 - p.tk0

  if p.reg != 5
    # xi1 = map u from [-1,1] → [ck0,ck1]; xi2 = map v → [tk0,tk1]
    αu = hc / 2
    αv = ht / 2
    βu = p.ck0 + αu
    βv = p.tk0 + αv
    c₁₁ = d.R1 * d.Cθ
    c₁₂ = d.R1 * d.Sθ
    c₂₁ = d.R2 * d.Cθ
    c₂₂ = -d.R2 * d.Sθ
    # Corner line (Xx,Xy) for this region
    # Iterate once over a single index range I that is 
    # valid for all four arrays.
    @inbounds for I in eachindex(u, v, Zx, Zy, DJ)

      û = muladd(αu, u[I], βu)
      v̂ = muladd(αv, v[I], βv)

      vt = π * (v̂ / 2 + (2 * p.reg - 3) / 4)
      svt, cvt = sincos(vt)

      # fixed edge point for given v̂
      if p.reg == 1
        Xx = d.A + d.L2 * d.Cθ / 2 - d.L1 * d.Sθ * (2 * v̂ - 1) / 2
        Xy = d.B + d.L2 * d.Sθ / 2 + d.L1 * d.Cθ * (2 * v̂ - 1) / 2

        dux = d.L1 * (2 * v̂ - 1) * d.Sθ / 2 - d.L2 * d.Cθ / 2 + c₁₁ * cvt + c₂₂ * svt
        dvx = (û - 1) * d.L1 * d.Sθ + π * û * (-c₁₁ * svt + c₂₂ * cvt) / 2
        duy = d.L1 * (1 - 2 * v̂) * d.Cθ / 2 - d.L2 * d.Sθ / 2 + c₁₂ * cvt + c₂₁ * svt
        dvy = (1 - û) * d.L1 * d.Cθ + π * û * (-c₁₂ * svt + c₂₁ * cvt) / 2

      elseif p.reg == 2
        Xx = d.A + d.L2 * d.Cθ * (1 - 2 * v̂) / 2 - d.L1 * d.Sθ / 2
        Xy = d.B + d.L2 * d.Sθ * (1 - 2 * v̂) / 2 + d.L1 * d.Cθ / 2

        dux = d.L2 * (2 * v̂ - 1) * d.Cθ / 2 + d.L1 * d.Sθ / 2 + c₁₁ * cvt + c₂₂ * svt
        dvx = (û - 1) * d.L2 * d.Cθ + π * û * (-c₁₁ * svt + c₂₂ * cvt) / 2
        duy = d.L2 * (2 * v̂ - 1) * d.Sθ / 2 - d.L1 * d.Cθ / 2 + c₁₂ * cvt + c₂₁ * svt
        dvy = (û - 1) * d.L2 * d.Sθ + π * û * (-c₁₂ * svt + c₂₁ * cvt) / 2

      elseif p.reg == 3
        Xx = d.A - d.L2 * d.Cθ / 2 + d.L1 * d.Sθ * (2 * v̂ - 1) / 2
        Xy = d.B - d.L2 * d.Sθ / 2 - d.L1 * d.Cθ * (2 * v̂ - 1) / 2

        dux = -d.L1 * (2 * v̂ - 1) * d.Sθ / 2 + d.L2 * d.Cθ / 2 + c₁₁ * cvt + c₂₂ * svt
        dvx = (1 - û) * d.L1 * d.Sθ + π * û * (-c₁₁ * svt + c₂₂ * cvt) / 2
        duy = -d.L1 * (1 - 2 * v̂) * d.Cθ / 2 + d.L2 * d.Sθ / 2 + c₁₂ * cvt + c₂₁ * svt
        dvy = (û - 1) * d.L1 * d.Cθ + π * û * (-c₁₂ * svt + c₂₁ * cvt) / 2

      else # p.reg == 4
        Xx = d.A - d.L2 * d.Cθ * (1 - 2 * v̂) / 2 + d.L1 * d.Sθ / 2
        Xy = d.B - d.L2 * d.Sθ * (1 - 2 * v̂) / 2 - d.L1 * d.Cθ / 2

        dux = -d.L2 * (2 * v̂ - 1) * d.Cθ / 2 - d.L1 * d.Sθ / 2 + c₁₁ * cvt + c₂₂ * svt
        dvx = (1 - û) * d.L2 * d.Cθ + π * û * (-c₁₁ * svt + c₂₂ * cvt) / 2
        duy = -d.L2 * (2 * v̂ - 1) * d.Sθ / 2 + d.L1 * d.Cθ / 2 + c₁₂ * cvt + c₂₁ * svt
        dvy = (1 - û) * d.L2 * d.Sθ + π * û * (-c₁₂ * svt + c₂₁ * cvt) / 2

      end


      #d.A + d.R1 * d.Cθ * c - d.R2 * d.Sθ * s
      Yx = muladd(c₂₂, svt, muladd(c₁₁, cvt, d.A)) 
      #d.B + d.R1 * d.Sθ * c + d.R2 * d.Cθ * s
      Yy = muladd(c₂₁, svt, muladd(c₁₂, cvt, d.B))  

      # blend: (1 - û) * X + û * Y
      Zx[I] = muladd(û, (Yx - Xx), Xx)
      Zy[I] = muladd(û, (Yy - Xy), Xy)

      DJ[I] = αu * αv * abs(dux * dvy - dvx * duy)
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
        Zx[I] = d.A + d.L2*(xi1 - 0.5)*d.Cθ - d.L1*(xi2 - 0.5)*d.Sθ
        Zy[I] = d.B + d.L2*(xi1 - 0.5)*d.Sθ + d.L1*(xi2 - 0.5)*d.Cθ
      end

      fill!(DJ, d.L1 * d.L2 * αu * αv)
    else
      αu = hc / 2
      αv = ht / 2
      βu = p.ck0 + αu
      βv = p.tk0 + αv
      @inbounds  for I in eachindex(u, v, Zx, Zy)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)
        Zx[I] = d.A + d.L2*(xi1 - 0.5)*d.Cθ - d.L1*(xi2 - 0.5)*d.Sθ
        Zy[I] = d.B + d.L2*(xi1 - 0.5)*d.Sθ + d.L1*(xi2 - 0.5)*d.Cθ
      end

      fill!(DJ, d.L1 * d.L2 * αu * αv)
    end
  end

  return nothing
end

function chk_map(d::ellipse)
  # --- Test function f(x,y) ---
  #f!(F, x, y) = fill!(F, 1.0)
  #Iex =  π * d.R1 * d.R2 

  n = 32

  f!(F, x, y) = @. F = x^2 + y^2
  Iex = π * d.R1 * d.R2 * ( (d.R1^2 + d.R2^2) / 4 + d.A^2 + d.B^2 )

  # Chebyshev nodes/weights (once)
  z = cospi.((2 .* (1:n) .- 1) ./ (2n))
  fw = Subroutines.getF1W(n)

  # meshgrid (once)
  zx = repeat(z', n, 1)   # n×n
  zy = repeat(z, 1, n)    # n×n

  # --- preallocate buffers reused across patches ---
  Zx = similar(zx)        # mapped x
  Zy = similar(zy)        # mapped y
  DJ = similar(zx)        # |det J|
  F  = similar(zx)        # integrand buffer f(Zx,Zy)

  I = 0.0
  @inbounds for k in 1:d.Npat
    # map & jacobian in-place
    mapxy_Dmap!(Zx, Zy, DJ, d, zx, zy, k)

    # no temporaries
    f!(F, Zx, Zy)

    # elementwise multiply into F to avoid extra array
    @. F = F * DJ

    I += dot(fw, F, fw)      # accumulate scalar
  end

  # --- report (unchanged) ---
  if abs(Iex - I) < 5e-14
    println("Okay!")
  else
    println("Bug!")
  end
  println("Integral approx: \n", I)
  println("Integral exact : \n", Iex)
  println("Difference     : \n", abs(Iex - I))
  return nothing
end

"""
    jinvmap!(d::ellipse, u::Float64, v::Float64, r::Int) -> tuple

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
function jinvmap(d::ellipse, u::Float64, v::Float64, r::Int)

  û = (u + 1) / 2
  v̂ = (v + 1) / 2

  vt = π * (v̂ / 2 + (2 * r - 3) / 4)
  svt, cvt = sincos(vt)

  if r == 1
    dux = d.L1 * (2 * v̂ - 1) * d.Sθ / 2 - d.L2 * d.Cθ / 2 + d.R1 * d.Cθ * cvt - d.R2 * d.Sθ * svt
    dvx = (û - 1) * d.L1 * d.Sθ - π * û * (d.R1 * d.Cθ * svt + d.R2 * d.Sθ * cvt) / 2
    duy = d.L1 * (1 - 2 * v̂) * d.Cθ / 2 - d.L2 * d.Sθ / 2 + d.R1 * d.Sθ * cvt + d.R2 * d.Cθ * svt
    dvy = (1 - û) * d.L1 * d.Cθ + π * û * (-d.R1 * d.Sθ * svt + d.R2 * d.Cθ * cvt) / 2

  elseif r == 2
    dux = d.L2 * (2 * v̂ - 1) * d.Cθ / 2 + d.L1 * d.Sθ / 2 + d.R1 * d.Cθ * cvt - d.R2 * d.Sθ * svt
    dvx = (û - 1) * d.L2 * d.Cθ - π * û * (d.R1 * d.Cθ * svt + d.R2 * d.Sθ * cvt) / 2
    duy = d.L2 * (2 * v̂ - 1) * d.Sθ / 2 - d.L1 * d.Cθ / 2 + d.R1 * d.Sθ * cvt + d.R2 * d.Cθ * svt
    dvy = (û - 1) * d.L2 * d.Sθ + π * û * (-d.R1 * d.Sθ * svt + d.R2 * d.Cθ * cvt) / 2

  elseif r == 3
    dux = -d.L1 * (2 * v̂ - 1) * d.Sθ / 2 + d.L2 * d.Cθ / 2 + d.R1 * d.Cθ * cvt - d.R2 * d.Sθ * svt
    dvx = (1 - û) * d.L1 * d.Sθ - π * û * (d.R1 * d.Cθ * svt + d.R2 * d.Sθ * cvt) / 2
    duy = -d.L1 * (1 - 2 * v̂) * d.Cθ / 2 + d.L2 * d.Sθ / 2 + d.R1 * d.Sθ * cvt + d.R2 * d.Cθ * svt
    dvy = (û - 1) * d.L1 * d.Cθ + π * û * (-d.R1 * d.Sθ * svt + d.R2 * d.Cθ * cvt) / 2

  else # r == 4
    dux = -d.L2 * (2 * v̂ - 1) * d.Cθ / 2 - d.L1 * d.Sθ / 2 + d.R1 * d.Cθ * cvt - d.R2 * d.Sθ * svt
    dvx = (1 - û) * d.L2 * d.Cθ - π * û * (d.R1 * d.Cθ * svt + d.R2 * d.Sθ * cvt) / 2
    duy = -d.L2 * (2 * v̂ - 1) * d.Sθ / 2 + d.L1 * d.Cθ / 2 + d.R1 * d.Sθ * cvt + d.R2 * d.Cθ * svt
    dvy = (1 - û) * d.L2 * d.Sθ + π * û * (-d.R1 * d.Sθ * svt + d.R2 * d.Cθ * cvt) / 2
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
  mapinv(d::ellipse, u::Float64, v::Float64, k::Int) -> Tuple{Float64,Float64}

Inverse mapping: given a physical point `(u,v)` on patch `k`, return
the reference coordinates `[t; s]` in `[-1,1]^2` such that `mapxy(d, t, s, k) = (u,v)`.

Strategy:
1. 2D Newton using `jinvmap`. If it fails to converge, fall back to a 1D root
   in `t ∈ [-30,30]` (`ApproxFun`) and reconstruct `s`.
"""
@inline function Xx(s::Float64, d::ellipse, r::Int)

  if r == 1
    return d.A + d.L2 * d.Cθ / 2 - d.L1 * d.Sθ * (2 * s - 1) / 2

  elseif r == 2
    return d.A + d.L2 * d.Cθ * (1 - 2 * s) / 2 - d.L1 * d.Sθ / 2

  elseif r == 3
    return d.A - d.L2 * d.Cθ / 2 + d.L1 * d.Sθ * (2 * s - 1) / 2

  else # r == 4
    return d.A - d.L2 * d.Cθ * (1 - 2 * s) / 2 + d.L1 * d.Sθ / 2

  end

end

@inline function Xy(s::Float64, d::ellipse, r::Int)
  if r == 1
    return d.B + d.L2 * d.Sθ / 2 + d.L1 * d.Cθ * (2 * s - 1) / 2

  elseif r == 2
    return d.B + d.L2 * d.Sθ * (1 - 2 * s) / 2 + d.L1 * d.Cθ / 2

  elseif r == 3
    return d.B - d.L2 * d.Sθ / 2 - d.L1 * d.Cθ * (2 * s - 1) / 2

  else # r == 4
    return d.B - d.L2 * d.Sθ * (1 - 2 * s) / 2 - d.L1 * d.Cθ / 2

  end
end

@inline function Yx(s::Float64, d::ellipse, r::Int)
  th = π * (s / 2 + (2 * r - 3) / 4)

  s, c = sincos(th)

  return d.A + d.R1 * d.Cθ * c - d.R2 * d.Sθ * s
end

@inline function Yy(s::Float64, d::ellipse, r::Int)
  th = π * (s / 2 + (2 * r - 3) / 4)

  s, c = sincos(th)

  return d.B + d.R1 * d.Sθ * c + d.R2 * d.Cθ * s
end

function fill_FTable!(tbl::FTable, d::ellipse, r::Int)
  vmin, vmax, N = tbl.vmin, tbl.vmax, tbl.N
  h = (vmax - vmin) / (N - 1)

  P1, P2, P3 = tbl.P1, tbl.P2, tbl.P3

  if r == 1
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A + d.L2 * d.Cθ / 2 - d.L1 * d.Sθ * (2 * v̂ - 1) / 2
      Xy = d.B + d.L2 * d.Sθ / 2 + d.L1 * d.Cθ * (2 * v̂ - 1) / 2
      s, c = sincos(π * (v̂ / 2 - 1 / 4))
      Yx = d.A + d.R1 * d.Cθ * c - d.R2 * d.Sθ * s
      Yy = d.B + d.R1 * d.Sθ * c + d.R2 * d.Cθ * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx
      P3[i] = Xy * Yx - Xx * Yy
    end

  elseif r == 2
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A + d.L2 * d.Cθ * (1 - 2 * v̂) / 2 - d.L1 * d.Sθ / 2
      Xy = d.B + d.L2 * d.Sθ * (1 - 2 * v̂) / 2 + d.L1 * d.Cθ / 2
      s, c = sincos(π * (v̂ / 2 + 1 / 4))
      Yx = d.A + d.R1 * d.Cθ * c - d.R2 * d.Sθ * s
      Yy = d.B + d.R1 * d.Sθ * c + d.R2 * d.Cθ * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx
      P3[i] = Xy * Yx - Xx * Yy
    end

  elseif r == 3
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A - d.L2 * d.Cθ / 2 + d.L1 * d.Sθ * (2 * v̂ - 1) / 2
      Xy = d.B - d.L2 * d.Sθ / 2 - d.L1 * d.Cθ * (2 * v̂ - 1) / 2
      s, c = sincos(π * (v̂ / 2 + 3 / 4))
      Yx = d.A + d.R1 * d.Cθ * c - d.R2 * d.Sθ * s
      Yy = d.B + d.R1 * d.Sθ * c + d.R2 * d.Cθ * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx
      P3[i] = Xy * Yx - Xx * Yy
    end

  else  # r == 4
    @inbounds for i in 1:N
      v̂ = vmin + (i - 1) * h
      Xx = d.A - d.L2 * d.Cθ * (1 - 2 * v̂) / 2 + d.L1 * d.Sθ / 2
      Xy = d.B - d.L2 * d.Sθ * (1 - 2 * v̂) / 2 - d.L1 * d.Cθ / 2
      s, c = sincos(π * (v̂ / 2 + 5 / 4))
      Yx = d.A + d.R1 * d.Cθ * c - d.R2 * d.Sθ * s
      Yy = d.B + d.R1 * d.Sθ * c + d.R2 * d.Cθ * s
      P1[i] = Yy - Xy
      P2[i] = Yx - Xx
      P3[i] = Xy * Yx - Xx * Yy
    end
  end

  tbl.reg = r
  return tbl
end

@inline function f1I(t::Float64, s::Float64,
  d::ellipse, u::Float64, v::Float64, r::Int)
  ŝ = (s + 1) / 2
  t̂ = (t + 1) / 2

  Xxv, Yxv = Xx(ŝ, d, r),  Yx(ŝ, d, r)

  return (1 - t̂) * Xxv + t̂ * Yxv - u
end

@inline function f2I(t::Float64, s::Float64,
  d::ellipse, u::Float64, v::Float64, r::Int)
    ŝ = (s + 1) / 2
    t̂ = (t + 1) / 2

    Xyv, Yyv = Xy(ŝ, d, r), Yy(ŝ, d, r)

    return (1 - t̂) * Xyv + t̂ * Yyv - v
end

@inline function JinvI(t::Float64, s::Float64,
  d::ellipse, u::Float64, v::Float64, r::Int)
  return jinvmap(d, t, s, r)  # returns J11,J12,J21,J22
end

#Used in find_roots! function
@inline function f_cont(v̂::Float64,
  d::ellipse, u::Float64, v::Float64, r::Int)

  if r == 1
    Xxv = d.A + d.L2 * d.Cθ / 2 - d.L1 * d.Sθ * (2 * v̂ - 1) / 2
    Xyv = d.B + d.L2 * d.Sθ / 2 + d.L1 * d.Cθ * (2 * v̂ - 1) / 2
  elseif r == 2
    Xxv = d.A + d.L2 * d.Cθ * (1 - 2 * v̂) / 2 - d.L1 * d.Sθ / 2
    Xyv = d.B + d.L2 * d.Sθ * (1 - 2 * v̂) / 2 + d.L1 * d.Cθ / 2
  elseif r == 3
    Xxv = d.A - d.L2 * d.Cθ / 2 + d.L1 * d.Sθ * (2 * v̂ - 1) / 2
    Xyv = d.B - d.L2 * d.Sθ / 2 - d.L1 * d.Cθ * (2 * v̂ - 1) / 2
  else
    Xxv = d.A - d.L2 * d.Cθ * (1 - 2 * v̂) / 2 + d.L1 * d.Sθ / 2
    Xyv = d.B - d.L2 * d.Sθ * (1 - 2 * v̂) / 2 - d.L1 * d.Cθ / 2
  end

  th = π * (v̂ / 2 + (2 * r - 3) / 4)

  s, c = sincos(th)

  Yxv = d.A + d.R1 * d.Cθ * c - d.R2 * d.Sθ * s
  Yyv = d.B + d.R1 * d.Sθ * c + d.R2 * d.Cθ * s

  return u * (Yyv - Xyv) - v * (Yxv - Xxv) + (Xyv*Yxv - Xxv*Yyv)
end

function mapinv(tbl::FTable, d::ellipse, u::Float64, 
  v::Float64, k::Int)::Tuple{Float64,Float64}

  p = d.pths[k]

  r = p.reg

  # central rectangle (affine inverse; watch axis swap per your map)
  if r == 5
    if d.L2 >= d.L1
      Z1 = xi_inv((((u - d.A)*d.Cθ + (v - d.B)*d.Sθ + d.L2 / 2) / d.L2), p.tk0, p.tk1)
      Z2 = xi_inv((((v - d.B)*d.Cθ - (u - d.A)*d.Sθ + d.L1 / 2) / d.L1), p.ck0, p.ck1)
    else
      Z1 = xi_inv((((u - d.A)*d.Cθ + (v - d.B)*d.Sθ + d.L2 / 2) / d.L2), p.ck0, p.ck1)
      Z2 = xi_inv((((v - d.B)*d.Cθ - (u - d.A)*d.Sθ + d.L1 / 2) / d.L1), p.tk0, p.tk1)
    end
    return Z1, Z2
  end

  # ----- curved patches reg = 1..4 -----
  # --- Stage 1: 2D Newton via Subroutines ---
  tN, sN = newtonR2D(f1I, f2I, JinvI,
    0.0, 0.0, 4, d, u, v, r; tol=5e-14)

  if tN !== :max

    x = xi_inv((1 + tN) / 2, p.ck0, p.ck1)
    y = xi_inv((1 + sN) / 2, p.tk0, p.tk1)

    return x, y

  end

  #--- Stage 2: BIS method inversion ---

  # ensure the table is filled for this region
  if tbl.reg != r
    display("Wrong table!")
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
    ptconv(d::ellipse,  t1::Float64, t2::Float64, idx::Int, ptdest::String)

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
function ptconv(d::ellipse, t1::Float64, t2::Float64, idx::Int, ptdest::String)
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
    mapinv2(d::ellipse, t1::Float64, t2::Float64, k2::Int, k::Int) -> Tuple{Float64,Float64}

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
function mapinv2(d::ellipse, t1::Float64, t2::Float64, k2::Int, k::Int)::Tuple{Float64,Float64}
    p = d.pths[k]  # Fields reg, ck0, ck1, tk0, tk1

    # Convert the given (t1,t2,k2) to the "regional" normalized coords of the *point*
    # Expecting something that returns two values in [-1,1]:
    #   tr1 ≈ u_ref, tr2 ≈ v_ref  (for the *global/regional* parameterization)
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
- `t` may be a scalar `Float64` or any `StridedArray{Float64}`; the return has the same shape.
"""
function dfunc(d::ellipse, k::Int, t::Float64, s::Float64)::Float64

  p = d.pths[k]

  if p.reg == 5
    return 1.0
  end

  hc = p.ck1 - p.ck0
  val = (1.0 - p.ck1) + hc * t / 2
  exp = s ≥ 0.5 ? (s - 1) : s
  return val^exp

end

function dfunc!(out::StridedArray{Float64}, d::ellipse, k::Int, t::StridedArray{Float64}, s::Float64)

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

"""
  diff_map!(out::Matrix{Float64}, Zx::Matrix{Float64}, Zy::Matrix{Float64}, 
  DJ::Matrix{Float64}, d::ellipse, u::Float64, v::Float64,
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
  d::ellipse, u::Float64, v::Float64,
  u2::Matrix{Float64}, v2::Matrix{Float64},
  du::Vector{Float64}, dv::Vector{Float64}, k::Int;
  tol::Float64=1e-5)

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
        Dx = (cDx * d.Cθ) * du[i] - (cDy * d.Sθ) * dj
        Dy = (cDx * d.Sθ) * du[i] + (cDy * d.Cθ) * dj
        out[i, j] = hypot(Dx, Dy)
      end
    end

    return nothing
  end

  # affine maps: [-1,1] → [ck0,ck1] × [tk0,tk1]
  û = muladd(αc, u, p.ck0 + αc)
  v̂ = muladd(αt, v, p.tk0 + αt)

  # angle/trig
  vt = muladd(π / 2, v̂, (2 * p.reg - 3) * (π / 4))
  svt, cvt = sincos(vt)

  # common combos 
  A = d.R1 * d.Cθ * svt + d.R2 * d.Sθ * cvt                  # (R1*ct*svt + R2*st*cvt)
  B = d.R1 * d.Cθ * cvt - d.R2 * d.Sθ * svt                  # (R1*ct*cvt - R2*st*svt)
  C = -d.R1 * d.Sθ * svt + d.R2 * d.Cθ * cvt                 # (-R1*st*svt + R2*ct*cvt)
  D = d.R1 * d.Sθ * cvt + d.R2 * d.Cθ * svt                  # (R1*st*cvt + R2*ct*svt)

  # folded pi/alpha_t powers
  πa = π * αt
  πa2 = πa * πa          # π^2 * αt^2
  πa3 = πa2 * πa         # π^3 * αt^3
  πa4 = πa2 * πa2        # π^4 * αt^4

  if p.reg == 1
    dux = αc * (d.L1 * (2 * v̂ - 1) * d.Sθ / 2 - d.L2 * d.Cθ / 2 + B)
    dvx = αt * ((û - 1) * d.L1 * d.Sθ - (π / 2) * û * A)
    duvx = αt * αc * (d.L1 * d.Sθ - (π / 2) * A)

    dv2x = -(πa2 / 4) * û * B
    duv2x = -(πa2 / 4) * αc * B
    dv3x = (πa3 / 8) * û * A
    duv3x = (πa3 / 8) * αc * A
    dv4x = (πa4 / 16) * û * B

    duy = αc * (d.L1 * (1 - 2 * v̂) * d.Cθ / 2 - d.L2 * d.Sθ / 2 + D)
    dvy = αt * ((1 - û) * d.L1 * d.Cθ + (π / 2) * û * C)
    duvy = αt * αc * (-d.L1 * d.Cθ + (π / 2) * C)

    dv2y = -(πa2 / 4) * û * D
    duv2y = -(πa2 / 4) * αc * D
    dv3y = -(πa3 / 8) * û * C
    duv3y = -(πa3 / 8) * αc * C
    dv4y = (πa4 / 16) * û * D

  elseif p.reg == 2
    dux = αc * (d.L2 * (2 * v̂ - 1) * d.Cθ / 2 + d.L1 * d.Sθ / 2 + B)
    dvx = αt * ((û - 1) * d.L2 * d.Cθ - (π / 2) * û * A)
    duvx = αt * αc * (d.L2 * d.Cθ - (π / 2) * A)

    dv2x = -(πa2 / 4) * û * B
    duv2x = -(πa2 / 4) * αc * B
    dv3x = (πa3 / 8) * û * A
    duv3x = (πa3 / 8) * αc * A
    dv4x = (πa4 / 16) * û * B

    duy = αc * (d.L2 * (2 * v̂ - 1) * d.Sθ / 2 - d.L1 * d.Cθ / 2 + D)
    dvy = αt * ((û - 1) * d.L2 * d.Sθ + (π / 2) * û * C)
    duvy = αt * αc * (d.L2 * d.Sθ + (π / 2) * C)

    dv2y = -(πa2 / 4) * û * D
    duv2y = -(πa2 / 4) * αc * D
    dv3y = -(πa3 / 8) * û * C
    duv3y = -(πa3 / 8) * αc * C
    dv4y = (πa4 / 16) * û * D

  elseif p.reg == 3
    dux = αc * (-d.L1 * (2 * v̂ - 1) * d.Sθ / 2 + d.L2 * d.Cθ / 2 + B)
    dvx = αt * (-(û - 1) * d.L1 * d.Sθ - (π / 2) * û * A)
    duvx = αt * αc * (-d.L1 * d.Sθ - (π / 2) * A)

    dv2x = -(πa2 / 4) * û * B
    duv2x = -(πa2 / 4) * αc * B
    dv3x = (πa3 / 8) * û * A
    duv3x = (πa3 / 8) * αc * A
    dv4x = (πa4 / 16) * û * B

    duy = αc * (-d.L1 * (1 - 2 * v̂) * d.Cθ / 2 + d.L2 * d.Sθ / 2 + D)
    dvy = αt * (-(1 - û) * d.L1 * d.Cθ + (π / 2) * û * C)
    duvy = αt * αc * (d.L1 * d.Cθ + (π / 2) * C)

    dv2y = -(πa2 / 4) * û * D
    duv2y = -(πa2 / 4) * αc * D
    dv3y = -(πa3 / 8) * û * C
    duv3y = -(πa3 / 8) * αc * C
    dv4y = (πa4 / 16) * û * D

  else # p.reg == 4
    dux = αc * (-d.L2 * (2 * v̂ - 1) * d.Cθ / 2 - d.L1 * d.Sθ / 2 + B)
    dvx = αt * (-(û - 1) * d.L2 * d.Cθ - (π / 2) * û * A)
    duvx = αt * αc * (-d.L2 * d.Cθ - (π / 2) * A)

    dv2x = -(πa2 / 4) * û * B
    duv2x = -(πa2 / 4) * αc * B
    dv3x = (πa3 / 8) * û * A
    duv3x = (πa3 / 8) * αc * A
    dv4x = (πa4 / 16) * û * B

    duy = αc * (-d.L2 * (2 * v̂ - 1) * d.Sθ / 2 + d.L1 * d.Cθ / 2 + D)
    dvy = αt * (-(û - 1) * d.L2 * d.Sθ + (π / 2) * û * C)
    duvy = αt * αc * (-d.L2 * d.Sθ + (π / 2) * C)

    dv2y = -(πa2 / 4) * û * D
    duv2y = -(πa2 / 4) * αc * D
    dv3y = -(πa3 / 8) * û * C
    duv3y = -(πa3 / 8) * αc * C
    dv4y = (πa4 / 16) * û * D
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
            d::ellipse, u::Float64, v::Float64,
            u2::Matrix{Float64}, v2::Matrix{Float64}, r::Matrix{Float64},
            du::Vector{Float64}, dv::Vector{Float64}, k::Int;
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
  d::ellipse, u::Float64, v::Float64,
  u2::Matrix{Float64}, v2::Matrix{Float64}, r::Matrix{Float64},
  du::Vector{Float64}, dv::Vector{Float64}, k::Int;
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
        Dx = (cDx * d.Cθ) * dui - (cDy * d.Sθ) * dvi
        Dy = (cDx * d.Sθ) * dui + (cDy * d.Cθ) * dvi
        hD = hypot(Dx, Dy)
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
  vt = muladd(π / 2, v̂, (2 * p.reg - 3) * (π / 4))
  svt, cvt = sincos(vt)

  # common combos 
  A = d.R1 * d.Cθ * svt + d.R2 * d.Sθ * cvt                  # (R1*ct*svt + R2*st*cvt)
  B = d.R1 * d.Cθ * cvt - d.R2 * d.Sθ * svt                  # (R1*ct*cvt - R2*st*svt)
  C = -d.R1 * d.Sθ * svt + d.R2 * d.Cθ * cvt                 # (-R1*st*svt + R2*ct*cvt)
  D = d.R1 * d.Sθ * cvt + d.R2 * d.Cθ * svt                  # (R1*st*cvt + R2*ct*svt)

  # folded pi/alpha_t powers
  πa = π * αt
  πa2 = πa * πa          # π^2 * αt^2
  πa3 = πa2 * πa         # π^3 * αt^3
  πa4 = πa2 * πa2        # π^4 * αt^4

  if p.reg == 1
    dux = αc * (d.L1 * (2 * v̂ - 1) * d.Sθ / 2 - d.L2 * d.Cθ / 2 + B)
    dvx = αt * ((û - 1) * d.L1 * d.Sθ - (π / 2) * û * A)
    duvx = αt * αc * (d.L1 * d.Sθ - (π / 2) * A)

    dv2x = -(πa2 / 4) * û * B
    duv2x = -(πa2 / 4) * αc * B
    dv3x = (πa3 / 8) * û * A
    duv3x = (πa3 / 8) * αc * A
    dv4x = (πa4 / 16) * û * B

    duy = αc * (d.L1 * (1 - 2 * v̂) * d.Cθ / 2 - d.L2 * d.Sθ / 2 + D)
    dvy = αt * ((1 - û) * d.L1 * d.Cθ + (π / 2) * û * C)
    duvy = αt * αc * (-d.L1 * d.Cθ + (π / 2) * C)

    dv2y = -(πa2 / 4) * û * D
    duv2y = -(πa2 / 4) * αc * D
    dv3y = -(πa3 / 8) * û * C
    duv3y = -(πa3 / 8) * αc * C
    dv4y = (πa4 / 16) * û * D

  elseif p.reg == 2
    dux = αc * (d.L2 * (2 * v̂ - 1) * d.Cθ / 2 + d.L1 * d.Sθ / 2 + B)
    dvx = αt * ((û - 1) * d.L2 * d.Cθ - (π / 2) * û * A)
    duvx = αt * αc * (d.L2 * d.Cθ - (π / 2) * A)

    dv2x = -(πa2 / 4) * û * B
    duv2x = -(πa2 / 4) * αc * B
    dv3x = (πa3 / 8) * û * A
    duv3x = (πa3 / 8) * αc * A
    dv4x = (πa4 / 16) * û * B

    duy = αc * (d.L2 * (2 * v̂ - 1) * d.Sθ / 2 - d.L1 * d.Cθ / 2 + D)
    dvy = αt * ((û - 1) * d.L2 * d.Sθ + (π / 2) * û * C)
    duvy = αt * αc * (d.L2 * d.Sθ + (π / 2) * C)

    dv2y = -(πa2 / 4) * û * D
    duv2y = -(πa2 / 4) * αc * D
    dv3y = -(πa3 / 8) * û * C
    duv3y = -(πa3 / 8) * αc * C
    dv4y = (πa4 / 16) * û * D

  elseif p.reg == 3
    dux = αc * (-d.L1 * (2 * v̂ - 1) * d.Sθ / 2 + d.L2 * d.Cθ / 2 + B)
    dvx = αt * (-(û - 1) * d.L1 * d.Sθ - (π / 2) * û * A)
    duvx = αt * αc * (-d.L1 * d.Sθ - (π / 2) * A)

    dv2x = -(πa2 / 4) * û * B
    duv2x = -(πa2 / 4) * αc * B
    dv3x = (πa3 / 8) * û * A
    duv3x = (πa3 / 8) * αc * A
    dv4x = (πa4 / 16) * û * B

    duy = αc * (-d.L1 * (1 - 2 * v̂) * d.Cθ / 2 + d.L2 * d.Sθ / 2 + D)
    dvy = αt * (-(1 - û) * d.L1 * d.Cθ + (π / 2) * û * C)
    duvy = αt * αc * (d.L1 * d.Cθ + (π / 2) * C)

    dv2y = -(πa2 / 4) * û * D
    duv2y = -(πa2 / 4) * αc * D
    dv3y = -(πa3 / 8) * û * C
    duv3y = -(πa3 / 8) * αc * C
    dv4y = (πa4 / 16) * û * D

  else # p.reg == 4
    dux = αc * (-d.L2 * (2 * v̂ - 1) * d.Cθ / 2 - d.L1 * d.Sθ / 2 + B)
    dvx = αt * (-(û - 1) * d.L2 * d.Cθ - (π / 2) * û * A)
    duvx = αt * αc * (-d.L2 * d.Cθ - (π / 2) * A)

    dv2x = -(πa2 / 4) * û * B
    duv2x = -(πa2 / 4) * αc * B
    dv3x = (πa3 / 8) * û * A
    duv3x = (πa3 / 8) * αc * A
    dv4x = (πa4 / 16) * û * B

    duy = αc * (-d.L2 * (2 * v̂ - 1) * d.Sθ / 2 + d.L1 * d.Cθ / 2 + D)
    dvy = αt * (-(û - 1) * d.L2 * d.Sθ + (π / 2) * û * C)
    duvy = αt * αc * (-d.L2 * d.Sθ + (π / 2) * C)

    dv2y = -(πa2 / 4) * û * D
    duv2y = -(πa2 / 4) * αc * D
    dv3y = -(πa3 / 8) * û * C
    duv3y = -(πa3 / 8) * αc * C
    dv4y = (πa4 / 16) * û * D
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

  return nothing
end


"""
    Dwall(d::ellipse, u::Float64, v::Float64, k::Int) -> dvx, dvy

Derivative of the wall curve w.r.t. `v` at the singular side `u = ±1`
for the `k`-th (non-quadrilateral) patch.

Returns a Tuple `dvx, dvy`.
"""
function Dwall(d::ellipse, u::Float64, v::Float64, k::Int)::Tuple{Float64, Float64}

   
  p = d.pths[k]
  hc = p.ck1 - p.ck0
  ht = p.tk1 - p.tk0
  reg = p.reg

  αu = hc / 2
  αv = ht / 2
  βu = p.ck0 + αu      # û = αu * u + βu
  βv = p.tk0 + αv      # v̂ = αv * v + βv

  û = muladd(αu, u, βu)
  v̂ = muladd(αv, v, βv)

  vt = π * (v̂ / 2 + (2 * reg - 3) / 4)
  svt, cvt = sincos(vt)

  if reg == 1
    dvx = (û - 1) * d.L1 * d.Sθ - π * û * (d.R1 * d.Cθ * svt + d.R2 * d.Sθ * cvt) / 2
    dvy = (1 - û) * d.L1 * d.Cθ + π * û * (-d.R1 * d.Sθ * svt + d.R2 * d.Cθ * cvt) / 2

  elseif reg == 2
    dvx = (û - 1) * d.L2 * d.Cθ - π * û * (d.R1 * d.Cθ * svt + d.R2 * d.Sθ * cvt) / 2
    dvy = (û - 1) * d.L2 * d.Sθ + π * û * (-d.R1 * d.Sθ * svt + d.R2 * d.Cθ * cvt) / 2

  elseif reg == 3
    dvx = (1 - û) * d.L1 * d.Sθ - π * û * (d.R1 * d.Cθ * svt + d.R2 * d.Sθ * cvt) / 2
    dvy = (û - 1) * d.L1 * d.Cθ + π * û * (-d.R1 * d.Sθ * svt + d.R2 * d.Cθ * cvt) / 2

  else # reg == 4
    dvx = (1 - û) * d.L2 * d.Cθ - π * û * (d.R1 * d.Cθ * svt + d.R2 * d.Sθ * cvt) / 2
    dvy = (1 - û) * d.L2 * d.Sθ + π * û * (-d.R1 * d.Sθ * svt + d.R2 * d.Cθ * cvt) / 2
  end

  return dvx, dvy
end


"""
    boundquad!(P::SubArray{Float64}, d::ellipse, k::Int) -> nothing

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
function boundquad!(P::SubArray{Float64}, d::ellipse, k::Int)

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
    refine!(d::ellipse, Nc::Int, Nt::Int, K::AbstractVector{<:Int})

Subdivide each patch in `K` into `Nc * Nt` subpatches (split in `c` and `t`),
update `d.pths`, `d.Npat`, recompute `d.kd`, and rebuild `d.Qpts` / `d.Qptsbd`.

Notes
- Patches are re-sorted lexicographically by `(reg, ck0, ck1, tk0, tk1)`
- `Qpts` uses `boundquad(d,k)`; `Qptsbd` uses `boundquadbd(d,k)` for patches in `d.kd`.
"""
function refine!(d::ellipse, Nc::Int, Nt::Int, K::Vector{Int})
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

function Base.show(io::IO, d::ellipse)
  println(io, "ellipse with properties:")
  println(io, "  (A , B )  = (", d.A,",", d.B,")")
  println(io, "  (R₁, R₂)  = (", d.R1,",", d.R2,")")
  println(io, "  (L₁, L₂)  = (", d.L1,",", d.L2,")")
  println(io, "  (cos,sin) = (", d.Cθ,",", d.Sθ,")")
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
