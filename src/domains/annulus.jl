"""
Mutable struct annulus

A  — x-coordinate of annulus center (Float64)

B  — y-coordinate of annulus center (Float64)

R₁₁ — Radius of γ₁ in x axis (Float64)

R₁₂ — Radius of γ₁ in y axis (Float64)

R₂₁ — Radius of γ₂ in x axis (Float64)

R₂₂ — Radius of γ₂ in y axis (Float64)

nₕ — number of holes (nonnegative Int) 

θ₁ — Angle of rotation of the ellipse γ₁ (Float64)

θ₂ — Angle of rotation of the ellipse γ₂ (Float64)

T₁ — Four parametric points in [0,2π) on γ₁. (Vector{Float64})

T₂ — Four parametric points in [0,2π) on γ₂. (Vector{Float64})

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
struct regionparams_annulus
  Cθ₁::Float64
  Sθ₁::Float64
  Cθ₂::Float64
  Sθ₂::Float64
end

mutable struct annulus <: abstractdomain

  A::Float64
  B::Float64
  R11::Float64
  R12::Float64
  R21::Float64
  R22::Float64
  tht1::Float64
  tht2::Float64
  T1::Vector{Float64}
  T2::Vector{Float64}
  nh::Int
  kd::Vector{Int}
  Npat::Int
  pths::Vector{Patch}   # 5 * Npat
  Qpts::Matrix{Float64}   # 8 * Npat
  Qptsbd::Matrix{Float64} # 8 * Npat
  RP::regionparams_annulus

  function annulus(; b,
    a=nothing, A=nothing, B=nothing,
    R11=nothing, R12=nothing, 
    R21=nothing, R22=nothing, 
    tht1=nothing, tht2=nothing, 
    ck=nothing, tk=nothing)
    #Construct an intsance for this structure

    """
    Required:
      - b :: Vector{Int} of length 8 (Col), positive entries
    
    Optional (all keywords):
      - a  :: Vector{Int} length 8 (defaults to ceil.(b/3))
      - A,B:: Center of the annulus (defaults (0,0))
      - R11 :: Float64 (default 1.0)
      - R12 :: Float64 (default 1.0)
      - R21 :: Float64 (default 0.5)
      - R22 :: Float64 (default 0.5)
      - θ₁ :: Float64 (default 0.0)
      - θ₂ :: Float64 (default 0.0)
      - ck :: vector{vector{Float64}} (default equispaced per a)
      - tk :: vector{vector{Float64}} (default equispaced per b)
    """
    @assert isa(b, AbstractVector{Int}) && length(b) == 8

    #The annulus has obviously one hole
    nh = 1

    a = isnothing(a) ? ceil.(Int, b ./ 3) : collect(Int, a)

    #By default, annulus at origin
    A = something(A, 0.0)
    B = something(B, 0.0)

    R11 = something(R11, 1.0)
    R12 = something(R12, 1.0)
    R21 = something(R21, 0.5)
    R22 = something(R22, 0.5)

    tht1 = something(tht1, 0.0)
    tht2 = something(tht2, 0.0)

    #Sum up all the subpatches
    Npat = dot(a, b)

    s₁, c₁ = sincos(tht1)
    s₂, c₂ = sincos(tht2)

    RP= regionparams_annulus(c₁, s₁, c₂, s₂)

    #--------- Computing T1 and T2 ---------

    if tht1 == tht2 #both the curves are rotated the same amount

      T1 = π .* [0, 1 / 2, 1, 3 / 2]

      T2 = π .* [0, 1 / 2, 1, 3 / 2]

    elseif R11 == R12 #γ₁ is a circle

      T1 = tht2 .+ π .* [0, 1 / 2, 1, 3 / 2]

      T2 = π .* [0, 1 / 2, 1, 3 / 2]

    elseif R21 == R22 #γ₂ is a circle

      T2 = tht1 .+ π .* [0, 1 / 2, 1, 3 / 2]

      T1 = π .* [0, 1 / 2, 1, 3 / 2]

    else
      T1 = Vector{Float64}(undef, 4)
      T2 = Vector{Float64}(undef, 4)

      # Unit directions for ellipse 
      e1x, e1y = c₁, s₁
      e2x, e2y = -s₁, c₁

      # ---- γ1 and derivatives (scalar, alloc-free) ----
      @inline g1x(t::Float64) = A + R11 * cos(t) * c₁ - R12 * sin(t) * s₁
      @inline g1y(t::Float64) = B + R11 * cos(t) * s₁ + R12 * sin(t) * c₁

      @inline dg1x(t::Float64) = -R11 * sin(t) * c₁ - R12 * cos(t) * s₁
      @inline dg1y(t::Float64) = -R11 * sin(t) * s₁ + R12 * cos(t) * c₁

      # ddg1 = -(g1 - center)  
      @inline ddg1x(t::Float64) = - R11 * cos(t) * c₁ + R12 * sin(t) * s₁
      @inline ddg1y(t::Float64) = - R11 * cos(t) * s₁ - R12 * sin(t) * c₁

      # ---- γ2 ----
      @inline g2x(t::Float64) = A + R21 * cos(t) * c₂ - R22 * sin(t) * s₂
      @inline g2y(t::Float64) = B + R21 * cos(t) * s₂ + R22 * sin(t) * c₂

      # Distance from point (x,y) to γ1 via "sample + Newton" (inlined dist_gam)
      # Returns (Dmin, idx) but we only need Dmin for GSS. 
      @inline function dist_to_g1(x::Float64, y::Float64)::Float64
        L = 9

        dt = 2.0 * π / (L - 1)

        Dmin = Inf

        @inbounds for i in 1:(L-1)
          t0 = (i - 1) * dt

          # Define h(t) and dh(t) for this fixed (x,y) WITHOUT allocating vectors
          # h(t)  = (g1(t)-X)·g1'(t)
          # dh(t) = (g1(t)-X)·g1''(t) + ||g1'(t)||^2
          @inline h(t::Float64) = (g1x(t) - x) * dg1x(t) + (g1y(t) - y) * dg1y(t)
          @inline dh(t::Float64) = (g1x(t) - x) * ddg1x(t) + (g1y(t) - y) * ddg1y(t) +
                           dg1x(t) * dg1x(t) + dg1y(t) * dg1y(t)

          t = newtonR1D(h, dh, t0, 32; tol=1e-5)

          if isfinite(t) && (t >= 0.0) && (t <= 2.0 * π)
            dx = g1x(t) - x
            dy = g1y(t) - y
            d = sqrt(dx * dx + dy * dy)
            if d < Dmin
              Dmin = d
            end
          end
        end

        return Dmin
      end
      
      fpi = Float64(pi)

      # dm(x) = distance from point on γ2 at parameter x to γ1
      @inline dm(x::Float64) = dist_to_g1(g2x(x), g2y(x))

      # Find two "farthest" points on γ2 relative to γ1, one in each half interval
      #There are some allocations present inside GSS, but they are minimal. 
      _, t1 = GSS(x -> dm(x), 0.0, fpi - 1e-14, 1e-5)
      _, t2 = GSS(x -> dm(x), fpi, 2.0 * fpi - 1e-14, 1e-5)

      # Build T2 around t1/t2 
      α = π / 4.0
      T2[1] = t1 - α
      T2[2] = t1 + α
      T2[3] = t2 - α
      T2[4] = t2 + α

      # For each T2[i], shoot a ray from γ₂(T2[i]) in normal direction to γ₂.
      # direction v = [dg2y(to), -dg2x(to)] intersect with ellipse 1 (Outer)
      @inbounds for i in 1:4
        s, c = sincos(T2[i])

        ux = A + R21 * c * c₂ - R22 * s * s₂ #γ₂x(T₂[i])
        uy = B + R21 * c * s₂ + R22 * s * c₂ #γ₂y(T₂[i])

        vx = -R21 * s * s₂ + R22 * c * c₂  # dγ₂y(T₂[i])
        vy =  R21 * s * c₂ + R22 * c * s₂  #-dγ₂x(T₂[i])

        Et, _ = ellipserinter(A, B, R11, R12, e1x, e1y, e2x, e2y, ux, uy, vx, vy)

        T1[i] = Et
      end

      # unwrap to enforce increasing order 
      @inbounds for i in 1:3
        if T1[i+1] < T1[i]
          T1[i+1] += 2.0 * pi
        end
        if T2[i+1] < T2[i]
          T2[i+1] += 2.0 * pi
        end
      end

    end

    # ck, tk as vector-of-vectors with *no padding*:
    #   ck[k] has length a[k] + 1
    #   tk[k] has length b[k] + 1
    ck = something(ck, [(0:a[k]) ./ a[k] for k in 1:8])
    tk = something(tk, [(0:b[k]) ./ b[k] for k in 1:8])

    # Preallocate pths and kd
    pths = Vector{Patch}(undef, Npat)

    #Starting with the assumption that all interioir 
    #patches touch the boundary
    kd = Vector{Int}(undef, Npat)   # upper bound length
    nkd = 0                         # actual count

    # Fill pths and kd in a single pass
    idx = 1
    @inbounds for k in 1:8
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

          if p.ck1 == 1.0
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

    d = new(A, B, R11, R12, R21, R22, tht1, tht2, T1, T2, nh, kd, Npat, pths, Qpts, Qptsbd, RP)
    
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

function mapx(d::annulus, u::Float64, v::Float64, k::Int)::Float64

  # p denotes information of kth patch
  p = d.pths[k]

  # Difference in c values
  hc = p.ck1 - p.ck0

  xi1 = hc * u / 2 + (p.ck1 + p.ck0) / 2

  # Difference in t values
  ht = p.tk1 - p.tk0

  xi2 = ht * v / 2 + (p.tk1 + p.tk0) / 2

  c₁ = d.RP.Cθ₁
  s₁ = d.RP.Sθ₁
  c₂ = d.RP.Cθ₂
  s₂ = d.RP.Sθ₂

  if p.reg == 1 

    th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
    th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

    th1 = muladd(Δth1, xi2, th₀1) 
    th5 = muladd(Δth5, xi2, th₀5) 

    st1, ct1 = sincos(th1)
    st5, ct5 = sincos(th5)

    y1 = d.A + d.R11 * ct1 * c₁ - d.R12 * st1 * s₁

    y5 = d.A + d.R21 * ct5 * c₂ - d.R22 * st5 * s₂

    return muladd((y1 - y5) / 2, xi1, (y1 + y5) / 2)

  elseif p.reg == 5

    th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
    th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

    th1 = muladd(Δth1, xi2, th₀1)
    th5 = muladd(Δth5, xi2, th₀5)

    st1, ct1 = sincos(th1)
    st5, ct5 = sincos(th5)

    y1 = d.A + d.R11 * ct1 * c₁ - d.R12 * st1 * s₁

    y5 = d.A + d.R21 * ct5 * c₂ - d.R22 * st5 * s₂

    return muladd((y5 - y1) / 2, xi1, (y5 + y1) / 2)

  elseif p.reg == 2

    th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
    th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

    th2 = muladd(Δth2, xi2, th₀2)
    th6 = muladd(Δth6, xi2, th₀6)

    st2, ct2 = sincos(th2)
    st6, ct6 = sincos(th6)

    y2 = d.A + d.R11 * ct2 * c₁ - d.R12 * st2 * s₁

    y6 = d.A + d.R21 * ct6 * c₂ - d.R22 * st6 * s₂

    return muladd((y2 - y6) / 2, xi1, (y2 + y6) / 2)

  elseif p.reg == 6

    th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
    th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

    th2 = muladd(Δth2, xi2, th₀2)
    th6 = muladd(Δth6, xi2, th₀6)

    st2, ct2 = sincos(th2)
    st6, ct6 = sincos(th6)

    y2 = d.A + d.R11 * ct2 * c₁ - d.R12 * st2 * s₁

    y6 = d.A + d.R21 * ct6 * c₂ - d.R22 * st6 * s₂

    return muladd((y6 - y2) / 2, xi1, (y6 + y2) / 2)

  elseif p.reg == 3

    th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
    th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

    th3 = muladd(Δth3, xi2, th₀3)
    th7 = muladd(Δth7, xi2, th₀7)

    st3, ct3 = sincos(th3)
    st7, ct7 = sincos(th7)

    y3 = d.A + d.R11 * ct3 * c₁ - d.R12 * st3 * s₁

    y7 = d.A + d.R21 * ct7 * c₂ - d.R22 * st7 * s₂

    return muladd((y3 - y7) / 2, xi1, (y3 + y7) / 2)

  elseif p.reg == 7

    th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
    th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

    th3 = muladd(Δth3, xi2, th₀3)
    th7 = muladd(Δth7, xi2, th₀7)

    st3, ct3 = sincos(th3)
    st7, ct7 = sincos(th7)

    y3 = d.A + d.R11 * ct3 * c₁ - d.R12 * st3 * s₁

    y7 = d.A + d.R21 * ct7 * c₂ - d.R22 * st7 * s₂

    return muladd((y7 - y3) / 2, xi1, (y7 + y3) / 2)

  elseif p.reg == 4

    th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
    th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

    th4 = muladd(Δth4, xi2, th₀4)
    th8 = muladd(Δth8, xi2, th₀8)

    st4, ct4 = sincos(th4)
    st8, ct8 = sincos(th8)

    y4 = d.A + d.R11 * ct4 * c₁ - d.R12 * st4 * s₁

    y8 = d.A + d.R21 * ct8 * c₂ - d.R22 * st8 * s₂

    return muladd((y4 - y8) / 2, xi1, (y4 + y8) / 2)

  else

    th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
    th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

    th4 = muladd(Δth4, xi2, th₀4)
    th8 = muladd(Δth8, xi2, th₀8)

    st4, ct4 = sincos(th4)
    st8, ct8 = sincos(th8)

    y4 = d.A + d.R11 * ct4 * c₁ - d.R12 * st4 * s₁

    y8 = d.A + d.R21 * ct8 * c₂ - d.R22 * st8 * s₂

    return muladd((y8 - y4) / 2, xi1, (y8 + y4) / 2)

  end

end

function mapy(d::annulus, u::Float64, v::Float64, k::Int)::Float64

  # p denotes information of kth patch
  p = d.pths[k]

  # Difference in c values
  hc = p.ck1 - p.ck0

  xi1 = hc * u / 2 + (p.ck1 + p.ck0) / 2

  # Difference in t values
  ht = p.tk1 - p.tk0

  xi2 = ht * v / 2 + (p.tk1 + p.tk0) / 2

  c₁ = d.RP.Cθ₁
  s₁ = d.RP.Sθ₁
  c₂ = d.RP.Cθ₂
  s₂ = d.RP.Sθ₂

  if p.reg == 1 

    th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
    th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

    th1 = muladd(Δth1, xi2, th₀1) 
    th5 = muladd(Δth5, xi2, th₀5) 

    st1, ct1 = sincos(th1)
    st5, ct5 = sincos(th5)

    y1 = d.B + d.R11 * ct1 * s₁ + d.R12 * st1 * c₁

    y5 = d.B + d.R21 * ct5 * s₂ + d.R22 * st5 * c₂

    return muladd((y1 - y5) / 2, xi1, (y1 + y5) / 2)

  elseif p.reg == 5

    th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
    th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

    th1 = muladd(Δth1, xi2, th₀1)
    th5 = muladd(Δth5, xi2, th₀5)

    st1, ct1 = sincos(th1)
    st5, ct5 = sincos(th5)

    y1 = d.B + d.R11 * ct1 * s₁ + d.R12 * st1 * c₁

    y5 = d.B + d.R21 * ct5 * s₂ + d.R22 * st5 * c₂

    return muladd((y5 - y1) / 2, xi1, (y5 + y1) / 2)

  elseif p.reg == 2

    th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
    th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

    th2 = muladd(Δth2, xi2, th₀2)
    th6 = muladd(Δth6, xi2, th₀6)

    st2, ct2 = sincos(th2)
    st6, ct6 = sincos(th6)

    y2 = d.B + d.R11 * ct2 * s₁ + d.R12 * st2 * c₁

    y6 = d.B + d.R21 * ct6 * s₂ + d.R22 * st6 * c₂

    return muladd((y2 - y6) / 2, xi1, (y2 + y6) / 2)

  elseif p.reg == 6

    th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
    th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

    th2 = muladd(Δth2, xi2, th₀2)
    th6 = muladd(Δth6, xi2, th₀6)

    st2, ct2 = sincos(th2)
    st6, ct6 = sincos(th6)

    y2 = d.B + d.R11 * ct2 * s₁ + d.R12 * st2 * c₁

    y6 = d.B + d.R21 * ct6 * s₂ + d.R22 * st6 * c₂

    return muladd((y6 - y2) / 2, xi1, (y6 + y2) / 2)

  elseif p.reg == 3

    th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
    th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

    th3 = muladd(Δth3, xi2, th₀3)
    th7 = muladd(Δth7, xi2, th₀7)

    st3, ct3 = sincos(th3)
    st7, ct7 = sincos(th7)

    y3 = d.B + d.R11 * ct3 * s₁ + d.R12 * st3 * c₁

    y7 = d.B + d.R21 * ct7 * s₂ + d.R22 * st7 * c₂

    return muladd((y3 - y7) / 2, xi1, (y3 + y7) / 2)

  elseif p.reg == 7

    th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
    th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

    th3 = muladd(Δth3, xi2, th₀3)
    th7 = muladd(Δth7, xi2, th₀7)

    st3, ct3 = sincos(th3)
    st7, ct7 = sincos(th7)

    y3 = d.B + d.R11 * ct3 * s₁ + d.R12 * st3 * c₁

    y7 = d.B + d.R21 * ct7 * s₂ + d.R22 * st7 * c₂

    return muladd((y7 - y3) / 2, xi1, (y7 + y3) / 2)

  elseif p.reg == 4

    th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
    th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

    th4 = muladd(Δth4, xi2, th₀4)
    th8 = muladd(Δth8, xi2, th₀8)

    st4, ct4 = sincos(th4)
    st8, ct8 = sincos(th8)

    y4 = d.B + d.R11 * ct4 * s₁ + d.R12 * st4 * c₁

    y8 = d.B + d.R21 * ct8 * s₂ + d.R22 * st8 * c₂

    return muladd((y4 - y8) / 2, xi1, (y4 + y8) / 2)

  else

    th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
    th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

    th4 = muladd(Δth4, xi2, th₀4)
    th8 = muladd(Δth8, xi2, th₀8)

    st4, ct4 = sincos(th4)
    st8, ct8 = sincos(th8)

    y4 = d.B + d.R11 * ct4 * s₁ + d.R12 * st4 * c₁

    y8 = d.B + d.R21 * ct8 * s₂ + d.R22 * st8 * c₂

    return muladd((y8 - y4) / 2, xi1, (y8 + y4) / 2)

  end

end

function mapxy(d::annulus, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64}

  @inbounds begin
    # p denotes information of kth patch
    p = d.pths[k]

    # Difference in c values
    hc = p.ck1 - p.ck0

    xi1 = hc * u / 2 + (p.ck1 + p.ck0) / 2

    # Difference in t values
    ht = p.tk1 - p.tk0

    xi2 = ht * v / 2 + (p.tk1 + p.tk0) / 2

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    if p.reg == 1

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      th1 = muladd(Δth1, xi2, th₀1)
      th5 = muladd(Δth5, xi2, th₀5)

      st1, ct1 = sincos(th1)
      st5, ct5 = sincos(th5)

      y1x = d.A + d.R11 * ct1 * c₁ - d.R12 * st1 * s₁

      y5x = d.A + d.R21 * ct5 * c₂ - d.R22 * st5 * s₂

      y1y = d.B + d.R11 * ct1 * s₁ + d.R12 * st1 * c₁

      y5y = d.B + d.R21 * ct5 * s₂ + d.R22 * st5 * c₂

      Zx = muladd((y1x - y5x) / 2, xi1, (y1x + y5x) / 2)

      Zy = muladd((y1y - y5y) / 2, xi1, (y1y + y5y) / 2)

      return Zx, Zy

    elseif p.reg == 5

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      th1 = muladd(Δth1, xi2, th₀1)
      th5 = muladd(Δth5, xi2, th₀5)

      st1, ct1 = sincos(th1)
      st5, ct5 = sincos(th5)

      y1x = d.A + d.R11 * ct1 * c₁ - d.R12 * st1 * s₁

      y5x = d.A + d.R21 * ct5 * c₂ - d.R22 * st5 * s₂

      y1y = d.B + d.R11 * ct1 * s₁ + d.R12 * st1 * c₁

      y5y = d.B + d.R21 * ct5 * s₂ + d.R22 * st5 * c₂

      Zx = muladd((y5x - y1x) / 2, xi1, (y5x + y1x) / 2)

      Zy = muladd((y5y - y1y) / 2, xi1, (y5y + y1y) / 2)

      return Zx, Zy

    elseif p.reg == 2

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      th2 = muladd(Δth2, xi2, th₀2)
      th6 = muladd(Δth6, xi2, th₀6)

      st2, ct2 = sincos(th2)
      st6, ct6 = sincos(th6)

      y2x = d.A + d.R11 * ct2 * c₁ - d.R12 * st2 * s₁

      y6x = d.A + d.R21 * ct6 * c₂ - d.R22 * st6 * s₂

      y2y = d.B + d.R11 * ct2 * s₁ + d.R12 * st2 * c₁

      y6y = d.B + d.R21 * ct6 * s₂ + d.R22 * st6 * c₂

      Zx = muladd((y2x - y6x) / 2, xi1, (y2x + y6x) / 2)

      Zy = muladd((y2y - y6y) / 2, xi1, (y2y + y6y) / 2)

      return Zx, Zy

    elseif p.reg == 6

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      th2 = muladd(Δth2, xi2, th₀2)
      th6 = muladd(Δth6, xi2, th₀6)

      st2, ct2 = sincos(th2)
      st6, ct6 = sincos(th6)

      y2x = d.A + d.R11 * ct2 * c₁ - d.R12 * st2 * s₁

      y6x = d.A + d.R21 * ct6 * c₂ - d.R22 * st6 * s₂

      y2y = d.B + d.R11 * ct2 * s₁ + d.R12 * st2 * c₁

      y6y = d.B + d.R21 * ct6 * s₂ + d.R22 * st6 * c₂

      Zx = muladd((y6x - y2x) / 2, xi1, (y6x + y2x) / 2)

      Zy = muladd((y6y - y2y) / 2, xi1, (y6y + y2y) / 2)

      return Zx, Zy

    elseif p.reg == 3

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      th3 = muladd(Δth3, xi2, th₀3)
      th7 = muladd(Δth7, xi2, th₀7)

      st3, ct3 = sincos(th3)
      st7, ct7 = sincos(th7)

      y3x = d.A + d.R11 * ct3 * c₁ - d.R12 * st3 * s₁

      y7x = d.A + d.R21 * ct7 * c₂ - d.R22 * st7 * s₂

      y3y = d.B + d.R11 * ct3 * s₁ + d.R12 * st3 * c₁

      y7y = d.B + d.R21 * ct7 * s₂ + d.R22 * st7 * c₂

      Zx = muladd((y3x - y7x) / 2, xi1, (y3x + y7x) / 2)

      Zy = muladd((y3y - y7y) / 2, xi1, (y3y + y7y) / 2)

      return Zx, Zy

    elseif p.reg == 7

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      th3 = muladd(Δth3, xi2, th₀3)
      th7 = muladd(Δth7, xi2, th₀7)

      st3, ct3 = sincos(th3)
      st7, ct7 = sincos(th7)

      y3x = d.A + d.R11 * ct3 * c₁ - d.R12 * st3 * s₁

      y7x = d.A + d.R21 * ct7 * c₂ - d.R22 * st7 * s₂

      y3y = d.B + d.R11 * ct3 * s₁ + d.R12 * st3 * c₁

      y7y = d.B + d.R21 * ct7 * s₂ + d.R22 * st7 * c₂

      Zx = muladd((y7x - y3x) / 2, xi1, (y7x + y3x) / 2)

      Zy = muladd((y7y - y3y) / 2, xi1, (y7y + y3y) / 2)

      return Zx, Zy

    elseif p.reg == 4

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      th4 = muladd(Δth4, xi2, th₀4)
      th8 = muladd(Δth8, xi2, th₀8)

      st4, ct4 = sincos(th4)
      st8, ct8 = sincos(th8)

      y4x = d.A + d.R11 * ct4 * c₁ - d.R12 * st4 * s₁

      y8x = d.A + d.R21 * ct8 * c₂ - d.R22 * st8 * s₂

      y4y = d.B + d.R11 * ct4 * s₁ + d.R12 * st4 * c₁

      y8y = d.B + d.R21 * ct8 * s₂ + d.R22 * st8 * c₂

      Zx = muladd((y4x - y8x) / 2, xi1, (y4x + y8x) / 2)

      Zy = muladd((y4y - y8y) / 2, xi1, (y4y + y8y) / 2)

      return Zx, Zy

    else

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      th4 = muladd(Δth4, xi2, th₀4)
      th8 = muladd(Δth8, xi2, th₀8)

      st4, ct4 = sincos(th4)
      st8, ct8 = sincos(th8)

      y4x = d.A + d.R11 * ct4 * c₁ - d.R12 * st4 * s₁

      y8x = d.A + d.R21 * ct8 * c₂ - d.R22 * st8 * s₂

      y4y = d.B + d.R11 * ct4 * s₁ + d.R12 * st4 * c₁

      y8y = d.B + d.R21 * ct8 * s₂ + d.R22 * st8 * c₂

      Zx = muladd((y8x - y4x) / 2, xi1, (y8x + y4x) / 2)

      Zy = muladd((y8y - y4y) / 2, xi1, (y8y + y4y) / 2)

      return Zx, Zy

    end
  end
end

function mapxy!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
  d::annulus, u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)
  @inbounds begin
    # p denotes information of kth patch
    p = d.pths[k]

    hc = p.ck1 - p.ck0
    ht = p.tk1 - p.tk0
    αu = 0.5 * hc
    αv = 0.5 * ht
    βu = p.ck0 + αu
    βv = p.tk0 + αv

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11s1, R11c1 = d.R11 * s₁, d.R11 * c₁
    R12s1, R12c1 = d.R12 * s₁, d.R12 * c₁

    R21s2, R21c2 = d.R21 * s₂, d.R21 * c₂
    R22s2, R22c2 = d.R22 * s₂, d.R22 * c₂

    if p.reg == 1

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      @inbounds for I in eachindex(Zx, Zy, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th1 = muladd(Δth1, xi2, th₀1)
        th5 = muladd(Δth5, xi2, th₀5)

        st1, ct1 = sincos(th1)
        st5, ct5 = sincos(th5)

        y1x = muladd(-st1, R12s1, muladd(ct1,  R11c1, d.A))
        y5x = muladd(-st5, R22s2, muladd(ct5,  R21c2, d.A))

        y1y = muladd(st1, R12c1, muladd(ct1, R11s1, d.B))
        y5y = muladd(st5, R22c2, muladd(ct5, R21s2, d.B))

        t = muladd(0.5, xi1, 0.5)

        Zx[I] = muladd(t, y1x - y5x, y5x)
        Zy[I] = muladd(t, y1y - y5y, y5y)
      end

    elseif p.reg == 5

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      @inbounds for I in eachindex(Zx, Zy, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th1 = muladd(Δth1, xi2, th₀1)
        th5 = muladd(Δth5, xi2, th₀5)

        st1, ct1 = sincos(th1)
        st5, ct5 = sincos(th5)

        y1x = muladd(-st1, R12s1, muladd(ct1,  R11c1, d.A))
        y5x = muladd(-st5, R22s2, muladd(ct5,  R21c2, d.A))

        y1y = muladd(st1, R12c1, muladd(ct1, R11s1, d.B))
        y5y = muladd(st5, R22c2, muladd(ct5, R21s2, d.B))

        t = muladd(0.5, xi1, 0.5)

        Zx[I] = muladd(t, y5x - y1x, y1x)
        Zy[I] = muladd(t, y5y - y1y, y1y)
      end

    elseif p.reg == 2

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      @inbounds for I in eachindex(Zx, Zy, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th2 = muladd(Δth2, xi2, th₀2)
        th6 = muladd(Δth6, xi2, th₀6)

        st2, ct2 = sincos(th2)
        st6, ct6 = sincos(th6)

        y2x = muladd(-st2, R12s1, muladd(ct2,  R11c1, d.A))
        y6x = muladd(-st6, R22s2, muladd(ct6,  R21c2, d.A))

        y2y = muladd(st2, R12c1, muladd(ct2, R11s1, d.B))
        y6y = muladd(st6, R22c2, muladd(ct6, R21s2, d.B))

        t = muladd(0.5, xi1, 0.5)

        Zx[I] = muladd(t, y2x - y6x, y6x)
        Zy[I] = muladd(t, y2y - y6y, y6y)
      end

    elseif p.reg == 6

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      @inbounds for I in eachindex(Zx, Zy, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th2 = muladd(Δth2, xi2, th₀2)
        th6 = muladd(Δth6, xi2, th₀6)

        st2, ct2 = sincos(th2)
        st6, ct6 = sincos(th6)

        y2x = muladd(-st2, R12s1, muladd(ct2,  R11c1, d.A))
        y6x = muladd(-st6, R22s2, muladd(ct6,  R21c2, d.A))

        y2y = muladd(st2, R12c1, muladd(ct2, R11s1, d.B))
        y6y = muladd(st6, R22c2, muladd(ct6, R21s2, d.B))

        t = muladd(0.5, xi1, 0.5)

        Zx[I] = muladd(t, y6x - y2x, y2x)
        Zy[I] = muladd(t, y6y - y2y, y2y)
      end

    elseif p.reg == 3

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      @inbounds for I in eachindex(Zx, Zy, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th3 = muladd(Δth3, xi2, th₀3)
        th7 = muladd(Δth7, xi2, th₀7)

        st3, ct3 = sincos(th3)
        st7, ct7 = sincos(th7)

        y3x = muladd(-st3, R12s1, muladd(ct3,  R11c1, d.A))
        y7x = muladd(-st7, R22s2, muladd(ct7,  R21c2, d.A))

        y3y = muladd(st3, R12c1, muladd(ct3, R11s1, d.B))
        y7y = muladd(st7, R22c2, muladd(ct7, R21s2, d.B))

        t = muladd(0.5, xi1, 0.5)

        Zx[I] = muladd(t, y3x - y7x, y7x)
        Zy[I] = muladd(t, y3y - y7y, y7y)
      end

    elseif p.reg == 7

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      @inbounds for I in eachindex(Zx, Zy, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th3 = muladd(Δth3, xi2, th₀3)
        th7 = muladd(Δth7, xi2, th₀7)

        st3, ct3 = sincos(th3)
        st7, ct7 = sincos(th7)

        y3x = muladd(-st3, R12s1, muladd(ct3,  R11c1, d.A))
        y7x = muladd(-st7, R22s2, muladd(ct7,  R21c2, d.A))

        y3y = muladd(st3, R12c1, muladd(ct3, R11s1, d.B))
        y7y = muladd(st7, R22c2, muladd(ct7, R21s2, d.B))

        t = muladd(0.5, xi1, 0.5)

        Zx[I] = muladd(t, y7x - y3x, y3x)
        Zy[I] = muladd(t, y7y - y3y, y3y)
      end

    elseif p.reg == 4

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      @inbounds for I in eachindex(Zx, Zy, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th4 = muladd(Δth4, xi2, th₀4)
        th8 = muladd(Δth8, xi2, th₀8)

        st4, ct4 = sincos(th4)
        st8, ct8 = sincos(th8)

        y4x = muladd(-st4, R12s1, muladd(ct4,  R11c1, d.A))
        y8x = muladd(-st8, R22s2, muladd(ct8,  R21c2, d.A))

        y4y = muladd(st4, R12c1, muladd(ct4, R11s1, d.B))
        y8y = muladd(st8, R22c2, muladd(ct8, R21s2, d.B))

        t = muladd(0.5, xi1, 0.5)

        Zx[I] = muladd(t, y4x - y8x, y8x)
        Zy[I] = muladd(t, y4y - y8y, y8y)
      end

    else

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      @inbounds for I in eachindex(Zx, Zy, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th4 = muladd(Δth4, xi2, th₀4)
        th8 = muladd(Δth8, xi2, th₀8)

        st4, ct4 = sincos(th4)
        st8, ct8 = sincos(th8)

        y4x = muladd(-st4, R12s1, muladd(ct4,  R11c1, d.A))
        y8x = muladd(-st8, R22s2, muladd(ct8,  R21c2, d.A))

        y4y = muladd(st4, R12c1, muladd(ct4, R11s1, d.B))
        y8y = muladd(st8, R22c2, muladd(ct8, R21s2, d.B))

        t = muladd(0.5, xi1, 0.5)

        Zx[I] = muladd(t, y8x - y4x, y4x)
        Zy[I] = muladd(t, y8y - y4y, y4y)
      end

    end
  end

  return nothing
end

function draw(d::annulus, flag=nothing)
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

  # define 8 colors
  colors = (RGBf(1.0000, 0.0000, 0.0000),  # Red
            RGBf(1.0000, 0.5000, 0.2500),  # Orange
            RGBf(0.5800, 0.0000, 0.8200),  # Purple
            RGBf(0.0000, 1.0000, 0.0000),  # Green
            RGBf(1.0000, 0.0000, 0.0000),  # Red
            RGBf(1.0000, 0.5000, 0.2500),  # Orange
            RGBf(0.5800, 0.0000, 0.8200),  # Purple
            RGBf(0.0000, 1.0000, 0.0000),  # Green
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

  display(GLMakie.Screen(),fig)
  return (fig, ax)
  #savefig(plt, "myplot.svg")
end

#-----------------------
"""
  gamx(d::annulus, t::Float64, k::Int)
  gamx!(out::StridedArray{Float64}, d::annulus, t::StridedArray{Float64}, k)
First (x) coordinate of the boundary parametrization γ(t).

- `gamx(d, t, k)` computes the **right boundary** of patch `k`, where `t ∈ [-1,1]`.

`t` may be a scalar or an array; the return has the same shape.
"""
function gamx(d::annulus, t::Float64, k::Int)::Float64
  @inbounds begin
    p = d.pths[k]

    # Map t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
    xi2 = muladd(0.5 * (p.tk1 - p.tk0), t, 0.5 * (p.tk0 + p.tk1))

    c1 = d.RP.Cθ₁
    s1 = d.RP.Sθ₁
    c2 = d.RP.Cθ₂
    s2 = d.RP.Sθ₂

    R11c1 = d.R11 * c1
    R12s1 = d.R12 * s1
    R21c2 = d.R21 * c2
    R22s2 = d.R22 * s2

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(-st, R12s1, muladd(ct, R11c1, d.A))

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(-st, R22s2, muladd(ct, R21c2, d.A))

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(-st, R12s1, muladd(ct, R11c1, d.A))

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(-st, R22s2, muladd(ct, R21c2, d.A))

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(-st, R12s1, muladd(ct, R11c1, d.A))

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(-st, R22s2, muladd(ct, R21c2, d.A))

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(-st, R12s1, muladd(ct, R11c1, d.A))

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(-st, R22s2, muladd(ct, R21c2, d.A))

    else
      throw(ArgumentError("gamx is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end
end

function gamx!(out::StridedArray{Float64}, d::annulus, t::StridedArray{Float64}, k::Int)
  @inbounds begin
    p = d.pths[k]

    ht = p.tk1 - p.tk0
    αv = 0.5 * ht
    βv = p.tk0 + αv

    c1 = d.RP.Cθ₁
    s1 = d.RP.Sθ₁
    c2 = d.RP.Cθ₂
    s2 = d.RP.Sθ₂

    R11c1 = d.R11 * c1
    R12s1 = d.R12 * s1
    R21c2 = d.R21 * c2
    R22s2 = d.R22 * s2

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
      end

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
      end

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
      end

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
      end

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
      end

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
      end

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
      end

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
      end

    else
      throw(ArgumentError("gamx! is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end

  return nothing
end

"""
    gamy(d::annulus, t::Float64, k::Int)
    gamy!(out::StridedArray{Float64}, d::annulus, t::StridedArray{Float64}, k)

Second (y) coordinate of the boundary parametrization γ(t).

- `gamy(d, t, k)` computes the right boundary of patch `k`, where `t ∈ [-1,1]`.

`t` can be a scalar `Float64` or an array; the return has the same shape.
"""
function gamy(d::annulus, t::Float64, k::Int)::Float64
  @inbounds begin
    p = d.pths[k]

    # Map t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
    xi2 = muladd(0.5 * (p.tk1 - p.tk0), t, 0.5 * (p.tk0 + p.tk1))

    c1 = d.RP.Cθ₁
    s1 = d.RP.Sθ₁
    c2 = d.RP.Cθ₂
    s2 = d.RP.Sθ₂

    R11s1 = d.R11 * s1
    R12c1 = d.R12 * c1
    R21s2 = d.R21 * s2
    R22c2 = d.R22 * c2

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(st, R12c1, muladd(ct, R11s1, d.B))

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(st, R22c2, muladd(ct, R21s2, d.B))

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(st, R12c1, muladd(ct, R11s1, d.B))

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(st, R22c2, muladd(ct, R21s2, d.B))

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(st, R12c1, muladd(ct, R11s1, d.B))

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(st, R22c2, muladd(ct, R21s2, d.B))

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(st, R12c1, muladd(ct, R11s1, d.B))

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      return muladd(st, R22c2, muladd(ct, R21s2, d.B))

    else
      throw(ArgumentError("gamy is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end
end

function gamy!(out::StridedArray{Float64}, d::annulus, t::StridedArray{Float64}, k::Int)
  @inbounds begin
    p = d.pths[k]

    # Map reference t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
    αv = 0.5 * (p.tk1 - p.tk0)
    βv = 0.5 * (p.tk0 + p.tk1)

    c1 = d.RP.Cθ₁
    s1 = d.RP.Sθ₁
    c2 = d.RP.Cθ₂
    s2 = d.RP.Sθ₂

    R11s1 = d.R11 * s1
    R12c1 = d.R12 * c1
    R21s2 = d.R21 * s2
    R22c2 = d.R22 * c2

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(st, R12c1, muladd(ct, R11s1, d.B))
      end

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(st, R22c2, muladd(ct, R21s2, d.B))
      end

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(st, R12c1, muladd(ct, R11s1, d.B))
      end

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(st, R22c2, muladd(ct, R21s2, d.B))
      end

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(st, R12c1, muladd(ct, R11s1, d.B))
      end

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(st, R22c2, muladd(ct, R21s2, d.B))
      end

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(st, R12c1, muladd(ct, R11s1, d.B))
      end

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
      @inbounds for I in eachindex(out, t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[I] = muladd(st, R22c2, muladd(ct, R21s2, d.B))
      end

    else
      throw(ArgumentError("gamy! is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end

  return nothing
end

function gam(d::annulus, t::Float64, k::Int)::Tuple{Float64,Float64}
  @inbounds begin
    # p denotes information of kth patch
    p = d.pths[k]

    # Map reference t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
    xi2 = muladd(0.5 * (p.tk1 - p.tk0), t, 0.5 * (p.tk0 + p.tk1))

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11s1, R11c1 = d.R11 * s₁, d.R11 * c₁
    R12s1, R12c1 = d.R12 * s₁, d.R12 * c₁

    R21s2, R21c2 = d.R21 * s₂, d.R21 * c₂
    R22s2, R22c2 = d.R22 * s₂, d.R22 * c₂

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      x = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
      y = muladd(st, R12c1, muladd(ct, R11s1, d.B))
      return x, y

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      x = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
      y = muladd(st, R22c2, muladd(ct, R21s2, d.B))
      return x, y

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      x = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
      y = muladd(st, R12c1, muladd(ct, R11s1, d.B))
      return x, y

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      x = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
      y = muladd(st, R22c2, muladd(ct, R21s2, d.B))
      return x, y

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      x = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
      y = muladd(st, R12c1, muladd(ct, R11s1, d.B))
      return x, y

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      x = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
      y = muladd(st, R22c2, muladd(ct, R21s2, d.B))
      return x, y

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      x = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
      y = muladd(st, R12c1, muladd(ct, R11s1, d.B))
      return x, y

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      x = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
      y = muladd(st, R22c2, muladd(ct, R21s2, d.B))
      return x, y

    else
      throw(ArgumentError("gam is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end
end

function gam!(out::Vector{Float64}, d::annulus, t::Float64, k::Int)
  @inbounds begin
    # p denotes information of kth patch
    p = d.pths[k]

    # Map reference t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
    xi2 = muladd(0.5 * (p.tk1 - p.tk0), t, 0.5 * (p.tk0 + p.tk1))

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11s1, R11c1 = d.R11 * s₁, d.R11 * c₁
    R12s1, R12c1 = d.R12 * s₁, d.R12 * c₁

    R21s2, R21c2 = d.R21 * s₂, d.R21 * c₂
    R22s2, R22c2 = d.R22 * s₂, d.R22 * c₂

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      out[1] = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
      out[2] = muladd(st, R12c1, muladd(ct, R11s1, d.B))

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      out[1] = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
      out[2] = muladd(st, R22c2, muladd(ct, R21s2, d.B))

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      out[1] = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
      out[2] = muladd(st, R12c1, muladd(ct, R11s1, d.B))

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      out[1] = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
      out[2] = muladd(st, R22c2, muladd(ct, R21s2, d.B))

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      out[1] = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
      out[2] = muladd(st, R12c1, muladd(ct, R11s1, d.B))

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      out[1] = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
      out[2] = muladd(st, R22c2, muladd(ct, R21s2, d.B))

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      out[1] = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
      out[2] = muladd(st, R12c1, muladd(ct, R11s1, d.B))

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
      th = muladd(Δth, xi2, th0)
      st, ct = sincos(th)
      out[1] = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
      out[2] = muladd(st, R22c2, muladd(ct, R21s2, d.B))

    else
      throw(ArgumentError("gam! is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end

  return nothing
end

function gam!(out::Matrix{Float64}, d::annulus, t::Vector{Float64}, k::Int)
  @inbounds begin
    p = d.pths[k]
    #@assert size(out,1) == 2
    #@assert size(out,2) == length(t)

    # Map reference t ∈ [-1,1] → ξ₂ ∈ [tk0, tk1]
    αv = 0.5 * (p.tk1 - p.tk0)
    βv = 0.5 * (p.tk0 + p.tk1)

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11s1, R11c1 = d.R11 * s₁, d.R11 * c₁
    R12s1, R12c1 = d.R12 * s₁, d.R12 * c₁

    R21s2, R21c2 = d.R21 * s₂, d.R21 * c₂
    R22s2, R22c2 = d.R22 * s₂, d.R22 * c₂

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      @inbounds for I in eachindex(t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[1, I] = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
        out[2, I] = muladd(st, R12c1, muladd(ct, R11s1, d.B))
      end

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      @inbounds for I in eachindex(t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[1, I] = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
        out[2, I] = muladd(st, R22c2, muladd(ct, R21s2, d.B))
      end

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      @inbounds for I in eachindex(t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[1, I] = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
        out[2, I] = muladd(st, R12c1, muladd(ct, R11s1, d.B))
      end

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      @inbounds for I in eachindex(t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[1, I] = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
        out[2, I] = muladd(st, R22c2, muladd(ct, R21s2, d.B))
      end

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      @inbounds for I in eachindex(t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[1, I] = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
        out[2, I] = muladd(st, R12c1, muladd(ct, R11s1, d.B))
      end

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      @inbounds for I in eachindex(t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[1, I] = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
        out[2, I] = muladd(st, R22c2, muladd(ct, R21s2, d.B))
      end

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      @inbounds for I in eachindex(t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[1, I] = muladd(-st, R12s1, muladd(ct, R11c1, d.A))
        out[2, I] = muladd(st, R12c1, muladd(ct, R11s1, d.B))
      end

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
      @inbounds for I in eachindex(t)
        xi2 = muladd(αv, t[I], βv)
        th = muladd(Δth, xi2, th0)
        st, ct = sincos(th)
        out[1, I] = muladd(-st, R22s2, muladd(ct, R21c2, d.A))
        out[2, I] = muladd(st, R22c2, muladd(ct, R21s2, d.B))
      end

    else
      throw(ArgumentError("gam! is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end

  return nothing
end

function drawbd(d::annulus)
  # Number of sampling points
  L = 33
  t = collect(range(-1, 1, L))

  # define 8 colors
  clr = (RGBf(1.0000, 0.0000, 0.0000),  # Red
    RGBf(1.0000, 0.5000, 0.2500),  # Orange
    RGBf(0.5800, 0.0000, 0.8200),  # Purple
    RGBf(0.0000, 1.0000, 0.0000),  # Green
    RGBf(1.0000, 0.0000, 0.0000),  # Red
    RGBf(1.0000, 0.5000, 0.2500),  # Orange
    RGBf(0.5800, 0.0000, 0.8200),  # Purple
    RGBf(0.0000, 1.0000, 0.0000),  # Green
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
    dgamx(d::annulus, t::Float64, k::Int) -> Float64

Derivative of the first coordinate of the boundary parametrization.

- With `k`: derivative of the **right boundary** of patch `k`
  (i.e., `∂/∂t mapx(d, 1, t, k)` for `t ∈ [-1,1]`).
"""
function dgamx(d::annulus, t::Float64, k::Int)::Float64
  @inbounds begin
    p = d.pths[k]

    ht = p.tk1 - p.tk0
    αv = 0.5 * ht
    βv = 0.5 * (p.tk0 + p.tk1)

    τ = muladd(αv, t, βv)

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11c1 = d.R11 * c₁
    R12s1 = d.R12 * s₁
    R21c2 = d.R21 * c₂
    R22s2 = d.R22 * s₂

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return (-Δth) * muladd(st, R11c1, ct * R12s1) * αv

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return (-Δth) * muladd(st, R21c2, ct * R22s2) * αv

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return (-Δth) * muladd(st, R11c1, ct * R12s1) * αv

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return (-Δth) * muladd(st, R21c2, ct * R22s2) * αv

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return (-Δth) * muladd(st, R11c1, ct * R12s1) * αv

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return (-Δth) * muladd(st, R21c2, ct * R22s2) * αv

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return (-Δth) * muladd(st, R11c1, ct * R12s1) * αv

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return (-Δth) * muladd(st, R21c2, ct * R22s2) * αv

    else
      throw(ArgumentError("dgamx is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end
end

function dgamx!(out::StridedArray{Float64}, d::annulus, t::StridedArray{Float64}, k::Int)
  @inbounds begin
    p = d.pths[k]

    ht = p.tk1 - p.tk0
    αv = 0.5 * ht
    βv = 0.5 * (p.tk0 + p.tk1)

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11c1 = d.R11 * c₁
    R12s1 = d.R12 * s₁
    R21c2 = d.R21 * c₂
    R22s2 = d.R22 * s₂

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      fac = (-Δth) * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(st, R11c1, ct * R12s1)
      end

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      fac = (-Δth) * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(st, R21c2, ct * R22s2)
      end

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      fac = (-Δth) * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(st, R11c1, ct * R12s1)
      end

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      fac = (-Δth) * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(st, R21c2, ct * R22s2)
      end

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      fac = (-Δth) * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(st, R11c1, ct * R12s1)
      end

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      fac = (-Δth) * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(st, R21c2, ct * R22s2)
      end

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      fac = (-Δth) * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(st, R11c1, ct * R12s1)
      end

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
      fac = (-Δth) * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(st, R21c2, ct * R22s2)
      end

    else
      throw(ArgumentError("dgamx! is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end
  return nothing
end

"""
    dgamy(d::annulus, t::Float64, k::Int) -> Float64

Derivative of the second coordinate of the boundary parametrization γ(t).

- With `k`: derivative of the **right boundary** of patch `k`
  (`∂/∂t mapy(d, 1, t, k)` for `t ∈ [-1,1]`).
"""
function dgamy(d::annulus, t::Float64, k::Int)::Float64
  @inbounds begin
    p = d.pths[k]

    ht = p.tk1 - p.tk0
    αv = 0.5 * ht
    βv = 0.5 * (p.tk0 + p.tk1)

    τ = muladd(αv, t, βv)

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11s1 = d.R11 * s₁
    R12c1 = d.R12 * c₁
    R21s2 = d.R21 * s₂
    R22c2 = d.R22 * c₂

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return Δth * muladd(-st, R11s1, ct * R12c1) * αv

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return Δth * muladd(-st, R21s2, ct * R22c2) * αv

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return Δth * muladd(-st, R11s1, ct * R12c1) * αv

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return Δth * muladd(-st, R21s2, ct * R22c2) * αv

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return Δth * muladd(-st, R11s1, ct * R12c1) * αv

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return Δth * muladd(-st, R21s2, ct * R22c2) * αv

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return Δth * muladd(-st, R11s1, ct * R12c1) * αv

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      return Δth * muladd(-st, R21s2, ct * R22c2) * αv

    else
      throw(ArgumentError("dgamy is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end
end

function dgamy!(out::StridedArray{Float64}, d::annulus, t::StridedArray{Float64}, k::Int)
  @inbounds begin
    p = d.pths[k]

    ht = p.tk1 - p.tk0
    αv = 0.5 * ht
    βv = 0.5 * (p.tk0 + p.tk1)

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11s1 = d.R11 * s₁
    R12c1 = d.R12 * c₁
    R21s2 = d.R21 * s₂
    R22c2 = d.R22 * c₂

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      fac = Δth * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(-st, R11s1, ct * R12c1)
      end

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      fac = Δth * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(-st, R21s2, ct * R22c2)
      end

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      fac = Δth * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(-st, R11s1, ct * R12c1)
      end

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      fac = Δth * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(-st, R21s2, ct * R22c2)
      end

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      fac = Δth * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(-st, R11s1, ct * R12c1)
      end

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      fac = Δth * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(-st, R21s2, ct * R22c2)
      end

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      fac = Δth * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(-st, R11s1, ct * R12c1)
      end

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
      fac = Δth * αv
      for I in eachindex(out, t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[I] = fac * muladd(-st, R21s2, ct * R22c2)
      end

    else
      throw(ArgumentError("dgamy! is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end
  return nothing
end

function dgam(d::annulus, t::Float64, k::Int)::Tuple{Float64,Float64}
  @inbounds begin
    p = d.pths[k]

    ht = p.tk1 - p.tk0
    αv = 0.5 * ht
    βv = 0.5 * (p.tk0 + p.tk1)

    τ = muladd(αv, t, βv)

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11c1 = d.R11 * c₁
    R12s1 = d.R12 * s₁
    R11s1 = d.R11 * s₁
    R12c1 = d.R12 * c₁

    R21c2 = d.R21 * c₂
    R22s2 = d.R22 * s₂
    R21s2 = d.R21 * s₂
    R22c2 = d.R22 * c₂

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      dx = (-Δth) * muladd(st, R11c1, ct * R12s1) * αv
      dy = (Δth) * muladd(-st, R11s1, ct * R12c1) * αv
      return dx, dy

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      dx = (-Δth) * muladd(st, R21c2, ct * R22s2) * αv
      dy = (Δth) * muladd(-st, R21s2, ct * R22c2) * αv
      return dx, dy

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      dx = (-Δth) * muladd(st, R11c1, ct * R12s1) * αv
      dy = (Δth) * muladd(-st, R11s1, ct * R12c1) * αv
      return dx, dy

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      dx = (-Δth) * muladd(st, R21c2, ct * R22s2) * αv
      dy = (Δth) * muladd(-st, R21s2, ct * R22c2) * αv
      return dx, dy

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      dx = (-Δth) * muladd(st, R11c1, ct * R12s1) * αv
      dy = (Δth) * muladd(-st, R11s1, ct * R12c1) * αv
      return dx, dy

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      dx = (-Δth) * muladd(st, R21c2, ct * R22s2) * αv
      dy = (Δth) * muladd(-st, R21s2, ct * R22c2) * αv
      return dx, dy

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      dx = (-Δth) * muladd(st, R11c1, ct * R12s1) * αv
      dy = (Δth) * muladd(-st, R11s1, ct * R12c1) * αv
      return dx, dy

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      dx = (-Δth) * muladd(st, R21c2, ct * R22s2) * αv
      dy = (Δth) * muladd(-st, R21s2, ct * R22c2) * αv
      return dx, dy

    else
      throw(ArgumentError("dgam is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end
end

function dgam!(out::Vector{Float64}, d::annulus, t::Float64, k::Int)
  @inbounds begin
    p = d.pths[k]

    ht = p.tk1 - p.tk0
    αv = 0.5 * ht
    βv = 0.5 * (p.tk0 + p.tk1)

    τ = muladd(αv, t, βv)

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11c1 = d.R11 * c₁
    R12s1 = d.R12 * s₁
    R11s1 = d.R11 * s₁
    R12c1 = d.R12 * c₁

    R21c2 = d.R21 * c₂
    R22s2 = d.R22 * s₂
    R21s2 = d.R21 * s₂
    R22c2 = d.R22 * c₂

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      out[1] = (-Δth) * muladd(st, R11c1, ct * R12s1) * αv
      out[2] = (Δth) * muladd(-st, R11s1, ct * R12c1) * αv

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      out[1] = (-Δth) * muladd(st, R21c2, ct * R22s2) * αv
      out[2] = (Δth) * muladd(-st, R21s2, ct * R22c2) * αv

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      out[1] = (-Δth) * muladd(st, R11c1, ct * R12s1) * αv
      out[2] = (Δth) * muladd(-st, R11s1, ct * R12c1) * αv

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      out[1] = (-Δth) * muladd(st, R21c2, ct * R22s2) * αv
      out[2] = (Δth) * muladd(-st, R21s2, ct * R22c2) * αv

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      out[1] = (-Δth) * muladd(st, R11c1, ct * R12s1) * αv
      out[2] = (Δth) * muladd(-st, R11s1, ct * R12c1) * αv

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      out[1] = (-Δth) * muladd(st, R21c2, ct * R22s2) * αv
      out[2] = (Δth) * muladd(-st, R21s2, ct * R22c2) * αv

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      out[1] = (-Δth) * muladd(st, R11c1, ct * R12s1) * αv
      out[2] = (Δth) * muladd(-st, R11s1, ct * R12c1) * αv

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)
      out[1] = (-Δth) * muladd(st, R21c2, ct * R22s2) * αv
      out[2] = (Δth) * muladd(-st, R21s2, ct * R22c2) * αv

    else
      throw(ArgumentError("dgam is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end

  return nothing
end

function dgam!(out::Matrix{Float64}, d::annulus, t::Vector{Float64}, k::Int)
  @inbounds begin
    p = d.pths[k]

    ht = p.tk1 - p.tk0
    αv = 0.5 * ht
    βv = 0.5 * (p.tk0 + p.tk1)

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11c1 = d.R11 * c₁
    R12s1 = d.R12 * s₁
    R11s1 = d.R11 * s₁
    R12c1 = d.R12 * c₁

    R21c2 = d.R21 * c₂
    R22s2 = d.R22 * s₂
    R21s2 = d.R21 * s₂
    R22c2 = d.R22 * c₂

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      facx = (-Δth) * αv
      facy = (Δth) * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = facx * muladd(st, R11c1, ct * R12s1)
        out[2, I] = facy * muladd(-st, R11s1, ct * R12c1)
      end

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      facx = (-Δth) * αv
      facy = (Δth) * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = facx * muladd(st, R21c2, ct * R22s2)
        out[2, I] = facy * muladd(-st, R21s2, ct * R22c2)
      end

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      facx = (-Δth) * αv
      facy = (Δth) * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = facx * muladd(st, R11c1, ct * R12s1)
        out[2, I] = facy * muladd(-st, R11s1, ct * R12c1)
      end

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      facx = (-Δth) * αv
      facy = (Δth) * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = facx * muladd(st, R21c2, ct * R22s2)
        out[2, I] = facy * muladd(-st, R21s2, ct * R22c2)
      end

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      facx = (-Δth) * αv
      facy = (Δth) * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = facx * muladd(st, R11c1, ct * R12s1)
        out[2, I] = facy * muladd(-st, R11s1, ct * R12c1)
      end

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      facx = (-Δth) * αv
      facy = (Δth) * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = facx * muladd(st, R21c2, ct * R22s2)
        out[2, I] = facy * muladd(-st, R21s2, ct * R22c2)
      end

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      facx = (-Δth) * αv
      facy = (Δth) * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = facx * muladd(st, R11c1, ct * R12s1)
        out[2, I] = facy * muladd(-st, R11s1, ct * R12c1)
      end

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
      facx = (-Δth) * αv
      facy = (Δth) * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = facx * muladd(st, R21c2, ct * R22s2)
        out[2, I] = facy * muladd(-st, R21s2, ct * R22c2)
      end

    else
      throw(ArgumentError("dgam! is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end
  return nothing
end

function gamp!(out::Matrix{Float64}, d::annulus, t::Vector{Float64}, k::Int)
  # out[1,:] = dy, out[2,:] = -dx
  @inbounds begin
    p = d.pths[k]

    ht = p.tk1 - p.tk0
    αv = 0.5 * ht
    βv = 0.5 * (p.tk0 + p.tk1)

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11c1 = d.R11 * c₁
    R12s1 = d.R12 * s₁
    R11s1 = d.R11 * s₁
    R12c1 = d.R12 * c₁

    R21c2 = d.R21 * c₂
    R22s2 = d.R22 * s₂
    R21s2 = d.R21 * s₂
    R22c2 = d.R22 * c₂

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = fac * muladd(-st, R11s1, ct * R12c1)  # dy
        out[2, I] = fac * muladd(st, R11c1, ct * R12s1)  # -dx
      end

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = fac * muladd(-st, R21s2, ct * R22c2)
        out[2, I] = fac * muladd(st, R21c2, ct * R22s2)
      end

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = fac * muladd(-st, R11s1, ct * R12c1)
        out[2, I] = fac * muladd(st, R11c1, ct * R12s1)
      end

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = fac * muladd(-st, R21s2, ct * R22c2)
        out[2, I] = fac * muladd(st, R21c2, ct * R22s2)
      end

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = fac * muladd(-st, R11s1, ct * R12c1)
        out[2, I] = fac * muladd(st, R11c1, ct * R12s1)
      end

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = fac * muladd(-st, R21s2, ct * R22c2)
        out[2, I] = fac * muladd(st, R21c2, ct * R22s2)
      end

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = fac * muladd(-st, R11s1, ct * R12c1)
        out[2, I] = fac * muladd(st, R11c1, ct * R12s1)
      end

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)
        out[1, I] = fac * muladd(-st, R21s2, ct * R22c2)
        out[2, I] = fac * muladd(st, R21c2, ct * R22s2)
      end

    else
      throw(ArgumentError("gamp! is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end
  return nothing
end

function nu!(out::Vector{Float64}, d::annulus, t::Float64, k::Int)
  @inbounds begin
    p = d.pths[k]

    ht = p.tk1 - p.tk0
    αv = 0.5 * ht
    βv = 0.5 * (p.tk0 + p.tk1)

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11c1 = d.R11 * c₁
    R12s1 = d.R12 * s₁
    R11s1 = d.R11 * s₁
    R12c1 = d.R12 * c₁

    R21c2 = d.R21 * c₂
    R22s2 = d.R22 * s₂
    R21s2 = d.R21 * s₂
    R22c2 = d.R22 * c₂

    τ = muladd(αv, t, βv)

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      fac = Δth * αv
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)

      dy = fac * muladd(-st, R11s1, ct * R12c1)
      dx = fac * muladd(-st, R11c1, ct * R12s1)

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      fac = Δth * αv
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)

      dy = fac * muladd(-st, R21s2, ct * R22c2)
      dx = fac * muladd(-st, R21c2, ct * R22s2)

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      fac = Δth * αv
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)

      dy = fac * muladd(-st, R11s1, ct * R12c1)
      dx = fac * muladd(-st, R11c1, ct * R12s1)

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      fac = Δth * αv
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)

      dy = fac * muladd(-st, R21s2, ct * R22c2)
      dx = fac * muladd(-st, R21c2, ct * R22s2)

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      fac = Δth * αv
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)

      dy = fac * muladd(-st, R11s1, ct * R12c1)
      dx = fac * muladd(-st, R11c1, ct * R12s1)

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      fac = Δth * αv
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)

      dy = fac * muladd(-st, R21s2, ct * R22c2)
      dx = fac * muladd(-st, R21c2, ct * R22s2)

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2π
      fac = Δth * αv
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)

      dy = fac * muladd(-st, R11s1, ct * R12c1)
      dx = fac * muladd(-st, R11c1, ct * R12s1)

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2π
      fac = Δth * αv
      th = muladd(Δth, τ, th0)
      st, ct = sincos(th)

      dy = fac * muladd(-st, R21s2, ct * R22c2)
      dx = fac * muladd(-st, R21c2, ct * R22s2)

    else
      throw(ArgumentError("nu! is defined only for regions 1–8; got reg=$(p.reg)"))
    end

    S = hypot(dx, dy)
    out[1] = dy / S
    out[2] = dx / S
  end
  return nothing
end

function nu!(out::Matrix{Float64}, d::annulus, t::Vector{Float64}, k::Int)
  @inbounds begin
    p = d.pths[k]

    ht = p.tk1 - p.tk0
    αv = 0.5 * ht
    βv = 0.5 * (p.tk0 + p.tk1)

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11c1 = d.R11 * c₁
    R12s1 = d.R12 * s₁
    R11s1 = d.R11 * s₁
    R12c1 = d.R12 * c₁

    R21c2 = d.R21 * c₂
    R22s2 = d.R22 * s₂
    R21s2 = d.R21 * s₂
    R22c2 = d.R22 * c₂

    if p.reg == 1
      th0, Δth = d.T1[1], d.T1[2] - d.T1[1]
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)

        dy = fac * muladd(-st, R11s1, ct * R12c1)
        dx = fac * muladd(-st, R11c1, ct * R12s1)

        S = hypot(dx, dy)
        out[1, I] = dy / S
        out[2, I] = dx / S
      end

    elseif p.reg == 5
      th0, Δth = d.T2[1], d.T2[2] - d.T2[1]
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)

        dy = fac * muladd(-st, R21s2, ct * R22c2)
        dx = fac * muladd(-st, R21c2, ct * R22s2)

        S = hypot(dx, dy)
        out[1, I] = dy / S
        out[2, I] = dx / S
      end

    elseif p.reg == 2
      th0, Δth = d.T1[2], d.T1[3] - d.T1[2]
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)

        dy = fac * muladd(-st, R11s1, ct * R12c1)
        dx = fac * muladd(-st, R11c1, ct * R12s1)

        S = hypot(dx, dy)
        out[1, I] = dy / S
        out[2, I] = dx / S
      end

    elseif p.reg == 6
      th0, Δth = d.T2[2], d.T2[3] - d.T2[2]
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)

        dy = fac * muladd(-st, R21s2, ct * R22c2)
        dx = fac * muladd(-st, R21c2, ct * R22s2)

        S = hypot(dx, dy)
        out[1, I] = dy / S
        out[2, I] = dx / S
      end

    elseif p.reg == 3
      th0, Δth = d.T1[3], d.T1[4] - d.T1[3]
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)

        dy = fac * muladd(-st, R11s1, ct * R12c1)
        dx = fac * muladd(-st, R11c1, ct * R12s1)

        S = hypot(dx, dy)
        out[1, I] = dy / S
        out[2, I] = dx / S
      end

    elseif p.reg == 7
      th0, Δth = d.T2[3], d.T2[4] - d.T2[3]
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)

        dy = fac * muladd(-st, R21s2, ct * R22c2)
        dx = fac * muladd(-st, R21c2, ct * R22s2)

        S = hypot(dx, dy)
        out[1, I] = dy / S
        out[2, I] = dx / S
      end

    elseif p.reg == 4
      th0, Δth = d.T1[4], d.T1[1] - d.T1[4] + 2π
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)

        dy = fac * muladd(-st, R11s1, ct * R12c1)
        dx = fac * muladd(-st, R11c1, ct * R12s1)

        S = hypot(dx, dy)
        out[1, I] = dy / S
        out[2, I] = dx / S
      end

    elseif p.reg == 8
      th0, Δth = d.T2[4], d.T2[1] - d.T2[4] + 2π
      fac = Δth * αv
      for I in eachindex(t)
        τ = muladd(αv, t[I], βv)
        th = muladd(Δth, τ, th0)
        st, ct = sincos(th)

        dy = fac * muladd(-st, R21s2, ct * R22c2)
        dx = fac * muladd(-st, R21c2, ct * R22s2)

        S = hypot(dx, dy)
        out[1, I] = dy / S
        out[2, I] = dx / S
      end

    else
      throw(ArgumentError("nu! is defined only for regions 1–8; got reg=$(p.reg)"))
    end
  end

  return nothing
end

#The following function is specific to domains with nₕ >= 1
function bdno(d::annulus, k::Int)::Int
    r = d.pths[k].reg
    return (r <= 4 ? 1 : 2)
end

"""
    DLP(d::annulus, t::Float64, l::Int, tau::StridedArray{Float64}, k::Int)

Double-layer kernel on the boundary:
K(t, τ) = ((γ_k(τ) - γ_l(t)) ⋅ γ'_k(τ)) / ‖γ_k(τ) - γ_l(t)‖² for k ≠ l,
and for k == l the limiting value is taken for patch k.
The array method returns an array with the same size as `tau`.
"""
function DLP!(out::StridedArray{Float64}, d::annulus, t::Float64, l::Int,
  tau::StridedArray{Float64}, k::Int, x::Vector{Float64},
  G::Vector{Float64}, GP::Vector{Float64})

  if k != l
    gam!(x, d, t, l)

    @inbounds for i in eachindex(tau)
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

    αv = 0.5 * ht
    βv = 0.5 * (p.tk0 + p.tk1)

    reg = p.reg

    @inbounds begin
      if reg == 1
        ak = d.T1[2] - d.T1[1]
        bk = d.T1[1]
      elseif reg == 2
        ak = d.T1[3] - d.T1[2]
        bk = d.T1[2]
      elseif reg == 3
        ak = d.T1[4] - d.T1[3]
        bk = d.T1[3]
      elseif reg == 4
        ak = d.T1[1] - d.T1[4] + 2.0 * pi
        bk = d.T1[4]
      elseif reg == 5
        ak = d.T2[2] - d.T2[1]
        bk = d.T2[1]
      elseif reg == 6
        ak = d.T2[3] - d.T2[2]
        bk = d.T2[2]
      elseif reg == 7
        ak = d.T2[4] - d.T2[3]
        bk = d.T2[3]
      elseif reg == 8
        ak = d.T2[1] - d.T2[4] + 2.0 * pi
        bk = d.T2[4]
      else
        throw(ArgumentError("DLP diagonal formula expects reg ∈ 1:8; got reg=$reg"))
      end
    end


    # With xi(t) = αv*t + βv:
    # tt = ak*(xi(t) + xi(tau))/2 + bk
    #    = ak*( (αv*t + βv) + (αv*τ + βv) )/2 + bk
    #    = ak*( αv*(t+τ) + 2βv )/2 + bk
    #    = (ak*αv/2)*(t+τ) + ak*βv + bk = c1*(t+τ) + c0
    c1 = 0.5 * ak * αv
    c0 = ak * βv + bk

    if reg <= 4
      pref = 0.25 * ht * ak * d.R12 * d.R11

      for i in eachindex(out, tau)
        τ = tau[i]

        st, ct = sincos(muladd(c1, t + τ, c0))

        den = (d.R11 * st)^2 + (d.R12 * ct)^2

        out[i] = pref / den
      end

    else
      pref = -0.25 * ht * ak * d.R22 * d.R21

      for i in eachindex(out, tau)
        τ = tau[i]

        st, ct = sincos(muladd(c1, t + τ, c0))

        den = (d.R21 * st)^2 + (d.R22 * ct)^2

        out[i] = pref / den
      end
    end
  end

  return nothing
end

#-----------------------

function Dmap!(out::StridedArray{Float64}, d::annulus,
  u::StridedArray{Float64}, v::StridedArray{Float64}, k::Int)
  @inbounds begin
    # p denotes information of kth patch
    p = d.pths[k]

    hc = p.ck1 - p.ck0
    ht = p.tk1 - p.tk0
    αu = 0.5 * hc
    αv = 0.5 * ht
    αuv= αu * αv
    βu = p.ck0 + αu
    βv = p.tk0 + αv

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11s1, R11c1 = d.R11 * s₁, d.R11 * c₁
    R12s1, R12c1 = d.R12 * s₁, d.R12 * c₁

    R21s2, R21c2 = d.R21 * s₂, d.R21 * c₂
    R22s2, R22c2 = d.R22 * s₂, d.R22 * c₂

    if p.reg == 1

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      @inbounds for I in eachindex(out, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th1 = muladd(Δth1, xi2, th₀1)
        th5 = muladd(Δth5, xi2, th₀5)

        st1, ct1 = sincos(th1)
        st5, ct5 = sincos(th5)

        y1x = muladd(-st1, R12s1, muladd(ct1,  R11c1, d.A))
        y1y = muladd(st1, R12c1, muladd(ct1, R11s1, d.B))

        dy1x= Δth1 * (-R11c1 * st1 - R12s1 * ct1 )
        dy1y= Δth1 * (-R11s1 * st1 + R12c1 * ct1 )

        y5x = muladd(-st5, R22s2, muladd(ct5,  R21c2, d.A))
        y5y = muladd(st5, R22c2, muladd(ct5, R21s2, d.B))

        dy5x= Δth5 * (-R21c2 * st5 - R22s2 * ct5 )
        dy5y= Δth5 * (-R21s2 * st5 + R22c2 * ct5 )

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        #Zx[I] = muladd(tp, y1x - y5x, y5x)
        #Zy[I] = muladd(tp, y1y - y5y, y5y)

        dux = (y1x - y5x) / 2
        duy = (y1y - y5y) / 2
        dvx = tp * dy1x + tm * dy5x
        dvy = tp * dy1y + tm * dy5y

        out[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    elseif p.reg == 5

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      @inbounds for I in eachindex(out, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th1 = muladd(Δth1, xi2, th₀1)
        th5 = muladd(Δth5, xi2, th₀5)

        st1, ct1 = sincos(th1)
        st5, ct5 = sincos(th5)

        y1x = muladd(-st1, R12s1, muladd(ct1,  R11c1, d.A))
        y1y = muladd(st1, R12c1, muladd(ct1, R11s1, d.B))

        dy1x= Δth1 * (-R11c1 * st1 - R12s1 * ct1 )
        dy1y= Δth1 * (-R11s1 * st1 + R12c1 * ct1 )

        y5x = muladd(-st5, R22s2, muladd(ct5,  R21c2, d.A))
        y5y = muladd(st5, R22c2, muladd(ct5, R21s2, d.B))

        dy5x= Δth5 * (-R21c2 * st5 - R22s2 * ct5 )
        dy5y= Δth5 * (-R21s2 * st5 + R22c2 * ct5 )

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        #Zx[I] = muladd(tp, y5x - y1x, y1x)
        #Zy[I] = muladd(tp, y5y - y1y, y1y)

        dux = (y5x - y1x) / 2
        duy = (y5y - y1y) / 2
        dvx = tp * dy5x + tm * dy1x
        dvy = tp * dy5y + tm * dy1y

        out[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    elseif p.reg == 2

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      @inbounds for I in eachindex(out, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th2 = muladd(Δth2, xi2, th₀2)
        th6 = muladd(Δth6, xi2, th₀6)

        st2, ct2 = sincos(th2)
        st6, ct6 = sincos(th6)

        y2x = muladd(-st2, R12s1, muladd(ct2,  R11c1, d.A))
        y2y = muladd(st2, R12c1, muladd(ct2, R11s1, d.B))

        dy2x= Δth2 * (-R11c1 * st2 - R12s1 * ct2 )
        dy2y= Δth2 * (-R11s1 * st2 + R12c1 * ct2 )

        y6x = muladd(-st6, R22s2, muladd(ct6,  R21c2, d.A))
        y6y = muladd(st6, R22c2, muladd(ct6, R21s2, d.B))

        dy6x= Δth6 * (-R21c2 * st6 - R22s2 * ct6 )
        dy6y= Δth6 * (-R21s2 * st6 + R22c2 * ct6 )

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        #Zx[I] = muladd(tp, y2x - y6x, y6x)
        #Zy[I] = muladd(tp, y2y - y6y, y6y)

        dux = (y2x - y6x) / 2
        duy = (y2y - y6y) / 2
        dvx = tp * dy2x + tm * dy6x
        dvy = tp * dy2y + tm * dy6y

        out[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    elseif p.reg == 6

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      @inbounds for I in eachindex(out, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th2 = muladd(Δth2, xi2, th₀2)
        th6 = muladd(Δth6, xi2, th₀6)

        st2, ct2 = sincos(th2)
        st6, ct6 = sincos(th6)

        y2x = muladd(-st2, R12s1, muladd(ct2,  R11c1, d.A))
        y2y = muladd(st2, R12c1, muladd(ct2, R11s1, d.B))

        dy2x= Δth2 * (-R11c1 * st2 - R12s1 * ct2 )
        dy2y= Δth2 * (-R11s1 * st2 + R12c1 * ct2 )

        y6x = muladd(-st6, R22s2, muladd(ct6,  R21c2, d.A))
        y6y = muladd(st6, R22c2, muladd(ct6, R21s2, d.B))

        dy6x= Δth6 * (-R21c2 * st6 - R22s2 * ct6 )
        dy6y= Δth6 * (-R21s2 * st6 + R22c2 * ct6 )

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        #Zx[I] = muladd(tp, y6x - y2x, y2x)
        #Zy[I] = muladd(tp, y6y - y2y, y2y)

        dux = (y6x - y2x) / 2
        duy = (y6y - y2y) / 2
        dvx = tp * dy6x + tm * dy2x
        dvy = tp * dy6y + tm * dy2y

        out[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    elseif p.reg == 3

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      @inbounds for I in eachindex(out, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th3 = muladd(Δth3, xi2, th₀3)
        th7 = muladd(Δth7, xi2, th₀7)

        st3, ct3 = sincos(th3)
        st7, ct7 = sincos(th7)

        y3x = muladd(-st3, R12s1, muladd(ct3,  R11c1, d.A))
        y3y = muladd(st3, R12c1, muladd(ct3, R11s1, d.B))

        dy3x= Δth3 * (-R11c1 * st3 - R12s1 * ct3 )
        dy3y= Δth3 * (-R11s1 * st3 + R12c1 * ct3 )

        y7x = muladd(-st7, R22s2, muladd(ct7,  R21c2, d.A))
        y7y = muladd(st7, R22c2, muladd(ct7, R21s2, d.B))

        dy7x= Δth7 * (-R21c2 * st7 - R22s2 * ct7 )
        dy7y= Δth7 * (-R21s2 * st7 + R22c2 * ct7 )

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        #Zx[I] = muladd(tp, y3x - y7x, y7x)
        #Zy[I] = muladd(tp, y3y - y7y, y7y)

        dux = (y3x - y7x) / 2
        duy = (y3y - y7y) / 2
        dvx = tp * dy3x + tm * dy7x
        dvy = tp * dy3y + tm * dy7y

        out[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    elseif p.reg == 7

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      @inbounds for I in eachindex(out, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th3 = muladd(Δth3, xi2, th₀3)
        th7 = muladd(Δth7, xi2, th₀7)

        st3, ct3 = sincos(th3)
        st7, ct7 = sincos(th7)

        y3x = muladd(-st3, R12s1, muladd(ct3,  R11c1, d.A))
        y3y = muladd(st3, R12c1, muladd(ct3, R11s1, d.B))

        dy3x= Δth3 * (-R11c1 * st3 - R12s1 * ct3 )
        dy3y= Δth3 * (-R11s1 * st3 + R12c1 * ct3 )

        y7x = muladd(-st7, R22s2, muladd(ct7,  R21c2, d.A))
        y7y = muladd(st7, R22c2, muladd(ct7, R21s2, d.B))

        dy7x= Δth7 * (-R21c2 * st7 - R22s2 * ct7 )
        dy7y= Δth7 * (-R21s2 * st7 + R22c2 * ct7 )

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        #Zx[I] = muladd(tp, y7x - y3x, y3x)
        #Zy[I] = muladd(tp, y7y - y3y, y3y)

        dux = (y7x - y3x) / 2
        duy = (y7y - y3y) / 2
        dvx = tp * dy7x + tm * dy3x
        dvy = tp * dy7y + tm * dy3y

        out[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    elseif p.reg == 4

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      @inbounds for I in eachindex(out, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th4 = muladd(Δth4, xi2, th₀4)
        th8 = muladd(Δth8, xi2, th₀8)

        st4, ct4 = sincos(th4)
        st8, ct8 = sincos(th8)

        y4x = muladd(-st4, R12s1, muladd(ct4,  R11c1, d.A))
        y4y = muladd(st4, R12c1, muladd(ct4, R11s1, d.B))

        dy4x= Δth4 * (-R11c1 * st4 - R12s1 * ct4 )
        dy4y= Δth4 * (-R11s1 * st4 + R12c1 * ct4 )

        y8x = muladd(-st8, R22s2, muladd(ct8,  R21c2, d.A))
        y8y = muladd(st8, R22c2, muladd(ct8, R21s2, d.B))

        dy8x= Δth8 * (-R21c2 * st8 - R22s2 * ct8 )
        dy8y= Δth8 * (-R21s2 * st8 + R22c2 * ct8 )

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        #Zx[I] = muladd(tp, y4x - y8x, y8x)
        #Zy[I] = muladd(tp, y4y - y8y, y8y)

        dux = (y4x - y8x) / 2
        duy = (y4y - y8y) / 2
        dvx = tp * dy4x + tm * dy8x
        dvy = tp * dy4y + tm * dy8y

        out[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    elseif p.reg == 8

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      @inbounds for I in eachindex(out, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th4 = muladd(Δth4, xi2, th₀4)
        th8 = muladd(Δth8, xi2, th₀8)

        st4, ct4 = sincos(th4)
        st8, ct8 = sincos(th8)

        y4x = muladd(-st4, R12s1, muladd(ct4,  R11c1, d.A))
        y4y = muladd(st4, R12c1, muladd(ct4, R11s1, d.B))

        dy4x= Δth4 * (-R11c1 * st4 - R12s1 * ct4 )
        dy4y= Δth4 * (-R11s1 * st4 + R12c1 * ct4 )

        y8x = muladd(-st8, R22s2, muladd(ct8,  R21c2, d.A))
        y8y = muladd(st8, R22c2, muladd(ct8, R21s2, d.B))

        dy8x= Δth8 * (-R21c2 * st8 - R22s2 * ct8 )
        dy8y= Δth8 * (-R21s2 * st8 + R22c2 * ct8 )

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        #Zx[I] = muladd(tp, y8x - y4x, y4x)
        #Zy[I] = muladd(tp, y8y - y4y, y4y)

        dux = (y8x - y4x) / 2
        duy = (y8y - y4y) / 2
        dvx = tp * dy8x + tm * dy4x
        dvy = tp * dy8y + tm * dy4y

        out[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    else
      throw(ArgumentError("Dmap! is defined only for regions 1–8; got reg=$(p.reg)"))

    end
  end

  return nothing
end

"""
A combination of mapxy! and Dmap! function (No allocations!)
The purpose of this function is to reduce computations  
related to cosine and sine's.  
"""
function mapxy_Dmap!(Zx::StridedArray{Float64}, Zy::StridedArray{Float64},
  DJ::StridedArray{Float64}, d::annulus, u::StridedArray{Float64}, 
  v::StridedArray{Float64}, k::Int)
  @inbounds begin
    # p denotes information of kth patch
    p = d.pths[k]

    hc = p.ck1 - p.ck0
    ht = p.tk1 - p.tk0
    αu = 0.5 * hc
    αv = 0.5 * ht
    αuv = αu * αv
    βu = p.ck0 + αu
    βv = p.tk0 + αv

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11s1, R11c1 = d.R11 * s₁, d.R11 * c₁
    R12s1, R12c1 = d.R12 * s₁, d.R12 * c₁

    R21s2, R21c2 = d.R21 * s₂, d.R21 * c₂
    R22s2, R22c2 = d.R22 * s₂, d.R22 * c₂

    if p.reg == 1

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th1 = muladd(Δth1, xi2, th₀1)
        th5 = muladd(Δth5, xi2, th₀5)

        st1, ct1 = sincos(th1)
        st5, ct5 = sincos(th5)

        y1x = muladd(-st1, R12s1, muladd(ct1, R11c1, d.A))
        y1y = muladd(st1, R12c1, muladd(ct1, R11s1, d.B))

        dy1x = Δth1 * (-R11c1 * st1 - R12s1 * ct1)
        dy1y = Δth1 * (-R11s1 * st1 + R12c1 * ct1)

        y5x = muladd(-st5, R22s2, muladd(ct5, R21c2, d.A))
        y5y = muladd(st5, R22c2, muladd(ct5, R21s2, d.B))

        dy5x = Δth5 * (-R21c2 * st5 - R22s2 * ct5)
        dy5y = Δth5 * (-R21s2 * st5 + R22c2 * ct5)

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        y1_y5x = y1x - y5x
        y1_y5y = y1y - y5y

        Zx[I] = muladd(tp, y1_y5x, y5x)
        Zy[I] = muladd(tp, y1_y5y, y5y)

        dux = y1_y5x / 2
        duy = y1_y5y / 2
        dvx = tp * dy1x + tm * dy5x
        dvy = tp * dy1y + tm * dy5y

        DJ[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    elseif p.reg == 5

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th1 = muladd(Δth1, xi2, th₀1)
        th5 = muladd(Δth5, xi2, th₀5)

        st1, ct1 = sincos(th1)
        st5, ct5 = sincos(th5)

        y1x = muladd(-st1, R12s1, muladd(ct1, R11c1, d.A))
        y1y = muladd(st1, R12c1, muladd(ct1, R11s1, d.B))

        dy1x = Δth1 * (-R11c1 * st1 - R12s1 * ct1)
        dy1y = Δth1 * (-R11s1 * st1 + R12c1 * ct1)

        y5x = muladd(-st5, R22s2, muladd(ct5, R21c2, d.A))
        y5y = muladd(st5, R22c2, muladd(ct5, R21s2, d.B))

        dy5x = Δth5 * (-R21c2 * st5 - R22s2 * ct5)
        dy5y = Δth5 * (-R21s2 * st5 + R22c2 * ct5)

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        y5_y1x = y5x - y1x
        y5_y1y = y5y - y1y

        Zx[I] = muladd(tp, y5_y1x, y1x)
        Zy[I] = muladd(tp, y5_y1y, y1y)

        dux = y5_y1x / 2
        duy = y5_y1y / 2
        dvx = tp * dy5x + tm * dy1x
        dvy = tp * dy5y + tm * dy1y

        DJ[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    elseif p.reg == 2

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th2 = muladd(Δth2, xi2, th₀2)
        th6 = muladd(Δth6, xi2, th₀6)

        st2, ct2 = sincos(th2)
        st6, ct6 = sincos(th6)

        y2x = muladd(-st2, R12s1, muladd(ct2, R11c1, d.A))
        y2y = muladd(st2, R12c1, muladd(ct2, R11s1, d.B))

        dy2x = Δth2 * (-R11c1 * st2 - R12s1 * ct2)
        dy2y = Δth2 * (-R11s1 * st2 + R12c1 * ct2)

        y6x = muladd(-st6, R22s2, muladd(ct6, R21c2, d.A))
        y6y = muladd(st6, R22c2, muladd(ct6, R21s2, d.B))

        dy6x = Δth6 * (-R21c2 * st6 - R22s2 * ct6)
        dy6y = Δth6 * (-R21s2 * st6 + R22c2 * ct6)

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        y2_y6x = y2x - y6x
        y2_y6y = y2y - y6y

        Zx[I] = muladd(tp, y2_y6x, y6x)
        Zy[I] = muladd(tp, y2_y6y, y6y)

        dux = y2_y6x / 2
        duy = y2_y6y / 2
        dvx = tp * dy2x + tm * dy6x
        dvy = tp * dy2y + tm * dy6y

        DJ[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    elseif p.reg == 6

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th2 = muladd(Δth2, xi2, th₀2)
        th6 = muladd(Δth6, xi2, th₀6)

        st2, ct2 = sincos(th2)
        st6, ct6 = sincos(th6)

        y2x = muladd(-st2, R12s1, muladd(ct2, R11c1, d.A))
        y2y = muladd(st2, R12c1, muladd(ct2, R11s1, d.B))

        dy2x = Δth2 * (-R11c1 * st2 - R12s1 * ct2)
        dy2y = Δth2 * (-R11s1 * st2 + R12c1 * ct2)

        y6x = muladd(-st6, R22s2, muladd(ct6, R21c2, d.A))
        y6y = muladd(st6, R22c2, muladd(ct6, R21s2, d.B))

        dy6x = Δth6 * (-R21c2 * st6 - R22s2 * ct6)
        dy6y = Δth6 * (-R21s2 * st6 + R22c2 * ct6)

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        y6_y2x = y6x - y2x
        y6_y2y = y6y - y2y

        Zx[I] = muladd(tp, y6_y2x, y2x)
        Zy[I] = muladd(tp, y6_y2y, y2y)

        dux = y6_y2x / 2
        duy = y6_y2y / 2
        dvx = tp * dy6x + tm * dy2x
        dvy = tp * dy6y + tm * dy2y

        DJ[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    elseif p.reg == 3

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th3 = muladd(Δth3, xi2, th₀3)
        th7 = muladd(Δth7, xi2, th₀7)

        st3, ct3 = sincos(th3)
        st7, ct7 = sincos(th7)

        y3x = muladd(-st3, R12s1, muladd(ct3, R11c1, d.A))
        y3y = muladd(st3, R12c1, muladd(ct3, R11s1, d.B))

        dy3x = Δth3 * (-R11c1 * st3 - R12s1 * ct3)
        dy3y = Δth3 * (-R11s1 * st3 + R12c1 * ct3)

        y7x = muladd(-st7, R22s2, muladd(ct7, R21c2, d.A))
        y7y = muladd(st7, R22c2, muladd(ct7, R21s2, d.B))

        dy7x = Δth7 * (-R21c2 * st7 - R22s2 * ct7)
        dy7y = Δth7 * (-R21s2 * st7 + R22c2 * ct7)

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        y3_y7x = y3x - y7x
        y3_y7y = y3y - y7y

        Zx[I] = muladd(tp, y3_y7x, y7x)
        Zy[I] = muladd(tp, y3_y7y, y7y)

        dux = y3_y7x / 2
        duy = y3_y7y / 2
        dvx = tp * dy3x + tm * dy7x
        dvy = tp * dy3y + tm * dy7y

        DJ[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    elseif p.reg == 7

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th3 = muladd(Δth3, xi2, th₀3)
        th7 = muladd(Δth7, xi2, th₀7)

        st3, ct3 = sincos(th3)
        st7, ct7 = sincos(th7)

        y3x = muladd(-st3, R12s1, muladd(ct3, R11c1, d.A))
        y3y = muladd(st3, R12c1, muladd(ct3, R11s1, d.B))

        dy3x = Δth3 * (-R11c1 * st3 - R12s1 * ct3)
        dy3y = Δth3 * (-R11s1 * st3 + R12c1 * ct3)

        y7x = muladd(-st7, R22s2, muladd(ct7, R21c2, d.A))
        y7y = muladd(st7, R22c2, muladd(ct7, R21s2, d.B))

        dy7x = Δth7 * (-R21c2 * st7 - R22s2 * ct7)
        dy7y = Δth7 * (-R21s2 * st7 + R22c2 * ct7)

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        y7_y3x = y7x - y3x
        y7_y3y = y7y - y3y

        Zx[I] = muladd(tp, y7_y3x, y3x)
        Zy[I] = muladd(tp, y7_y3y, y3y)

        dux = y7_y3x / 2
        duy = y7_y3y / 2
        dvx = tp * dy7x + tm * dy3x
        dvy = tp * dy7y + tm * dy3y

        DJ[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    elseif p.reg == 4

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th4 = muladd(Δth4, xi2, th₀4)
        th8 = muladd(Δth8, xi2, th₀8)

        st4, ct4 = sincos(th4)
        st8, ct8 = sincos(th8)

        y4x = muladd(-st4, R12s1, muladd(ct4, R11c1, d.A))
        y4y = muladd(st4, R12c1, muladd(ct4, R11s1, d.B))

        dy4x = Δth4 * (-R11c1 * st4 - R12s1 * ct4)
        dy4y = Δth4 * (-R11s1 * st4 + R12c1 * ct4)

        y8x = muladd(-st8, R22s2, muladd(ct8, R21c2, d.A))
        y8y = muladd(st8, R22c2, muladd(ct8, R21s2, d.B))

        dy8x = Δth8 * (-R21c2 * st8 - R22s2 * ct8)
        dy8y = Δth8 * (-R21s2 * st8 + R22c2 * ct8)

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        y4_y8x = y4x - y8x
        y4_y8y = y4y - y8y

        Zx[I] = muladd(tp, y4_y8x, y8x)
        Zy[I] = muladd(tp, y4_y8y, y8y)

        dux = y4_y8x / 2
        duy = y4_y8y / 2
        dvx = tp * dy4x + tm * dy8x
        dvy = tp * dy4y + tm * dy8y

        DJ[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    elseif p.reg == 8

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      @inbounds for I in eachindex(Zx, Zy, DJ, u, v)
        xi1 = muladd(αu, u[I], βu)
        xi2 = muladd(αv, v[I], βv)

        th4 = muladd(Δth4, xi2, th₀4)
        th8 = muladd(Δth8, xi2, th₀8)

        st4, ct4 = sincos(th4)
        st8, ct8 = sincos(th8)

        y4x = muladd(-st4, R12s1, muladd(ct4, R11c1, d.A))
        y4y = muladd(st4, R12c1, muladd(ct4, R11s1, d.B))

        dy4x = Δth4 * (-R11c1 * st4 - R12s1 * ct4)
        dy4y = Δth4 * (-R11s1 * st4 + R12c1 * ct4)

        y8x = muladd(-st8, R22s2, muladd(ct8, R21c2, d.A))
        y8y = muladd(st8, R22c2, muladd(ct8, R21s2, d.B))

        dy8x = Δth8 * (-R21c2 * st8 - R22s2 * ct8)
        dy8y = Δth8 * (-R21s2 * st8 + R22c2 * ct8)

        tp = muladd(0.5, xi1, 0.5)
        tm = muladd(0.5, -xi1, 0.5)

        y8_y4x = y8x - y4x
        y8_y4y = y8y - y4y

        Zx[I] = muladd(tp, y8_y4x, y4x)
        Zy[I] = muladd(tp, y8_y4y, y4y)

        dux = y8_y4x / 2
        duy = y8_y4y / 2
        dvx = tp * dy8x + tm * dy4x
        dvy = tp * dy8y + tm * dy4y

        DJ[I] = αuv * abs(dux * dvy - dvx * duy)
      end

    else
      throw(ArgumentError("Dmap! is defined only for regions 1–8; got reg=$(p.reg)"))

    end
  end

  return nothing
end

function chk_map(d::annulus)
  # --- Test function f(x,y) ---
  #f!(F, x, y) = fill!(F, 1.0)
  #Iex =  π * (d.R11 * d.R12 -  d.R21 * d.R22)

  n = 32

  f!(F, x, y) = @. F = x^2 + y^2
  Iex = pi * (d.R11 * d.R12 * (d.R11^2 + d.R12^2) / 4
              -d.R21 * d.R22 * (d.R21^2 + d.R22^2) / 4
              +(d.R11 * d.R12 - d.R21 * d.R22) * (d.A^2 + d.B^2))

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
    jinvmap(d::annulus, u::Float64, v::Float64, r::Int) -> tuple

Inverse Jacobian of the mapping at a single reference point `(u,v) ∈ [-1,1]^2`
for the given patch. Returns a 4 Tuple which are components of matrix `J⁻¹`.
Only the region of the patch is needed as inverse are only computed over them.
ht = hc = 1
"""
# Return (J11, J12, J21, J22) of the inverse Jacobian
function jinvmap(d::annulus, u::Float64, v::Float64, r::Int)

  @inbounds begin
    xi1 = muladd(0.5, u, 0.5)
    xi2 = muladd(0.5, v, 0.5)

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11s1, R11c1 = d.R11 * s₁, d.R11 * c₁
    R12s1, R12c1 = d.R12 * s₁, d.R12 * c₁

    R21s2, R21c2 = d.R21 * s₂, d.R21 * c₂
    R22s2, R22c2 = d.R22 * s₂, d.R22 * c₂

    if r == 1

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      th1 = muladd(Δth1, xi2, th₀1)
      th5 = muladd(Δth5, xi2, th₀5)

      st1, ct1 = sincos(th1)
      st5, ct5 = sincos(th5)

      y1x = muladd(-st1, R12s1, muladd(ct1, R11c1, d.A))
      y1y = muladd(st1, R12c1, muladd(ct1, R11s1, d.B))

      dy1x = Δth1 * (-R11c1 * st1 - R12s1 * ct1)
      dy1y = Δth1 * (-R11s1 * st1 + R12c1 * ct1)

      y5x = muladd(-st5, R22s2, muladd(ct5, R21c2, d.A))
      y5y = muladd(st5, R22c2, muladd(ct5, R21s2, d.B))

      dy5x = Δth5 * (-R21c2 * st5 - R22s2 * ct5)
      dy5y = Δth5 * (-R21s2 * st5 + R22c2 * ct5)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dux = (y1x - y5x) / 2
      duy = (y1y - y5y) / 2
      dvx = tp * dy1x + tm * dy5x
      dvy = tp * dy1y + tm * dy5y

    elseif r == 5

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      th1 = muladd(Δth1, xi2, th₀1)
      th5 = muladd(Δth5, xi2, th₀5)

      st1, ct1 = sincos(th1)
      st5, ct5 = sincos(th5)

      y1x = muladd(-st1, R12s1, muladd(ct1, R11c1, d.A))
      y1y = muladd(st1, R12c1, muladd(ct1, R11s1, d.B))

      dy1x = Δth1 * (-R11c1 * st1 - R12s1 * ct1)
      dy1y = Δth1 * (-R11s1 * st1 + R12c1 * ct1)

      y5x = muladd(-st5, R22s2, muladd(ct5, R21c2, d.A))
      y5y = muladd(st5, R22c2, muladd(ct5, R21s2, d.B))

      dy5x = Δth5 * (-R21c2 * st5 - R22s2 * ct5)
      dy5y = Δth5 * (-R21s2 * st5 + R22c2 * ct5)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dux = (y5x - y1x) / 2
      duy = (y5y - y1y) / 2
      dvx = tp * dy5x + tm * dy1x
      dvy = tp * dy5y + tm * dy1y

    elseif r == 2

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      th2 = muladd(Δth2, xi2, th₀2)
      th6 = muladd(Δth6, xi2, th₀6)

      st2, ct2 = sincos(th2)
      st6, ct6 = sincos(th6)

      y2x = muladd(-st2, R12s1, muladd(ct2, R11c1, d.A))
      y2y = muladd(st2, R12c1, muladd(ct2, R11s1, d.B))

      dy2x = Δth2 * (-R11c1 * st2 - R12s1 * ct2)
      dy2y = Δth2 * (-R11s1 * st2 + R12c1 * ct2)

      y6x = muladd(-st6, R22s2, muladd(ct6, R21c2, d.A))
      y6y = muladd(st6, R22c2, muladd(ct6, R21s2, d.B))

      dy6x = Δth6 * (-R21c2 * st6 - R22s2 * ct6)
      dy6y = Δth6 * (-R21s2 * st6 + R22c2 * ct6)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dux = (y2x - y6x) / 2
      duy = (y2y - y6y) / 2
      dvx = tp * dy2x + tm * dy6x
      dvy = tp * dy2y + tm * dy6y

    elseif r == 6

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      th2 = muladd(Δth2, xi2, th₀2)
      th6 = muladd(Δth6, xi2, th₀6)

      st2, ct2 = sincos(th2)
      st6, ct6 = sincos(th6)

      y2x = muladd(-st2, R12s1, muladd(ct2, R11c1, d.A))
      y2y = muladd(st2, R12c1, muladd(ct2, R11s1, d.B))

      dy2x = Δth2 * (-R11c1 * st2 - R12s1 * ct2)
      dy2y = Δth2 * (-R11s1 * st2 + R12c1 * ct2)

      y6x = muladd(-st6, R22s2, muladd(ct6, R21c2, d.A))
      y6y = muladd(st6, R22c2, muladd(ct6, R21s2, d.B))

      dy6x = Δth6 * (-R21c2 * st6 - R22s2 * ct6)
      dy6y = Δth6 * (-R21s2 * st6 + R22c2 * ct6)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dux = (y6x - y2x) / 2
      duy = (y6y - y2y) / 2
      dvx = tp * dy6x + tm * dy2x
      dvy = tp * dy6y + tm * dy2y

    elseif r == 3

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      th3 = muladd(Δth3, xi2, th₀3)
      th7 = muladd(Δth7, xi2, th₀7)

      st3, ct3 = sincos(th3)
      st7, ct7 = sincos(th7)

      y3x = muladd(-st3, R12s1, muladd(ct3, R11c1, d.A))
      y3y = muladd(st3, R12c1, muladd(ct3, R11s1, d.B))

      dy3x = Δth3 * (-R11c1 * st3 - R12s1 * ct3)
      dy3y = Δth3 * (-R11s1 * st3 + R12c1 * ct3)

      y7x = muladd(-st7, R22s2, muladd(ct7, R21c2, d.A))
      y7y = muladd(st7, R22c2, muladd(ct7, R21s2, d.B))

      dy7x = Δth7 * (-R21c2 * st7 - R22s2 * ct7)
      dy7y = Δth7 * (-R21s2 * st7 + R22c2 * ct7)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dux = (y3x - y7x) / 2
      duy = (y3y - y7y) / 2
      dvx = tp * dy3x + tm * dy7x
      dvy = tp * dy3y + tm * dy7y

    elseif r == 7

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      th3 = muladd(Δth3, xi2, th₀3)
      th7 = muladd(Δth7, xi2, th₀7)

      st3, ct3 = sincos(th3)
      st7, ct7 = sincos(th7)

      y3x = muladd(-st3, R12s1, muladd(ct3, R11c1, d.A))
      y3y = muladd(st3, R12c1, muladd(ct3, R11s1, d.B))

      dy3x = Δth3 * (-R11c1 * st3 - R12s1 * ct3)
      dy3y = Δth3 * (-R11s1 * st3 + R12c1 * ct3)

      y7x = muladd(-st7, R22s2, muladd(ct7, R21c2, d.A))
      y7y = muladd(st7, R22c2, muladd(ct7, R21s2, d.B))

      dy7x = Δth7 * (-R21c2 * st7 - R22s2 * ct7)
      dy7y = Δth7 * (-R21s2 * st7 + R22c2 * ct7)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dux = (y7x - y3x) / 2
      duy = (y7y - y3y) / 2
      dvx = tp * dy7x + tm * dy3x
      dvy = tp * dy7y + tm * dy3y

    elseif r == 4

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      th4 = muladd(Δth4, xi2, th₀4)
      th8 = muladd(Δth8, xi2, th₀8)

      st4, ct4 = sincos(th4)
      st8, ct8 = sincos(th8)

      y4x = muladd(-st4, R12s1, muladd(ct4, R11c1, d.A))
      y4y = muladd(st4, R12c1, muladd(ct4, R11s1, d.B))

      dy4x = Δth4 * (-R11c1 * st4 - R12s1 * ct4)
      dy4y = Δth4 * (-R11s1 * st4 + R12c1 * ct4)

      y8x = muladd(-st8, R22s2, muladd(ct8, R21c2, d.A))
      y8y = muladd(st8, R22c2, muladd(ct8, R21s2, d.B))

      dy8x = Δth8 * (-R21c2 * st8 - R22s2 * ct8)
      dy8y = Δth8 * (-R21s2 * st8 + R22c2 * ct8)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dux = (y4x - y8x) / 2
      duy = (y4y - y8y) / 2
      dvx = tp * dy4x + tm * dy8x
      dvy = tp * dy4y + tm * dy8y

    elseif r == 8

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      th4 = muladd(Δth4, xi2, th₀4)
      th8 = muladd(Δth8, xi2, th₀8)

      st4, ct4 = sincos(th4)
      st8, ct8 = sincos(th8)

      y4x = muladd(-st4, R12s1, muladd(ct4, R11c1, d.A))
      y4y = muladd(st4, R12c1, muladd(ct4, R11s1, d.B))

      dy4x = Δth4 * (-R11c1 * st4 - R12s1 * ct4)
      dy4y = Δth4 * (-R11s1 * st4 + R12c1 * ct4)

      y8x = muladd(-st8, R22s2, muladd(ct8, R21c2, d.A))
      y8y = muladd(st8, R22c2, muladd(ct8, R21s2, d.B))

      dy8x = Δth8 * (-R21c2 * st8 - R22s2 * ct8)
      dy8y = Δth8 * (-R21s2 * st8 + R22c2 * ct8)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dux = (y8x - y4x) / 2
      duy = (y8y - y4y) / 2
      dvx = tp * dy8x + tm * dy4x
      dvy = tp * dy8y + tm * dy4y

    else
      throw(ArgumentError("jinvmap is defined only for regions 1–8; got reg=$r"))

    end
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
  mapinv(d::annulus, u::Float64, v::Float64, k::Int) -> Tuple{Float64,Float64}

Inverse mapping: given a physical point `(u,v)` on patch `k`, return
the reference coordinates `[t; s]` in `[-1,1]^2` such that `mapxy(d, t, s, k) = (u,v)`.
"""
function Xx(s::Float64, d::annulus, r::Int)
  @inbounds begin

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    if r == 1

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      th1 = muladd(Δth1, s, th₀1)
      th5 = muladd(Δth5, s, th₀5)

      st1, ct1 = sincos(th1)
      st5, ct5 = sincos(th5)

      y1x = d.A + d.R11 * ct1 * c₁ - d.R12 * st1 * s₁
      y5x = d.A + d.R21 * ct5 * c₂ - d.R22 * st5 * s₂

      return (y1x + y5x) / 2

    elseif r == 5

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      th1 = muladd(Δth1, s, th₀1)
      th5 = muladd(Δth5, s, th₀5)

      st1, ct1 = sincos(th1)
      st5, ct5 = sincos(th5)

      y1x = d.A + d.R11 * ct1 * c₁ - d.R12 * st1 * s₁
      y5x = d.A + d.R21 * ct5 * c₂ - d.R22 * st5 * s₂

      return (y5x + y1x) / 2

    elseif r == 2

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      th2 = muladd(Δth2, s, th₀2)
      th6 = muladd(Δth6, s, th₀6)

      st2, ct2 = sincos(th2)
      st6, ct6 = sincos(th6)

      y2x = d.A + d.R11 * ct2 * c₁ - d.R12 * st2 * s₁
      y6x = d.A + d.R21 * ct6 * c₂ - d.R22 * st6 * s₂

      return (y2x + y6x) / 2

    elseif r == 6

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      th2 = muladd(Δth2, s, th₀2)
      th6 = muladd(Δth6, s, th₀6)

      st2, ct2 = sincos(th2)
      st6, ct6 = sincos(th6)

      y2x = d.A + d.R11 * ct2 * c₁ - d.R12 * st2 * s₁
      y6x = d.A + d.R21 * ct6 * c₂ - d.R22 * st6 * s₂

      return (y6x + y2x) / 2

    elseif r == 3

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      th3 = muladd(Δth3, s, th₀3)
      th7 = muladd(Δth7, s, th₀7)

      st3, ct3 = sincos(th3)
      st7, ct7 = sincos(th7)

      y3x = d.A + d.R11 * ct3 * c₁ - d.R12 * st3 * s₁
      y7x = d.A + d.R21 * ct7 * c₂ - d.R22 * st7 * s₂

      return (y3x + y7x) / 2

    elseif r == 7

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      th3 = muladd(Δth3, s, th₀3)
      th7 = muladd(Δth7, s, th₀7)

      st3, ct3 = sincos(th3)
      st7, ct7 = sincos(th7)

      y3x = d.A + d.R11 * ct3 * c₁ - d.R12 * st3 * s₁
      y7x = d.A + d.R21 * ct7 * c₂ - d.R22 * st7 * s₂

      return (y7x + y3x) / 2

    elseif r == 4

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      th4 = muladd(Δth4, s, th₀4)
      th8 = muladd(Δth8, s, th₀8)

      st4, ct4 = sincos(th4)
      st8, ct8 = sincos(th8)

      y4x = d.A + d.R11 * ct4 * c₁ - d.R12 * st4 * s₁
      y8x = d.A + d.R21 * ct8 * c₂ - d.R22 * st8 * s₂

      return (y4x + y8x) / 2

    else

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      th4 = muladd(Δth4, s, th₀4)
      th8 = muladd(Δth8, s, th₀8)

      st4, ct4 = sincos(th4)
      st8, ct8 = sincos(th8)

      y4x = d.A + d.R11 * ct4 * c₁ - d.R12 * st4 * s₁
      y8x = d.A + d.R21 * ct8 * c₂ - d.R22 * st8 * s₂

      return (y8x + y4x) / 2

    end
  end
end

function Xy(s::Float64, d::annulus, r::Int)
  @inbounds begin

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    if r == 1

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      th1 = muladd(Δth1, s, th₀1)
      th5 = muladd(Δth5, s, th₀5)

      st1, ct1 = sincos(th1)
      st5, ct5 = sincos(th5)

      y1y = d.B + d.R11 * ct1 * s₁ + d.R12 * st1 * c₁
      y5y = d.B + d.R21 * ct5 * s₂ + d.R22 * st5 * c₂

      return (y1y + y5y) / 2

    elseif r == 5

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      th1 = muladd(Δth1, s, th₀1)
      th5 = muladd(Δth5, s, th₀5)

      st1, ct1 = sincos(th1)
      st5, ct5 = sincos(th5)

      y1y = d.B + d.R11 * ct1 * s₁ + d.R12 * st1 * c₁
      y5y = d.B + d.R21 * ct5 * s₂ + d.R22 * st5 * c₂

      return (y5y + y1y) / 2

    elseif r == 2

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      th2 = muladd(Δth2, s, th₀2)
      th6 = muladd(Δth6, s, th₀6)

      st2, ct2 = sincos(th2)
      st6, ct6 = sincos(th6)

      y2y = d.B + d.R11 * ct2 * s₁ + d.R12 * st2 * c₁
      y6y = d.B + d.R21 * ct6 * s₂ + d.R22 * st6 * c₂

      return (y2y + y6y) / 2

    elseif r == 6

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      th2 = muladd(Δth2, s, th₀2)
      th6 = muladd(Δth6, s, th₀6)

      st2, ct2 = sincos(th2)
      st6, ct6 = sincos(th6)

      y2y = d.B + d.R11 * ct2 * s₁ + d.R12 * st2 * c₁
      y6y = d.B + d.R21 * ct6 * s₂ + d.R22 * st6 * c₂

      return (y6y + y2y) / 2

    elseif r == 3

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      th3 = muladd(Δth3, s, th₀3)
      th7 = muladd(Δth7, s, th₀7)

      st3, ct3 = sincos(th3)
      st7, ct7 = sincos(th7)

      y3y = d.B + d.R11 * ct3 * s₁ + d.R12 * st3 * c₁
      y7y = d.B + d.R21 * ct7 * s₂ + d.R22 * st7 * c₂

      return (y3y + y7y) / 2

    elseif r == 7

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      th3 = muladd(Δth3, s, th₀3)
      th7 = muladd(Δth7, s, th₀7)

      st3, ct3 = sincos(th3)
      st7, ct7 = sincos(th7)

      y3y = d.B + d.R11 * ct3 * s₁ + d.R12 * st3 * c₁
      y7y = d.B + d.R21 * ct7 * s₂ + d.R22 * st7 * c₂

      return (y7y + y3y) / 2

    elseif r == 4

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      th4 = muladd(Δth4, s, th₀4)
      th8 = muladd(Δth8, s, th₀8)

      st4, ct4 = sincos(th4)
      st8, ct8 = sincos(th8)

      y4y = d.B + d.R11 * ct4 * s₁ + d.R12 * st4 * c₁
      y8y = d.B + d.R21 * ct8 * s₂ + d.R22 * st8 * c₂

      return (y4y + y8y) / 2

    else

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      th4 = muladd(Δth4, s, th₀4)
      th8 = muladd(Δth8, s, th₀8)

      st4, ct4 = sincos(th4)
      st8, ct8 = sincos(th8)

      y4y = d.B + d.R11 * ct4 * s₁ + d.R12 * st4 * c₁
      y8y = d.B + d.R21 * ct8 * s₂ + d.R22 * st8 * c₂

      return (y8y + y4y) / 2

    end
  end
end

function Yx(s::Float64, d::annulus, r::Int)
  @inbounds begin

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    if r == 1

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]

      th1 = muladd(Δth1, s, th₀1)

      st1, ct1 = sincos(th1)

      return d.A + d.R11 * ct1 * c₁ - d.R12 * st1 * s₁

    elseif r == 5

      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      th5 = muladd(Δth5, s, th₀5)

      st5, ct5 = sincos(th5)

      return d.A + d.R21 * ct5 * c₂ - d.R22 * st5 * s₂

    elseif r == 2

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]

      th2 = muladd(Δth2, s, th₀2)

      st2, ct2 = sincos(th2)

      return d.A + d.R11 * ct2 * c₁ - d.R12 * st2 * s₁

    elseif r == 6

      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      th6 = muladd(Δth6, s, th₀6)

      st6, ct6 = sincos(th6)

      return d.A + d.R21 * ct6 * c₂ - d.R22 * st6 * s₂

    elseif r == 3

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]

      th3 = muladd(Δth3, s, th₀3)

      st3, ct3 = sincos(th3)

      return d.A + d.R11 * ct3 * c₁ - d.R12 * st3 * s₁

    elseif r == 7

      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      th7 = muladd(Δth7, s, th₀7)

      st7, ct7 = sincos(th7)

      return d.A + d.R21 * ct7 * c₂ - d.R22 * st7 * s₂

    elseif r == 4

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi

      th4 = muladd(Δth4, s, th₀4)

      st4, ct4 = sincos(th4)

      return d.A + d.R11 * ct4 * c₁ - d.R12 * st4 * s₁

    else

      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      th8 = muladd(Δth8, s, th₀8)

      st8, ct8 = sincos(th8)

      return d.A + d.R21 * ct8 * c₂ - d.R22 * st8 * s₂

    end
  end
end

function Yy(s::Float64, d::annulus, r::Int)

  @inbounds begin

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    if r == 1

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]

      th1 = muladd(Δth1, s, th₀1)

      st1, ct1 = sincos(th1)

      return d.B + d.R11 * ct1 * s₁ + d.R12 * st1 * c₁

    elseif r == 5

      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      th5 = muladd(Δth5, s, th₀5)

      st5, ct5 = sincos(th5)

      return d.B + d.R21 * ct5 * s₂ + d.R22 * st5 * c₂

    elseif r == 2

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]

      th2 = muladd(Δth2, s, th₀2)

      st2, ct2 = sincos(th2)

      return d.B + d.R11 * ct2 * s₁ + d.R12 * st2 * c₁

    elseif r == 6

      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      th6 = muladd(Δth6, s, th₀6)

      st6, ct6 = sincos(th6)

      return d.B + d.R21 * ct6 * s₂ + d.R22 * st6 * c₂

    elseif r == 3

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]

      th3 = muladd(Δth3, s, th₀3)

      st3, ct3 = sincos(th3)

      return d.B + d.R11 * ct3 * s₁ + d.R12 * st3 * c₁

    elseif r == 7

      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      th7 = muladd(Δth7, s, th₀7)

      st7, ct7 = sincos(th7)

      return d.B + d.R21 * ct7 * s₂ + d.R22 * st7 * c₂

    elseif r == 4

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi

      th4 = muladd(Δth4, s, th₀4)

      st4, ct4 = sincos(th4)

      return d.B + d.R11 * ct4 * s₁ + d.R12 * st4 * c₁

    else

      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      th8 = muladd(Δth8, s, th₀8)

      st8, ct8 = sincos(th8)

      return d.B + d.R21 * ct8 * s₂ + d.R22 * st8 * c₂

    end
  end
end

function fill_FTable!(tbl::FTable, d::annulus, r::Int)
  vmin, vmax, N = tbl.vmin, tbl.vmax, tbl.N
  h = (vmax - vmin) / (N - 1)

  P1, P2, P3 = tbl.P1, tbl.P2, tbl.P3

  @inbounds begin

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    if r == 1
      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      @inbounds for i in 1:N
        v̂ = vmin + (i - 1) * h
        th1 = muladd(Δth1, v̂, th₀1)
        th5 = muladd(Δth5, v̂, th₀5)

        st1, ct1 = sincos(th1)
        st5, ct5 = sincos(th5)

        y1x = d.A + d.R11 * ct1 * c₁ - d.R12 * st1 * s₁
        y1y = d.B + d.R11 * ct1 * s₁ + d.R12 * st1 * c₁
        y5x = d.A + d.R21 * ct5 * c₂ - d.R22 * st5 * s₂
        y5y = d.B + d.R21 * ct5 * s₂ + d.R22 * st5 * c₂

        Xx = (y1x + y5x) / 2
        Xy = (y1y + y5y) / 2

        Yx = y1x
        Yy = y1y

        P1[i] = Yy - Xy
        P2[i] = Yx - Xx
        P3[i] = Xy * Yx - Xx * Yy
      end

    elseif r == 5
      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      @inbounds for i in 1:N
        v̂ = vmin + (i - 1) * h
        th1 = muladd(Δth1, v̂, th₀1)
        th5 = muladd(Δth5, v̂, th₀5)

        st1, ct1 = sincos(th1)
        st5, ct5 = sincos(th5)

        y1x = d.A + d.R11 * ct1 * c₁ - d.R12 * st1 * s₁
        y1y = d.B + d.R11 * ct1 * s₁ + d.R12 * st1 * c₁
        y5x = d.A + d.R21 * ct5 * c₂ - d.R22 * st5 * s₂
        y5y = d.B + d.R21 * ct5 * s₂ + d.R22 * st5 * c₂

        Xx = (y5x + y1x) / 2
        Xy = (y5y + y1y) / 2

        Yx = y5x
        Yy = y5y

        P1[i] = Yy - Xy
        P2[i] = Yx - Xx
        P3[i] = Xy * Yx - Xx * Yy
      end

    elseif r == 2
      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      @inbounds for i in 1:N
        v̂ = vmin + (i - 1) * h

        th2 = muladd(Δth2, v̂, th₀2)
        th6 = muladd(Δth6, v̂, th₀6)

        st2, ct2 = sincos(th2)
        st6, ct6 = sincos(th6)

        y2x = d.A + d.R11 * ct2 * c₁ - d.R12 * st2 * s₁
        y2y = d.B + d.R11 * ct2 * s₁ + d.R12 * st2 * c₁
        y6x = d.A + d.R21 * ct6 * c₂ - d.R22 * st6 * s₂
        y6y = d.B + d.R21 * ct6 * s₂ + d.R22 * st6 * c₂

        Xx = (y2x + y6x) / 2
        Xy = (y2y + y6y) / 2

        Yx = y2x
        Yy = y2y

        P1[i] = Yy - Xy
        P2[i] = Yx - Xx
        P3[i] = Xy * Yx - Xx * Yy
      end

    elseif r == 6
      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      @inbounds for i in 1:N
        v̂ = vmin + (i - 1) * h

        th2 = muladd(Δth2, v̂, th₀2)
        th6 = muladd(Δth6, v̂, th₀6)

        st2, ct2 = sincos(th2)
        st6, ct6 = sincos(th6)

        y2x = d.A + d.R11 * ct2 * c₁ - d.R12 * st2 * s₁
        y2y = d.B + d.R11 * ct2 * s₁ + d.R12 * st2 * c₁
        y6x = d.A + d.R21 * ct6 * c₂ - d.R22 * st6 * s₂
        y6y = d.B + d.R21 * ct6 * s₂ + d.R22 * st6 * c₂

        Xx = (y6x + y2x) / 2
        Xy = (y6y + y2y) / 2

        Yx = y6x
        Yy = y6y

        P1[i] = Yy - Xy
        P2[i] = Yx - Xx
        P3[i] = Xy * Yx - Xx * Yy
      end

    elseif r == 3
      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      @inbounds for i in 1:N
        v̂ = vmin + (i - 1) * h

        th3 = muladd(Δth3, v̂, th₀3)
        th7 = muladd(Δth7, v̂, th₀7)

        st3, ct3 = sincos(th3)
        st7, ct7 = sincos(th7)

        y3x = d.A + d.R11 * ct3 * c₁ - d.R12 * st3 * s₁
        y3y = d.B + d.R11 * ct3 * s₁ + d.R12 * st3 * c₁
        y7x = d.A + d.R21 * ct7 * c₂ - d.R22 * st7 * s₂
        y7y = d.B + d.R21 * ct7 * s₂ + d.R22 * st7 * c₂

        Xx = (y3x + y7x) / 2
        Xy = (y3y + y7y) / 2

        Yx = y3x
        Yy = y3y

        P1[i] = Yy - Xy
        P2[i] = Yx - Xx
        P3[i] = Xy * Yx - Xx * Yy
      end

    elseif r == 7
      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      @inbounds for i in 1:N
        v̂ = vmin + (i - 1) * h

        th3 = muladd(Δth3, v̂, th₀3)
        th7 = muladd(Δth7, v̂, th₀7)

        st3, ct3 = sincos(th3)
        st7, ct7 = sincos(th7)

        y3x = d.A + d.R11 * ct3 * c₁ - d.R12 * st3 * s₁
        y3y = d.B + d.R11 * ct3 * s₁ + d.R12 * st3 * c₁
        y7x = d.A + d.R21 * ct7 * c₂ - d.R22 * st7 * s₂
        y7y = d.B + d.R21 * ct7 * s₂ + d.R22 * st7 * c₂

        Xx = (y7x + y3x) / 2
        Xy = (y7y + y3y) / 2

        Yx = y7x
        Yy = y7y

        P1[i] = Yy - Xy
        P2[i] = Yx - Xx
        P3[i] = Xy * Yx - Xx * Yy
      end
    elseif r == 4
      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      @inbounds for i in 1:N
        v̂ = vmin + (i - 1) * h

        th4 = muladd(Δth4, v̂, th₀4)
        th8 = muladd(Δth8, v̂, th₀8)

        st4, ct4 = sincos(th4)
        st8, ct8 = sincos(th8)

        y4x = d.A + d.R11 * ct4 * c₁ - d.R12 * st4 * s₁
        y4y = d.B + d.R11 * ct4 * s₁ + d.R12 * st4 * c₁
        y8x = d.A + d.R21 * ct8 * c₂ - d.R22 * st8 * s₂
        y8y = d.B + d.R21 * ct8 * s₂ + d.R22 * st8 * c₂

        Xx = (y8x + y4x) / 2
        Xy = (y8y + y4y) / 2

        Yx = y4x
        Yy = y4y

        P1[i] = Yy - Xy
        P2[i] = Yx - Xx
        P3[i] = Xy * Yx - Xx * Yy
      end

    elseif r == 8
      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      @inbounds for i in 1:N
        v̂ = vmin + (i - 1) * h

        th4 = muladd(Δth4, v̂, th₀4)
        th8 = muladd(Δth8, v̂, th₀8)

        st4, ct4 = sincos(th4)
        st8, ct8 = sincos(th8)

        y4x = d.A + d.R11 * ct4 * c₁ - d.R12 * st4 * s₁
        y4y = d.B + d.R11 * ct4 * s₁ + d.R12 * st4 * c₁
        y8x = d.A + d.R21 * ct8 * c₂ - d.R22 * st8 * s₂
        y8y = d.B + d.R21 * ct8 * s₂ + d.R22 * st8 * c₂

        Xx = (y4x + y8x) / 2
        Xy = (y4y + y8y) / 2

        Yx = y8x
        Yy = y8y

        P1[i] = Yy - Xy
        P2[i] = Yx - Xx
        P3[i] = Xy * Yx - Xx * Yy
      end

    else
      throw(ArgumentError("fill_FTable! is defined only for regions 1–8; got reg=$(p.reg)"))

    end
  end

  tbl.reg = r
  return tbl
end

@inline function f1I(t::Float64, s::Float64,
  d::annulus, u::Float64, v::Float64, r::Int)
  ŝ = (s + 1) / 2
  t̂ = (t + 1) / 2

  Xxv, Yxv = Xx(ŝ, d, r),  Yx(ŝ, d, r)

  return (1 - t̂) * Xxv + t̂ * Yxv - u
end

@inline function f2I(t::Float64, s::Float64,
  d::annulus, u::Float64, v::Float64, r::Int)
    ŝ = (s + 1) / 2
    t̂ = (t + 1) / 2

    Xyv, Yyv = Xy(ŝ, d, r), Yy(ŝ, d, r)

    return (1 - t̂) * Xyv + t̂ * Yyv - v
end

@inline function JinvI(t::Float64, s::Float64,
  d::annulus, u::Float64, v::Float64, r::Int)
  return jinvmap(d, t, s, r)  # returns J11,J12,J21,J22
end

@inline function f_cont(v̂::Float64,
  d::annulus, u::Float64, v::Float64, r::Int)

  @inbounds begin

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    if r == 1

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      th1 = muladd(Δth1, v̂, th₀1)
      th5 = muladd(Δth5, v̂, th₀5)

      st1, ct1 = sincos(th1)
      st5, ct5 = sincos(th5)

      y1x = d.A + d.R11 * ct1 * c₁ - d.R12 * st1 * s₁
      y1y = d.B + d.R11 * ct1 * s₁ + d.R12 * st1 * c₁
      y5x = d.A + d.R21 * ct5 * c₂ - d.R22 * st5 * s₂
      y5y = d.B + d.R21 * ct5 * s₂ + d.R22 * st5 * c₂

      Xxv = (y1x + y5x) / 2
      Xyv = (y1y + y5y) / 2

      Yxv = y1x
      Yyv = y1y

    elseif r == 5

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      th1 = muladd(Δth1, v̂, th₀1)
      th5 = muladd(Δth5, v̂, th₀5)

      st1, ct1 = sincos(th1)
      st5, ct5 = sincos(th5)

      y1x = d.A + d.R11 * ct1 * c₁ - d.R12 * st1 * s₁
      y1y = d.B + d.R11 * ct1 * s₁ + d.R12 * st1 * c₁
      y5x = d.A + d.R21 * ct5 * c₂ - d.R22 * st5 * s₂
      y5y = d.B + d.R21 * ct5 * s₂ + d.R22 * st5 * c₂

      Xxv = (y5x + y1x) / 2
      Xyv = (y5y + y1y) / 2

      Yxv = y5x
      Yyv = y5y

    elseif r == 2

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      th2 = muladd(Δth2, v̂, th₀2)
      th6 = muladd(Δth6, v̂, th₀6)

      st2, ct2 = sincos(th2)
      st6, ct6 = sincos(th6)

      y2x = d.A + d.R11 * ct2 * c₁ - d.R12 * st2 * s₁
      y2y = d.B + d.R11 * ct2 * s₁ + d.R12 * st2 * c₁
      y6x = d.A + d.R21 * ct6 * c₂ - d.R22 * st6 * s₂
      y6y = d.B + d.R21 * ct6 * s₂ + d.R22 * st6 * c₂

      Xxv = (y2x + y6x) / 2
      Xyv = (y2y + y6y) / 2

      Yxv = y2x
      Yyv = y2y

    elseif r == 6

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      th2 = muladd(Δth2, v̂, th₀2)
      th6 = muladd(Δth6, v̂, th₀6)

      st2, ct2 = sincos(th2)
      st6, ct6 = sincos(th6)

      y2x = d.A + d.R11 * ct2 * c₁ - d.R12 * st2 * s₁
      y2y = d.B + d.R11 * ct2 * s₁ + d.R12 * st2 * c₁
      y6x = d.A + d.R21 * ct6 * c₂ - d.R22 * st6 * s₂
      y6y = d.B + d.R21 * ct6 * s₂ + d.R22 * st6 * c₂

      Xxv = (y6x + y2x) / 2
      Xyv = (y6y + y2y) / 2

      Yxv = y6x
      Yyv = y6y

    elseif r == 3

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      th3 = muladd(Δth3, v̂, th₀3)
      th7 = muladd(Δth7, v̂, th₀7)

      st3, ct3 = sincos(th3)
      st7, ct7 = sincos(th7)

      y3x = d.A + d.R11 * ct3 * c₁ - d.R12 * st3 * s₁
      y3y = d.B + d.R11 * ct3 * s₁ + d.R12 * st3 * c₁
      y7x = d.A + d.R21 * ct7 * c₂ - d.R22 * st7 * s₂
      y7y = d.B + d.R21 * ct7 * s₂ + d.R22 * st7 * c₂

      Xxv = (y3x + y7x) / 2
      Xyv = (y3y + y7y) / 2

      Yxv = y3x
      Yyv = y3y

    elseif r == 7

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      th3 = muladd(Δth3, v̂, th₀3)
      th7 = muladd(Δth7, v̂, th₀7)

      st3, ct3 = sincos(th3)
      st7, ct7 = sincos(th7)

      y3x = d.A + d.R11 * ct3 * c₁ - d.R12 * st3 * s₁
      y3y = d.B + d.R11 * ct3 * s₁ + d.R12 * st3 * c₁
      y7x = d.A + d.R21 * ct7 * c₂ - d.R22 * st7 * s₂
      y7y = d.B + d.R21 * ct7 * s₂ + d.R22 * st7 * c₂

      Xxv = (y7x + y3x) / 2
      Xyv = (y7y + y3y) / 2

      Yxv = y7x
      Yyv = y7y

    elseif r == 4

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      th4 = muladd(Δth4, v̂, th₀4)
      th8 = muladd(Δth8, v̂, th₀8)

      st4, ct4 = sincos(th4)
      st8, ct8 = sincos(th8)

      y4x = d.A + d.R11 * ct4 * c₁ - d.R12 * st4 * s₁
      y4y = d.B + d.R11 * ct4 * s₁ + d.R12 * st4 * c₁
      y8x = d.A + d.R21 * ct8 * c₂ - d.R22 * st8 * s₂
      y8y = d.B + d.R21 * ct8 * s₂ + d.R22 * st8 * c₂
      
      Xxv = (y4x + y8x) / 2
      Xyv = (y4y + y8y) / 2

      Yxv = y4x
      Yyv = y4y

    else

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      th4 = muladd(Δth4, v̂, th₀4)
      th8 = muladd(Δth8, v̂, th₀8)

      st4, ct4 = sincos(th4)
      st8, ct8 = sincos(th8)

      y4x = d.A + d.R11 * ct4 * c₁ - d.R12 * st4 * s₁
      y4y = d.B + d.R11 * ct4 * s₁ + d.R12 * st4 * c₁
      y8x = d.A + d.R21 * ct8 * c₂ - d.R22 * st8 * s₂
      y8y = d.B + d.R21 * ct8 * s₂ + d.R22 * st8 * c₂

      Yxv = y8x
      Yyv = y8y

      Xxv = (y8x + y4x) / 2
      Xyv = (y8y + y4y) / 2

    end
  end

  return u * (Yyv - Xyv) - v * (Yxv - Xxv) + (Xyv*Yxv - Xxv*Yyv)
end

function mapinv(tbl::FTable, d::annulus, u::Float64, 
  v::Float64, k::Int)::Tuple{Float64,Float64}

  p = d.pths[k]

  r = p.reg

  # ----- curved patches reg = 1..8 -----
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
    ptconv(d::annulus, t1::Float64, t2::Float64, idx::Int, ptdest::String)

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
function ptconv(d::annulus, t1::Float64, t2::Float64, idx::Int, ptdest::String)
  #@assert length(t) == 2 "t must be length-2 vector [t1,t2]"
  if ptdest == "to_pth"
    r = idx

    for k in 1:d.Npat
      p = d.pths[k]
      p.reg == r || continue

      t1k = xi_inv((t1 + 1) / 2, p.ck0, p.ck1)
      t2k = xi_inv((t2 + 1) / 2, p.tk0, p.tk1)

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

    t1r = (p.ck1 - p.ck0) * t1 + (p.ck1 + p.ck0) - 1
    t2r = (p.tk1 - p.tk0) * t2 + (p.tk1 + p.tk0) - 1

    return t1r, t2r, r

  else
    error("ptconv: ptdest must be \"to_pth\" or \"to_reg\"")
  end
end

"""
    mapinv2(d::annulus, t1::Float64, t2::Float64, k2::Int, k::Int) -> Tuple{Float64,Float64}

Given a point `(u,v) = τₖ₂(t1,t2)` on patch `k2`, return its reference
coordinates on patch `k` in `[-1,1]^2`.

This differs from `mapinv` because we already know `(u,v)` comes from `(t1,t2)` on `k2`.
"""
function mapinv2(d::annulus, t1::Float64, t2::Float64, k2::Int, k::Int)::Tuple{Float64,Float64}
  p = d.pths[k]  # Fields reg, ck0, ck1, tk0, tk1

  # Convert the given (t1,t2,k2) to the "regional" normalized coords of the *point*
  # Expecting something that returns two values in [-1,1]:
  # tr1 ≈ u_ref, tr2 ≈ v_ref  (for the *global/regional* parameterization)
  tr1, tr2, _ = ptconv(d, t1, t2, k2, "to_reg")

  # Map from [-1,1] -> [0,1]
  û = (tr1 + 1) / 2
  v̂ = (tr2 + 1) / 2

  t1k = xi_inv(û, p.ck0, p.ck1)
  t2k = xi_inv(v̂, p.tk0, p.tk1)

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
function dfunc(d::annulus, k::Int, t::Float64, s::Float64)::Float64

  p = d.pths[k]
  hc = p.ck1 - p.ck0
  val = (1.0 - p.ck1) + hc * t / 2
  exp = s ≥ 0.5 ? (s - 1) : s
  return val^exp

end

function dfunc!(out::StridedArray{Float64}, d::annulus, k::Int, t::StridedArray{Float64}, s::Float64)

  p = d.pths[k]
  exp = s ≥ 0.5 ? (s - 1) : s
  αc = (p.ck1 - p.ck0) / 2
  βc = 1.0 - p.ck1
  @inbounds for i in eachindex(t)
    out[i] = (muladd(αc, t[i], βc))^exp
  end

  return nothing

end

"""
  diff_map!(out::Matrix{Float64}, Zx::Matrix{Float64}, Zy::Matrix{Float64}, 
  DJ::Matrix{Float64}, d::annulus, u::Float64, v::Float64,
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
  d::annulus, u::Float64, v::Float64,
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
  reg = p.reg

  # affine scalings for reference -> (ck0,ck1)/(tk0,tk1)
  hc = p.ck1 - p.ck0
  ht = p.tk1 - p.tk0
  αu = 0.5 * hc
  αv = 0.5 * ht
  αuv= αu * αv
  βu = p.ck0 + αu
  βv = p.tk0 + αv

  mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)

  # --- Precomputed geometry params
  c₁ = d.RP.Cθ₁
  s₁ = d.RP.Sθ₁
  c₂ = d.RP.Cθ₂
  s₂ = d.RP.Sθ₂

  # ------------------------
  # General case reg ∈ 1:8
  # ------------------------
  R11s1, R11c1 = d.R11 * s₁, d.R11 * c₁
  R12s1, R12c1 = d.R12 * s₁, d.R12 * c₁

  R21s2, R21c2 = d.R21 * s₂, d.R21 * c₂
  R22s2, R22c2 = d.R22 * s₂, d.R22 * c₂

  mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)

  # τ(u,v) scalar
  tux, tvy = mapxy(d, u, v, k)

  
  xi1 = muladd(αu, u, βu)
  xi2 = muladd(αv, v, βv)

  #a is in the set {1,2,3,4}
  #b is in the set {5,6,7,8}
  if reg == 1 || reg == 5
    th0a, Δtha = d.T1[1], d.T1[2] - d.T1[1]
    th0b, Δthb = d.T2[1], d.T2[2] - d.T2[1]
  elseif reg == 2 || reg == 6
    th0a, Δtha = d.T1[2], d.T1[3] - d.T1[2]
    th0b, Δthb = d.T2[2], d.T2[3] - d.T2[2]
  elseif reg == 3 || reg == 7
    th0a, Δtha = d.T1[3], d.T1[4] - d.T1[3]
    th0b, Δthb = d.T2[3], d.T2[4] - d.T2[3]
  elseif reg == 4 || reg == 8
    th0a, Δtha = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
    th0b, Δthb = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
  else
    throw(ArgumentError("diff_map! Taylor fixup expects reg in 1–8; got reg=$reg"))
  end

  # angles at xi2
  tha = muladd(Δtha, xi2, th0a)
  thb = muladd(Δthb, xi2, th0b)
  sta, cta = sincos(tha)
  stb, ctb = sincos(thb)

  # curve a: uses (R11,R12,θ₁), curve b: uses (R21,R22,θ₂)
  # x components
  yax = muladd(-sta, R12s1, muladd(cta, R11c1, d.A))
  ybx = muladd(-stb, R22s2, muladd(ctb, R21c2, d.A))
  # y components
  yay = muladd(sta, R12c1, muladd(cta, R11s1, d.B))
  yby = muladd(stb, R22c2, muladd(ctb, R21s2, d.B))

  # first derivatives w.r.t. xi2 (v-mapped)
  dyax = Δtha * (-R11c1 * sta - R12s1 * cta)
  dyay = Δtha * (-R11s1 * sta + R12c1 * cta)

  dybx = Δthb * (-R21c2 * stb - R22s2 * ctb)
  dyby = Δthb * (-R21s2 * stb + R22c2 * ctb)

  # Higher derivatives: exactly as MATLAB
  # d2 = -(Δth)^2 * (y - center), 
  # d3 = -(Δth)^2 * d1,
  # d4 = -(Δth)^2 * d2
  d2yax = -Δtha^2 * (yax - d.A)
  d2yay = -Δtha^2 * (yay - d.B)
  d3yax = -Δtha^2 * dyax
  d3yay = -Δtha^2 * dyay
  d4yax = -Δtha^2 * d2yax
  d4yay = -Δtha^2 * d2yay

  d2ybx = -Δthb^2 * (ybx - d.A)
  d2yby = -Δthb^2 * (yby - d.B)
  d3ybx = -Δthb^2 * dybx
  d3yby = -Δthb^2 * dyby
  d4ybx = -Δthb^2 * d2ybx
  d4yby = -Δthb^2 * d2yby

  if reg <= 4
    dux = αu * (yax - ybx) / 2
    duvx = αuv * (dyax - dybx) / 2
    duv2x = αv * αuv * (d2yax - d2ybx) / 2
    duv3x = αv^2 * αuv * (d3yax - d3ybx) / 2

    dvx = αv * ((1 + xi1) * dyax + (1 - xi1) * dybx) / 2
    dv2x = αv^2 * ((1 + xi1) * d2yax + (1 - xi1) * d2ybx) / 2
    dv3x = αv^3 * ((1 + xi1) * d3yax + (1 - xi1) * d3ybx) / 2
    dv4x = αv^4 * ((1 + xi1) * d4yax + (1 - xi1) * d4ybx) / 2

    duy = αu * (yay - yby) / 2
    duvy = αuv * (dyay - dyby) / 2
    duv2y = αv * αuv * (d2yay - d2yby) / 2
    duv3y = αv^2 * αuv * (d3yay - d3yby) / 2

    dvy = αv * ((1 + xi1) * dyay + (1 - xi1) * dyby) / 2
    dv2y = αv^2 * ((1 + xi1) * d2yay + (1 - xi1) * d2yby) / 2
    dv3y = αv^3 * ((1 + xi1) * d3yay + (1 - xi1) * d3yby) / 2
    dv4y = αv^4 * ((1 + xi1) * d4yay + (1 - xi1) * d4yby) / 2

  else
    dux = αu * (ybx - yax) / 2
    duvx = αuv * (dybx - dyax) / 2
    duv2x = αv * αuv * (d2ybx - d2yax) / 2
    duv3x = αv^2 * αuv * (d3ybx - d3yax) / 2

    dvx = αv * ((1 + xi1) * dybx + (1 - xi1) * dyax) / 2
    dv2x = αv^2 * ((1 + xi1) * d2ybx + (1 - xi1) * d2yax) / 2
    dv3x = αv^3 * ((1 + xi1) * d3ybx + (1 - xi1) * d3yax) / 2
    dv4x = αv^4 * ((1 + xi1) * d4ybx + (1 - xi1) * d4yax) / 2

    duy = αu * (yby - yay) / 2
    duvy = αuv * (dyby - dyay) / 2
    duv2y = αv * αuv * (d2yby - d2yay) / 2
    duv3y = αv^2 * αuv * (d3yby - d3yay) / 2

    dvy = αv * ((1 + xi1) * dyby + (1 - xi1) * dyay) / 2
    dv2y = αv^2 * ((1 + xi1) * d2yby + (1 - xi1) * d2yay) / 2
    dv3y = αv^3 * ((1 + xi1) * d3yby + (1 - xi1) * d3yay) / 2
    dv4y = αv^4 * ((1 + xi1) * d4yby + (1 - xi1) * d4yay) / 2
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
  DJ::StridedArray{Float64}, d::annulus, u::Float64, v::Float64,
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
with affine mappings are not updated! (no affine map in annulus)
But the Jacobian DJ is always updated.
"""
function diff_rmap!(out::Matrix{Float64},
  Zx::Matrix{Float64}, Zy::Matrix{Float64}, DJ::StridedArray{Float64},
  d::annulus, u::Float64, v::Float64,
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
  hc = p.ck1 - p.ck0
  ht = p.tk1 - p.tk0
  αu = 0.5 * hc
  αv = 0.5 * ht
  αuv = αu * αv
  βu = p.ck0 + αu
  βv = p.tk0 + αv

  mapxy_Dmap!(Zx, Zy, DJ, d, u2, v2, k)

  # --- Precomputed geometry params
  c₁ = d.RP.Cθ₁
  s₁ = d.RP.Sθ₁
  c₂ = d.RP.Cθ₂
  s₂ = d.RP.Sθ₂

  # ------------------------
  # General case reg ∈ 1:8
  # ------------------------
  R11s1, R11c1 = d.R11 * s₁, d.R11 * c₁
  R12s1, R12c1 = d.R12 * s₁, d.R12 * c₁

  R21s2, R21c2 = d.R21 * s₂, d.R21 * c₂
  R22s2, R22c2 = d.R22 * s₂, d.R22 * c₂

  # τ(u,v) scalar
  tux, tvy = mapxy(d, u, v, k)

  # mapped scalar (u,v) to (xi1,xi2)
  xi1 = muladd(αu, u, βu)
  xi2 = muladd(αv, v, βv)

  #a is in the set {1,2,3,4}
  #b is in the set {5,6,7,8}
  if reg == 1 || reg == 5
    th0a, Δtha = d.T1[1], d.T1[2] - d.T1[1]
    th0b, Δthb = d.T2[1], d.T2[2] - d.T2[1]
  elseif reg == 2 || reg == 6
    th0a, Δtha = d.T1[2], d.T1[3] - d.T1[2]
    th0b, Δthb = d.T2[2], d.T2[3] - d.T2[2]
  elseif reg == 3 || reg == 7
    th0a, Δtha = d.T1[3], d.T1[4] - d.T1[3]
    th0b, Δthb = d.T2[3], d.T2[4] - d.T2[3]
  elseif reg == 4 || reg == 8
    th0a, Δtha = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
    th0b, Δthb = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi
  else
    throw(ArgumentError("diff_map! Taylor fixup expects reg in 1–8; got reg=$reg"))
  end

  # angles at xi2
  tha = muladd(Δtha, xi2, th0a)
  thb = muladd(Δthb, xi2, th0b)
  sta, cta = sincos(tha)
  stb, ctb = sincos(thb)

  # curve a: uses (R11,R12,θ₁), curve b: uses (R21,R22,θ₂)
  # x components
  yax = muladd(-sta, R12s1, muladd(cta, R11c1, d.A))
  ybx = muladd(-stb, R22s2, muladd(ctb, R21c2, d.A))
  # y components
  yay = muladd(sta, R12c1, muladd(cta, R11s1, d.B))
  yby = muladd(stb, R22c2, muladd(ctb, R21s2, d.B))

  # first derivatives w.r.t. xi2 (v-mapped)
  dyax = Δtha * (-R11c1 * sta - R12s1 * cta)
  dyay = Δtha * (-R11s1 * sta + R12c1 * cta)

  dybx = Δthb * (-R21c2 * stb - R22s2 * ctb)
  dyby = Δthb * (-R21s2 * stb + R22c2 * ctb)

  # Higher derivatives: exactly as MATLAB
  # d2 = -(Δth)^2 * (y - center), 
  # d3 = -(Δth)^2 * d1,
  # d4 = -(Δth)^2 * d2
  d2yax = -Δtha^2 * (yax - d.A)
  d2yay = -Δtha^2 * (yay - d.B)
  d3yax = -Δtha^2 * dyax
  d3yay = -Δtha^2 * dyay
  d4yax = -Δtha^2 * d2yax
  d4yay = -Δtha^2 * d2yay

  d2ybx = -Δthb^2 * (ybx - d.A)
  d2yby = -Δthb^2 * (yby - d.B)
  d3ybx = -Δthb^2 * dybx
  d3yby = -Δthb^2 * dyby
  d4ybx = -Δthb^2 * d2ybx
  d4yby = -Δthb^2 * d2yby

  if reg <= 4
    dux = αu * (yax - ybx) / 2
    duvx = αuv * (dyax - dybx) / 2
    duv2x = αv * αuv * (d2yax - d2ybx) / 2
    duv3x = αv^2 * αuv * (d3yax - d3ybx) / 2

    dvx = αv * ((1 + xi1) * dyax + (1 - xi1) * dybx) / 2
    dv2x = αv^2 * ((1 + xi1) * d2yax + (1 - xi1) * d2ybx) / 2
    dv3x = αv^3 * ((1 + xi1) * d3yax + (1 - xi1) * d3ybx) / 2
    dv4x = αv^4 * ((1 + xi1) * d4yax + (1 - xi1) * d4ybx) / 2

    duy = αu * (yay - yby) / 2
    duvy = αuv * (dyay - dyby) / 2
    duv2y = αv * αuv * (d2yay - d2yby) / 2
    duv3y = αv^2 * αuv * (d3yay - d3yby) / 2

    dvy = αv * ((1 + xi1) * dyay + (1 - xi1) * dyby) / 2
    dv2y = αv^2 * ((1 + xi1) * d2yay + (1 - xi1) * d2yby) / 2
    dv3y = αv^3 * ((1 + xi1) * d3yay + (1 - xi1) * d3yby) / 2
    dv4y = αv^4 * ((1 + xi1) * d4yay + (1 - xi1) * d4yby) / 2

  else
    dux = αu * (ybx - yax) / 2
    duvx = αuv * (dybx - dyax) / 2
    duv2x = αv * αuv * (d2ybx - d2yax) / 2
    duv3x = αv^2 * αuv * (d3ybx - d3yax) / 2

    dvx = αv * ((1 + xi1) * dybx + (1 - xi1) * dyax) / 2
    dv2x = αv^2 * ((1 + xi1) * d2ybx + (1 - xi1) * d2yax) / 2
    dv3x = αv^3 * ((1 + xi1) * d3ybx + (1 - xi1) * d3yax) / 2
    dv4x = αv^4 * ((1 + xi1) * d4ybx + (1 - xi1) * d4yax) / 2

    duy = αu * (yby - yay) / 2
    duvy = αuv * (dyby - dyay) / 2
    duv2y = αv * αuv * (d2yby - d2yay) / 2
    duv3y = αv^2 * αuv * (d3yby - d3yay) / 2

    dvy = αv * ((1 + xi1) * dyby + (1 - xi1) * dyay) / 2
    dv2y = αv^2 * ((1 + xi1) * d2yby + (1 - xi1) * d2yay) / 2
    dv3y = αv^3 * ((1 + xi1) * d3yby + (1 - xi1) * d3yay) / 2
    dv4y = αv^4 * ((1 + xi1) * d4yby + (1 - xi1) * d4yay) / 2
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
    Dwall(d::annulus, u::Float64, v::Float64, k::Int) -> dvx, dvy

Derivative of the wall curve w.r.t. `v` at the singular side `u = ±1`
for the `k`-th (non-quadrilateral) patch.
Returns a Tuple `dvx, dvy`.
"""
function Dwall(d::annulus, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64}
  @inbounds begin
    # p denotes information of kth patch
    p = d.pths[k]

    hc = p.ck1 - p.ck0
    ht = p.tk1 - p.tk0
    αu = 0.5 * hc
    αv = 0.5 * ht
    βu = p.ck0 + αu
    βv = p.tk0 + αv

    c₁ = d.RP.Cθ₁
    s₁ = d.RP.Sθ₁
    c₂ = d.RP.Cθ₂
    s₂ = d.RP.Sθ₂

    R11s1, R11c1 = d.R11 * s₁, d.R11 * c₁
    R12s1, R12c1 = d.R12 * s₁, d.R12 * c₁

    R21s2, R21c2 = d.R21 * s₂, d.R21 * c₂
    R22s2, R22c2 = d.R22 * s₂, d.R22 * c₂

    if p.reg == 1

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      xi1 = muladd(αu, u, βu)
      xi2 = muladd(αv, v, βv)

      th1 = muladd(Δth1, xi2, th₀1)
      th5 = muladd(Δth5, xi2, th₀5)

      st1, ct1 = sincos(th1)
      st5, ct5 = sincos(th5)

      dy1x = Δth1 * (-R11c1 * st1 - R12s1 * ct1)
      dy1y = Δth1 * (-R11s1 * st1 + R12c1 * ct1)

      dy5x = Δth5 * (-R21c2 * st5 - R22s2 * ct5)
      dy5y = Δth5 * (-R21s2 * st5 + R22c2 * ct5)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dvx = tp * dy1x + tm * dy5x
      dvy = tp * dy1y + tm * dy5y

    elseif p.reg == 5

      th₀1, Δth1 = d.T1[1], d.T1[2] - d.T1[1]
      th₀5, Δth5 = d.T2[1], d.T2[2] - d.T2[1]

      xi1 = muladd(αu, u, βu)
      xi2 = muladd(αv, v, βv)

      th1 = muladd(Δth1, xi2, th₀1)
      th5 = muladd(Δth5, xi2, th₀5)

      st1, ct1 = sincos(th1)
      st5, ct5 = sincos(th5)

      dy1x = Δth1 * (-R11c1 * st1 - R12s1 * ct1)
      dy1y = Δth1 * (-R11s1 * st1 + R12c1 * ct1)

      dy5x = Δth5 * (-R21c2 * st5 - R22s2 * ct5)
      dy5y = Δth5 * (-R21s2 * st5 + R22c2 * ct5)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dvx = tp * dy5x + tm * dy1x
      dvy = tp * dy5y + tm * dy1y

    elseif p.reg == 2

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      xi1 = muladd(αu, u, βu)
      xi2 = muladd(αv, v, βv)

      th2 = muladd(Δth2, xi2, th₀2)
      th6 = muladd(Δth6, xi2, th₀6)

      st2, ct2 = sincos(th2)
      st6, ct6 = sincos(th6)

      dy2x = Δth2 * (-R11c1 * st2 - R12s1 * ct2)
      dy2y = Δth2 * (-R11s1 * st2 + R12c1 * ct2)

      dy6x = Δth6 * (-R21c2 * st6 - R22s2 * ct6)
      dy6y = Δth6 * (-R21s2 * st6 + R22c2 * ct6)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dvx = tp * dy2x + tm * dy6x
      dvy = tp * dy2y + tm * dy6y

    elseif p.reg == 6

      th₀2, Δth2 = d.T1[2], d.T1[3] - d.T1[2]
      th₀6, Δth6 = d.T2[2], d.T2[3] - d.T2[2]

      xi1 = muladd(αu, u, βu)
      xi2 = muladd(αv, v, βv)

      th2 = muladd(Δth2, xi2, th₀2)
      th6 = muladd(Δth6, xi2, th₀6)

      st2, ct2 = sincos(th2)
      st6, ct6 = sincos(th6)

      dy2x = Δth2 * (-R11c1 * st2 - R12s1 * ct2)
      dy2y = Δth2 * (-R11s1 * st2 + R12c1 * ct2)

      dy6x = Δth6 * (-R21c2 * st6 - R22s2 * ct6)
      dy6y = Δth6 * (-R21s2 * st6 + R22c2 * ct6)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dvx = tp * dy6x + tm * dy2x
      dvy = tp * dy6y + tm * dy2y

    elseif p.reg == 3

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      xi1 = muladd(αu, u, βu)
      xi2 = muladd(αv, v, βv)

      th3 = muladd(Δth3, xi2, th₀3)
      th7 = muladd(Δth7, xi2, th₀7)

      st3, ct3 = sincos(th3)
      st7, ct7 = sincos(th7)

      dy3x = Δth3 * (-R11c1 * st3 - R12s1 * ct3)
      dy3y = Δth3 * (-R11s1 * st3 + R12c1 * ct3)

      dy7x = Δth7 * (-R21c2 * st7 - R22s2 * ct7)
      dy7y = Δth7 * (-R21s2 * st7 + R22c2 * ct7)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dvx = tp * dy3x + tm * dy7x
      dvy = tp * dy3y + tm * dy7y

    elseif p.reg == 7

      th₀3, Δth3 = d.T1[3], d.T1[4] - d.T1[3]
      th₀7, Δth7 = d.T2[3], d.T2[4] - d.T2[3]

      xi1 = muladd(αu, u, βu)
      xi2 = muladd(αv, v, βv)

      th3 = muladd(Δth3, xi2, th₀3)
      th7 = muladd(Δth7, xi2, th₀7)

      st3, ct3 = sincos(th3)
      st7, ct7 = sincos(th7)

      dy3x = Δth3 * (-R11c1 * st3 - R12s1 * ct3)
      dy3y = Δth3 * (-R11s1 * st3 + R12c1 * ct3)

      dy7x = Δth7 * (-R21c2 * st7 - R22s2 * ct7)
      dy7y = Δth7 * (-R21s2 * st7 + R22c2 * ct7)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dvx = tp * dy7x + tm * dy3x
      dvy = tp * dy7y + tm * dy3y

    elseif p.reg == 4

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      xi1 = muladd(αu, u, βu)
      xi2 = muladd(αv, v, βv)

      th4 = muladd(Δth4, xi2, th₀4)
      th8 = muladd(Δth8, xi2, th₀8)

      st4, ct4 = sincos(th4)
      st8, ct8 = sincos(th8)

      dy4x = Δth4 * (-R11c1 * st4 - R12s1 * ct4)
      dy4y = Δth4 * (-R11s1 * st4 + R12c1 * ct4)

      dy8x = Δth8 * (-R21c2 * st8 - R22s2 * ct8)
      dy8y = Δth8 * (-R21s2 * st8 + R22c2 * ct8)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dvx = tp * dy4x + tm * dy8x
      dvy = tp * dy4y + tm * dy8y

    elseif p.reg == 8

      th₀4, Δth4 = d.T1[4], d.T1[1] - d.T1[4] + 2.0 * pi
      th₀8, Δth8 = d.T2[4], d.T2[1] - d.T2[4] + 2.0 * pi

      xi1 = muladd(αu, u, βu)
      xi2 = muladd(αv, v, βv)

      th4 = muladd(Δth4, xi2, th₀4)
      th8 = muladd(Δth8, xi2, th₀8)

      st4, ct4 = sincos(th4)
      st8, ct8 = sincos(th8)

      dy4x = Δth4 * (-R11c1 * st4 - R12s1 * ct4)
      dy4y = Δth4 * (-R11s1 * st4 + R12c1 * ct4)

      dy8x = Δth8 * (-R21c2 * st8 - R22s2 * ct8)
      dy8y = Δth8 * (-R21s2 * st8 + R22c2 * ct8)

      tp = muladd(0.5, xi1, 0.5)
      tm = muladd(0.5, -xi1, 0.5)

      dvx = tp * dy8x + tm * dy4x
      dvy = tp * dy8y + tm * dy4y

    else
      throw(ArgumentError("Dwall is defined only for regions 1–8; got reg=$(p.reg)"))

    end
  end

  return αv * dvx, αv * dvy
end

"""
    boundquad!(P::SubArray{Float64}, d::annulus, k::Int) -> nothing

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
function boundquad!(P::SubArray{Float64}, d::annulus, k::Int)

  p = d.pths[k]

  # Four corners (tuples immediately unpacked into scalars)
  L1p1x, L1p1y = mapm1(d, -1.0, k)   # (u=-1, v=-1)
  L1p2x, L1p2y = mapp1(d, -1.0, k)   # (u=+1, v=-1)
  L2p1x, L2p1y = mapm1(d, 1.0, k)   # (u=-1, v=+1)
  L2p2x, L2p2y = mapp1(d, 1.0, k)   # (u=+1, v=+1)

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
    refine!(d::annulus, Nc::Int, Nt::Int, K::AbstractVector{<:Int})

Subdivide each patch in `K` into `Nc * Nt` subpatches (split in `c` and `t`),
update `d.pths`, `d.Npat`, recompute `d.kd`, and rebuild `d.Qpts` / `d.Qptsbd`.

Notes
- Patches are re-sorted lexicographically by `(reg, ck0, ck1, tk0, tk1)`
- `Qpts` uses `boundquad(d,k)`; `Qptsbd` uses `boundquadbd(d,k)` for patches in `d.kd`.
"""
function refine!(d::annulus, Nc::Int, Nt::Int, K::Vector{Int})
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
    d.kd = [k for k in 1:d.Npat if d.pths[k].ck1 == 1.0]

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

function Base.show(io::IO, d::annulus)
  println(io, "annulus with properties:")
  println(io, " (A , B )    = (", d.A,",", d.B,")")
  println(io, " (R₁₁, R₁₂)  = (", d.R11,",", d.R12,")")
  println(io, " (R₂₁, R₂₂)  = (", d.R21,",", d.R22,")")
  println(io, " (θ₁ , θ₂ )  = (", d.tht1,",", d.tht2,")")
  #println(io, " (cos₁,sin₁) = (", d.RP.Cθ₁,",", d.RP.Sθ₁,")")
  #println(io, " (cos₂,sin₂) = (", d.RP.Cθ₂,",", d.RP.Sθ₂,")")
  println(io, "  No. of holes= ", d.nh)

  if length(d.kd) <= 6
    L = "Int[" * join(d.kd, ' ') * "]"
  else
    head = join(d.kd[1:3], ' ')
    tail = join(d.kd[end-3+1:end], ' ')
    L = "Int[$head … $tail]"
  end
  
  print(io, "  T₁ = [")
  @inbounds for i in eachindex(d.T1)
    i > firstindex(d.T1) && print(io, ", ")
    @printf(io, "%.3f", d.T1[i])
  end
  println(io, "]")

  print(io, "  T₂ = [")
  @inbounds for i in eachindex(d.T2)
    i > firstindex(d.T2) && print(io, ", ")
    @printf(io, "%.3f", d.T2[i])
  end
  println(io, "]")

  println(io, "  kd     = ", L)
  println(io, "  Npat   = ", d.Npat)
  println(io, "  pths   = Vector{Patch} (", length(d.pths), " patches)")
  println(io, "  Qpts   = ", size(d.Qpts, 1), "×", size(d.Qpts, 2), " Matrix")
  println(io, "  Qptsbd = ", size(d.Qptsbd, 1), "×", size(d.Qptsbd, 2), " Matrix")

end
