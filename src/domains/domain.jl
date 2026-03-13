struct Patch
  reg::Int
  ck0::Float64
  ck1::Float64
  tk0::Float64
  tk1::Float64
end

mutable struct FTable
    vmin::Float64
    vmax::Float64
    N::Int
    P1::Vector{Float64}
    P2::Vector{Float64}
    P3::Vector{Float64}
    #roots in mapinv function
    rmi::Vector{Float64}
    zxi::Vector{Float64}
    reg::Int    
end

function inFTable(N::Int=100_001)
    vmin = -30.0
    vmax =  30.0

    P1   = Vector{Float64}(undef, N)
    P2   = Vector{Float64}(undef, N)
    P3   = Vector{Float64}(undef, N)
    rmi  = Vector{Float64}(undef, 256)
    zxi  = Vector{Float64}(undef, 256)
    # reg=0 = "uninitialized"
    return FTable(vmin, vmax, N, P1, P2, P3, rmi, zxi, 0) 
end

@inline function xi_inv(x::Float64, a::Float64, b::Float64)::Float64
    half = 0.5*(b - a)
    mid  = a + half
    return (x - mid)/half
end

include(joinpath(@__DIR__, "disc.jl")); 
include(joinpath(@__DIR__, "ellipse.jl")); 
include(joinpath(@__DIR__, "kite.jl")); 
include(joinpath(@__DIR__, "annulus.jl")); 

"""
    find_roots!(roots, tbl, d, p, u, v;
                maxiter=48, tol=1e-12)

- `roots`  : preallocated Vector{Float64} to store found roots
- returns  : number of roots written into `roots[1:n]`
"""
function find_roots!(roots::Vector{Float64},
  tbl::FTable, d::abstractdomain, u::Float64, v::Float64, 
  r::Int; maxi::Int=128, tol::Float64=1e-16)

  vmin, vmax, N = tbl.vmin, tbl.vmax, tbl.N
  P1, P2, P3 = tbl.P1, tbl.P2, tbl.P3
  h = (vmax - vmin) / (N - 1)

  nroots = 0

  # F at first grid point
  v_prev = vmin

  F_prev = u * P1[1] - v * P2[1] + P3[1]

  @inbounds for i in 2:N
    v_curr = vmin + (i - 1) * h
    F_curr = u * P1[i] - v * P2[i] + P3[i]

    # Check sign change or exact zero
    if abs(F_curr) < 1e-14
      # root exactly on grid node
        nroots += 1
        roots[nroots] = v_curr
    elseif F_prev * F_curr < 0.0
      # bracket [v_prev, v_curr]
        z = Bis(f_cont, v_prev, v_curr, maxi,
          d, u, v, r; tol=tol)
        nroots += 1
        roots[nroots] = z
    end

    v_prev = v_curr
    F_prev = F_curr
  end

  return nroots
end

# Non-mutating convenience 
function mapxy(d::abstractdomain, u::AbstractArray, v::AbstractArray, k::Int)
    Zx = similar(u, Float64)
    Zy = similar(v, Float64)
    mapxy!(Zx, Zy, d, u, v, k)
    return Zx, Zy
end

"""
    gam(d::abstractdomain, t::Float64, k::Int)
    gam(d::abstractdomain, t::Float64)
    gam(d::abstractdomain, t::AbstractArray, k::Int)
    gam(d::abstractdomain, t::AbstractArray)

Combine `gamx` and `gamy` vertically to form the 2D parametrization
γ(t) = [γₓ(t); γᵧ(t)].

- With `k`: gives the boundary parametrization of patch `k`.
- Without `k`: gives the global boundary.

`t` can be a scalar or an array; the output is a 2×N array if `t` is an array,
or a 2×1 vector if `t` is a scalar.
"""

function gamy(d::abstractdomain, t::AbstractArray, k::Int)
    Z = similar(t)
    gamy!(Z, d, t, k)
    return Z
end

function gamx(d::abstractdomain, t::AbstractArray, k::Int)
    Z = similar(t)
    gamx!(Z, d, t, k)
    return Z
end

# Patch-specific: array
function gam(d::abstractdomain, t::Vector{Float64}, k::Int)
  Z = Matrix{Float64}(undef,2,length(t))
  gam!(Z, d, t, k)
  return Z 
end

# --- dgam: derivative of γ(t) = [γx(t); γy(t)] ---

function dgamx(d::abstractdomain, t::AbstractArray, k::Int)
    Z = similar(t)
    dgamx!(Z,d,t,k)
    return Z
end

function dgamy(d::abstractdomain, t::AbstractArray, k::Int)
    Z = similar(t)
    dgamy!(Z,d,t,k)
    return Z
end

function dgam(d::abstractdomain, t::Vector{Float64}, k::Int)
  Z = Matrix{Float64}(undef, 2, length(t))
  dgam!(Z, d, t, k)
  return Z
end

# --- gamp: a perpendicular to γ′(t) (swap + negate x) ---

function gamp(d::abstractdomain, t::Float64, k::Int)::Tuple{Float64,Float64}
    return dgamy(d, t, k), -dgamx(d, t, k)
end

function gamp(d::abstractdomain, t::Vector{Float64}, k::Int)
  Z = Matrix{Float64}(undef, 2, length(t))
  gamp!(Z, d, t, k)
  return Z
end

function gamp!(out::Vector{Float64}, d::abstractdomain, t::Float64, k::Int)

  dgam!(out, d, t, k)

  out1 = out[1]
  out[1] = out[2]
  out[2] = -out1

  return nothing
end

# --- Unit normal on the right boundary of patch k ---
function nu(d::abstractdomain, t::Float64, k::Int)::Vector{Float64}
    Z = Vector{Float64}(undef,2)
    nu!(Z, d, t, k)
    return Z
end

function nu(d::abstractdomain, t::Vector{Float64}, k::Int) 
  Z = Matrix{Float64}(undef,2,length(t))
  nu!(Z, d, t, k)
  return Z
end

function DLP(d::abstractdomain, t::Float64, l::Int, tau::AbstractArray, k::Int)
    Z = similar(tau)
    x = Vector{Float64}(undef, 2)
    G = similar(x)
    GP = similar(x)
    DLP!(Z, d, t, l, tau, k, x, G, GP)
    return Z
end

function Dmap(d::abstractdomain, u::AbstractArray{Float64}, v::AbstractArray{Float64}, k::Int)

  @assert size(u) == size(v)
  
  Z = similar(u, Float64)

  Dmap!(Z, d, u, v, k)

  return Z

end


"""
    knbh!(kdx, d, l, dnbh, loc1, loc2)

Fill `kdx` with the `dnbh`-nearest boundary neighbors of boundary patch `l`
(centered on `l`), using preallocated workspace `loc1`, `loc2` for the gamma
coordinates.

Requirements:
- `length(kdx) == 2*dnbh + 1`
- `length(loc1) == length(loc2) == 2`
"""
function knbh!(kdx::Vector{Int}, d::abstractdomain, l::Int, dnbh::Int,
  γ1::Vector{Float64}, γ2::Vector{Float64}; tol::Float64=1e-14)

    mid = dnbh + 1
    kdx[mid] = l
    left  = l
    right = l

    # grow neighbors step by step
    for n in 1:dnbh
        # ---- extend to the right (use upper endpoint of current 'right') ----
        gam!(γ1, d, 1.0, right)  # upper end of `right`, matches neighbor's -1

        found_right = 0
        for k in d.kd
            gam!(γ2, d, -1.0, k)  # lower end of candidate k
            if abs(γ1[1] - γ2[1]) <= tol && abs(γ1[2] - γ2[2]) <= tol
                found_right = k
                break
            end
        end
        right = found_right
        kdx[mid + n] = right

        # ---- extend to the left (use lower endpoint of current 'left') ----
        gam!(γ1, d, -1.0, left)   # lower end of `left`, matches neighbor's +1

        found_left = 0
        for k in d.kd
            gam!(γ2, d, 1.0, k)  # upper end of candidate k
            if abs(γ1[1] - γ2[1]) <= tol && abs(γ1[2] - γ2[2]) <= tol
                found_left = k
                break
            end
        end
        left = found_left
        kdx[mid - n] = left
    end

    return kdx
end

function knbh(d::abstractdomain, l::Int, dnbh::Int)
    kdx  = Vector{Int}(undef, 2*dnbh + 1)
    #Two temporary allocations
    γ1 = Vector{Float64}(undef, 2)
    γ2 = Vector{Float64}(undef, 2)
    knbh!(kdx, d, l, dnbh, γ1, γ2)
    return kdx
end


"""
The following two functions are used to minimise load of isnsp function,
function in boundquad and boundquadbd.
"""
# u = -1 endpoint
function mapm1(d::D, v::Float64, k::Int)::Tuple{Float64, Float64} where {D<:abstractdomain}
  return mapxy(d, -1.0, v, k)
end

# u = +1 endpoint
function mapp1(d::D, v::Float64, k::Int)::Tuple{Float64, Float64} where {D<:abstractdomain}
  return mapxy(d,  1.0, v, k)
end

"""
    isnsp(d::abstractdomain, u::Float64, v::Float64, k::Int, del::Float64)::Bool

Return true iff point `(u,v)` is within distance `del` of 
the `k`-th patch boundary. Point `(u,v)` is assumed outside
of the `k`-th patch 
"""

#Two functions below are used in domprop
@inline function dist_gam(t::Float64, d::D, u::Float64, v::Float64, k::Int)::Float64 where {D<:abstractdomain}
    γx, γy = gam(d, t, k)
    return hypot(γx - u, γy - v)
end

@inline function dist_gamBis(t::Float64, d::D, u::Float64, v::Float64, k::Int)::Float64 where {D<:abstractdomain}
    γx, γy = gam(d, t, k)
    dγx, dγy = dgam(d, t, k)

    return dγx*(γx - u) + dγy*(γy - v)
end

@inline function dist_mapm1(t::Float64, d::D, u::Float64, v::Float64, k::Int)::Float64 where {D<:abstractdomain}
  x, y = mapm1(d, t, k)
  return hypot(u - x, v - y)
end

@inline function dist_mapp1(t::Float64, d::D, u::Float64, v::Float64, k::Int)::Float64 where {D<:abstractdomain}
  x, y = mapp1(d, t, k)
  return hypot(u - x, v - y)
end

function isnsp(d::D, u::Float64, v::Float64, k::Int, del::Float64)::Bool where {D<:abstractdomain}

  # 1) quick corner checks
  x, y = mapm1(d, -1.0, k)                    # (-1,-1)
  if hypot(u - x, v - y) ≤ del
    return true
  end

  x, y = mapm1(d, 1.0, k)                    # (-1, 1)
  if hypot(u - x, v - y) ≤ del
    return true
  end

  x, y = mapp1(d, 1.0, k)                    # ( 1, 1)
  if hypot(u - x, v - y) ≤ del
    return true
  end

  x, y = mapp1(d, -1.0, k)                    # ( 1,-1)
  if hypot(u - x, v - y) ≤ del
    return true
  end

  # 2) midpoints (cheap early exits)
  x, y = mapm1(d, 0.0, k)                  # (-1,0)
  if hypot(u - x, v - y) ≤ del
    return true
  end

  x, y = mapp1(d, 0.0, k)                  # ( 1,0)
  if hypot(u - x, v - y) ≤ del
    return true
  end

  x, y = mapm1(d, -0.5, k)                  # (-1,-1/2)
  if hypot(u - x, v - y) ≤ del
    return true
  end

  x, y = mapm1(d, 0.5, k)                  # (-1, 1/2)
  if hypot(u - x, v - y) ≤ del
    return true
  end

  # 3) straight walls (v = ±1): point-to-segment distance
  Ax, Ay = mapm1(d, -1.0, k)
  Bx, By = mapp1(d, -1.0, k)   # v=-1
  if dlineseg2D(u, v, Ax, Ay, Bx, By) ≤ del
    return true
  end

  Ax, Ay = mapm1(d, 1.0, k)
  Bx, By = mapp1(d, 1.0, k)   # v=+1
  if dlineseg2D(u, v, Ax, Ay, Bx, By) ≤ del
    return true
  end

  # 4) curved walls (u = ±1): minimize distance along v ∈ [-1,1] via GSS
  tol = 1e-10

  d1, _ = GSS(dist_mapm1, -1.0, 1.0, tol, d, u, v, k)        

  if d1 ≤ del
    return true
  end

  d3, _ = GSS(dist_mapp1, -1.0, 1.0, tol, d, u, v, k)

  if d3 ≤ del
    return true
  end

  return false

end

function isnspbd(d::D, u::Float64, v::Float64, k::Int, del::Float64)::Bool where {D<:abstractdomain}

  # quick corner checks
  # gam(-1.0,k) = mapp1(d, -1.0, k)  for k in kd  
  x, y = gam(d, -1.0, k)                    # at t = -1.0 corner
  if hypot(u - x, v - y) ≤ del
    return true
  end

  # gam(1.0,k) = mapp1(d, 1.0, k)  for k in kd  
  x, y = gam(d, 1.0, k)                    # at t = 1.0 corner
  if hypot(u - x, v - y) ≤ del
    return true
  end

  # Some other checks is in the middle
  x, y = gam(d, -0.5, k)                   
  if hypot(u - x, v - y) ≤ del
    return true
  end

  x, y = gam(d, 0.0, k)                   
  if hypot(u - x, v - y) ≤ del
    return true
  end

  x, y = gam(d, 0.5, k)                   
  if hypot(u - x, v - y) ≤ del
    return true
  end

  tol = 1e-10

  d1, _ = GSS(dist_gam, -1.0, 1.0, tol, d, u, v, k)

  if d1 ≤ del
    return true
  end

  return false

end

function projbd(d::D, u::Float64, v::Float64, k::Int; tol::Float64 = 1e-14)::Float64 where {D<:abstractdomain}

  d1, t = GSS(dist_gam, -1.0, 1.0, tol, d, u, v, k)

  x, y = gam(d, 1.0, k)

  if abs( hypot(u - x, v - y) - d1 ) < tol
    return 1.0
  end

  x, y = gam(d, -1.0, k)

  if abs( hypot(u - x, v - y) - d1 ) < tol
    return -1.0
  end

  return t

end

# in your abstractdomain type: kd::Set{Int}
isbdpatch(d::abstractdomain, k::Int)::Bool = k in d.kd

"""
Draw Bounding quadrilateral for listed K patches
  use like : d = abstractdomain(b = [3,3,3,3,3]) for 
             drawboundquad(d,0.1,1:6:30) 
"""
function drawboundquad(d::abstractdomain, δ::Float64, K::Vector{Int}, flag=nothing)

  fig, ax = draw(d, flag);

  Pe = Vector{Float64}(undef,8)

  for p in K
      @views P = d.Qpts[:,p]
      extendqua!(Pe, P, δ)
      # polyline order: P1→P2→P3→P4→P1 uses (1,3,5,7) for x and (2,4,6,8) for y
      lines!(ax,P[[1, 3, 5, 7, 1]], P[[2, 4, 6, 8, 2]], color=:black, linewidth=2)
      lines!(ax,Pe[[1, 3, 5, 7, 1]], Pe[[2, 4, 6, 8, 2]], color=:gold, linewidth=2)
  end

  display(fig)
  return (fig,ax)
  #savefig(plt, "myplot.svg")
end

"""
    boundquadbd(d::abstractdomain, k::Int) -> Vector{Float64}

Bounding quadrilateral for the **boundary side** of boundary patch `k`.
`k` is expected to be in `d.kd`. Returns an 8-vector
`[x1,y1,x2,y2,x3,y3,x4,y4]` in counter-clockwise order
"""
@inline function SVum1(z1::Float64, z2::Float64)::Tuple{Bool,Float64}
    (z1 <= 0.0 && z2 <= 0.0) || return (false, 0.0)
    s = abs(z1)
    t = abs(z2)
    return true, (s >= t ? s : t)      # max(abs(z))
end

#Score Valid function for u=+1
@inline function SVup1(z1::Float64, z2::Float64)::Tuple{Bool,Float64}
    (z1 >= 1.0 && z2 >= 1.0) || return (false, 0.0)
    s = abs(z1)
    t = abs(z2)
    return true, (s >= t ? s : t)      # max(abs(z))
end

# U: valid hits with both params in [0,1]; score = 1 - min(z)
@inline function SVUm(z1::Float64, z2::Float64)::Tuple{Bool,Float64}
    (0.0 <= z1 <= 1.0 && 0.0 <= z2 <= 1.0) || return (false, 0.0)
    m = z1 <= z2 ? z1 : z2              # min(z1,z2)
    return true, 1.0 - m
end

# V: far hits with both params ≥ 1; score = max(z)
@inline function SVVm(z1::Float64, z2::Float64)::Tuple{Bool,Float64}
    (z1 >= 1.0 && z2 >= 1.0) || return (false, 0.0)
    m = z1 >= z2 ? z1 : z2              # max(z1,z2)
    return true, m
end

function boundquadbd!(P::SubArray{Float64}, d::abstractdomain, k::Int)

  # Four corners, unpack to scalars
  L1p1x, L1p1y = mapm1(d, -1.0, k)   # (u=-1, v=-1)
  L1p2x, L1p2y = mapp1(d, -1.0, k)   # (u=+1, v=-1)
  L2p1x, L2p1y = mapm1(d, 1.0, k)   # (u=-1, v=+1)
  L2p2x, L2p2y = mapp1(d, 1.0, k)   # (u=+1, v=+1)

  # boundary side u = +1

  # sample N-1 interior points: t ∈ (-1,1), exclude endpoints
  N = 101

  # U: valid hits     (0≤z≤1), score = 1 - min(z)
  best_U_score = Inf
  best_U_idx = 0

  # V: far hits       (z≥1),   score = max(z)
  best_V_score = Inf
  best_V_idx = 0

  # Loop over interior t points without allocating tpts array
  @inbounds for i in 1:(N-1)
    t = -1.0 + 2.0 * (i / N)          # same as your tpts formula

    # tau(1.0, t): point on boundary side u=+1
    Cx, Cy = mapp1(d, t, k)
    dvx, dvy = Dwall(d, 1.0, t, k)

    # z1,z2 are the intersection parameters
    z1, z2 = tanginterp(
      L1p1x, L1p1y, L1p2x, L1p2y,
      Cx, Cy, dvx, dvy,
      L2p1x, L2p1y, L2p2x, L2p2y)

    # U-set scoring
    okU, sU = SVUm(z1, z2)
    if okU && sU < best_U_score
      best_U_score = sU
      best_U_idx = i
    end

    # V-set scoring
    okV, sV = SVVm(z1, z2)
    if okV && sV < best_V_score
      best_V_score = sV
      best_V_idx = i
    end
  end

  # ----- pick best in U -----
  if best_U_idx != 0
    to = -1.0 + 2.0 * (best_U_idx / N)
    Cx, Cy = mapp1(d, to, k)
    dvx, dvy = Dwall(d, 1.0, to, k)

    P1x, P1y, P2x, P2y = tanginterx(
      L1p1x, L1p1y, L1p2x, L1p2y,
      Cx, Cy, dvx, dvy,
      L2p1x, L2p1y, L2p2x, L2p2y)
  else
    P1x, P1y = L1p2x, L1p2y
    P2x, P2y = L2p2x, L2p2y
  end

  # ----- pick best in V -----
  if best_V_idx != 0
    to = -1.0 + 2.0 * (best_V_idx / N)
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

  # assemble (CCW) as P1, P3, P4, P2
  P[1] = P1x
  P[2] = P1y
  P[3] = P3x
  P[4] = P3y
  P[5] = P4x
  P[6] = P4y
  P[7] = P2x
  P[8] = P2y

  return nothing   # or `return P` if you prefer
end

"""
Draw Bounding quadrilateral for for all the boundaries
  use like : d = abstractdomain(...) for 
             drawboundquad(d) 
"""
function drawboundquadbd(d::abstractdomain)

  fig, ax = drawbd(d);

  Pe = Vector{Float64}(undef,8)

  for ℓ in eachindex(d.kd)
    @views  P = d.Qptsbd[:, ℓ]
    extendqua!(Pe, P, 0.05)
    # polyline order: P1→P2→P3→P4→P1 uses (1,3,5,7) for x and (2,4,6,8) for y
    lines!(ax, P[[1, 3, 5, 7, 1]], P[[2, 4, 6, 8, 2]], color=:black, linewidth=2)
    lines!(ax, Pe[[1, 3, 5, 7, 1]], Pe[[2, 4, 6, 8, 2]], color=:gold, linewidth=2)
  end

  display(fig)
  return (fig, ax)
  #savefig(plt, "myplot.svg")
end