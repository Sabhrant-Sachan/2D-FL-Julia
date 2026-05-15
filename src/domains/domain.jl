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

abstract type abstractdomain end

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
function gamp(d::abstractdomain, t::Vector{Float64}, k::Int)
  Z = Matrix{Float64}(undef, 2, length(t))
  gamp!(Z, d, t, k)
  return Z
end

function gamp!(out::Vector{Float64}, d::abstractdomain, t::Float64, k::Int)
  @inbounds begin
    gp1, gp2 = gamp(d, t, k)

    out[1] = gp1
    out[2] = gp2
  end

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

function nu!(out::Vector{Float64}, d::abstractdomain, t::Float64, k::Int)
    @inbounds begin
        gp1, gp2 = gamp(d, t, k)

        S = hypot(gp1, gp2)

        out[1] = gp1 / S
        out[2] = gp2 / S
    end

    return nothing
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

function draw_geom(d::abstractdomain, colors; flag=nothing, L::Int=33, show::Bool=true)

    t = collect(range(-1.0, 1.0, length=L))

    tx = similar(t)
    ty = similar(t)
    xedge = similar(t)
    yedge = similar(t)

    fig = Figure(size=(650, 650))
    ax = Axis(fig[1, 1];
        aspect=DataAspect(),
        xlabel=latexstring("\$x\$"),
        ylabel=latexstring("\$y\$"),
        xlabelsize=22,
        ylabelsize=22,
        xticklabelsize=18,
        yticklabelsize=18)

    label_patches = !isnothing(flag)

    xs = Float64[]
    ys = Float64[]

    nreg = length(colors)

    @inbounds for reg in 1:nreg

        empty!(xs)
        empty!(ys)

        npreg = 0
        for p in 1:d.Npat
            npreg += d.pths[p].reg == reg
        end

        needed = npreg * 4 * (L + 1)
        sizehint!(xs, needed)
        sizehint!(ys, needed)

        for p in 1:d.Npat

            d.pths[p].reg == reg || continue

            # left edge
            fill!(tx, -1.0)
            copyto!(ty, t)
            mapxy!(xedge, yedge, d, tx, ty, p)
            append!(xs, xedge); append!(ys, yedge)
            push!(xs, NaN); push!(ys, NaN)

            # right edge
            fill!(tx, 1.0)
            copyto!(ty, t)
            mapxy!(xedge, yedge, d, tx, ty, p)
            append!(xs, xedge); append!(ys, yedge)
            push!(xs, NaN); push!(ys, NaN)

            # bottom edge
            copyto!(tx, t)
            fill!(ty, -1.0)
            mapxy!(xedge, yedge, d, tx, ty, p)
            append!(xs, xedge); append!(ys, yedge)
            push!(xs, NaN); push!(ys, NaN)

            # top edge
            copyto!(tx, t)
            fill!(ty, 1.0)
            mapxy!(xedge, yedge, d, tx, ty, p)
            append!(xs, xedge); append!(ys, yedge)
            push!(xs, NaN); push!(ys, NaN)

            if label_patches
                cx, cy = mapxy(d, 0.0, 0.0, p)
                text!(ax, cx, cy;
                    text=latexstring("\$$(p)\$"),
                    align=(:center, :center),
                    fontsize=16,
                    color=:black)
            end
        end

        if !isempty(xs)
            lines!(ax, xs, ys; color=colors[reg], linewidth=2)
        end
    end

    show && display(GLMakie.Screen(), fig)

    return fig, ax
end

function drawbd_geom(d::abstractdomain, colors; flag = true, L::Int = 33, show::Bool = true)

    t = collect(range(-1.0, 1.0, length = L))

    Tx = similar(t)
    Ty = similar(t)

    fig = Figure(size = (650, 650))

    ax = Axis(fig[1, 1];
        aspect = DataAspect(),
        xlabel = latexstring("\$x\$"),
        ylabel = latexstring("\$y\$"),
        xlabelsize = 22,
        ylabelsize = 22,
        xticklabelsize = 18,
        yticklabelsize = 18)

    label_patches = !isnothing(flag) && flag != 0

    # Group boundary curves by region/color.
    nreg = length(colors)

    xs = Float64[]
    ys = Float64[]

    @inbounds for reg in 1:nreg

        empty!(xs)
        empty!(ys)

        # Count boundary patches in this region for sizehint.
        npreg = 0
        for p in d.kd
            npreg += d.pths[p].reg == reg
        end

        needed = npreg * (L + 1)
        sizehint!(xs, needed)
        sizehint!(ys, needed)

        for p in d.kd
            d.pths[p].reg == reg || continue

            gamx!(Tx, d, t, p)
            gamy!(Ty, d, t, p)

            append!(xs, Tx)
            append!(ys, Ty)

            # Prevent Makie from connecting distinct boundary panels.
            push!(xs, NaN)
            push!(ys, NaN)

            if label_patches
                cx, cy = mapxy(d, 0.0, 0.0, p)
                text!(ax, cx, cy;
                    text = latexstring("\$$(p)\$"),
                    align = (:center, :center),
                    fontsize = 16,
                    color = :black)
            end
        end

        if !isempty(xs)
            lines!(ax, xs, ys;
                color = colors[reg],
                linewidth = 2)
        end
    end

    show && display(GLMakie.Screen(), fig)

    return fig, ax
end

function chkmap_geom(d::D, Iex::Float64; n::Int = 32, tol::Float64 = 5e-14) where {D<:abstractdomain}

    # Test function f(x,y) = x^2 + y^2
    f!(Fv, x, y) = @. Fv = x^2 + y^2

    # Chebyshev nodes and Fejér weights
    z = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        z[i] = cospi((2*i - 1) / (2n))
    end

    fw = Subroutines.getF1W(n)

    # Tensor grid
    zx = repeat(z', n, 1)   # n × n
    zy = repeat(z, 1, n)    # n × n

    # Reusable buffers
    Zx = similar(zx)
    Zy = similar(zy)
    DJ = similar(zx)
    Fv = similar(zx)

    I = 0.0

    @inbounds for k in 1:d.Npat
        mapxy_Dmap!(Zx, Zy, DJ, d, zx, zy, k)

        f!(Fv, Zx, Zy)

        @. Fv = Fv * DJ

        I += dot(fw, Fv, fw)
    end

    err = abs(Iex - I)

    if err < tol
        println("Okay!")
    else
        println("Bug!")
    end

    println("Integral approx: \n", I)
    println("Integral exact : \n", Iex)
    println("Difference     : \n", err)

    return err
end

"""
    GSS(f, a, b, tol) -> fmin, x

Golden Section Search to **minimize** a unimodal function `f` on `[a,b]`.

Returns:
- `fmin` — minimum value (≈ `f(x)`)
- `x`    — argmin in `[a,b]` (accuracy ~ `tol`)
- `iter` — iterations taken

Notes:
- To maximize `g`, call `GSS(x -> -g(x), a, b, tol=...)`.
- Assumes a single interior minimum on `[a,b]`.
- args... means “Collect all remaining positional arguments 
into a tuple named args"". 
Example:

function foo(a, b, args...)
    @show a b args
end
foo(1, 2, 3, 4, 5)

Output: 

a = 1
b = 2
args = (3, 4, 5)
"""
function GSS(f::F, a::Float64, b::Float64, tol::Float64, 
    d::D, u::Float64, v::Float64, k::Int)::Tuple{Float64,Float64} where {F,D<:abstractdomain}

    τ = (sqrt(5) - 1) / 2
    ρ = 1 - τ
    h = b - a

    x1 = a + ρ * h
    x2 = a + τ * h
    f1 = f(x1, d, u, v, k)
    f2 = f(x2, d, u, v, k)

    L = τ
    iter = 0
    while (b - a) > tol && iter < 128
        if f1 > f2
            a = x1
            x1 = x2
            f1 = f2
            x2 = a + τ * L * h
            f2 = f(x2, d, u, v, k)
        else
            b = x2
            x2 = x1
            f2 = f1
            x1 = a + ρ * L * h
            f1 = f(x1, d, u, v, k)
        end
        L *= τ
        iter += 1
    end

    return (f1 <= f2) ? (f1, x1) : (f2, x2)
end


"""
    Bis(f, a, b, maxi; tol=1e-15) -> Float64

Bisection method on interval `[a,b]` starting from midpoint 

Arguments:
- `f`     :: Function   — scalar function
- `a,b`   :: Float64       — bracketing interval
- `maxi`  :: Int    — maximum iterations
- `tol`   :: Float64       — tolerance on `abs(f(c))` (default `1e-15`)

Returns: approximation of a root in `[a,b]`.
"""
function Bis(f::F, a::Float64, b::Float64, maxi::Int, 
    d::D, u::Float64, v::Float64, r::Int; tol::Float64 = 1e-15) where {F,D<:abstractdomain}

    left  = a
    right = b
    fL    = f(a, d, u, v, r)
    fR    = f(b, d, u, v, r)

    for _ in 1:maxi
        mid = (left + right) / 2
        fM  = f(mid, d, u, v, r)

        if abs(fM) < tol
            return mid
        end

        if fL * fM < 0
            right = mid
            fR    = fM
        else
            left = mid
            fL   = fM
        end
    end

    return (left + right) / 2
end

function Bis(f::F,a::Float64, b::Float64, c::Float64, maxi::Int,
    d::D,u::Float64,v::Float64,k::Int;tol::Float64 = 1e-15) where {F,D<:abstractdomain}

    iter = 0
    fc = f(c, d, u, v, k)

    while abs(fc) > tol && iter < maxi

        fa = f(a, d, u, v, k)

        if fc * fa < 0
            b = c
        else
            a = c
        end

        c  = (a + b) / 2
        fc = f(c, d, u, v, k)

        iter += 1
    end

    return c
end

"""
    NewtonR2D(f1, f2, Jfinv, x0, y0; maxi, tol=1e-14)

Solve the 2D nonlinear system `f(x,y) = (0,0)` by Newton's method,
where `f(x,y) = (f1(x,y), f2(x,y))` and `Jfinv(x,y)` returns the **inverse**
Jacobian at `(x,y)` as a 2×2 matrix.

Returns a tuple `(x, y)` if converged within `maxi` iterations.
If the maximum iteration count is reached, returns `(:max, :max)`.

# Arguments
- `f1, f2  :: Function` — scalar-valued functions of two Float64s
- `Jfinv   :: Function` — returns a 4 tuple evaluated at `(x,y)`
- `t0, s0  :: Float64`     — initial guess
- `maxi    :: Int`      — max iterations
- `tol     :: Float64`     — convergence tolerance on component residuals 
                          (default 1e-14)

# Notes: Convergence test : it checks `|f1(x,y)|`and `|f2(x,y)|` are less
than the given tolerance.
"""
function newtonR2D(f1::F1, f2::F2, Jinv::FJ, t0::Float64, s0::Float64, maxi::Int,
                   d::D, u::Float64, v::Float64, r::Int; tol=1e-14) where {F1,F2,FJ,D<:abstractdomain}
    t = t0
    s = s0
    iter = 0
    while iter < maxi
        F1val = f1(t, s, d, u, v, r)
        F2val = f2(t, s, d, u, v, r)
        J11, J12, J21, J22 = Jinv(t, s, d, u, v, r) 

        # solve J * Δ = -F
        detJ = J11*J22 - J12*J21
        Δt   = (-F1val*J22 + F2val*J12) / detJ
        Δs   = (-F2val*J11 + F1val*J21) / detJ

        if abs(detJ) < 1e-15 
            return :max, :max  
        end
        
        t += Δt
        s += Δs

        if max(abs(Δt), abs(Δs)) < tol
            return t, s
        end
        iter += 1
    end

    return :max, :max  
end
