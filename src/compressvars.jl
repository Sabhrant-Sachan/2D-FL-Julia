function compress_vars(d::abstractdomain, N::Int, s::Float64,
    p::Int, dₙₕ::Int, δclsbd::Float64)
    # ----------------- BASIC CONSTANTS -------------------
    # Number of points per patch
    Np = N * N

    # Constant Cs
    Cs = -sinpi(s) * (gamma(s) / (π * 2^(1 - s)))^2

    # ----------------- REGULAR INTEGRALS -------------------
    # Define nr, number of nodes for regular integral, with weights
    nr = 32
    fwr = getF1W(nr)

    # Number of regular nodes per patch
    nrp = nr * nr

    # Store regular points
    zr = Vector{Float64}(undef, nr)
    zx = Matrix{Float64}(undef, nr, nr)
    zt = Matrix{Float64}(undef, nr, nr)
    zy = Matrix{Float64}(undef, nr, nr)

    @inbounds for i in 1:nr
        ph = π * (2 * i - 1) / (2 * nr)
        cs = cos(ph)
        ss = 2.0 * sin(ph/2.0)^2 #1-cs
        zr[i] = cs
        @inbounds for j in 1:nr
            # "meshgrid(zr)" :-> zx, zy of size (nr, nr)
            zy[j, i] = cs   # rows are zr
            zx[i, j] = cs   # columns are zr
            zt[i, j] = ss   # 1-zx
        end
    end

    # Store idct values
    t₂ = Vector{Float64}(undef, nr)
    idctrg = Matrix{Float64}(undef, nr, N)

    @inline function gs(t::Float64, N::Int)::Float64
        th = 0.5 * t                 # t/2
        return sin(N * th) * cos((N - 1) * th) / sin(th)
    end

    @inbounds for i in 1:nr
        t₂[i] = π * (2 * i - 1) / (2 * nr)
    end

    @inbounds for j in 1:N
        t₁ = π * (2 * j - 1) / (2 * N)
        @inbounds for i in 1:nr
            u = t₁ + t₂[i]
            v = t₁ - t₂[i]

            idctrg[i, j] = (gs(u, N) + gs(v, N) - 1.0)/ N
        end
    end

    # ----------------- BOUNDARY PATCH INTEGRALS ------------
    # For boundary patches, which have singularity on boundary
    nbd = 128
    fwbd = getF1W(nbd)

    # z1 = cospi((2*(1:nbd)-1)/(4*nbd)).^2
    z1 = Vector{Float64}(undef, nbd)
    # We need yvec = 1 - 2*w(-z1).
    wz1 = similar(z1)             
    yvec= similar(z1)
    # [zy2, zx2] = meshgrid(zr, 1 - 2*w(-z1));
    # In MATLAB [X,Y] = meshgrid(x,y) implies
    #  X: rows = x, Y: cols = y
    zy2 = Matrix{Float64}(undef, nbd, nr)
    zx2 = Matrix{Float64}(undef, nbd, nr)
    Tz1 = Matrix{Float64}(undef, nbd, N)
    # Tz: we want nr×N like ChebyTN(N,zr)
    # Compute nr×N and then transpose once.
    Tz = Matrix{Float64}(undef, nr, N)
    gam = Vector{Float64}(undef, N)
    coeffs = Matrix{Float64}(undef, N, N)

    @inbounds for i in 1:nbd
        z1[i] = cos(π * (2 * i - 1) / (4 * nbd))^2
    end

    wfunc!(wz1, p, z1; α=-1.0, β=1.0, γ=0.0)

    @. yvec = 1 - 2 * wz1

    @inbounds for j in 1:nr
        @inbounds for i in 1:nbd
            zy2[i, j] = zr[j]     # each row of zy2 is zr
            zx2[i, j] = yvec[i]   # each column of zx2 is yvec
        end
    end

    ChebyTN!(Tz1, N, yvec)          
    ChebyTN!(Tz, N, zr)         

    gam[1] = 1.0
    @inbounds for k in 2:N
        gam[k] = 2.0
    end

    @inbounds for j in 1:N
        c = π * (2j - 1) / (2N)
        @inbounds for k in 1:N
            coeffs[k, j] = gam[k] * cos(c * (k - 1)) / N
        end
    end

    # ----------------- BOUNDARY INTEGRAL -------------------
    # For the boundary integral on the boundary curves
    N_intbd = 128
    N_nearbd = 768      # N_nearbd = 768;

    N₁ = N_intbd
    N₂ = N_nearbd

    fw₁ = getF1W(N₁)
    fw₂ = getF1W(N₂)

    # Compute boundary parametrization on N₂ points
    y₁ = Vector{Float64}(undef, N₁)
    y₂ = Vector{Float64}(undef, N₂)
    t₁ = Vector{Float64}(undef, N₁)
    t₂ = Vector{Float64}(undef, N₂)
    idctbd₁ = Matrix{Float64}(undef, N₁, N)
    idctbd₂ = Matrix{Float64}(undef, N₂, N)

    Mbd = length(d.kd)

    @inbounds for j in 1:N₁
        t₁[j] = π * (2 * j - 1) / (2 * N₁)
        y₁[j] = cos(t₁[j])
    end

    @inbounds for j in 1:N₂
        t₂[j] = π * (2 * j - 1) / (2 * N₂)
        y₂[j] = cos(t₂[j])
    end

    @inbounds for j in 1:N
        tp = π * (2j - 1) / (2N)
        @inbounds for i in 1:N₁
            u = tp + t₁[i]
            v = tp - t₁[i]
            idctbd₁[i, j] = (gs(u, N) + gs(v, N) - 1.0) / N
        end
        @inbounds for i in 1:N₂
            u = tp + t₂[i]
            v = tp - t₂[i]
            idctbd₂[i, j] = (gs(u, N) + gs(v, N) - 1.0) / N
        end
    end

    # Storage for gam!(...) or gamp!(...) 
    gamk₁ = Matrix{Float64}(undef, 2, N₁)
    gamp₁ = Matrix{Float64}(undef, 2, N₁)
    gamk₂ = Matrix{Float64}(undef, 2, N₂)
    gamp₂ = Matrix{Float64}(undef, 2, N₂)
    kbd₁  = Vector{Float64}(undef, N₁)
    kbd₂  = Vector{Float64}(undef, N₂)
    γx    = Vector{Float64}(undef, 2)
    μ₀    = Vector{Float64}(undef, 2)
    γt1   = Vector{Float64}(undef, 2)
    γt2   = Vector{Float64}(undef, 2)
    
    #For interpolating 
    Lᵢₙ = 7
    xvals = zeros(Float64,Lᵢₙ)
    # δclsbd = 0.01 and δclsbd/2 = 0.005. 
    for j in 2:Lᵢₙ
        xvals[j] = δclsbd * j / 2 
    end
    #2nd value slighty bigger than δclsbd
    #comment this line if you want to match matlab
    #xvals[2] += 1e-10 #<--- I don't do this in matlab
    fvals = Vector{Float64}(undef, Lᵢₙ)
    γxbd  = Vector{Float64}(undef, 2)
    kdx   = Vector{Int}(undef, 2*dₙₕ + 1)
    Cbd = Vector{Float64}(undef, N₂)
    Tbd = Vector{Float64}(undef, N₂)

    pCbd = FFTW.plan_r2r!(Cbd, FFTW.REDFT10; flags=FFTW.ESTIMATE)

    # ------------ TEMPORARY VARS  ------------
    # ------------ Reg Integration ------------
    CN = Matrix{Float64}(undef, N, N)
    Zx = Matrix{Float64}(undef, nr, nr)
    Zy = Matrix{Float64}(undef, nr, nr)
    DJ = Matrix{Float64}(undef, nr, nr)
    Df = Matrix{Float64}(undef, nr, nr)

    Ker = Matrix{Float64}(undef, nr, nr)  # kernel at each row x
    KIr = Matrix{Float64}(undef, nr, nr)  # Ker .* Ir
    Ur = Matrix{Float64}(undef, nr, nr)  # outer product g1*g2'
    # ------------ Bd Integration ------------
    Zx₂ = Matrix{Float64}(undef, nbd, nr)
    Zy₂ = Matrix{Float64}(undef, nbd, nr)
    DJ₂ = Matrix{Float64}(undef, nbd, nr)
    mfw = Vector{Float64}(undef, nbd)

    if s>=0.5
        qw2func!(mfw, p, z1, s)

        @. mfw = fwbd * mfw

    else
        dwfunc!(mfw, p, z1)

        tmpbd = Vector{Float64}(undef, nbd)

        wfunc!(tmpbd, p, z1; α=-1.0)

        @. mfw = fwbd * (tmpbd)^s * mfw

    end

    Ker₂= Matrix{Float64}(undef, nbd, nr)  # kernel at each row x
    KIbd= Matrix{Float64}(undef, nbd, nr)  # Ker₂ .* Ubd .* DJ₂
    Ubd = Matrix{Float64}(undef, nbd, nr)  # outer product g1*g2'
    Ubd1= Vector{Float64}(undef, nbd)      # To store g1
    Ubd2= Vector{Float64}(undef, nr)  # To store g2

    CT = Matrix{Float64}(undef, N, N)

    @inbounds for j in 1:N
        # t_j = cospi((2j - 1)/(2N))  # not needed
        @inbounds for q in 0:N-1
            CT[q+1, j] = cospi(q * (2j - 1) / (2N))
        end
    end
    # ------------ PACK EVERYTHING ------------

    #using "NamedTuples" in Julia (These are immutable)
    IV1 = (N=N, Np=Np, Cs=Cs, M=d.Npat, Mbd = Mbd)

    IVr = (nr=nr, fwr=fwr, nrp=nrp, zx=zx, zt=zt, zy=zy, Df=Df, idctrg=idctrg)

    IVbdth = ( zx2=zx2, zy2=zy2, Tz1=Tz1, Tz=Tz, coeffs=coeffs)
    
    IVbd = (N₁=N₁, N₂=N₂, fw₁=fw₁, fw₂=fw₂, idctbd₁=idctbd₁, idctbd₂=idctbd₂, dₙₕ=dₙₕ)

    IVt = (Zx=Zx, Zy=Zy, DJ=DJ, Ker=Ker, KIr=KIr, Ur=Ur, CN=CN)

    IVbt1 = (Zx₂=Zx₂, Zy₂=Zy₂, DJ₂=DJ₂, Ker₂=Ker₂, KIbd=KIbd)

    IVbt2 = (Ubd=Ubd, Ubd1=Ubd1, Ubd2=Ubd2, mfw=mfw, CT=CT)

    IVbdt1= (y₁=y₁, y₂=y₂, gamk₁=gamk₁, gamp₁=gamp₁, gamk₂=gamk₂, gamp₂=gamp₂)

    IVbdt2= (kbd₁=kbd₁, kbd₂=kbd₂, γx=γx, μ₀=μ₀, γt1=γt1, γt2=γt2)

    IVbdt3= (Lᵢₙ=Lᵢₙ, xvals=xvals, fvals=fvals, γxbd=γxbd, kdx=kdx, Cbd=Cbd, pCbd=pCbd, Tbd=Tbd)
    
    IV = (IV1=IV1, IVr=IVr, IVbdth=IVbdth, IVbd=IVbd, IVt=IVt, 
    IVbt1=IVbt1, IVbt2=IVbt2, IVbdt1=IVbdt1, IVbdt2=IVbdt2, IVbdt3=IVbdt3)

    return IV
end

using FFTW

"""
compress_vars_2(d, N, s, p, dₙₕ, δclsbd)

Builds a NamedTuple IV containing only what is needed for Ax! (matvec),
with preallocated scratch arrays and FFTW r2r plans to avoid allocations.
"""
function compress_vars_2(d::abstractdomain, N::Int, s::Float64,
    p::Int, dₙₕ::Int, δclsbd::Float64)

    # ----------------- BASIC CONSTANTS -------------------
    Np  = N * N
    M   = d.Npat
    Mbd = length(d.kd)

    Cs = -sinpi(s) * (gamma(s) / (π * 2^(1 - s)))^2

    # ----------------- REGULAR INTEGRALS GRID -------------------
    nr  = 32
    nrp = nr * nr
    fwr = getF1W(nr)                    # must be pre-existing

    # regular Chebyshev nodes (nr)
    zr = Vector{Float64}(undef, nr)
    for i in 1:nr
        zr[i] = cospi((2i - 1) / (2nr))
    end

    # meshgrid(zr)
    zx = Matrix{Float64}(undef, nr, nr)
    zy = Matrix{Float64}(undef, nr, nr)
    for j in 1:nr
        cj = zr[j]
        @inbounds for i in 1:nr
            zx[i,j] = zr[i]
            zy[i,j] = cj
        end
    end

    # ----------------- BOUNDARY INTEGRAL GRIDS -------------------
    N₁ = 128
    N₂ = 768
    fw₁ = getF1W(N₁)
    fw₂ = getF1W(N₂)

    y₁ = Vector{Float64}(undef, N₁)
    y₂ = Vector{Float64}(undef, N₂)
    for j in 1:N₁
        y₁[j] = cospi((2j - 1) / (2N₁))
    end
    for j in 1:N₂
        y₂[j] = cospi((2j - 1) / (2N₂))
    end

    # ----------------- Scratch vectors for this Ax! -------------------
    Lp   = M * Np
    Ltot = Lp + Mbd * N

    chebcoef  = Vector{Float64}(undef, Ltot)
    ufin      = Vector{Float64}(undef, M * nrp)

    zeta1     = Vector{Float64}(undef, Mbd * N₁)
    zeta2     = Vector{Float64}(undef, Mbd * N₂)
    zeta2coef = Vector{Float64}(undef, Mbd * N₂)

    # ----------------- FFTW buffers + plans -------------------
    #
    # MATLAB dct == DCT-II  => FFTW.REDFT10
    # MATLAB idct == DCT-III => FFTW.REDFT01
    #
    # 2D transforms on N×N buffer:
    FV = Matrix{Float64}(undef, N, N)   # function values
    CN = Matrix{Float64}(undef, N, N)   # coeffs (in-place target)

    # Plans: apply DCT-II in-place along dim 2 then dim 1
    p_dct2_dim2 = plan_r2r!(CN, FFTW.REDFT10, 2; flags=FFTW.ESTIMATE)
    p_dct2_dim1 = plan_r2r!(CN, FFTW.REDFT10, 1; flags=FFTW.ESTIMATE)

    # For refinement: need IDCT-III from N×N coeffs to nr×nr values.
    # We’ll do it as two 1D IDCTs using intermediate buffers:
    # - first IDCT along dim2 into a nr×N buffer
    # - then IDCT along dim1 into nr×nr
    #
    # Buffers:
    TMP_nrN = Matrix{Float64}(undef, nr, N)
    U_nrnr  = Matrix{Float64}(undef, nr, nr)

    # Plans for 1D IDCT-III on length=nr and length=N etc.
    # But FFTW plans are tied to array sizes. We use r2r! on full matrices:
    p_idct_nr_dim2 = plan_r2r!(TMP_nrN, FFTW.REDFT01, 2; flags=FFTW.ESTIMATE)  # acts on each row length N (careful)
    # For simplicity/robustness, we instead use *matrix-based ChebyTN eval* for refinement in Julia later
    # OR build proper plans per dimension lengths. If you want pure FFTW, we need dedicated buffers sized (nr, nr) etc.

    # Practical approach (less finicky): do refinement using Chebyshev evaluation matrices you precompute once:
    # Build TnrN = ChebyTN(n, zr) etc. This is deterministic and allocation-free with mul!.
    # You already did this style for bdop; here we do similar:
    TnrN = Matrix{Float64}(undef, nr, N)
    ChebyTN!(TnrN, N, zr)   # must exist (fills TnrN[i,j] = T_{j-1}(zr[i]))

    # For boundary ζ upsampling (N -> N₁/N₂):
    TN1N = Matrix{Float64}(undef, N₁, N)
    TN2N = Matrix{Float64}(undef, N₂, N)
    ChebyTN!(TN1N, N, y₁)
    ChebyTN!(TN2N, N, y₂)

    # 1D DCT-II buffer for boundary coeffs:
    cN   = Vector{Float64}(undef, N)
    p_dct1_N = plan_r2r!(cN, FFTW.REDFT10; flags=FFTW.ESTIMATE)

    # --- Boundary geometry scratch (used in bd integral + regular integrals) ---
    gamk₁  = Matrix{Float64}(undef, 2, N₁)
    gamp₁  = Matrix{Float64}(undef, 2, N₁)
    gamk₂  = Matrix{Float64}(undef, 2, N₂)
    gamp₂  = Matrix{Float64}(undef, 2, N₂)

    # kernel dot helper buffers
    kbd₁ = Vector{Float64}(undef, N₁)
    kbd₂ = Vector{Float64}(undef, N₂)

    γx   = Vector{Float64}(undef, 2)
    μ₀   = Vector{Float64}(undef, 2)
    γxbd = Vector{Float64}(undef, 2)
    νvec = Vector{Float64}(undef, 2)

    # interpolation scratch
    Lᵢₙ = 7
    xvals = Vector{Float64}(undef, Lᵢₙ)
    xvals[1] = 0.0
    for i in 2:Lᵢₙ
        xvals[i] = δclsbd + (i-2) * (δclsbd/2)   # matches your MATLAB [0, del, del+0.005, ...]
    end
    fvals = Vector{Float64}(undef, Lᵢₙ)

    kdx = Vector{Int}(undef, 2dₙₕ + 1)

    # regular integral scratch
    phix = Matrix{Float64}(undef, nr, nr)
    phiy = Matrix{Float64}(undef, nr, nr)
    DJ   = Matrix{Float64}(undef, nr, nr)
    Df   = Matrix{Float64}(undef, nr, nr)
    Ker  = Matrix{Float64}(undef, nr, nr)

    # boundary-patch regularization scratch (if you keep that part in Ax! later)
    # add as needed.

    IV1 = (N=N, Np=Np, Cs=Cs, M=M, Mbd=Mbd, nr=nr, nrp=nrp, N₁=N₁, N₂=N₂, dₙₕ=dₙₕ, δclsbd=δclsbd)

    IVx = (
        chebcoef=chebcoef, ufin=ufin,
        zeta1=zeta1, zeta2=zeta2, zeta2coef=zeta2coef
    )

    IVt = (
        FV=FV, CN=CN, TMP_nrN=TMP_nrN, U_nrnr=U_nrnr,
        TnrN=TnrN, TN1N=TN1N, TN2N=TN2N,
        cN=cN, p_dct1_N=p_dct1_N,
        p_dct2_dim2=p_dct2_dim2, p_dct2_dim1=p_dct2_dim1
    )

    IVbdt = (
        y₁=y₁, y₂=y₂, fw₁=fw₁, fw₂=fw₂,
        gamk₁=gamk₁, gamp₁=gamp₁, gamk₂=gamk₂, gamp₂=gamp₂,
        kbd₁=kbd₁, kbd₂=kbd₂,
        γx=γx, μ₀=μ₀, γxbd=γxbd, νvec=νvec,
        Lᵢₙ=Lᵢₙ, xvals=xvals, fvals=fvals, kdx=kdx
    )

    IVreg = (
        fwr=fwr, zx=zx, zy=zy,
        phix=phix, phiy=phiy, DJ=DJ, Df=Df, Ker=Ker
    )

    return (IV1=IV1, IVx=IVx, IVt=IVt, IVbdt=IVbdt, IVreg=IVreg)
end

#Unpacking is also simple:
# Unpack the four structs
#(; IV1, IVr, IVbdth, IVbd) = IV
# Unpack fields from IV1
#(; N, Np, Cs, M, Mbd) = IV1
#NOTE: this is named-field unpacking
# and not positional unpacking. For 
# instance, if nt = (a = 10, b = 20, c = 30)   
# then (; a, c) = nt gives in scope a = 10 and 
#c = 30, b is not given. 
