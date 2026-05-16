"""
If matrix_form=true, compress_vars builds a NamedTuple IV which is needed for 
functions Axbdop.jl, Axbdpth.jl, and Axintpth.jl (these functions help build 
the matrix A) with preallocated scratch arrays to avoid allocations.

If matrix_form=false, compress_vars builds a NamedTuple IV which is needed for 
functions Axb.jl (For a matrix free iterative solver) with preallocated scratch 
arrays to avoid allocations.
"""
function compress_vars(d::abstractdomain, dp::domprop, s::Float64,
    p::Int; matrix_form=true)
    # ----------------- BASIC CONSTANTS -------------------
    # Number of points per patch
    N = dp.N
    Np = N * N
    M = d.Npat
    Mbd = length(d.kd)
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

    # ----------------- Boundary DLP quadrature -------------------
    nr_bdy = 64
    ns_near = 128

    fwr_bdy = getF1W(nr_bdy)
    fw_near = getF1W(ns_near)

    zr_bdy = Vector{Float64}(undef, nr_bdy)

    @inbounds for j in 1:nr_bdy
        zr_bdy[j] = cospi((2 * j - 1) / (2 * nr_bdy))
    end

    z_near1 = Vector{Float64}(undef, ns_near)
    z_near2 = Vector{Float64}(undef, ns_near)

    @inbounds for i in 1:ns_near
        c = π * (2 * i - 1) / (2 * ns_near)
        z_near1[i] = cos(c / 2)^2
        z_near2[i] = sin(c / 2)^2
    end

    # Dynamic split points for interior projected parameter t0.
    y₁ = Vector{Float64}(undef, ns_near)
    y₂ = Vector{Float64}(undef, ns_near)

    # Fixed endpoint split nodes:
    # t0 = +1 uses yt1, t0 = -1 uses yt2.
    yt1 = Vector{Float64}(undef, ns_near)
    yt2 = Vector{Float64}(undef, ns_near)

    wmz₁ = Vector{Float64}(undef, ns_near)
    wmz₂ = Vector{Float64}(undef, ns_near)
    dwz₁ = Vector{Float64}(undef, ns_near)
    dwz₂ = Vector{Float64}(undef, ns_near)

    wfunc!(wmz₁, p, z_near1; α=-1.0)
    wfunc!(wmz₂, p, z_near2; α=-1.0)

    @. yt1 = 1.0 - 2.0 * wmz₁
    @. yt2 = -1.0 + 2.0 * wmz₂

    dwfunc!(dwz₁, p, z_near1)
    dwfunc!(dwz₂, p, z_near2)

    # Kernel-weight buffers.
    kbd_bdy = Vector{Float64}(undef, nr_bdy)
    kbd₁ = Vector{Float64}(undef, ns_near)
    kbd₂ = Vector{Float64}(undef, ns_near)

    # Geometry buffers for regular boundary quadrature.
    gamkbdy = Matrix{Float64}(undef, 2, nr_bdy)
    gampkbdy = Matrix{Float64}(undef, 2, nr_bdy)

    # Geometry buffers for dynamic split.
    gamks₁ = Matrix{Float64}(undef, 2, ns_near)
    gamks₂ = Matrix{Float64}(undef, 2, ns_near)
    gamperpks₁ = Matrix{Float64}(undef, 2, ns_near)
    gamperpks₂ = Matrix{Float64}(undef, 2, ns_near)

    # Geometry buffers for endpoint split.
    gamkt1 = Matrix{Float64}(undef, 2, ns_near)
    gamkt2 = Matrix{Float64}(undef, 2, ns_near)
    gamperpkt1 = Matrix{Float64}(undef, 2, ns_near)
    gamperpkt2 = Matrix{Float64}(undef, 2, ns_near)

    μ₀ = Vector{Float64}(undef, 2)
    γt1 = Vector{Float64}(undef, 2)
    γt2 = Vector{Float64}(undef, 2)

    fvals = Vector{Float64}(undef, dp.Lᵢₙ)
    γxbd = Vector{Float64}(undef, 2)
    Tbd = Vector{Float64}(undef, N)

    inv2π = 1 / (2π)

    BD_FAR = UInt8(0)
    BD_NEAR = UInt8(1)
    BD_INTP = UInt8(2)

    # Chebyshev polynomials at fixed endpoint split nodes.
    Tbdt1 = Matrix{Float64}(undef, ns_near, N)
    Tbdt2 = Matrix{Float64}(undef, ns_near, N)

    ChebyTN!(Tbdt1, N, yt1)
    ChebyTN!(Tbdt2, N, yt2)

    # Dynamic split basis values.
    Tbdt₀ = Matrix{Float64}(undef, ns_near, N)

    # ------------ TEMPORARY VARS  ------------
    # ------------ Reg Integration ------------
    Zx = Matrix{Float64}(undef, nr, nr)
    Zy = Matrix{Float64}(undef, nr, nr)
    DJ = Matrix{Float64}(undef, nr, nr)
    Df = Matrix{Float64}(undef, nr, nr)

    KIr = Matrix{Float64}(undef, nr, nr)  # Ker .* Ir

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

    KIbd= Matrix{Float64}(undef, nbd, nr)  # Ker₂ .* Ubd .* DJ₂
    CT = Matrix{Float64}(undef, N, N)

    @inbounds for j in 1:N
        # t_j = cospi((2j - 1)/(2N))  # not needed
        @inbounds for q in 0:N-1
            CT[q+1, j] = cospi(q * (2j - 1) / (2N))
        end
    end


    if matrix_form

        # ------ matrix-form-only coefficient machinery ------
        # Store idct values
        t₂ = Vector{Float64}(undef, nr)
        idctrg = Matrix{Float64}(undef, nr, N)

        @inbounds for i in 1:nr
            t₂[i] = π * (2 * i - 1) / (2 * nr)
        end

        @inline function gs(t::Float64, N::Int)::Float64
            th = 0.5 * t                 # t/2
            return sin(N * th) * cos((N - 1) * th) / sin(th)
        end

        @inbounds for j in 1:N
            t₁ = π * (2 * j - 1) / (2 * N)
            for i in 1:nr
                u = t₁ + t₂[i]
                v = t₁ - t₂[i]
                idctrg[i, j] = (gs(u, N) + gs(v, N) - 1.0) / N
            end
        end

        gamvals = Vector{Float64}(undef, N)
        coeffs = Matrix{Float64}(undef, N, N)
        coeffs_sum = Vector{Float64}(undef, N)

        gamvals[1] = 1.0
        @inbounds for k in 2:N
            gamvals[k] = 2.0
        end

        @inbounds for j in 1:N
            c = π * (2j - 1) / (2N)
            @inbounds for k in 1:N
                coeffs[k, j] = gamvals[k] * cos(c * (k - 1)) / N
            end
        end

        @inbounds for j in 1:N
            @views cc = coeffs[:, j]
            coeffs_sum[j] = sum(cc)
        end

        B1bd = Matrix{Float64}(undef, nbd, N)
        B2bd = Matrix{Float64}(undef, nr, N)

        mul!(B1bd, Tz1, coeffs)   # nbd × N
        mul!(B2bd, Tz, coeffs)   # nr × N

        Gbdtmp = Matrix{Float64}(undef, N, nr)
        Gbd = Matrix{Float64}(undef, N, N)
        Wbd = Matrix{Float64}(undef, nbd, nr)

        Grgtmp = Matrix{Float64}(undef, N, nr)
        Grg = Matrix{Float64}(undef, N, N)
        Wrg = Matrix{Float64}(undef, nr, nr)

        # Boundary basis values on regular boundary nodes.
        γx = Vector{Float64}(undef, 2)
        tbd = Vector{Float64}(undef, nr_bdy)
        idctbd = Matrix{Float64}(undef, nr_bdy, N)

        @inbounds for i in 1:nr_bdy
            tbd[i] = π * (2 * i - 1) / (2 * nr_bdy)
        end

        @inbounds for j in 1:N
            tp = π * (2 * j - 1) / (2 * N)

            @inbounds for i in 1:nr_bdy
                u = tp + tbd[i]
                v = tp - tbd[i]
                idctbd[i, j] = (gs(u, N) + gs(v, N) - 1.0) / N
            end
        end

        # Basis values at fixed endpoint split nodes.
        Bbdt1 = Matrix{Float64}(undef, ns_near, N)
        Bbdt2 = Matrix{Float64}(undef, ns_near, N)

        mul!(Bbdt1, Tbdt1, coeffs)
        mul!(Bbdt2, Tbdt2, coeffs)

        Bbdt₀ = Matrix{Float64}(undef, ns_near, N)

        Tbd = Vector{Float64}(undef, N)

        #using "NamedTuples" in Julia (These are immutable)
        # ------------ PACK EVERYTHING ------------
        # ============================================================
        # Matrix-form IV Used by:
        #   Axbdop!, Axbdpth!, Axintpth!
        # I do not pack matrix-free-only buffers or FFTW matrices here.
        # ============================================================
        IV1 = (N=N, Np=Np, Cs=Cs, M=M, Mbd=Mbd, nbd=nbd)

        IVr = (nr=nr, fwr=fwr, nrp=nrp, zx=zx, zt=zt, zy=zy, Df=Df, idctrg=idctrg)

        IVbdth = (zx2=zx2, zy2=zy2, coeffs=coeffs, coeffs_sum=coeffs_sum)

        IVt = (Zx=Zx, Zy=Zy, DJ=DJ, KIr=KIr, Wrg=Wrg, Grgtmp=Grgtmp, Grg=Grg)

        IVbt1 = (Zx₂=Zx₂, Zy₂=Zy₂, DJ₂=DJ₂, KIbd=KIbd)

        IVbt2 = (mfw=mfw, CT=CT, B1bd=B1bd, B2bd=B2bd, Gbdtmp=Gbdtmp, Gbd=Gbd, Wbd=Wbd)
        
        #Storing boundary variables in the following 4 structs
        #IVbd, IVbdt1, IVbdt2, and IVbdt3.
        IVbd = (
            nr_bdy=nr_bdy,
            ns_near=ns_near,
            fwr_bdy=fwr_bdy,
            fw_near=fw_near,
            zr_bdy=zr_bdy,
            idctbd=idctbd,
            inv2π=inv2π,
            BD_FAR=BD_FAR,
            BD_NEAR=BD_NEAR,
            BD_INTP=BD_INTP)

        IVbdt1 = (
            y₁=y₁,y₂=y₂,
            yt1=yt1,yt2=yt2,
            wmz₁=wmz₁,wmz₂=wmz₂,
            dwz₁=dwz₁,dwz₂=dwz₂,
            gamkbdy=gamkbdy,
            gampkbdy=gampkbdy,
            gamks₁=gamks₁,
            gamks₂=gamks₂,
            gamperpks₁=gamperpks₁,
            gamperpks₂=gamperpks₂,
            gamkt1=gamkt1,
            gamkt2=gamkt2,
            gamperpkt1=gamperpkt1,
            gamperpkt2=gamperpkt2)

        IVbdt2 = (
            kbd_bdy=kbd_bdy,
            kbd₁=kbd₁,
            kbd₂=kbd₂,
            γx=γx,μ₀=μ₀,
            γt1=γt1,γt2=γt2)

        IVbdt3 = (
            fvals=fvals,
            γxbd=γxbd,
            Tbd=Tbd,
            Tbdt1=Tbdt1,
            Tbdt2=Tbdt2,
            Tbdt₀=Tbdt₀,
            Bbdt1=Bbdt1,
            Bbdt2=Bbdt2,
            Bbdt₀=Bbdt₀)

        IV = (IV1=IV1, IVr=IVr, IVbdth=IVbdth, IVbd=IVbd, IVt=IVt, IVbt1=IVbt1,
            IVbt2=IVbt2, IVbdt1=IVbdt1, IVbdt2=IVbdt2, IVbdt3=IVbdt3)

    else
        # ------ matrix-free-only coefficient machinery ------

        chebcoef = Vector{Float64}(undef, M * Np + Mbd * N)
        ufin = Vector{Float64}(undef, M * nrp)

        # Boundary density values:
        # ζ_bdy     : ζ evaluated on the regular boundary Chebyshev grid zr_bdy
        # ζ_neart1  : ζ evaluated on fixed endpoint-smoothed nodes yt1
        # ζ_neart2  : ζ evaluated on fixed endpoint-smoothed nodes yt2
        #
        # Dynamic near nodes y₁/y₂ are handled inside Ax! by Chebyshev
        # interpolation using Tbdt₀ and the stored boundary coefficients.
        ζ_bdy = Vector{Float64}(undef, Mbd * nr_bdy)
        ζ_neart1 = Vector{Float64}(undef, Mbd * ns_near)
        ζ_neart2 = Vector{Float64}(undef, Mbd * ns_near)

        UV = Matrix{Float64}(undef, N, N)
        CN = Matrix{Float64}(undef, N, N)
        CNnr = Matrix{Float64}(undef, N, nr)
        UFV = Matrix{Float64}(undef, nr, nr)

        Ker = Matrix{Float64}(undef, nr, nr)
        Ur = Matrix{Float64}(undef, nr, nr)

        Ker₂ = Matrix{Float64}(undef, nbd, nr)
        Ubd = Matrix{Float64}(undef, nbd, nr)

        p_dct2_dim2 = FFTW.plan_r2r!(UV, FFTW.REDFT10, 2; flags=FFTW.MEASURE)
        p_dct2_dim1 = FFTW.plan_r2r!(UV, FFTW.REDFT10, 1; flags=FFTW.MEASURE)

        p_dct3_dim2 = FFTW.plan_r2r!(UFV, FFTW.REDFT01, 2; flags=FFTW.MEASURE)
        p_dct3_dim1 = FFTW.plan_r2r!(UFV, FFTW.REDFT01, 1; flags=FFTW.MEASURE)

        # Boundary DCT work vectors.
        ζv = Vector{Float64}(undef, N)
        ζfv_bdy = Vector{Float64}(undef, nr_bdy)

        p_dct2_N = FFTW.plan_r2r!(ζv, FFTW.REDFT10; flags=FFTW.MEASURE)
        p_dct3_bdy = FFTW.plan_r2r!(ζfv_bdy, FFTW.REDFT01; flags=FFTW.MEASURE)

        TzT = transpose(Tz)

        # Boundary Chebyshev values at one projected point tₛ.
        # Used for ζ(tₛ)/2 in the DLP jump term.
        Tbd = Vector{Float64}(undef, N)

        # ============================================================
        # Matrix-free IV used by Ax!
        # ============================================================

        IV1 = (
            N=N,
            Np=Np,
            Cs=Cs,
            M=M,
            Mbd=Mbd
        )

        IVr = (
            nr=nr,
            fwr=fwr,
            nrp=nrp,
            zx=zx,
            zt=zt,
            zy=zy,
            Df=Df
        )

        IVbdth = (
            zx2=zx2,
            zy2=zy2,
            Tz1=Tz1
        )

        IVt = (
            Zx=Zx,
            Zy=Zy,
            DJ=DJ,
            Ker=Ker,
            KIr=KIr,
            Ur=Ur,
            CN=CN
        )

        IVbt1 = (
            Zx₂=Zx₂,
            Zy₂=Zy₂,
            DJ₂=DJ₂,
            Ker₂=Ker₂,
            KIbd=KIbd
        )

        IVbt2 = (
            Ubd=Ubd,
            mfw=mfw,
            CT=CT
        )

        IVbd = (
            nr_bdy=nr_bdy,
            ns_near=ns_near,
            fwr_bdy=fwr_bdy,
            fw_near=fw_near,
            zr_bdy=zr_bdy,
            inv2π=inv2π,
            BD_FAR=BD_FAR,
            BD_NEAR=BD_NEAR,
            BD_INTP=BD_INTP
        )

        IVbdt1 = (
            y₁=y₁,
            y₂=y₂,
            yt1=yt1,
            yt2=yt2,
            wmz₁=wmz₁,
            wmz₂=wmz₂,
            dwz₁=dwz₁,
            dwz₂=dwz₂,
            gamkbdy=gamkbdy,
            gampkbdy=gampkbdy,
            gamks₁=gamks₁,
            gamks₂=gamks₂,
            gamperpks₁=gamperpks₁,
            gamperpks₂=gamperpks₂,
            gamkt1=gamkt1,
            gamkt2=gamkt2,
            gamperpkt1=gamperpkt1,
            gamperpkt2=gamperpkt2
        )

        IVbdt2 = (
            kbd_bdy=kbd_bdy,
            kbd₁=kbd₁,
            kbd₂=kbd₂,
            μ₀=μ₀,
            γt1=γt1,
            γt2=γt2
        )

        IVbdt3 = (
            fvals=fvals,
            γxbd=γxbd,
            Tbd=Tbd,
            Tbdt₀=Tbdt₀,
            Tbdt1=Tbdt1,
            Tbdt2=Tbdt2
        )

        IVAf = (
            chebcoef=chebcoef,
            ufin=ufin,
            ζ_bdy=ζ_bdy,
            ζ_neart1=ζ_neart1,
            ζ_neart2=ζ_neart2,
            UV=UV,
            UFV=UFV,
            ζv=ζv,
            ζfv_bdy=ζfv_bdy,
            CNnr=CNnr,
            TzT=TzT
        )

        IVAdct = (
            p_dct2_dim2=p_dct2_dim2,
            p_dct2_dim1=p_dct2_dim1,
            p_dct3_dim2=p_dct3_dim2,
            p_dct3_dim1=p_dct3_dim1,
            p_dct2_N=p_dct2_N,
            p_dct3_bdy=p_dct3_bdy
        )

        IV = (
            IV1=IV1,
            IVr=IVr,
            IVbdth=IVbdth,
            IVbd=IVbd,
            IVt=IVt,
            IVbt1=IVbt1,
            IVbt2=IVbt2,
            IVbdt1=IVbdt1,
            IVbdt2=IVbdt2,
            IVbdt3=IVbdt3,
            IVAf=IVAf,
            IVAdct=IVAdct
        )
    end

    return IV

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
