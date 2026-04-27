
using Revise, Dates, FractionalLaplace2D

using FractionalLaplace2D.FLdata

dobenchmark, docondnum = true, true

open("solve_outputs075.txt", "w") do io
    # =============================================
    s, p = 0.75, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)

    δ, δclsbd = 0.1, 0.01

    dom = FractionalLaplace2D.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [32, 32, 32, 32, 64, 64, 64, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(2, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=false)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        #------------- Matrix free --------

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=true)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FractionalLaplace2D.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    Anₚᵣ = [128, 128]

    AN = [10, 12]

    for i in 1:2

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=false)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        #------------- Matrix free --------

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=true)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FractionalLaplace2D.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    Anₚᵣ = [128, 128]

    AN = [10, 12]

    for i in 1:2

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=false)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        #------------- Matrix free --------

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=true)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FractionalLaplace2D.disc(b=[3, 3, 3, 3, 4], a=[2, 2, 2, 2, 3], L1=0.8, L2=0.8)

    Anₚᵣ = [128, 256]

    AN = [10, 12]

    for i in 1:2

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=false)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        #------------- Matrix free --------

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=true)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

    flush(io)

end

# ================== 2nd file ==================

open("solve_outputs025.txt", "w") do io
    
    # =============================================
    s, p = 0.25, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)

    δ, δclsbd = 0.1, 0.01

    dom = FractionalLaplace2D.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [32, 32, 32, 32, 64, 64, 64, 128, 128]

    AN = [3, 4, 5, 6, 7, 8, 9, 10, 12]

    f!, uex, fv = makediscfuex(2, s)

    for i in 1:9

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=false)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        #------------- Matrix free --------

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=true)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FractionalLaplace2D.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    Anₚᵣ = [128, 128]

    AN = [10, 12]

    for i in 1:2

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=false)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        #------------- Matrix free --------

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=true)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FractionalLaplace2D.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    Anₚᵣ = [128, 128]

    AN = [10, 12]

    for i in 1:2

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=false)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        #------------- Matrix free --------

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=true)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FractionalLaplace2D.disc(b=[3, 3, 3, 3, 4], a=[2, 2, 2, 2, 3], L1=0.8, L2=0.8)

    Anₚᵣ = [128, 256]

    AN = [10, 12]

    for i in 1:2

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=false)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        #------------- Matrix free --------

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=true)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

    flush(io)

end

#flush(io) forces Julia to write any buffered output to the actual file immediately.