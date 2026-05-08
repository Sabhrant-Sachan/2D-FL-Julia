using Revise, Dates, FL2D

using FL2D.FLdata

dobenchmark, docondnum = true, true

# Direct vs Matrix Free solver:

# ================== 1st file ==================
open("solve_outputs075.txt", "w") do io
    # =============================================
    s, p = 0.75, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)

    δ, δclsbd = 0.1, 0.01

    dom = FL2D.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

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
    dom = FL2D.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

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
    dom = FL2D.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

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
    dom = FL2D.disc(b=[3, 3, 3, 3, 4], a=[2, 2, 2, 2, 3], L1=0.8, L2=0.8)

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
open("solve_outputs05.txt", "w") do io
    
    # =============================================
    s, p = 0.5, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)

    δ, δclsbd = 0.1, 0.01

    dom = FL2D.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

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
    dom = FL2D.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

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
    dom = FL2D.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

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
    dom = FL2D.disc(b=[3, 3, 3, 3, 4], a=[2, 2, 2, 2, 3], L1=0.8, L2=0.8)

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

# ================== 3rd file ==================
open("solve_outputs025.txt", "w") do io
    
    # =============================================
    s, p = 0.25, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)

    δ, δclsbd = 0.1, 0.01

    dom = FL2D.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

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
    dom = FL2D.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    Anₚᵣ = [128, 128]

    AN = [10, 12]

    for i in 1:2

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=false)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

        #------------- Matrix free --------

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=true)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    Anₚᵣ = [128, 128]

    AN = [10, 12]

    for i in 1:2

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=false)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

        #------------- Matrix free --------

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=true)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D.disc(b=[3, 3, 3, 3, 4], a=[2, 2, 2, 2, 3], L1=0.8, L2=0.8)

    Anₚᵣ = [128, 256]

    AN = [10, 12]

    for i in 1:2

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=false)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

        #------------- Matrix free --------

        opts = Options(; plot=false, solver=:direct, cond_num=docondnum, benchmark=dobenchmark, matrixfree=true)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    println(io, "==End of Disc==")

    flush(io)

end

# =============================================
# =============================================
#Direct solver for many value of s

open("solve_outputs9999_direct.txt", "w") do io
    # =============================================
    s, p = 0.9999, 5000

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)

    δ, δclsbd = 0.1, 0.01

    dom = FL2D.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [128, 128]

    AN = [10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:2

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    nₚᵣ, N = 256, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    nₚᵣ, N = 256, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    nₚᵣ, N = 512, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    nₚᵣ, N = 512, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    nₚᵣ, N = 1024, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    println(io, "==End of Disc==")

end

open("solve_outputs0999_direct.txt", "w") do io
    # =============================================
    s, p = 0.999, 500

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)

    δ, δclsbd = 0.1, 0.01

    dom = FL2D.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [128, 128]

    AN = [10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:2

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    nₚᵣ, N = 256, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    nₚᵣ, N = 256, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    nₚᵣ, N = 512, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    nₚᵣ, N = 512, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    println(io, "==End of Disc==")

end

open("solve_outputs0990_direct.txt", "w") do io
    # =============================================
    s, p = 0.99, 50

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)

    δ, δclsbd = 0.1, 0.01

    dom = FL2D.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [128, 128]

    AN = [10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:2

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    nₚᵣ, N = 256, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    nₚᵣ, N = 256, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    nₚᵣ, N = 512, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    nₚᵣ, N = 512, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    println(io, "==End of Disc==")

end

open("solve_outputs0900_direct.txt", "w") do io
    # =============================================
    s, p = 0.9, 5

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)

    δ, δclsbd = 0.1, 0.01

    dom = FL2D.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    AN = [10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:2

        nₚᵣ, N = 128, AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    nₚᵣ, N = 128, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    nₚᵣ, N = 128, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    nₚᵣ, N = 256, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    nₚᵣ, N = 512, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    println(io, "==End of Disc==")

end

open("solve_outputs0750_direct.txt", "w") do io
    # =============================================
    s, p = 0.75, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)

    δ, δclsbd = 0.1, 0.01

    dom = FL2D.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [128, 128]

    AN = [10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:2

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    nₚᵣ, N = 128, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    nₚᵣ, N = 128, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    nₚᵣ, N = 256, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    nₚᵣ, N = 256, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    println(io, "==End of Disc==")

    # =============================================
    # dom = FL2D.disc(b=[6, 6, 6, 6, 6], a=[4, 4, 4, 4, 5], L1=0.8, L2=0.8)
    # nₚᵣ, N = 256, 12 ; dₙₕ = 1 :--> this yields 2.30e-12 Max rel error

end

open("solve_outputs0500_direct.txt", "w") do io
    # =============================================
    s, p = 0.5, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)

    δ, δclsbd = 0.1, 0.01

    dom = FL2D.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

    Anₚᵣ = [128, 128]

    AN = [10, 12]

    f!, uex, fv = makediscfuex(5, s)

    for i in 1:2

        nₚᵣ, N = Anₚᵣ[i], AN[i]

        prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
            dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

        opts = Options(; plot=false, solver=:direct)

        core_res = solveFL(prob; opts=opts)

        println(io, SolveView(prob, opts, core_res))

        flush(io)

    end

    # =============================================
    dom = FL2D.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)

    nₚᵣ, N = 128, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8)

    nₚᵣ, N = 128, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

    nₚᵣ, N = 256, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=1, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    nₚᵣ, N = 256, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    # =============================================
    dom = FL2D.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8)

    nₚᵣ, N = 512, 12

    prob = Problem(; N=N, nₚᵣ=nₚᵣ, s=s, p=p, δ=δ, δclsbd=δclsbd,
        dₙₕ=2, (f!)=f!, uex=uex, dom=dom)

    opts = Options(; plot=false, solver=:direct)

    core_res = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core_res))

    flush(io)

    println(io, "==End of Disc==")

end

#For small s, refer to other solver

#flush(io) forces Julia to write any buffered output to the actual file immediately.