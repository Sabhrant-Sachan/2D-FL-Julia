using Revise, Dates, FL2D

using FL2D.FLdata

function getfield_default(x, name::Symbol, default)
    return name in propertynames(x) ? getproperty(x, name) : default
end

function run_one_squircle!(io, case; s::Float64, p::Int, nJac::Int=5)

    f!, uex, _ = makesquirclefuex(nJac, s, 4)

    # Per-case controls, with defaults
    δ        = getfield_default(case, :δ, 0.1)
    δ_near   = getfield_default(case, :δ_near, 0.15)
    δ_intp   = getfield_default(case, :δ_intp, 5e-3)

    solver     = getfield_default(case, :solver, :direct)
    matrixfree = getfield_default(case, :matrixfree, false)
    s_small    = getfield_default(case, :s_small, false)
    benchmark  = getfield_default(case, :benchmark, false)
    cond_num   = getfield_default(case, :cond_num, false)

    prob = Problem(;
        N = case.N,
        nₚᵣ = case.nₚᵣ,
        s = s, p = p, δ = δ,
        δ_near = δ_near,
        δ_intp = δ_intp,
        f! = f!, uex = uex,
        dom = case.dom)

    opts = Options(;
        plot = false,
        solver = solver,
        cond_num = cond_num,
        benchmark = benchmark,
        matrixfree = matrixfree,
        s_small = s_small)

    core = solveFL(prob; opts=opts)

    println(io, SolveView(prob, opts, core))
    flush(io)

    return nothing
end

open("test.txt", "w") do io

    s, p = 0.2, 5

    println(io, "s = ", s)
    println(io, "p = ", p)
    flush(io)

    case = (dom=FL2D.squircle(b=[1, 1, 1, 1, 1], P=4.0),
            N=10, nₚᵣ=128, δ=0.1, δ_near=0.15, δ_intp=5e-3, 
            s_small = false, benchmark = false, cond_num = false)

    run_one_squircle!(io, case; s=s, p=p, nJac=5)

    flush(io)
end

