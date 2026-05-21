using Revise, Dates, FL2D

using FL2D.FLdata

function getfield_default(x, name::Symbol, default)
    return name in propertynames(x) ? getproperty(x, name) : default
end

function run_one!(io, case; s::Float64, p::Int, nJac::Int=5)

    f!, uex, _ = makediscfuex(nJac, s)

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

open("test_temp.txt", "w") do io

    s, p = 0.999, 500

    println(io, "s = ", s)
    println(io, "p = ", p)
    flush(io)

    case = (dom=FL2D.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8),
            N=12, nₚᵣ=128, δ=0.1, δ_near=0.15, δ_intp=5e-3, 
            s_small = false, benchmark = false, cond_num = false)

    run_one!(io, case; s=s, p=p, nJac=5)

    flush(io)
end

open("test.txt", "w") do io

    s, p = 0.999, 500

    println(io, "s = ", s)
    println(io, "p = ", p)
    flush(io)

    case = (dom=FL2D.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8),
            N=12, nₚᵣ=512, δ=0.1, δ_near=0.15, δ_intp=5e-3, 
            s_small = false, benchmark = false, cond_num = false)

    run_one!(io, case; s=s, p=p, nJac=5)

    flush(io)
end

# ------------------------------------------------------------
# Shared domains
# ------------------------------------------------------------
d1 = FL2D.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8);
d2 = FL2D.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8);
d3 = FL2D.disc(b=[3, 3, 3, 3, 3], L1=0.8, L2=0.8);
d5 = FL2D.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8);
d6 = FL2D.disc(b=[6, 6, 6, 6, 6], a=[3, 3, 3, 3, 5], L1=0.8, L2=0.8);

# ------------------------------------------------------------
# Direct vs Matrix Free solver:
case_specs = [
    (dom=d1, N=3,  nₚᵣ=32), (dom=d1, N=4,  nₚᵣ=32), (dom=d1, N=5,  nₚᵣ=32),
    (dom=d1, N=6,  nₚᵣ=64), (dom=d1, N=7,  nₚᵣ=64), (dom=d1, N=8,  nₚᵣ=64),

    (dom=d1, N=9,  nₚᵣ=128), (dom=d1, N=10, nₚᵣ=128), (dom=d1, N=12, nₚᵣ=128),
    (dom=d2, N=10, nₚᵣ=128), (dom=d2, N=12, nₚᵣ=128), (dom=d2, N=12, nₚᵣ=256),
    (dom=d3, N=10, nₚᵣ=128), (dom=d3, N=12, nₚᵣ=128), (dom=d3, N=12, nₚᵣ=256),
    (dom=d5, N=12, nₚᵣ=128), (dom=d5, N=12, nₚᵣ=256), (dom=d5, N=12, nₚᵣ=512), 
    (dom=d6, N=12, nₚᵣ=128), (dom=d6, N=12, nₚᵣ=256), (dom=d6, N=12, nₚᵣ=512)
];

direct_opts     = (solver=:direct, matrixfree=false, benchmark=true, cond_num=true);
matrixfree_opts = (solver=:gmres,  matrixfree=true,  benchmark=true, cond_num=false);


direct_cases      = [(; c..., direct_opts...) for c in case_specs];
matrixfree_cases  = [(; c..., matrixfree_opts...) for c in case_specs];


for (s, p, tag) in [
    (0.90, 5, "0900"),
    (0.75, 4, "0750"),
    (0.50, 4, "0500"),
    (0.25, 4, "0250")]

    open("solve_outputs$(tag)_direct_matrix.txt", "w") do io
        println(io, "Run started: ", Dates.now())
        println(io, "s = ", s)
        println(io, "p = ", p)
        flush(io)

        for case in direct_cases
            run_one!(io, case; s=s, p=p, nJac=5)
        end

        println(io, "\n==End of Disc: direct matrix==")
        println(io, "Run finished: ", Dates.now())
        flush(io)
    end

    open("solve_outputs$(tag)_matrixfree.txt", "w") do io
        println(io, "Run started: ", Dates.now())
        println(io, "s = ", s)
        println(io, "p = ", p)
        flush(io)

        for case in matrixfree_cases
            run_one!(io, case; s=s, p=p, nJac=5)
        end

        println(io, "\n==End of Disc: matrix free==")
        println(io, "Run finished: ", Dates.now())
        flush(io)
    end
end

#Small s algorithm comparison
small_true_opts = (solver=:direct, matrixfree=false, benchmark=true, cond_num=true, s_small=true)
small_false_opts = (solver=:direct, matrixfree=false, benchmark=true, cond_num=true, s_small=false)

small_true_cases  = [(; c..., small_true_opts...) for c in case_specs]
small_false_cases = [(; c..., small_false_opts...) for c in case_specs]

for (s, p, tag) in [
    (0.50, 4, "0500"),
    (0.25, 4, "0250"),
    (0.10, 4, "0100"),]

    open("solve_outputs$(tag)_s_small_false.txt", "w") do io
        println(io, "Run started: ", Dates.now())
        println(io, "s = ", s)
        println(io, "p = ", p)
        flush(io)

        for case in small_false_cases
            run_one!(io, case; s=s, p=p, nJac=5)
        end

        println(io, "\n==End of Disc: s_small=false==")
        println(io, "Run finished: ", Dates.now())
        flush(io)
    end

    open("solve_outputs$(tag)_s_small_true.txt", "w") do io
        println(io, "Run started: ", Dates.now())
        println(io, "s = ", s)
        println(io, "p = ", p)
        flush(io)

        for case in small_true_cases
            run_one!(io, case; s=s, p=p, nJac=5)
        end

        println(io, "\n==End of Disc: s_small=true==")
        println(io, "Run finished: ", Dates.now())
        flush(io)
    end
end

# =============================================
#Direct solver for many value of s

open("solve_outputs9999_direct.txt", "w") do io

    s, p = 0.9999, 5000

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    println(io, "p = ", p)
    flush(io)

    cases = [
        (dom=d1, N=10, nₚᵣ=128),
        (dom=d1, N=12, nₚᵣ=128),

        (dom=d2, N=12, nₚᵣ=128),
        (dom=d2, N=12, nₚᵣ=256),

        (dom=d3, N=12, nₚᵣ=128),
        (dom=d3, N=12, nₚᵣ=256),

        (dom=d5, N=12, nₚᵣ=128),
        (dom=d5, N=12, nₚᵣ=256),
        (dom=d5, N=12, nₚᵣ=512),
        (dom=d5, N=12, nₚᵣ=1024),

        (dom=d6, N=12, nₚᵣ=128),
        (dom=d6, N=12, nₚᵣ=256),
        (dom=d6, N=12, nₚᵣ=512),
        (dom=d6, N=12, nₚᵣ=1024),
    ]

    for case in cases
        run_one!(io, case; s=s, p=p, nJac=5)
    end

    println(io, "\n==End of Disc==")
    println(io, "Run finished: ", Dates.now())
    flush(io)
end

open("solve_outputs0999_direct.txt", "w") do io

    s, p = 0.999, 500

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    println(io, "p = ", p)
    flush(io)

    cases = [
        (dom=d1, N=10, nₚᵣ=128),
        (dom=d1, N=12, nₚᵣ=128),

        (dom=d2, N=12, nₚᵣ=128),
        (dom=d2, N=12, nₚᵣ=256),

        (dom=d3, N=12, nₚᵣ=128),
        (dom=d3, N=12, nₚᵣ=256),

        (dom=d5, N=12, nₚᵣ=128),
        (dom=d5, N=12, nₚᵣ=256),
        (dom=d5, N=12, nₚᵣ=512),
        (dom=d5, N=12, nₚᵣ=1024),

        (dom=d6, N=12, nₚᵣ=128),
        (dom=d6, N=12, nₚᵣ=256),
        (dom=d6, N=12, nₚᵣ=512),
        (dom=d6, N=12, nₚᵣ=1024),
    ]

    for case in cases
        run_one!(io, case; s=s, p=p, nJac=5)
    end

    println(io, "\n==End of Disc==")
    println(io, "Run finished: ", Dates.now())
    flush(io)
end

open("solve_outputs0990_direct.txt", "w") do io

    s, p = 0.99, 50

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    println(io, "p = ", p)
    flush(io)

    cases = [
        (dom=d1, N=10, nₚᵣ=128),
        (dom=d1, N=12, nₚᵣ=128),

        (dom=d2, N=12, nₚᵣ=128),
        (dom=d2, N=12, nₚᵣ=256),

        (dom=d3, N=12, nₚᵣ=128),
        (dom=d3, N=12, nₚᵣ=256),

        (dom=d5, N=12, nₚᵣ=128),
        (dom=d5, N=12, nₚᵣ=256),
        (dom=d5, N=12, nₚᵣ=512),

        (dom=d6, N=12, nₚᵣ=128),
        (dom=d6, N=12, nₚᵣ=256),
        (dom=d6, N=12, nₚᵣ=512),
    ]

    for case in cases
        run_one!(io, case; s=s, p=p, nJac=5)
    end

    println(io, "\n==End of Disc==")
    println(io, "Run finished: ", Dates.now())
    flush(io)
end

open("solve_outputs0900_direct.txt", "w") do io

    s, p = 0.9, 5

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    println(io, "p = ", p)
    flush(io)

    cases = [
        (dom=d1, N=10, nₚᵣ=128),
        (dom=d1, N=12, nₚᵣ=128),

        (dom=d2, N=12, nₚᵣ=128),

        (dom=d3, N=12, nₚᵣ=128),

        (dom=d5, N=12, nₚᵣ=128),
        (dom=d5, N=12, nₚᵣ=256),
        (dom=d5, N=12, nₚᵣ=512),

        (dom=d6, N=12, nₚᵣ=256),
        (dom=d6, N=12, nₚᵣ=512),
    ]

    for case in cases
        run_one!(io, case; s=s, p=p, nJac=5)
    end

    println(io, "\n==End of Disc==")
    println(io, "Run finished: ", Dates.now())
    flush(io)
end

open("solve_outputs0750_direct.txt", "w") do io

    s, p = 0.75, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    println(io, "p = ", p)
    flush(io)

    cases = [
        (dom=d1, N=10, nₚᵣ=128),
        (dom=d1, N=12, nₚᵣ=128),

        (dom=d2, N=12, nₚᵣ=128),

        (dom=d3, N=12, nₚᵣ=128),

        (dom=d5, N=12, nₚᵣ=256),

        (dom=d6, N=12, nₚᵣ=256),
    ]

    for case in cases
        run_one!(io, case; s=s, p=p, nJac=5)
    end

    println(io, "\n==End of Disc==")
    println(io, "Run finished: ", Dates.now())
    flush(io)
end

open("solve_outputs0500_direct.txt", "w") do io

    s, p = 0.5, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    println(io, "p = ", p)
    flush(io)

    cases = [
        (dom=d1, N=10, nₚᵣ=128),
        (dom=d1, N=12, nₚᵣ=128),

        (dom=d2, N=12, nₚᵣ=128),

        (dom=d3, N=12, nₚᵣ=128),

        (dom=d5, N=12, nₚᵣ=128),
        (dom=d5, N=12, nₚᵣ=256),

        (dom=d6, N=12, nₚᵣ=128),
        (dom=d6, N=12, nₚᵣ=256),
        (dom=d6, N=12, nₚᵣ=512),
    ]

    for case in cases
        run_one!(io, case; s=s, p=p, nJac=5)
    end

    println(io, "\n==End of Disc==")
    println(io, "Run finished: ", Dates.now())
    flush(io)
end

open("solve_outputs0250_direct.txt", "w") do io

    s, p = 0.25, 4

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    println(io, "p = ", p)
    flush(io)

    cases = [
        (dom=d1, N=10, nₚᵣ=128),
        (dom=d1, N=12, nₚᵣ=128),

        (dom=d2, N=12, nₚᵣ=128),

        (dom=d3, N=12, nₚᵣ=128),

        (dom=d5, N=12, nₚᵣ=128),
        (dom=d5, N=12, nₚᵣ=256),

        (dom=d6, N=12, nₚᵣ=128),
        (dom=d6, N=12, nₚᵣ=256),
        (dom=d6, N=12, nₚᵣ=512),
    ]

    for case in cases
        run_one!(io, case; s=s, p=p, nJac=5)
    end

    println(io, "\n==End of Disc==")
    println(io, "Run finished: ", Dates.now())
    flush(io)
end

open("solve_outputs0100_direct.txt", "w") do io

    s, p = 0.1, 10

    println(io, "Run started: ", Dates.now())
    println(io, "s = ", s)
    println(io, "p = ", p)
    flush(io)

    cases = [
        (dom=d1, N=10, nₚᵣ=128),
        (dom=d1, N=12, nₚᵣ=128),

        (dom=d2, N=12, nₚᵣ=128),

        (dom=d3, N=12, nₚᵣ=128),

        (dom=d5, N=12, nₚᵣ=128),
        (dom=d5, N=12, nₚᵣ=256),

        (dom=d6, N=12, nₚᵣ=128),
        (dom=d6, N=12, nₚᵣ=256),
        (dom=d6, N=12, nₚᵣ=512),
    ]

    for case in cases
        run_one!(io, case; s=s, p=p, nJac=5)
    end

    println(io, "\n==End of Disc==")
    println(io, "Run finished: ", Dates.now())
    flush(io)
end

#From this point on, s_small is true
case_specs = [
    (dom=d1, N=3,  nₚᵣ=32), (dom=d1, N=4,  nₚᵣ=32), (dom=d1, N=5,  nₚᵣ=32),
    (dom=d1, N=6,  nₚᵣ=64), (dom=d1, N=7,  nₚᵣ=64), (dom=d1, N=8,  nₚᵣ=64),

    (dom=d1, N=9,  nₚᵣ=128), (dom=d1, N=10, nₚᵣ=128), (dom=d1, N=12, nₚᵣ=128),
    (dom=d2, N=10, nₚᵣ=128), (dom=d2, N=12, nₚᵣ=128), (dom=d2, N=12, nₚᵣ=256),
    (dom=d3, N=10, nₚᵣ=128), (dom=d3, N=12, nₚᵣ=128), (dom=d3, N=12, nₚᵣ=256),
    (dom=d5, N=12, nₚᵣ=128), (dom=d5, N=12, nₚᵣ=256),
];


#Small s algorithm comparison
small_true_opts = (solver=:direct, matrixfree=false, benchmark=false, cond_num=true, s_small=true);
small_true_cases  = [(; c..., small_true_opts...) for c in case_specs];

for (s, p, tag) in [
    (0.50, 4, "0500"),
    (0.25, 4, "0250"),
    (0.10, 4, "0100"),
    (0.01, 4, "0010"),
    (1e-3, 4, "1e-3"),
    (1e-4, 4, "1e-4"),
    (1e-7, 4, "1e-7"),
    (1e-16, 4, "1e-16"),]

    open("solve_outputs$(tag)_s_small_true.txt", "w") do io
        println(io, "Run started: ", Dates.now())
        println(io, "s = ", s)
        println(io, "p = ", p)
        flush(io)

        for case in small_true_cases
            run_one!(io, case; s=s, p=p, nJac=5)
        end

        println(io, "\n==End of Disc: s_small=true==")
        println(io, "Run finished: ", Dates.now())
        flush(io)
    end
end

#flush(io) forces Julia to write any buffered output to the actual file immediately.