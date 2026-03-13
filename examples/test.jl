using GLMakie, LaTeXStrings, ColorSchemes, LinearAlgebra

using FFTW

using Subroutines, Random, IterativeSolvers

using BenchmarkTools, Profile, ProfileView

# For gamma function 
using SpecialFunctions: gamma

# brings Patch and Disc into this scope
include(joinpath(@__DIR__, "domains/domain.jl"));
include(joinpath(@__DIR__, "initialvars/compressvars.jl")); 
include(joinpath(@__DIR__, "domprop.jl"));
include(joinpath(@__DIR__, "precompsH.jl"));
include(joinpath(@__DIR__, "precompsL.jl"));
include(joinpath(@__DIR__, "bvec.jl"));
include(joinpath(@__DIR__, "matrixvecprod/Axintpth.jl"));
include(joinpath(@__DIR__, "matrixvecprod/Axbdpth.jl"));
include(joinpath(@__DIR__, "matrixvecprod/Axbdop.jl"));

print("\033c")

d = disc(b =  [3,3,3,3,3])   

d = ellipse(b =  [3,3,3,3,3]) 

@btime disc(b =  [3,3,3,3,3]);

@benchmark disc(b = [3,3,3,3,3])


dp = domprop(12,0.1,0.01,d)

@btime domprop(12,0.1,0.01,$d)

chkinvpts(dp,d)

@btime chkinvpts(dp,d)

@btime domprop(12,0.1,0.01,d)


d = disc(b =  [1,1,1,1,1])   

dp = domprop(12,0.1,0.01,d)

s = 0.9
p = 5
# Compute Precomputations 
if s >= 0.5
    IntS = precompsH(d, dp, s, p);
else
    IntS = precompsL(d, dp, s, p);
end

using MAT
f = matopen("I_ref.mat"); 
I_ref = read(f, "Int");
close(f)

err = abs.(IntS .- I_ref);
for j in 1273:48+1273
err_absj = maximum(abs.(IntS[:,j] .- I_ref[:,j]))
display(err_absj)
end

err_abs = maximum(abs.(IntS .- I_ref))
err_rel = err_abs / maximum(abs.(I_ref))

function f!(F, x, y)
    fill!(F, 1.0)
    return nothing
end

# Compute the right hand vector
b = bvec(d, dp, s, f!)

# Forming the Matrix A
# Compress <-> decompress Vars 
#Check everything below this once!
δclsbd = 0.01 
dₙₕ = 1
IV = compress_vars(d, dp.N, s, p, dₙₕ, δclsbd);

(; N, Np, M, Mbd) = IV.IV1;

Lₚₙ = M * Np
Lₚ = Lₚₙ + Mbd * N

#A = Matrix{Float64}(undef, Lₚ, Lₚ)

A = zeros(Float64, Lₚ, Lₚ);

function constuctA!(A, IntS, d, dp, s, IV)

    (; N, Np, M, Mbd) = IV.IV1;

    Lₚₙ = M * Np
    for k in 1:M
        @views v = A[:, (1+Np*(k-1)):Np*k]

        if k in d.kd
            Axbdpth!(v, k, IntS, d, dp, s, IV)
        else
            Axintpth!(v, k, IntS, d, dp, s, IV)
        end
    end

    for k = 1:Mbd
        @views v = A[:, Lₚₙ+N*(k-1)+1:Lₚₙ+N*k]
        Axbdop!(v, k, d, dp, s, IV)
    end

    return nothing
end

constuctA!(A, IntS, d, dp, s, IV)

@btime constuctA!($A, $IntS, $d, $dp, $s, $IV)

@btime gmres($A, $b; reltol = 0.0, abstol = 1e-15, log = true)

using MAT
Af = matopen("A_ref.mat"); 
Aref = read(Af, "A");
close(Af)

for i in 1:Np*M
err = abs.(Aref[:,i] .- A[:,i]);
err_abs = maximum(err);
display(err_abs)
end

for i in Np*M+1:Lₚ
#i = Np*M+2
err = abs.(Aref[:,i] .- A[:,i])
err_abs = maximum(err);
display(err_abs)
end

IntS = precompsH(d,dp,0.75,4);

using MAT
f = matopen("I_ref.mat"); 
I_ref = read(f, "Int");
close(f)

err = abs.(IntS .- I_ref);
err_abs = maximum(abs.(IntS .- I_ref))
err_rel = err_abs / maximum(abs.(I_ref))

@btime precompsH($d, $dp, 0.75, 4);

@btime compress_vars(d, dp.N, 0.75, 4, 1);

Profile.clear()
@profile precompsH(d, dp, 0.75, 4)
ProfileView.view()

@benchmark precompsH($d,$dp,0.75,4) evals=1 samples=10 seconds=1000

@btime precompsH($d,$dp,0.75,4)

fig = Figure()
ax = Axis(fig[1,1])

heatmap!(ax, err, colormap = :viridis)


fig = Figure(size=(650, 650))
ax = Axis(fig[1, 1];
    #aspect=DataAspect(),
    xlabel=latexstring("\$x\$"),
    ylabel=latexstring("\$y\$"),
    xlabelsize=22,      # label font sizes
    ylabelsize=22,
    xticklabelsize=18,  # tick label sizes
    yticklabelsize=18,
)
#hidespines!(ax, :top, :right)
ax.xlabel = "x"
ax.ylabel = "y"

t = range(-1,1,1001) |> collect;

s = 0.9 

p = 500

lines!(t,wfunc(p,t))

Z = qw3func(p,t,s);

Z2= dwfunc(p,t)./(wfunc(p,-t).^(2*s-1).*wfunc(p,t).^(1-s));

abs.(Z .- Z2)

sum(map(x -> isnan(x) || !isfinite(x),Z))

dp = domprop(7,0.1,0.01,d)

plotdp(dp,d; label=:int)

#plotns(dp,d,1)

# fig, _ = plotu(dp, d, 0.75; interpolate=true)

# save("myplot.png",fig; px_per_unit=10)

#@btime disc(b = [3,3,3,3,3]);

function build_dp()
    d = disc(b = [7,7,7,7,7])  
    return domprop(12,0.1,0.01,d)
end

@btime build_dp(); 

#t1 = time()

d = disc(b = [3,3,3,3,3])    

#d = disc(b = [7,7,7,7,7])  

dp = domprop(10,0.1,0.01,d)

include(joinpath(@__DIR__, "bvec.jl"));
function f!(F, x, y)
    #fill!(F, 1.0)
    @. F = sqrt(x^2 + y^2)
    return nothing
end

b = bvec(d,dp,0.75,f!)

#er = similar(b);

for i in 1:M*N*N
    display(i - (ceil(Int, i/(N*N))-1)*N*N)
end

for i in M*N*N+1:M*N*N+Mbd*N
    k₀ = ceil(Int, (i-M*N*N)/N )
    l = d.kd[k₀]
    display(l)
    #display(i - M*N*N - (k₀-1)*N)
end

d = disc(b = [1,1,1,1,1])     

dp = domprop(4,0.1,0.01,d)

@btime bvec($d,$dp,0.75,$f!);

# #@code_warntype bvec(d, dp, 0.75, f!)

plotfunc(dp,d,b; interpolate=true)

er = similar(b);

for i in eachindex(b)

    er[i] = abs.( b[i] - ( ((dp.tgtpts[i].x)^2+(dp.tgtpts[i].y)^2)^(3/2)-1)/9 )

end

plotfunc(dp,d,er; interpolate=false)

maximum(er)

N = 8;

u, v = 0.0, 0.0;

# z= collect(range(-1,1,N));

# u2 = repeat(z, 1, N);
# v2 = repeat(z', N, 1);

# k = 1

# Zx = similar(u2);
# Zy = similar(v2);

# FXY = similar(Zx)

# t1 = copy(z);
# t2 = copy(z);

# DJ = similar(t1)

# @btime Dmap!($DJ, $d, $t1, $t2, 3);

# @btime f!($FXY, $Zx, $Zy);

z= collect(range(-1,1,N));

u2 = repeat(z, 1, N);
v2 = repeat(z', N, 1);

k = 1

Zx = similar(u2);
Zy = similar(v2);

out=similar(u2);

du = cospi.((2*(1:N).-1)./(2*N));
dv = cospi.((2*(1:N).-1)./(2*N));

diff_map!(out,Zx,Zy,d,u,v,u2,v2,du,dv,k)


display(out)

@btime diff_map!($out,$Zx,$Zy,$d,$u,$v,$u2,$v2,$du,$dv,3)

# u = collect(range(-1,1,512));
# v = collect(range(-1,1,512));

# D = Dmap(d, u, v, 1);
# De = similar(u);
# Dmap!(De,d, u, v, 1);

# maximum(abs.(D.-De))
# function D1()
#     for _ in 1:100_000
#     Dmap(d, u, v, 1)
#     end
# end


# function D2()
#     D = similar(u);
#     for _ in 1:100_000
#     Dmap!(D,d, u, v, 1)
#     end
# end

#@btime D1(); 

#@btime D2(); 

# dp = domprop(10,0.1,0.01,d)

# f(x,y) = x.^2 + cos.(10 .* y.^2)

# plotfunc(dp, d, f)

# fv = f(getfield.(dp.tgtpts, :x), getfield.(dp.tgtpts, :y));

# plotfunc(dp, d, fv)

# plotfunc(dp, d, fv; interpolate=true)

# y = chkinvpts(dp,d);

# u(x,y) = (1 .- x.^2 .- y.^2).^(3/4) 

# plotu(dp, d, u)

# maximum(y)

# fig , _ = drawboundquad(d,0.1,1:6:30,1)

#plot(1:5104,err)

#elapsed_time = time() - t1

#plotdp(dp,d)

#println("Elapsed time: ", elapsed_time, " seconds")
