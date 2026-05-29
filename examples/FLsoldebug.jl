using Revise, FL2D

import FL2D.FLdata as FLdata

d = FL2D.peanut(b=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

d = FL2D.peanut(b=[6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]);

FL2D.refine!(d, 2, 2, [67,96,187,216])



d = FL2D.kite(b=[4, 5, 9, 9, 5, 4, 5, 5, 5, 5, 4, 4],
a=[4, 3, 7, 7, 3, 4, 3, 6, 6, 3, 3, 3]);

FL2D.refine!(d, 1, 2, [9, 12, 13, 16, 181, 184, 185, 188])

d = FL2D.kite(b=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])


d = FL2D.bean(b=[4, 5, 5, 5, 5, 4, 5, 9, 9, 5, 4, 4],
a=[4, 3, 7, 7, 3, 4, 3, 6, 6, 3, 4, 4]);

FL2D.refine!(d, 1, 2, [9, 12, 13, 16, 125, 128, 129, 132])

d = FL2D.bean(b=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])



d =  FL2D.disc(b=[6, 6, 6, 6, 5],a=[4, 4, 4, 4, 4], L1=0.8, L2=0.8)

d = FL2D.disc(b=[2, 2, 2, 2, 2], L1=0.8, L2=0.8)


d =  FL2D.squircle(b=[6, 6, 6, 6, 5],a=[4, 4, 4, 4, 4], L1=0.8, L2=0.8, P=4.0)

d = FL2D.squircle(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8, P=4.0)



d = FL2D.ellipse(
      b=[2, 2, 2, 2, 2],
      a=[1, 1, 1, 1, 1],
      R1=1.0, R2=2.0, L1=2.0, L2=0.8)

d = FL2D.ellipse(b=[6, 6, 6, 6, 6],a=[3, 3, 3, 3, 4],R1=1,R2=2,L1=2,L2=0.8)



d = FL2D.star(b=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],P=4,tht1=0.24,tht2=0.225,L1=0.5,L2=0.5)

d = FL2D.star(b=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

d = FL2D.star(b=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],P=6,tht1=0.14,tht2=0.15,L1=0.3,L2=0.52)


d = FL2D.ellipseNh(b=[1, 1, 1, 1, 1, 1, 1, 1])


dp = FL2D.domprop(12, 0.15, 8e-3, 0.15, 5e-3, d; Lᵢₙ=5)

FL2D.memory_report(dp)

fig, ax = FL2D.plotns(dp, d, 1)

FL2D.drawextregion(fig, ax, d, 1, 1.3)

FL2D.chk_map(d);

FL2D.draw(d,1)

FL2D.drawbd(d)

FL2D.chkinvpts(dp,d; flag=true)

FL2D.testDLP(d, dp; nr=64)

f!, uex, fv = FLdata.makesquirclefuex(5, 0.1, 4);

f!, uex, fv = FLdata.makediscfuex(5, 0.75);

f!, uex, fv = FLdata.makepeanutfuex(0, 0.75);

f!, uex, fv = FLdata.makebeanfuex(5, 0.75);

f!, uex, fv = FLdata.makestarfuex(0, 0.75);

FL2D.plotfunc(dp, d, f!)

bex = FL2D.bvec(d, dp, 0.75, f!);

FL2D.plotfunc(dp,d,bex)

for n in [16,32]
   b = FL2D.bvec(d, dp, 0.75, f!; n=n);
   display(maximum(abs.(b.-bex)))
end  

s=0.25
p=4 

IntS = FL2D.precompsL(d, dp, s, p; n=128);

IV = FL2D.compress_vars(d, dp, s, p; matrix_form=true);

(; N, Np, M, Mbd, Ni) = IV.IV1;

Lp = Ni + Mbd * N;

A = zeros(Float64, Lp, Lp);

for k in 1:M
   @views v = A[:, (1+Np*(k-1)):(Np*k)]
   if !(k in d.kd)
      FL2D.Axintpth!(v, k, IntS, d, dp, s, IV)
   end
end

for k in 1:M
   @views v = A[:, (1+Np*(k-1)):(Np*k)]
   if k in d.kd
      FL2D.Axbdpth!(v, k, IntS, d, dp, s, IV)
   end
end



for k in 1:Mbd
   @views v = A[:, Ni+N*(k-1)+1:Ni+N*k]
   Axbdop!(v, k, d, dp, s, IV)
end

bex = FL2D.bvec(d, dp, 0.75, f!, n=256); 
for n in [32, 64, 128] 
   b = FL2D.bvec(d, dp, 0.75, f!; n=n); 
   display(maximum(abs.(b.-bex))) 
end

d = FL2D.disc(b=[5, 5, 5, 5, 5], a=[3, 3, 3, 3, 4], L1=0.8, L2=0.8)

d = FL2D.disc(b=[1, 1, 1, 1, 1], L1=0.8, L2=0.8)

dp = FL2D.domprop(12, 0.15, 8e-3, 0.15, 5e-3, d; Lᵢₙ=5)

s, p, n = 0.1, 4, 128

IntSex = FL2D.precompsLs(d, dp, s, p; n=512);

for n in [16,32,64,128]
   IntS = FL2D.precompsLs(d, dp, s, p; n=n);
   display(maximum(abs.(IntS .- IntSex)))
end

IntS = FL2D.precompsLs(d, dp, s, p; n=256);
display(maximum(abs.(IntS .- IntSex)))

for k in 1:d.Npat
    S = dp.pthgo[k] + dp.N^2
    E = dp.pthgo[k+1] - 1
    println(k, " ", maximum(abs.(IntS[:, S:E] .- IntSex[:, S:E])))
end

# Singular/self blocks only
for k in 1:d.Npat
    S = dp.pthgo[k]
    E = dp.pthgo[k] + dp.N^2 - 1
    abs_err = maximum(abs.(IntS[:, S:E] .- IntSex[:, S:E]))
    rel_err = abs_err / maximum(abs.(IntSex[:, S:E]))
    println(k, " abs=", abs_err, " rel=", rel_err)
end

# Apply log10 transformation. 
# We add 1e-16 to safeguard against log10(0) if your matrix contains exact zeros.
log_errors = log10.(abs.(IntS .- IntSex)/ maximum(abs.(IntSex)) .+ 1e-16);

# Define axes using the matrix indices
# In Makie, dimension 1 (rows) maps to the X-axis, dimension 2 (columns) maps to the Y-axis
row_indices = 1:size(IntS, 1);
col_indices = 1:size(IntS, 2);


fig = Figure(size = (1300, 600), fontsize = 16)

# Left Column: 2D Heatmap Axis mapped to element indices
ax2d = Axis(fig[1, 1], 
    title = "2D Error Heatmap (By Index)", 
    xlabel = "Matrix Row Index", 
    ylabel = "Matrix Column Index"
)

# Right Column: 3D Surface Axis mapped to element indices
ax3d = Axis3(fig[1, 3], 
    title = "3D Error Landscape (By Index)", 
    xlabel = "Matrix Row Index", 
    ylabel = "Matrix Column Index", 
    zlabel = "log10(Error)",
    aspect = (1, 1, 0.8), 
    perspectiveness = 0.5 
)


# ==========================================
# 3. Plotting using Indices
# ==========================================
cmap = :inferno  # Perceptually uniform colormap perfect for scale visualization

# Pass the index ranges as the X and Y coordinates
hm = heatmap!(ax2d, row_indices, col_indices, log_errors, colormap = cmap)
sf = surface!(ax3d, row_indices, col_indices, log_errors, colormap = cmap)


# ==========================================
# 4. Colorbar and Render
# ==========================================
cb = Colorbar(fig[1, 2], hm, 
    label = "log10(Absolute Error)", 
    width = 20, 
    tickalign = 1
)

# Display the interactive window
display(fig)
