"""
Mutable struct squircle

A  — x-coordinate of squircle center (Float64)

B  — y-coordinate of squircle center (Float64)

P  — A parameter bigger than 2. By default 4

θ₁ — A parameter used to constuct the patches in squircle. (Float64)

R₁ — Scale of the squircle in x axis (Float64)

R₂ — Scale of the squircle in y axis (Float64)

L₁ — Length of the y axis side of the rectangle (Float64)

L₂ — Length of the x axis side of the rectangle (Float64)

Together, L1×L2 are the dimensions of the inscribed rectangle. 

nₕ — number of holes (nonnegative Int) 

kd — indices of patches touching boundary (Vector{Int})  

Nₚₐₜ — total number of patches (≥ 1)  (Int)

pths — Npat stuctures of type Patch. The first
       row `reg` represents region in which the patch is present.
       2,3 row represents partition ck values.
       4,5 row represents partition tk values.
       With this five data values, the patch and the mapping
       is uniquely found and can be used in further 
       calculations. Definitions of ck and tk given in
       constructor of class. It is a matrix of size 5*Npat. 

Qpts — 8*Npat matrix (filled by bd_quad)  

Qptsbd — 8*Npat matrix (filled by bd_quadbd)

"""