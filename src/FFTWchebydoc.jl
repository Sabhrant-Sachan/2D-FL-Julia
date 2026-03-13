using FFTW

# ============================================================
# 2D Chebyshev interpolation using FFTW DCTs in Julia
#
# Goal:
#   1. Sample a function f(x,y) on an N×N Chebyshev grid
#   2. Convert those sampled values to Chebyshev coefficients
#   3. Zero-pad the coefficient matrix to a larger size N2×N2
#   4. Transform back to values on the finer Chebyshev grid
#
# Important idea:
#   Zero-padding the Chebyshev coefficients does NOT add new
#   information. It simply evaluates the same interpolating
#   polynomial on a finer grid.
#
# Recall: 
# https://www.fftw.org/fftw3_doc/1d-Real_002deven-DFTs-_0028DCTs_0029.html
#   FFTW's REDFT10 is the unnormalized DCT-II
#   FFTW's REDFT01 is the unnormalized DCT-III
#
# For our Chebyshev convention:
#
#   Forward (values -> coefficients):
#       c = (1/N) * DCT-II(values)
#       c[1] /= 2
#
#   Inverse (coefficients -> values):
#       c[1] *= 2
#       values = 0.5 * DCT-III(c)
#
# In 2D, we apply the same logic one dimension at a time.
# ============================================================

# Number of Chebyshev points in each direction for the original grid
N = 32

# Number of Chebyshev points in each direction for the finer grid
# We will zero-pad the coefficient matrix from N×N to N2×N2.
N2 = 64

# Number of total grid points in the original grid
Np = N * N


# ------------------------------------------------------------
# Test function
# ------------------------------------------------------------

# A smooth analytic function.
# Chebyshev coefficients for such functions decay very quickly,
# which is why interpolation works so well with modest N.
#f(x::Float64, y::Float64) = exp(sin(x) * cos(y))

# test with a polynomial:
 f(x::Float64, y::Float64) = x^2 + y^2 
#  = T_0(x)T_0(y) + 1/2(T_2(x)T_0(y)+T_0(x)T_2(y))
#
# For a low-degree polynomial, the Chebyshev interpolant becomes
# exact once N is large enough to represent all needed modes.


# ------------------------------------------------------------
# Chebyshev nodes of the first kind (Gauss-Chebyshev points)
# ------------------------------------------------------------
#
# For DCT-II based Chebyshev interpolation, the grid points are:
#
#   x_i = cos( (2i-1)π / (2N) ),   i = 1,...,N
#
# In Julia, cospi(z) computes cos(π*z), which is usually more
# accurate than writing cos(pi*z).
# ------------------------------------------------------------

y = Vector{Float64}(undef, N)

for i = 1:N
    y[i] = cospi((2*i - 1) / (2*N))
end


# ------------------------------------------------------------
# Sample the function on the N×N Chebyshev grid
# ------------------------------------------------------------
# See : https://prnt.sc/_EK7Tjm_T1dM
# FV[i,j] stores f(x_j, y_i)
#
# So:
#   - columns correspond to x-direction
#   - rows correspond to y-direction
# ------------------------------------------------------------

FV = Matrix{Float64}(undef, N, N)

for i = 1:N
    for j = 1:N
        FV[i, j] = f(y[j], y[i])
    end
end


# ------------------------------------------------------------
# Forward transform: values -> Chebyshev coefficients
# ------------------------------------------------------------
#
# Step 1:
#   Apply DCT-II along columns (dim = 2), then scale by 1/N,
#   then halve the first coefficient in that transformed direction.
#
# Step 2:
#   Apply DCT-II along rows (dim = 1), then scale by 1/N,
#   then halve the first coefficient in that transformed direction.
#
# After both steps, CN[k,l] is the 2D Chebyshev coefficient matrix.
# ------------------------------------------------------------

# Coefficient matrix
CN = similar(FV)

# Transform along dim = 2 (across columns)
CN = (1 / N) .* FFTW.r2r(FV, FFTW.REDFT10, 2)

# In Chebyshev normalization, the first mode in that direction
# gets an extra factor of 1/2.
CN[:, 1] ./= 2

# Transform along dim = 1 (across rows)
CN = (1 / N) .* FFTW.r2r(CN, FFTW.REDFT10, 1)

# Again halve the first mode in this direction
CN[1, :] ./= 2


# ------------------------------------------------------------
# Optional: pack coefficients into a vector
# ------------------------------------------------------------
#
# Sometimes PDE / spectral codes store coefficients in a long vector.
# Julia uses column-major ordering, so vec(CN) stacks columns.
# ------------------------------------------------------------

CNV = Vector{Float64}(undef, Np)
CNV .= vec(CN)


# ------------------------------------------------------------
# Inverse transform on the SAME grid: coefficients -> values
# ------------------------------------------------------------
# (This is just a correctness check.)
# To invert the forward transform:
#
#   1. Undo the "halve first coefficient" in reverse order
#   2. Apply DCT-III (REDFT01)
#   3. Multiply by 0.5 after each inverse transform
#
# Reverse order matters because the forward transform was applied
# one dimension at a time.
# ------------------------------------------------------------

G = copy(CN)

G[:, 1] .*= 2

# Inverse of DCT-II along dim = 2 is DCT-III, with factor 0.5
G = 0.5 .* FFTW.r2r(G, FFTW.REDFT01, 2)

# Now undo the scaling from the dim = 1 stage
G[1, :] .*= 2

# Invert along dim = 1
G = 0.5 .* FFTW.r2r(G, FFTW.REDFT01, 1)

# Reconstruction error on the original grid
println("Max reconstruction error on original grid = ", maximum(abs.(G .- FV)))


# ------------------------------------------------------------
# Zero-pad the coefficient matrix to a larger size N2×N2
# ------------------------------------------------------------
#
# This creates a larger Chebyshev expansion where:
#
#   - all low modes (the original coefficients) are preserved
#   - all higher modes are set to zero
#
# This does NOT create new information.
# It only evaluates the same truncated polynomial on a finer grid.
# ------------------------------------------------------------

CNp = zeros(Float64, N2, N2)

# Copy original coefficients into the upper-left corner.
# This is the correct place because coefficient indices simply
# represent mode numbers.
CNp[1:N, 1:N] .= CN


# ------------------------------------------------------------
# Inverse transform on the finer N2×N2 grid
# ------------------------------------------------------------
#
# Same inverse logic as before, but now applied at size N2.
# The output is the interpolating polynomial evaluated on the
# finer Chebyshev grid.
# ------------------------------------------------------------

Gfine = copy(CNp)

# Undo the first-mode scaling for dim = 2
Gfine[:, 1] .*= 2

# Invert along dim = 2
Gfine = 0.5 .* FFTW.r2r(Gfine, FFTW.REDFT01, 2)

# Undo the first-mode scaling for dim = 1
Gfine[1, :] .*= 2

# Invert along dim = 1
Gfine = 0.5 .* FFTW.r2r(Gfine, FFTW.REDFT01, 1)


# ------------------------------------------------------------
# Build the exact function values on the finer N2×N2 grid
# ------------------------------------------------------------
# (To measure interpolation error.)
# ------------------------------------------------------------

y2 = Vector{Float64}(undef, N2)

for i = 1:N2
    y2[i] = cospi((2*i - 1) / (2*N2))
end

FV2 = Matrix{Float64}(undef, N2, N2)

for i = 1:N2
    for j = 1:N2
        FV2[i, j] = f(y2[j], y2[i])
    end
end

# ------------------------------------------------------------
# Compare interpolated values against exact function values
# ------------------------------------------------------------
#
# For a polynomial of sufficiently low degree, this will be
# near machine precision.For a non-polynomial analytic function, 
# this measures the truncation/interpolation error of the N×N 
# Chebyshev expansion, evaluated on the finer N2×N2 grid.
# ------------------------------------------------------------

println("Max interpolation error on finer grid = ", maximum(abs.(Gfine .- FV2)))
