# Program to generate the reflection coefficient plot in the real frequency space
L = 10000.
H = 800.
h = 200.
nev = 10
N = 5

ω = 2π/100;
Hmat, F, Ref, RefModes, RefDiff, X, U, Lc = solveIceVibration(L, h, H, nev, N, ω);

# Run command without @show in the REPL
@show plotIce(X,U,ω,[-2.2,2.2],Ref)


# Real space coarse2fine
H_new, F_new, λ_new, ω_new = coarse2fine(2π/400, 2π/20, 10 ,200, L, h, H, nev, N)

# Run command without @show in the REPL
@show plotMode(ω_new, λ_new, 4)

# Complex space coarse2fine
H_new, F_new, λ_new, ω_new = coarse2fine(2π/400, 2π/20, -0.01, 0.01, 3 ,200, L, h, H, nev, N)
