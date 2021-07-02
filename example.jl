# Program to generate the reflection coefficient plot in the real frequency space
L = 10000.
H = 800.
h = 200.
nev = 10
N = 5

ω = 2π/100;
Ω₁ = solveIceVibration(L, h, H, nev, N, ω);
H₁, F₁, Ref₁, RefModes₁, RefDiff₁, X, U, Lc = Ω₁
"""
Run command in the REPL
    plotIce(X,U,ω,[-2.2,2.2],Ref)
"""

# Real space coarse2fine (10 × 1) ω-Line => (200 × 1) ω-Line
a = 2π/300
b = 2π/30
Ω₂ = coarse2fine(a, b, 10 ,200, L, h, H, nev, N)
H₂, F₂, λ₂, Ref₂, ω₂ = Ω₂
"""
Run commands in the REPL
    plotMode(ω₂, λ₂, 4)
    plotMode(ω₂, λ₂, [1,5,10])
    plotRefCoeff(ω₂, Ref₂)
"""

# Complex space coarse2fine (5 × 5) ω-Grid => (200 × 200) ω-Grid
a = 2π/300
b = 2π/30
c=-0.02
d=0.02
Ω₃ = coarse2fine(a, b, c, d, 5 ,200, L, h, H, nev, N)
H₃, F₃, λ₃, Ref₃, ω₃ = Ω₃

"""
Run command without @show in the REPL
    heatmap(angle.(Ref₃), colormap=:ascii)
    plotComplexRefCoeff(a, b, c, d, Ref₃)
"""
