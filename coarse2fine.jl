function coarse2fine(a, b, npts, npts_new, L, h, H, nev, N)
    ω = LinRange(a,b,npts)
    ω_new = LinRange(a,b,npts_new)

    H_old = zeros(ComplexF64, nev^2, npts)
    F_old = zeros(ComplexF64, nev, npts)
    for i ∈ 1:npts
        Hω,Fω,Refω,RefModes_ω,RefDiff_ω, X_ω, U_ω, Lc = solveIceVibration(L, h, H, nev, N, ω[i]);
        #print("\n")
        H_old[:,i] = collect(Iterators.flatten(Hω))
        F_old[:,i] = Fω
    end

    H_new = zeros(ComplexF64, nev^2, npts_new)
    F_new = zeros(ComplexF64, nev, npts_new)
    λ_new = zeros(ComplexF64, nev, npts_new)
    nodes = (collect(ω),)

    for i ∈ 1:nev^2
        itp = Interpolations.interpolate(nodes, H_old[i,:], Gridded(Linear()))
        H_new[i,:] = itp.(ω_new)
        if(i ≤ nev)
            itp = Interpolations.interpolate(nodes, F_old[i,:], Gridded(Linear()))
            F_new[i,:] = itp.(ω_new)
        end
    end

    for i ∈ 1:npts_new
        λ_new[:,i] = reshape(H_new[:,i],(nev,nev))\F_new[:,i]
    end

    return H_new, F_new, λ_new, ω_new

end
