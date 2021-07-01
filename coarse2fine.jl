function coarse2fine(a, b, npts, npts_new, L, h, H, nev, N)
    ω = LinRange(a,b,npts)
    ω_new = LinRange(a,b,npts_new)

    H_old = zeros(ComplexF64, nev^2, npts)
    F_old = zeros(ComplexF64, nev, npts)
    print("Total of "*string(npts)*" FEM problems on the Real-ω Line\n")
    print("Solving problem ... ")
    for i ∈ 1:npts
        print(string(i)*"..")
        Hω,Fω,Refω,RefModes_ω,RefDiff_ω, X_ω, U_ω, Lc = solveIceVibration(L, h, H, nev, N, ω[i]);
        #print("\n")
        H_old[:,i] = collect(Iterators.flatten(Hω))
        F_old[:,i] = Fω
    end

    H_new = zeros(ComplexF64, nev^2, npts_new)
    F_new = zeros(ComplexF64, nev, npts_new)
    λ_new = zeros(ComplexF64, nev, npts_new)
    nodes = (collect(ω),)

    print("\nInterpolating on a "*string(npts_new)" Real-ω Line\n")
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

function coarse2fine(a, b, c, d, npts, npts_new, L, h, H, nev, N)
    ω = LinRange(a,b,npts)' .* ones(npts) + (LinRange(c,d,npts)' .* ones(npts))*im
    H_old = zeros(ComplexF64, nev^2, npts^2)
    F_old = zeros(ComplexF64, nev, npts^2)
    print("Total of "*string(npts^2)*" FEM problems on the Complex ω-grid\n")
    print("Solving problem ... ")
    for i ∈ 1:npts^2
        print(string(i)*"..")
        Hω,Fω,Refω,RefModes_ω,RefDiff_ω, X_ω, U_ω, Lc = solveIceVibration(L, h, H, nev, N, ω[i]);
        H_old[:,i] = collect(Iterators.flatten(Hω))
        F_old[:,i] = Fω
    end

    # Interpolate to new Frequency Space.
    H_new = zeros(ComplexF64, nev^2, npts_new^2)
    F_new = zeros(ComplexF64, nev, npts_new^2)
    λ_new = zeros(ComplexF64, nev, npts_new^2)
    nodes = (LinRange(a,b,npts), LinRange(c,d,npts))

    print("\nInterpolating on a "*string(npts_new)*"×"*string(npts_new)*" Complex ω-Grid\n")
    for i ∈ 1:nev^2
        itp = Interpolations.interpolate(nodes, reshape(H_old[i,:],(npts,npts)), Gridded(Linear()))
        X = itp(LinRange(a,b,npts_new), LinRange(c,d,npts_new))
        H_new[i,:] = reshape(X,(npts_new^2,))
        if(i ≤ nev)
            itp = Interpolations.interpolate(nodes, reshape(F_old[i,:],(npts,npts)), Gridded(Linear()))
            X = itp(LinRange(a,b,npts_new), LinRange(c,d,npts_new))
            F_new[i,:] = reshape(X,(npts_new^2,))
        end
    end

    for i ∈ 1:npts_new^2
        λ_new[:,i] = reshape(H_new[:,i],(nev,nev))\F_new[:,i]
    end

    return H_new, F_new, λ_new, LinRange(a,b,npts_new)' .* ones(npts_new) + (LinRange(c,d,npts_new)' .* ones(npts_new))*im

end
