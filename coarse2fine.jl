function coarse2fine(a, b, npts, npts_new, L, h, H, nev, N)
    ω = LinRange(a,b,npts)
    ω_new = LinRange(a,b,npts_new)
    H_old = zeros(ComplexF64, nev^2, npts)
    F_old = zeros(ComplexF64, nev, npts)
    RefDiff_ω_old = zeros(ComplexF64, npts)
    RefModes_ω_old = zeros(ComplexF64, nev, npts)

    print("Total of "*string(npts)*" FEM problems on the Real-ω Line\n")
    print("Solving problem ... ")
    for i ∈ 1:npts
        print(string(i)*"..")
        Hω,Fω,Refω,RefModes_ω,RefDiff_ω, X_ω, U_ω, Lc = solveIceVibration(L, h, H, nev, N, ω[i]);
        H_old[:,i] = collect(Iterators.flatten(Hω))
        F_old[:,i] = Fω
        RefDiff_ω_old[i] = RefDiff_ω
        RefModes_ω_old[:,i] = RefModes_ω
    end

    H_new = zeros(ComplexF64, nev^2, npts_new)
    F_new = zeros(ComplexF64, nev, npts_new)
    λ_new = zeros(ComplexF64, nev, npts_new)
    RefDiff_ω_new = zeros(ComplexF64, npts_new)
    RefModes_ω_new = zeros(ComplexF64, nev, npts_new)
    nodes = (collect(ω),)
    print("\nInterpolating on a "*string(npts_new)*" Real-ω Line\n")

    # Interpolate Diffraction coefficient
    itp_ref = Interpolations.interpolate(nodes, RefDiff_ω_old, Gridded(Linear()))
    RefDiff_ω_new = itp_ref.(ω_new)

    for i ∈ 1:nev^2
        itp = Interpolations.interpolate(nodes, H_old[i,:], Gridded(Linear()))
        H_new[i,:] = itp.(ω_new)
        if(i ≤ nev)
            itp = Interpolations.interpolate(nodes, F_old[i,:], Gridded(Linear()))
            F_new[i,:] = itp.(ω_new)
            # Radiation Reflection coefficients
            itp_ref = Interpolations.interpolate(nodes, RefModes_ω_old[i,:], Gridded(Linear()))
            RefModes_ω_new[i,:] = itp_ref.(ω_new)
        end
    end

    Refω_new = zeros(ComplexF64, npts_new)
    for i ∈ 1:npts_new
        λ_new[:,i] = reshape(H_new[:,i],(nev,nev))\F_new[:,i]
        Refω_new[i] = RefDiff_ω_new[i] + conj(RefModes_ω_new[:,i])⋅λ_new[:,i]
    end

    return H_new, F_new, λ_new, Refω_new, ω_new

end

function coarse2fine(a, b, c, d, npts, npts_new, L, h, H, nev, N)
    ω = LinRange(a,b,npts)' .* ones(npts) + (ones(npts)' .* LinRange(c,d,npts))*im
    H_old = zeros(ComplexF64, nev^2, npts^2)
    F_old = zeros(ComplexF64, nev, npts^2)
    RefDiff_ω_old = zeros(ComplexF64, npts^2)
    RefModes_ω_old = zeros(ComplexF64, nev, npts^2)

    print("Total of "*string(npts^2)*" FEM problems on the Complex ω-grid\n")
    print("Solving problem ... ")
    for i ∈ 1:npts^2
        print(string(i)*"..")
        Hω,Fω,Refω,RefModes_ω,RefDiff_ω, X_ω, U_ω, Lc = solveIceVibration(L, h, H, nev, N, ω[i]);
        H_old[:,i] = collect(Iterators.flatten(Hω))
        F_old[:,i] = Fω
        RefDiff_ω_old[i] = RefDiff_ω
        RefModes_ω_old[:,i] = RefModes_ω
    end

    # Interpolate to new Frequency Space.
    H_new = zeros(ComplexF64, nev^2, npts_new^2)
    F_new = zeros(ComplexF64, nev, npts_new^2)
    λ_new = zeros(ComplexF64, nev, npts_new^2)
    RefModes_ω_new = zeros(ComplexF64, nev, npts_new^2)
    nodes = (LinRange(a,b,npts), LinRange(c,d,npts))

    print("\nInterpolating on a "*string(npts_new)*"×"*string(npts_new)*" Complex ω-Grid\n")

    # Interpolate Diffraction coefficient
    itp_ref = Interpolations.interpolate(nodes, reshape(RefDiff_ω_old,(npts,npts)), Gridded(Linear()))
    X = itp_ref(LinRange(a,b,npts_new), LinRange(c,d,npts_new))
    RefDiff_ω_new = reshape(X,(npts_new^2,))

    for i ∈ 1:nev^2
        itp = Interpolations.interpolate(nodes, reshape(H_old[i,:],(npts,npts)), Gridded(Linear()))
        X = itp(LinRange(a,b,npts_new), LinRange(c,d,npts_new))
        H_new[i,:] = reshape(X,(npts_new^2,))
        if(i ≤ nev)
            itp = Interpolations.interpolate(nodes, reshape(F_old[i,:],(npts,npts)), Gridded(Linear()))
            X = itp(LinRange(a,b,npts_new), LinRange(c,d,npts_new))
            F_new[i,:] = reshape(X,(npts_new^2,))
            # Radiation Reflection coefficients
            itp_ref = Interpolations.interpolate(nodes, reshape(RefModes_ω_old[i,:],(npts,npts)), Gridded(Linear()))
            X = itp_ref(LinRange(a,b,npts_new), LinRange(c,d,npts_new))
            RefModes_ω_new[i,:] = reshape(X,(npts_new^2,))
        end
    end

    Refω_new = zeros(ComplexF64, npts_new^2)
    for i ∈ 1:npts_new^2
        λ_new[:,i] = reshape(H_new[:,i],(nev,nev))\F_new[:,i]
        Refω_new[i] = RefDiff_ω_new[i] + conj(RefModes_ω_new[:,i])⋅λ_new[:,i]
    end

    ω_new = LinRange(a,b,npts_new)' .* ones(npts_new) + (ones(npts_new)' .* LinRange(c,d,npts_new))*im
    Refω_new = reshape(Refω_new, (npts_new,npts_new))
    return H_new, F_new, λ_new, Refω_new, ω_new

end
