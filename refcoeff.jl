function getRefCoeff(phi, NModes, k, kd, HH, dd, Γ, Ap)
    aa = zeros(Complex{Float64}, NModes+1, 1)

    dΓ = Measure(Γ,2)
    for i ∈ 1:NModes+1
        τ(x) = cos(kd[i]*(x[2]+HH))/cos(kd[i]*(HH-dd))
        aa[i] = sum(∫(τ*phi)*dΓ)
    end

    A, M, f, g = getMAT(k, kd, HH, dd, NModes, Ap)
    Mt = transpose(M)
    T = inv(Mt)*A*inv(M)

    bb = T*aa + inv(Mt)*g - T*f
    c = inv(A)*(Mt*bb - g)
    Ref = c/Ap
end

function getRefModes(phi, NModes, k, kd, HH, dd, Γ, Ap)
    aa = zeros(Complex{Float64}, NModes+1, 1)

    dΓ = Measure(Γ,2)
    for i ∈ 1:NModes+1
        τ(x) = cos(kd[i]*(x[2]+HH))/cos(kd[i]*(HH-dd))
        aa[i] = sum(∫(τ*phi)*dΓ)
    end

    A, M, f, g = getMAT(k, kd, HH, dd, NModes, Ap)
    Mt = transpose(M)
    T = inv(Mt)*A*inv(M)

    bb = T*aa
    c = inv(A)*Mt*bb
    Ref = c/Ap
end
