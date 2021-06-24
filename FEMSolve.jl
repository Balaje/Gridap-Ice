# Solve the ice-shelf 2d eigenvalue problem. [ε is a function to store the symmetric gradient]
function solveEigen(iceModel,Vh,Vh0,order,N)
    λ=3500;
    ν=0.3;
    μ=λ*(1-2*ν)/(2*ν);
    ρᵢ=1;
    # The stress tensor
    σ(ε)=λ*tr(ε)*one(ε) + 2*μ*ε;
    Ω=Triangulation(iceModel);
    dΩ=Measure(Ω,2*order);
    # Weak form
    a(u,v) = ∫(ε(v)⊙(σ∘ε(u)))*dΩ;
    m(u,v) = ∫(ρᵢ*u⋅v)*dΩ;
    l(v) = 0;
    # Get the matrices
    opK = AffineFEOperator(a,l,Vh,Vh0);
    opM = AffineFEOperator(m,l,Vh,Vh0);
    K=opK.op.matrix;
    M=opM.op.matrix;
    # Use Arpack to solve the Eigenvalue problem
    ξ,Vec = eigs(K,M; nev=N, which=:LM);
    return ξ,Vec;
end


# Solve the velocity potentials.
function getLaplaceMatEB(Ω, Γ₃, Vh, Vh0, QΦ, χ, μₘ, L, ω)
    # Measure of the domains
    dΩ=Measure(Ω,2);
    dΓ₃=Measure(Γ₃,6); #Interface boundary.

    η(x) = ω*((cos(L*μₘ) + cosh(L*μₘ))*(sin(μₘ*x[1]) + sinh(μₘ*x[1]))-
              (sin(L*μₘ) + sinh(L*μₘ))*(cos(μₘ*x[1]) + cosh(μₘ*x[1])))/
              (cos(L*μₘ) + cosh(L*μₘ));

    a(u,v) = ∫( ∇(v)⊙∇(u) )*dΩ;
    b(v) = ∫(η*v)*dΓ₃;
    op=AffineFEOperator(a,b,Vh,Vh0);
    K=op.op.matrix+QΦ;
    #print(norm(op.op.vector),"\n")
    f=-(1im)*op.op.vector-χ[1,:];
    return K,f,op
end

# Function to build the reduced system
function buildReducedSystem(μ, ϕ₀, ϕⱼ, α, β, γ, Γ, L, ω, V)
    dΓ=Measure(Γ,6);
    nev=length(μ);

    B=zeros(ComplexF64, nev, nev);
    K=zeros(ComplexF64, nev, nev);
    AB=zeros(ComplexF64, nev, nev);
    F=zeros(ComplexF64,nev,1);

    φ₀=FEFunction(V, real(ϕ₀))+1im*FEFunction(V, imag(ϕ₀));
    for i=1:nev
        μₘ=μ[i];
        η(x)=((cos(L*μₘ) + cosh(L*μₘ))*(sin(μₘ*x[1]) + sinh(μₘ*x[1]))-
              (sin(L*μₘ) + sinh(L*μₘ))*(cos(μₘ*x[1]) + cosh(μₘ*x[1])))/
              (cos(L*μₘ) + cosh(L*μₘ));

        B[i,i]=(1-γ*α)*(cosh(L*μₘ)*sin(L*μₘ) - cos(L*μₘ)*sinh(L*μₘ) -
                        (cos(L*μₘ)^2*sinh(2*L*μₘ))/2 + (cosh(L*μₘ)^2*sin(2*L*μₘ))/2
                        - L*μₘ*cos(L*μₘ)^2 + L*μₘ*cosh(L*μₘ)^2 + 2*L*μₘ*sin(L*μₘ)*
                        sinh(L*μₘ))/(μₘ*(cos(L*μₘ) + cosh(L*μₘ))^2);
        K[i,i]=β*μₘ^4*B[i,i];
        F[i]=(1im*ω/10)*sum(∫(η*φ₀)*dΓ);

        φₖ=FEFunction(V,real(ϕⱼ[:,i])) + 1im*FEFunction(V,imag(ϕⱼ[:,i]));
        for j=1:nev
            μₘ=μ[j];
            ξ(x)=((cos(L*μₘ) + cosh(L*μₘ))*(sin(μₘ*x[1]) + sinh(μₘ*x[1]))-
                  (sin(L*μₘ) + sinh(L*μₘ))*(cos(μₘ*x[1]) + cosh(μₘ*x[1])))/
                  (cos(L*μₘ) + cosh(L*μₘ));

            AB[i,j]=-(1im)*(ω/10)*sum(∫(ξ*φₖ)*dΓ);
        end
    end
    H=K+B+AB;
    λ=H\F;
    return λ, K, B, AB, F;
end
