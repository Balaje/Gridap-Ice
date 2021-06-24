# Some parameters
#ω=2*π/200; # 40s incident wave.
#N=5; # Modal expansion in the ocean
#nev=20; #Number of eigenvalues
#L=10000; #Shelf length
#h=200; #Shelf thickness
#d=0.9*h; #Submergence.
#H=800; #Ocean depth
function solveIceVibration(L=10000, h=200, H=800, nev=10, N=5, ω=2*π/200)
    d=0.9*h
    E=2.0e9;
    ρᵢ=922.5;
    ρₗ=1025;
    ν=0.33;
    EI=E*h^3/(12*(1-ν^2));
    Lc=(EI/(ρₗ*9.8))^0.25;
    tc=sqrt(ρₗ*Lc^6/(EI*H));
    LL=L/Lc;
    HH=H/Lc;
    hh=h/Lc;
    dd=d/Lc;
    α=HH*(ω*tc)^2;
    Ap=(9.8/(1im*ω));

    # Solve the dispersion equation
    k=dispersionfreesurface(α, N, HH);
    k[1]=-k[1];
    kd=dispersionfreesurface(α, N, HH-dd);
    kd[1]=-kd[1];
    #print("Solved dispersion equations\n")

    partition=(100,20);

    # Build model for the ice-shelf.
    #iceDomain=(0,L,-d,h-d);
    #iceModel=CartesianDiscreteModel(iceDomain,partition);
    #iceLabels=get_face_labeling(iceModel);
    #Ωs=Triangulation(iceModel); #Build the triangulation
    #add_tag_from_tags!(iceLabels,"neumannIce",[1,2,5])
    #add_tag_from_tags!(iceLabels,"dirichletIce",[4,8])
    #Γs₁=BoundaryTriangulation(iceModel,iceLabels,tags="neumannIce");
    #nΓs₁=get_normal_vector(Γs₁);# Get the normal
    #dΓs₁=Measure(Γs₁,2);
    #reffe=ReferenceFE(lagrangian,VectorValue{2,Float64},1);
    #Vsₕ=TestFESpace(iceModel,reffe,conformity=:H1,dirichlet_tags="dirichletIce",
    #                dirichlet_masks=(true,true)); #Test space for Ice

    ## Build model for the cavity region
    cavDomain=(0,LL,-HH,-dd);
    cavModel=CartesianDiscreteModel(cavDomain,partition);
    cavModel=simplexify(cavModel)
    cavLabels=get_face_labeling(cavModel);
    Ωf=Triangulation(cavModel); #Build the triangulation
    add_tag_from_tags!(cavLabels,"neumannIce",[6]);
    add_tag_from_tags!(cavLabels,"NonLocal",[3,7]);
    Γf₃=BoundaryTriangulation(cavModel,cavLabels,tags="neumannIce"); # Shelf/Cavity
    Γf₄=BoundaryTriangulation(cavModel,cavLabels,tags="NonLocal"); # Non-Local
    reffe=ReferenceFE(lagrangian,Float64,1);
    Vfₕ=TestFESpace(cavModel,reffe,conformity=:H1); #Test space for cavity.

    # ------ Get the non-local boundary condition
    Qϕ,χ=getMQχ(k, kd, HH, dd, N, Ap, cavModel, Γf₄, Vfₕ, Vfₕ);
    #print("Done computing non-local boundary condition\n")
    # -----------------------------------------

    # First attempt at solving the Elasticity Eigenvalue problem
    # Question on interpolation?
    #ξ,Vec=FEMSolvers.solveEigen(iceModel, Vsₕ, Vsₕ, 1, 10); # Returns the raw vectors.
    #uh=FEFunction(Vsₕ, real(Vec[:,1]))
    #writevtk(Ωs,"results",cellfields=["uh"=>uh]); #To visualize the solution.

    # ------- Solve for the velocity potentials
    #Diffraction Potential (with a rigid shelf)
    K,f,op=getLaplaceMatEB(Ωf, Γf₃, Vfₕ, Vfₕ, Qϕ, χ, 0, LL, 0)
    ϕ₀=K\f;
    # Compute diffraction refcoeffs
    real_uh=FEFunction(Vfₕ,real(ϕ₀[:,1]))
    imag_uh=FEFunction(Vfₕ,imag(ϕ₀[:,1]))
    real_Ref = getRefCoeff(real_uh, N, k, kd, HH, dd, Γf₄, Ap)
    imag_Ref = getRefCoeff(imag_uh, N, k, kd, HH, dd, Γf₄, Ap)
    RefDiff = zeros(ComplexF64, 2, 1)
    RefDiff[1] = real_Ref[1]
    RefDiff[2] = imag_Ref[1]

    # Radiation potential from the Eigenmodes
    μ=solveEigenEB(nev, LL);# Obtain the Eigenvalues for the beam equation
    ϕₖ=zeros(ComplexF64,length(χ),nev)
    RefModes=zeros(ComplexF64, 2, nev)
    for m=1:nev
        Km,fm,op1=getLaplaceMatEB(Ωf, Γf₃, Vfₕ, Vfₕ, Qϕ, 0*χ, μ[m], LL, ω*Lc)
        ϕₖ[:,m]=Km\fm; # Using the linear algebra package (raw vector)

        # Compute modal refcoeffs
        POTₖ=ϕₖ[:,m]
        real_uh=FEFunction(Vfₕ,real(POTₖ[:,1]))
        imag_uh=FEFunction(Vfₕ,imag(POTₖ[:,1]))
        real_Ref = getRefModes(real_uh, N, k, kd, HH, dd, Γf₄, Ap)
        imag_Ref = getRefModes(imag_uh, N, k, kd, HH, dd, Γf₄, Ap)
        RefModes[1,m] = real_Ref[1]
        RefModes[2,m] = imag_Ref[1]
    end
    #print("Done computing potentials\n")
    # -----------------------------------------

    # ------ Build and solve the reduced system
    λ,K,B,AB,F=buildReducedSystem(μ, ϕ₀, ϕₖ, α, 1, dd, Γf₃, LL, ω, Vfₕ);
    #print("Done computing coefficients\n")
    # ----------------------------------

    ## Construct the displacement
    function u(x, μ, L, λ)
        nev=length(μ);
        ξ=0;
        for m=1:nev
            μₘ=μ[m];
            ηₘ = ((cos(L*μₘ) + cosh(L*μₘ))*(sin(μₘ*x) + sinh(μₘ*x))-
                (sin(L*μₘ) + sinh(L*μₘ))*(cos(μₘ*x) + cosh(μₘ*x)))/
                (cos(L*μₘ) + cosh(L*μₘ));
            ξ = ξ+λ[m]*ηₘ;
        end
        return ξ;
    end

    X=collect(range(0, LL, length=200));
    U=zeros(ComplexF64,length(X),1);
    for m=1:length(U)
        U[m]=u(X[m], μ, LL, λ);
    end

    # Construct the velocity potential
    POT=ϕ₀+ϕₖ*λ
    #Compute the reflection coefficient
    real_uh=FEFunction(Vfₕ,real(POT[:,1]))
    imag_uh=FEFunction(Vfₕ,imag(POT[:,1]))
    #writevtk(Ωf,"results",cellfields=["uh"=>uh]); #To visualize the solution.
    real_Ref = getRefCoeff(real_uh, N, k, kd, HH, dd, Γf₄, Ap)
    imag_Ref = getRefCoeff(imag_uh, N, k, kd, HH, dd, Γf₄, Ap)
    Ref = real_Ref .+ 1im*(imag_Ref .+ 1) # Don't know how but works 😉

    H = K+AB+B
    return H, F, Ref, RefModes, RefDiff, X, U, Lc
end

function buildRefCoeff(RefDiff, RefModes, λ)
    X1=RefDiff + RefModes*λ
    (2 .+ X1)[1]
end


function plotIce(X,U,ω,ylim=[-2,2],Ref=:none)
    plt = lineplot(X*Lc,real(U[:,1]), width = 40, xlim = [minimum(X*Lc), maximum(X*Lc)],
                   ylim=ylim,
                   xlabel="x in m",
                   ylabel="η(x,ω) in m",
                   title="Displacement of Ice for incident T = "*string(round(2π/ω))*" s",
                   name="Real part",
                   border=:ascii)
    lineplot!(plt,X*Lc,imag(U[:,1]),color=:red, name="Imaginary Part")
    if(Ref==:none || typeof(Ref)!=Matrix{ComplexF64})
        return plt
    else
        annotate!(plt, :r, 3, "R(ω) = "*string(Ref[1]))
        annotate!(plt, :r, 4, "|R(ω)| = "*string(round(abs(Ref[1]))))
    end
    plt
end

function plotMode(ω, λ, N)
    plt=lineplot(ω, abs.(λ[N[1],:]), width = 80,
                 xlim = [minimum(ω), maximum(ω)],
                 xlabel="ω in s⁻¹",
                 ylabel="|λ|",
                 name="Euler Bernoulli Mode Number "*string(N[1])*" vs ω",
                 border=:ascii)
    if(length(N)>1)
        for n ∈ 2:length(N)
            i = N[n]
            lineplot!(plt, collect(ω), abs.(λ[i,:]),
                         name="Euler Bernoulli Mode Number "*string(i)*" vs ω")
        end
    end
    plt

end
