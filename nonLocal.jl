# Function to get the associated matrices on the semi--infinite domain
function innerproduct(k, kappa, H, d)
    if(abs(k-kappa)>=1e-7)
        return ( (kappa*sin(kappa*(H-d))*cos(k*(H-d)) - k*cos(kappa*(H-d))*sin(k*(H-d)))/(kappa^2-k^2) );
    else
        return ( (H-d)/2 + sin(2*k*(H-d))/(4*k) );
    end
end

# Get Matrices in the open ocean.
function getMAT(k, kd, H, d, NModes, Ap)
    A = zeros(Complex{Float64},NModes+1,NModes+1);
    M = zeros(Complex{Float64},NModes+1,NModes+1);
    f = zeros(Complex{Float64},NModes+1,1);
    g = zeros(Complex{Float64},NModes+1,1);

    for i=1:NModes+1
        A[i,i]=0.5*(cos(k[i]*H)*sin(k[i]*H) + k[i]*H)/(cos(k[i]*H))^2
        f[i]=Ap*innerproduct(k[1], kd[i], H, d)/(cos(kd[i]*(H-d))*cos(k[1]*H));
        for j=1:NModes+1
            M[i,j]=innerproduct(k[j], kd[i], H, d)/(cos(kd[i]*(H-d))*cos(k[j]*H));
        end
    end

    g[1]=-Ap*A[1,1];
    return A, M, f, g;
end

# Get the non-local matrix on the boundary
function getMQχ(k, kd, H, d, NModes, Ap, model, Γ, V, V0)
    A,M,f,g=getMAT(k,kd,H,d,NModes,Ap);
    ndofs=num_free_dofs(V); #Get the number of dofs in the domain.
    pp=zeros(ComplexF64,NModes+1,ndofs);
    for m=1:NModes+1
        τ(x)=cos(kd[m]*(x[2]+H))/cos(kd[m]*(H-d));
        dΓ=Measure(Γ,2);
        a(u,v)=∫(u*v)*dΓ; #Dummy bilinear form ?
        b(v)=∫(τ*v)*dΓ;
        op=AffineFEOperator(a,b,V,V0);
        pp[m,:]=op.op.vector;
    end
    # Get the matrix corresponding to Qϕ
    Mt=transpose(M);
    T=inv(Mt)*A*inv(M);
    tP=T*pp;
    Qϕ=transpose(pp)*tP;

    c=inv(Mt)*g-T*f;
    χ=transpose(c)*pp;

    return Qϕ,χ;
end
