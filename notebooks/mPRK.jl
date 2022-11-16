using LinearAlgebra

function mPRK22(prod_dest,rhs, tspan, u0::AbstractArray{T}, α) where {T}
    if α<1/2
        error("Negative RK coefficients")
    end
    c=T[0,α,1] #c,1
    A=T[[0 0]; [α 0]; [1-1/2/α 1/2/α]] #[A;b]
    ex=T[[0 0 0];[1 0 0];[1-1/α 1/α 0]]
    S = 3 # stages+1
    dim=length(u0)
    Nt=length(tspan)
    U=zeros(T, dim,Nt)
    y=zeros(T, S,dim)
    p=zeros(T, S,dim,dim)
    d=zeros(T, S,dim,dim)
    rhs_p =zeros(T, S,dim)
    dens=zeros(T, dim)
    U[:,1]=u0
    for it=2:Nt
        dt=tspan[it]-tspan[it-1]
        for k=1:S
            rhs_p[k,:] = U[:,it-1]
            for i=1:dim
                dens[i]=1
                for z=1:k-1
                    dens[i]=dens[i]*y[z,i]^ex[k,z]
                end
            end
            MM=Matrix{T}(I, dim, dim)
            for z=1:k-1
                rhs_p[k,:] = rhs_p[k,:] +dt*A[k,z]*rhs(y[z,:])
                for i=1:dim
                    for j=1:dim
                        MM[i,i]=MM[i,i]+A[k,z]*dt*d[z,i,j]/dens[i]
                        MM[i,j]=MM[i,j]-A[k,z]*dt*p[z,i,j]/dens[j]
                    end
                end
            end
            y[k,:]=MM\rhs_p[k,:]
            p[k,:,:],d[k,:,:]=prod_dest(y[k,:])
        end
        U[:,it]=y[S,:]
    end
    
    return tspan, U
end


function mPRK32(prod_dest,rhs, tspan, u0)
    c=[0,1,0.5,1] #c,1
    A=[[0 0 0]; [1 0 0]; [0.25 0.25 0]; [1/6 1/6 2/3]] #[A;b]
    S = 4 # stages+1
    dim=length(u0)
    Nt=length(tspan)
    U=zeros(dim,Nt)
    y=zeros(S,dim)
    p=zeros(S,dim,dim)
    d=zeros(S,dim,dim)
    rhs_p =zeros(S,dim)
    dens=zeros(dim)
    U[:,1]=u0
    for it=2:Nt
        dt=tspan[it]-tspan[it-1]
        for k=1:S
            rhs_p[k,:] = U[:,it-1]
            for i=1:dim
                if k==1
                    dens[i]=1
                else
                    dens[i]=y[minimum([k-1,2]),i]
                end
            end
            MM=Matrix{Float64}(I, dim, dim)
            for z=1:k-1
                rhs_p[k,:] = rhs_p[k,:] +dt*A[k,z]*rhs(y[z,:])
                for i=1:dim
                    for j=1:dim
                        MM[i,i]=MM[i,i]+A[k,z]*dt*d[z,i,j]/dens[i]
                        MM[i,j]=MM[i,j]-A[k,z]*dt*p[z,i,j]/dens[j]
                    end
                end
            end
            y[k,:]=MM\rhs_p[k,:]
            p[k,:,:],d[k,:,:]=prod_dest(y[k,:])
        end
        U[:,it]=y[S,:]
    end
    
    return tspan, U
end


function mPRKSO2(prod_dest,rhs,tspan, u0,α,β)
#    if α<0 || α>1 || β<0 || α*β+1/2/β>1 
#        error("Negative RK coefficients")
#    end
    c=[0,α,1] #c,1
    A=[[0 0]; [1 0]; [1-α α]] #[A;b]
    b = [[0 0]; [β 0]; [1-1/2/β-β*α 1/2/β]]
    s=(1-α*β+α*β^2)/β/(1-α*β)
    ex=[[0 0 0];[1 0 0];[1-s s 0]]
    S = 3 # stages+1
    dim=length(u0)
    Nt=length(tspan)
    U=zeros(dim,Nt)
    y=zeros(S,dim)
    p=zeros(S,dim,dim)
    d=zeros(S,dim,dim)
    rhs_p =zeros(S,dim)
    dens=zeros(dim)
    U[:,1]=u0
    for it=2:Nt
        dt=tspan[it]-tspan[it-1]
        y[1,:]=U[:,it-1]
        p[1,:,:],d[1,:,:]=prod_dest(y[1,:])
        rhs_p[1,:] = zeros(dim)
        for k=2:S
        rhs_p[k,:] = zeros(dim)
            for i=1:dim
                dens[i]=1
                for z=1:k-1
                    dens[i]=dens[i]*y[z,i]^ex[k,z]
                end
            end
            MM=Matrix{Float64}(I, dim, dim)
            for z=1:k-1
                rhs_p[k,:] = rhs_p[k,:] +dt*b[k,z]*rhs(y[z,:])
                for i=1:dim
                    for j=1:dim
                        if b[k,z]>=0
                            MM[i,i]=MM[i,i]+b[k,z]*dt*d[z,i,j]/dens[i]
                            MM[i,j]=MM[i,j]-b[k,z]*dt*p[z,i,j]/dens[j]
                        else
                            MM[i,i]=MM[i,i]-b[k,z]*dt*p[z,i,j]/dens[i]
                            MM[i,j]=MM[i,j]+b[k,z]*dt*d[z,i,j]/dens[j]
                        end
                    end
                end
            end
            RHS=A[:,:]*y[1:S-1,:]
            RHS[k,:] = RHS[k,:]+rhs_p[k,:]
            y[k,:]=MM\RHS[k,:]
            p[k,:,:],d[k,:,:]=prod_dest(y[k,:])
        end
        U[:,it]=y[S,:]
    end
    
    return tspan, U
end


function mPRK3(prod_dest,rhs,tspan, u0,α,β)
    α0=1/6*(3+(3-2*sqrt(2))^(1/3)+(3+2*sqrt(2))^(1/3))
#    if α<1/3  || β<0
#        error("Negative RK coefficients")
#    elseif α<2/3 && (β<2/3 || β >3*α*(1-α))
#        error("Negative RK coefficients")
#    elseif α>=2/3 && α<α0 && (β>2/3 || β <3*α*(1-α))
#        error("Negative RK coefficients")
#    elseif α>=α0 && (β>2/3 || β <(3*α-2)/(6*α-3))
#        error("Negative RK coefficients")
#    end
    
    c=[0,α,β, 1]
    A=[ [0 0 0]; [α 0 0]; [(3*α*β*(1-α)-β^2)/α/(2-3*α)  β*(β-α)/α/(2-3*α) 0]; [1+(2-3*(α+β))/6/α/β  (3*β-2)/(6*α*(β-α)) (2-3*α)/(6*β*(β-α))]]
    if maximum(A)<0 || maximum(c)<0
        error("Negative RK coefficients")
    end

    S = 4 # stages+1

    pp=3*A[2,1]*(A[3,1]+A[3,2])*A[S,3]
    q=A[2,1]
    beta=zeros(2)
    beta[2]=1/2/A[2,1]
    beta[1]=1-beta[2]
    #println("beta $(beta)")
    
    dim=length(u0)
    Nt=length(tspan)
    U=zeros(dim,Nt)
    y=zeros(S,dim)
    p=zeros(S,dim,dim)
    d=zeros(S,dim,dim)
    rhs_p =zeros(S,dim)
    dens=zeros(dim)
    U[:,1]=u0
    for it=2:Nt
        dt=tspan[it]-tspan[it-1]
        # u1
        y[1,:]=U[:,it-1]
        p[1,:,:],d[1,:,:]=prod_dest(y[1,:])
        rhs_p[1,:] = rhs(y[1,:])

        # u2
        k=2
        rhs_p[k,:] = zeros(dim)
        MM=Matrix{Float64}(I, dim, dim)
        for z=1:k-1
            rhs_p[k,:] = rhs_p[k,:] +dt*A[k,z]*rhs(y[z,:])
            for i=1:dim
                for j=1:dim
                    if A[k,z]>=0
                        MM[i,i]=MM[i,i]+A[k,z]*dt*d[z,i,j]/y[1,i]
                        MM[i,j]=MM[i,j]-A[k,z]*dt*p[z,i,j]/y[1,j]
                    else
                        MM[i,i]=MM[i,i]-A[k,z]*dt*p[z,i,j]/y[1,i]
                        MM[i,j]=MM[i,j]+A[k,z]*dt*d[z,i,j]/y[1,j]
                    end
                end
            end
        end
        RHS= y[1,:]+rhs_p[k,:]
        y[k,:]=MM\RHS
        p[k,:,:],d[k,:,:]=prod_dest(y[k,:])
        rhs_p[k,:] = rhs(y[k,:])

        # u3
        k=3
        rhs_p[k,:] = zeros(dim)
        rho=zeros(dim)

        for i=1:dim
            rho[i]=y[2,i]^(1/pp)*y[1,i]^(1-1/pp)
        end

        MM=Matrix{Float64}(I, dim, dim)
        for z=1:k-1
            rhs_p[k,:] = rhs_p[k,:] +dt*A[k,z]*rhs(y[z,:])
            for i=1:dim
                for j=1:dim
                    if A[k,z]>=0
                        MM[i,i]=MM[i,i]+A[k,z]*dt*d[z,i,j]/rho[i]
                        MM[i,j]=MM[i,j]-A[k,z]*dt*p[z,i,j]/rho[j]
                    else
                        MM[i,i]=MM[i,i]-A[k,z]*dt*p[z,i,j]/rho[i]
                        MM[i,j]=MM[i,j]+A[k,z]*dt*d[z,i,j]/rho[j]
                    end
                end
            end
        end
        RHS = y[1,:]+rhs_p[k,:]
        y[k,:]=MM\RHS
        p[k,:,:],d[k,:,:]=prod_dest(y[k,:])
        rhs_p[k,:] = rhs(y[k,:])

        # u4
        k=4
                        
        mu=zeros(dim)
        as = zeros(dim)

        for i=1:dim
            mu[i]=U[i,it-1].^(1-1/q).*y[2,i].^(1/q)
        end
        
        rrr=zeros(dim)
        MM=Matrix{Float64}(I, dim, dim)
        rrr = y[1,:] + dt*beta[1]*rhs(y[1,:])+dt*beta[2]*rhs(y[2,:])
        for i = 1:dim
            for j=1:dim
                for k=1:2
                    if beta[k]>=0
                        MM[i,j] = MM[i,j] -dt*(beta[k]*p[k,i,j])/mu[j]
                        MM[i,i] = MM[i,i] +dt*(beta[k]*d[k,i,j])/mu[i]
                    else
                        MM[i,j] = MM[i,j] +dt*(beta[k]*d[k,i,j])/mu[j]
                        MM[i,i] = MM[i,i] -dt*(beta[k]*p[k,i,j])/mu[i]
                    end
                end
            end
        end
        σ=MM\rrr

        rhs_p[k,:] = zeros(dim)
        MM=Matrix{Float64}(I, dim, dim)
        for z=1:k-1
            rhs_p[k,:] = rhs_p[k,:] +dt*A[k,z]*rhs(y[z,:])
            for i=1:dim
                for j=1:dim
                    if A[k,z]>=0
                        MM[i,i]=MM[i,i]+A[k,z]*dt*d[z,i,j]/σ[i]
                        MM[i,j]=MM[i,j]-A[k,z]*dt*p[z,i,j]/σ[j]
                    else
                        MM[i,i]=MM[i,i]-A[k,z]*dt*p[z,i,j]/σ[i]
                        MM[i,j]=MM[i,j]+A[k,z]*dt*d[z,i,j]/σ[j]
                    end
                end
            end
        end
        RHS = y[1,:]+rhs_p[k,:]
        y[k,:]=MM\RHS
        p[k,:,:],d[k,:,:]=prod_dest(y[k,:])
        U[:,it]=y[S,:]
    end
    return tspan, U
end




function mPRK3gam(prod_dest,rhs,tspan, u0, gamma)
#    if α<1/3  || β<0
#        error("Negative RK coefficients")
#    elseif α<2/3 && (β<2/3 || β >3*α*(1-α))
#        error("Negative RK coefficients")
#    elseif α>=2/3 && α<α0 && (β>2/3 || β <3*α*(1-α))
#        error("Negative RK coefficients")
#    elseif α>=α0 && (β>2/3 || β <(3*α-2)/(6*α-3))
#        error("Negative RK coefficients")
#    end
    
    c=[0,2.0/3.0,2.0/3.0, 1]
    A=[ [0 0 0]; [2.0/3.0 0 0]; [2.0/3.0-1.0/4.0/gamma  1.0/4.0/gamma 0]; [1.0/4.0  3.0/4.0-gamma gamma]]
    if maximum(A)<0 || maximum(c)<0
        println("Negative RK coefficients")
    end

    S = 4 # stages+1

    pp=3*A[2,1]*(A[3,1]+A[3,2])*A[S,3]
    q=A[2,1]
    beta=zeros(2)
    beta[2]=1/2/A[2,1]
    beta[1]=1-beta[2]
    #println("beta $(beta)")
    
    dim=length(u0)
    Nt=length(tspan)
    U=zeros(dim,Nt)
    y=zeros(S,dim)
    p=zeros(S,dim,dim)
    d=zeros(S,dim,dim)
    rhs_p =zeros(S,dim)
    dens=zeros(dim)
    U[:,1]=u0
    for it=2:Nt
        dt=tspan[it]-tspan[it-1]
        # u1
        y[1,:]=U[:,it-1]
        p[1,:,:],d[1,:,:]=prod_dest(y[1,:])
        rhs_p[1,:] = rhs(y[1,:])

        # u2
        k=2
        rhs_p[k,:] = zeros(dim)
        MM=Matrix{Float64}(I, dim, dim)
        for z=1:k-1
            rhs_p[k,:] = rhs_p[k,:] +dt*A[k,z]*rhs(y[z,:])
            for i=1:dim
                for j=1:dim
                    if A[k,z]>=0
                        MM[i,i]=MM[i,i]+A[k,z]*dt*d[z,i,j]/y[1,i]
                        MM[i,j]=MM[i,j]-A[k,z]*dt*p[z,i,j]/y[1,j]
                    else
                        MM[i,i]=MM[i,i]-A[k,z]*dt*p[z,i,j]/y[1,i]
                        MM[i,j]=MM[i,j]+A[k,z]*dt*d[z,i,j]/y[1,j]
                    end
                end
            end
        end
        RHS= y[1,:]+rhs_p[k,:]
        y[k,:]=MM\RHS
        p[k,:,:],d[k,:,:]=prod_dest(y[k,:])
        rhs_p[k,:] = rhs(y[k,:])

        # u3
        k=3
        rhs_p[k,:] = zeros(dim)
        rho=zeros(dim)

        for i=1:dim
            rho[i]=y[2,i]^(1/pp)*y[1,i]^(1-1/pp)
        end

        MM=Matrix{Float64}(I, dim, dim)
        for z=1:k-1
            rhs_p[k,:] = rhs_p[k,:] +dt*A[k,z]*rhs(y[z,:])
            for i=1:dim
                for j=1:dim
                    if A[k,z]>=0
                        MM[i,i]=MM[i,i]+A[k,z]*dt*d[z,i,j]/rho[i]
                        MM[i,j]=MM[i,j]-A[k,z]*dt*p[z,i,j]/rho[j]
                    else
                        MM[i,i]=MM[i,i]-A[k,z]*dt*p[z,i,j]/rho[i]
                        MM[i,j]=MM[i,j]+A[k,z]*dt*d[z,i,j]/rho[j]
                    end
                end
            end
        end
        RHS = y[1,:]+rhs_p[k,:]
        y[k,:]=MM\RHS
        p[k,:,:],d[k,:,:]=prod_dest(y[k,:])
        rhs_p[k,:] = rhs(y[k,:])

        # u4
        k=4
                        
        mu=zeros(dim)
        as = zeros(dim)

        for i=1:dim
            mu[i]=U[i,it-1].^(1-1/q).*y[2,i].^(1/q)
        end
        
        rrr=zeros(dim)
        MM=Matrix{Float64}(I, dim, dim)
        rrr = y[1,:] + dt*beta[1]*rhs(y[1,:])+dt*beta[2]*rhs(y[2,:])
        for i = 1:dim
            for j=1:dim
                for k=1:2
                    if beta[k]>=0
                        MM[i,j] = MM[i,j] -dt*(beta[k]*p[k,i,j])/mu[j]
                        MM[i,i] = MM[i,i] +dt*(beta[k]*d[k,i,j])/mu[i]
                    else
                        MM[i,j] = MM[i,j] +dt*(beta[k]*d[k,i,j])/mu[j]
                        MM[i,i] = MM[i,i] -dt*(beta[k]*p[k,i,j])/mu[i]
                    end
                end
            end
        end
        σ=MM\rrr

        rhs_p[k,:] = zeros(dim)
        MM=Matrix{Float64}(I, dim, dim)
        for z=1:k-1
            rhs_p[k,:] = rhs_p[k,:] +dt*A[k,z]*rhs(y[z,:])
            for i=1:dim
                for j=1:dim
                    if A[k,z]>=0
                        MM[i,i]=MM[i,i]+A[k,z]*dt*d[z,i,j]/σ[i]
                        MM[i,j]=MM[i,j]-A[k,z]*dt*p[z,i,j]/σ[j]
                    else
                        MM[i,i]=MM[i,i]-A[k,z]*dt*p[z,i,j]/σ[i]
                        MM[i,j]=MM[i,j]+A[k,z]*dt*d[z,i,j]/σ[j]
                    end
                end
            end
        end
        RHS = y[1,:]+rhs_p[k,:]
        y[k,:]=MM\RHS
        p[k,:,:],d[k,:,:]=prod_dest(y[k,:])
        U[:,it]=y[S,:]
    end
    return tspan, U
end


function mPRKSO3(prod_dest,rhs,tspan, u0)
    A=[[0 0 0]; [1 0 0]; [9.2600312554031827E-01 7.3996874459681783E-02 0]; 
       [7.0439040373427619E-01 2.0662904223744017E-10 2.9560959605909481E-01]] #[A;b]
    b = [[0 0 0 0]; [4.7620819268131703E-01 0 0 0]; [7.7545442722396801E-02 5.9197500149679749E-01 0 0];
         [2.0044747790361456E-01 6.8214380786704851E-10 5.9121918658514827E-01 0]]

    n1=2.569046025732011E-01; n2 = 7.430953974267989E-01; z = 6.288938077828750E-01;	
    η1 = 3.777285888379173E-02; η2 = 1/3; η3 = 1.868649805549811E-01; η4=2.224876040351123; s = 5.721964308755304

    S = 4 # stages+1
    dim=length(u0)
    Nt=length(tspan)
    U=zeros(dim,Nt)
    y=zeros(S,dim)
    p=zeros(S,dim,dim)
    d=zeros(S,dim,dim)
    rhs_p =zeros(S,dim)
    dens=zeros(dim)
    U[:,1]=u0
    for it=2:Nt
        dt=tspan[it]-tspan[it-1]
        # u0
        y[1,:]=U[:,it-1]
        p[1,:,:],d[1,:,:]=prod_dest(y[1,:])
        rhs_p[1,:] = rhs(y[1,:])

        # u1
        k=2
        rhs_p[k,:] = zeros(dim)
        MM=Matrix{Float64}(I, dim, dim)
        for z=1:k-1
            rhs_p[k,:] = rhs_p[k,:] +dt*b[k,z]*rhs(y[z,:])
            for i=1:dim
                for j=1:dim
                    MM[i,i]=MM[i,i]+b[k,z]*dt*d[z,i,j]/y[1,i]
                    MM[i,j]=MM[i,j]-b[k,z]*dt*p[z,i,j]/y[1,j]
                end
            end
        end
        RHS=A[:,:]*y[1:S-1,:]
        RHS[k,:] = RHS[k,:]+rhs_p[k,:]
        y[k,:]=MM\RHS[k,:]
        p[k,:,:],d[k,:,:]=prod_dest(y[k,:])
        rhs_p[k,:] = rhs(y[k,:])

        # u2
        k=3
        rhs_p[k,:] = zeros(dim)
        rho=zeros(dim)

        for i=1:dim
            rho[i]=n1*y[2,i]+n2*y[1,i]*(y[2,i]/y[1,i])^2
        end
    
        MM=Matrix{Float64}(I, dim, dim)
        for z=1:k-1
            rhs_p[k,:] = rhs_p[k,:] +dt*b[k,z]*rhs(y[z,:])
            for i=1:dim
                for j=1:dim
                    MM[i,i]=MM[i,i]+b[k,z]*dt*d[z,i,j]/rho[i]
                    MM[i,j]=MM[i,j]-b[k,z]*dt*p[z,i,j]/rho[j]
                end
            end
        end
        RHS=A[:,:]*y[1:S-1,:]
        RHS[k,:] = RHS[k,:]+rhs_p[k,:]
        y[k,:]=MM\RHS[k,:]
        p[k,:,:],d[k,:,:]=prod_dest(y[k,:])
        rhs_p[k,:] = rhs(y[k,:])

        # u3
        k=4
        rhs_p[k,:] = zeros(dim)
        mu=zeros(dim)
        as = zeros(dim)

        for i=1:dim
            mu[i]=y[1,i]*(y[2,i]/y[1,i])^s
        end

        rrr=zeros(dim)
        MM=Matrix{Float64}(I, dim, dim)
        for i = 1:dim
            rrr[i] = η1*y[1,i]+η2*y[2,i]+dt*η3*(η1+η2)*rhs_p[1,i]+dt*η4*(η1+η2)*rhs_p[2,i]
            for j=1:dim
                MM[i,j] = MM[i,j] -dt*(η3*p[1,i,j]+η4*p[2,i,j])/mu[j]
                MM[i,i] = MM[i,i] +dt*(η3*d[1,i,j]+η4*d[2,i,j])/mu[i]
            end
        end
        as=MM\rrr

        σ=zeros(dim)
        for i=1:dim
            σ[i]=as[i]+z*y[1,i]*y[3,i]/rho[i]
        end

        MM=Matrix{Float64}(I, dim, dim)
        for z=1:k-1
            rhs_p[k,:] = rhs_p[k,:] +dt*b[k,z]*rhs(y[z,:])
            for i=1:dim
                for j=1:dim
                    MM[i,i]=MM[i,i]+b[k,z]*dt*d[z,i,j]/σ[i]
                    MM[i,j]=MM[i,j]-b[k,z]*dt*p[z,i,j]/σ[j]
                end
            end
        end
        RHS=A[:,:]*y[1:S-1,:]
        RHS[k,:] = RHS[k,:]+rhs_p[k,:]
        y[k,:]=MM\RHS[k,:]
        p[k,:,:],d[k,:,:]=prod_dest(y[k,:])
        U[:,it]=y[S,:]
    end
    return tspan, U
end



function SIRK2(f,g,tspan, u0)
    S = 4 # stages+1

    dim=length(u0)
    Nt=length(tspan)
    U=zeros(dim,Nt)
    y=zeros(S,dim)
    ff=zeros(S,dim)
    gg=zeros(S,dim)
    U[:,1]=u0
    for it=2:Nt
        dt=tspan[it]-tspan[it-1]
        # u0
        y[1,:]=U[:,it-1]
        ff[1,:]=f(y[1,:])
        gg[1,:]=g(y[1,:])

        # u1
        k=2
        y[k,:] = (U[:,it-1] .+ dt*ff[1,:])./ (1 .-dt.*gg[1,:])
        ff[k,:]=f(y[k,:])
        gg[k,:]=g(y[k,:])
        
        # u2
        k = 3
        y[k,:] = 0.5*U[:,it-1] .+ 0.5*(y[2,:] .+ dt*ff[2,:])./ (1 .-dt.*gg[2,:])
        ff[k,:]=f(y[k,:])
        gg[k,:]=g(y[k,:])
        
        # u^{n+1}
        k = 4
        y[k,:] = (y[3,:] .- dt^2*ff[3,:].*gg[3,:])./ (1 .+dt^2 .*gg[3,:].*gg[3,:])
        
        U[:,it]=y[S,:]
    end
    return tspan, U
end




function SIRK3(f,g,tspan, u0)
    S = 5 # stages+1

    dim=length(u0)
    Nt=length(tspan)
    U=zeros(dim,Nt)
    y=zeros(S,dim)
    ff=zeros(S,dim)
    gg=zeros(S,dim)
    U[:,1]=u0
    for it=2:Nt
        dt=tspan[it]-tspan[it-1]
        # u0
        y[1,:]=U[:,it-1]
        ff[1,:]=f(y[1,:])
        gg[1,:]=g(y[1,:])

        # u1
        k=2
        y[k,:] = (U[:,it-1] .+ dt*ff[1,:])./ (1 .-dt.*gg[1,:])
        ff[k,:]=f(y[k,:])
        gg[k,:]=g(y[k,:])
        
        # u2
        k = 3
        y[k,:] = 3.0/4.0*U[:,it-1].+ 1.0/4.0*(y[2,:] .+ dt*ff[2,:])./ (1 .-dt.*gg[2,:])
        ff[k,:]=f(y[k,:])
        gg[k,:]=g(y[k,:])
        
        # u3
        k = 4
        y[k,:]=1.0/3.0*U[:,it-1].+2.0/3.0.*(y[3,:] .+ dt*ff[3,:])./(1 .-dt.*gg[3,:])
        ff[k,:]=f(y[k,:])
        gg[k,:]=g(y[k,:])
        
        # u^{n+1}
        k = 5
        y[k,:] = (y[4,:] .- dt^2*ff[4,:].*gg[4,:])./ (1 .+dt^2 .*gg[4,:].*gg[4,:])
        
        U[:,it]=y[S,:]
    end
    return tspan, U
end



