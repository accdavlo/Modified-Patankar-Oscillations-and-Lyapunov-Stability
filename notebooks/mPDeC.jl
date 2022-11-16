using FastGaussQuadrature, Compat
using LinearAlgebra

"""Different models"""

#include("Linear_model.jl")
#include("Robertson.jl")
#include("Nonlinear_model.jl")
#include("Nonlinear_model2.jl")

#include("Lagrange_polynomials.jl")
#include("Error_module.jl")


function lagrange_basis(X,t,i)
    idxs = eachindex(X)
    if t in X
        if t == X[i]
            return 1
        else
            return 0
        end
    end
    return prod((t-X[j])/(X[i]-X[j]) for j in idxs if j != i)
end


function get_nodes(order,nodes_type)
    if nodes_type=="equispaced"
        nodes= range(0,stop=1, length=order)
        w = 1/order*ones(order)
    elseif nodes_type=="gaussLobatto"
        nodes, w = gausslobatto(order)
        w = map( x-> x*0.5, w)
        nodes = nodes*0.5.+0.5
    elseif nodes_type=="gaussLegendre"
        nodes, w = gausslegendre(order)
        w = map( x-> x*0.5, w)
        nodes = nodes*0.5.+0.5
    end
    
    return nodes, w
end


""" Computing the matrix of the values T_k^m = int_{t^0}^{t^m} varphi_k"""
function compute_theta_DeC(order, nodes_type)
    
    nodes, w = get_nodes(order,nodes_type)
    int_nodes, int_w = get_nodes(order,"gaussLobatto")
    # generate theta coefficients 
    theta = zeros(order,order)
    for m in 1:order
        nodes_m = int_nodes*(nodes[m])
        w_m = int_w*(nodes[m])
        for r in 1:order
            theta[r,m] = sum([lagrange_basis(nodes,nodes_m[k],r)*w_m[k] for k in 1:order])
        end
    end
    return theta
    
end


#Dec_correction-
#Func is the function of the ode
#N_step -> Number of subtimesteps
#Corr-> Number of Corrections
#y_0 Initial condition
# t_0 Start value
# t_end End value
# Maybe logical variable for checking for the different solvers or the distributaiton
# of subtimesteps up to some order

""" Dec_correction constant time step"""
function dec_correction( prod_dest, rhs, tspan, y_0, M_sub::Int, K_corr::Int,
    distribution)

   N_time=length(tspan)
   dim=length(y_0)
   U=zeros(dim, N_time)
   u_p=zeros(dim, M_sub+1)
   u_a=zeros(dim, M_sub+1)
   prod_p=zeros(dim,dim,M_sub+1)
   dest_p=zeros(dim,dim,M_sub+1)
   rhs_p =zeros(dim,M_sub+1)
   Theta=compute_theta_DeC(M_sub+1,distribution)
   U[:,1]=y_0
   for it=2: N_time
        delta_t=(tspan[it]-tspan[it-1])
        for m=1:M_sub+1
         u_a[:,m]=U[:,it-1]
         u_p[:,m]=U[:,it-1]
        end
        for k=2:K_corr+1
            u_p=copy(u_a)
            for r=1:M_sub+1
                prod_p[:,:,r], dest_p[:,:,r]=prod_dest(u_p[:,r])
                rhs_p[:,r] = rhs(u_p[:,r])
            end
            for m=2:M_sub+1
            u_a[:,m]=patanker_type_dec(prod_p, dest_p, rhs_p, delta_t, m, M_sub, Theta, u_p, dim)
            end
        end
        U[:,it]=u_a[:,M_sub+1]
    end
    return tspan, U
end



"""prod_p and dest_p dim x dim
delta_t is timestep length
m_substep is the subtimesteps for which we are solving the system
M_sub is the maximum number of subtimesteps (0,..., M)
theta coefficient vector matrix  in  M_sub+1 x M_sub+1, first index is the index we are looping on,
the second is the one we are solving for, theta_r,m= int_t0^tm phi_r
u_p is the previous correction solution for all subtimesteps in dim x (M_sub+1) """
function patanker_type_dec(prod_p, dest_p, rhs_p, delta_t,  m_substep::Int, M_sub, theta, u_p, dim)

mass=Matrix{Float64}(I, dim, dim);
#println(mass)
RHS = u_p[:,1]
for i=1:dim
    for r=1: M_sub+1
				RHS[i]=RHS[i] +delta_t*theta[r,m_substep]*rhs_p[i,r]
        if theta[r,m_substep]>0
            for j= 1: dim
              mass[i,j]=mass[i,j]- delta_t*theta[r,m_substep]*(prod_p[i,j,r]/u_p[j,m_substep])
              mass[i,i]=mass[i,i]+ delta_t*theta[r,m_substep]*(dest_p[i,j,r]/u_p[i,m_substep])
            end
        elseif theta[r,m_substep]<0
            for j= 1: dim
              mass[i,i]=mass[i,i]- delta_t*theta[r,m_substep]*(prod_p[i,j,r]/u_p[i,m_substep])
              mass[i,j]=mass[i,j]+ delta_t*theta[r,m_substep]*(dest_p[i,j,r]/u_p[j,m_substep])
            end
        end
    end
end
return mass\RHS
end

##dim=2;
##M_sub=4;
##u_p=rand(dim, M_sub+1 );
##prod_p=dest_p= zeros(dim,dim,M_sub+1)
##for r=1:M_sub+1
##    prod_p[:,:,r],dest_p[:,:,r]=prod_dest(u_p[:,1])
##end

##delta_t=0.1
##m_substep=2
##theta=-1.0 .+2.0.*rand(M_sub+1,M_sub+1)

##patanker_type_dec(prod_p,dest_p, delta_t, m_substep, M_sub, theta, u_p, dim)


#####Convergence test
###u0=[0.9;0.1]
###dim = length(u0)
###T_final = 1.75
###n_ns = 5
###Ns = 2 .^range(3,stop=2+n_ns)
###dts = T_final./Ns
###U_exact = exact_solution(u0,T_final)
###plot()
###for order = 2:6
###  Us = zeros(dim, n_ns)
###  for k=1:n_ns
###     t,U=dec_correction(order,order+1,u0,0,T_final,Ns[k],true)
###     Us[:,k]=U[:,end]
###  end


###  plot_convergence(dts,Us,U_exact,"dec $(order)", order)
###end
###plot!()



"""
Later Version for different quadrature rules, for later worth it to check!

"""

function quadrature_distribution_name_to_type(name::Symbol)
    if name == :Periodic
        #Going in the slope
    elseif name == :Lobatto
        gausslobatto
    elseif name == :Chebyshev
        gausschebyshev
    else
        error("Quadrature `$name` unknown.")
    end
end
