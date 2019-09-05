
module escape_num
export escape_exit_num
using DifferentialEquations
using Roots
using Printf
# using MATLAB
# using ProgressMeter

include("is_it.jl")

Hamil(XX,YY,QQ,PP)=( (QQ-XX)^2+(PP-YY)^2 )*( (QQ+XX)^2+(PP+YY)^2 )/((PP^4+2*PP^2*(XX^2-1)+(1+XX^2)^2 )*(QQ^4+2*QQ^2*(YY^2+1)+(YY^2-1)^2 ))
H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))
ODE2(z,w)=conj(  im * z.*( 1 ./(w.^2-z.^2)+1 ./(1+z.^2) ))
ODE1(z,w)=conj(  im * w.*( 1 ./(z.^2-w.^2)+1 ./(1+w.^2) ))

function Eq_of_M(du,u,p,t)
    du[1]=real(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[2]=imag(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[3]=real(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[4]=imag(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[5]=0
    return
end

function Yfind(Q,P,H)
    Y_find(y)=Hamil(0,y,Q,P)-H
    try
        Y=find_zeros(Y_find,0,10, maxeval=100,maxfnevals=300,tol=1e-15)
    catch
        Y=zeros(0)
    end
end



# function condition_escapes(u,t,p,integrator)
# end



function escape_exit_num(mesh_list,t_end,Energy_A,max_hit)
    Q=mesh_list[1]; P=mesh_list[2];
    H=(2*Energy_A)^2
    tol_dist=1e-5
    Y=Yfind(Q,P,H)
    p=zeros(2)
    barrier=5
    p[1]=max_hit
    p[2]=barrier

    condition_max_hits(u,t,integrator)= u[5]>p[1] || maximum([abs(u[1]),abs(u[2]),abs(u[3]),abs(u[4])])>p[2]
    affect_stop!(integrator) = terminate!(integrator)

    function condition_hits_PSS(u,t,integrator) # Event when event_f(u,t) == 0
       u[1]
    end

    function affect_update_iterator!(integrator)
        integrator.u[5]=integrator.u[5]+1
    end

    callback_max_hits=DiscreteCallback(condition_max_hits,affect_stop!)
    callback_hits_PSS=ContinuousCallback(condition_hits_PSS, affect_update_iterator!,nothing)
    cb=CallbackSet(callback_hits_PSS, callback_max_hits)

    if ~isempty(Y)
        u0=zeros(5)
        u0[1]=0 #X
        u0[2]=P #P
        u0[3]=Q #Q
        u0[4]=Y[1] #Y
        u0[5]=0 #Y
        
        prob = ODEProblem(Eq_of_M,u0,(0., t_end),p)
        # sol=solve(prob,RK4(),maxiters=1e20, reltol=1e-6,abstol=1e-8,callback=cb)
        sol=solve(prob, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb)
        uf=zeros(5)
        uf[1]=sol[1,end] #X
        uf[2]=sol[2,end] #P
        uf[3]=sol[3,end] #Q
        uf[4]=sol[4,end] #Y
        uf[5]=0
        dH=abs(H_test(uf)-H)
        if sol.u[end][5]==5
            return t_end
        else
            return
            sol.t[end]
        end
        # if dH<1e-5
        #     if is_it(uf,H,tol_dist)
        #         return sol.u[end][5]
        #     else
        #         return -1
        #     end
        # end
    end
    return 0
end


end
