
module escape_num
export escape_exit_num
using DifferentialEquations
using Roots
using Printf
# using MATLAB
# using ProgressMeter

include("is_it.jl")
include("leap_frog_definitions.jl")

condition_max_hits(u,t,p,integrator)= u[5]>p
affect_stop!(integrator) = terminate!(integrator)

function condition_hits_PSS(u,t,integrator) # Event when event_f(u,t) == 0
    u[1]
end
function affect_update_iterator!(integrator)
    integrator.u[5]=integrator.u[5]+1
end

callback_max_hits=DiscreteCallback(condition_max_hits,affect_stop!)
callback_hits_PSS_=ContinuousCallback(condition_hits_PSS, affect_update_iterator!,nothing)
cb=CallbackSet(callback_hits_PSS, callback_max_hits)


function escape_exit_num(mesh_list,t_end,Energy_A,max_hit)
    Q=mesh_list[1]; P=mesh_list[2];
    H=(2*Energy_A)^2
    tol_dist=1e-5
    Y=Yfind(Q,P,H)
    p=max_hit
    if ~isempty(Y)
        u0=zeros(5)
        u0[1]=0 #X
        u0[2]=P #P
        u0[3]=Q #Q
        u0[4]=Y[1] #Y
        u0[5]=0 #Y
        
        prob = ODEProblem(Eq_of_M,u0,(0., t_end),max_hit,p)
        # sol=solve(prob,RK4(),maxiters=1e20, reltol=1e-6,abstol=1e-8,callback=cb)
        sol=solve(prob, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb)
        uf=zeros(5)
        uf[1]=sol[1,end] #X
        uf[2]=sol[2,end] #P
        uf[3]=sol[3,end] #Q
        uf[4]=sol[4,end] #Y
        uf[5]=0
        dH=abs(H_test(uf)-H)
        if dH<1e-5
            if is_it(uf,H,tol_dist)
                return sol.u[end][5]
            else
                return -1
            end
        end
    end
    return -2
end


end
