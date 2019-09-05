
module escape_num
export escape_exit_num
using DifferentialEquations
using Roots
using Printf
# using MATLAB
# using ProgressMeter


include("leap_frog_definitions.jl")

max_hit=20
barrier=10
condition_max_hits(u,t,integrator)= u[5]>max_hit || maximum([abs(u[1]),abs(u[2]),abs(u[3]),abs(u[4])])>barrier
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

function escape_exit_num(mesh_list,t_end,Energy_A,max_hit)
    Q=mesh_list[2]; P=mesh_list[1];
    H=(2*Energy_A)^2
    tol_dist=1e-5
    Y=Yfind(Q,P,H)




    if ~isempty(Y)
        u0=zeros(5)
        u0[1]=0 #X
        u0[2]=P #P
        u0[3]=Q #Q
        u0[4]=Y #Y
        u0[5]=0 #Y

        prob = ODEProblem(Eq_of_M,u0,(0., t_end))
        # sol=solve(prob,RK4(),maxiters=1e20, reltol=1e-6,abstol=1e-8,callback=cb)
        sol=solve(prob, RK4(),maxiters=1e20, reltol=1e-6,abstol=1e-8,callback=cb,save_start=true,save_end=true,save_everystep=false)
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
            return sol.t[end]

        end
        # if dH<1e-5
        #     if is_it(uf,H,tol_dist)
        #         return sol.u[end][5]
        #     else
        #         return -1
        #     end
        # end
    else
        return 0.0
    end
end


end
