
module escape_num_dimer
export escape_exit_num_dimer
using DifferentialEquations
using Roots
using Printf
# using MATLAB
# using ProgressMeter


include("leap_frog_definitions.jl")
max_hit=500
barrier=10


condition_max_hits(u,t,integrator)= u[6]>max_hit || maximum([abs(u[1]),abs(u[2]),abs(u[4]),abs(u[5])])>barrier
affect_stop!(integrator) = terminate!(integrator)
function condition_hits_PSS(u,t,integrator) # Event when event_f(u,t) == 0
   u[1]
end
function affect_update_iterator!(integrator)
    integrator.u[6]=integrator.u[6]+1
end


callback_max_hits=DiscreteCallback(condition_max_hits,affect_stop!)
callback_hits_PSS=ContinuousCallback(condition_hits_PSS, affect_update_iterator!,nothing)
cb=CallbackSet(callback_hits_PSS, callback_max_hits)

function escape_exit_num_dimer(mesh_list,t_end,H)
    Q2=mesh_list[2]; P2=mesh_list[1];
    P1=P1_find_dimer(Q2,P2,H)
    if isempty(P1)
        P1=P1_find_dimer_second(Q2,P2,H)
        # t_end=1
    end
    if ~isempty(P1)
        # println(Q1)
        q0,p0=[zeros(3) for i in 1:2]
        q0[1]=0
        q0[2]=Q2
        q0[3]=0
        p0[1]=P1[1]
        p0[2]=P2
        p0[3]=0
       #constructor for ODE
        prob= HamiltonianProblem{true}(Hamiltonian_Dimer, q0, p0, (0., t_end));
        #solve ode , save_everystep=false is important to prevent sol to include all points, not just event points
        sol=solve(prob, RK4(),maxiters=1e20, reltol=1e-8,callback=cb,save_start=true,save_end=true,save_everystep=false)
        if sol.u[end][6]==max_hit
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
