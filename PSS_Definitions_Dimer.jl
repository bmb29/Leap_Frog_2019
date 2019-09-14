include("leap_frog_definitions.jl")

function PSS_function(Q2,P2, H,  t_end)
    #for a given Q,P,H with X=0
    Q1=Q1_find_dimer(Q2,P2,H)
    if ~isempty(Q2)
        println(Q2)
        q0,p0=[zeros(2) for i in 1:2]
        q0[1]=Q1
        q0[2]=Q2
        p0[1]=0
        p0[2]=P2
        #constructor for ODE
        prob= HamiltonianProblem{true}(Hamiltonian_Dimer, q0, p0, (0., t_end));
        #solve ode , save_everystep=false is important to prevent sol to include all points, not just event points
        # sol=solve(prob, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb,save_start=true,save_end=true,save_everystep=false)
        sol=solve(prob, RK4(),maxiters=1e20, reltol=1e-6,abstol=1e-6,callback=cb,save_start=true,save_end=true,save_everystep=false)

        #output Q, P and dH
        return sol[:,2:end-1][2,:],sol[:,2:end-1][4,:]
    else
    #need to return 3 values, dH=1 flags that there is no Y
        return 0, 0
    end
end
