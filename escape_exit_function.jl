
include("YfindNP.jl")
include("is_it_exit.jl")
include("is_it.jl")


Hamil(XX,YY,QQ,PP)=( (QQ-XX)^2+(PP-YY)^2 )*( (QQ+XX)^2+(PP+YY)^2 )/((PP^4+2*PP^2*(XX^2-1)+(1+XX^2)^2 )*(QQ^4+2*QQ^2*(YY^2+1)+(YY^2-1)^2 ))
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
condition(u,t,integrator)= u[5]>max_hit
function condition2(u,t,integrator) # Event when event_f(u,t) == 0
u[1]
end
affect!(integrator) = terminate!(integrator)
function affect2!(integrator)
integrator.u[5]=integrator.u[5]+1
end
cb2 = ContinuousCallback(condition2,affect2!,nothing)
cb1 = DiscreteCallback(condition,affect!)

cb=CallbackSet(cb2,cb1)
H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))


function escape_exit_function(Q,P, Energy, t_end, max_hit)

    tol_dist=1e-5
    H=(2*Energy)^2

    Q_0=0
    P_0=0

    Q_1=0
    P_1=0

    Q_2=0
    P_2=0

    Q_3=0
    P_3=0


    Y=YfindNP(Q,P,H)
    if ~isempty(Y)
        u0=zeros(5)
        u0[1]=0 #X
        u0[2]=P #P
        u0[3]=Q #Q
        u0[4]=Y[1] #Y
        u0[5]=0 #Y
        # counter variable
        #are the initial conditions right?
        # prob_test= ODEProblem(Eq_of_M,u0,(0., 10.))
        # sol_test=solve(prob_test, RK4(),maxiters=1e10 )
        # S=sign(sol_test.u[2][1])
        if true
            prob = ODEProblem(Eq_of_M,u0,(0., t_end))
            # sol=solve(prob,RK4(),maxiters=1e20, reltol=1e-6,abstol=1e-8,callback=cb)
            sol=solve(prob, Vern9(),maxiters=1e10, reltol=1e-8,abstol=1e-10,callback=cb)
            uf=zeros(5)
            uf[1]=sol[1,end] #X
            uf[2]=sol[2,end] #P
            uf[3]=sol[3,end] #Q
            uf[4]=sol[4,end] #Y
            uf[5]=0
            dH=abs(H_test(uf)-H)
            #
            if dH<1e-5
                B=hcat(sol.u[end-10],sol.u[end])'
                if is_it_exit(sol.u[end],H,tol_dist)
                    return sol.u[end][5]
                else
                    return -1
                end
            end
        end
    end
    return -2
end
