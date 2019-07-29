using DelimitedFiles
using DifferentialEquations
using ProgressMeter
include("PSS_vern9.jl")
 using PyCall
pygui(:qt)
using PyPlot
pygui(true)
max_hit=2
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


Energy=.25
t_end=100.
n_iter_P=500
n_iter_Q=501

# ArrP=range(1.1,stop=2,length=n_iter_P)
# ArrQ=range(1e-10,stop=1,length=n_iter_Q)
# ArrP=range(0,stop=.1,length=n_iter_P)
# ArrQ=range(1e-10,stop=.1,length=n_iter_Q)
ArrP=zeros(n_iter_P)
ArrQ=zeros(n_iter_Q)
Q_bound=zeros(0)
P_bound=zeros(0)

Q_0=zeros(0)
P_0=zeros(0)

Q_1=zeros(0)
P_1=zeros(0)

Q_2=zeros(0)
P_2=zeros(0)

Q_3=zeros(0)
P_3=zeros(0)

Q_4=zeros(0)
P_4=zeros(0)

Q_5=zeros(0)
P_5=zeros(0)

@showprogress 1 "Computing..." for j=1:n_iter_P
    for k=1:n_iter_Q
        # ArrQ[k]=.1*rand()
        # ArrP[j]=.01*rand()
        K=PSS_vern9(ArrQ[k],ArrP[j], Energy, t_end)
        if K==0
            push!(Q_0,ArrQ[k])
            push!(P_0,ArrP[j])
        elseif K==1
            push!(Q_1,ArrQ[k])
            push!(P_1,ArrP[j])
        elseif K==2
            push!(Q_2,ArrQ[k])
            push!(P_2,ArrP[j])
        elseif K==3
            push!(Q_3,ArrQ[k])
            push!(P_3,ArrP[j])
        elseif K==4
            push!(Q_4,ArrQ[k])
            push!(P_4,ArrP[j])

        elseif K==-1
            push!(Q_bound,ArrQ[k])
            push!(P_bound,ArrP[j])
        end
    end
end
