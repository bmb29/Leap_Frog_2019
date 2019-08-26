using Distributed
# @everywhere using DelimitedFiles
@everywhere using DifferentialEquations
@everywhere using ProgressMeter
@everywhere using MATLAB
@everywhere using Printf


include("escape_exit_function_parallel.jl")
Energy=.25
t_end=1e3
max_hit=5
@everywhere Hamil(XX,YY,QQ,PP)=( (QQ-XX)^2+(PP-YY)^2 )*( (QQ+XX)^2+(PP+YY)^2 )/((PP^4+2*PP^2*(XX^2-1)+(1+XX^2)^2 )*(QQ^4+2*QQ^2*(YY^2+1)+(YY^2-1)^2 ))
@everywhere ODE2(z,w)=conj(  im * z.*( 1 ./(w.^2-z.^2)+1 ./(1+z.^2) ))
@everywhere ODE1(z,w)=conj(  im * w.*( 1 ./(z.^2-w.^2)+1 ./(1+w.^2) ))
@everywhere function Eq_of_M(du,u,p,t)
    du[1]=real(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[2]=imag(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[3]=real(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[4]=imag(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[5]=0
    return
end
@everywhere condition(u,t,integrator)= u[5]>max_hit
@everywhere function condition2(u,t,integrator) # Event when event_f(u,t) == 0
    u[1]
end
@everywhere affect!(integrator) = terminate!(integrator)
@everywhere function affect2!(integrator)
    integrator.u[5]=integrator.u[5]+1
end
@everywhere cb2 = ContinuousCallback(condition2,affect2!,nothing)
@everywhere cb1 = DiscreteCallback(condition,affect!)

@everywhere cb=CallbackSet(cb2,cb1)
@everywhere H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))
@everywhere ODE1(z,w)=conj(  im * w.*( 1 ./(z.^2-w.^2)+1 ./(1+w.^2) ))
@everywhere ODE2(z,w)=conj(  im * z.*( 1 ./(w.^2-z.^2)+1 ./(1+z.^2) ))
@everywhere function Eq_of_M(du,u,t)
    du[1]=real(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[2]=imag(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[3]=real(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[4]=imag(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[5]=0
    return
end

#initialize lists
# Q_bound=zeros(0);P_bound=zeros(0)
Q_0=zeros(0);P_0=zeros(0)
Q_1=zeros(0);P_1=zeros(0)
Q_2=zeros(0);P_2=zeros(0)
Q_3=zeros(0);P_3=zeros(0)
Q_4=zeros(0);P_4=zeros(0)
Q_5=zeros(0);P_5=zeros(0)
width=2.5; height=1.5
n_iter_Q=20;n_iter_P=40;N=n_iter_P*n_iter_Q;

ArrP=range(-width,stop=width,length=n_iter_P)
ArrQ=range(-height,stop=height,length=n_iter_Q)
mesh=[(Q,P) for Q in ArrQ, P in ArrP]
@everywhere mesh_list=reshape(mesh,1,:)
@everywhere t_end=t_end*ones(N);
@everywhere Energy=Energy*ones(N)
# location="/mnt/bdd38f66-9ece-451a-b915-952523c139d2/Escape/"
# h=replace(@sprintf("%.13f",Energy),"."=>"_")
# file_name=location*"Escape_"*h*".fig"
file_name="test"


T=@showprogress pmap(escape_exit_function_parallel,mesh_list,t_end, Energy)
for i=1:N
    Q,P=mesh_list[i]
    if T[i]==0
        push!(Q_0,Q)
        push!(P_0,P)
    elseif T[i]==1
        push!(Q_1,Q)
        push!(P_1,P)
    elseif T[i]==2
        push!(Q_2,Q)
        push!(P_2,P)
    elseif T[i]==3
        push!(Q_3,Q)
        push!(P_3,P)
    elseif T[i]==4
        push!(Q_4,Q)
        push!(P_4,P)
    elseif T[i]==5
        push!(Q_5,Q)
        push!(P_5,P)
    end
end
mat"figure(); hold on;"
mat"axis([ -$width,$width,-$height,$height ])"
mat"plot($P_0,$Q_0 ,'b.','MarkerSize',3)"
mat"plot($P_1,$Q_1 ,'r.','MarkerSize',3)"
mat"plot($P_2,$Q_2 ,'g.','MarkerSize',3)"
mat"plot($P_4,$Q_4 ,'c.','MarkerSize',3)"
mat"plot($P_5,$Q_5 ,'y.','MarkerSize',3)"
mat"savefig($file_name)"
