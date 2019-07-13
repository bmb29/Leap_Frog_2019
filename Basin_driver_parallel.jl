using Distributed
@everywhere using DelimitedFiles
@everywhere using DifferentialEquations
@everywhere using ProgressMeter

include("BASIN_PSS_function_parallel.jl")
@everywhere using PyCall
pygui(:qt)
@everywhere using PyPlot
pygui(true)

@everywhere max_hit=30
@everywhere Energy=.25
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
@everywhere cb2 = ContinuousCallback(condition2,affect2!,rootfind=false)
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

# ArrQ[k]=.1*rand()
# ArrP[j]=.01*rand()
@everywhere Q_start=-1.5
@everywhere Q_end=1.5
@everywhere P_start=-2.5
@everywhere P_end=2.5
@everywhere n_iter_P=2000
@everywhere n_iter_Q=2001
@everywhere N=n_iter_P*n_iter_Q

@everywhere delta_P=(P_end-P_start)/n_iter_P
@everywhere delta_Q=(Q_end-Q_start)/n_iter_Q

@everywhere ArrP=P_start: delta_P:P_end
@everywhere ArrQ=Q_start: delta_Q: Q_end

T=@showprogress pmap(BASIN_PSS_function_parallel,1:N)
for i=1:N
    Q=ArrQ[Int(ceil(i/n_iter_Q))]
    P=ArrP[mod(i,n_iter_P)+1]
    if floor(T[i])==0
        push!(Q_0,Q)
        push!(P_0,P)
    elseif floor(T[i])==1
        push!(Q_1,Q)
        push!(P_1,P)
    elseif floor(T[i])==2
        push!(Q_2,Q)
        push!(P_2,P)
    elseif floor(T[i])==3
        push!(Q_3,Q)
        push!(P_3,P)
    elseif floor(T[i])==4
        push!(Q_4,Q)
        push!(P_4,P)
    elseif floor(T[i])==5
        push!(Q_5,Q)
        push!(P_5,P)
    elseif T[i]==-1

        push!(Q_bound,Q)
        push!(P_bound,P)
    end
end

figure()
# axis([ -.6, .6,-1.1, 1.1])
axis([ -2.5, 2.5,-1.5, 1.5])

plot(P_bound,Q_bound ,color="k",".",markersize=2, markeredgewidth=.1)

T0=hcat(Q_0,P_0)
plot(P_0,Q_0 ,color="r",".",markersize=2, markeredgewidth=.1)

T1=hcat(Q_1,P_1)
plot(P_1,Q_1 ,color="b",".",markersize=2, markeredgewidth=.1)

T2=hcat(Q_2,P_2)
plot(P_2,Q_2 ,color="g",".",markersize=2, markeredgewidth=.1)

T3=hcat(Q_3,P_3)
plot(P_3,Q_3 ,color="m",".",markersize=2, markeredgewidth=.1)

T4=hcat(Q_4,P_4)
plot(P_4,Q_4 ,color="c",".",markersize=2, markeredgewidth=.1)

T5=hcat(Q_5,P_5)
plot(P_5,Q_5 ,color="y",".",markersize=2, markeredgewidth=.1)
#
# outfile_0 = "ParVern9_outfile0_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
# outfile_1 = "ParVern9_outfile1_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
# outfile_2 = "ParVern9_outfile2_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
# outfile_3 = "ParVern9_outfile3_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
# outfile_4 = "ParVern9_outfile4_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
# outfile_5 = "ParVern9_outfile5_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
# outfile_bound = "ParVern9_outfile1_bound"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
#
#
# a = open(outfile_0, "w")
# b = open(outfile_1, "w")
# c = open(outfile_2, "w")
# d = open(outfile_3, "w")
# e = open(outfile_4, "w")
# f = open(outfile_5, "w")
# g = open(outfile_bound, "w")
#
# writedlm(outfile_0, T0)
# writedlm(outfile_1, T1)
# writedlm(outfile_2, T2)
# writedlm(outfile_3, T3)
# writedlm(outfile_4, T4)
# writedlm(outfile_5, T5)
# writedlm(outfile_bound, Tbound)
#
#
# close(a)
# close(b)
# close(c)
# close(d)
# close(e)
# close(f)
# close(g)
