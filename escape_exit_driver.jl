using DifferentialEquations
using ProgressMeter
using Printf
include("escape_exit_function.jl")

max_hit=5
Energy=.2

Energy=.15
t_end=100.
n_iter_P=500
n_iter_Q=501
width=.03
height=.4
ArrP=range(-width,stop=width,length=n_iter_P)
ArrQ=range(-height,stop=height,length=n_iter_Q)
location="~/Desktop/"
h=replace(@sprintf("%.13f",Energy),"."=>"_")
file_name=location*"Escape_"*h*".fig"
println(file_name)


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
        K=escape_exit_function(ArrQ[k],ArrP[j], Energy, t_end, max_hit)
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



#
# mat"axis([ -$width,$width,-$height,$height ])"
mat"""
figure; hold on;
axis([ -$width,$width,-$height,$height ])
plot($P0,$Q0 ,'r.','MarkerSize',1)
plot($P1,$Q1 ,'b.','MarkerSize',1)
plot($P2,$Q2 ,'g.','MarkerSize',1)
plot($P3,$Q3 ,'m.','MarkerSize',1)
plot($P4,$Q4 ,'c.','MarkerSize',1)
plot($P5,$Q5 ,'y.','MarkerSize',1)
"""
mat"savefig($file_name),close"
