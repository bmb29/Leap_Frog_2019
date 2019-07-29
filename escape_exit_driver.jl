using DifferentialEquations
using ProgressMeter
using Printf
using MATLAB
include("escape_exit_function.jl")

max_hit=5
Energy=.25
t_end=1000.
n_iter_P=500
n_iter_Q=501
width=.99
height=1
ArrP=range(-width,stop=width,length=n_iter_P)
ArrQ=range(-height,stop=height,length=n_iter_Q)
location="~/Desktop/"
h=replace(@sprintf("%.13f",Energy),"."=>"_")
file_name=location*"Escape_"*h*".fig"
println(file_name)

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

@showprogress 1 "Computing..." for j=1:n_iter_P
    for k=1:n_iter_Q
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
        #
        # elseif K==-1
        #     push!(Q_bound,ArrQ[k])
        #     push!(P_bound,ArrP[j])
        end
    end
end



#
# mat"axis([ -$width,$width,-$height,$height ])"
mat"figure; hold on;"
mat"axis([ -$width,$width,-$height,$height ])"
# mat"plot($P_bound,$Q_bound ,'k.','MarkerSize',10)"
mat"plot($P_0,$Q_0 ,'b.','MarkerSize',2)"
mat"plot($P_1,$Q_1 ,'r.','MarkerSize',2)"
mat"plot($P_2,$Q_2 ,'g.','MarkerSize',2)"
mat"plot($P_4,$Q_4 ,'c.','MarkerSize',2)"
# mat"plot($P_5,$Q_5 ,'y.','MarkerSize',2)"
mat"savefig($file_name)"
