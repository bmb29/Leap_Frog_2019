include("PSS_function.jl")
using ProgressMeter
using MATLAB


 # @showprogress 1 "Computing..."
H=[.09,.095,.1,.115,.11]
max_hit=50
for Energy in H
mat"figure; hold on;"
location="~/MATLAB-Drive/PSS_LEAP_FROG_CENTER/"

h=replace(string(Energy),"."=>"_")
file_name=location*h*".fig"
println(file_name)

println(max_hit)
# n_iter_Q=100;
# Q_start=-1
# Q_end=1
# n_iter_P=401;
# P_start=-2
# P_end=2
# t_end=1000

n_iter_Q=30;
Q_start=-.25
Q_end=.25
n_iter_P=31;
P_start=-.015
P_end=.015
t_end=1e3
#create an empty list to store dH
H_differences=zeros(0)
#julia's version of linspace
ArrP=range(P_start,stop=P_end,length=n_iter_P)
ArrQ=range(Q_start,stop=Q_end,length=n_iter_Q)

#defining colors for PSS
COLOR=["b", "g", "c", "m", "y", "k", "r"]
# COLOR=['b','g','c','m','y','k','r']
# axis([ -.003,.003,-.08, .08])
mat"axis([ -.03,.03,-.4,.4 ])"
#
# mat"axis([ -3,3,-1.5, 1.5])"
# @showprogress 1 "Computing..."
for j=1:n_iter_P
    for k=1:n_iter_Q
        Q,P,dH=PSS_function(ArrQ[k], ArrP[j], Energy, t_end, max_hit)
        if dH<5e-12 #make sure dH is ok

            #push new dH to list
            push!(H_differences,dH)
            #choose new color
            current_color=COLOR[mod(k,length(COLOR))+1]

            #take all symmetric combinations
            Q_PSS=vcat(Q)
            P_PSS=vcat(P)
            # Qy=A[2,:];
            # Px=A[4,:];

            # mat"plot($P_PSS,$Q_PSS,'.','MarkerSize',.1,'color',$current_color); hold on;"
            mat"plot($P_PSS,$Q_PSS,'k.','MarkerSize',.1); hold on;"

            # plot(P_PSS,Q_PSS,color=current_color,".",markersize=1, markeredgewidth=.1)
        end
    end
end
mat"savefig($file_name),close"
end
