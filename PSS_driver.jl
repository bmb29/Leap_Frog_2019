include("PSS_function.jl")
using ProgressMeter
using MATLAB
    mat"figure; hold on;"

Energy=.1251
t_end=1e3


n_iter_Q=10;
Q_start=0
Q_end=.5
n_iter_P=10;
P_start=1e-5
P_end=.9
#create an empty list to store dH
H_differences=zeros(0)
#julia's version of linspace
ArrP=range(P_start,stop=P_end,length=n_iter_P)
ArrQ=range(Q_start,stop=Q_end,length=n_iter_Q)

#defining colors for PSS
COLOR=["b", "g", "c", "m", "y", "k", "r"]
# COLOR=['b','g','c','m','y','k','r']
title("P-Q Plane for PSS X=0")
# axis([ -.003,.003,-.08, .08])

mat"axis([ -2.5,2.5,-1, 1])"
# @showprogress 1 "Computing..."
 @showprogress 1 "Computing..." for j=1:n_iter_P
    for k=1:n_iter_Q
        Q,P,dH=PSS_function(ArrQ[k], ArrP[j], Energy, t_end, max_hit)
        if dH<1e-5 #make sure dH is ok

            #push new dH to list
            push!(H_differences,dH)
            #choose new color
            current_color=COLOR[mod(k,length(COLOR))+1]

            #take all symmetric combinations
            Q_PSS=vcat(Q,Q,-Q,-Q)
            P_PSS=vcat(P,-P,P,-P)
            # Qy=A[2,:];
            # Px=A[4,:];

            mat"plot($P_PSS,$Q_PSS,'.','MarkerSize',.1,'color',$current_color); hold on;"
            # plot(P_PSS,Q_PSS,color=current_color,".",markersize=1, markeredgewidth=.1)
        end
    end
end
