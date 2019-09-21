using DifferentialEquations
using ProgressMeter
using PyCall
# using MATLAB
using Printf
using BSON: @save, @load
# pygui(:qt)
using PyPlot
include("PSS_Definitions_Dimer_X.jl")
# pygui(true);

# H=range(.191,stop=0.25,length=20)

# @showprogress 1 "Computing..." for Energy in H
Energy=.185
location="/mnt/bdd38f66-9ece-451a-b915-952523c139d2/"
h=replace(@sprintf("%.15f",Energy),"."=>"_")
h_BSON=h*".bson"
h_title=@sprintf("%.8f",Energy)
file_name=location*"DIMER_PSS_X"*h*".fig"
# mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
# mat"xlabel('Q')"
# mat"ylabel('P')"
# mat"title($h_title)"
# mat"axis([-3, 3, -3, 3])"
t_end=1e3;
N_iter_Q=100;#50
Q_start=-3
Q_end=3
N_iter_P=50
P_start=-1.395
P_end=1.4
ArrP=range(P_start,stop=P_end,length=N_iter_P)
ArrQ=range(Q_start,stop=Q_end,length=N_iter_Q)
Brandons_Colors=["#393b79" ,"#5254a3","#6b6ecf","#9c9ede" ,"#637939","#8ca252" ,"#b5cf6b" ,"#cedb9c" ,"#8c6d31","#bd9e39" ,"#e7ba52","#e7cb94","#843c39","#ad494a" ,"#d6616b","#e7969c" ,"#7b4173" ,"#a55194","#ce6dbd" ,"#de9ed6"];
SAVE_DATA=Dict()
# figure()
@showprogress for k=1:N_iter_Q
    for j=1:N_iter_P
        Q_PSS,P_PSS=PSS_function(ArrQ[k], ArrP[j], Energy, t_end);
        SAVE_DATA[(ArrQ[k],ArrP[j])]=Q_PSS, P_PSS
        # current_color=Brandons_Colors[mod(k,length(Brandons_Colors))+1]
    end
end
# @showprogress for k=1:N_iter_P
#     Q_PSS,P_PSS=PSS_function(0,ArrP[k], Energy, t_end);
#     current_color=Brandons_Colors[mod(k,length(Brandons_Colors))+1]
#     SAVE_DATA[k+N_iter_Q]=Q_PSS, P_PSS
#     if Q_PSS!=0
#         # mat"plot($Q_PSS,$P_PSS,'.','MarkerSize',3,'color',$current_color); hold on;"
#         # mat"plot(-$Q_PSS,-$P_PSS,'.','MarkerSize',3,'color',$current_color); hold on;"
#     end
# end

@save h_BSON SAVE_DATA
mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
mat"xlabel('Q')"
mat"ylabel('P')"
mat"title($h_title)"
mat"axis([-3, 3, -3, 3])"
for k=1:N_iter_Q
    for j=1:N_iter_P
        Q_PSS,P_PSS= SAVE_DATA[(ArrQ[k],ArrP[j])]
        SAVE_DATA[(ArrQ[k],ArrP[j])]=Q_PSS, P_PSS
        current_color=Brandons_Colors[mod(k*j,length(Brandons_Colors))+1]
        mat"plot($Q_PSS,$P_PSS,'.','MarkerSize',3,'color',$current_color); hold on;"
        # mat"plot(-$Q_PSS,-$P_PSS,'.','MarkerSize',3,'color',$current_color); hold on;"
    end
end
mat"savefig($file_name)"
# end
