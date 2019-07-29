using ProgressMeter
using MATLAB
using Printf
include("PSS_function.jl")
# H=range(.12,stop=0.144281,length=211);
# H=[0.112012,.125789, .135]
H=range(.1,stop=0.15,length=10);
max_hit=100
width=.03
height=.4
h_start=replace(@sprintf("%.3f",H[1]),"."=>"_")
h_end=replace(@sprintf("%.3f",H[end]),"."=>"_")

location="~/Desktop/PSS_MOVIE/"
video_name=location*"PSS_MOVIE"*h_start*"to"*h_end*".avi"
mat"writerObj = VideoWriter($video_name);"

mat"""
writerObj.FrameRate=15;
writerObj.Quality = 100;
open(writerObj);
figure('position',[0, 0, 1200, 900])
hold on
"""
mat"axis([-$width, $width, -$height, $height])"

# mat"""
# set(gca,'XTick',[])
# set(gca,'YTick',[])
# set(gca,'xticklabel',[])
# set(gca,'yticklabel',[])
# set(gca,'XColor','none')
# set(gca,'YColor','none')
# set(gca,'XColor','none')
# set(gca,'YColor','none')
# """

n_iter_Q=11;#50
Q_start=-.25
Q_end=.25
n_iter_P=10;#51
P_start=-.015
P_end=.015
t_end=1e3

H_differences=zeros(0)

ArrP=range(P_start,stop=P_end,length=n_iter_P)
ArrQ=range(Q_start,stop=Q_end,length=n_iter_Q)

#defining colors for PSS
COLOR=["b", "g", "c", "m", "y", "k", "r"]
k=1;
@showprogress 1 "Computing..."for Energy in H
    # h=replace(@sprintf("%.8f",Energy),"."=>"_")
    # file_name=location*"PSS_"*h*".fig"
    mat"""
    figure('position',[0, 0, 1200, 900])
    hold on
    """
    mat"axis([-$width, $width, -$height, $height])"
    for j=1:n_iter_P
        for k=1:n_iter_Q
            Q,P,dH=PSS_function(ArrQ[k], ArrP[j], Energy, t_end, max_hit)
            if dH<5e-5 #make sure dH is ok
                current_color=COLOR[mod(k,length(COLOR))+1]
                # mat"plot($P_PSS,$Q_PSS,'.','MarkerSize',.1,'color',$current_color); hold on;"
                mat"plot($P_PSS,$Q_PSS,'k.','MarkerSize',.1); hold on;"
            end
        end
    end
    # mat"savefig($file_name)"
    if k>1
        mat"thisFrame = getframe;"
        mat"writeVideo(writerObj, thisFrame);"
    end
    mat"cla reset"
    k=k+1
end
mat"close(writerObj)"
