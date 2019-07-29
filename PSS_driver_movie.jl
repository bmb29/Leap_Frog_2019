using ProgressMeter
using MATLAB
include("PSS_function.jl")


@showprogress 1 "Computing..." for Energy in H
    # H=range(.12,stop=0.144281,length=211);
    H=[0.112012,.125789, .135]
    H=.13
    max_hit=100
    width=.03
    height=.4
    location="~/Desktop/"
    h=replace(@sprintf("%.13f",Energy),"."=>"_")
    file_name=location*"PSS_"*h*".fig"
    video_name_fast=location*"PSS_"*h*"MOVIE.avi"

    mat"""
    writerObj = VideoWriter(video_name_fast);
    writerObj.FrameRate=15;
    writerObj.Quality = 100;
    open(writerObj);
    figure; hold on;
    set(gcf, 'Position',  [0, 0, 1500, 1000])
    figure('position',[0, 0, 1200, 900])
    hold on
    """
    mat"axis([-$width, $width, -$height, $height])"

    mat"""
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'XColor','none')
    set(gca,'YColor','none')
    set(gca,'XColor','none')
    set(gca,'YColor','none')
    """

    n_iter_Q=11;#50
    Q_start=-.25
    Q_end=.25
    n_iter_P=10;#51
    P_start=-.015
    P_end=.015
    t_end=1e4
    #create an empty list to store dH
    H_differences=zeros(0)
    #julia's version of linspace
    ArrP=range(P_start,stop=P_end,length=n_iter_P)
    ArrQ=range(Q_start,stop=Q_end,length=n_iter_Q)

    #defining colors for PSS
    COLOR=["b", "g", "c", "m", "y", "k", "r"]
    # COLOR=['b','g','c','m','y','k','r']

    #
    # mat"axis([ -3,3,-1.5, 1.5])"
    # @showprogress 1 "Computing..."
    for j=1:n_iter_P
        for k=1:n_iter_Q
            Q,P,dH=PSS_function(ArrQ[k], ArrP[j], Energy, t_end, max_hit)
            if dH<5e-5 #make sure dH is ok

                #push new dH to list
                # push!(H_differences,dH)
                #choose new color
                current_color=COLOR[mod(k,length(COLOR))+1]

                #take all symmetric combinations
                Q_PSS=vcat(Q)
                P_PSS=vcat(P)
                # Qy=A[2,:];
                # Px=A[4,:];
                mat"plot($P_PSS,$Q_PSS,'.','MarkerSize',.1,'color',$current_color); hold on;"
                # mat"plot($P_PSS,$Q_PSS,'k.','MarkerSize',.1); hold on;"
            end
        end
    end
    mat"savefig($file_name)"

    mat"""
    try
        thisFrame = getframe;
    end
    myMovie(frameIndex) = thisFrame;
    try
        writeVideo(writerObj, thisFrame);
    end
end
cla reset
"""
end
