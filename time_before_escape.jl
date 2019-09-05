

# @everywhere module_dir = "/Users/brandonbehring/Desktop/Leap_Frog_2019"
# @everywhere push!(LOAD_PATH, "/Users/brandonbehring/Desktop/Leap_Frog_2019")
# # @everywhere module_dir ="/Users/brandonbehring/Desktop/Leap_Frog_2019"
# # @everywhere push!(LOAD_PATH, $module_dir)
# # @everywhere thisDir = dirname(@__FILE__())
# # @everywhere any(path -> path == thisDir, LOAD_PATH) || push!(LOAD_PATH, thisDir)
# include("escape.jl")
# @everywhere using .escape

using ProgressMeter
using Printf
using MATLAB
include("escape_num.jl")
ty

@everywhere begin
    include("escape_num.jl")
    push!(LOAD_PATH, "/Users/brandonbehring/Desktop/Leap_Frog_2019")
    # push!(LOAD_PATH, "/home/brandon_behring/Desktop/Leap_Frog_2019")
    using ProgressMeter
    using Printf
    using MATLAB
    using .escape_num
    t_end = 1e3
    width = 1; height = 1
    n_iter_Q = 80;n_iter_P = 80;N = n_iter_P * n_iter_Q;
    ArrP = range(-width, stop = width, length = n_iter_P)
    ArrQ = range(-height, stop = height, length = n_iter_Q)
    mesh = [(Q, P) for Q in ArrQ, P in ArrP]
    mesh_list = reshape(mesh, 1, :)
    t_end = t_end * ones(n_iter_Q,n_iter_P);
    max_hits=5*ones(n_iter_Q,n_iter_P);
    H=.25
    # location = "/mnt/bdd38f66-9ece-451a-b915-952523c139d2/Escape/"
    location = "/Users/brandonbehring/Desktop/"
end


    @everywhere Energy = H* ones(n_iter_Q,n_iter_P);
    @everywhere h = replace(@sprintf("%.13f",H), "." => "_")
    @everywhere file_name = location * "Escape_" * h * ".fig"
    println(file_name)
    num_until_exit = @showprogress map(escape_num.escape_exit_num, mesh, t_end, Energy,max_hits)

    # @everywhere begin
    #     Q_0 = zeros(0);P_0 = zeros(0)
    #     Q_1 = zeros(0);P_1 = zeros(0)
    #     Q_2 = zeros(0);P_2 = zeros(0)
    #     Q_3 = zeros(0);P_3 = zeros(0)
    #     Q_4 = zeros(0);P_4 = zeros(0)
    #     Q_5 = zeros(0);P_5 = zeros(0)
    # end 
    # for i = 1:N
    #     for 
    #     Q, P = mesh_list[i]
    #     if num_until_exit[i] == 0
    #         push!(Q_0, Q)
    #         push!(P_0, P)
    #     elseif num_until_exit[i] == 1
    #         push!(Q_1, Q)
    #         push!(P_1, P)
    #     elseif num_until_exit[i] == 2
    #         push!(Q_2, Q)
    #         push!(P_2, P)
    #     elseif num_until_exit[i] == 3
    #         push!(Q_3, Q)
    #         push!(P_3, P)
    #     elseif num_until_exit[i] == 4
    #         push!(Q_4, Q)
    #         push!(P_4, P)
    #     elseif num_until_exit[i] == 5
    #         push!(Q_5, Q)
    #         push!(P_5, P)
    #     end
    # end
    m2=transpose(mesh)
    mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
    # mat"plot($P_0,$Q_0 ,'b.','MarkerSize',3)"
    # mat"plot($P_1,$Q_1 ,'r.','MarkerSize',3)"
    # mat"plot($P_2,$Q_2 ,'g.','MarkerSize',3)"
    # mat"plot($P_4,$Q_4 ,'c.','MarkerSize',3)"
    # mat"plot($P_5,$Q_5 ,'y.','MarkerSize',3)"
    mat"surf($mesh, $mesh' , $num_until_exit)"
    # mat"axis([ -$width,$width,-$height,$height ])"
    mat"axis([ -1,1,-1,1])"
    mat"savefig($file_name)"

