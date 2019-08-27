using Distributed
addprocs(3)
# module_dir ="/Users/brandonbehring/Desktop/Leap_Frog_2019"
module_dir ="/Users/brandonbehring/Desktop/Leap_Frog_2019"
@everywhere thisDir = dirname(@__FILE__())
@everywhere any(path -> path == thisDir, LOAD_PATH) || push!(LOAD_PATH, thisDir)
# @everywhere push!(LOAD_PATH, $module_dir)
@everywhere using par_module
@everywhere N=1000
@everywhere E=ones(N)*.5
@everywhere A=1:N;
@everywhere B=A+ones(N)
pmap(par_module.f,A,B,E)
