

using DifferentialEquations
using MATLAB

# using PyCall
using ProgressMeter
# pygui(:qt)
# using PyPlot
# pygui(true)


function condition(u,t,integrator) # Event when event_f(u,t) == 0
 u[3]
end

function affect!(integrator)
end
#cb1 = SavingCallback((u,t,integrator)->(u[4],u[2]), saved_values)

cb2= ContinuousCallback(condition, affect!)
Energy=0.112012;
H=(2*Energy)^2
h=sqrt(H)
n_iter=20
Y_start=8e-5
Y_end=1.5e-2
t_end=10000
Y=range(Y_start,stop=Y_end,length=n_iter)
P=zeros(n_iter)
PofY_find(YY)=sqrt((h+YY.^2-h*YY.^2 )./(1+h-h*YY.^2))
H1(p,q,t)= ( (q[2]-q[1])^2+(p[2]-p[1])^2 )*( (q[2]+q[1])^2+(p[2]+p[1])^2 )
H2(p,q,t)= (p[2]^4+2*p[2]^2*(q[1]^2-1)+(1+q[1]^2)^2 )*(q[2]^4+2*q[2]^2*(p[1]^2+1)+(p[1]^2-1)^2 )
HAM(p,q,t)=H1(p,q,t)/H2(p,q,t)
PofY_find(YY)=sqrt((h+YY.^2-h*YY.^2 )./(1+h-h*YY.^2))
using MATLAB
    mat"figure; hold on;"
 @showprogress 1 "Computing..." for k=1:n_iter
# k=1
    q0,p0=[zeros(2) for i in 1:2]
    P[k]=PofY_find(Y[k])
    q0[1]=0
    p0[1]=P[k]
    q0[2]=0
    p0[2]=Y[k]
    prob= HamiltonianProblem{true}(HAM, p0, q0, (0., t_end))
    t,A=solve(prob, Vern9(),reltol=1e-13,abstol=1e-20,maxiters=1e15, callback=cb2,save_end=false,save_everystep=false)
    Qy=A[2,:];
    Px=A[4,:];

    mat"plot($Qy,$Px,'.',MarkerSize,1); hold on;"
#     plot!(A[2,:],-A[4,:],seriestype=:scatter,marker=:circle,msw = 0,ms=.1,seriescolor = :black,legend=false)
#     plot!(-A[2,:],A[4,:],seriestype=:scatter,marker=:circle,msw = 0,ms=.1,seriescolor = :black,legend=false)
#     plot(-A[2,:],-A[4,:],seriestype=:scatter,marker=:circle,msw = 0,ms=.1,seriescolor = :black,legend=false)

end
# QList=zeros(N)
# PList=zeros(N)
#
#
# for k in 1:N
#     QList[k],PList[k]=saved_values.saveval[k]
# end
#     plot(QList,PList,seriestype=:scatter,marker=:circle,msw = 0,ms=.1)
