using DifferentialEquations
using Plots
using LinearAlgebra
H=range(0.01,stop=.49,length=100)
FLO=zeros(length(H))
two_mu=zeros(length(H))
u0=zeros(4)
function DeQ(du,u,p,t)
du[1]=-u[2]*(1+u[1]^2) / ( (u[1]^2+u[2]^2)*(1-u[2]^2));
du[2]=u[1]*(1-u[2]^2)/ ( (u[1]^2+u[2]^2)*(1+u[1]^2));
du[3]= 2*u[1]*u[2]/(u[1]^2+u[2]^2)^2 *u[3]-1/(u[2]^2+u[1]^2)^2* (1-u[2]^2)*(3*u[1]^4+u[2]^2*u[1]^2-u[2]^2+u[1]^2) /( (1+u[1]^2)^2 )*u[4];
du[4]=-(1/(u[1]^2+u[2]^2)^2* (1+u[1]^2)*(3*u[2]^4+u[1]^2*u[2]^2+u[1]^2-u[2]^2)/( (1-u[2]^2)^2) )*u[3]-(2*u[1]*u[2]/(u[1]^2+u[2]^2)^2)*u[4];
end
condition(u,t,integrator)= u[1]
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition,affect!)


for k=1:length(H)
    h=1/(2*H[k])
    Y0=1.0/sqrt(h[k]+1)

    tspan0=(0.,1e5)
    u0=[0.,Y0,0.,0.]
    prob = ODEProblem(DeQ,u0,tspan0)
    t,A=solve(prob, Vern9(),reltol=1e-13,abstol=1e-20,maxiters=1e15,callback=cb)
    FinalT=t.t[end]


    tspan=(0.,FinalT)
    u0=[0.,Y0,1.,0.]
    prob = ODEProblem(DeQ,u0,tspan)
    t,A1=solve(prob, Vern9(),reltol=1e-20,abstol=1e-12,maxiters=1e15)

    tspan=(0.,FinalT)
    u0=[0.,Y0,0.,1.]
    prob = ODEProblem(DeQ,u0,tspan)
    t,A2=solve(prob, Vern9(),reltol=1e-20,abstol=1e-12,maxiters=1e15)

    #Create Monodromy matrix
    a=A1[3,end]
    b=A1[4,end]
    c=A2[3,end]
    d=A2[4,end]
    M=[a c; b d]
    FLOQUET_MULTIPLIER= eigvals(M)
    FLOQUET_EXPONENT= abs(log(Complex(FLOQUET_MULTIPLIER[1])))
    FLO[k]=FLOQUET_EXPONENT
    two_mu[k]=2*FLOQUET_EXPONENT/FinalT
end


 plot(h,FLO)
 plot!(h,two_mu)
