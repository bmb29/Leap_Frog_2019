# using IntervalArithmetic
using DifferentialEquations
using LinearAlgebra
using TaylorIntegration

function Flo(du,u,p,t)
    h=BigFloat(".125")
    A11=(-1) *(1+4. *h^2+4. *h*cos(2. *t)).^(-1/2).*sin(2. *t)
    A12=(1+4. *h^2+(1+4. *h^2+4. *h*cos(2. *t)).^(1/2)+2. *h*cos(2. *t).*(2+(1+4. *h^2+ 4. *h*cos(2. *t)).^(1/2))).^(-1) *((-1) *cos(2. *t).*(1+(-8).*h^2+(1+4. *h^2+4. *h*cos(2. *t)).^(1/2))+h*(1+(-4).*h+16. *h^2+2. *(1+ 4. *h^2+4. *h*cos(2. *t)).^(1/2)+(-8).*h*(1+4. *h^2+4. *h*cos(2. * t)).^(1/2)+(-1) *cos(4. *t)))
    A21=(1+4. *h^2+(1+4. *h^2+4. *h*cos(2. *t) ).^(1/2)+2. *h*cos(2. *t).*(2+(1+4. *h^2+4. *h*cos(2. *t)).^(1/2))).^(-1) *(h+4. *h^2. *(1+4. *h)+2. *h*(1+4. *h).*(1+4. *h^2+4. *h*cos(2. *t)).^(1/2)+(-1) *cos(2. *t).*(1+(-8).*h^2+(1+4. *h^2+4. *h*cos( 2. *t)).^(1/2))+(-1) *h*cos(4. *t))
    A22=(1+4. *h^2+4. *h*cos(2. *t)).^( .-1/2).*sin(2. *t)

    du[1]=A11*u[1]+A12*u[2]
    du[2]=A21*u[1]+A22*u[2]
end
#

T_end=BigFloat(pi)
tspan=(BigFloat("0"),T_end)
tol=BigFloat("1e-100")

u0=[BigFloat("1"),BigFloat("0")]
prob = ODEProblem(Flo,u0,tspan)

t,u=taylorinteg(Flo,u0, tspan[1],tspan[2],50,1e-80, maxsteps=1_000_000_000)
# t,u1=solve(prob, Feagin14(),reltol=relT,abstol=absT,maxiters=1e20)
diff=sqrt((u[end,1]-u0[1])^2+(u[end,2]-u0[2])^2)
# diff=sqrt((u1[1,end]-u1[1,1])^2+(u1[2,end]-u1[2,1])^2)
println(diff)
