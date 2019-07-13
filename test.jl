include("YfindNP.jl")
Energy=.25
H=(2*Energy)^2

Q_start=0
Q_end=2
P_start=0
P_end=2

n_iter_P=5
n_iter_Q=5

delta_P=(P_end-P_start)/n_iter_P
delta_Q=(Q_end-Q_start)/n_iter_Q

ArrP= P_start: delta_P: P_end
ArrQ= Q_start: delta_Q: Q_end
Ylist=zeros(0)
for j=0:20
    for k=0:20
        Q=.1*j
        P=.1*k
        println(YfindNP(Q,P,H))
        println("\n")
    end
end
