using Distributed
addprocs()

@everywhere function f(k,b,c)
     workernum = myid() - 1
     sleep(workernum)
     println("job $k")
    return b*[k k^2; k^2 k]+c*[1 1; 1 1]

end




n_iter_Q=10;#50
Q_start=-.25
Q_end=.25
n_iter_P=10;#51
P_start=-.015
P_end=.015


ArrP=range(P_start,stop=P_end,length=n_iter_P)
ArrQ=range(Q_start,stop=Q_end,length=n_iter_Q)
E=ones(6)*.5
A=[1,2,3,4,5,6]; B=[2,3,4,5,6,7]
pmap(f,A,B,E)
