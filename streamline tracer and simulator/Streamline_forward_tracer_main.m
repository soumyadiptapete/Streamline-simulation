global i1 j1 dx dy stx sty tofs actofs actofp_center  tof_curr1;

actofp_center=0;
beta0_s=[0.2 0.4 0.6 0.8];
for j1=1:4
    
n=2;
m=1;
alpha0_s=0;
i1=1;
stx(j1,i1)= sum(dx(m,1:n-1));
sty(j1,i1)=(dy(m,n))*beta0_s(j1);
tofs(j1)=0;
actofs(j1)=0; 
tof_curr1(j1,m,n)=0;
actofp_center(j1,m,n)=0;
singlecell_forward_tracer(alpha0_s,beta0_s(j1),m,n,qx,qy);
a1(j1)=i1;

end

alpha0_s=[0.8 0.6 0.4 0.2];

for j1=5:8
    
    n=1;
    m=2;
    beta0_s=0;
    i1=1;
    stx(j1,i1)=alpha0_s(j1-4)*dx(m-1,n);
    sty(j1,i1)=sum(dy(1:m-1,n));
    tofs(j1)=0;
    actofs(j1)=0;
    tof_curr1(j1,m,n)=0;
    actofp_center(j1,m,n)=0;
    singlecell_forward_tracer(alpha0_s(j1-4),beta0_s,m,n,qx,qy);
   
    a1(j1)=i1;
end