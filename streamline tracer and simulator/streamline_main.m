global nx ny  i1 j1 dx dy stx sty tofs actofs;
beta0_s=0.8;
for j1=1:4
n=nx-1;
m=ny;
alpha0_s=1;
i1=1;
stx(j1,i1)= sum(dx(m,1:n));
sty(j1,i1)=sum(dy(1:m-1,n))+(dy(m,n))*beta0_s;
tofs(j1)=0;
actofs(j1)=0; 
singlecell_streamline_tracer(alpha0_s,beta0_s,m,n,-qx,-qy);
a1(j1)=i1;
beta0_s=beta0_s-0.2;
end

alpha0_s=0.2;

for j1=5:8

    n=nx;
    m=ny-1;
    beta0_s=1;
    i1=1;
    stx(j1,i1)=sum(dx(m,1:n))-(1-alpha0_s)*(dx(m,n));
    sty(j1,i1)=sum(dy(1:m,n));
    tofs(j1)=0;
    actofs(j1)=0;
    singlecell_streamline_tracer(alpha0_s,beta0_s,m,n,-qx,-qy);
    alpha0_s=alpha0_s+0.2;
    a1(j1)=i1;
end
% plotting streamlines

axis([0 sum(dx(1,1:nx)) 0 sum(dy(1:ny,1))]);
plot(stx(1,1:a1(1)),sty(1,1:a1(1)),stx(2,1:a1(2)),sty(2,1:a1(2)),stx(3,1:a1(3)),sty(3,1:a1(3)),stx(4,1:a1(4)),sty(4,1:a1(4)),stx(5,1:a1(5)),sty(5,1:a1(5)),stx(6,1:a1(6)),sty(6,1:a1(6)),stx(7,1:a1(7)),sty(7,1:a1(7)),stx(8,1:a1(8)),sty(8,1:a1(8)));

%plotting time of flight for each streamline

plot(1:8,actofs);

%rate for each streamline


  