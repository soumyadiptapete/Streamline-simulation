clear all;

L=1;
dtc=0.001;
dtg=0.001;
dx=0.01;
nx=L/dx;
t=400;
Sw_conv=zeros(200,nx+2);
Sw_diff=zeros(200,nx+2);
Sw=zeros(200,nx);

% initial conditions
Swc=0.2;
Sor=0.2;
phi=1;
M=3;
alphag=20;
Sw_conv(1,1:nx+2)=Swc;


n=1;
% operator splitting solution
for timestep=1:t
 
for i=2:nx
    
%convective step
    S=(Sw_conv(n,1)-Swc)/(1-Swc-Sor);
    Fw(n,1)= (M*S^2)/(M*S^2+(1-S)^2);
    Sw_conv(n+1,1)=Sw_conv(n,1)+ (dtc/phi)*(1-Fw(n,1))*(1/dx);
    S=(Sw_conv(n,i)-Swc)/(1-Swc-Sor);
    Fw(n,i)= (M*S^2)/(M*S^2+(1-S)^2);
Sw_conv(n+1,i)=Sw_conv(n,i)+ (dtc/phi)*(Fw(n,i-1)-Fw(n,i))*(1/dx);
end

Sw_conv(n+1,nx+1:nx+2)=Swc;

if n==1 
Sw_diff(n,1:nx+2)=Sw_conv(n+1,1:nx+2);% initial condition for gravity step
end

for i=2:nx
%diffusive step

S1=(Sw_conv(n+1,2)-Swc)/(1-Swc-Sor);
S2=(Sw_conv(n+1,1)-Swc)/(1-Swc-Sor);
Gw(n,1)=(M*S1^2*(1-S2)^2)/(M*S1^2+(1-S2)^2);
Sw_diff(n+1,1)=Sw_conv(n+1,1)+(dtg/phi)*(Gw(n,1))*(alphag/dx);
S1=(Sw_conv(n+1,i+1)-Swc)/(1-Swc-Sor);
S2=(Sw_conv(n+1,i)-Swc)/(1-Swc-Sor);
Gw(n,i)=(M*S1^2*(1-S2)^2)/(M*S1^2+(1-S2)^2);

Sw_diff(n+1,i)=Sw_conv(n+1,i)+(dtg/phi)*(Gw(n,i)-Gw(n,i-1))*(alphag/dx);
end


n=n+1;
Sw_conv(n,1:nx+2)=Sw_diff(n,1:nx+2); 
end

% convective solution
Sw_conv=zeros(400,nx+2);
Sw_conv(1,1:nx+2)=Swc;
n=1;
for timestep=1:t
  
for i=2:nx

    S=(Sw_conv(n,1)-Swc)/(1-Swc-Sor);
    Fw(n,1)= (M*S^2)/(M*S^2+(1-S)^2);
    Sw_conv(n+1,1)=Sw_conv(n,1)+ (dtc/phi)*(1-Fw(n,1))*(1/dx);
    S=(Sw_conv(n,i)-Swc)/(1-Swc-Sor);
    Fw(n,i)= (M*S^2)/(M*S^2+(1-S)^2);
Sw_conv(n+1,i)=Sw_conv(n,i)+ (dtc/phi)*(Fw(n,i-1)-Fw(n,i))*(1/dx);
Sw_conv(n+1,nx+1:nx+2)=Swc;

end
n=n+1;
end

% Without operator splitting
L=1;
dtc=0.001;
dtg=0.001;
dx=0.01;
nx=L/dx;
t=400;
Sw=zeros(200,nx+2);

% initial conditions
Swc=0.2;
Sor=0.2;
phi=1;


Sw(1,1:nx+2)=Swc;

n=1;

for timestep=1:t
 
for i=2:nx
    
%convective step
    S=(Sw(n,1)-Swc)/(1-Swc-Sor);
    Fw(n,1)= (M*S^2)/(M*S^2+(1-S)^2);
    S1=(Sw(n,2)-Swc)/(1-Swc-Sor);
    S2=(Sw(n,1)-Swc)/(1-Swc-Sor);
    if S1==0 && S2==1
        Gw(n,1)=0;
    else
    Gw(n,1)=(M*S1^2*(1-S2)^2)/(M*S1^2+(1-S2)^2);
    end
    Sw(n+1,1)=Sw(n,1)+ (dtc/phi)*(1-Fw(n,1))*(1/dx)+(dtg/phi)*(Gw(n,1))*(alphag/dx);
    S=(Sw(n,i)-Swc)/(1-Swc-Sor);
    Fw(n,i)= (M*S^2)/(M*S^2+(1-S)^2);
    S1=(Sw(n,i+1)-Swc)/(1-Swc-Sor);
    S2=(Sw(n,i)-Swc)/(1-Swc-Sor);
    if S1==0 && S2==1
        Gw(n,i)=0;
    else
    Gw(n,i)=(M*S1^2*(1-S2)^2)/(M*S1^2+(1-S2)^2);
    end
    Sw(n+1,i)=Sw(n,i)+ (dtc/phi)*(Fw(n,i-1)-Fw(n,i))*(1/dx)+(dtg/phi)*(Gw(n,i)-Gw(n,i-1))*(alphag/dx);
    Sw(n+1,nx+1:nx+2)=Swc;
end
n=n+1;
end


figure();
x=[0.01:0.01:1];
a=['Plot for M=',num2str(M),' and Gravity Number =',num2str(alphag),' at t=0.2 secs'];
plot(x,Sw_diff(200,1:nx),x,Sw(200,1:nx),x,Sw_conv(200,1:nx));
axis([0 1 0 1]);
title(a);
legend('Operator Splitting','No Operator Splitting','Convective Solution');
