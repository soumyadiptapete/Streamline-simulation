
L=1;
dt=0.001;
dx=0.01;
nx=L/dx;
c=dt/dx;
Swc=0.2;
Sor=0.2;
Sw1=Swc*ones(60,nx+1);
Fw1=zeros(60,nx+1);
M=3;

%upstream
n=2;
while n<600
Sw1(n,1)=Sw1(n-1,1)+(-dt/dx)*Fw1(n-1,1)+(dt/dx);
S=(Sw1(n,1)-Swc)/(1-Swc-Sor);
Fw1(n,1)=M*S^2/(M*S^2+(1-S)^2);
for i=2:nx
    Sw1(n,i)=Sw1(n-1,i)+(-dt/dx)*Fw1(n-1,i)+(dt/dx)*Fw1(n-1,i-1);
    S=(Sw1(n,i)-Swc)/(1-Swc-Sor);
    Fw1(n,i)=M*S^2/(M*S^2+(1-S)^2);
end
% if(abs(C1(n,10)-C1(n-1,10))>t || C1(n,10)==0)
n=n+1;
% else 
% break;
% end
end
n1=n;

% Midpoint
n=2;
Sw2=Swc*ones(60,nx+1);
Fw2=zeros(60,nx+1);
while n<600
%     C2(n,1)=(1+(c*(1-c)/2))*C2(n-1,1)-(c*(1-c)/2)*C2(n-1,2)+c;
    Sw2(n,1)=Sw2(n-1,1)+(-c^2)*Fw2(n-1,1)-(c*(1-c)/2)*Fw2(n-1,2)+(c*(1+c)/2);
    S=(Sw2(n,1)-Swc)/(1-Swc-Sor);
    Fw2(n,1)=M*S^2/(M*S^2+(1-S)^2);
    for i=2:nx
        Sw2(n,i)=Sw2(n-1,i)+(-c^2)*Fw2(n-1,i)-(c*(1-c)/2)*Fw2(n-1,i+1)+(c*(1+c)/2)*Fw2(n-1,i-1);
        S=(Sw2(n,i)-Swc)/(1-Swc-Sor);
        Fw2(n,i)=M*S^2/(M*S^2+(1-S)^2);
    end
%     if(abs(C2(n,10)-C2(n-1,10))>t || C2(n,10)==0)
    n=n+1;
%     else
%         break;
%     end
end
n2=n;

% 2 point upstream
n=2;
 Sw3=Swc*ones(60,nx+1);
 Fw3=zeros(60,nx+1);
 while n<600
   Sw3(n,1)=Sw3(n-1,1)+(-c)*Fw3(n-1,1)+c-(c*(1-c)/2)*(Fw3(n-1,1)-2+1);
   Sw3(n,2)=Sw3(n-1,2)+(-c)*Fw3(n-1,2)+c*Fw3(n-1,1)-(c*(1-c)/2)*(Fw3(n-1,2)-2*Fw3(n-1,1)+1);
   S=(Sw3(n,1)-Swc)/(1-Swc-Sor);
   Fw3(n,1)=M*S^2/(M*S^2+(1-S)^2);
   S=(Sw3(n,2)-Swc)/(1-Swc-Sor);
   Fw3(n,2)=M*S^2/(M*S^2+(1-S)^2);

     for i=3:nx
         Sw3(n,i)=Sw3(n-1,i)+(-c)*Fw3(n-1,i)+c*Fw3(n-1,i-1)-(c*(1-c)/2)*(Fw3(n-1,i)-2*Fw3(n-1,i-1)+Fw3(n-1,i-2));
         S=(Sw3(n,i)-Swc)/(1-Swc-Sor);
         Fw3(n,i)=M*S^2/(M*S^2+(1-S)^2);
     end
%      if(abs(C3(n,10)-C3(n-1,10))>t || C3(n,10)==0)
     n=n+1;
%      else
%          break;
%      end
 end
 n3=n;
 
 %leonard scheme
 n=2;
 Sw4=Swc*ones(60,nx+1);
 Fw4=zeros(60,nx+1);
 while n<600
     Sw4(n,1)=Sw4(n-1,1)+(-c+(c*(1-c)^2)/2)*Fw4(n-1,1)+(1+c*(1-c)/2)*c+(c*(c-2)*(1-c)*Fw4(n-1,2)/6)+(c*(c-1)*(1+c)/6);
     Sw4(n,2)=Sw4(n-1,2)+(-c+(c*(1-c)^2)/2)*Fw4(n-1,2)+(1+c*(1-c)/2)*c*Fw4(n-1,1)+(c*(c-2)*(1-c)*Fw4(n-1,3)/6)+(c*(c-1)*(1+c)/6);
     S=(Sw4(n,1)-Swc)/(1-Swc-Sor);
     Fw4(n,1)=M*S^2/(M*S^2+(1-S)^2);
     S=(Sw4(n,2)-Swc)/(1-Swc-Sor);
     Fw4(n,2)=M*S^2/(M*S^2+(1-S)^2);
      
     for i=3:nx
         Sw4(n,i)=Sw4(n-1,i)+(-c+(c*(1-c)^2)/2)*Fw4(n-1,i)+(1+c*(1-c)/2)*c*Fw4(n-1,i-1)+(c*(c-2)*(1-c)*Fw4(n-1,i+1)/6)+(c*(c-1)*(1+c)*Fw4(n-1,i-2)/6);
         S=(Sw4(n,i)-Swc)/(1-Swc-Sor);
         Fw4(n,i)=M*S^2/(M*S^2+(1-S)^2);
     end
%      if(abs(C4(n,10)-C4(n-1,10))>t || C4(n,10)==0)
     n=n+1;
%      else
%          break;
%      end
 end
 n4=n;
 axis([0 1 0 1]);
 plot(1:100,Sw1(200,1:100),1:100,Sw2(200,1:100),1:100,Sw3(200,1:100),1:100,Sw4(200,1:100));

