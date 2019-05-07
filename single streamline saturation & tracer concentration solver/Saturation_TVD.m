
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
    r1=(Fw2(n-1,1)-1)/(Fw2(n-1,2)-Fw2(n-1,1));
    r2=0;
    fl2=0;
    if r1>0
        if 2*r1<=1
            fl1=2*r1;
        else
            fl1=1;
        end
    else
        fl1=0;
    end
        
    Sw2(n,1)=Sw2(n-1,1)+ c*(-1+(fl1+fl2)*(1-c)/2)*Fw2(n-1,1)-c*((1-c)/2)*fl1*Fw2(n-1,2)+(c-fl2*c*(1-c)/2);
    S=(Sw2(n,1)-Swc)/(1-Swc-Sor);
    Fw2(n,1)=M*S^2/(M*S^2+(1-S)^2);
    for i=2:nx
        r1=(Fw2(n-1,i)-Fw2(n-1,i-1))/(Fw2(n-1,i+1)-Fw2(n-1,i));
        if i==2
        r2=(Fw2(n-1,1)-1)/(Fw2(n-1,2)-Fw2(n-1,1));
        else
          r2=(Fw2(n-1,i-1)-Fw2(n-1,i-2))/(Fw2(n-1,i)-Fw2(n-1,i-1));
        end
    if r1>0
        if 2*r1<=1
            fl1=2*r1;
        else
            fl1=1;
        end
    else
        fl1=0;
    end
    if r2>0
        if 2*r2<=1
            fl2=2*r2;
        else
            fl2=1;
        end
    else
        fl2=0;
    end
        Sw2(n,i)=Sw2(n-1,i)+ c*(-1+(fl1+fl2)*(1-c)/2)*Fw2(n-1,i)-c*((1-c)/2)*fl1*Fw2(n-1,i+1)+(c-fl2*c*(1-c)/2)*Fw2(n-1,i-1);
        S=(Sw2(n,i)-Swc)/(1-Swc-Sor);
        Fw2(n,i)=M*S^2/(M*S^2+(1-S)^2);
    end
n=n+1;
end
n2=n;

% 2 point upstream
n=2;
Sw3=Swc*ones(60,nx+1);
Fw3=zeros(60,nx+1);
 while n<600
    r1=(Fw3(n-1,1)-1)/(Fw3(n-1,2)-Fw3(n-1,1));
    r2=0;
    fl2=0;
    if r1>0
        if r1<=2
            fl1=r1;
        else
            fl1=2;
        end
    else
        fl1=0;
    end    
    Sw3(n,1)=Sw3(n-1,1)+ c*(-1+(fl1+fl2)*(1-c)/2)*Fw3(n-1,1)-c*((1-c)/2)*fl1*Fw3(n-1,2)+(c-fl2*c*(1-c)/2);
    S=(Sw3(n,1)-Swc)/(1-Swc-Sor);
    Fw3(n,1)=M*S^2/(M*S^2+(1-S)^2);
    for i=2:nx
        r1=(Fw3(n-1,i)-Fw3(n-1,i-1))/(Fw3(n-1,i+1)-Fw3(n-1,i));
        if i==2
        r2=(Fw3(n-1,1)-1)/(Fw3(n-1,2)-Fw3(n-1,1));
        else
          r2=(Fw3(n-1,i-1)-Fw3(n-1,i-2))/(Fw3(n-1,i)-Fw3(n-1,i-1));
        end
    if r1>0
        if r1<=2
            fl1=r1;
        else
            fl1=1;
        end
    else
        fl1=0;
    end
    if r2>0
        if r2<=2
            fl2=r2;
        else
            fl2=1;
        end
    else
        fl2=0;
    end
        Sw3(n,i)=Sw3(n-1,i)+ c*(-1+(fl1+fl2)*(1-c)/2)*Fw3(n-1,i)-c*((1-c)/2)*fl1*Fw3(n-1,i+1)+(c-fl2*c*(1-c)/2)*Fw3(n-1,i-1);
        S=(Sw3(n,i)-Swc)/(1-Swc-Sor);
        Fw3(n,i)=M*S^2/(M*S^2+(1-S)^2);
    end
     n=n+1;
 end
 n3=n;
 
%  %leonard scheme
n=2;
Sw4=Swc*ones(60,nx+1);
Fw4=zeros(60,nx+1);
 while n<600
    r1=(Fw4(n-1,1)-1)/(Fw4(n-1,2)-Fw4(n-1,1));
    r2=0;
    fl2=(2-c)/3;
    fl1=((2-c)/3)+((1+c)/3)*r1;
    if r1>0
        if fl1/r1>=0 && fl1<=2 && fl1/r1<=2
            fl1=fl1;
        else
            fl1=0;
        end
    else
        fl1=0;
    end
    Sw4(n,1)=Sw4(n-1,1)+ c*(-1+(fl1+fl2)*(1-c)/2)*Fw4(n-1,1)-c*((1-c)/2)*fl1*Fw4(n-1,2)+(c-fl2*c*(1-c)/2);
    S=(Sw4(n,1)-Swc)/(1-Swc-Sor);
    Fw4(n,1)=M*S^2/(M*S^2+(1-S)^2);
    for i=2:nx
        r1=(Fw4(n-1,i)-Fw4(n-1,i-1))/(Fw4(n-1,i+1)-Fw4(n-1,i));
        if i==2
        r2=(Fw4(n-1,1)-1)/(Fw4(n-1,2)-Fw4(n-1,1));
        else
          r2=(Fw4(n-1,i-1)-Fw4(n-1,i-2))/(Fw4(n-1,i)-Fw4(n-1,i-1));
        end
        fl1=((2-c)/3)+((1+c)/3)*r1;
        fl2=((2-c)/3)+((1+c)/3)*r2;
    if r1>0
        if fl1/r1>=0 && fl1<=2 && fl1/r1<=2
            fl1=fl1;
        else
            fl1=0;
        end
    else
        fl1=0;
    end
    if r2>0
        if fl2/r2>=0 && fl2<=2 && fl2/r2<=2
            fl2=fl2;
        else
            fl2=0;
        end
    else
        fl2=0;
    end
        Sw4(n,i)=Sw4(n-1,i)+ c*(-1+(fl1+fl2)*(1-c)/2)*Fw4(n-1,i)-c*((1-c)/2)*fl1*Fw4(n-1,i+1)+(c-fl2*c*(1-c)/2)*Fw4(n-1,i-1);
        S=(Sw4(n,i)-Swc)/(1-Swc-Sor);
        Fw4(n,i)=M*S^2/(M*S^2+(1-S)^2);
    end
     n=n+1;
 end
 n4=n;
 axis([0 1 0 1]);
 plot(1:100,Sw1(300,1:100),1:100,Sw2(300,1:100),1:100,Sw3(300,1:100),1:100,Sw4(300,1:100));


