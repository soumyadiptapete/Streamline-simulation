
L=1;
dt=0.001;
dx=0.01;
nx=L/dx;
c=dt/dx;
C1=zeros(60,nx+1);
%upstream
n=2;
while n<600
C1(n,1)=(1-dt/dx)*C1(n-1,1)+(dt/dx);
for i=2:nx
    C1(n,i)=(1-dt/dx)*C1(n-1,i)+(dt/dx)*C1(n-1,i-1);
end
n=n+1;
end
n1=n;

% Midpoint
n=2;
C2=zeros(60,nx+1);
while n<600
    r1=(C2(n-1,1)-1)/(C2(n-1,2)-C2(n-1,1));
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
        
    C2(n,1)=C2(n-1,1)+ c*(-1+(fl1+fl2)*(1-c)/2)*C2(n-1,1)-c*((1-c)/2)*fl1*C2(n-1,2)+(c-fl2*c*(1-c)/2);
    for i=2:nx
        r1=(C2(n-1,i)-C2(n-1,i-1))/(C2(n-1,i+1)-C2(n-1,i));
        if i==2
        r2=(C2(n-1,1)-1)/(C2(n-1,2)-C2(n-1,1));
        else
          r2=(C2(n-1,i-1)-C2(n-1,i-2))/(C2(n-1,i)-C2(n-1,i-1));
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
    
    C2(n,i)=C2(n-1,i)+ c*(-1+(fl1+fl2)*(1-c)/2)*C2(n-1,i)-c*((1-c)/2)*fl1*C2(n-1,i+1)+(c-fl2*c*(1-c)/2)*C2(n-1,i-1);
    end

    n=n+1;

end
n2=n;

% 2 point upstream
n=2;
 C3=zeros(60,nx+1);
 while n<600
     r1=(C3(n-1,1)-1)/(C3(n-1,2)-C3(n-1,1));
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
        
    C3(n,1)=C3(n-1,1)+ c*(-1+(fl1+fl2)*(1-c)/2)*C3(n-1,1)-c*((1-c)/2)*fl1*C3(n-1,2)+(c-fl2*c*(1-c)/2);
    for i=2:nx
        r1=(C3(n-1,i)-C3(n-1,i-1))/(C3(n-1,i+1)-C3(n-1,i));
        if i==2
        r2=(C3(n-1,1)-1)/(C3(n-1,2)-C3(n-1,1));
        else
          r2=(C3(n-1,i-1)-C3(n-1,i-2))/(C3(n-1,i)-C3(n-1,i-1));
        end
    if r1>0
        if r1<=2
            fl1=r1;
        else
            fl1=2;
        end
    else
        fl1=0;
    end
     if r2>0
        if r2<=2
            fl2=r2;
        else
            fl2=2;
        end
    else
        fl2=0;
    end
    
    C3(n,i)=C3(n-1,i)+ c*(-1+(fl1+fl2)*(1-c)/2)*C3(n-1,i)-c*((1-c)/2)*fl1*C3(n-1,i+1)+(c-fl2*c*(1-c)/2)*C3(n-1,i-1);
    end
n=n+1;
 end
 n3=n;
 
%  %leonard scheme
n=2;
 C4=zeros(60,nx+1);
 while n<600
     r1=(C4(n-1,1)-1)/(C4(n-1,2)-C4(n-1,1));
     fl1=((2-c)/3)+((1+c)/3)*r1;
    r2=0;
    fl2=0;
    if r1>0
        if fl1/r1>0 && fl1<=2
            fl1=fl1;
        else
            fl1=0;
        end
    else
        fl1=0;
    end
        
    C4(n,1)=C4(n-1,1)+ c*(-1+(fl1+fl2)*(1-c)/2)*C4(n-1,1)-c*((1-c)/2)*fl1*C4(n-1,2)+(c-fl2*c*(1-c)/2);
    for i=2:nx
        r1=(C4(n-1,i)-C4(n-1,i-1))/(C4(n-1,i+1)-C4(n-1,i));
        if i==2
        r2=(C4(n-1,1)-1)/(C4(n-1,2)-C4(n-1,1));
        else
          r2=(C4(n-1,i-1)-C4(n-1,i-2))/(C4(n-1,i)-C4(n-1,i-1));
        end
        
         fl1=((2-c)/3)+((1+c)/3)*r1;
         fl2=((2-c)/3)+((1+c)/3)*r2;
    if r1>0
        if fl1/r1>0 && fl1<=2
            fl1=fl1;
        else
            fl1=0;
        end
    else
        fl1=0;
    end
     if r2>0
        if fl2/r2>0 && fl2<=2
            fl2=fl2;
        else
            fl2=0;
        end
    else
        fl2=0;
    end
    
    C4(n,i)=C4(n-1,i)+ c*(-1+(fl1+fl2)*(1-c)/2)*C4(n-1,i)-c*((1-c)/2)*fl1*C4(n-1,i+1)+(c-fl2*c*(1-c)/2)*C4(n-1,i-1);
    end
n=n+1;
 end
 n4=n;
 plot(1:100,C1(400,1:100),1:100,C2(400,1:100),1:100,C3(400,1:100),1:100,C4(400,1:100));
 
