nx_grid=41;
ny_grid=41;
dtau=0.6;
dx=dtau;
dt=0.05;
time=300;%time at which saturation is to be calculated
timesteps=time/dt;

L=(floor(max(max(actofi))/dtau))*dtau+2*dtau;

nx=L/dx;
c=dt/dx;
M=3;
n=2;
Swc=0.2;
Sor=0.2;
Sw3=Swc*ones(60,nx+1);
Fw3=zeros(60,nx+1);
 while n<=timesteps
     disp('timesteps left');
     disp(timesteps-n);
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
 
 Sw=Sw3(timesteps,1:nx);
 
 for i=1:ny_grid
     for j=1:nx_grid
         if actofi(i,j)==0
             Sw_avg(i,j)=0;
         else
         Sw_avg(i,j)= (Sw(floor(actofi(i,j)/dtau))+Sw(floor(actofi(i,j)/dtau)+1))/2;
         end
     end
 end
 
         
         
         
 
