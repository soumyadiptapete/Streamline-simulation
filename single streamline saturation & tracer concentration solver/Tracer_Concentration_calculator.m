
L=1;
dt=0.001;
dx=0.01;
nx=L/dx;
t=1e-12;
c=dt/dx;
C1=zeros(60,nx+1);
%upstream
n=2;
while n<600
C1(n,1)=(1-dt/dx)*C1(n-1,1)+(dt/dx);
for i=2:nx
    C1(n,i)=(1-dt/dx)*C1(n-1,i)+(dt/dx)*C1(n-1,i-1);
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
C2=zeros(60,nx+1);
while n<600
%     C2(n,1)=(1+(c*(1-c)/2))*C2(n-1,1)-(c*(1-c)/2)*C2(n-1,2)+c;
    C2(n,1)=(1-c^2)*C2(n-1,1)-(c*(1-c)/2)*C2(n-1,2)+(c*(1+c)/2);
    for i=2:nx
        C2(n,i)=(1-c^2)*C2(n-1,i)-(c*(1-c)/2)*C2(n-1,i+1)+(c*(1+c)/2)*C2(n-1,i-1);
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
 C3=zeros(60,nx+1);
 while n<600
    C3(n,1)=(1-c)*C3(n-1,1)+c-(c*(1-c)/2)*(C3(n-1,1)-2+1);
   C3(n,2)=(1-c)*C3(n-1,2)+c*C3(n-1,1)-(c*(1-c)/2)*(C3(n-1,2)-2*C3(n-1,1)+1);

     for i=3:nx
         C3(n,i)=(1-c)*C3(n-1,i)+c*C3(n-1,i-1)-(c*(1-c)/2)*(C3(n-1,i)-2*C3(n-1,i-1)+C3(n-1,i-2));
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
 C4=zeros(60,nx+1);
 while n<600
     C4(n,1)=(1-c+(c*(1-c)^2)/2)*C4(n-1,1)+(1+c*(1-c)/2)*c+(c*(c-2)*(1-c)*C4(n-1,2)/6)+(c*(c-1)*(1+c)/6);
     C4(n,2)=(1-c+(c*(1-c)^2)/2)*C4(n-1,2)+(1+c*(1-c)/2)*c*C4(n-1,1)+(c*(c-2)*(1-c)*C4(n-1,3)/6)+(c*(c-1)*(1+c)/6);
      
     for i=3:nx
         C4(n,i)=(1-c+(c*(1-c)^2)/2)*C4(n-1,i)+(1+c*(1-c)/2)*c*C4(n-1,i-1)+(c*(c-2)*(1-c)*C4(n-1,i+1)/6)+(c*(c-1)*(1+c)*C4(n-1,i-2)/6);
     end
%      if(abs(C4(n,10)-C4(n-1,10))>t || C4(n,10)==0)
     n=n+1;
%      else
%          break;
%      end
 end
 n4=n;
 plot(1:100,C1(400,1:100),1:100,C2(400,1:100),1:100,C3(400,1:100),1:100,C4(400,1:100));

