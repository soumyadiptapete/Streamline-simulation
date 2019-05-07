uw=1;
uo=3;
Swc=0.25;
Sor=0.25;
a=0.3/uw;
b=0.9/uo;
Sw=0.8;%uncomment the line if using for computation
% Sw=0.2:0.001:1;
S=(Sw-Swc)/(1-Swc-Sor);
%1st correlation
% Fw=(a.*(S.^2))./(a.*(S.^2)+b.*((1-S).^3));
% dFw= (((2.*a.*S)./(a.*(S.^2)+b.*((1-S).^3)))+((a.*(S.^2)).*(3*b.*((1-S).^2)-2.*a.*S))./(a.*(S.^2)+b.*((1-S).^3)).^2)./(1-Swc-Sor);
%2nd correlation.comment out 1st correlation if using 2nd one
% M=3;
% Fw=(M.*S.^2)./((M.*S.^2)+(1-S).^2);
% dFw=(((2.*M.*S)./((M.*S.^2)+(1-S).^2))-((M.*S.^2).*(2.*M.*S-2.*(1-S)))./((M.*S.^2)+(1-S).^2).^2).*(1/(1-Swc-Sor));

% figure();
% plot(Sw,Fw);
% figure();
% plot(Sw,dFw);
% single value of tau/t to saturation
% 
% x=tofs_bl./t;
% n=0;
% 
%     Sw=0.8;
% while(Sw>0.2);
%    S=(Sw-Swc)/(1-Swc-Sor);
%   dFw= (((2.*a.*S)./(a.*(S.^2)+b.*((1-S).^3)))+((a.*(S.^2)).*(3*b.*((1-S).^2)-2.*a.*S))./(a.*(S.^2)+b.*((1-S).^3)).^2)./(1-Swc-Sor);
% dFw=(((2.*M.*S)./((M.*S.^2)+(1-S).^2))-((M.*S.^2).*(2.*M.*S-2.*(1-S)))./((M.*S.^2)+(1-S).^2).^2).*(1/(1-Swc-Sor));
%   if(x>=2)
%       y=0.2;
%       z=0;
%   elseif x==0
%       y=0;
%       z=0;
%   else
%   if abs(x-dFw)<=0.001
%       y=Sw;
%       z=(a.*(S.^2))./(a.*(S.^2)+b.*((1-S).^3));
%       z=(M.*S.^2)./((M.*S.^2)+(1-S).^2);
%       break;
%   end
%   end
%   Sw=Sw-0.0001;
%   n=n+1;
% end
%     
% A 1-d array of values of tau/t to saturation values
% x=tofs_bl./t;
% n=0;
% for i=1:size(x)
% Sw=0.8;
% while(Sw>0.2);
%    S=(Sw-Swc)/(1-Swc-Sor);
%   dFw= (((2.*a.*S)./(a.*(S.^2)+b.*((1-S).^3)))+((a.*(S.^2)).*(3*b.*((1-S).^2)-2.*a.*S))./(a.*(S.^2)+b.*((1-S).^3)).^2)./(1-Swc-Sor);
% % dFw=(((2.*M.*S)./((M.*S.^2)+(1-S).^2))-((M.*S.^2).*(2.*M.*S-2.*(1-S)))./((M.*S.^2)+(1-S).^2).^2).*(1/(1-Swc-Sor));
%   if(x(i)>=2)
%       y(i)=0.2;
%       z(i)=0;
%   elseif x(i)==0
%       y(i)=0;
%       z(i)=0;
%   else
%   if abs(x(i)-dFw)<=0.001
%       y(i)=Sw;
%       z(i)=(a.*(S.^2))./(a.*(S.^2)+b.*((1-S).^3));
% %       z(i)=(M.*S.^2)./((M.*S.^2)+(1-S).^2);
%       break;
%   end
%   end
%   Sw=Sw-0.0001;
%   
% end
% end
% 
% % A 2-d array of values of tau/t to saturation values
x(41,41)=0;
y(41,41)=0;
z(41,41)=0;
x=tofs_bl./t;
[r,c]=size(x);
for i=1:r
    for j=1:c
     Sw=0.8;
while(Sw>0.2);
   S=(Sw-Swc)/(1-Swc-Sor);
%   dFw= (((2.*a.*S)./(a.*(S.^2)+b.*((1-S).^3)))+((a.*(S.^2)).*(3*b.*((1-S).^2)-2.*a.*S))./(a.*(S.^2)+b.*((1-S).^3)).^2)./(1-Swc-Sor);
  dFw= (((2*a*S)/(a*(S^2)+b*((1-S)^3)))+((a*(S^2))*(3*b*((1-S)^2)-2*a*S))/(a*(S^2)+b*((1-S)^3))^2)/(1-Swc-Sor);
% dFw=(((2.*M.*S)./((M.*S.^2)+(1-S).^2))-((M.*S.^2).*(2.*M.*S-2.*(1-S)))./((M.*S.^2)+(1-S).^2).^2).*(1/(1-Swc-Sor));
  if(x(i,j)>=2)
      y(i,j)=0.2;
%       z(i,j)=0;
  elseif x(i,j)==0
      y(i,j)=0;
%       z(i,j)=0;
  else
  if abs(x(i,j)-dFw)<=0.001
      y(i,j)=Sw;
%       z(i,j)=(a.*(S.^2))./(a.*(S.^2)+b.*((1-S).^3));
%       z(i,j)=(M.*S.^2)./((M.*S.^2)+(1-S).^2);
      break;
  end
  end
  Sw=Sw-0.0001;
  
end   
end
end

    