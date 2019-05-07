function [ Sw_d Fw_d ] = buckley_leverett( x_d )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Sor=0.2;
Swc=0.2;
M=0.3;
Sw=0.8;
uw=1;
uo=1;
a=0.3/uw;
b=0.9/uo;

while(Sw>0.2)
 
 S=(Sw-Swc)/(1-Swc-Sor);
%  Fw=(M*S^2)/(M*S^2+(1-S)^2);
% dFw=(2*M*S)/(M*S^2+(1-S)^2)-(M*S^2*(2*M*S-2*(1-S))/(M*S^2+(1-S)^2)^2);
Fw=(a*S^2/(a*S^2+b*(1-S)^3));
dFw=(2*a*S/(a*S^2+b*(1-S)^3))-(a*S^2*(2*a*S-3*b*(1-S)^2))/(a*S^2+b*(1-S)^3)^2;
    if x_d>=2
        Sw_d=0.2;
        Fw_d=0;
    elseif x_d==0
        Sw_d=0;
        Fw_d=0;
    else
        if abs(x_d-dFw)<=0.001
            Sw_d=Sw;
            Fw_d=Fw;
            break;
        end
    end
    Sw=Sw-0.0001;
end
        

end

