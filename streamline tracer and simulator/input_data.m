global nx ny dx dy dz kx ky phi uo uw bo bw Sw Krw Kro dv;

nx=input('please enter number of cells in x direction\n');
ny=input('please enter number of cells in y direction\n');
% grid spacing
for i= 1:ny
    for j=1:nx
        dx(i,j)=25;
        dy(i,j)=25;
        dz(i,j)=10;
        kx(i,j)=100;
        ky(i,j)=100;
        phi(i,j)=0.3;
        uo(i,j)=3;
        uw(i,j)=1;
        bo(i,j)=1;
        bw(i,j)=1;
        Sw(i,j)=0.25;
%         Krw(i,j)=(5/3)*(Sw(i,j)-0.2);% corey relative perms with Swc=0.2 and Sor=0.2. At the end points relative perm is assumed to be 1.
%         Kro(i,j)=(4/3)-(5/3)*Sw(i,j);
        dv(i,j)=dx(i,j)*dy(i,j)*dz(i,j);
    end
end

% permeability alteration
for j=16:26
    kx(21,j)=1;
    ky(21,j)=1;
end

% rel perm calculation if the below correlation is given

Swc=0.25;
Sor=0.25;
%read saturation values
%calculate the rel perms
for i=1:ny
    for j=1:nx
        S=(Sw(i,j)-Swc)/(1-Swc-Sor);
        Krw(i,j)=0.3*S^2;
        Kro(i,j)=0.9*(1-S)^3;
    end
end
% well information

% rel perm calculation if the below correlation is given
%  Swc=0.2;
% Sor=0.2;
% %read saturation values
% %calculate the rel perms
% for i=1:ny
%     for j=1:nx
%         S=(Sw(i,j)-Swc)/(1-Swc-Sor);
%         Krw(i,j)=(1-S)^2;
%         Kro(i,j)=M*(uw(i,j)/uo(i,j))*S^2;
%     end
% end
n=input('please enter number of wells\n');
b1=1;
c1=1;
for a1=1:n
    Wellj=input('please enter X-index of the cell which has the well\n');
    Welli=input('please enter Y-index of the cell which has the well\n');
    i=input('please enter 0 if the well is rate constraint or enter 1 if the well is pressure constraint\n');
    if i==0
        Wx_r(b1)=Wellj;
        Wy_r(b1)=Welli;
        q_d= input(' Please enter the rate of the well in STB/D. Positive if injector and Negative if producer\n');
        q(Welli,Wellj)=5.615*q_d;%convert rate to scf/d
        b1=b1+1;
    elseif i==1
        Wx_p(c1)=Wellj;
        Wy_p(c1)=Welli;
        Pwf(Welli,Wellj)=input('please enter the bottomhole flowing pressure in psi\n');
        rw=input('please enter the wellbore radius in ft.\n');
        C(Welli,Wellj)=(0.0397*kx(Welli,Wellj)*dz(Welli,Wellj)/log(0.2*dx(Welli,Wellj)/rw))*((Kro(Welli,Wellj)/(uo(Welli,Wellj)*bo(Welli,Wellj)))+(Krw(Welli,Wellj)/(uw(Welli,Wellj)*bw(Welli,Wellj))));
        c1=c1+1;
    else
        disp('you did not enter the correct number\n');
    end
end

        
