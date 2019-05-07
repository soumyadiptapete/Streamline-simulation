%no flow boundaries on the sides of the grid
Tx(1:ny,1)=0;
Tx(1:ny,nx+1)=0;
Ty(1,1:nx)=0;
Ty(ny+1,1:nx)=0;
% transmissibility calculation for intercell boundaries
for i=1:ny
    for j=2:nx
        Uo=(uo(i,j)+uo(i,j-1))/2;
        Bo=(bo(i,j)+bo(i,j-1))/2;
        Uw=(uw(i,j)+uw(i,j-1))/2;
        Bw=(bw(i,j)+bw(i,j-1))/2;
        Kxavg=(dx(i,j)+dx(i,j-1))*(((dx(i,j)/kx(i,j))+(dx(i,j-1)/kx(i,j-1)))^-1);
        Tx(i,j)=(0.01266*dy(i,j)*dz(i,j)/(dx(i,j)+dx(i,j-1)))*Kxavg*((Kro(i,j-1))/(Uo*Bo) + (Krw(i,j-1))/(Uw*Bw)); % flow is assumed from left to right
        
    end
end

for i=2:ny
    for j=1:nx
        Uo=(uo(i,j)+uo(i-1,j))/2;
        Bo=(bo(i,j)+bo(i-1,j))/2;
        Uw=(uw(i,j)+uw(i-1,j))/2;
        Bw=(bw(i,j)+bw(i-1,j))/2;
        Kyavg=(dy(i,j)+dx(i-1,j))*(((dy(i,j)/ky(i,j))+(dy(i-1,j)/ky(i-1,j)))^-1);
        Ty(i,j)=(0.01266*dx(i,j)*dz(i,j)/(dy(i,j)+dy(i-1,j)))*Kyavg*((Kro(i-1,j))/(Uo*Bo) + (Krw(i-1,j))/(Uw*Bw)); % flow is assumed from bottom to top
    end
end

    