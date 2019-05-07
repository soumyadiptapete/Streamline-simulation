% Velocity x direction

vx(1:ny,1)=0;
vx(1:ny,nx+1)=0;
for i=1:ny
    for j=2:nx
        Uo=(uo(i,j)+uo(i,j-1))/2;
        Bo=(bo(i,j)+bo(i,j-1))/2;
        Uw=(uw(i,j)+uw(i,j-1))/2;
        Bw=(bw(i,j)+bw(i,j-1))/2;
        Kxavg=(dx(i,j)+dx(i,j-1))*(((dx(i,j)/kx(i,j))+(dx(i,j-1)/kx(i,j-1)))^-1);
        %total particle velocity calculation in x direction
        vx(i,j)=-0.01266*(P(i,j)-P(i,j-1))*Kxavg*((Kro(i,j-1))/(Uo*Bo) + (Krw(i,j-1))/(Uw*Bw))/(dx(i,j)+dx(i,j-1));%ft/day
    end
end

%velocity in y direction
vy(1,1:nx)=0;
vy(ny+1,1:nx)=0;
for i= 2:ny
    for j=1:nx
        Uo=(uo(i,j)+uo(i-1,j))/2;
        Bo=(bo(i,j)+bo(i-1,j))/2;
        Uw=(uw(i,j)+uw(i-1,j))/2;
        Bw=(bw(i,j)+bw(i-1,j))/2;
        Kyavg=(dy(i,j)+dy(i-1,j))*(((dy(i,j)/ky(i,j))+(dy(i-1,j)/ky(i-1,j)))^-1);
        % total particle velocity calculation in y direction
        vy(i,j)=-0.01266*(P(i,j)-P(i-1,j))*Kyavg*((Kro(i-1,j))/(Uo*Bo) + (Krw(i-1,j))/(Uw*Bw))/(dy(i,j)+dy(i-1,j));%ft/day
    end
end

    % rates in x and y direction
qx(1:ny,1)=0;
qx(1:ny,nx+1)=0;
for i=1:ny
    for j=2:nx
        qx(i,j)=(vx(i,j)*dy(i,j)*dz(i,j)); % q= scf/day
    end
end

qy(1,1:nx)=0;
qy(ny+1,1:nx)=0;
for i= 2:ny
    for j=1:nx
        qy(i,j)=(vy(i,j)*dx(i,j)*dz(i,j)); % q= scf/day
    end
end


