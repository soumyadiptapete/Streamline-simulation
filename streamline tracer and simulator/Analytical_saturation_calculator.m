t=400;
x=actofi./t;
Sw=zeros(ny,nx);
Fw=zeros(ny,nx);


%calculating saturation at each cell

for i=1:ny
    for j=1:nx
        [Sw(i,j),Fw(i,j)]=buckley_leverett(x(i,j));
    end
end

contourf(Sw,50);