run input_data.m
run Transmissibility.m
run Diagonal.m


%transmissibility matrix calculation

T= diag(g,0)+diag(b,-1)+diag(c,1)+diag(d,-nx)+diag(e,nx);

% rhs column matrix
a=1;
for i= 1:ny
    for j= 1:nx
        z(a)=0;
        [j1, k1]=size(Wx_r);
        for x= 1:k1
            if j==Wx_r(x) && i==Wy_r(x)
                z(a)=q(i,j);
            end
        end
        [j2, k2]=size(Wx_p);
           for y=1:k2
            if j==Wx_p(y) && i==Wy_p(y)
                z(a)=C(i,j)*Pwf(i,j);
            end
           end
               
        a=a+1;
    end
end
                
% Pressure calculation
z=z';
P=T\z;
P=reshape(P,nx,ny);
P=P';

global e;
e=0.0000001;
run Velocities_rates.m
% run streamline_main.m
run Tof.m
