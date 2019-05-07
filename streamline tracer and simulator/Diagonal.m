%j-1 diagonal

a=1;
for i=2: ny
    for j= 1: nx
        d(a)=-Ty(i,j);
        a=a+1;
    end
end

%  j+1 diagonal
a=1;
for i= 1: ny-1
    for j=1: nx
        e(a)=-Ty(i+1,j);
        a=a+1;
    end
end
        
% i-1 diagonal
a=1;
for i= 1: ny
    for j= 1: nx
        l(a)= -Tx(i,j);
        a=a+1;
    end
end
[j1, k1]=size(l);
b=l(2:k1);

%i+1 diagonal
a=1;
for i= 1:ny
    for j= 1:nx
        m(a)=-Tx(i,j+1);
        a=a+1;
    end
end
[j1, k1]=size(m);
c=m(1:k1-1);

% main diagonal

a=1;
for i= 1:ny
    for j= 1:nx
        g(a)= Tx(i,j+1)+Tx(i,j)+Ty(i+1,j)+Ty(i,j);
        
        for z= 1:size(Wx_p)
            if j==Wx_p(z) && i==Wy_p(z)
                g(a)= Tx(i,j+1)+Tx(i,j)+Ty(i+1,j)+Ty(i,j)+C(i,j);
            end
        end
            a=a+1;
      end
    end
    
                
                

        
    