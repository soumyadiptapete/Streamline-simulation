i=1;
j=0;
for i=1:ny
    actofi1(i,:)=actofi(ny-j,:);
    j=j+1;
end
image(actofi1);