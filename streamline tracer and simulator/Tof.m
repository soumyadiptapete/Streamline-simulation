% start the outer loop for selecting the initial alpha0 and beta0 values
% at centre of each cell


global p q tofp actofp tofi actofi;


for p= 1:ny

    for q= 1:nx
        tofp(p,q)=0;
        actofp(p,q)=0;
        alpha0_s=0.5;
        beta0_s=0.5;
        % select cell indices
        n=q;
        m=p;
        singlecell_tofp_calculator(alpha0_s,beta0_s,m,n,qx,qy);
    end
  
end

        
 for p= 1:ny

    for q= 1:nx
        tofi(p,q)=0;
        actofi(p,q)=0;
        alpha0_s=0.5;
        beta0_s=0.5;
        % select cell indices
        n=q;
        m=p;
        singlecell_tofi_calculator(alpha0_s,beta0_s,m,n,-qx,-qy);
    end
  
 end          
 
TOF=actofp+actofi;
image(TOF,'CDataMapping','scaled');

        