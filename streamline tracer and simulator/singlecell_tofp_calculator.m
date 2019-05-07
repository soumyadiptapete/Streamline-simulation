function singlecell_tofp_calculator(alpha0_d,beta0_d,m_d,n_d,qx_d,qy_d )
global p q tofp actofp phi dv nx ny e;
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
        qx_d1=qx_d;
        qy_d1=qy_d;
        alpha0=alpha0_d;
        beta0=beta0_d;
        testa=alpha0_d;
        testb=beta0_d;
 if(m_d<=ny && m_d>=1 && n_d<=nx && n_d>=1)      
        q11=qx_d1(m_d,n_d);
        q12=qx_d1(m_d,n_d+1);  
        q21=qy_d1(m_d,n_d); 
        q22=qy_d1(m_d+1,n_d);
       
        
        dq1=q12-q11;
        dq2=q22-q21;

        if dq1~=0
            argalpha1=(q11)/(q11+alpha0*dq1);
            testalpha1=argalpha1;
            argalpha2=(q11+dq1)/(q11+alpha0*dq1);
            testalpha2=argalpha2;
        end
        if dq2~=0
            argbeta1=(q21)/(q21+beta0*dq2);
            argbeta2=(q21+dq2)/(q21+beta0*dq2);
        end
           
end
        %test boundary conditions
while(m_d<=ny && m_d>=1 && n_d<=nx && n_d>=1)
        % straight line conditions
        if abs(dq1)<=e && abs(dq2)>e
           t11= -alpha0/q11;
           t12=(1-alpha0)/q11;
           
           if argbeta1<=0
               t21=100000;
           else
               t21=(1/dq2)*log(argbeta1);
           end
           
           if argbeta2<=0
               t22=100000;
           else
               t22=(1/dq2)*log(argbeta2);
           end
           t=[t11 t12 t21 t22];
           tof_curr=tofc(t);
           alpha=alpha0+q11*tof_curr;
           beta=beta0*exp(dq2*tof_curr)+q21*((exp(dq2*tof_curr)-1)/(dq2));
           tofp(p,q)= tofp(p,q)+tof_curr;
           actofp(p,q)=actofp(p,q)+(tof_curr*phi(p,q)*dv(p,q));
           [m_d2,n_d2,alpha0,beta0]=Exit_face(tof_curr,t11,t12,t21,t22,alpha,beta,m_d,n_d);
           singlecell_tofp_calculator(alpha0,beta0,m_d2,n_d2,qx_d1,qy_d1);
           break;
           
        elseif(abs(dq2)<=e && abs(dq1)>e)
             t21= -beta0/q21;
           t22=(1-beta0)/q21;
           
           if argalpha1<=0
               t11=100000;
           else
               t11=(1/dq1)*log(argalpha1);
           end
           
           if argalpha2<=0
               t12=100000;
           else
               t12=(1/dq1)*log(argalpha2);
           end
           t=[t11 t12 t21 t22];
           tof_curr=tofc(t);
           beta=beta0+q21*tof_curr;
           alpha=alpha0*exp(dq1*tof_curr)+q11*((exp(dq1*tof_curr)-1)/(dq1));
           tofp(p,q)= tofp(p,q)+tof_curr;
           actofp(p,q)=actofp(p,q)+(tof_curr*phi(p,q)*dv(p,q));
           [m_d2,n_d2,alpha0,beta0]=Exit_face(tof_curr,t11,t12,t21,t22,alpha,beta,m_d,n_d);
           singlecell_tofp_calculator(alpha0,beta0,m_d2,n_d2,qx_d1,qy_d1);
           break;
           
        elseif(abs(dq1)<=e && abs(dq2)<=e)
            
           t11= -alpha0/q11;
           t12=(1-alpha0)/q11;
           t21= -beta0/q21;
           t22=(1-beta0)/q21;
           t=[t11 t12 t21 t22];
           tof_curr=tofc(t);
           beta=beta0+q21*tof_curr;
           alpha=alpha0+q11*tof_curr;
           tofp(p,q)= tofp(p,q)+tof_curr;
           actofp(p,q)=actofp(p,q)+(tof_curr*phi(p,q)*dv(p,q));
           [m_d2,n_d2,alpha0,beta0]=Exit_face(tof_curr,t11,t12,t21,t22,alpha,beta,m_d,n_d);
           singlecell_tofp_calculator(alpha0,beta0,m_d2,n_d2,qx_d1,qy_d1);
           break;
           
        % producer cell
        elseif(q11>=0 && q12<=0 && q21>=0 && q22<=0)
            tof_curr=0;
            tofp(p,q)=tofp(p,q)+tof_curr;
            actofp(p,q)=actofp(p,q)+ tof_curr;
            break;
            
        %injector cell
        elseif(q11<=0 && q12>=0 && q21<=0 && q22>=0)
            tof_curr=0;
            tofp(p,q)=tofp(p,q)+tof_curr;
            actofp(p,q)=actofp(p,q)+ tof_curr;
            alpha0=0;
            beta0=0;
            m_d2=m_d+1;
            n_d2=n_d+1;
            singlecell_tofp_calculator(alpha0,beta0,m_d2,n_d2,qx_d1,qy_d1);
            break;
            
         % stagnation point   
        elseif(q11==0 && q12==0 && q21==0 && q22==0)
            tof_curr=0;
            tofp(p,q)=tofp(p,q)+tof_curr;
            actofp(p,q)=actofp(p,q)+ tof_curr;
            break;
            
            
        else
            if(argalpha1<=0)
                t11=100000;
            else
                t11=(1/dq1)*log(argalpha1);
            end
            if(argalpha2<=0)
                t12=100000;
            else
                t12=(1/dq1)*log(argalpha2);
            end
           if argbeta1<=0
               t21=100000;
           else
               t21=(1/dq2)*log(argbeta1);
           end
           
           if argbeta2<=0
               t22=100000;
           else
               t22=(1/dq2)*log(argbeta2);
           end
           t=[t11 t12 t21 t22];
           tof_curr=tofc(t); 
           beta=beta0*exp(dq2*tof_curr)+q21*((exp(dq2*tof_curr)-1)/(dq2));
           alpha=alpha0*exp(dq1*tof_curr)+q11*((exp(dq1*tof_curr)-1)/(dq1));
           tofp(p,q)= tofp(p,q)+tof_curr;
           actofp(p,q)=actofp(p,q)+(tof_curr*phi(p,q)*dv(p,q));
           [m_d2,n_d2,alpha0,beta0]=Exit_face(tof_curr,t11,t12,t21,t22,alpha,beta,m_d,n_d);
           singlecell_tofp_calculator(alpha0,beta0,m_d2,n_d2,qx_d1,qy_d1);
           break;
        
        end
end

