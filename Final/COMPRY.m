function[z] = COMPRY(ik,a,alpha_0,alpha,option,mop)

%% orders of bessel functions is 0 for heave and 1 for roll/pitch

if(option==1)
    
    if(ik==1)
        
        z = besselh(mop,1,alpha_0*a); %/besselk(0,alpha_0*a)
        
%         z = besselk(0,alpha_0*a);
        
    else
        
        z = besselk(mop,alpha(ik-1,3)*a);
%         z = 1.0;
    end

%         z = 1.0;
    
elseif(option==2)
    
    if(ik == 1)
        
        z = dbesselh(mop,1,alpha_0*a);  %
%         alpha_0_new = alpha_0;
%         z = (besselk(mop,alpha_0_new*a)*dbesselh(mop,1,alpha_0*a) - ...
%             (alpha_0_new/alpha_0)*besselh(mop,1,alpha_0*a)*dbesselk(mop,alpha_0_new*a))/...
%             besselk(mop,alpha_0_new*a);

    else
        
        z = dbesselk(mop,alpha(ik-1,3)*a);    %
%         mk = alpha(ik-1,3);
%         z = besselk(mop,mk*a)*dbesselk(mop,mk*a) - 
        
    end

    
end

