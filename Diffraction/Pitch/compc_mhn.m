function[z] = compc_mhn(Nn,mn,ik2,dn,option)

lambda = ik2*pi/dn;

if(option==1)
    
    if(lambda == 0)
        z = sinh(mn*dn)/mn;
        z = z*Nn^-0.5;
        
%         disp('active..');
    else
    
    z = mn*sinh(mn*dn);
    z = z*(-1.0)^ik2;
    z = z*Nn^-0.5;
    z = z/(lambda^2 + mn^2);
    
%     disp('active..')
    
    end
    
elseif(option==2)
    
    if(lambda==0)
        z = sin(mn*dn)/mn;
        z = z*Nn^-0.5;
    else
    
    z = mn*sin(mn*dn);
    z = z/(mn^2 - lambda^2);
    z = z*(-1.0)^ik2;
    z = z*Nn^-0.5;
    end
    
end