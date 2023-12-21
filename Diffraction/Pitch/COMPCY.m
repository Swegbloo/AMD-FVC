function[z] = COMPCY(dn,ik,in,alpha_0,N0,alpha)

if(ik==1)
    arg = alpha_0*dn;
    z = 2.0*sinh(arg)/arg;
    z = z*(-1)^in;
    
    z2 = 1.0 + (in*pi/arg)^2;
    z2 = z2*N0^0.5;
    z = z/z2;

%     z = 2.0*sqrt(N0)*((-1.0)^in)*sinh(alpha_0*dn)/...
%         (dn*alpha_0*(1.0+(in*pi/(alpha_0*dn))^2)); 

    
%     lambda = in*pi/dn;
%     z = ((-1.0)^in)*sinh(alpha_0*dn)/(alpha_0*sqrt(N0)*(1.0+(lambda/alpha_0)^2));
    
else
    
%     arg = alpha(ik-1,3)*dn + in*pi; 
%     
%     z = sin(arg)/arg;
%     
%     arg = in*pi -alpha(ik-1,3)*dn ;
%     
%     z = z + sin(arg)/arg;
%     
%     z = z*alpha(ik-1,4)^-0.5;
%     
%     z2 = alpha(ik-1,3)*dn/(alpha(ik-1,3)*dn + in*pi);
%     
%     z = z*z2;

    mk = alpha(ik-1,3);nk = alpha(ik-1,4);
%     
    arg1 = mk*dn - in*pi;
    arg2 = mk*dn + in*pi;
    
    z = 2.0*sin(arg1)/arg1; %
    z = z*mk*dn/arg2;
    z = z*nk^-0.5;

%     lambda = in*pi/dn;
%     z = ((-1.0)^in)*sin(mk*dn)/(mk*sqrt(nk)*(1.0-(lambda/mk)^2));
    
end