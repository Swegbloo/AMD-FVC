function out = A_star(a,d,k,m)
out = 0;

if k==0
    %out = -1/(2*d*N_const(m,k)^0.5)*((d^2*sinh(m*d)/m(k+1))-(2*d*cosh(m*d)/m^2)+(2*sinh(m*d)/m^3))+0.75*a^2/N_const(m,k)^0.5*sinh(m*d)/(2*d*m)+1/N_const(m,k)^0.5*((d*sinh(m*d)/m)-(sinh(m)/m)-(cosh(m*d)/m^2)+cosh(m)/m^2);
    out = (16.0*m*d*cosh(m*d) + 8.0*m^2*d*sinh(m) - 8.0*m*d*cosh(m)...
            -(12.0*m^2*d^2 - 3.0*m^2*a^2 + 8.0)*sinh(m*d))/(8.0*m^3*d*sqrt(N_const(m,k)));
     % disp(out);
elseif k>0
    % arg = m*d;
    % Nk = 0.5*(1+sin(arg)/arg);
    % out = (a/(2.0*d))*(1.0/m)*Nk^-0.5;
    % out = -1.0*out*sin(m*d);
    out = (-16.0*m*d*cos(m*d) + 8.0*m^2*d*sin(m) + 8.0*m*d*cos(m)...
            -(12.0*m^2*d^2 - 3.0*m^2*a^2 - 8.0)*sin(m*d))/(8.0*m^3*d*sqrt(N_const(m,k)));
    % disp(out);
    % out = -1/(2*d*N_const(m,k)^0.5)*((d^2*sin(m*d)/m)+(2*d*cos(m*d)/m^2)-(2*sin(m*d)/m^3))+0.75*a^2/N_const(m,k)^0.5*sin(m*d)/(2*d*m)+1/N_const(m,k)^0.5*((d*sin(m*d)/m)-(sin(m)/m)+(cos(m*d)/m^2)-cos(m)/m^2);
end
end