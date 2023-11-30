function out = A_star(a,d,k,m)
out = 0;
if k==0
    out = -1/(2*d*N_const(m,k)^0.5)*((d^2*sinh(m*d)/m(k+1))-(2*d*cosh(m*d)/m^2)+(2*sinh(m*d)/m^3))+0.75*a^2/N_const(m,k)^0.5*sinh(m*d)/(2*d*m)+1/N_const(m,k)^0.5*((d*sinh(m*d)/m)-(sinh(m)/m)-(cosh(m*d)/m^2)+cosh(m)/m^2);
elseif k>0
    out = -1/(2*d*N_const(m,k)^0.5)*((d^2*sin(m*d)/m)+(2*d*cos(m*d)/m^2)-(2*sin(m*d)/m^3))+0.75*a^2/N_const(m,k)^0.5*sin(m*d)/(2*d*m)+1/N_const(m,k)^0.5*((d*sin(m*d)/m)-(sin(m)/m)+(cos(m*d)/m^2)-cos(m)/m^2);
end
end