function out = Z_star(m,k,d)
out = 0;
if k == 0
    out = (1/N_const(m,k))*(sinh(m)/m-d*sinh(m*d)/m-cosh(m)/m^2+cosh(m*d)/m^2);
elseif k>0
    out = (1/N_const(m,k))*(sin(m)/m-d*sin(m*d)/m+cos(m)/m^2-cos(m*d)/m^2);
end
end