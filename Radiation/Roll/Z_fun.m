function out = Z_fun(n,~)
if n == 0
    out = (cosh(m*d) + m*sinh(m) - cosh(m) - dn*m*sinh(m*d))/(sqrt(N_const(m,k))*m^2);
    % disp(out);
elseif n>0
    out = (cos(m) + m*sin(m) - d*m*sin(m*d) - cos(m*d))/(sqrt(N_const(m,k))*m^2);
    % disp(out);
end
end