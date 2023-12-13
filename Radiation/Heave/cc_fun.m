function out = cc_fun(n,m,d,k)
if k == 0
    out = (2*(-1)^(n)*sinh(m*d)/(m*d))/((1+((n*pi)/(m*d))^2)*sqrt(N_const(m,k)));
elseif k >= 1
    out = (2*sin(m*d-n*pi)*m*d)/((m*d-n*pi)*(m*d+n*pi)*sqrt(N_const(m,k)));
end
end
