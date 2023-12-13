function out = E_new_fun(m,k,d)
out = 0;
if k == 0
    out = N_const(m(k+1),k)^(-0.5)*((d-1)*sinh(m(k+1)*d)/m(k+1)+(cosh(m(k+1))-cosh(m(k+1)*d))/m(k+1)^2);
elseif k>0
    out = N_const(m(k+1),k)^(-0.5)*((d-1)*sin(m(k+1)*d)/m(k+1)-(cos(k+1)-cos(m(k+1)*d))/m(k+1)^2);
end
end