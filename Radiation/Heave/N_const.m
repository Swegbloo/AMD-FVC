function out = N_const(m,k)
if k == 0
    out = 0.5*(1+sinh(2*m)/(2*m));
else
    out = 0.5*(1+sin(2*m)/(2*m));
end
