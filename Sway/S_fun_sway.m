function out = S_fun_sway(a,d,k,m)
if k == 0
    out = d/(4*a);
elseif k>= 0
    out = 0.5*k*pi*dbesseli(1,m*a)/besseli(1,m*a);
end