function out = R_fun(k,r,m) %inspection needed
if k == 0
    out = besselh(0,m*r);
elseif k >= 1
    out = besselk(0,m*r);
end

