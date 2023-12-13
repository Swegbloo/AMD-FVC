function out = R_prime(m,r,k)
if k==0
%     out = -m*besselh(1,m*r);
    out = dbesselh(0,1,m*r);
elseif k>=1
%     out = -m*besselk(1,m*r);
    out = dbesselk(0,m*r);
end
end
