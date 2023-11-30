function out = R_ratio(m,r,k)
out = 0;
if k==0
%     out = -m*besselh(1,m*r);
    out = besselh(1,m*r)/(dbesselh(1,m*r));
elseif k>=1
%     out = -m*besselk(1,m*r);
    out = besselk(1,m*r)/(dbesselk(1,m*r));
end
end
