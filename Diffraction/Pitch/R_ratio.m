function out = R_ratio(p,m,r,k)
out = 0;
    if k==0
    %     out = -m*besselh(1,m*r);
        out = besselh(p,1,m*r)/(m*dbesselh(p,1,m*r));
    elseif k>=1
    %     out = -m*besselk(1,m*r);
        out = besselk(p,m*r)/(m*dbesselk(p,m*r));
    end
end

