function out = alpha_star(n,d,a)
if n~=0
    out = (2*d/(n*pi)^2)*(-1)^n;
elseif n==0
    out = d/3-a^2/(2*d);
end
end