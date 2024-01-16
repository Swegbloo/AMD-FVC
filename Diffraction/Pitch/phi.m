function out = phi(p,n,a,d)
out = 0;
if n == 0
    out = a^p;
elseif n>0
    out = besseli(p,n*pi*a/d);
end
end