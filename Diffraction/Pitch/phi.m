function out = psi(p,n,a,d)
out = 0;
if n == 0
    out = a^p;
elseif n>0
    out = (n*pi/d)*besseli(p,n*pi*a/d);
end
end
