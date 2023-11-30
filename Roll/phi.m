function out = phi(n,a,d)
out = 0;
if n == 0
    out = a;
elseif n>0
    out = besseli(1,n*pi*a/d);
end
end