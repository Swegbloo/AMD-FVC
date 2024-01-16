function out = phi_prime(p,n,a,d)
out = 0;
if n == 0
    out = p*r^(p-1);
elseif n>0
    out = (n*pi/d)*0.5*(besseli(p-1,n*pi*a/d)+besseli(p+1,n*pi*a/d));
end
end