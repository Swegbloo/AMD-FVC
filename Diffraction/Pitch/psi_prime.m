function out = psi_prime(p,n,a,d)
out = 0;
if n == 0
    out = p*a^(p-1);
elseif n>0
    out = (n*pi/d)^2*dbesseli(p,n*pi*a/d);
end
end
