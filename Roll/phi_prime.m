function out = phi_prime(n,a,d)
out = 0;
if n == 0
    out = 1;
elseif n>0
    out = dbesseli(1,n*pi*a/d);
end
end