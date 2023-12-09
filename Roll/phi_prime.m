function out = phi_prime(n,a,d)
out = 0;
if n == 0
    out = 1;
elseif n>0
    out = n*pi*dbesseli(1,n*pi*a/d)/d;
end
end