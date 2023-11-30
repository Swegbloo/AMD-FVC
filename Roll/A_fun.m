function out = A_fun(a,d,k,m,n,alp)
out = 0;
for i = 0:n
    out = out + alp(i+1)*phi_prime(i,a,d)*cc_fun(n,m(i+1),d,k)/phi(i,a,d);
end
out = out + A_star(a,d,k,m(k+1));
end