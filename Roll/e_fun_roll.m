function out = e_fun_roll(k,n,j,a,d,m,~)
out = 0;
if n == 0
    eps = 1;
elseif n > 0
    eps = 2;
end
for i = 0:k
    out = out + (eps/d)*phi_prime(i,a,d)*R_ratio(m,a,i)*cc_fun(j,m,d,i)*cc_fun(n,m,d,k)/phi(i,a,d);
end
end