function out = e_fun_roll(p,k,n,j,a,d,m,~)
out = 0;
if n == 0
    eps = 1;
elseif n > 0
    eps = 2;
end
for i = 0:k
    out = out + (eps/d)*phi_prime(j,a,d)*R_ratio(m(i+1),a,i)*cc_fun(j,m(i+1),d,i)*cc_fun(n,m(i+1),d,i)/phi(j,a,d);
end
end