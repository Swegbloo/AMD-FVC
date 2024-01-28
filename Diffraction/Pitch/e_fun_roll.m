function out = e_fun_roll(p,j,n,l,a,d,m)
out = 0;
if l == 0
    eps = 1;
elseif l > 0
    eps = 2;
end
for i = 0:j
    out = out + (eps/d)*phi_prime(p,l,a,d)*R_ratio(p,m(i+1),a,i)*cc_fun(l,m(i+1),d,i)*cc_fun(n,m(i+1),d,i)/phi(p,l,a,d);
end
end