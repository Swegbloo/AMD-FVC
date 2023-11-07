function out = g_fun(k,n,m,d,l,a)
out = 0;
for i = 1:k
    out = out + (cc_fun(n,m(i,l),d,i-1)*R_fun(i-1,a,m(i,l))*A_star(m(i,l),d,a,i-1))/(m(i,l)*R_prime(m(i,l),a,i-1));
end
out = out - alpha_star(n,d,a);
end
