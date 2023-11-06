function out = e_fun(k,n,j,a,d,m,l)
out = 0;
for i = 1:k
    out = out + (cc_fun(n,m(i,l),d,i-1)*cc_fun(j,m(i,l),d,i-1)*R_fun(i-1,a,m(i,l)))/(m(i,l)*R_prime(m(i,l),a,i-1));
end
out = out*s_fun(j,a,d);
end