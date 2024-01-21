function out = A_fun(p,a,d,k,m,n,alp)
out = 0;
% disp(alp)
% disp(n)
for i = 0:n
    out = out + alp(i+1)*phi_prime(p,i,a,d)*cc_fun(n,m(i+1),d,k)/phi(p,i,a,d);
end
out = out + A_star(a,d,k,m(k+1));
end