function out = A_fun(a,d,k,m,n,alp)
out = 0;
for i = 0:n
    
    disp(alp(i+1))
    disp(phi_prime(i,a,d))
    disp(cc_fun(i,m(k+1),d,k))
    disp(phi(i,a,d))
    return
    out = out + alp(i+1)*phi_prime(i,a,d)*cc_fun(i,m(k+1),d,k)/phi(i,a,d);
end
out = out + A_star(a,d,k,m(k+1));
disp(out)
end