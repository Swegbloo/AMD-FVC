function out = A_fun(a,d,p,j,m,n,alp)
out = 0;
if n == 0
    eps = 1;
else
    eps = 2;
% disp(alp)
% disp(n)
% for i = 0:n
%      out = out + alp(i+1)*phi_prime(p,i,a,d)*cc_fun(n,m(i+1),d,k)/phi(p,i,a,d);
% end
% out = out + A_star(a,d,k,m(k+1));
% end
for i = 0:n
    out = out + eps/(pi*d)*phi_prime(p,i,a,d)/phi(p,i,a,d)*cc_fun(i,m,d,j)

