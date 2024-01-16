function out = g_fun_roll(k,n,m,d,a)
out = 0;
eps=zeros(n);
if n == 0
    eps = 1;
elseif n>0
    eps = 2;
end
for i = 0:k
out = out + A_star(a,d,i,m(i+1))*R_ratio(m(i+1),a,i)*cc_fun(n,m(i+1),d,i);
%disp(A_star(a,d,i,m))
%disp(R_ratio(m,a,i))
%disp(cc_fun(n,m,d,i))

end
out = (eps/d)*out - alpha_star(n,d,a);

end
% if n == 0
% disp();
% end