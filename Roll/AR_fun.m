function out = AR_fun(a,d,k,m,n,alp)
out = 0;

%disp(alpha_matrix);

for i = 0:n
    out = out + S_fun_roll(a,d,i,m(i+1))*cc_fun(i,m(i+1),d,k)*alp(i+1,1);
end

out = out + A_star(a,d,k,m(k+1));
out = out*R_ratio(m(k+1),a,k)/m(k+1);
end