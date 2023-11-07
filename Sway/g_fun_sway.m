function out = g_fun_sway(k,n,m,d,a)
out = 0;
for i = 1:k+1
    out = out + (cc_fun(n,m(i),d,i-1)*R_ratio(m(i),a,i-1)*A_star(a,d,i-1,m(i)))/(m(i));
end
end