function out = e_fun_sway(k,n,j,a,d,m,~)
out = 0;
for i = 1:k+1
    out = out + (cc_fun(n,m(i),d,i-1)*cc_fun(j,m(i),d,i-1)*R_ratio(m(i),a,i-1))/(m(i));
end
out = out*S_fun_sway(a,d,j,m(i));
% if j == 0
% disp(S_fun_sway(a,d,j,m));
% end
end