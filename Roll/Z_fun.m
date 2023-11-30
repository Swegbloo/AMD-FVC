function out = Z_fun(n,a)
if n == 0
    out = 1;
elseif n>0
    out = dbesseli(1,n*pi*a/d);
end
end