function out = S_fun_roll(a,d,n,~)
out = 0;
if n == 0
    out = 1/(a*d);
elseif n>0
    out = (2*n*pi/d)*dbesseli(1,n*pi*a/d)/(d*besseli(1,n*pi*a/d));
end